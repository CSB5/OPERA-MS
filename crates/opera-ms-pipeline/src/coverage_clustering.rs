//! Species and strain sets (stage 6, first part): finds the species whose clusters add up to more
//! than 1 Mbp, writes `reference_length.dat` (species, reference genome
//! length, length of its contigs, reference genome) and, per species, the
//! files the strain clustering reads (`matrix`, `contigs_window_cov`,
//! `pairedEdges_i<n>/pairedEdges_i<n>`, `contigs.fa`, `wrong`); then the
//! clusters and contigs of the other species (`clusters_single_strain`,
//! `single_strain_species.fa` in the reference clustering directory).
//!
//! Species and genomes are taken in sorted order. The reference genome of a
//! species is the one most of its clusters have as first genome in
//! `cluster_species.dat` (ties: first in sorted order); plasmid genomes are
//! listed in `plasmids.txt`; contigs shorter than 500 bp, without coverage
//! windows, are left out of the strain sets and assembled with the other
//! clusters.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Read, Write};

use anyhow::{Context, Result, bail, ensure};

use crate::edge_filter::perl_number as perl_value;
use crate::files::perl_number as perl_format;
use crate::ref_clustering::species_name;

const MIN_CONTIG_LENGTH: f64 = 500.0;
const MIN_SPECIES_LENGTH: f64 = 1_000_000.0;
const DEFAULT_GENOME_LENGTH: &str = "3500000";
const EDGE_SETS: usize = 6;

/// Arguments of the coverage clustering.
pub struct Args<'a> {
    pub inter_dir: &'a str,
    pub out_dir: &'a str,
    pub contig_file: &'a str,
    pub ref_dir: &'a str,
    pub opera_ms_db: &'a str,
}

fn read(path: &str) -> Result<String> {
    fs::read_to_string(path).with_context(|| format!("cannot read {path}"))
}

fn write(path: &str, text: &str) -> Result<()> {
    fs::write(path, text).with_context(|| format!("cannot write {path}"))
}

/// Perl `chop`: removes the last character, whatever it is.
fn chop(text: &str) -> &str {
    let mut chars = text.chars();
    chars.next_back();
    chars.as_str()
}

/// The lines of a file without their `\n` (Perl `chomp`).
fn lines(text: &str) -> impl Iterator<Item = &str> {
    text.split_inclusive('\n')
        .map(|l| l.strip_suffix('\n').unwrap_or(l))
}

/// The first `coverage_estimation/contigs_*` file (Perl `glob`, sorted).
fn windows_file(inter_dir: &str) -> Result<String> {
    let dir = format!("{inter_dir}/coverage_estimation");
    let mut found = Vec::new();
    for entry in fs::read_dir(&dir).with_context(|| format!("cannot read {dir}"))? {
        let name = entry
            .with_context(|| format!("cannot read {dir}"))?
            .file_name();
        let name = name.to_string_lossy().into_owned();
        if name.starts_with("contigs_") {
            found.push(name);
        }
    }
    found.sort();
    match found.first() {
        Some(name) => Ok(format!("{dir}/{name}")),
        None => bail!("no {dir}/contigs_* file"),
    }
}

/// A contig of at least 500 bp: coverage and length (as written).
struct Contig {
    cov: f64,
    length: f64,
}

/// The contig lines of the windows file (each followed by two lines).
fn read_windows(text: &str) -> Result<HashMap<String, Contig>> {
    let mut contigs = HashMap::new();
    let mut it = text.split_inclusive('\n');
    let header = chop(it.next().unwrap_or(""));
    let window_size = perl_value(header.split(' ').nth(3).unwrap_or(""));
    while let Some(line) = it.next() {
        let line = line.strip_suffix('\n').unwrap_or(line);
        let fields: Vec<&str> = line.split('\t').collect();
        let get = |i: usize| fields.get(i).copied().unwrap_or("");
        let read_count = perl_value(chop(it.next().unwrap_or("")));
        it.next();
        let length = perl_value(get(1));
        if length >= MIN_CONTIG_LENGTH {
            let denominator = perl_value(get(4)) * window_size;
            ensure!(
                denominator != 0.0,
                "contig {}: no window (Illegal division by zero)",
                get(0)
            );
            contigs.insert(
                get(0).to_string(),
                Contig {
                    cov: read_count / denominator,
                    length,
                },
            );
        }
    }
    Ok(contigs)
}

/// `genomes_length.txt`: genome -> length (as written).
fn genome_lengths(text: &str) -> HashMap<String, String> {
    let mut lengths = HashMap::new();
    for line in text.split_inclusive('\n') {
        let fields: Vec<&str> = chop(line).split('\t').collect();
        lengths.insert(
            fields.first().copied().unwrap_or("").to_string(),
            fields.get(2).copied().unwrap_or("").to_string(),
        );
    }
    lengths
}

/// The first genome listed for each cluster of `cluster_species.dat`, in
/// file order (clusters without a genome are skipped).
fn first_genomes(text: &str) -> Vec<(String, String)> {
    let mut result = Vec::new();
    let mut it = lines(text);
    let after_marker = |l: &str| l.find('>').map(|i| l[i + 1..].to_string());
    while let Some(line) = it.next() {
        let Some(mut cluster) = after_marker(line) else {
            continue;
        };
        let mut strain = it.next().unwrap_or("");
        while let Some(c) = after_marker(strain) {
            cluster = c;
            strain = it.next().unwrap_or("");
        }
        if !strain.is_empty() {
            result.push((cluster, strain.to_string()));
        }
    }
    result
}

/// Calls `f(header, sequence, is_last)` for each record of a FASTA file the
/// way OPERA-MS 0.9 read them: a header is any line containing `>`, the
/// sequence lines are joined; the last call is always made, with `None`
/// as header if the file has no header.
fn for_each_record(
    path: &str,
    mut f: impl FnMut(Option<&str>, &str, bool) -> Result<()>,
) -> Result<()> {
    let mut reader =
        BufReader::new(File::open(path).with_context(|| format!("File {path} not found"))?);
    let mut header: Option<String> = None;
    let mut seq = String::new();
    let mut buf = Vec::new();
    loop {
        buf.clear();
        if reader
            .read_until(b'\n', &mut buf)
            .with_context(|| format!("cannot read {path}"))?
            == 0
        {
            break;
        }
        let line = String::from_utf8_lossy(&buf);
        let line = line.strip_suffix('\n').unwrap_or(&line);
        if line.contains('>') {
            if let Some(h) = &header {
                f(Some(h), &seq, false)?;
            }
            header = Some(line.to_string());
            seq.clear();
        } else {
            seq.push_str(line);
        }
    }
    f(header.as_deref(), &seq, true)
}

/// `substr((split(/\s/, $header))[0], 1)`: the contig name of a header.
fn record_name(header: &str) -> String {
    header
        .split(char::is_whitespace)
        .next()
        .unwrap_or("")
        .chars()
        .skip(1)
        .collect()
}

/// The edges between the contigs of a matrix, given as (contig, coverage)
/// lines.
fn edges_from_matrix(
    inter_dir: &str,
    matrix: &[(String, String)],
    out_dir: &str,
    contig_file: &str,
) -> Result<()> {
    let species: HashSet<&str> = matrix
        .iter()
        .filter(|(_, cov)| {
            let cov = perl_value(cov);
            cov > 0.0 && cov < 1000.0
        })
        .map(|(c, _)| c.as_str())
        .collect();
    let mut sets: [String; EDGE_SETS] = Default::default();
    for i in 0..EDGE_SETS {
        let dir = format!("{out_dir}/pairedEdges_i{i}");
        fs::create_dir_all(&dir).with_context(|| format!("cannot create {dir}"))?;
    }

    let windows = read(&windows_file(inter_dir)?)?;
    let mut window_cov = String::new();
    let mut it = windows.split_inclusive('\n');
    window_cov.push_str(it.next().unwrap_or(""));
    while let Some(line) = it.next() {
        let body = line.strip_suffix('\n').unwrap_or(line);
        let contig = body.split('\t').next().unwrap_or("");
        let (l1, l2) = (it.next().unwrap_or(""), it.next().unwrap_or(""));
        if species.contains(contig) {
            window_cov.push_str(body);
            window_cov.push('\n');
            window_cov.push_str(l1);
            window_cov.push_str(l2);
        }
    }
    write(&format!("{out_dir}/contigs_window_cov"), &window_cov)?;

    let mut wrong = String::new();
    let edges = read(&format!("{inter_dir}/read_mapping/pairedEdges"))?;
    for line in lines(&edges) {
        let fields: Vec<&str> = line.split('\t').collect();
        let get = |i: usize| fields.get(i).copied().unwrap_or("");
        if perl_value(get(6)) < 2.0 {
            continue;
        }
        if species.contains(get(0)) && species.contains(get(2)) {
            let length = perl_value(get(4));
            if let Some(set) = [300.0, 1000.0, 2000.0, 5000.0, 15000.0, 40000.0]
                .iter()
                .position(|&limit| length < limit)
            {
                sets[set].push_str(line);
                sets[set].push('\n');
            }
        } else {
            wrong.push_str(line);
            wrong.push('\n');
        }
    }
    for (i, text) in sets.iter().enumerate() {
        write(
            &format!("{out_dir}/pairedEdges_i{i}/pairedEdges_i{i}"),
            text,
        )?;
    }
    write(&format!("{out_dir}/wrong"), &wrong)?;

    let path = format!("{out_dir}/contigs.fa");
    let mut out = std::io::BufWriter::new(
        File::create(&path).with_context(|| format!("cannot write {path}"))?,
    );
    for_each_record(contig_file, |header, seq, _| {
        if let Some(header) = header
            && species.contains(record_name(header).as_str())
        {
            write!(out, "{header}\n{seq}\n").with_context(|| format!("cannot write {path}"))?;
        }
        Ok(())
    })?;
    out.flush().with_context(|| format!("cannot write {path}"))
}

/// The species and strain sets.
pub fn coverage_clustering(a: &Args<'_>) -> Result<()> {
    let log_path = format!("{}/coverage_clustering.err", a.out_dir);
    let mut log = File::create(&log_path).with_context(|| format!("cannot write {log_path}"))?;
    let mut note = |text: String| -> Result<()> {
        log.write_all(text.as_bytes())
            .with_context(|| format!("cannot write {log_path}"))
    };
    let mut lengths_text = String::new();
    File::open(format!("{}/genomes_length.txt", a.opera_ms_db))
        .and_then(|mut f| f.read_to_string(&mut lengths_text))
        .with_context(|| format!("cannot read {}/genomes_length.txt", a.opera_ms_db))?;
    let genome_length = genome_lengths(&lengths_text);
    let plasmids_text =
        fs::read_to_string(format!("{}/plasmids.txt", a.opera_ms_db)).unwrap_or_default();
    let plasmids: HashSet<&str> = plasmids_text.lines().collect();
    let contigs = read_windows(&read(&windows_file(a.inter_dir)?)?)?;

    let ref_clusters_path = format!("{}/clusters_seq_similarity", a.ref_dir);
    let ref_clusters = read(&ref_clusters_path)?;
    let mut cluster_lengths: HashMap<&str, f64> = HashMap::new();
    for line in lines(&ref_clusters) {
        let mut fields = line.split('\t');
        let contig = fields.next().unwrap_or("");
        let cluster = fields.next().unwrap_or("");
        *cluster_lengths.entry(cluster).or_insert(0.0) +=
            contigs.get(contig).map_or(0.0, |c| c.length);
    }

    // votes (every genome gets 1), species length and clusters
    let mut genomes: BTreeMap<String, BTreeMap<String, usize>> = BTreeMap::new();
    let mut species_length: HashMap<String, f64> = HashMap::new();
    let mut species_clusters: HashMap<String, HashSet<String>> = HashMap::new();
    let species_file = read(&format!("{}/cluster_species.dat", a.ref_dir))?;
    for (cluster, strain) in first_genomes(&species_file) {
        let species = species_name(&strain).to_string();
        let Some(length) = cluster_lengths.get(cluster.as_str()) else {
            continue;
        };
        let plasmid = if strain.contains("XXX") {
            false
        } else {
            // the genomes whose first header names a plasmid, listed by
            // references::prepare
            plasmids.contains(strain.as_str())
        };
        if !plasmid {
            *genomes
                .entry(species.clone())
                .or_default()
                .entry(strain.clone())
                .or_insert(0) += 1;
        }
        *species_length.entry(species.clone()).or_insert(0.0) += length;
        species_clusters.entry(species).or_default().insert(cluster);
    }

    let mut reference_length = String::new();
    let mut to_analyze: BTreeSet<String> = BTreeSet::new();
    for (species, strains) in &genomes {
        // the most voted genome; ties: the first in sorted order
        let best = strains
            .iter()
            .filter(|(g, _)| !g.contains("XXX"))
            .fold(None::<(&String, usize)>, |acc, (g, &n)| match acc {
                Some((_, m)) if m >= n => acc,
                _ => Some((g, n)),
            })
            .map_or("NA", |(g, _)| g.as_str());
        let seq_length = if best == "NA" {
            DEFAULT_GENOME_LENGTH
        } else {
            note(format!(" *** get_ref_genome_length\t{best}\n"))?;
            genome_length.get(best).map_or("", String::as_str)
        };
        let length = species_length.get(species).copied().unwrap_or(0.0);
        let line = format!("{species}\t{seq_length}\t{}\t{best}\n", perl_format(length));
        if length > MIN_SPECIES_LENGTH {
            reference_length.push_str(&line);
            to_analyze.insert(species.clone());
        } else {
            note(line)?;
        }
    }
    write(
        &format!("{}/reference_length.dat", a.out_dir),
        &reference_length,
    )?;

    let mut strain_contigs: HashSet<&str> = HashSet::new();
    for species in &to_analyze {
        let dir = format!("{}/{species}", a.out_dir);
        fs::create_dir_all(&dir).with_context(|| format!("cannot create {dir}"))?;
        let clusters = species_clusters.get(species);
        let mut matrix = Vec::new();
        for line in lines(&ref_clusters) {
            let mut fields = line.split('\t');
            let contig = fields.next().unwrap_or("");
            let cluster = fields.next().unwrap_or("");
            if clusters.is_some_and(|c| c.contains(cluster)) {
                // contigs without coverage windows (< 500 bp) stay out of
                // the strain sets
                if let Some(c) = contigs.get(contig) {
                    matrix.push((contig.to_string(), perl_format(c.cov)));
                    strain_contigs.insert(contig);
                }
            }
        }
        let text: String = matrix
            .iter()
            .map(|(c, cov)| format!("{c}\t{cov}\n"))
            .collect();
        write(&format!("{dir}/matrix"), &text)?;
        edges_from_matrix(a.inter_dir, &matrix, &dir, a.contig_file)
            .with_context(|| format!("Error during get_edges_from_matrix for {species}"))?;
    }

    let single: String = lines(&ref_clusters)
        .filter(|l| !strain_contigs.contains(l.split('\t').next().unwrap_or("")))
        .map(|l| format!("{l}\n"))
        .collect();
    write(&format!("{}/clusters_single_strain", a.ref_dir), &single)?;

    let path = format!("{}/single_strain_species.fa", a.ref_dir);
    let mut out = std::io::BufWriter::new(
        File::create(&path).with_context(|| format!("cannot write {path}"))?,
    );
    for_each_record(a.contig_file, |header, seq, _last| {
        let name = record_name(header.unwrap_or(""));
        if !strain_contigs.contains(name.as_str()) {
            write!(out, "{}\n{seq}\n", header.unwrap_or(""))
                .with_context(|| format!("cannot write {path}"))?;
        }
        Ok(())
    })?;
    out.flush().with_context(|| format!("cannot write {path}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn first_genome_of_each_cluster() {
        let text = ">0S\ng/A__1\ng/A__2\n>1\n>2\ng/B__1\n>3\n";
        assert_eq!(
            first_genomes(text),
            [
                ("0S".to_string(), "g/A__1".to_string()),
                ("2".to_string(), "g/B__1".to_string())
            ]
        );
    }

    #[test]
    fn windows_skip_short_contigs() {
        let text = "1 500 80 340\na\t1000\tx\tx\t2\n680\n1 1\nb\t400\tx\tx\t1\n9\n1\n";
        let contigs = read_windows(text).unwrap();
        assert_eq!(contigs["a"].cov, 1.0);
        assert!(!contigs.contains_key("b"));
    }

    #[test]
    fn records_and_names() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("c.fa");
        fs::write(&path, ">a x\nAC\nGT\n>b\nTT\n").unwrap();
        let mut seen = Vec::new();
        for_each_record(path.to_str().unwrap(), |h, s, last| {
            seen.push((h.map(record_name), s.to_string(), last));
            Ok(())
        })
        .unwrap();
        assert_eq!(
            seen,
            [
                (Some("a".to_string()), "ACGT".to_string(), false),
                (Some("b".to_string()), "TT".to_string(), true)
            ]
        );
    }
}
