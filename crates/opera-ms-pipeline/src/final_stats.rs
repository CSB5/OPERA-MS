//! Assembly statistics: `contig_info.txt`, `scaffolds.scaf.cname` and
//! `cluster_info.txt`.
//!
//! A scaffold whose name is neither `opera_scaffold_<n>` nor
//! `strain<s>_opera_scaffold_<n>` gets a new number, above every number used
//! so far. An empty assembly or cluster gets "NA" for its Nx/Lx and longest
//! contig; so does a cluster without a contig of 500 bp or more, or a
//! scaffold of contig length 0. Clusters with the same N50 are listed by
//! name.

use std::collections::{BTreeMap, HashMap};
use std::fmt::Write as _;
use std::fs;
use std::path::Path;

use anyhow::{Context, Result};

use crate::edge_filter::perl_number as perl_value;
use crate::files::perl_number as perl_format;

fn read(path: &str) -> Result<String> {
    fs::read_to_string(path).with_context(|| format!("cannot read {path}"))
}

fn write(path: &str, text: &str) -> Result<()> {
    fs::write(path, text).with_context(|| format!("cannot write {path}"))
}

/// The text without its last character (Perl's `chop`).
fn chop(text: &str) -> &str {
    let mut chars = text.chars();
    chars.next_back();
    chars.as_str()
}

fn lines(text: &str) -> impl Iterator<Item = &str> {
    text.split_inclusive('\n')
        .map(|l| l.strip_suffix('\n').unwrap_or(l))
}

/// Arguments of `generate_assembly_stats`.
pub struct AssemblyStatsArgs<'a> {
    pub scaffold_file: &'a str,
    /// "NULL" for none.
    pub cluster_file: &'a str,
    pub read_size: &'a str,
    pub coverage_file: &'a str,
    pub long_read_coverage_file: &'a str,
    /// "NULL" for none.
    pub species_file: &'a str,
    pub strain_dir: &'a str,
    pub scaffold_seq_file: &'a str,
    pub outfile: &'a str,
}

#[derive(Default, Clone, Copy)]
struct Coverage {
    size: Option<f64>,
    short: Option<f64>,
    long: Option<f64>,
}

/// `read_coverage_file` + `read_long_read_coverage_file`.
fn read_coverages(windows: &str, long: &str, read_size: f64) -> HashMap<String, Coverage> {
    let mut covs: HashMap<String, Coverage> = HashMap::new();
    let mut it = windows.split_inclusive('\n');
    it.next();
    while let Some(line) = it.next() {
        let line = line.strip_suffix('\n').unwrap_or(line);
        let f: Vec<&str> = line.split('\t').collect();
        let read_count = perl_value(chop(it.next().unwrap_or("")));
        it.next();
        covs.insert(
            f.first().copied().unwrap_or("").to_string(),
            Coverage {
                size: Some(perl_value(f.get(1).copied().unwrap_or(""))),
                short: Some(read_count * (2.0 * read_size)),
                long: Some(0.0),
            },
        );
    }
    for line in long.split_inclusive('\n') {
        let f: Vec<&str> = chop(line).split('\t').collect();
        let get = |i: usize| f.get(i).copied().unwrap_or("");
        covs.entry(get(0).to_string()).or_default().long =
            Some(perl_value(get(1)) * perl_value(get(2)));
    }
    covs
}

/// `read_species_file`: the first genome of each cluster.
fn read_species(text: &str) -> HashMap<String, String> {
    let mut species = HashMap::new();
    let (mut cluster, mut best) = (String::new(), String::new());
    for line in lines(text) {
        if let Some(i) = line.find('>') {
            if !cluster.is_empty() && !best.is_empty() {
                species.insert(cluster.clone(), best.clone());
            }
            cluster = line[i + 1..].to_string();
            best.clear();
        } else if best.is_empty() {
            best = line.to_string();
        }
    }
    // the last cluster
    if !cluster.is_empty() && !best.is_empty() {
        species.insert(cluster, best);
    }
    species
}

/// Name and cluster of each scaffold (`analyse_scaff_seq`): returns the
/// map, the filled lengths and the `.cname` text.
/// Each scaffold's (contig name, cluster), its filled length, the `.cname` text.
type ScaffoldNames = (
    HashMap<String, (String, String)>,
    HashMap<String, i64>,
    String,
);

fn scaffold_names(scaf: &str, seqs: &str) -> ScaffoldNames {
    let mut names = HashMap::new();
    let mut lengths = HashMap::new();
    let mut cname = String::new();
    let mut seq_lines = seqs.split_inclusive('\n');
    let (mut id, mut cluster) = (String::new(), String::new());
    // the highest number used so far (for scaffolds of an unknown name)
    let mut max_id = 0.0_f64;
    for line in lines(scaf) {
        if !line.contains('>') {
            continue;
        }
        let scaffold = line.split('\t').next().unwrap_or("");
        let mut known = false;
        if let Some(i) = scaffold.find(">opera_scaffold_") {
            id = scaffold[i + ">opera_scaffold_".len()..].to_string();
            cluster = "NA".to_string();
            known = true;
        }
        if let Some(i) = scaffold.find(">strain") {
            let rest = &scaffold[i + ">strain".len()..];
            if let Some(j) = rest.rfind("_opera_scaffold_") {
                id = perl_format(perl_value(&id) + 1.0);
                cluster = format!("cluster_{}", &rest[..j]);
                known = true;
            }
        }
        if !known {
            // a new number
            id = perl_format(max_id + 1.0);
            cluster = "NA".to_string();
        }
        max_id = max_id.max(perl_value(&id));
        let name = format!("opera_contig_{id}");
        seq_lines.next();
        let length = seq_lines.next().map_or(-1, |s| s.len() as i64 - 1);
        lengths.insert(scaffold.to_string(), length);
        let _ = writeln!(cname, "{}\t{name}", scaffold.get(1..).unwrap_or(""));
        names.insert(scaffold.to_string(), (name, cluster.clone()));
    }
    (names, lengths, cname)
}

/// The number of `STRAIN_*/contigs.fa` of a species (`ls -l ... | wc -l`),
/// "NA" without a strain directory.
fn nb_strains(strain_dir: &str, species: &str) -> Result<String> {
    let dir = format!("{strain_dir}/{species}");
    if !Path::new(&dir).is_dir() {
        return Ok("NA".to_string());
    }
    let mut count = 0;
    for entry in fs::read_dir(&dir).with_context(|| format!("cannot list {dir}"))? {
        let name = entry
            .with_context(|| format!("cannot list {dir}"))?
            .file_name();
        let name = name.to_string_lossy();
        if name.starts_with("STRAIN_") && Path::new(&format!("{dir}/{name}/contigs.fa")).exists() {
            count += 1;
        }
    }
    Ok(count.to_string())
}

/// `contig_info.txt` and `scaffolds.scaf.cname`.
pub fn generate_assembly_stats(a: &AssemblyStatsArgs<'_>) -> Result<()> {
    let covs = read_coverages(
        &read(a.coverage_file)?,
        &read(a.long_read_coverage_file)?,
        perl_value(a.read_size),
    );
    let mut clusters: HashMap<String, String> = HashMap::new();
    if a.cluster_file != "NULL" {
        for line in lines(&read(a.cluster_file)?) {
            let mut f = line.split('\t');
            clusters.insert(
                f.next().unwrap_or("").to_string(),
                f.next().unwrap_or("").to_string(),
            );
        }
    }
    let species = if a.species_file != "NULL" {
        read_species(&read(a.species_file)?)
    } else {
        HashMap::new()
    };
    let scaf = read(a.scaffold_file)?;
    let (names, lengths, cname) = scaffold_names(&scaf, &read(a.scaffold_seq_file)?);
    write(&format!("{}.cname", a.scaffold_file), &cname)?;

    let mut out = String::from(
        "CONTIG_ID\tLENGTH\tSHORT_READ_COV\tLONG_READ_COV\tCLUSTER\tSPECIES\tNB_STRAIN\tREFERENCE_GENOME\n",
    );
    let mut strain_counts: HashMap<String, String> = HashMap::new();
    let mut scaffold: Option<String> = None;
    let mut cluster: Option<String> = None;
    let mut sums = [0.0_f64; 3];
    let mut flush = |scaffold: &str,
                     cluster: Option<&String>,
                     sums: [f64; 3],
                     out: &mut String|
     -> Result<()> {
        let id = scaffold.split('\t').next().unwrap_or("");
        let (name, cluster_id) = names
            .get(id)
            .map_or(("", ""), |(n, c)| (n.as_str(), c.as_str()));
        let length = lengths.get(id).map_or(String::new(), |l| l.to_string());
        out.push_str(name);
        match cluster {
            Some(_) if sums[0] == 0.0 => {
                // "NA" instead of a division by zero
                let _ = write!(out, "\t{length}\tNA\tNA\t");
            }
            Some(_) => {
                let _ = write!(
                    out,
                    "\t{length}\t{:.2}\t{:.2}\t",
                    sums[1] / sums[0],
                    sums[2] / sums[0]
                );
            }
            None => {
                let _ = write!(out, "\t{length}\tNA\tNA\t");
            }
        }
        match cluster.and_then(|c| species.get(c)) {
            Some(genome) => {
                let file = genome.split('/').rfind(|s| !s.is_empty()).unwrap_or("");
                let mut parts = file.split("__");
                let species_name = parts.next().unwrap_or("").to_string();
                let reference = parts.next().unwrap_or("");
                if !strain_counts.contains_key(&species_name) {
                    let n = nb_strains(a.strain_dir, &species_name)?;
                    strain_counts.insert(species_name.clone(), n);
                }
                let _ = write!(
                    out,
                    "{cluster_id}\t{species_name}\t{}\t{reference}",
                    strain_counts[&species_name]
                );
            }
            None => out.push_str("NA\tNA\tNA\tNA"),
        }
        out.push('\n');
        Ok(())
    };
    for line in lines(&scaf) {
        if line.contains('>') {
            if let Some(s) = scaffold.as_deref().filter(|s| !s.is_empty()) {
                flush(s, cluster.as_ref(), sums, &mut out)?;
            }
            scaffold = Some(line.to_string());
            cluster = None;
            sums = [0.0; 3];
        } else {
            let contig = line.split('\t').next().unwrap_or("");
            let Some(c) = clusters.get(contig) else {
                continue;
            };
            cluster = Some(c.clone());
            if let Some(cov) = covs.get(contig) {
                for (sum, value) in sums.iter_mut().zip([cov.size, cov.short, cov.long]) {
                    *sum += value.unwrap_or(0.0);
                }
            }
        }
    }
    flush(
        scaffold.as_deref().unwrap_or(""),
        cluster.as_ref(),
        sums,
        &mut out,
    )?;
    write(a.outfile, &out)
}

// ---- cluster_info ---------------------------------------------------------

const COLUMNS: [&str; 16] = [
    "CLUSTER_ID",
    "SHORT_READ_COV",
    "LONG_READ_COV",
    "SPECIES",
    "SIMILARITY",
    "CLUSTER_SIZE",
    "NB_CONTIG",
    "LONGEST_CONTIG",
    "N50",
    "L50",
    "N90",
    "L90",
    "SHORT_READ_NB_CONTIG",
    "SHORT_READ_LONGEST_CONTIG",
    "SHORT_READ_N50",
    "SHORT_READ_L50",
];

#[derive(Default)]
struct Bin {
    fields: HashMap<&'static str, String>,
    contig_sizes: Vec<String>,
    short_sizes: Vec<String>,
    cluster_size: f64,
    short_assembly_size: f64,
    short_cov: f64,
    long_cov: f64,
}

impl Bin {
    /// `init_bin`.
    fn new() -> Self {
        let mut b = Self::default();
        for key in [
            "SPECIES",
            "SIMILARITY",
            "LONGEST_CONTIG",
            "NB_CONTIG",
            "N50",
            "L50",
            "N90",
            "L90",
            "SHORT_READ_LONGEST_CONTIG",
            "SHORT_READ_NB_CONTIG",
            "SHORT_READ_N50",
            "SHORT_READ_L50",
        ] {
            b.fields.insert(key, "NA".to_string());
        }
        b
    }
}

/// `compute_Nx`: the position of the contig reaching x% of the length
/// (-1 when the length is 0).
fn compute_nx(x: f64, length: f64, sorted: &[String]) -> i64 {
    let target = length * (x / 100.0);
    let (mut total, mut i) = (0.0, 0usize);
    while total < target {
        total += sorted.get(i).map_or(0.0, |v| perl_value(v));
        i += 1;
        if i > sorted.len() && total < target {
            break;
        }
    }
    i as i64 - 1
}

/// Lx from an Nx position: "NA" for an empty assembly (-1).
fn lx(i: i64) -> String {
    if i < 0 {
        "NA".to_string()
    } else {
        (i + 1).to_string()
    }
}

/// The value at an Nx position: "NA" for -1 (an empty assembly) and outside
/// the list.
fn perl_index(values: &[String], i: i64) -> String {
    usize::try_from(i)
        .ok()
        .and_then(|i| values.get(i).cloned())
        .unwrap_or_else(|| "NA".to_string())
}

fn sorted_desc(values: &[String]) -> Vec<String> {
    let mut v = values.to_vec();
    v.sort_by(|a, b| perl_value(b).total_cmp(&perl_value(a)));
    v
}

/// GNU `sort -k3,3 -n` (C locale): field 3 as a plain number (no
/// exponent), then the whole line.
fn sort_numeric_field3(text: &str) -> Vec<&str> {
    let key = |l: &str| -> f64 {
        let field = l.split_whitespace().nth(2).unwrap_or("");
        let mut end = 0;
        for (i, c) in field.char_indices() {
            if c.is_ascii_digit() || c == '.' || (i == 0 && c == '-') {
                end = i + 1;
            } else {
                break;
            }
        }
        field[..end].parse().unwrap_or(0.0)
    };
    let mut v: Vec<(f64, &str)> = lines(text).map(|l| (key(l), l)).collect();
    v.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.as_bytes().cmp(b.1.as_bytes())));
    v.into_iter().map(|(_, l)| l).collect()
}

/// The text of
/// `cluster_info.txt`.
pub fn cluster_info(assembly_dir: &str, bin_dir: &str) -> Result<String> {
    // analyze_scaffold: the short-read contig lengths of each OPERA-MS contig
    let scaffold_file = format!("{assembly_dir}/intermediate_files/opera_long_read/scaffolds.scaf");
    let mut contig_names: HashMap<String, String> = HashMap::new();
    for line in read(&format!("{scaffold_file}.cname"))?.split_inclusive('\n') {
        let mut f = chop(line).split('\t');
        contig_names.insert(
            f.next().unwrap_or("").to_string(),
            f.next().unwrap_or("").to_string(),
        );
    }
    let mut short_lengths: HashMap<String, Vec<String>> = HashMap::new();
    let mut current = String::new();
    for line in read(&scaffold_file)?.split_inclusive('\n') {
        let f: Vec<&str> = chop(line).split('\t').collect();
        let first = f.first().copied().unwrap_or("");
        if let Some(i) = first.find('>') {
            current = contig_names
                .get(&first[i + 1..])
                .cloned()
                .unwrap_or_default();
            short_lengths.insert(current.clone(), Vec::new());
        } else {
            short_lengths
                .entry(current.clone())
                .or_default()
                .push(f.get(2).copied().unwrap_or("").to_string());
        }
    }

    // analyze_bin
    let mut bins: BTreeMap<String, Bin> = BTreeMap::new();
    let mut contig_to_bin: HashMap<String, String> = HashMap::new();
    for entry in fs::read_dir(bin_dir).with_context(|| format!("cannot list {bin_dir}"))? {
        let file = entry
            .with_context(|| format!("cannot list {bin_dir}"))?
            .file_name()
            .to_string_lossy()
            .into_owned();
        let Some(pos) = file.rfind(".fa").filter(|p| *p > 0) else {
            continue;
        };
        let id = file[..pos].to_string();
        bins.insert(id.clone(), Bin::new());
        for line in read(&format!("{bin_dir}/{file}"))?.split_inclusive('\n') {
            let name = chop(line);
            if let Some(i) = name.find('>').filter(|i| name.len() > i + 1) {
                contig_to_bin.insert(name[i + 1..].to_string(), id.clone());
            }
        }
    }

    // analyze_contig_info
    let info = read(&format!("{assembly_dir}/contig_info.txt"))?;
    for line in info.split_inclusive('\n').skip(1) {
        let f: Vec<&str> = chop(line).split('\t').collect();
        let get = |i: usize| f.get(i).copied().unwrap_or("");
        let contig = get(0);
        let size = perl_value(get(1));
        if size < 500.0 {
            continue;
        }
        let bin_name = contig_to_bin
            .get(contig)
            .cloned()
            .unwrap_or_else(|| "UN_BINNED".to_string());
        let bin = bins.entry(bin_name).or_default();
        bin.contig_sizes.push(get(1).to_string());
        bin.cluster_size += size;
        let shorts = short_lengths.get(contig).cloned().unwrap_or_default();
        bin.short_assembly_size += shorts.iter().map(|v| perl_value(v)).sum::<f64>();
        bin.short_sizes.extend(shorts);
        bin.short_cov += perl_value(get(2)) * size;
        bin.long_cov += perl_value(get(3)) * size;
    }
    for bin in bins.values_mut() {
        let f = &mut bin.fields;
        if bin.cluster_size == 0.0 {
            // no contig of 500 bp or more: "NA" values
            // instead of a division by zero
            for key in [
                "SHORT_READ_COV",
                "LONG_READ_COV",
                "LONGEST_CONTIG",
                "N50",
                "L50",
                "N90",
                "L90",
                "SHORT_READ_LONGEST_CONTIG",
                "SHORT_READ_N50",
                "SHORT_READ_L50",
            ] {
                f.insert(key, "NA".to_string());
            }
            f.insert("CLUSTER_SIZE", "0".to_string());
            f.insert("NB_CONTIG", bin.contig_sizes.len().to_string());
            f.insert("SHORT_READ_NB_CONTIG", bin.short_sizes.len().to_string());
            continue;
        }
        f.insert(
            "SHORT_READ_COV",
            format!("{:.2}", bin.short_cov / bin.cluster_size),
        );
        f.insert(
            "LONG_READ_COV",
            format!("{:.2}", bin.long_cov / bin.cluster_size),
        );
        f.insert("CLUSTER_SIZE", perl_format(bin.cluster_size));
        let sorted = sorted_desc(&bin.contig_sizes);
        let (n50, n90) = (
            compute_nx(50.0, bin.cluster_size, &sorted),
            compute_nx(90.0, bin.cluster_size, &sorted),
        );
        f.insert("LONGEST_CONTIG", perl_index(&sorted, 0));
        f.insert("NB_CONTIG", sorted.len().to_string());
        f.insert("N50", perl_index(&sorted, n50));
        f.insert("L50", lx(n50));
        f.insert("N90", perl_index(&sorted, n90));
        f.insert("L90", lx(n90));
        let sorted = sorted_desc(&bin.short_sizes);
        let n50 = compute_nx(50.0, bin.short_assembly_size, &sorted);
        f.insert("SHORT_READ_LONGEST_CONTIG", perl_index(&sorted, 0));
        f.insert("SHORT_READ_NB_CONTIG", sorted.len().to_string());
        f.insert("SHORT_READ_N50", perl_index(&sorted, n50));
        f.insert("SHORT_READ_L50", lx(n50));
    }

    // analyze_mash: the closest genome of each cluster
    let mash = read(&format!("{bin_dir}/mash.dist"))?;
    for line in sort_numeric_field3(&mash) {
        let f: Vec<&str> = line.split('\t').collect();
        let get = |i: usize| f.get(i).copied().unwrap_or("");
        let file = get(1).split('/').rfind(|s| !s.is_empty()).unwrap_or("");
        let parts: Vec<&str> = file.split('.').collect();
        let id = parts[..parts.len().saturating_sub(1)].join(".");
        // an unknown cluster is created here, printed with empty fields
        let bin = bins.entry(id).or_default();
        if bin.fields.get("SPECIES").map(String::as_str) != Some("NA") {
            continue;
        }
        // the accession for `<accession>.fna.gz` genomes
        let species = crate::ref_clustering::species_name(get(0)).to_string();
        if species.contains("multispecies") {
            continue;
        }
        bin.fields.insert("SPECIES", species);
        if get(2) != "NA" {
            bin.fields
                .insert("SIMILARITY", format!("{:.6}", 1.0 - perl_value(get(2))));
        }
    }

    // by decreasing N50, ties by name
    let unclustered = bins.entry("unclustered".to_string()).or_default();
    unclustered.fields.insert("SPECIES", "NA".to_string());
    unclustered.fields.insert("SIMILARITY", "NA".to_string());
    let mut order: Vec<(&String, &Bin)> = bins.iter().collect();
    order.sort_by(|(_, a), (_, b)| {
        let n50 = |b: &Bin| b.fields.get("N50").map_or(0.0, |v| perl_value(v));
        n50(b).total_cmp(&n50(a))
    });
    let mut out = COLUMNS.join("\t") + "\n";
    for (name, bin) in order {
        out.push_str(name);
        for column in &COLUMNS[1..] {
            out.push('\t');
            out.push_str(bin.fields.get(column).map_or("", String::as_str));
        }
        out.push('\n');
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn species_file_keeps_the_last_cluster() {
        let s = read_species(">0S\ng/a__x\ng/b__y\n>3\ng/c__z\n");
        assert_eq!(s["0S"], "g/a__x");
        // the last cluster
        assert_eq!(s["3"], "g/c__z");
    }

    #[test]
    fn nx_positions() {
        let sizes: Vec<String> = ["50", "30", "20"].iter().map(|s| s.to_string()).collect();
        assert_eq!(compute_nx(50.0, 100.0, &sizes), 0);
        assert_eq!(compute_nx(90.0, 100.0, &sizes), 2);
        assert_eq!(compute_nx(50.0, 0.0, &sizes), -1);
        // an empty assembly (-1): "NA"
        assert_eq!(perl_index(&sizes, -1), "NA");
        assert_eq!(perl_index(&[], 0), "NA");
        assert_eq!((lx(-1), lx(2)), ("NA".to_string(), "3".to_string()));
    }

    #[test]
    fn scaffold_numbering() {
        let scaf = ">opera_scaffold_1\tlength: 5\n c\n>opera_scaffold_2\tl\n>strain1_opera_scaffold_1\tl\n";
        let seqs = ">a\nACGT\n>b\nAC\n>c\nA\n";
        let (names, lengths, cname) = scaffold_names(scaf, seqs);
        assert_eq!(
            names[">strain1_opera_scaffold_1"],
            ("opera_contig_3".to_string(), "cluster_1".to_string())
        );
        assert_eq!(lengths[">opera_scaffold_2"], 2);
        assert!(cname.starts_with("opera_scaffold_1\topera_contig_1\n"));
    }
}
