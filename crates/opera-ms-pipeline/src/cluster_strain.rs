//! Strain clustering (stage 6, per species): finds the
//! coverage modes of the species' windows, assigns contigs to strains, and
//! writes each strain's OPERA-LG input (`STRAIN_<n>/opera.config`,
//! `contigs.fa`, `pairedEdges_i<n>`), `strain_cluster.dat`,
//! `window_distibution.dat` and `excluded_contigs.fa`. OPERA-LG itself is
//! run by `run.rs` (timeout and fallback).
//!
//! Density modes and negative binomial quantiles come from
//! [`crate::rstats`]; values pass through R's printing (15 significant
//! digits) as in OPERA-MS 0.9, which called R.
//!
//! A mode is accepted when the windows of the contigs whose median falls
//! in its interval cover more than 500 kb (100 kb once a strain exists)
//! with a mean within 30% of the estimate, then kept if its re-estimated
//! interval still holds more than 500 kb of contigs for the first strain,
//! 100 kb for the next ones; when no mode is accepted every contig goes to
//! strain 1. Contigs are taken in sorted order.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::fs;
use std::path::Path;

use anyhow::{Context, Result, bail, ensure};

use crate::edge_filter::perl_number as perl_value;
use crate::files::perl_number as perl_format;
use crate::{rstats, stats};

const CLUSTER_SIZE_THRESHOLD_LOW: f64 = 100_000.0;
const CLUSTER_SIZE_THRESHOLD_HIGH: f64 = 500_000.0;
const MAX_CONTIG_COV: f64 = 100_000.0;
const EDGE_MEANS: [u32; 6] = [300, 1000, 2000, 5000, 15000, 40000];

/// The distance set of an edge line (`pairedEdges` format): the first whose
/// mean is above the edge's distance (field 5), `None` from 40 kb.
pub fn edge_set(line: &str) -> Option<usize> {
    let length = perl_value(line.split('\t').nth(4).unwrap_or(""));
    EDGE_MEANS.iter().position(|&mean| length < f64::from(mean))
}

/// The edge files of a strain directory, given to OPERA-LG (they can be
/// rebuilt from the species' `pairedEdges.gz`).
fn strain_edge_files(dir: &str) -> Vec<String> {
    (0..EDGE_MEANS.len())
        .map(|i| format!("{dir}/pairedEdges_i{i}"))
        .collect()
}

/// A contig of the species (`contigs_window_cov`).
#[derive(Debug, Clone)]
pub struct Contig {
    length: f64,
    /// The window read counts as written.
    windows: Vec<String>,
    window_values: Vec<f64>,
    mean_read_count: f64,
    median_read_count: f64,
    /// 0: not assigned.
    strain: usize,
}

/// The contigs of a species, by name.
#[derive(Clone)]
pub struct Species {
    pub contigs: BTreeMap<String, Contig>,
    pub window_size: f64,
}

/// The text without its last character (Perl's `chop`).
fn chop(text: &str) -> &str {
    let mut chars = text.chars();
    chars.next_back();
    chars.as_str()
}

/// `split(/ /, $text)`: Perl drops the trailing empty fields.
fn split_space(text: &str) -> Vec<&str> {
    let mut fields: Vec<&str> = text.split(' ').collect();
    while fields.last().is_some_and(|f| f.is_empty()) {
        fields.pop();
    }
    fields
}

/// `contigs_window_cov`: contig line, read count line, windows line.
pub fn read_species(text: &str) -> Species {
    let mut it = text.split_inclusive('\n');
    let header = chop(it.next().unwrap_or(""));
    let window_size = perl_value(split_space(header).get(3).copied().unwrap_or(""));
    let mut contigs = BTreeMap::new();
    while let Some(line) = it.next() {
        let line = line.strip_suffix('\n').unwrap_or(line);
        let fields: Vec<&str> = line.split('\t').collect();
        it.next(); // read count
        let windows: Vec<String> = split_space(chop(it.next().unwrap_or("")))
            .into_iter()
            .map(str::to_string)
            .collect();
        let window_values: Vec<f64> = windows.iter().map(|w| perl_value(w)).collect();
        contigs.insert(
            fields.first().copied().unwrap_or("").to_string(),
            Contig {
                length: perl_value(fields.get(1).copied().unwrap_or("")),
                mean_read_count: stats::mean(&window_values).trunc(),
                median_read_count: stats::median(&window_values),
                windows,
                window_values,
                strain: 0,
            },
        );
    }
    Species {
        contigs,
        window_size,
    }
}

/// The values the density is computed on: the windows of the contigs by
/// decreasing mean read count (ties by name), without the contigs above
/// 100,000.
pub fn mode_input(species: &Species) -> Vec<f64> {
    let mut contigs: Vec<&Contig> = species.contigs.values().collect();
    contigs.sort_by(|a, b| b.mean_read_count.total_cmp(&a.mean_read_count));
    contigs
        .into_iter()
        .filter(|c| c.mean_read_count <= MAX_CONTIG_COV)
        .flat_map(|c| c.window_values.iter().copied())
        .collect()
}

/// The tokens of the mode search: the number of windows, then the modes
/// as R prints them (`numeric(0)` when none).
pub fn mode_tokens(values: &[f64]) -> Result<Vec<String>> {
    ensure!(
        values.len() >= 2,
        "need at least 2 windows to estimate the coverage density (R: bw.nrd0)"
    );
    let modes = rstats::density_modes(values);
    // the number of windows is not a coverage mode
    let mut tokens = Vec::new();
    if modes.is_empty() {
        tokens.push("numeric(0)".to_string());
    } else {
        tokens.extend(rstats::r_printed(&modes));
    }
    Ok(tokens)
}

/// A number printed with `%.15g` and read back.
fn through_perl(value: f64) -> f64 {
    perl_format(value).parse().unwrap_or(value)
}

/// `compute_confidance_interval`: the 2% and 98% negative binomial
/// quantiles.
fn confidence_interval(mean: f64, dispersion: f64) -> (f64, f64) {
    let (size, mu) = (through_perl(dispersion), through_perl(mean));
    (
        rstats::qnbinom(0.02, size, mu),
        rstats::qnbinom(0.98, size, mu),
    )
}

fn estimate_mean(mode: f64, dispersion: f64) -> f64 {
    mode * (dispersion / (dispersion - 1.0))
}

/// The mode search and strain assignment (`while($mode_accepted == -1)`):
/// returns the strain list; `log` gets the messages.
pub fn assign_strains(
    species: &mut Species,
    tokens: &[String],
    dispersion_text: &str,
    log: &mut String,
) -> Result<Vec<usize>> {
    let dispersion = perl_value(dispersion_text);
    let _ = writeln!(log, " *** MODE DETECTED {}", tokens.join(" "));
    let mut accepted: i64 = -1;
    let mut strains = Vec::new();
    for token in tokens.iter().rev() {
        if token.contains('[') {
            continue;
        }
        let mode = perl_value(token);
        let strain_mean = estimate_mean(mode, dispersion);
        let (lo, hi) = confidence_interval(strain_mean, dispersion);
        let _ = writeln!(
            log,
            "\nMODE_INFO\t{token}\t{dispersion_text}\t{}\t[{} {}]",
            perl_format(strain_mean),
            perl_format(lo),
            perl_format(hi)
        );
        let mut in_dist = Vec::new();
        for contig in species.contigs.values() {
            let median = contig.median_read_count;
            if contig.strain != 0 || !(lo < median && median < hi) {
                continue;
            }
            in_dist.extend(
                contig
                    .window_values
                    .iter()
                    .copied()
                    .filter(|&c| c < median * 3.0),
            );
        }
        if in_dist.is_empty() {
            continue;
        }
        let size = in_dist.len() as f64 * species.window_size;
        let mean = in_dist.iter().sum::<f64>() / in_dist.len() as f64;
        let variance = stats::variance(&in_dist);
        let ratio = mean / strain_mean;
        let _ = writeln!(
            log,
            "MODE_EVAL\t{token}\t{dispersion_text}\t{}\t{}\t{}\t{}\t{}",
            perl_format(strain_mean),
            perl_format(size),
            perl_format(mean),
            perl_format(ratio),
            stats::basic_string(variance)
        );
        let big_enough = (accepted != -1 && size > CLUSTER_SIZE_THRESHOLD_LOW)
            || size > CLUSTER_SIZE_THRESHOLD_HIGH;
        if 0.70 < ratio && ratio < 1.3 && big_enough {
            if accepted == -1 {
                accepted = 0;
            }
            accepted += 1;
            let id = usize::try_from(accepted).unwrap_or(0);
            if select_mode_contigs(species, id, token, mean, variance, log)? {
                strains.push(id);
            } else {
                if accepted == 1 {
                    accepted = 0;
                }
                accepted -= 1;
            }
        }
    }
    if accepted == -1 {
        for contig in species.contigs.values_mut() {
            contig.strain = 1;
        }
        strains.push(1);
    }
    let _ = writeln!(log, " *** End mode detection");
    Ok(strains)
}

/// `select_mode_contigs`: re-estimates the dispersion and the interval,
/// assigns the contigs; undone when they cover 100 kb or less.
fn select_mode_contigs(
    species: &mut Species,
    id: usize,
    mode_text: &str,
    mean: f64,
    variance: f64,
    log: &mut String,
) -> Result<bool> {
    let excess = variance - mean;
    if excess == 0.0 {
        bail!("strain {id}: window variance equals the mean (division by zero)");
    }
    let dispersion = if excess < 0.0 {
        500.0
    } else {
        mean * mean / excess
    };
    let strain_mean = estimate_mean(perl_value(mode_text), dispersion);
    let (lo, hi) = confidence_interval(strain_mean, dispersion);
    let _ = writeln!(
        log,
        "MODE_IMPROVEMENT\t{mode_text}\t{}\t{}\t[{} {}]",
        perl_format(strain_mean),
        perl_format(dispersion),
        perl_format(lo),
        perl_format(hi)
    );
    let mut size = 0.0;
    for contig in species.contigs.values_mut() {
        let median = contig.median_read_count;
        if contig.strain == 0 && lo < median && median < hi {
            contig.strain = id;
            size += contig.length;
        }
    }
    // 500 kb for the first strain, 100 kb for the next ones
    let threshold = if id == 1 {
        CLUSTER_SIZE_THRESHOLD_HIGH
    } else {
        CLUSTER_SIZE_THRESHOLD_LOW
    };
    if size > threshold {
        return Ok(true);
    }
    for contig in species.contigs.values_mut() {
        if contig.strain == id {
            contig.strain = 0;
        }
    }
    let _ = writeln!(log, "MODE_REJECTED\t{mode_text}\t{}", perl_format(size));
    Ok(false)
}

fn read_if_exists(path: &str) -> Result<String> {
    if Path::new(path).exists() {
        fs::read_to_string(path).with_context(|| format!("cannot read {path}"))
    } else {
        Ok(String::new())
    }
}

fn write(path: &str, text: &str) -> Result<()> {
    fs::write(path, text).with_context(|| format!("cannot write {path}"))
}

/// `extract_edge`: the edges of a strain (both contigs in it, or one in it
/// and the other in a strain with a smaller ID).
fn strain_edges(species: &Species, edges: &str, id: usize) -> String {
    let mut out = String::new();
    for line in edges.split_inclusive('\n') {
        let fields: Vec<&str> = line.split('\t').collect();
        let strain_of = |i: usize| {
            fields
                .get(i)
                .and_then(|c| species.contigs.get(*c))
                .map(|c| c.strain)
        };
        let (Some(s1), Some(s2)) = (strain_of(0), strain_of(2)) else {
            continue;
        };
        let both = s1 == id && s2 == id;
        let rescued = s1 * s2 != 0 && (s1 == id || s2 == id) && (s1 < id || s2 < id);
        if both || rescued {
            out.push_str(line);
        }
    }
    out
}

/// The first word after the first `>` of a header line, if the line has one.
fn header_name(line: &str) -> Option<&str> {
    let rest = &line[line.find('>')? + 1..];
    let rest = rest.strip_suffix('\n').unwrap_or(rest);
    Some(rest.split(char::is_whitespace).next().unwrap_or(""))
}

/// Writes the strain directories and the species files (everything the
/// strain clustering writes before running OPERA-LG, and `excluded_contigs.fa`).
pub fn write_strains(species_dir: &str, species: &Species, strains: &[usize]) -> Result<()> {
    let contigs_fa = read_if_exists(&format!("{species_dir}/contigs.fa"))?;
    let species_edges = crate::gz::read(&format!("{species_dir}/pairedEdges"))?;
    let mut edge_sets: [String; EDGE_MEANS.len()] = Default::default();
    for line in species_edges.split_inclusive('\n') {
        if let Some(set) = edge_set(line) {
            edge_sets[set].push_str(line);
        }
    }
    for &id in strains {
        let dir = format!("{species_dir}/STRAIN_{id}");
        fs::create_dir_all(&dir).with_context(|| format!("cannot create {dir}"))?;
        let mut config = format!(
            "output_folder={dir}\ncontig_file={dir}/contigs.fa\nkeep_repeat=yes\nfilter_repeat=no\n\
             cluster_threshold=1\ncluster_increased_step=5\nkmer=60\n"
        );
        let files = strain_edge_files(&dir);
        for ((edge_file, mean), edges) in files.iter().zip(EDGE_MEANS).zip(&edge_sets) {
            write(edge_file, &strain_edges(species, edges, id))?;
            let _ = write!(
                config,
                "[LIB]\nmap_type=opera\nmap_file={edge_file}\nlib_mean={mean}\nlib_std={}\n",
                perl_format(f64::from(mean) / 10.0)
            );
        }
        write(&format!("{dir}/opera.config"), &config)?;

        // strain contigs: a header and the line after it
        let mut out = String::new();
        let mut lines = contigs_fa.split_inclusive('\n');
        while let Some(line) = lines.next() {
            let Some(name) = header_name(line) else {
                continue;
            };
            let seq = lines.next().unwrap_or("");
            if species.contigs.get(name).is_some_and(|c| c.strain == id) {
                let _ = write!(out, ">{name}\n{seq}");
            }
        }
        write(&format!("{dir}/contigs.fa"), &out)?;
    }

    let (mut clusters, mut windows) = (String::new(), String::new());
    for (name, contig) in &species.contigs {
        let _ = writeln!(
            clusters,
            "{name}\t{}\t{}\t0",
            perl_format(contig.mean_read_count),
            contig.strain
        );
        windows.push_str(&contig.windows.join("\n"));
        windows.push('\n');
    }
    write(&format!("{species_dir}/strain_cluster.dat"), &clusters)?;
    write(&format!("{species_dir}/window_distibution.dat"), &windows)?;

    // get_contig_sequence: the contigs left without a strain
    let mut excluded = String::new();
    let mut print = false;
    for line in contigs_fa.split_inclusive('\n') {
        match header_name(line) {
            Some(name) => {
                print = species.contigs.get(name).is_some_and(|c| c.strain == 0);
                if print {
                    let _ = writeln!(excluded, ">{name}");
                }
            }
            None if print => excluded.push_str(line),
            None => {}
        }
    }
    write(&format!("{species_dir}/excluded_contigs.fa"), &excluded)
}

/// The strain clustering up to OPERA-LG: returns the strain IDs; the log is
/// written to `log_path`.
pub fn cluster_strain(species_dir: &str, dispersion: &str, log_path: &str) -> Result<Vec<usize>> {
    let path = format!("{species_dir}/contigs_window_cov");
    let mut species =
        read_species(&fs::read_to_string(&path).with_context(|| format!("cannot read {path}"))?);
    let mut log = String::from(" *** End window construction\n");
    let tokens = mode_tokens(&mode_input(&species))?;
    let result = assign_strains(&mut species, &tokens, dispersion, &mut log)
        .and_then(|strains| write_strains(species_dir, &species, &strains).map(|()| strains));
    write(log_path, &log)?;
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    fn species(spec: &[(&str, f64, &[f64])]) -> Species {
        let mut text = String::from("x x x 100\n");
        for (name, length, windows) in spec {
            let w: Vec<String> = windows.iter().map(|v| v.to_string()).collect();
            text.push_str(&format!("{name}\t{length}\tx\tx\t1\n0\n{} \n", w.join(" ")));
        }
        read_species(&text)
    }

    #[test]
    fn reads_windows_mean_and_median() {
        let s = species(&[("a", 1000.0, &[10.0, 20.0, 31.0, 40.0])]);
        let a = &s.contigs["a"];
        assert_eq!(a.mean_read_count, 25.0);
        assert_eq!(a.median_read_count, 25.5);
        assert_eq!(a.windows, ["10", "20", "31", "40"]);
        assert_eq!(s.window_size, 100.0);
    }

    #[test]
    fn no_mode_accepted_gives_one_strain() {
        let mut s = species(&[("a", 1000.0, &[10.0, 12.0]), ("b", 1000.0, &[11.0, 9.0])]);
        let tokens = vec!["4".to_string(), "10.5".to_string()];
        let mut log = String::new();
        let strains = assign_strains(&mut s, &tokens, "80", &mut log).unwrap();
        assert_eq!(strains, [1]);
        assert!(s.contigs.values().all(|c| c.strain == 1));
    }

    #[test]
    fn edges_of_a_strain() {
        let mut s = species(&[("a", 1.0, &[1.0]), ("b", 1.0, &[1.0]), ("c", 1.0, &[1.0])]);
        s.contigs.get_mut("a").unwrap().strain = 1;
        s.contigs.get_mut("b").unwrap().strain = 2;
        s.contigs.get_mut("c").unwrap().strain = 2;
        let edges = "a\t+\tb\t+\t1\t1\t5\nb\t+\tc\t+\t1\t1\t5\nb\t+\tx\t+\t1\t1\t5\n";
        assert_eq!(strain_edges(&s, edges, 1), "");
        // strain 2: b-c, and a-b rescued (a is in the smaller strain 1)
        assert_eq!(
            strain_edges(&s, edges, 2),
            "a\t+\tb\t+\t1\t1\t5\nb\t+\tc\t+\t1\t1\t5\n"
        );
    }

    /// Opt-in: `STRAIN_CASES=<file from tests/r_cases/strain_cases.R>`: the
    /// strain assignment with R's modes vs with ours.
    #[test]
    fn against_r_modes() {
        let Ok(path) = std::env::var("STRAIN_CASES") else {
            return;
        };
        let text = std::fs::read_to_string(path).unwrap();
        let (mut cases, mut same_modes, mut same_strains, mut multi) = (0, 0, 0, 0);
        let mut differences = Vec::new();
        let mut dispersion = String::new();
        let mut contigs = String::new();
        for line in text.lines() {
            let fields: Vec<&str> = line.split('\t').collect();
            match fields[0] {
                "case" => {
                    dispersion = fields[1].to_string();
                    contigs = String::from("x x x 340\n");
                }
                "contig" => {
                    contigs.push_str(&format!(
                        "{}\t{}\tx\tx\t1\n0\n{}\n",
                        fields[1], fields[2], fields[3]
                    ));
                }
                "rmodes" => {
                    cases += 1;
                    let species = read_species(&contigs);
                    let r_tokens: Vec<String> = fields[1].split(' ').map(str::to_string).collect();
                    let mine = mode_tokens(&mode_input(&species)).unwrap();
                    let r_values: Vec<&String> =
                        r_tokens.iter().filter(|t| !t.contains('[')).collect();
                    if r_values.iter().copied().eq(mine.iter()) {
                        same_modes += 1;
                    }
                    let assign = |tokens: &[String]| {
                        let mut s = species.clone();
                        let mut log = String::new();
                        let strains = assign_strains(&mut s, tokens, &dispersion, &mut log);
                        strains
                            .map(|strains| {
                                let ids: Vec<usize> =
                                    s.contigs.values().map(|c| c.strain).collect();
                                (strains, ids)
                            })
                            .map_err(|e| e.to_string())
                    };
                    let (a, b) = (assign(&r_tokens), assign(&mine));
                    if a.as_ref().is_ok_and(|(strains, _)| strains.len() > 1) {
                        multi += 1;
                    }
                    if a.as_ref().is_ok_and(|(strains, _)| strains.len() > 1) {
                        multi += 1;
                    }
                    if a == b {
                        same_strains += 1;
                    } else {
                        let summary = |r: &Result<(Vec<usize>, Vec<usize>), String>| match r {
                            Ok((strains, ids)) => {
                                let counts: Vec<String> = strains
                                    .iter()
                                    .map(|s| {
                                        format!("{s}:{}", ids.iter().filter(|i| *i == s).count())
                                    })
                                    .collect();
                                format!(
                                    "{} (unassigned {})",
                                    counts.join(" "),
                                    ids.iter().filter(|i| **i == 0).count()
                                )
                            }
                            Err(e) => e.clone(),
                        };
                        differences.push(format!(
                            "case {cases}: R modes -> {}; ours -> {}\n  R    {}\n  ours {}",
                            summary(&a),
                            summary(&b),
                            fields[1].chars().take(200).collect::<String>(),
                            mine.join(" ").chars().take(200).collect::<String>()
                        ));
                    }
                }
                _ => {}
            }
        }
        for d in differences.iter().take(15) {
            opera_ms_core::note!("{d}");
        }
        opera_ms_core::note!(
            "{cases} species ({multi} with several strains): same modes {same_modes}, same strain assignment {same_strains}"
        );
    }
}
