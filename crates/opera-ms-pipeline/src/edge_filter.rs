//! Edge filter: removes long-read edges
//! whose support is low compared with the coverage of their contigs
//! (log(support / (min coverage + 0.0001)) below median - 1.5 * stddev,
//! or support 1), in each of the 6 edge sets of the long-read mapping.
//!
//! Each edge line loses its last character (`chop`)
//! and gets a newline back; edges with a contig missing from the coverage
//! file are dropped from both outputs; both output files are created even
//! when there's no edge. Errors stop the run (a missing edge file, a
//! malformed coverage file, a division by zero or log of 0).

use std::collections::HashMap;
use std::fs;

use anyhow::{Context, Result, bail, ensure};

use crate::stats;

const EDGE_STD_THRESHOLD: f64 = 1.5;
const EDGE_SETS: usize = 6;

/// Numeric value of a string as Perl reads it (leading number, else 0).
pub fn perl_number(text: &str) -> f64 {
    let text = text.trim_start();
    let bytes = text.as_bytes();
    let mut end = 0;
    let mut seen_digit = false;
    let mut seen_dot = false;
    let mut seen_exp = false;
    let mut i = 0;
    if i < bytes.len() && (bytes[i] == b'-' || bytes[i] == b'+') {
        i += 1;
    }
    while i < bytes.len() {
        let c = bytes[i];
        if c.is_ascii_digit() {
            seen_digit = true;
            end = i + 1;
        } else if c == b'.' && !seen_dot && !seen_exp {
            seen_dot = true;
        } else if (c == b'e' || c == b'E') && seen_digit && !seen_exp {
            seen_exp = true;
            if i + 1 < bytes.len() && (bytes[i + 1] == b'-' || bytes[i + 1] == b'+') {
                i += 1;
            }
        } else {
            break;
        }
        i += 1;
    }
    if !seen_digit {
        return 0.0;
    }
    text[..end].parse().unwrap_or(0.0)
}

/// `chop`: removes the last character.
fn chop(text: &str) -> &str {
    let mut chars = text.chars();
    chars.next_back();
    chars.as_str()
}

/// Coverage per contig: read count / (number of windows * window size).
pub fn contig_coverage(cov_file: &str) -> Result<HashMap<String, f64>> {
    let mut lines = cov_file.split_inclusive('\n');
    let header = chop(lines.next().context("empty coverage file")?);
    let window_size = perl_number(header.split(' ').nth(3).unwrap_or(""));
    let mut coverage = HashMap::new();
    while let Some(line) = lines.next() {
        let line = line.strip_suffix('\n').unwrap_or(line);
        let fields: Vec<&str> = line.split('\t').collect();
        let contig = fields.first().copied().unwrap_or("");
        let nb_window = perl_number(fields.get(4).copied().unwrap_or(""));
        let read_count = perl_number(chop(lines.next().unwrap_or("")));
        lines.next();
        let denominator = nb_window * window_size;
        ensure!(
            denominator != 0.0,
            "coverage file: contig {contig} has no window (division by zero)"
        );
        coverage.insert(contig.to_string(), read_count / denominator);
    }
    Ok(coverage)
}

fn ratio(support: f64, cov1: f64, cov2: f64) -> Result<f64> {
    let cov = if cov1 > cov2 { cov2 } else { cov1 };
    let value = support / (cov + 0.0001);
    if value <= 0.0 {
        bail!("Can't take log of {value} (edge support {support})");
    }
    Ok(value.ln())
}

/// Filters one edge set: returns (kept, below threshold) file contents.
pub fn filter_edges(edges_text: &str, coverage: &HashMap<String, f64>) -> Result<(String, String)> {
    let mut edges: Vec<(&str, f64, f64, f64)> = Vec::new();
    let mut ratios = Vec::new();
    for line in edges_text.split_inclusive('\n') {
        let edge = chop(line);
        let fields: Vec<&str> = edge.split('\t').collect();
        let (c1, c2) = (
            fields.first().copied().unwrap_or(""),
            fields.get(2).copied().unwrap_or(""),
        );
        let (Some(&cov1), Some(&cov2)) = (coverage.get(c1), coverage.get(c2)) else {
            opera_ms_core::note!("{c1} {c2} DOESN'T EXIST");
            continue;
        };
        let support = perl_number(fields.get(6).copied().unwrap_or(""));
        if support != 1.0 {
            ratios.push(ratio(support, cov1, cov2)?);
        }
        edges.push((edge, support, cov1, cov2));
    }
    let (mut kept, mut below) = (String::new(), String::new());
    if edges.is_empty() {
        return Ok((kept, below));
    }
    let cutoff = stats::median(&ratios) - EDGE_STD_THRESHOLD * stats::stddev(&ratios);
    for (edge, support, cov1, cov2) in edges {
        let out = if ratio(support, cov1, cov2)? > cutoff && support > 1.0 {
            &mut kept
        } else {
            &mut below
        };
        out.push_str(edge);
        out.push('\n');
    }
    Ok((kept, below))
}

/// Filters the edge sets of `mapping_dir` with the coverage of `cov_file`.
pub fn filter_edge_coverage_ratio(mapping_dir: &str, cov_file: &str) -> Result<()> {
    let text = fs::read_to_string(cov_file)
        .with_context(|| format!("The file {cov_file} does not exists"))?;
    let coverage = contig_coverage(&text).with_context(|| format!("in {cov_file}"))?;
    let mut above = String::new();
    for i in 0..EDGE_SETS {
        let old = format!("{mapping_dir}/prefiltered_pairedEdges_i{i}");
        let filtered = format!("{mapping_dir}/pairedEdges_i{i}");
        let below_file = format!("{mapping_dir}/below_thresh_pairedEdges_i{i}");
        fs::rename(&filtered, &old).with_context(|| format!("cannot move {filtered} to {old}"))?;
        let edges = fs::read_to_string(&old).with_context(|| format!("cannot read {old}"))?;
        let (kept, below) = filter_edges(&edges, &coverage).with_context(|| format!("in {old}"))?;
        fs::write(&filtered, &kept).with_context(|| format!("cannot write {filtered}"))?;
        fs::write(&below_file, &below).with_context(|| format!("cannot write {below_file}"))?;
        above.push_str(&kept);
    }
    let all = format!("{mapping_dir}/pairedEdges_above_thresh");
    fs::write(&all, above).with_context(|| format!("cannot write {all}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coverage_from_sigma_contigs_file() {
        let text =
            "1 500 80 340\nk99_1\t11099\t109\t10988\t32\n841\n21 20\nk99_2\t600\tx\tx\t1\n0\n1\n";
        let cov = contig_coverage(text).unwrap();
        assert_eq!(cov["k99_1"], 841.0 / (32.0 * 340.0));
        assert_eq!(cov["k99_2"], 0.0);
    }

    #[test]
    fn edges_split_by_threshold() {
        let cov: HashMap<String, f64> = [("a", 10.0), ("b", 20.0), ("c", 1.0)]
            .map(|(k, v)| (k.to_string(), v))
            .into_iter()
            .collect();
        // supports: a-b 50 (high), a-c 2 (low ratio?), b-c 1 (always below), x-a missing
        let edges = "a\t+\tb\t+\t100\t10\t50\n\
                     a\t+\tc\t+\t100\t10\t2\n\
                     b\t+\tc\t-\t100\t10\t1\n\
                     x\t+\ta\t+\t100\t10\t9\n";
        let (kept, below) = filter_edges(edges, &cov).unwrap();
        // ratios: ln(50/10.0001), ln(2/1.0001); median = mean of both,
        // stddev = half their difference: cutoff = mean - 0.75 * diff
        let r1 = (50.0f64 / 10.0001).ln();
        let r2 = (2.0f64 / 1.0001).ln();
        let cutoff = stats::median(&[r1, r2]) - 1.5 * stats::stddev(&[r1, r2]);
        assert!(r2 > cutoff && r1 > cutoff);
        assert_eq!(kept, "a\t+\tb\t+\t100\t10\t50\na\t+\tc\t+\t100\t10\t2\n");
        assert_eq!(below, "b\t+\tc\t-\t100\t10\t1\n");
    }

    #[test]
    fn perl_numbers() {
        assert_eq!(perl_number("32"), 32.0);
        assert_eq!(perl_number(" 12abc"), 12.0);
        assert_eq!(perl_number("1.5e3x"), 1500.0);
        assert_eq!(perl_number("abc"), 0.0);
        assert_eq!(perl_number("-2"), -2.0);
        assert_eq!(chop("ab\n"), "ab");
        assert_eq!(chop("ab"), "a");
    }
}
