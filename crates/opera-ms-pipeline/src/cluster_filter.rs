//! Repeat contigs of the clusters: finds repeat contigs (coverage
//! above `threshold` x the median coverage of their cluster) and removes,
//! from the 6 edge sets, the edges touching a repeat or a contig outside
//! the clusters. Writes `NO_REPEAT/filtered_pairedEdges_i<i>`,
//! `repeat_contigs_info.txt` and `clusters_cov`.
//!
//! The first line of each edge file is copied as is; medians are printed
//! the way Perl's `Statistics::Basic` prints them (2 decimals, thousands
//! separators); a contig missing from the coverage file counts as coverage
//! 0 and length 0. The two diagnostic files list the clusters by cluster ID.

use std::collections::{HashMap, HashSet};
use std::fs;

use anyhow::{Context, Result, ensure};

use crate::edge_filter::perl_number as perl_value;
use crate::files::perl_number as perl_format;
use crate::stats;

const EDGE_SETS: usize = 6;

fn chop(text: &str) -> &str {
    let mut chars = text.chars();
    chars.next_back();
    chars.as_str()
}

/// `split(/\s+/, $_)` with Perl's trailing-empty-field removal.
fn split_whitespace_perl(line: &str) -> Vec<&str> {
    let mut fields: Vec<&str> = line.split(|c: char| c.is_whitespace()).collect();
    // consecutive whitespace is one separator, except a leading empty field
    let leading = fields.first().is_some_and(|f| f.is_empty());
    fields.retain(|f| !f.is_empty());
    if leading && !fields.is_empty() {
        fields.insert(0, "");
    }
    fields
}

struct Contig<'a> {
    cov: f64,
    length: f64,
    length_text: &'a str,
}

struct Cluster<'a> {
    length: f64,
    contigs: Vec<(&'a str, f64)>,
}

/// What the filter computes from its two input files.
pub struct Filter<'a> {
    contigs: HashMap<&'a str, Contig<'a>>,
    clusters: Vec<(&'a str, Cluster<'a>)>,
    in_cluster: HashSet<&'a str>,
    repeats: HashSet<&'a str>,
}

impl<'a> Filter<'a> {
    pub fn new(coverage_file: &'a str, cluster_file: &'a str, threshold: f64) -> Result<Self> {
        let mut contigs = HashMap::new();
        let mut lines = coverage_file.split_inclusive('\n');
        let header = chop(lines.next().context("empty coverage file")?);
        let window_size = perl_value(header.split(' ').nth(3).unwrap_or(""));
        while let Some(line) = lines.next() {
            let line = line.strip_suffix('\n').unwrap_or(line);
            let fields: Vec<&str> = line.split('\t').collect();
            let contig = fields.first().copied().unwrap_or("");
            let length_text = fields.get(1).copied().unwrap_or("");
            let nb_window = perl_value(fields.get(4).copied().unwrap_or(""));
            let read_count = perl_value(chop(lines.next().unwrap_or("")));
            lines.next();
            let denominator = nb_window * window_size;
            ensure!(
                denominator != 0.0,
                "contig {contig}: no window (division by zero)"
            );
            contigs.insert(
                contig,
                Contig {
                    cov: read_count / denominator,
                    length: perl_value(length_text),
                    length_text,
                },
            );
        }

        let mut index: HashMap<&str, usize> = HashMap::new();
        let mut clusters: Vec<(&str, Cluster)> = Vec::new();
        let mut in_cluster = HashSet::new();
        for line in cluster_file.lines() {
            let fields: Vec<&str> = line.split('\t').collect();
            let contig = fields.first().copied().unwrap_or("");
            let cluster = fields.get(1).copied().unwrap_or("");
            let (cov, length) = contigs
                .get(contig)
                .map_or((0.0, 0.0), |c| (c.cov, c.length));
            let i = *index.entry(cluster).or_insert_with(|| {
                clusters.push((
                    cluster,
                    Cluster {
                        length: 0.0,
                        contigs: Vec::new(),
                    },
                ));
                clusters.len() - 1
            });
            clusters[i].1.contigs.push((contig, cov));
            clusters[i].1.length += length;
            in_cluster.insert(contig);
        }
        // clusters by ID (numeric when possible)
        clusters.sort_by(
            |(a, _), (b, _)| match (a.parse::<u64>(), b.parse::<u64>()) {
                (Ok(x), Ok(y)) => x.cmp(&y),
                _ => a.cmp(b),
            },
        );

        let mut repeats = HashSet::new();
        for (_, cluster) in &clusters {
            let covs: Vec<f64> = cluster.contigs.iter().map(|(_, cov)| *cov).collect();
            let median = stats::median(&covs);
            for (contig, cov) in &cluster.contigs {
                if *cov > median * threshold {
                    repeats.insert(*contig);
                }
            }
        }
        Ok(Self {
            contigs,
            clusters,
            in_cluster,
            repeats,
        })
    }

    pub fn nb_repeats(&self) -> usize {
        self.repeats.len()
    }

    /// One edge file filtered (its first line is always kept).
    pub fn filter_edges(&self, edges: &str) -> String {
        let mut out = String::new();
        for (i, line) in edges.split_inclusive('\n').enumerate() {
            if i == 0 {
                out.push_str(line);
                continue;
            }
            let fields = split_whitespace_perl(line);
            if fields.len() < 3 {
                continue;
            }
            let (c1, c2) = (fields[0], fields[2]);
            let drop = self.repeats.contains(c1)
                || self.repeats.contains(c2)
                || !self.in_cluster.contains(c1)
                || !self.in_cluster.contains(c2);
            if !drop {
                out.push_str(line);
            }
        }
        out
    }

    /// `repeat_contigs_info.txt` and `clusters_cov`.
    pub fn reports(&self) -> (String, String) {
        let (mut repeat_info, mut clusters_cov) = (String::new(), String::new());
        for (id, cluster) in &self.clusters {
            let covs: Vec<f64> = cluster.contigs.iter().map(|(_, cov)| *cov).collect();
            let median = stats::basic_string(stats::median(&covs));
            repeat_info.push_str(&format!(
                "{id} {median} {} {}\n",
                perl_format(cluster.length),
                cluster.contigs.len()
            ));
            clusters_cov.push_str(&format!("{id} {median}\n"));
            for (contig, _) in &cluster.contigs {
                if self.repeats.contains(contig) {
                    let (cov, length) =
                        self.contigs.get(contig).map_or(("".to_string(), ""), |c| {
                            (perl_format(c.cov), c.length_text)
                        });
                    repeat_info.push_str(&format!("\t{contig}\t{cov}\t{length}\n"));
                }
            }
        }
        (repeat_info, clusters_cov)
    }
}

/// The repeat filter of the edge sets, written to `dir_res`.
pub fn filter_cluster_coverage(
    cov_file: &str,
    cluster_file: &str,
    edge_file_dir: &str,
    threshold: &str,
    dir_res: &str,
) -> Result<()> {
    let coverage =
        fs::read_to_string(cov_file).with_context(|| format!("cannot read {cov_file}"))?;
    let clusters =
        fs::read_to_string(cluster_file).with_context(|| format!("cannot read {cluster_file}"))?;
    let filter = Filter::new(&coverage, &clusters, perl_value(threshold))
        .with_context(|| format!("in {cov_file}"))?;
    opera_ms_core::note!(" *** {} repeats filtered", filter.nb_repeats());
    fs::create_dir_all(dir_res).with_context(|| format!("cannot create {dir_res}"))?;
    for i in 0..EDGE_SETS {
        let name = format!("filtered_pairedEdges_i{i}");
        let input = format!("{edge_file_dir}/{name}");
        let edges = fs::read_to_string(&input).with_context(|| format!("cannot read {input}"))?;
        let output = format!("{dir_res}/{name}");
        fs::write(&output, filter.filter_edges(&edges))
            .with_context(|| format!("cannot write {output}"))?;
    }
    let (repeat_info, clusters_cov) = filter.reports();
    for (name, text) in [
        ("repeat_contigs_info.txt", repeat_info),
        ("clusters_cov", clusters_cov),
    ] {
        let path = format!("{dir_res}/{name}");
        fs::write(&path, text).with_context(|| format!("cannot write {path}"))?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    const COVERAGE: &str = "1 500 80 340\n\
        a\t1000\tx\tx\t1\n340\n1\n\
        b\t2000\tx\tx\t1\n340\n1\n\
        c\t3000\tx\tx\t1\n3400\n1\n\
        d\t4000\tx\tx\t1\n68\n1\n";

    #[test]
    fn repeats_edges_and_reports() {
        // cluster 2: a (cov 1), b (1), c (10) -> median 1, c is a repeat;
        // cluster 10: d (0.2); e is in no cluster
        let clusters = "a\t2\t0\t1\nc\t2\t0\t1\nb\t2\t0\t1\nd\t10\t0\t1\n";
        let f = Filter::new(COVERAGE, clusters, 1.5).unwrap();
        assert_eq!(f.nb_repeats(), 1);
        let edges = "a\t+\tb\t+\t1\t1\t5\na\t+\tb\t+\t1\t1\t5\na\t+\tc\t-\t1\t1\t5\nb\t+\te\t+\t1\t1\t5\n  \nd\t+\ta\t-\t1\t1\t5\n";
        // first line always kept, a-c (repeat) and b-e (not clustered) and the
        // short line dropped
        assert_eq!(
            f.filter_edges(edges),
            "a\t+\tb\t+\t1\t1\t5\na\t+\tb\t+\t1\t1\t5\nd\t+\ta\t-\t1\t1\t5\n"
        );
        let (info, cov) = f.reports();
        assert_eq!(info, "2 1 6000 3\n\tc\t10\t3000\n10 0.2 4000 1\n");
        assert_eq!(cov, "2 1\n10 0.2\n");
    }

    #[test]
    fn perl_split() {
        assert_eq!(split_whitespace_perl(" a\tb  c\n"), ["", "a", "b", "c"]);
        assert_eq!(split_whitespace_perl("a b\n"), ["a", "b"]);
    }
}
