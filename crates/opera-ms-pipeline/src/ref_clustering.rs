//! Reference-based clustering (stage 5): merges sigma clusters whose
//! contigs are linked by edges that a common reference genome confirms.
//!
//! 1. skani: one FASTA per cluster against the run's reference genomes; the
//!    top 5 genomes of each cluster.
//! 2. Edges between clusters (long-read edges, and short-read edges with
//!    |distance| <= 500 and support >= 5), sorted by the difference of the
//!    cluster coverages.
//! 3. minimap2 of each linked cluster pair's contigs on their common
//!    genomes; an edge is good when the contig positions on a genome agree
//!    with it.
//! 4. Greedy merging of clusters along the good edges into super clusters
//!    (when they share a species), then the edges kept for OPERA-LG.
//!
//! Files written from hash maps are sorted, so that the output doesn't
//! depend on hash order.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;
use std::sync::Mutex;
use std::time::Instant;

use anyhow::{Context, Result, bail, ensure};
use opera_ms_core::bench::StepRecorder;
use opera_ms_core::parallel::run_parallel;
use opera_ms_core::runner::{Cmd, Output, Pipeline};

use crate::edge_filter::perl_number as perl_value;
use crate::files;
use crate::long_read::general_number;

const TOP_PREDICTIONS: usize = 5;
const EDGE_SETS: usize = 6;

/// Arguments of the reference-based clustering.
pub struct Args<'a> {
    pub inter_dir: &'a str,
    pub ref_dir: &'a str,
    pub mapping_dir: &'a str,
    pub sigma_dir: &'a str,
    pub contig_file: &'a str,
    pub nproc: usize,
    /// The run's reference genomes (`references::prepare`): files
    /// `<accession>.fna.gz` and `genomes_length.txt`.
    pub genome_dir: &'a str,
    /// Directories of `skani` and `minimap2` (with a trailing slash).
    pub skani_dir: &'a str,
    pub minimap2_dir: &'a str,
}

/// One line of sigma's `clusters` file.
#[derive(Debug, Default)]
struct ContigInfo {
    cluster: Option<String>,
    mean: Option<String>,
}

fn field(fields: &[&str], i: usize) -> Option<String> {
    fields.get(i).map(|s| (*s).to_string())
}

/// `get_cluster_contigs_info`: contigs, and the number of contigs per
/// cluster.
fn read_clusters(text: &str) -> (HashMap<String, ContigInfo>, HashMap<String, i64>) {
    let mut contigs = HashMap::new();
    let mut counts: HashMap<String, i64> = HashMap::new();
    for line in text.split_inclusive('\n') {
        let line = line.strip_suffix('\n').unwrap_or(line);
        let fields: Vec<&str> = line.split('\t').collect();
        let cluster = field(&fields, 1);
        *counts
            .entry(cluster.clone().unwrap_or_default())
            .or_insert(0) += 1;
        contigs.insert(
            fields.first().copied().unwrap_or("").to_string(),
            ContigInfo {
                cluster,
                mean: field(&fields, 3),
            },
        );
    }
    (contigs, counts)
}

/// The first field of `split(/\s+/, $text)`: empty when the text starts
/// with a space.
fn first_word(text: &str) -> &str {
    text.split(char::is_whitespace).next().unwrap_or("")
}

/// The contig loop of `run_mash_on_clusters`: a cluster's file is written
/// when its last contig has been read (the last contig of the file
/// always writes its cluster's file).
fn cluster_fastas(
    contig_file: impl BufRead,
    contigs: &HashMap<String, ContigInfo>,
    mut counts: HashMap<String, i64>,
    mut on_written: impl FnMut(&str, &str) -> Result<()>,
) -> Result<HashMap<String, String>> {
    let mut contig_seqs: HashMap<String, String> = HashMap::new();
    let mut cluster_seqs: HashMap<String, String> = HashMap::new();
    let mut name: Option<String> = None;
    let mut cluster: Option<String> = None;
    let mut seq = String::new();
    let cluster_of = |name: &str| contigs.get(name).and_then(|c| c.cluster.clone());
    for line in contig_file.split(b'\n') {
        let line = line.context("cannot read the contig file")?;
        let line = String::from_utf8_lossy(&line);
        if line.is_empty() {
            continue;
        }
        let header = line.strip_prefix('>').filter(|rest| !rest.is_empty());
        let Some(rest) = header else {
            seq.push_str(&line);
            seq.push('\n');
            continue;
        };
        if let (Some(n), Some(c)) = (&name, &cluster) {
            contig_seqs.entry(n.clone()).or_default().push_str(&seq);
            let cluster_seq = cluster_seqs.entry(c.clone()).or_default();
            cluster_seq.push_str(&seq);
            let count = counts.entry(c.clone()).or_insert(0);
            *count -= 1;
            if *count == 0 {
                let text = std::mem::take(cluster_seq);
                on_written(c, &text)?;
            }
        }
        seq = format!("{line}\n");
        let new_name = first_word(rest).to_string();
        cluster = cluster_of(&new_name);
        name = Some(new_name);
    }
    if let Some(n) = &name {
        contig_seqs.entry(n.clone()).or_default().push_str(&seq);
        if let Some(c) = &cluster {
            let cluster_seq = cluster_seqs.entry(c.clone()).or_default();
            cluster_seq.push_str(&seq);
            let text = std::mem::take(cluster_seq);
            on_written(c, &text)?;
        }
    }
    Ok(contig_seqs)
}

/// The sequences of the contigs of an edge between two clusters.
#[derive(Default)]
struct EdgeSeq {
    seq: String,
    contigs: HashSet<String>,
}

/// `generate_edges_between_clusters`: `edges_between_clusters` and the
/// contig sequences of each cluster pair.
fn edges_between_clusters(
    long_edges: &str,
    short_edges: &str,
    contigs: &HashMap<String, ContigInfo>,
    contig_seqs: &HashMap<String, String>,
) -> (String, HashMap<String, EdgeSeq>) {
    let mut cov_diff: BTreeMap<String, f64> = BTreeMap::new();
    let mut edge_seqs: HashMap<String, EdgeSeq> = HashMap::new();
    let cluster_of = |c: &str| contigs.get(c).and_then(|i| i.cluster.as_deref());
    let mean_of = |c: &str| {
        contigs
            .get(c)
            .and_then(|i| i.mean.as_deref())
            .map_or(0.0, perl_value)
    };
    let mut add = |line: &str, contig1: &str, contig2: &str, c1: &str, c2: &str| {
        let pair = if perl_value(c1) > perl_value(c2) {
            format!("{c2}-vs-{c1}")
        } else {
            format!("{c1}-vs-{c2}")
        };
        let edge = edge_seqs.entry(pair).or_default();
        for contig in [contig1, contig2] {
            if edge.contigs.insert(contig.to_string()) {
                edge.seq
                    .push_str(contig_seqs.get(contig).map_or("", String::as_str));
            }
        }
        cov_diff.insert(
            line.to_string(),
            (mean_of(contig1) - mean_of(contig2)).abs(),
        );
    };
    for line in lines(long_edges) {
        let fields: Vec<&str> = line.split('\t').collect();
        let contig1 = fields.first().copied().unwrap_or("");
        let contig2 = fields.get(2).copied().unwrap_or("");
        // both contigs clustered
        let (Some(c1), Some(c2)) = (cluster_of(contig1), cluster_of(contig2)) else {
            continue;
        };
        if c1 == c2 {
            continue;
        }
        add(line, contig1, contig2, c1, c2);
    }
    for line in lines(short_edges) {
        let fields: Vec<&str> = line.split('\t').collect();
        let contig1 = fields.first().copied().unwrap_or("");
        let contig2 = fields.get(2).copied().unwrap_or("");
        let (Some(c1), Some(c2)) = (cluster_of(contig1), cluster_of(contig2)) else {
            continue;
        };
        if c1 == c2 {
            continue;
        }
        let overlap = perl_value(fields.get(4).copied().unwrap_or(""));
        let support = perl_value(fields.get(6).copied().unwrap_or(""));
        if overlap.abs() > 500.0 || support < 5.0 {
            continue;
        }
        add(line, contig1, contig2, c1, c2);
    }
    // sorted lines, then a stable sort on the coverage difference
    let mut edges: Vec<(&String, &f64)> = cov_diff.iter().collect();
    edges.sort_by(|(_, a), (_, b)| a.total_cmp(b));
    let text = edges.iter().map(|(line, _)| format!("{line}\n")).collect();
    (text, edge_seqs)
}

/// The lines of a file without their `\n`.
fn lines(text: &str) -> impl Iterator<Item = &str> {
    text.split_inclusive('\n')
        .map(|l| l.strip_suffix('\n').unwrap_or(l))
}

/// `sort -k3,3 -g` (C locale): the third field as a general number,
/// non-numbers first, then the whole line bytewise.
fn sort_by_distance(text: &str) -> Vec<&str> {
    let mut keyed: Vec<(Option<f64>, &str)> = lines(text)
        .map(|l| {
            let key = l.split_whitespace().nth(2).and_then(general_number);
            (key, l)
        })
        .collect();
    keyed.sort_by(|(ka, la), (kb, lb)| {
        let by_key = match (ka, kb) {
            (None, None) => std::cmp::Ordering::Equal,
            (None, Some(_)) => std::cmp::Ordering::Less,
            (Some(_), None) => std::cmp::Ordering::Greater,
            (Some(x), Some(y)) => x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal),
        };
        by_key.then_with(|| la.as_bytes().cmp(lb.as_bytes()))
    });
    keyed.into_iter().map(|(_, l)| l).collect()
}

/// `generate_clusters_species`: the top 5 genomes of each cluster, by
/// Mash distance.
fn cluster_species(mash_dist: &str) -> HashMap<String, Vec<String>> {
    let mut species: HashMap<String, Vec<String>> = HashMap::new();
    for line in sort_by_distance(mash_dist) {
        let fields: Vec<&str> = line.split('\t').collect();
        let query = fields.get(1).copied().unwrap_or("");
        let file = query.split('/').rfind(|s| !s.is_empty()).unwrap_or("");
        // only the extension
        let cluster = file.strip_suffix(".fa").unwrap_or(file).to_string();
        let list = species.entry(cluster).or_default();
        if list.len() < TOP_PREDICTIONS {
            list.push(fields.first().copied().unwrap_or("").to_string());
        }
    }
    species
}

/// The species label of a genome file: the name up to the first `__`
/// (`Species__Accession_genomic.fna.gz`), or the accession for
/// `<accession>.fna.gz`.
pub(crate) fn species_name(strain: &str) -> &str {
    let name = strain.split('/').rfind(|s| !s.is_empty()).unwrap_or("");
    let name = name.split("__").next().unwrap_or("");
    name.strip_suffix(".fna.gz").unwrap_or(name)
}

/// One nucmer run: a cluster pair's contigs on one common genome.
struct NucmerJob {
    pair: String,
    index: usize,
    genome: String,
}

/// What `run_nucmer` prepares: the pair FASTA files, the nucmer runs, and
/// the number of common genomes of each pair (`NUMB_SPECIES`).
#[derive(Default)]
struct NucmerPlan {
    fastas: Vec<(String, String)>,
    jobs: Vec<NucmerJob>,
    nb_genomes: HashMap<String, usize>,
}

fn nucmer_plan(
    edges: &str,
    contigs: &HashMap<String, ContigInfo>,
    species: &HashMap<String, Vec<String>>,
    edge_seqs: &HashMap<String, EdgeSeq>,
) -> NucmerPlan {
    let mut plan = NucmerPlan::default();
    let mut seen: HashSet<String> = HashSet::new();
    let cluster_of = |c: &str| {
        contigs
            .get(c)
            .and_then(|i| i.cluster.clone())
            .unwrap_or_default()
    };
    for line in lines(edges) {
        let fields: Vec<&str> = line.split('\t').collect();
        let c1 = cluster_of(fields.first().copied().unwrap_or(""));
        let c2 = cluster_of(fields.get(2).copied().unwrap_or(""));
        let pair = format!("{c1}-vs-{c2}");
        if !seen.contains(&pair) && !seen.contains(&format!("{c2}-vs-{c1}")) {
            // compare_cluster
            if let (Some(g1), Some(g2)) = (species.get(&c1), species.get(&c2)) {
                let seq = edge_seqs
                    .get(&pair)
                    .or_else(|| edge_seqs.get(&format!("{c2}-vs-{c1}")))
                    .map_or("", |e| e.seq.as_str());
                plan.fastas.push((pair.clone(), seq.to_string()));
                let mut index = 0;
                for genome in g1 {
                    for other in g2 {
                        if other == genome {
                            plan.jobs.push(NucmerJob {
                                pair: pair.clone(),
                                index,
                                genome: genome.clone(),
                            });
                            index += 1;
                        }
                    }
                }
                if index > 0 {
                    plan.nb_genomes.insert(pair.clone(), index);
                }
            }
        }
        seen.insert(pair);
    }
    plan
}

/// The useful columns of a `show-coords -lrcT` line.
struct Hit {
    start: f64,
    end: f64,
    ref_length: f64,
    percent: f64,
    contig: String,
    /// The genome sequence the contig is on (`TAGS` Q).
    target: String,
}

/// The hits of one genome, read once and shared.
type Coords = std::rc::Rc<Vec<Hit>>;

/// A `show-coords -lrcT` file: its hits (after the 4 header lines).
fn parse_coords(text: &str) -> Vec<Hit> {
    text.split_inclusive('\n')
        .skip(4)
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            let get = |i: usize| fields.get(i).copied().unwrap_or("");
            let (a, b) = (perl_value(get(2)), perl_value(get(3)));
            Hit {
                start: a.min(b),
                end: a.max(b),
                ref_length: perl_value(get(8)),
                percent: perl_value(get(9)),
                contig: get(11).to_string(),
                target: get(12).trim_end().to_string(),
            }
        })
        .collect()
}

/// `check_reference_position`: whether the distance of an edge agrees with
/// the contig positions on one of the pair's genomes. `coords(pair, i)`
/// gives the hits of the pair's i-th genome.
fn check_reference_position(
    contig1: &str,
    contig2: &str,
    length: f64,
    stdev: f64,
    pair: &str,
    nb_genomes: usize,
    coords: &mut impl FnMut(&str, usize) -> Result<Coords>,
) -> Result<bool> {
    // the hits of one genome at a time; two positions are compared only on
    // the same genome sequence, with that sequence's length
    let margin = 3.0 * stdev.abs();
    for i in 0..nb_genomes {
        let hits = coords(pair, i)?;
        let hits1: Vec<&Hit> = hits.iter().filter(|h| h.contig == contig1).collect();
        let hits2: Vec<&Hit> = hits.iter().filter(|h| h.contig == contig2).collect();
        if hits1.is_empty() || hits2.is_empty() {
            continue;
        }
        // repeat detection on the fraction of the contig mapped
        let percent1: f64 = hits1.iter().map(|h| h.percent).sum();
        let percent2: f64 = hits2.iter().map(|h| h.percent).sum();
        if percent1 > 150.0 || percent2 > 150.0 {
            continue;
        }
        for h1 in &hits1 {
            for h2 in hits2.iter().filter(|h| h.target == h1.target) {
                let ref_length = h1.ref_length;
                let (p1, p2) = ((h1.start, h1.end), (h2.start, h2.end));
                let (min, max, min_circ, max_circ) = if p2.1 > p1.1 {
                    (
                        p2.0 - p1.1 - margin,
                        p2.1 - p1.0 + margin,
                        p1.0 - p2.1 - margin + ref_length,
                        p1.1 - p2.0 + margin + ref_length,
                    )
                } else {
                    (
                        p1.0 - p2.1 - margin,
                        p1.1 - p2.0 + margin,
                        p2.0 - p1.1 - margin + ref_length,
                        p2.1 - p1.0 + margin + ref_length,
                    )
                };
                if (max > length && length > min) || (max_circ > length && min_circ < length) {
                    return Ok(true);
                }
            }
        }
    }
    Ok(false)
}

/// The super clusters made by `merge_clusters`.
#[derive(Default)]
struct Merge {
    super_of: HashMap<String, usize>,
    /// `None`: merged into another one.
    supers: Vec<Option<Vec<String>>>,
    species: HashMap<usize, Vec<String>>,
    strains: HashMap<usize, HashMap<String, u64>>,
}

impl Merge {
    /// `common_species`.
    fn common_species(&self, c1: &str, c2: &str, clusters: &HashMap<String, Vec<String>>) -> bool {
        let (Some(list1), Some(list2)) = (clusters.get(c1), clusters.get(c2)) else {
            return false;
        };
        let names = |c: &str, list: &Vec<String>| -> Vec<String> {
            match self.super_of.get(c) {
                Some(id) => self.species.get(id).cloned().unwrap_or_default(),
                None => list.iter().map(|s| species_name(s).to_string()).collect(),
            }
        };
        let (s1, s2) = (names(c1, list1), names(c2, list2));
        s1.iter().any(|a| s2.contains(a))
    }

    fn merge(
        &mut self,
        mut c1: String,
        mut c2: String,
        clusters: &HashMap<String, Vec<String>>,
        log: &mut impl Write,
    ) -> io::Result<()> {
        if !self.common_species(&c1, &c2, clusters) {
            return Ok(());
        }
        writeln!(log, " *** *** *** Clusters to merge {c1} {c2}")?;
        let empty = Vec::new();
        match (
            self.super_of.get(&c1).copied(),
            self.super_of.get(&c2).copied(),
        ) {
            (None, None) => {
                let id = self.supers.len();
                let (list1, list2) = (
                    clusters.get(&c1).unwrap_or(&empty),
                    clusters.get(&c2).unwrap_or(&empty),
                );
                let mut intersection = BTreeSet::new();
                for s1 in list1 {
                    for s2 in list2 {
                        if species_name(s1) == species_name(s2) {
                            intersection.insert(species_name(s1).to_string());
                        }
                    }
                }
                let mut strains: HashMap<String, u64> = HashMap::new();
                for s in list1 {
                    strains.insert(s.clone(), 1);
                }
                for s in list2 {
                    *strains.entry(s.clone()).or_insert(0) += 1;
                }
                self.super_of.insert(c1.clone(), id);
                self.super_of.insert(c2.clone(), id);
                self.supers.push(Some(vec![c1, c2]));
                self.strains.insert(id, strains);
                self.species.insert(id, intersection.into_iter().collect());
            }
            (Some(_), None) | (None, Some(_)) => {
                if self.super_of.contains_key(&c2) {
                    std::mem::swap(&mut c1, &mut c2);
                }
                let Some(&id) = self.super_of.get(&c1) else {
                    return Ok(());
                };
                self.super_of.insert(c2.clone(), id);
                if let Some(Some(members)) = self.supers.get_mut(id) {
                    members.push(c2.clone());
                }
                let list2 = clusters.get(&c2).unwrap_or(&empty);
                let mut intersection = BTreeSet::new();
                for s in self.species.get(&id).unwrap_or(&empty) {
                    for s2 in list2 {
                        if species_name(s2) == s {
                            intersection.insert(s.clone());
                        }
                    }
                }
                // tested against a copy made before the loop
                let before: HashSet<String> = self
                    .strains
                    .get(&id)
                    .map(|s| s.keys().cloned().collect())
                    .unwrap_or_default();
                let strains = self.strains.entry(id).or_default();
                for s in list2 {
                    if before.contains(s) {
                        *strains.entry(s.clone()).or_insert(0) += 1;
                    } else {
                        strains.insert(s.clone(), 1);
                    }
                }
                self.species.insert(id, intersection.into_iter().collect());
            }
            (Some(id1), Some(id2)) => {
                if id1 == id2 {
                    return Ok(());
                }
                writeln!(log, " *** Merge super cluster !!!")?;
                let moved = self
                    .supers
                    .get_mut(id2)
                    .and_then(Option::take)
                    .unwrap_or_default();
                for c in &moved {
                    self.super_of.insert(c.clone(), id1);
                }
                if let Some(Some(members)) = self.supers.get_mut(id1) {
                    members.extend(moved);
                }
                let mut intersection = BTreeSet::new();
                let species2 = self.species.remove(&id2).unwrap_or_default();
                for s in self.species.get(&id1).unwrap_or(&empty) {
                    if species2.contains(s) {
                        intersection.insert(s.clone());
                    }
                }
                let strains2 = self.strains.remove(&id2).unwrap_or_default();
                let strains1 = self.strains.entry(id1).or_default();
                for (s, n) in strains2 {
                    *strains1.entry(s).or_insert(0) += n;
                }
                self.species.insert(id1, intersection.into_iter().collect());
            }
        }
        Ok(())
    }
}

/// The edge set (`filtered_pairedEdges_i<n>`) of a distance, if any.
fn edge_set(distance: f64) -> Option<usize> {
    [300.0, 1000.0, 2000.0, 5000.0, 15000.0, 40000.0]
        .iter()
        .position(|&limit| distance < limit)
}

/// `recover_edges`: the good edges inside a super cluster, and the
/// above-threshold edges inside a cluster with support > 1.
fn recover_edges(
    good_edges: &str,
    above_threshold: &str,
    contigs: &HashMap<String, ContigInfo>,
    merge: &Merge,
    too_long: &mut usize,
) -> Result<[String; EDGE_SETS]> {
    let mut sets: [String; EDGE_SETS] = Default::default();
    let cluster_of = |c: &str| contigs.get(c).and_then(|i| i.cluster.as_deref());
    for line in lines(good_edges) {
        let fields: Vec<&str> = line.split('\t').collect();
        let get = |i: usize| fields.get(i).copied().unwrap_or("");
        let super1 = cluster_of(get(0)).and_then(|c| merge.super_of.get(c));
        let super2 = cluster_of(get(2)).and_then(|c| merge.super_of.get(c));
        let (Some(s1), Some(s2)) = (super1, super2) else {
            continue;
        };
        if s1 == s2
            && let Some(set) = edge_set(perl_value(get(4)))
        {
            sets[set].push_str(line);
            sets[set].push('\n');
        }
    }
    for line in lines(above_threshold) {
        let fields: Vec<&str> = line.split('\t').collect();
        let get = |i: usize| fields.get(i).copied().unwrap_or("");
        // two unclustered contigs aren't "the same cluster"; an edge of 40 kb
        // or more, outside the six edge sets, is skipped
        let same = match (cluster_of(get(0)), cluster_of(get(2))) {
            (Some(c1), Some(c2)) => c1 == c2,
            _ => false,
        };
        if same && perl_value(get(6)) > 1.0 {
            match edge_set(perl_value(get(4))) {
                Some(set) => {
                    sets[set].push_str(line);
                    sets[set].push('\n');
                }
                None => *too_long += 1,
            }
        }
    }
    Ok(sets)
}

/// `export_info`: `clusters_seq_similarity`, `super_cluster.dat` and
/// `cluster_species.dat`.
fn export_info(
    contigs: &HashMap<String, ContigInfo>,
    species: &HashMap<String, Vec<String>>,
    merge: &Merge,
) -> (String, String, String) {
    let mut clusters = String::new();
    let mut names: Vec<&String> = contigs.keys().collect();
    names.sort();
    for name in names {
        let info = &contigs[name.as_str()];
        let Some(cluster) = &info.cluster else {
            continue;
        };
        match merge.super_of.get(cluster) {
            Some(id) => clusters.push_str(&format!("{name}\t{id}S\t001100\t0.0088888\n")),
            None => clusters.push_str(&format!(
                "{name}\t{cluster}\t001100\t{}\n",
                info.mean.as_deref().unwrap_or("")
            )),
        }
    }
    let mut supers = String::new();
    for (i, members) in merge.supers.iter().enumerate() {
        if let Some(members) = members {
            supers.push_str(&format!("{i}S\t{}\n", members.join("\t")));
        }
    }
    let mut out = String::new();
    let mut ids: Vec<(String, usize)> = merge
        .strains
        .keys()
        .map(|id| (id.to_string(), *id))
        .collect();
    ids.sort();
    for (text, id) in ids {
        out.push_str(&format!(">{text}S\n"));
        let mut strains: Vec<(&String, &u64)> = merge.strains[&id].iter().collect();
        strains.sort_by(|(sa, na), (sb, nb)| nb.cmp(na).then_with(|| sa.cmp(sb)));
        for (strain, _) in strains {
            out.push_str(&format!("{strain}\n"));
        }
    }
    let mut clusters_sorted: Vec<(&String, &Vec<String>)> = species.iter().collect();
    clusters_sorted.sort();
    for (id, list) in clusters_sorted {
        out.push_str(&format!(">{id}\n"));
        for s in list {
            out.push_str(&format!("{s}\n"));
        }
    }
    (clusters, supers, out)
}

fn read(path: &str) -> Result<String> {
    fs::read_to_string(path).with_context(|| format!("cannot read {path}"))
}

fn write(path: &str, text: &str) -> Result<()> {
    fs::write(path, text).with_context(|| format!("cannot write {path}"))
}

fn remove_dir(path: &str) -> Result<()> {
    if Path::new(path).exists() {
        fs::remove_dir_all(path).with_context(|| format!("cannot remove {path}"))?;
    }
    Ok(())
}

fn time_log(inter_dir: &str, step: &str, start: Instant) -> Result<()> {
    let path = format!("{inter_dir}/time.log");
    let mut file = fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(&path)
        .with_context(|| format!("cannot open {path}"))?;
    file.write_all(files::time_log_line(step, start.elapsed().as_secs()).as_bytes())
        .with_context(|| format!("cannot write {path}"))
}

/// The single `lib_1_bundles/clusters_*` file (short-read edges).
fn short_edges_file(inter_dir: &str) -> Result<String> {
    let dir = format!("{inter_dir}/lib_1_bundles");
    let mut found = Vec::new();
    for entry in fs::read_dir(&dir).with_context(|| format!("cannot read {dir}"))? {
        let name = entry
            .with_context(|| format!("cannot read {dir}"))?
            .file_name();
        let name = name.to_string_lossy();
        if name.starts_with("clusters_") {
            found.push(format!("{dir}/{name}"));
        }
    }
    match found.as_slice() {
        [one] => Ok(one.clone()),
        [] => bail!("no {dir}/clusters_* file"),
        _ => bail!("several {dir}/clusters_* files (multiple samples are not supported)"),
    }
}

/// The species of the clusters, with skani: one FASTA per cluster, `skani dist` against the run's
/// reference genomes (the ones sylph found), at most `TOP_PREDICTIONS`
/// per cluster, written in the `mash_dist.dat` layout (distance =
/// 1 - ANI/100, genomes named `<accession>.fna.gz`). Returns the edge
/// sequences of the clusters.
fn run_skani(
    rec: &StepRecorder<'_>,
    a: &Args<'_>,
    contigs: &HashMap<String, ContigInfo>,
    counts: HashMap<String, i64>,
) -> Result<HashMap<String, String>> {
    let skani_dir = format!("{}/SKANI/", a.ref_dir);
    remove_dir(&skani_dir)?;
    let inter_fa_dir = format!("{skani_dir}/INTER_FILES/");
    fs::create_dir_all(&inter_fa_dir).with_context(|| format!("cannot create {inter_fa_dir}"))?;
    let mut queries: Vec<String> = Vec::new();
    let reader = BufReader::new(
        File::open(a.contig_file)
            .with_context(|| format!("Could not open contigs file {}", a.contig_file))?,
    );
    let contig_seqs = cluster_fastas(reader, contigs, counts, |cluster, text| {
        let path = format!("{inter_fa_dir}/{cluster}.fa");
        write(&path, text)?;
        if !queries.contains(&path) {
            queries.push(path);
        }
        Ok(())
    })?;
    ensure!(
        !queries.is_empty(),
        "no cluster sequence for skani (empty {})",
        a.contig_file
    );
    queries.sort();
    let genomes = crate::references::genome_files(&format!("{}/..", a.genome_dir))?;
    let dist = crate::references::skani_as_mash(
        rec,
        a.skani_dir,
        &genomes,
        &queries,
        TOP_PREDICTIONS,
        a.nproc,
        &format!("{skani_dir}/clusters"),
    )?;
    // genomes by file name
    let dist: String = dist
        .lines()
        .map(|l| {
            let (genome, rest) = l.split_once('\t').unwrap_or((l, ""));
            format!("{}\t{rest}\n", genome.rsplit('/').next().unwrap_or(genome))
        })
        .collect();
    write(&format!("{skani_dir}/mash_dist.dat"), &dist)?;
    remove_dir(&inter_fa_dir)?;
    Ok(contig_seqs)
}

/// The contig-on-genome alignments of the linked cluster pairs (was
/// The contigs of each cluster pair on their common genomes: minimap2
/// (`-cx asm5 -N 50`, secondary alignments kept for repeats), one run per
/// genome with every contig of the pairs that need
/// it (an alignment doesn't depend on the other contigs), then each
/// pair's file `<pair>_<i>.txt` from its own contigs' hits, in the
/// `show-coords -lrcT` columns that `parse_coords` reads.
fn run_minimap2(
    rec: &StepRecorder<'_>,
    a: &Args<'_>,
    nucmer_dir: &str,
    plan: &NucmerPlan,
) -> Result<()> {
    remove_dir(nucmer_dir)?;
    fs::create_dir_all(nucmer_dir).with_context(|| format!("cannot create {nucmer_dir}"))?;
    // the contigs of each pair (FASTA text), and of each genome
    let mut pair_contigs: HashMap<&str, Vec<String>> = HashMap::new();
    let mut records: HashMap<String, String> = HashMap::new();
    for (pair, seq) in &plan.fastas {
        write(&format!("{nucmer_dir}/{pair}.fa"), seq)?;
        let mut name = String::new();
        for line in seq.lines() {
            if let Some(h) = line.strip_prefix('>') {
                name = h.split_whitespace().next().unwrap_or("").to_string();
                pair_contigs.entry(pair).or_default().push(name.clone());
                records.entry(name.clone()).or_default();
            } else if let Some(r) = records.get_mut(&name) {
                if r.len() < line.len() || r.is_empty() {
                    // one record per contig (the same sequence in every pair)
                    *r = line.to_string();
                }
            }
        }
    }
    let mut by_genome: BTreeMap<&str, BTreeSet<&str>> = BTreeMap::new();
    for job in &plan.jobs {
        let set = by_genome.entry(job.genome.as_str()).or_default();
        for c in pair_contigs.get(job.pair.as_str()).into_iter().flatten() {
            set.insert(c.as_str());
        }
    }
    let hits: Mutex<HashMap<String, Vec<String>>> = Mutex::new(HashMap::new());
    let mut jobs = Vec::new();
    for (g, (genome, contig_set)) in by_genome.iter().enumerate() {
        let query = format!("{nucmer_dir}/genome_{g}_contigs.fa");
        let mut text = String::new();
        for c in contig_set {
            text.push_str(&format!(
                ">{c}\n{}\n",
                records.get(*c).map_or("", String::as_str)
            ));
        }
        write(&query, &text)?;
        let genome_file = genome_path(a.genome_dir, genome);
        let paf = format!("{nucmer_dir}/genome_{g}.paf");
        let pipeline = Pipeline::new(
            Cmd::new(format!("{}minimap2", a.minimap2_dir))
                .args(["-cx", "asm5", "-N", "50", "-t", "1"])
                .args([&genome_file, &query]),
        )
        .stdout(Output::File(paf.clone().into()))
        .stderr(Output::File(format!("{nucmer_dir}/genome_{g}.err").into()));
        let (hits, genome) = (&hits, genome.to_string());
        jobs.push(move || -> Result<()> {
            let report = pipeline.run().with_context(|| {
                format!("Error during minimap2 of the cluster contigs on {genome}")
            })?;
            rec.record(&report);
            let text = read(&paf)?;
            let mut rows: Vec<(String, String)> = Vec::new();
            for line in text.lines() {
                let f: Vec<&str> = line.split('\t').collect();
                if f.len() < 11 {
                    continue;
                }
                let num = |i: usize| f[i].parse::<u64>().unwrap_or(0);
                let (qlen, qs, qe, tlen, ts, te, m, b) = (
                    num(1),
                    num(2),
                    num(3),
                    num(6),
                    num(7),
                    num(8),
                    num(9),
                    num(10),
                );
                if qlen == 0 || tlen == 0 || b == 0 {
                    continue;
                }
                let (s2, e2) = if f[4] == "-" {
                    (te, ts + 1)
                } else {
                    (ts + 1, te)
                };
                // [S1] [E1] [S2] [E2] [LEN 1] [LEN 2] [% IDY] [LEN R] [LEN Q]
                // [COV R] [COV Q] [TAGS]: R = the contig, Q = the genome
                let row = format!(
                    "{}\t{qe}\t{s2}\t{e2}\t{}\t{}\t{:.2}\t{qlen}\t{tlen}\t{:.2}\t{:.2}\t{}\t{}\n",
                    qs + 1,
                    qe - qs,
                    te - ts,
                    100.0 * m as f64 / b as f64,
                    100.0 * (qe - qs) as f64 / qlen as f64,
                    100.0 * (te - ts) as f64 / tlen as f64,
                    f[0],
                    f[5]
                );
                rows.push((f[0].to_string(), row));
            }
            if let Ok(mut h) = hits.lock() {
                h.insert(
                    genome,
                    rows.into_iter()
                        .map(|(c, r)| format!("{c}\u{1}{r}"))
                        .collect(),
                );
            }
            Ok(())
        });
    }
    run_parallel(jobs, a.nproc.max(1)).map_err(|e| anyhow::anyhow!("{e}"))?;
    let hits = hits
        .into_inner()
        .map_err(|_| anyhow::anyhow!("lock poisoned"))?;
    let header = "minimap2\nNUCMER\n\n[S1]\t[E1]\t[S2]\t[E2]\t[LEN 1]\t[LEN 2]\t[% IDY]\t[LEN R]\t[LEN Q]\t[COV R]\t[COV Q]\t[TAGS]\n";
    for job in &plan.jobs {
        let wanted: HashSet<&str> = pair_contigs
            .get(job.pair.as_str())
            .map(|v| v.iter().map(String::as_str).collect())
            .unwrap_or_default();
        let mut text = header.to_string();
        for entry in hits.get(&job.genome).into_iter().flatten() {
            let (contig, row) = entry.split_once('\u{1}').unwrap_or(("", entry));
            if wanted.contains(contig) {
                text.push_str(row);
            }
        }
        write(
            &format!("{nucmer_dir}/{}_{}.txt", job.pair, job.index),
            &text,
        )?;
    }
    for g in 0..by_genome.len() {
        for name in [format!("genome_{g}_contigs.fa"), format!("genome_{g}.paf")] {
            let path = format!("{nucmer_dir}/{name}");
            fs::remove_file(&path).with_context(|| format!("cannot remove {path}"))?;
        }
    }
    Ok(())
}

/// `get_correct_edges`: the edges between clusters confirmed by a genome.
fn correct_edges(
    edges: &str,
    contigs: &HashMap<String, ContigInfo>,
    plan: &NucmerPlan,
    nucmer_dir: &str,
) -> Result<String> {
    let mut cache: HashMap<(String, usize), Coords> = HashMap::new();
    let mut coords = |pair: &str, i: usize| -> Result<Coords> {
        let key = (pair.to_string(), i);
        if let Some(hits) = cache.get(&key) {
            return Ok(hits.clone());
        }
        let path = format!("{nucmer_dir}/{pair}_{i}.txt");
        let hits = std::rc::Rc::new(parse_coords(
            &read(&path).context("Parsing problem during read rescue using NUCMER mapping")?,
        ));
        cache.insert(key, hits.clone());
        Ok(hits)
    };
    let mut good = String::new();
    for line in lines(edges) {
        let fields: Vec<&str> = line.split('\t').collect();
        let get = |i: usize| fields.get(i).copied().unwrap_or("");
        let (contig1, contig2) = (get(0), get(2));
        let cluster_of = |c: &str| {
            contigs
                .get(c)
                .map(|i| i.cluster.clone().unwrap_or_default())
        };
        let (Some(c1), Some(c2)) = (cluster_of(contig1), cluster_of(contig2)) else {
            continue;
        };
        let pair = [format!("{c1}-vs-{c2}"), format!("{c2}-vs-{c1}")]
            .into_iter()
            .find_map(|p| plan.nb_genomes.get(&p).map(|n| (p, *n)));
        let Some((pair, nb)) = pair else {
            continue;
        };
        let (overlap, stdev) = (perl_value(get(4)), perl_value(get(5)));
        if check_reference_position(contig1, contig2, overlap, stdev, &pair, nb, &mut coords)? {
            good.push_str(line);
            good.push('\n');
        }
    }
    Ok(good)
}

/// The reference-based clustering.
pub fn sequence_similarity_clustering(rec: &StepRecorder<'_>, a: &Args<'_>) -> Result<()> {
    let ref_dir = a.ref_dir;
    fs::create_dir_all(ref_dir).with_context(|| format!("cannot create {ref_dir}"))?;
    let log_path = format!("{ref_dir}/sequence_similarity_clustering.err");
    let log =
        Mutex::new(File::create(&log_path).with_context(|| format!("cannot write {log_path}"))?);
    let clusters_file = format!("{}/clusters", a.sigma_dir);
    let (contigs, counts) = read_clusters(&read(&clusters_file)?);

    let start = Instant::now();
    let contig_seqs = rec.substep("skani", |sub| run_skani(sub, a, &contigs, counts))?;
    let long_edges = read(&format!("{}/pairedEdges", a.mapping_dir))?;
    let short_edges = read(&short_edges_file(a.inter_dir)?)?;
    let (edges, edge_seqs) =
        edges_between_clusters(&long_edges, &short_edges, &contigs, &contig_seqs);
    drop(contig_seqs);
    write(&format!("{ref_dir}/edges_between_clusters"), &edges)?;
    time_log(a.inter_dir, "r_mash", start)?;
    let species = cluster_species(&read(&format!("{ref_dir}/SKANI/mash_dist.dat"))?);

    let start = Instant::now();
    let nucmer_dir = format!("{ref_dir}/NUCMER_OUT/");
    let plan = nucmer_plan(&edges, &contigs, &species, &edge_seqs);
    drop(edge_seqs);
    rec.substep("minimap2", |sub| run_minimap2(sub, a, &nucmer_dir, &plan))?;
    let good = correct_edges(&edges, &contigs, &plan, &nucmer_dir)?;
    write(&format!("{ref_dir}/edges_between_clusters_good"), &good)?;
    time_log(a.inter_dir, "r_mummer", start)?;

    let start = Instant::now();
    let mut merge = Merge::default();
    let mut log = log
        .into_inner()
        .map_err(|_| anyhow::anyhow!("log lock poisoned"))?;
    for line in lines(&good) {
        let fields: Vec<&str> = line.split('\t').collect();
        let cluster_of = |i: usize| {
            fields
                .get(i)
                .and_then(|c| contigs.get(*c))
                .and_then(|c| c.cluster.clone())
                .unwrap_or_default()
        };
        merge
            .merge(cluster_of(0), cluster_of(2), &species, &mut log)
            .with_context(|| format!("cannot write {log_path}"))?;
    }
    let above = read(&format!("{}/pairedEdges_above_thresh", a.mapping_dir))?;
    let mut too_long = 0;
    let sets = recover_edges(&good, &above, &contigs, &merge, &mut too_long)?;
    if too_long > 0 {
        opera_ms_core::note!(
            " *** {too_long} above-threshold edge(s) of 40 kb or more skipped (no edge set)"
        );
    }
    for (i, text) in sets.iter().enumerate() {
        write(&format!("{ref_dir}/filtered_pairedEdges_i{i}"), text)?;
    }
    let (clusters, supers, cluster_species) = export_info(&contigs, &species, &merge);
    write(&format!("{ref_dir}/clusters_seq_similarity"), &clusters)?;
    write(&format!("{ref_dir}/super_cluster.dat"), &supers)?;
    write(&format!("{ref_dir}/cluster_species.dat"), &cluster_species)?;
    time_log(a.inter_dir, "r_cluster_merging", start)
}

/// The file of a genome: `<db>/<name>`, without a leading `OPERA-MS-DB/`
/// (genome names of OPERA-MS 0.9 databases).
pub fn genome_path(db: &str, name: &str) -> String {
    let inside = name.strip_prefix("OPERA-MS-DB/").unwrap_or(name);
    format!("{}/{inside}", db.trim_end_matches('/'))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn contigs(spec: &[(&str, &str, &str)]) -> HashMap<String, ContigInfo> {
        let text: String = spec
            .iter()
            .map(|(c, cl, mean)| format!("{c}\t{cl}\t0\t{mean}\n"))
            .collect();
        read_clusters(&text).0
    }

    #[test]
    fn cluster_files_written_when_complete() {
        let info = contigs(&[("a", "1", "10"), ("b", "2", "20"), ("c", "1", "10")]);
        let (_, counts) = read_clusters("a\t1\t0\t10\nb\t2\t0\t20\nc\t1\t0\t10\n");
        let fasta = ">a x\nAC\nGT\n\n>u\nTT\n>b\nCC\n>c\nGG\n";
        let mut written = Vec::new();
        let seqs = cluster_fastas(fasta.as_bytes(), &info, counts, |c, t| {
            written.push((c.to_string(), t.to_string()));
            Ok(())
        })
        .unwrap();
        // cluster 2 (b) completes first; the last contig writes cluster 1
        assert_eq!(
            written,
            [
                ("2".to_string(), ">b\nCC\n".to_string()),
                ("1".to_string(), ">a x\nAC\nGT\n>c\nGG\n".to_string())
            ]
        );
        assert_eq!(seqs["a"], ">a x\nAC\nGT\n");
        assert!(!seqs.contains_key("u"));
    }

    #[test]
    fn edges_sorted_by_coverage_difference() {
        let info = contigs(&[("a", "1", "10"), ("b", "2", "20"), ("c", "3", "12")]);
        let seqs: HashMap<String, String> = [("a", ">a\n"), ("b", ">b\n"), ("c", ">c\n")]
            .iter()
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect();
        let long = "b\t+\ta\t-\t100\t10\t3\na\t+\tc\t+\t50\t10\t2\na\t+\tb\t+\t5\t1\t1\nx\t+\ta\t+\t1\t1\t1\nx\t+\ty\t+\t1\t1\t1\n";
        // short: |overlap| > 500 or support < 5 dropped
        let short = "a\t+\tc\t-\t-20\t5\t5\na\t+\tc\t-\t600\t5\t9\na\t+\tb\t-\t0\t5\t4\n";
        let (text, edge_seqs) = edges_between_clusters(long, short, &info, &seqs);
        assert_eq!(
            text,
            // x is unclustered: its edges are dropped
            "a\t+\tc\t+\t50\t10\t2\na\t+\tc\t-\t-20\t5\t5\na\t+\tb\t+\t5\t1\t1\nb\t+\ta\t-\t100\t10\t3\n"
        );
        // pair id ordered by cluster number; contigs in order of appearance
        assert_eq!(edge_seqs["1-vs-2"].seq, ">b\n>a\n");
        // no pair with an unclustered contig
        assert!(!edge_seqs.contains_key("-vs-1"));
    }

    #[test]
    fn top_species_by_distance() {
        let dist = "g/s1__x.gz\t/d//1.fa\t0.2\t0\t1/1000\n\
                    g/s2__y.gz\t/d//1.fa\t0.1\t0\t1/1000\n\
                    g/s3__z.gz\t/d//10.fa\t1e-2\t0\t1/1000\n";
        let species = cluster_species(dist);
        assert_eq!(species["1"], ["g/s2__y.gz", "g/s1__x.gz"]);
        assert_eq!(species["10"], ["g/s3__z.gz"]);
        assert_eq!(species_name("g/s1__x.gz"), "s1");
    }

    #[test]
    fn reference_positions() {
        let header = "a\nNUCMER\n\n[S1]\n";
        // contigs 2 kb apart on the genome (spanning 4 kb), stdev 10
        let coords = format!(
            "{header}1\t1000\t1000\t2000\t0\t0\t99\t1000\t100000\t100.0\t1\tc1\tg\n\
             1\t1000\t5000\t4000\t0\t0\t99\t1000\t100000\t100.0\t1\tc2\tg\n"
        );
        let hits = std::rc::Rc::new(parse_coords(&coords));
        let mut get = |_: &str, _: usize| Ok(hits.clone());
        assert!(check_reference_position("c1", "c2", 2000.0, 10.0, "p", 1, &mut get).unwrap());
        assert!(!check_reference_position("c1", "c2", 5000.0, 10.0, "p", 1, &mut get).unwrap());
        // circular: distance through the end of the genome
        assert!(check_reference_position("c1", "c2", 97500.0, 10.0, "p", 1, &mut get).unwrap());
        // a contig mapped more than 1.5 times is a repeat
        let repeat =
            format!("{coords}1\t1000\t9000\t10000\t0\t0\t99\t1000\t100000\t60.0\t1\tc1\tg\n");
        let hits = std::rc::Rc::new(parse_coords(&repeat));
        let mut get = |_: &str, _: usize| Ok(hits.clone());
        assert!(!check_reference_position("c1", "c2", 2000.0, 10.0, "p", 1, &mut get).unwrap());
    }

    #[test]
    fn greedy_merging_and_export() {
        let species: HashMap<String, Vec<String>> = [
            ("1", vec!["d/Ecoli__a", "d/Shig__b"]),
            ("2", vec!["d/Ecoli__c", "d/Shig__b"]),
            ("3", vec!["d/Ecoli__a"]),
            ("4", vec!["d/Kleb__k"]),
            ("5", vec!["d/Kleb__k"]),
        ]
        .iter()
        .map(|(k, v)| (k.to_string(), v.iter().map(|s| s.to_string()).collect()))
        .collect();
        let mut merge = Merge::default();
        let mut log = Vec::new();
        for (a, b) in [("1", "2"), ("3", "1"), ("4", "3"), ("4", "5"), ("2", "5")] {
            merge.merge(a.into(), b.into(), &species, &mut log).unwrap();
        }
        assert_eq!(
            merge.supers,
            [
                Some(vec!["1".into(), "2".into(), "3".into()]),
                Some(vec!["4".into(), "5".into()])
            ]
        );
        assert_eq!(merge.species[&0], ["Ecoli"]);
        let info = contigs(&[("b", "1", "10"), ("a", "4", "5"), ("c", "6", "7.5")]);
        let (clusters, supers, out) = export_info(&info, &species, &merge);
        assert_eq!(
            clusters,
            "a\t1S\t001100\t0.0088888\nb\t0S\t001100\t0.0088888\nc\t6\t001100\t7.5\n"
        );
        assert_eq!(supers, "0S\t1\t2\t3\n1S\t4\t5\n");
        assert!(out.starts_with(">0S\nd/Ecoli__a\nd/Shig__b\nd/Ecoli__c\n>1S\nd/Kleb__k\n>1\n"));
    }

    #[test]
    fn edges_recovered_by_distance() {
        let info = contigs(&[("a", "1", "1"), ("b", "2", "1"), ("c", "3", "1")]);
        let mut merge = Merge::default();
        merge.super_of.insert("1".into(), 0);
        merge.super_of.insert("2".into(), 0);
        let good = "a\t+\tb\t+\t250\t1\t1\na\t+\tc\t+\t250\t1\t1\na\t+\tb\t+\t50000\t1\t1\n";
        let above = "a\t+\ta\t+\t1500\t1\t2\na\t+\tb\t+\t1500\t1\t5\nx\t+\ty\t+\t-5\t1\t2\nc\t+\tc\t+\t1\t1\t1\n";
        let mut too_long = 0;
        let sets = recover_edges(good, above, &info, &merge, &mut too_long).unwrap();
        // x and y are unclustered: not "the same cluster"
        assert_eq!(sets[0], "a\t+\tb\t+\t250\t1\t1\n");
        assert_eq!(sets[2], "a\t+\ta\t+\t1500\t1\t2\n");
        // a 40 kb edge is skipped and counted
        let sets = recover_edges(
            "",
            "a\t+\ta\t+\t40000\t1\t2\n",
            &info,
            &merge,
            &mut too_long,
        )
        .unwrap();
        assert!(sets.iter().all(String::is_empty));
        assert_eq!(too_long, 1);
    }

    #[test]
    fn genome_path_inside_any_database_folder() {
        assert_eq!(
            genome_path("/db/legacy/", "OPERA-MS-DB/genomes_2/g.fna.gz"),
            "/db/legacy/genomes_2/g.fna.gz"
        );
        assert_eq!(
            genome_path("/x/OPERA-MS-DB", "genomes_1/g.fna.gz"),
            "/x/OPERA-MS-DB/genomes_1/g.fna.gz"
        );
    }
}
