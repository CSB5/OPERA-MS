//! Long-read edges: maps the long reads on the contigs (minimap2), keeps
//! the useful alignments of each read, derives the contig pairs (edges) linked by the reads, bundles them
//! by orientation and distance, and writes the OPERA-LG inputs:
//! `opera.map`, `opera.map.sort`, `opera.map.sort.status`, `pairedEdges`,
//! `pairedEdges_i0..5`, `edge_read_info.dat`, `opera.map.cov`, `config`,
//! `results/contigs`.
//!
//! The behaviour of OPERA-MS 0.9 is kept, including its quirks (see the
//! comments: the last read of the mapping is never used for edges, a small
//! contig of a new read joins the previous read's list, the score of the
//! previous read carries over, `check_overlap` receives the contig start as
//! its score). Hash iterations that reach the outputs use sorted keys. The
//! mapping is sorted exactly as GNU `sort -k1,1 -k<n>,<n>g` sorts in the C
//! locale.

use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs;
use std::path::Path;

use anyhow::{Context, Result, ensure};
use opera_ms_core::bench::StepRecorder;
use opera_ms_core::runner::{Cmd, Output, Pipeline};

use crate::edge_filter::perl_number as perl_value;
use crate::files::perl_number as perl_format;

const MIN_OPERA_CONTIG_SIZE: f64 = 500.0;
const MIN_CONTIG_SIZE: f64 = 500.0;
const MIN_CONTIG_SIZE_FOR_GAP_FILLING: f64 = 100.0;
const FRACTION: f64 = 0.9;
const FRACTION_IN_GAP: f64 = 0.8;
const MIN_ALIGNMENT_LENGTH: f64 = 400.0;
const OVERLAP: f64 = 200.0;
const LONG_READ_CLUSTER_THRESHOLD: u32 = 2;
const MEANS: [u32; 6] = [300, 1000, 2000, 5000, 15000, 40000];
const STDS: [u32; 6] = [30, 100, 200, 500, 1500, 4000];
/// The distance intervals of `pairedEdges_i0..5` (open intervals).
const INTERVALS: [(f64, f64); 6] = [
    (-200.0, 300.0),
    (300.0, 1000.0),
    (1000.0, 2000.0),
    (2000.0, 5000.0),
    (5000.0, 15000.0),
    (15000.0, 40000.0),
];

/// Arguments of the long-read step.
pub struct Args<'a> {
    pub contig_file: &'a str,
    pub kmer: &'a str,
    pub long_read_file: &'a str,
    /// `--output-directory` (a `/` is appended if missing).
    pub output_dir: &'a str,
    pub nproc: &'a str,
    /// `--minimap2`: minimap2's directory.
    pub minimap2_dir: &'a str,
}

/// A value of an alignment: its text as Perl would print it, and its
/// numeric value.
#[derive(Debug, Clone)]
struct Val {
    text: String,
    num: f64,
}

impl Val {
    fn text(text: &str) -> Self {
        Self {
            text: text.to_string(),
            num: perl_value(text),
        }
    }
    fn num(num: f64) -> Self {
        Self {
            text: perl_format(num),
            num,
        }
    }
}

/// `@data = ($rn, $cn, $ro, $co, $cs, $ce, $cl, $rs, $re, $rl, $score)`.
#[derive(Debug, Clone)]
struct Data([Val; 11]);

impl Data {
    fn get(&self, i: usize) -> &Val {
        &self.0[i]
    }
    /// `"@data"`: the values joined with spaces.
    fn spaced(&self) -> String {
        self.0
            .iter()
            .map(|v| v.text.as_str())
            .collect::<Vec<_>>()
            .join(" ")
    }
    fn tabbed(&self) -> String {
        self.0
            .iter()
            .map(|v| v.text.as_str())
            .collect::<Vec<_>>()
            .join("\t")
    }
}

const RN: usize = 0;
const CS: usize = 4;
const CE: usize = 5;
const CL: usize = 6;
const RS: usize = 7;
const RE: usize = 8;
const RL: usize = 9;
const SCORE: usize = 10;

/// `split /\s+/` of a chomped line (Perl keeps a leading empty field).
fn split_ws(line: &str) -> Vec<&str> {
    let leading = line.starts_with(|c: char| c.is_whitespace());
    let mut fields: Vec<&str> = line.split_whitespace().collect();
    if leading {
        fields.insert(0, "");
    }
    fields
}

/// One line of the mapping as `checkMapping` reads it.
fn parse_alignment(line: &str) -> Result<Option<Data>> {
    let f = split_ws(line);
    if f.len() < 10 {
        return Ok(None);
    }
    let t = |i: usize| Val::text(f.get(i).copied().unwrap_or(""));
    let data = {
        // ($rn, $rl, $rs, $re, $ro, $cn, $cl, $cs, $ce)
        let (rn, rl, rs, re, cn, cl) = (t(0), t(1), t(2), t(3), t(5), t(6));
        let (mut cs, mut ce) = (t(7), t(8));
        let strand = f[4];
        let (ro, co) = if strand == "-" {
            let (cs_rev, ce_rev) = (cl.num - ce.num, cl.num - cs.num);
            cs = Val::num(cs_rev);
            ce = Val::num(ce_rev);
            (Val::num(0.0), Val::num(1.0))
        } else if strand == "+" {
            (Val::num(0.0), Val::num(0.0))
        } else {
            (Val::text(strand), Val::num(0.0))
        };
        let span = ce.num - cs.num;
        ensure!(
            span != 0.0,
            "mapping line with an empty contig alignment: {line}"
        );
        let score = Val::num(perl_value(f[9]) / span);
        [rn, cn, ro, co, cs, ce, cl, rs, re, rl, score]
    };
    Ok(Some(Data(data)))
}

/// `getAlignmentType`.
fn alignment_type(d: &Data) -> &'static str {
    let (cs, ce, cl) = (d.get(CS).num, d.get(CE).num, d.get(CL).num);
    let (rs, re, rl) = (d.get(RS).num, d.get(RE).num, d.get(RL).num);
    if ce - cs >= FRACTION * cl || ce - cs >= cl - OVERLAP {
        "contig-contained"
    } else if re - rs >= FRACTION * rl || re - rs >= rl - OVERLAP {
        "read-contained"
    } else if ce >= cl - OVERLAP && rs <= OVERLAP {
        "contig-at-start"
    } else if cs <= OVERLAP && re >= rl - OVERLAP {
        "contig-at-end"
    } else {
        "partial-match"
    }
}

/// An alignment of `@alignments`: `"@data"`, plus the tag Perl appends.
#[derive(Debug, Clone)]
struct Alignment {
    text: String,
    small_contig: bool,
}

/// The fields `printEdges` reads from an alignment string
/// (`split(/ /, ...)`): the first ten (the score is never used).
struct Fields {
    rn: String,
    cn: String,
    ro: String,
    co: String,
    cs: String,
    ce: String,
    rs_t: String,
    re_t: String,
    cl: f64,
    rs: f64,
    re: f64,
    cs_n: f64,
    ce_n: f64,
}

impl Fields {
    fn of(text: &str) -> Self {
        let f: Vec<&str> = text.split(' ').collect();
        let s = |i: usize| f.get(i).copied().unwrap_or("").to_string();
        let n = |i: usize| perl_value(f.get(i).copied().unwrap_or(""));
        Self {
            rn: s(0),
            cn: s(1),
            ro: s(2),
            co: s(3),
            cs: s(4),
            ce: s(5),
            rs_t: s(7),
            re_t: s(8),
            cl: n(6),
            rs: n(7),
            re: n(8),
            cs_n: n(4),
            ce_n: n(5),
        }
    }
    fn ori(&self) -> &'static str {
        if perl_value(&self.co) == 0.0 {
            "+"
        } else {
            "-"
        }
    }
}

/// An entry of `@contig_for_gapfilling_list`: `($name, $ori, $start, $end, $score)`.
struct GapContig {
    name: String,
    ori: &'static str,
    end: f64,
    score: f64,
}

/// `check_overlap_contig_for_gapfilling($list, $start, $end, $score, $size)`.
/// As in OPERA-MS 0.9: the "previous size" is the last element of the
/// 5-element entry (its score) and the "previous score" its end.
fn check_overlap(list: &mut Vec<GapContig>, start: f64, score: f64, size: f64) -> bool {
    let Some(last) = list.last() else { return true };
    if last.end - OVERLAP > start {
        let (prev_size, prev_score) = (last.score, last.end);
        if size >= MIN_CONTIG_SIZE || (prev_size < MIN_CONTIG_SIZE && prev_score > score) {
            list.pop();
            true
        } else {
            false
        }
    } else {
        true
    }
}

/// `add_contig_for_gapfilling`.
#[allow(clippy::too_many_arguments)]
fn add_gap_contig(
    list: &mut Vec<GapContig>,
    name: &str,
    ori: &'static str,
    start: f64,
    end: f64,
    score: f64,
    size: f64,
    check: bool,
) {
    if !check || check_overlap(list, start, score, size) {
        list.push(GapContig {
            name: name.to_string(),
            ori,
            end,
            score,
        });
    }
}

#[derive(Default)]
struct EdgeInfo {
    edge: String,
    dist: Vec<f64>,
    reads: Vec<String>,
    coord1: Vec<String>,
    coord2: Vec<String>,
    coord1_on_read: Vec<String>,
    coord2_on_read: Vec<String>,
    gapfilling: Vec<String>,
}

/// Perl's `int()`: truncation toward zero.
fn int(x: f64) -> f64 {
    x.trunc()
}

/// `printEdges`: the edges between the alignments of one read.
fn print_edges(alignments: &[Alignment], edges: &mut BTreeMap<String, EdgeInfo>) {
    for (i, a1) in alignments.iter().enumerate() {
        if a1.small_contig {
            continue;
        }
        let f1 = Fields::of(&a1.text);
        let mut gap_list: Vec<GapContig> = Vec::new();
        for a2 in &alignments[i + 1..] {
            let f2 = Fields::of(&a2.text);
            let ori2 = f2.ori();
            if a2.small_contig {
                if f1.re - OVERLAP < f2.rs {
                    add_gap_contig(
                        &mut gap_list,
                        &f2.cn,
                        ori2,
                        f2.rs,
                        f2.re,
                        f2.cs_n,
                        f2.cl,
                        true,
                    );
                }
                continue;
            }
            let gap = f2.rs - f1.re;
            let mut distance = int(gap + gap.abs() * 0.09);
            distance += -f2.cs_n - (f1.cl - f1.ce_n);
            if f1.cn == f2.cn || f1.cl < MIN_OPERA_CONTIG_SIZE || f2.cl < MIN_OPERA_CONTIG_SIZE {
                continue;
            }
            let first_is_1 = f1.cn <= f2.cn;
            let edge_id = if first_is_1 {
                format!("{} {}", f1.cn, f2.cn)
            } else {
                format!("{} {}", f2.cn, f1.cn)
            };
            let info = edges.entry(edge_id).or_insert_with(|| EdgeInfo {
                edge: if first_is_1 {
                    format!("{}\t{}", f1.cn, f2.cn)
                } else {
                    format!("{}\t{}", f2.cn, f1.cn)
                },
                ..Default::default()
            });
            let c1 = format!("{}_{}_{}", f1.co, f1.cs, f1.ce);
            let c2 = format!("{}_{}_{}", f2.co, f2.cs, f2.ce);
            let r1 = format!("{}_{}_{}", f1.ro, f1.rs_t, f1.re_t);
            let r2 = format!("{}_{}_{}", f2.ro, f2.rs_t, f2.re_t);
            if first_is_1 {
                info.coord1.push(c1);
                info.coord2.push(c2);
                info.coord1_on_read.push(r1);
                info.coord2_on_read.push(r2);
            } else {
                info.coord1.push(c2);
                info.coord2.push(c1);
                info.coord1_on_read.push(r2);
                info.coord2_on_read.push(r1);
            }
            info.reads.push(f1.rn.clone());
            info.dist.push(distance);
            // as in OPERA-MS 0.9: called with the contig start as "score"
            check_overlap(&mut gap_list, f2.rs, f2.cs_n, f2.cl);
            let names: Vec<String> = gap_list
                .iter()
                .map(|c| format!("{}{}", c.ori, c.name))
                .collect();
            info.gapfilling.push(if names.is_empty() {
                "NONE".to_string()
            } else {
                names.join(",")
            });
            add_gap_contig(
                &mut gap_list,
                &f2.cn,
                ori2,
                f2.rs,
                f2.re,
                f2.cs_n,
                f2.cl,
                false,
            );
        }
    }
}

fn add_sd(d: f64) -> f64 {
    d + 6.0 * (d.abs() * 0.1 + 50.0)
}

fn remove_sd(d: f64) -> f64 {
    d - 6.0 * (d.abs() * 0.1 + 50.0)
}

/// `checkMapping`: returns (`opera.map.sort.status`, `pairedEdges`,
/// `edge_read_info.dat`).
pub fn check_mapping(map_sorted: &str) -> Result<(String, String, String)> {
    let mut status = String::new();
    let mut edges: BTreeMap<String, EdgeInfo> = BTreeMap::new();
    let mut alignments: Vec<Alignment> = Vec::new();
    let mut prev_rn = String::new();
    let mut prev_alignment_end = 0.0;
    let (mut current_score, mut previous_score) = (0.0, 0.0);
    for line in map_sorted.lines() {
        let Some(d) = parse_alignment(line)? else {
            continue;
        };
        let (cs, ce, cl) = (d.get(CS).num, d.get(CE).num, d.get(CL).num);
        let (rs, re) = (d.get(RS).num, d.get(RE).num);
        if cl < MIN_CONTIG_SIZE {
            status.push_str(&format!("{} | small-contig\n", d.tabbed()));
            if cl > MIN_CONTIG_SIZE_FOR_GAP_FILLING && ce - cs >= FRACTION_IN_GAP * cl {
                alignments.push(Alignment {
                    text: d.spaced(),
                    small_contig: true,
                });
            }
            continue;
        }
        if re - rs < MIN_ALIGNMENT_LENGTH {
            status.push_str(&format!("{} | small-alignment\n", d.tabbed()));
            continue;
        }
        let kind = alignment_type(&d);
        status.push_str(&format!("{} | {kind}", d.tabbed()));
        if kind == "partial-match" {
            status.push('\n');
            continue;
        }
        let rn = &d.get(RN).text;
        if *rn != prev_rn {
            if alignments.len() > 1 {
                print_edges(&alignments, &mut edges);
            }
            prev_rn = rn.clone();
            alignments = vec![Alignment {
                text: d.spaced(),
                small_contig: false,
            }];
            prev_alignment_end = re + cl - ce;
            // as in OPERA-MS 0.9: the last score of the previous read
            previous_score = current_score;
        } else {
            let curr_start = rs - cs;
            let curr_end = re + cl - ce;
            current_score = d.get(SCORE).num;
            if prev_alignment_end - OVERLAP < curr_start {
                alignments.push(Alignment {
                    text: d.spaced(),
                    small_contig: false,
                });
                prev_alignment_end = curr_end;
                previous_score = current_score;
            } else {
                if previous_score < current_score {
                    alignments.pop();
                    alignments.push(Alignment {
                        text: d.spaced(),
                        small_contig: false,
                    });
                    prev_alignment_end = curr_end;
                    previous_score = current_score;
                    status.push_str(" overlapped better");
                } else {
                    status.push_str(" overlapped");
                }
            }
        }
        status.push('\n');
    }
    // as in OPERA-MS 0.9: the alignments of the last read are not used

    let (mut paired, mut edge_read) = (String::new(), String::new());
    for info in edges.values_mut() {
        let mut edge_types = Vec::with_capacity(info.dist.len());
        let mut all: Vec<(f64, String)> = Vec::with_capacity(info.dist.len());
        for i in 0..info.dist.len() {
            let first = |s: &str| s.split('_').next().map(perl_value).unwrap_or(0.0);
            let second = |s: &str| s.split('_').nth(1).map(perl_value).unwrap_or(0.0);
            let mut ori_1 = if first(&info.coord1[i]) == 0.0 {
                "+"
            } else {
                "-"
            };
            let mut ori_2 = if first(&info.coord2[i]) == 0.0 {
                "+"
            } else {
                "-"
            };
            if second(&info.coord1_on_read[i]) > second(&info.coord2_on_read[i]) {
                ori_1 = if ori_1 == "+" { "-" } else { "+" };
                ori_2 = if ori_2 == "+" { "-" } else { "+" };
            }
            let edge = format!("{ori_1}:{ori_2}");
            all.push((info.dist[i], edge.clone()));
            edge_types.push(edge);
        }
        // a stable sort: equal distances keep their order
        all.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(Ordering::Equal));
        // edge type -> [first distance, distance sum, support]
        let mut support: BTreeMap<String, Vec<(f64, f64, u32)>> = BTreeMap::new();
        for (distance, edge) in &all {
            let groups = support.entry(edge.clone()).or_default();
            match groups
                .iter_mut()
                .find(|g| remove_sd(g.0) <= *distance && *distance <= add_sd(g.0))
            {
                Some(g) => {
                    g.1 += distance;
                    g.2 += 1;
                }
                None => groups.push((*distance, *distance, 1)),
            }
        }
        // the first key (sorted) is the starting point; strict < keeps it on ties
        let Some(first_edge) = support.keys().next().cloned() else {
            continue;
        };
        let (mut best_edge, mut best_id) = (first_edge.clone(), 0usize);
        let mut best_support = support[&first_edge][0].2;
        for (edge, groups) in &support {
            for (j, g) in groups.iter().enumerate() {
                if best_support < g.2 {
                    best_support = g.2;
                    best_id = j;
                    best_edge = edge.clone();
                }
            }
        }
        let g = support[&best_edge][best_id];
        let best_distance = g.1 / f64::from(g.2);
        let sd = int(best_distance.abs() * 0.1 + 50.0);
        let mut oris = best_edge.split(':');
        let (o1, o2) = (oris.next().unwrap_or(""), oris.next().unwrap_or(""));
        let mut ids = info.edge.split('\t');
        let (id1, id2) = (
            ids.next().unwrap_or("").to_string(),
            ids.next().unwrap_or("").to_string(),
        );
        info.edge = format!("{id1}\t{o1}\t{id2}\t{o2}");
        paired.push_str(&format!(
            "{}\t{}\t{}\t{best_support}\n",
            info.edge,
            perl_format(int(best_distance)),
            perl_format(sd)
        ));
        // keep only the reads supporting the chosen edge (the distance and
        // type lists aren't spliced, as in OPERA-MS 0.9)
        let mut offset = 0;
        for (d, edge_type) in info.dist.iter().zip(&edge_types) {
            let d = *d;
            if best_edge != *edge_type
                || !(remove_sd(best_distance) <= d && d <= add_sd(best_distance))
            {
                for list in [
                    &mut info.reads,
                    &mut info.coord1,
                    &mut info.coord2,
                    &mut info.coord1_on_read,
                    &mut info.coord2_on_read,
                    &mut info.gapfilling,
                ] {
                    if offset < list.len() {
                        list.remove(offset);
                    }
                }
            } else {
                offset += 1;
            }
        }
        edge_read.push_str(&format!(
            "{}\t{}\tCOORD_CONTIG_1:{}\tCOORD_CONTIG_2:{}\tCOORD_CONTIG_1_ON_READ:{}\t\
             COORD_CONTIG_2_ON_READ:{}\tCONTIG_FOR_GAPFILLING:{}\n",
            info.edge,
            info.reads.join(";"),
            info.coord1.join(";"),
            info.coord2.join(";"),
            info.coord1_on_read.join(";"),
            info.coord2_on_read.join(";"),
            info.gapfilling.join(";")
        ));
    }
    Ok((status, paired, edge_read))
}

/// `long_read_coverage_estimate`: mapped bases per contig length, from the
/// status file (without the `small-alignment` lines). Contigs sorted.
pub fn long_read_coverage(status: &str) -> String {
    let mut map_length: BTreeMap<&str, f64> = BTreeMap::new();
    let mut contig_length: HashMap<&str, &str> = HashMap::new();
    let mut on_read: HashSet<&str> = HashSet::new();
    let mut current_read = "";
    for line in status
        .split_inclusive('\n')
        .filter(|l| !l.contains("small-alignment"))
    {
        let f: Vec<&str> = line.split('\t').collect();
        let get = |i: usize| f.get(i).copied().unwrap_or("");
        let (read, contig) = (get(0), get(1));
        contig_length.insert(contig, get(6));
        if current_read != read {
            on_read.clear();
            current_read = read;
        }
        if on_read.insert(contig) {
            *map_length.entry(contig).or_insert(0.0) += perl_value(get(8)) - perl_value(get(7));
        }
    }
    map_length
        .iter()
        .map(|(c, len)| {
            let cl = contig_length.get(c).copied().unwrap_or("");
            format!("{c}\t{cl}\t{}\n", perl_format(len / perl_value(cl)))
        })
        .collect()
}

/// `extract_edge`: `pairedEdges` split by distance interval.
pub fn extract_edges(paired: &str) -> [String; 6] {
    let mut out: [String; 6] = Default::default();
    for line in paired.split_inclusive('\n') {
        let mut chars = line.chars();
        chars.next_back(); // chop
        let mut fields: Vec<&str> = chars.as_str().split('\t').collect();
        while fields.last().is_some_and(|f| f.is_empty()) {
            fields.pop();
        }
        let dist = perl_value(fields.get(4).copied().unwrap_or(""));
        let support = perl_value(fields.get(6).copied().unwrap_or(""));
        if let Some(i) = INTERVALS
            .iter()
            .position(|(lo, hi)| support >= 0.0 && *lo < dist && dist < *hi)
        {
            out[i].push_str(&fields.join("\t"));
            out[i].push('\n');
        }
    }
    out
}

/// `CreateConfigFile` (no short-read mapping: no BAM library).
pub fn opera_config(contig_file: &str, output_dir: &str, kmer: &str) -> String {
    let mut out = format!(
        "#\n# Essential Parameters\n#\n\n# Output folder for final results\n\
         output_folder={output_dir}/results\n# Contig file\ncontig_file={contig_file}\nkmer={kmer}\n\
         # Mapped read locations\n"
    );
    for i in 0..6 {
        out.push_str(&format!(
            "[LIB]\ncluster_threshold={LONG_READ_CLUSTER_THRESHOLD}\nmap_file={output_dir}/pairedEdges_i{i}\n\
             lib_mean={}\nlib_std={}\nmap_type=opera\n",
            MEANS[i], STDS[i]
        ));
    }
    out
}

/// GNU `sort -k1,1 -k<n>,<n>g` in the C locale: field 1 bytewise, field n
/// as a general number (non-numbers first), then the whole line bytewise.
pub fn sort_mapping(text: &str, numeric_field: usize) -> String {
    let key = |line: &str| -> (String, Option<f64>) {
        let fields: Vec<&str> = line.split(['\t', ' ']).filter(|f| !f.is_empty()).collect();
        let first = line.split(['\t', ' ']).next().unwrap_or("").to_string();
        let numeric = fields
            .get(numeric_field - 1)
            .and_then(|f| general_number(f));
        (first, numeric)
    };
    let mut lines: Vec<(&str, (String, Option<f64>))> = text
        .split_inclusive('\n')
        .map(|l| (l, key(l.trim_end_matches('\n'))))
        .collect();
    lines.sort_by(|(la, (fa, na)), (lb, (fb, nb))| {
        fa.as_bytes()
            .cmp(fb.as_bytes())
            .then_with(|| match (na, nb) {
                (None, None) => Ordering::Equal,
                (None, Some(_)) => Ordering::Less,
                (Some(_), None) => Ordering::Greater,
                (Some(x), Some(y)) => x.partial_cmp(y).unwrap_or(Ordering::Equal),
            })
            .then_with(|| {
                la.trim_end_matches('\n')
                    .as_bytes()
                    .cmp(lb.trim_end_matches('\n').as_bytes())
            })
    });
    lines.into_iter().map(|(l, _)| l).collect()
}

/// `strtold` on a field (leading number), `None` if it isn't a number.
pub(crate) fn general_number(field: &str) -> Option<f64> {
    let t = field.trim_start();
    let mut end = 0;
    for i in (1..=t.len()).rev() {
        if t.is_char_boundary(i) && t[..i].parse::<f64>().is_ok() {
            end = i;
            break;
        }
    }
    if end == 0 {
        None
    } else {
        t[..end].parse().ok()
    }
}

/// `cut -f1-21`.
pub fn cut_mapping(raw: &str) -> String {
    let mut out = String::with_capacity(raw.len());
    for line in raw.split_inclusive('\n') {
        let (body, newline) = match line.strip_suffix('\n') {
            Some(b) => (b, "\n"),
            None => (line, ""),
        };
        out.push_str(&body.split('\t').take(21).collect::<Vec<_>>().join("\t"));
        out.push_str(newline);
    }
    out
}

fn write(path: &str, text: &str) -> Result<()> {
    fs::write(path, text).with_context(|| format!("cannot write {path}"))
}

/// The long-read mapping and the OPERA-LG inputs.
pub fn opera_long_read(rec: &StepRecorder<'_>, a: &Args<'_>) -> Result<()> {
    let out = if a.output_dir.ends_with('/') {
        a.output_dir.to_string()
    } else {
        format!("{}/", a.output_dir)
    };
    fs::create_dir_all(&out).with_context(|| format!("cannot create {out}"))?;
    for (what, path) in [("contig", a.contig_file), ("long read", a.long_read_file)] {
        ensure!(
            Path::new(path).exists(),
            "Error: {path} - {what} file does not exist"
        );
    }
    let file = |name: &str| format!("{out}{name}");
    let sorted_path = file("opera.map.sort");
    if !Path::new(&sorted_path).exists() {
        let raw = file("opera.map.raw");
        let program = format!("{}/minimap2", a.minimap2_dir);
        let args = [
            "-t",
            a.nproc,
            "-w5",
            "-m0",
            "--cs=short",
            a.contig_file,
            a.long_read_file,
        ];
        let err = file("minimap2.err");
        // the mapping is sorted by read name, then by the read start
        let numeric_field = 3;
        let pipeline = Pipeline::new(Cmd::new(&program).args(args))
            .stdout(Output::File(raw.clone().into()))
            .stderr(Output::File(err.clone().into()));
        let report = pipeline.run().with_context(|| {
            format!("Error in the long-read mapping. Please see {err} for details")
        })?;
        rec.record(&report);
        let raw_text = fs::read_to_string(&raw).with_context(|| format!("cannot read {raw}"))?;
        let map = cut_mapping(&raw_text);
        drop(raw_text);
        // only the sorted mapping is kept
        fs::remove_file(&raw).with_context(|| format!("cannot remove {raw}"))?;
        write(&sorted_path, &sort_mapping(&map, numeric_field))?;
    }
    let sorted =
        fs::read_to_string(&sorted_path).with_context(|| format!("cannot read {sorted_path}"))?;
    let (status, paired, edge_read) = check_mapping(&sorted)?;
    drop(sorted);
    write(&file("opera.map.sort.status"), &status)?;
    write(&file("pairedEdges"), &paired)?;
    write(&file("edge_read_info.dat"), &edge_read)?;
    write(&file("opera.map.cov"), &long_read_coverage(&status))?;
    for (i, text) in extract_edges(&paired).iter().enumerate() {
        write(&file(&format!("pairedEdges_i{i}")), text)?;
    }
    write(&file("config"), &opera_config(a.contig_file, &out, a.kmer))?;
    let results = file("results");
    if Path::new(&results).is_dir() {
        fs::remove_dir_all(&results).with_context(|| format!("cannot remove {results}"))?;
    }
    fs::create_dir(&results).with_context(|| format!("cannot create {results}"))?;
    let link = format!("{results}/contigs");
    std::os::unix::fs::symlink("../opera.map.cov", &link)
        .with_context(|| format!("cannot link {link}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    const RUN: &str = "/data/ecoli_run/intermediate_files";

    #[test]
    fn config_matches_the_reference_run() {
        let fixture = std::fs::read_to_string(format!(
            "{}/../opera-ms-io/tests/fixtures/cpp_config/read_mapping.config",
            env!("CARGO_MANIFEST_DIR")
        ))
        .unwrap();
        let config = opera_config(
            &format!("{RUN}/user_assembly/contigs.fasta"),
            &format!("{RUN}/read_mapping/"),
            "60",
        );
        assert_eq!(config, fixture);
    }

    #[test]
    fn sort_like_gnu_sort() {
        let text = "b\tx\t10\n\
                    a\tx\t9\n\
                    a\tx\t10\n\
                    a\ty\t10\n\
                    a\tx\tz\n";
        assert_eq!(
            sort_mapping(text, 3),
            "a\tx\tz\na\tx\t9\na\tx\t10\na\ty\t10\nb\tx\t10\n"
        );
    }

    #[test]
    fn cut_fields() {
        let line = (1..=25)
            .map(|i| i.to_string())
            .collect::<Vec<_>>()
            .join("\t")
            + "\n";
        let cut = cut_mapping(&line);
        assert_eq!(cut.trim_end().split('\t').count(), 21);
    }

    #[test]
    fn minimap2_reverse_strand_and_score() {
        let d = parse_alignment("r1\t1000\t0\t900\t-\tc1\t5000\t100\t1000\t450\t900\t60")
            .unwrap()
            .unwrap();
        // cs = cl - ce = 4000, ce = cl - cs = 4900; score = 450 / 900
        assert_eq!(
            d.tabbed(),
            "r1\tc1\t0\t1\t4000\t4900\t5000\t0\t900\t1000\t0.5"
        );
        assert_eq!(alignment_type(&d), "read-contained");
    }

    #[test]
    fn two_contigs_on_two_reads_make_one_edge() {
        // reads r1 and r2 both link c1 (+) and c2 (+); a third read r3 is
        // the last read: never used (as in OPERA-MS 0.9)
        let map = "r1\t10000\t0\t3000\t+\tc1\t3000\t0\t3000\t2900\t3000\t60\n\
                   r1\t10000\t4000\t7000\t+\tc2\t3000\t0\t3000\t2900\t3000\t60\n\
                   r2\t10000\t0\t3000\t+\tc1\t3000\t0\t3000\t2900\t3000\t60\n\
                   r2\t10000\t4100\t7100\t+\tc2\t3000\t0\t3000\t2900\t3000\t60\n\
                   r3\t10000\t0\t3000\t+\tc1\t3000\t0\t3000\t2900\t3000\t60\n";
        let (status, paired, edge_read) = check_mapping(map).unwrap();
        assert_eq!(status.lines().count(), 5);
        assert!(status.starts_with(
            "r1\tc1\t0\t0\t0\t3000\t3000\t0\t3000\t10000\t0.966666666666667 | contig-contained\n"
        ));
        // only r1's edge is printed: r2's edges wait for the next read (r3)
        // ... and r2 is processed when r3 starts; r3 itself is never used
        // r1: distance int(1000 + 90) = 1090; r2: int(1100 + 99) = 1199
        assert_eq!(paired, "c1\t+\tc2\t+\t1144\t164\t2\n");
        assert!(edge_read.starts_with("c1\t+\tc2\t+\tr1;r2\tCOORD_CONTIG_1:0_0_3000;0_0_3000\t"));
        assert_eq!(long_read_coverage(&status), "c1\t3000\t3\nc2\t3000\t2\n");
        let split = extract_edges(&paired);
        assert_eq!(split[2], "c1\t+\tc2\t+\t1144\t164\t2\n");
        assert!(split[0].is_empty());
    }
}
