//! Gap filling (stage 8).
//!
//! 1. For the scaffolds with several contigs above the length threshold,
//!    the long reads of their edges are pooled (`READ_POOLS/`, a bounded
//!    number of files) and a pre-consensus is made from the contigs and the
//!    read sequence of each gap (`pre_consensus.fa`, all scaffolds;
//!    `extract_read_sequence`).
//! 2. minimap2 + racon polish each pre-consensus (`scaffold_racon`), in a
//!    directory per scaffold that only exists while its job runs.
//! 3. Two tilings: the contigs are mapped on the consensus with minimap2 and
//!    the contig sequences put back in place of the consensus (`remap`).
//! 4. `scaffoldSeq.fasta.filled`: filled scaffolds, the others unchanged.
//!
//! A gap whose read sequence is not found is filled with Ns. Directory
//! listings are used in file-system order; a read's edges in sorted order.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::fmt::Write as _;
use std::fs::{self, OpenOptions};
use std::io::{BufRead, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;
use std::sync::Mutex;
use std::time::Instant;

use anyhow::{Context, Result, bail, ensure};
use flate2::Compression;
use flate2::write::GzEncoder;
use opera_ms_core::bench::StepRecorder;
use opera_ms_core::parallel::run_parallel;
use opera_ms_core::runner::{Cmd, Output, Pipeline};

use crate::edge_filter::perl_number as perl_value;
use crate::files::perl_number as perl_format;

/// Arguments of the gap filling.
pub struct Args<'a> {
    pub working_dir: &'a str,
    pub contig_file: &'a str,
    pub read_file: &'a str,
    pub read_mapping_dir: &'a str,
    pub nproc: usize,
    pub scaffold_length_threshold: &'a str,
    pub racon_dir: &'a str,
    pub minimap2_dir: &'a str,
}

fn read(path: &str) -> Result<String> {
    fs::read_to_string(path).with_context(|| format!("cannot read {path}"))
}

fn write(path: &str, text: &str) -> Result<()> {
    fs::write(path, text).with_context(|| format!("cannot write {path}"))
}

fn append(path: &str, text: &[u8]) -> Result<()> {
    OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)
        .and_then(|mut f| f.write_all(text))
        .with_context(|| format!("cannot append to {path}"))
}

fn remove_if_exists(path: &str) -> Result<()> {
    match fs::symlink_metadata(path) {
        Ok(meta) if meta.is_dir() => fs::remove_dir_all(path),
        Ok(_) => fs::remove_file(path),
        Err(_) => return Ok(()),
    }
    .with_context(|| format!("cannot remove {path}"))
}

/// The line without its newline.
fn chop(text: &str) -> &str {
    text.strip_suffix('\n').unwrap_or(text)
}

/// Lines with their `\n`.
fn raw_lines(text: &str) -> impl Iterator<Item = &str> {
    text.split_inclusive('\n')
}

/// `split(/\s+/, $text)`: a leading empty field when the text starts with
/// whitespace, no trailing empty fields.
fn split_ws(text: &str) -> Vec<&str> {
    let leading = text.starts_with(char::is_whitespace);
    let mut fields: Vec<&str> = text.split_whitespace().collect();
    if leading && !fields.is_empty() {
        fields.insert(0, "");
    }
    fields
}

fn time_log(path: &str, step: &str, start: Instant) -> Result<()> {
    append(
        path,
        crate::files::time_log_line(step, start.elapsed().as_secs()).as_bytes(),
    )
}

/// Perl's `substr($s, $offset, $length)` with `$length >= 0` (`None`
/// where Perl gives undef).
fn perl_substr(s: &str, offset: i64, length: i64) -> Option<&str> {
    let n = s.len() as i64;
    let mut start = offset;
    if start < 0 {
        start += n;
        if start + length < 0 {
            return None;
        }
    } else if start > n {
        return None;
    }
    let end = (start + length).clamp(0, n);
    let start = start.clamp(0, n);
    s.get(start as usize..end.max(start) as usize)
}

/// Reverse complement, IUPAC codes in both cases.
fn reverse_complement(seq: &str) -> String {
    const FROM: &[u8] = b"ABCDGHKMNRSTUVWXYabcdghkmnrstuvwxy";
    const TO: &[u8] = b"TVGHCDMKNYSAABWXRtvghcdmknysaabwxr";
    seq.bytes()
        .rev()
        .map(|b| match FROM.iter().position(|&f| f == b) {
            Some(i) => TO[i] as char,
            None => b as char,
        })
        .collect()
}

// ---- split_scaffold_file ------------------------------------------------------

/// The scaffolds: order, the ones to fill (several contigs, longer than
/// the threshold).
struct Scaffolds {
    order: Vec<String>,
    multi: HashSet<String>,
    to_fill: Vec<String>,
}

fn split_scaffold_file(text: &str, threshold: f64) -> Scaffolds {
    let mut s = Scaffolds {
        order: Vec::new(),
        multi: HashSet::new(),
        to_fill: Vec::new(),
    };
    let (mut nb_contig, mut name, mut length) = (0usize, String::new(), 0.0);
    let close = |s: &mut Scaffolds, nb_contig: usize, name: &str, length: f64| {
        if nb_contig > 1 && length > threshold && s.multi.insert(name.to_string()) {
            s.to_fill.push(name.to_string());
        }
    };
    let mut started = false;
    for line in raw_lines(text) {
        let body = line.strip_suffix('\n').unwrap_or(line);
        match body.strip_prefix('>').filter(|r| !r.is_empty()) {
            Some(rest) => {
                if started && nb_contig > 0 {
                    close(&mut s, nb_contig, &name, length);
                }
                started = true;
                let fields = split_ws(rest);
                name = fields.first().copied().unwrap_or("").to_string();
                length = perl_value(fields.get(2).copied().unwrap_or(""));
                s.order.push(name.clone());
                nb_contig = 0;
            }
            None => nb_contig += 1,
        }
    }
    close(&mut s, nb_contig, &name, length);
    s
}

// ---- read pools and pre-consensus --------------------------------------------

struct EdgeInfo {
    scaff: String,
    /// The gap sequence; empty until the reference read is seen (empty at the
    /// end: the gap is filled with Ns).
    seq: String,
    ori: String,
}

/// `get_edge_ID`: contigs in sorted order, the orientation of the first.
fn edge_id(c1: &str, o1: &str, c2: &str, o2: &str) -> (String, String) {
    if c1 <= c2 {
        (format!("{c1}-vs-{c2}"), o1.to_string())
    } else {
        (format!("{c2}-vs-{c1}"), o2.to_string())
    }
}

/// The edges of the scaffolds (consecutive contigs).
fn scaffold_edges(scaf: &str) -> HashMap<String, EdgeInfo> {
    let mut edges = HashMap::new();
    let mut lines = raw_lines(scaf);
    let (mut c1, mut o1, mut scaff) = (String::new(), String::new(), String::new());
    while let Some(line) = lines.next() {
        let fields: Vec<&str> = chop(line).split('\t').collect();
        let first = fields.first().copied().unwrap_or("");
        let orientation = |f: &[&str]| if f.get(1) == Some(&"EB") { "-" } else { "+" }.to_string();
        if let Some(i) = first.find('>').filter(|i| first.len() > i + 1) {
            scaff = first[i + 1..].to_string();
            let next = lines.next().unwrap_or("");
            let f: Vec<&str> = chop(next).split('\t').collect();
            c1 = f.first().copied().unwrap_or("").to_string();
            o1 = orientation(&f);
        } else {
            let c2 = first.to_string();
            let o2 = orientation(&fields);
            let (id, ori) = edge_id(&c1, &o1, &c2, &o2);
            edges.insert(
                id,
                EdgeInfo {
                    scaff: scaff.clone(),
                    seq: String::new(),
                    ori,
                },
            );
            c1 = c2;
            o1 = o2;
        }
    }
    edges
}

/// A read's coordinates for an edge: start and end of the gap on the read,
/// reverse complement or not.
type ReadCoords = (f64, f64, i64);

fn coord_list(field: &str) -> Vec<&str> {
    field.split(':').nth(1).unwrap_or("").split(';').collect()
}

fn coord_values(entry: &str) -> Vec<f64> {
    entry.split('_').map(perl_value).collect()
}

/// The edges of `edge_read_info.dat` that are in a scaffold: the reads of
/// each edge and their gap coordinates; writes `contig_extention.log`.
fn read_edges(
    text: &str,
    edges: &HashMap<String, EdgeInfo>,
    contig_lengths: &HashMap<String, f64>,
    extension_log: &mut String,
) -> (
    HashMap<String, HashMap<String, ReadCoords>>,
    HashMap<String, String>,
) {
    let mut read_to_edge: HashMap<String, HashMap<String, ReadCoords>> = HashMap::new();
    let mut reference_reads = HashMap::new();
    for line in raw_lines(text) {
        let f: Vec<&str> = chop(line).split('\t').collect();
        let get = |i: usize| f.get(i).copied().unwrap_or("");
        let (c1, o1, c2, o2) = (get(0), get(1), get(2), get(3));
        let (id, ori) = edge_id(c1, o1, c2, o2);
        let Some(info) = edges.get(&id) else {
            continue;
        };
        let reverse_edge = i64::from(ori != info.ori);
        let reads: Vec<&str> = get(4).split(';').collect();
        let (cc1, cc2) = (coord_list(get(5)), coord_list(get(6)));
        let (cr1, cr2) = (coord_list(get(7)), coord_list(get(8)));
        let c1_r1 = coord_values(cc1.first().copied().unwrap_or(""));
        let c2_r1 = coord_values(cc2.first().copied().unwrap_or(""));
        let at = |v: &[f64], i: usize| v.get(i).copied().unwrap_or(0.0);
        let length = |c: &str| contig_lengths.get(c).copied().unwrap_or(0.0);
        reference_reads.insert(id.clone(), reads.first().copied().unwrap_or("").to_string());
        for (i, r) in reads.iter().enumerate() {
            let p1 = coord_values(cc1.get(i).copied().unwrap_or(""));
            let p2 = coord_values(cc2.get(i).copied().unwrap_or(""));
            let ori1 = if at(&p1, 0) == 0.0 { "+" } else { "-" };
            let ori2 = if at(&p2, 0) == 0.0 { "+" } else { "-" };
            let mut reverse = reverse_edge;
            if ori1 != o1 && ori2 != o2 {
                reverse = (reverse + 1) % 2;
            }
            let on1 = coord_values(cr1.get(i).copied().unwrap_or(""));
            let on2 = coord_values(cr2.get(i).copied().unwrap_or(""));
            let entry = read_to_edge.entry((*r).to_string()).or_default();
            if at(&on1, 2) < at(&on2, 2) {
                let end_extra = length(c1) - at(&c1_r1, 2);
                let start_extra = at(&c2_r1, 1);
                let _ = writeln!(
                    extension_log,
                    "XXX\t{r}\t{id}\t{}\t{}\t{}\t{}",
                    perl_format(at(&on1, 2)),
                    perl_format(end_extra),
                    perl_format(at(&on2, 1)),
                    perl_format(start_extra)
                );
                entry.insert(
                    id.clone(),
                    (at(&on1, 2) + end_extra, at(&on2, 1) - start_extra, reverse),
                );
            } else {
                let end_extra = length(c2) - at(&c2_r1, 2);
                let start_extra = at(&c1_r1, 1);
                let _ = writeln!(
                    extension_log,
                    "YYY\t{r}\t{id}\t{}\t{}\t{}\t{}",
                    perl_format(at(&on2, 2)),
                    perl_format(end_extra),
                    perl_format(at(&on1, 1)),
                    perl_format(start_extra)
                );
                entry.insert(
                    id.clone(),
                    (at(&on2, 2) + end_extra, at(&on1, 1) - start_extra, reverse),
                );
            }
        }
    }
    (read_to_edge, reference_reads)
}

/// A scaffold's read pool, gzipped (level 1) to keep the disk use well below
/// the long reads' uncompressed size; minimap2 and racon read it directly.
const POOL: &str = "POOL.fastq.gz";

/// The directory of the read pools, in the gap-filling directory.
const POOL_DIR: &str = "READ_POOLS";

/// The read pools of the scaffolds, written in batches: each batch of a
/// scaffold is one gzip member (compressed in parallel). The members go to
/// a bounded number of bucket files (`READ_POOLS/pool_<n>.fastq.gz`, a
/// scaffold always in the same bucket) instead of one file per scaffold;
/// a scaffold's pool is the concatenation of its members in order, the
/// bytes a pool file of its own would have.
struct Pools {
    dir: String,
    bucket_of: HashMap<String, usize>,
    /// The size of each bucket file.
    ends: Vec<u64>,
    /// Per scaffold, its members: (bucket, offset, length).
    members: HashMap<String, Vec<(usize, u64, u64)>>,
    buffers: HashMap<String, Vec<u8>>,
    size: usize,
    threads: usize,
}

impl Pools {
    const FLUSH_SIZE: usize = 256 << 20;
    const MAX_BUCKETS: usize = 256;

    /// The pools of `scaffolds` (the scaffolds to fill, in order), in
    /// `<ana_dir>/READ_POOLS`.
    fn new(ana_dir: &str, scaffolds: &[String], threads: usize) -> Result<Self> {
        let dir = format!("{ana_dir}/{POOL_DIR}");
        remove_if_exists(&dir)?;
        fs::create_dir_all(&dir).with_context(|| format!("cannot create {dir}"))?;
        let buckets = (threads.max(1) * 4)
            .min(Self::MAX_BUCKETS)
            .min(scaffolds.len())
            .max(1);
        let bucket_of = scaffolds
            .iter()
            .enumerate()
            .map(|(i, name)| (name.clone(), i % buckets))
            .collect();
        Ok(Self {
            dir,
            bucket_of,
            ends: vec![0; buckets],
            members: HashMap::new(),
            buffers: HashMap::new(),
            size: 0,
            threads,
        })
    }

    fn bucket_path(&self, bucket: usize) -> String {
        format!("{}/pool_{bucket}.fastq.gz", self.dir)
    }

    fn add(&mut self, scaff: &str, record: &str) -> Result<()> {
        self.buffers
            .entry(scaff.to_string())
            .or_default()
            .extend_from_slice(record.as_bytes());
        self.size += record.len();
        if self.size > Self::FLUSH_SIZE {
            self.flush()?;
        }
        Ok(())
    }

    /// Appends one member per buffered scaffold to its bucket.
    fn flush(&mut self) -> Result<()> {
        let mut per_bucket: Vec<Vec<(String, Vec<u8>)>> = vec![Vec::new(); self.ends.len()];
        for (scaff, data) in self.buffers.drain() {
            let bucket = *self
                .bucket_of
                .get(&scaff)
                .with_context(|| format!("no read pool for {scaff}"))?;
            per_bucket[bucket].push((scaff, data));
        }
        let jobs: Vec<_> = per_bucket
            .into_iter()
            .enumerate()
            .filter(|(_, list)| !list.is_empty())
            .map(|(bucket, mut list)| {
                list.sort_by(|a, b| a.0.cmp(&b.0));
                let path = self.bucket_path(bucket);
                move || -> Result<(usize, Vec<(String, u64)>)> {
                    let file = OpenOptions::new()
                        .create(true)
                        .append(true)
                        .open(&path)
                        .with_context(|| format!("cannot append to {path}"))?;
                    let mut out = BufWriter::new(file);
                    let mut lengths = Vec::with_capacity(list.len());
                    for (scaff, data) in list {
                        let member = gzip_member(&data)
                            .with_context(|| format!("cannot compress the reads of {scaff}"))?;
                        out.write_all(&member)
                            .with_context(|| format!("cannot append to {path}"))?;
                        lengths.push((scaff, member.len() as u64));
                    }
                    out.flush()
                        .with_context(|| format!("cannot append to {path}"))?;
                    Ok((bucket, lengths))
                }
            })
            .collect();
        let written = run_parallel(jobs, self.threads)
            .map_err(|e| anyhow::anyhow!("cannot write the read pools: {e}"))?;
        for (bucket, lengths) in written {
            for (scaff, length) in lengths {
                let offset = self.ends[bucket];
                self.members
                    .entry(scaff)
                    .or_default()
                    .push((bucket, offset, length));
                self.ends[bucket] += length;
            }
        }
        self.size = 0;
        Ok(())
    }

    /// The size of a scaffold's pool (0: no reads).
    fn pool_size(&self, scaff: &str) -> u64 {
        self.members
            .get(scaff)
            .map_or(0, |m| m.iter().map(|&(_, _, length)| length).sum())
    }

    /// Writes a scaffold's pool to `path`; false (nothing written) if it
    /// has no reads.
    fn write_pool(&self, scaff: &str, path: &str) -> Result<bool> {
        let Some(members) = self.members.get(scaff) else {
            return Ok(false);
        };
        let mut out =
            BufWriter::new(fs::File::create(path).with_context(|| format!("cannot write {path}"))?);
        for &(bucket, offset, length) in members {
            let source = self.bucket_path(bucket);
            let mut file =
                fs::File::open(&source).with_context(|| format!("cannot read {source}"))?;
            file.seek(SeekFrom::Start(offset))
                .with_context(|| format!("cannot read {source}"))?;
            let copied = std::io::copy(&mut file.take(length), &mut out)
                .with_context(|| format!("cannot copy {source} to {path}"))?;
            ensure!(copied == length, "{source}: truncated read pool of {scaff}");
        }
        out.flush()
            .with_context(|| format!("cannot write {path}"))?;
        Ok(true)
    }

    /// Removes the bucket files.
    fn remove(&self) -> Result<()> {
        remove_if_exists(&self.dir)
    }
}

/// `data` as one gzip member (fast compression).
fn gzip_member(data: &[u8]) -> std::io::Result<Vec<u8>> {
    let mut gz = GzEncoder::new(Vec::with_capacity(data.len() / 3), Compression::fast());
    gz.write_all(data)?;
    gz.finish()
}

/// The pre-consensus of the scaffolds to fill, in one FASTA file
/// (`pre_consensus.fa`) with the place of each record.
struct PreConsensus {
    path: String,
    records: HashMap<String, (u64, usize)>,
}

impl PreConsensus {
    /// The record of a scaffold, if it has one.
    fn get(&self, scaff: &str) -> Result<Option<Vec<u8>>> {
        let Some(&(offset, length)) = self.records.get(scaff) else {
            return Ok(None);
        };
        let mut file =
            fs::File::open(&self.path).with_context(|| format!("cannot read {}", self.path))?;
        let mut data = vec![0; length];
        file.seek(SeekFrom::Start(offset))
            .and_then(|_| file.read_exact(&mut data))
            .with_context(|| format!("cannot read {}", self.path))?;
        Ok(Some(data))
    }
}

/// The read pools and the pre-consensus of every scaffold to fill.
fn extract_read_sequence(
    edge_file: &str,
    contig_file: &str,
    opera_lr_dir: &str,
    read_file: &str,
    ana_dir: &str,
    to_fill: &[String],
    threads: usize,
) -> Result<(Pools, PreConsensus)> {
    if !Path::new(contig_file).exists() {
        bail!("CONTIG FILE NOT FOUND {contig_file}");
    }
    if !Path::new(read_file).exists() {
        bail!("READ FILE NOT FOUND {read_file}");
    }
    let is_fastq = {
        let mut first = String::new();
        opera_ms_io::seq::open_maybe_gzip(Path::new(read_file))
            .with_context(|| format!("cannot read {read_file}"))?
            .read_line(&mut first)
            .with_context(|| format!("cannot read {read_file}"))?;
        first.starts_with('@')
    };
    let log_dir = format!("{ana_dir}/LOG");
    remove_if_exists(&log_dir)?;
    fs::create_dir_all(&log_dir).with_context(|| format!("cannot create {log_dir}"))?;

    let scaffold_file = format!("{opera_lr_dir}/scaffolds.scaf");
    let scaf = read(&scaffold_file)?;
    let mut edges = scaffold_edges(&scaf);
    let mut contig_lengths = HashMap::new();
    for line in raw_lines(&read(&format!("{opera_lr_dir}/contigs"))?).skip(1) {
        let mut f = line.split('\t');
        contig_lengths.insert(
            f.next().unwrap_or("").to_string(),
            perl_value(f.next().unwrap_or("")),
        );
    }
    let mut extension_log = String::new();
    let (mut read_to_edge, reference_reads) = read_edges(
        &crate::gz::read(edge_file)?,
        &edges,
        &contig_lengths,
        &mut extension_log,
    );
    write(&format!("{ana_dir}/contig_extention.log"), &extension_log)?;

    let edge_err_path = format!("{ana_dir}/edge_err");
    write(&edge_err_path, "")?;
    let scaff_to_fill: HashSet<String> = to_fill.iter().cloned().collect();

    let mut edge_err = String::new();
    let mut pure_gap: HashMap<String, HashMap<String, f64>> = HashMap::new();
    let mut pools = Pools::new(ana_dir, to_fill, threads)?;
    let mut process = |name: &str,
                       seq: &str,
                       quality: &str,
                       edges: &mut HashMap<String, EdgeInfo>|
     -> Result<()> {
        let Some(list) = read_to_edge.remove(name) else {
            return Ok(());
        };
        let mut ids: Vec<&String> = list.keys().collect();
        ids.sort();
        let mut written = HashSet::new();
        for id in ids {
            let Some(info) = edges.get_mut(id) else {
                continue;
            };
            if !scaff_to_fill.contains(&info.scaff) {
                continue;
            }
            if written.insert(info.scaff.clone()) {
                pools.add(&info.scaff, &format!("@{name}\n{seq}\n+\n{quality}\n"))?;
            }
            if reference_reads.get(id).map(String::as_str) == Some(name) {
                let (start, end, reverse) = list[id];
                if end < start {
                    pure_gap
                        .entry(id.clone())
                        .or_default()
                        .insert(name.to_string(), start - end);
                    let _ = writeln!(
                        edge_err,
                        "{id}\t{name}\tNA\t\t{}\t{}\t{}\tOVERLAP",
                        perl_format(start),
                        perl_format(end),
                        perl_format(start - end)
                    );
                    info.seq = String::new();
                } else {
                    let gap = perl_substr(seq, start.trunc() as i64, (end - start).trunc() as i64)
                        .unwrap_or("");
                    info.seq = if reverse != 0 {
                        reverse_complement(gap)
                    } else {
                        gap.to_string()
                    };
                    let _ = writeln!(
                        edge_err,
                        "{id}\t{name}\t{reverse}\t{}\t{}\t{}\tGAP",
                        perl_format(start),
                        perl_format(end),
                        perl_format(end - start)
                    );
                }
            }
        }
        Ok(())
    };

    let mut reader = opera_ms_io::seq::open_maybe_gzip(Path::new(read_file))
        .with_context(|| format!("cannot read {read_file}"))?;
    let next_line = |reader: &mut Box<dyn BufRead + Send>| -> Result<Option<String>> {
        let mut buf = Vec::new();
        if reader
            .read_until(b'\n', &mut buf)
            .with_context(|| format!("cannot read {read_file}"))?
            == 0
        {
            return Ok(None);
        }
        Ok(Some(String::from_utf8_lossy(&buf).into_owned()))
    };
    let (mut name, mut seq, mut quality) = (String::new(), String::new(), String::new());
    while let Some(raw) = next_line(&mut reader)? {
        let line = chop(&raw);
        if line.starts_with('>') || is_fastq {
            if !seq.is_empty() {
                process(&name, &seq, &quality, &mut edges)?;
            }
            if is_fastq {
                let first = line.split(char::is_whitespace).next().unwrap_or("");
                name = first.get(1..).unwrap_or("").to_string();
                let s = next_line(&mut reader)?.unwrap_or_default();
                seq = s.strip_suffix('\n').unwrap_or(&s).to_string();
                next_line(&mut reader)?;
                let q = next_line(&mut reader)?.unwrap_or_default();
                quality = q.strip_suffix('\n').unwrap_or(&q).to_string();
            } else {
                name = line[1..]
                    .split(char::is_whitespace)
                    .next()
                    .unwrap_or("")
                    .to_string();
                seq.clear();
                quality.clear();
            }
        } else if !line.is_empty() {
            // FASTA sequence line (multi-line records)
            seq.push_str(line);
            quality = "I".repeat(seq.len());
        }
    }
    // the last read
    if !seq.is_empty() {
        process(&name, &seq, &quality, &mut edges)?;
    }
    pools.flush()?;
    write(&edge_err_path, &edge_err)?;

    let pre = write_filled_gap_scaffold(
        contig_file,
        &scaf,
        &edges,
        &pure_gap,
        &reference_reads,
        &scaff_to_fill,
        ana_dir,
    )?;
    Ok((pools, pre))
}

/// `write_filled_gap_scaffold`: `pre_consensus.fa` (the pre-consensus of
/// the scaffolds to fill), `gap_size.dat`.
fn write_filled_gap_scaffold(
    contig_file: &str,
    scaf: &str,
    edges: &HashMap<String, EdgeInfo>,
    pure_gap: &HashMap<String, HashMap<String, f64>>,
    reference_reads: &HashMap<String, String>,
    scaff_to_fill: &HashSet<String>,
    ana_dir: &str,
) -> Result<PreConsensus> {
    // contig name: the [\w|/.] characters after '>'
    let mut contigs: HashMap<String, String> = HashMap::new();
    let (mut name, mut seq) = (String::new(), String::new());
    for line in raw_lines(&read(contig_file)?) {
        let line = chop(line);
        if let Some(rest) = line.strip_prefix('>') {
            if !name.is_empty() {
                contigs.insert(std::mem::take(&mut name), std::mem::take(&mut seq));
            }
            name = rest
                .chars()
                .take_while(|c| c.is_ascii_alphanumeric() || matches!(c, '_' | '|' | '/' | '.'))
                .collect();
            seq.clear();
        } else {
            seq.push_str(line);
        }
    }
    contigs.insert(name, seq);
    let contig_seq = |c: &str, o: &str| {
        let s = contigs.get(c).map_or("", String::as_str);
        if o == "-" {
            reverse_complement(s)
        } else {
            s.to_string()
        }
    };

    let mut gap_size = String::new();
    let (mut scaff, mut scaff_seq) = (String::new(), String::new());
    let (mut c1, mut o1, mut c2, mut o2) =
        (String::new(), String::new(), String::new(), String::new());
    let mut gap_length = String::new();
    let mut trimming = 0.0;
    let mut lines = raw_lines(scaf);
    let pre_path = format!("{ana_dir}/pre_consensus.fa");
    let mut pre_file = BufWriter::new(
        fs::File::create(&pre_path).with_context(|| format!("cannot write {pre_path}"))?,
    );
    let mut pre = PreConsensus {
        path: pre_path.clone(),
        records: HashMap::new(),
    };
    let mut pre_end = 0u64;
    let mut write_pre = |scaff: &str, scaff_seq: &str, c2: &str, o2: &str| -> Result<()> {
        if scaff_to_fill.contains(scaff) && !c2.is_empty() {
            let record = format!(">{scaff}\n{scaff_seq}{}\n", contig_seq(c2, o2));
            pre_file
                .write_all(record.as_bytes())
                .with_context(|| format!("cannot write {pre_path}"))?;
            pre.records
                .insert(scaff.to_string(), (pre_end, record.len()));
            pre_end += record.len() as u64;
        }
        Ok(())
    };
    while let Some(line) = lines.next() {
        let fields: Vec<&str> = chop(line).split('\t').collect();
        let first = fields.first().copied().unwrap_or("");
        let orientation = |f: &[&str]| if f.get(1) == Some(&"EB") { "-" } else { "+" }.to_string();
        if let Some(i) = first.find('>').filter(|i| first.len() > i + 1) {
            write_pre(&scaff, &scaff_seq, &c2, &o2)?;
            scaff_seq.clear();
            scaff = first[i + 1..].to_string();
            trimming = 0.0;
            let next = lines.next().unwrap_or("");
            let f: Vec<&str> = chop(next).split('\t').collect();
            c1 = f.first().copied().unwrap_or("").to_string();
            o1 = orientation(&f);
            gap_length = f.get(3).copied().unwrap_or("").to_string();
            c2.clear();
            continue;
        }
        if !scaff_to_fill.contains(&scaff) {
            continue;
        }
        c2 = first.to_string();
        o2 = orientation(&fields);
        let mut c1_seq = contig_seq(&c1, &o1);
        if trimming != 0.0 && trimming < c1_seq.len() as f64 {
            c1_seq = perl_substr(&c1_seq, trimming.trunc() as i64, i64::MAX / 2)
                .unwrap_or("")
                .to_string();
        }
        scaff_seq.push_str(&c1_seq);
        let (id, _) = edge_id(&c1, "", &c2, "");
        // extract_gap_seq
        let stored = edges.get(&id).map(|e| e.seq.clone()).unwrap_or_default();
        let gap_seq = if !stored.is_empty() {
            stored
        } else {
            reference_reads
                .get(&id)
                .and_then(|r| pure_gap.get(&id).and_then(|m| m.get(r)))
                .map(|v| perl_format(*v))
                .unwrap_or_default()
        };
        let stats;
        if !gap_seq.is_empty() {
            if gap_seq
                .chars()
                .all(|c| c.is_ascii_digit() || matches!(c, ',' | '.' | 'E'))
            {
                trimming = perl_value(&gap_seq);
                stats = perl_format(-trimming);
            } else {
                stats = gap_seq.len().to_string();
                trimming = 0.0;
                scaff_seq.push_str(&gap_seq);
            }
        } else {
            stats = "NA".to_string();
            trimming = 0.0;
            let length = perl_value(&gap_length);
            let n = if length > 0.0 {
                length.trunc() as usize
            } else {
                3
            };
            scaff_seq.push_str(&"N".repeat(n));
        }
        let _ = writeln!(gap_size, "{id}\t{gap_length}\t{stats}\t");
        c1 = c2.clone();
        o1 = o2.clone();
        gap_length = fields.get(3).copied().unwrap_or("").to_string();
    }
    write_pre(&scaff, &scaff_seq, &c2, &o2)?;
    pre_file
        .flush()
        .with_context(|| format!("cannot write {}", pre.path))?;
    drop(pre_file);
    write(&format!("{ana_dir}/gap_size.dat"), &gap_size)?;
    Ok(pre)
}

// ---- tiling alignment -----------------------------------------------------

/// The tiling alignments with minimap2: every contig is mapped on the
/// consensus in one `minimap2 -cx asm20 -N 50` call (secondary alignments
/// kept: a repeated contig keeps all its placements), written as a
/// `show-coords -T -l` table (4 header lines, then `S1 E1 S2 E2 LEN1 LEN2
/// %IDY LENR LENQ TAGR TAGQ`) for `get_paf_file`. Filters: identity >= 95 %,
/// and per consensus sequence the alignments are taken by decreasing
/// identity, then number of matching bases (ranking by matches alone lets
/// longer contigs of a dominant strain take the regions of a minor one),
/// dropping one only when a single kept alignment covers at least 90 % of
/// its consensus interval: a contig inside a region already tiled by a
/// better one goes; partial overlaps of neighbouring contigs stay.
fn map_contigs_minimap2(
    rec: &StepRecorder<'_>,
    consensus: &str,
    contig_file: &str,
    dir_res: &str,
    final_res_file: &str,
    nproc: usize,
    minimap2_dir: &str,
) -> Result<()> {
    fs::create_dir_all(dir_res).with_context(|| format!("cannot create {dir_res}"))?;
    let paf_path = format!("{final_res_file}.minimap2.paf");
    let err_path = format!("{final_res_file}.minimap2.err");
    let report = Pipeline::new(
        Cmd::new(format!("{minimap2_dir}/minimap2"))
            .args(["-cx", "asm20", "-N", "50", "-t"])
            .arg(nproc.max(1).to_string())
            .args([consensus, contig_file]),
    )
    .stdout(Output::File(paf_path.clone().into()))
    .stderr(Output::File(err_path.clone().into()))
    .run()
    .with_context(|| format!("Error during the tiling alignment. Please see {err_path}"))?;
    rec.record(&report);
    let paf = read(&paf_path)?;
    let coords = coords_from_paf(&paf, consensus, contig_file);
    write(final_res_file, &coords)?;
    fs::remove_file(&paf_path).with_context(|| format!("cannot remove {paf_path}"))
}

/// One alignment of a contig on the consensus.
struct TileHit<'a> {
    s1: u64,
    e1: u64,
    s2: u64,
    e2: u64,
    identity: f64,
    matches: u64,
    tlen: &'a str,
    qlen: &'a str,
    target: &'a str,
    query: &'a str,
}

/// minimap2 PAF lines to the `show-coords -T -l` table (see
/// [`map_contigs_minimap2`]).
fn coords_from_paf(paf: &str, ref_file: &str, query_file: &str) -> String {
    let mut by_target: BTreeMap<&str, Vec<TileHit<'_>>> = BTreeMap::new();
    for line in paf.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 {
            continue;
        }
        let num = |i: usize| f[i].parse::<u64>().unwrap_or(0);
        let (qs, qe, ts, te, matches, block) = (num(2), num(3), num(7), num(8), num(9), num(10));
        if block == 0 {
            continue;
        }
        let identity = 100.0 * matches as f64 / block as f64;
        if identity < 95.0 {
            continue;
        }
        // show-coords: 1-based, reverse alignments with S2 > E2
        let (s2, e2) = if f[4] == "-" {
            (qe, qs + 1)
        } else {
            (qs + 1, qe)
        };
        by_target.entry(f[5]).or_default().push(TileHit {
            s1: ts + 1,
            e1: te,
            s2,
            e2,
            identity,
            matches,
            tlen: f[6],
            qlen: f[1],
            target: f[5],
            query: f[0],
        });
    }
    let mut out = format!(
        "{ref_file} {query_file}\nNUCMER\n\n[S1]\t[E1]\t[S2]\t[E2]\t[LEN 1]\t[LEN 2]\t[% IDY]\t[LEN R]\t[LEN Q]\t[TAGS]\n"
    );
    for (_, mut hits) in by_target {
        // best first: highest identity, then most matches (identity first,
        // so that a contig of the right strain wins over a longer one of a
        // related strain); ties in a fixed order (deterministic output)
        hits.sort_by(|a, b| {
            b.identity
                .total_cmp(&a.identity)
                .then(b.matches.cmp(&a.matches))
                .then(a.s1.cmp(&b.s1))
                .then(a.query.cmp(b.query))
        });
        let mut kept: Vec<TileHit<'_>> = Vec::new();
        for hit in hits {
            let span = hit.e1 + 1 - hit.s1;
            let contained = kept.iter().any(|k| {
                let overlap = (hit.e1.min(k.e1) + 1).saturating_sub(hit.s1.max(k.s1));
                overlap * 10 >= span * 9
            });
            if !contained {
                kept.push(hit);
            }
        }
        kept.sort_by(|a, b| a.s1.cmp(&b.s1).then(a.query.cmp(b.query)));
        for h in kept {
            let _ = writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{}",
                h.s1,
                h.e1,
                h.s2,
                h.e2,
                h.e1 + 1 - h.s1,
                h.s2.max(h.e2) + 1 - h.s2.min(h.e2),
                h.identity,
                h.tlen,
                h.qlen,
                h.target,
                h.query
            );
        }
    }
    out
}

// ---- get_paf_file / get_contig_in_gap ----------------------------------------

/// GNU `sort -k10,10 -k1,1n` in the C locale, then `grep -v '\['`.
fn sort_coords(text: &str) -> Vec<&str> {
    // a GNU field: the blanks before it and the non-blanks
    let field = |line: &str, n: usize| -> String {
        let bytes = line.as_bytes();
        let mut i = 0;
        let mut count = 0;
        while i < bytes.len() {
            let start = i;
            while i < bytes.len() && (bytes[i] == b' ' || bytes[i] == b'\t') {
                i += 1;
            }
            while i < bytes.len() && bytes[i] != b' ' && bytes[i] != b'\t' {
                i += 1;
            }
            count += 1;
            if count == n {
                return line[start..i].to_string();
            }
        }
        String::new()
    };
    let numeric = |line: &str| -> f64 {
        let t = line.trim_start_matches([' ', '\t']);
        let mut end = 0;
        for (i, c) in t.char_indices() {
            if c.is_ascii_digit() || (i == 0 && c == '-') || c == '.' {
                end = i + 1;
            } else {
                break;
            }
        }
        t[..end].parse().unwrap_or(0.0)
    };
    let mut lines: Vec<(String, f64, &str)> = text
        .lines()
        .map(|l| (field(l, 10), numeric(l), l))
        .collect();
    lines.sort_by(|a, b| {
        a.0.as_bytes()
            .cmp(b.0.as_bytes())
            .then(a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
            .then(a.2.as_bytes().cmp(b.2.as_bytes()))
    });
    lines
        .into_iter()
        .map(|(_, _, l)| l)
        .filter(|l| !l.contains('['))
        .collect()
}

/// `get_paf_file`: the contigs tiled on the consensus, as PAF lines;
/// `S1_gap_size.dat`.
fn get_paf_file(scaf: &str, coords: &str, only_scaffold_contigs: bool) -> (String, String) {
    let mut keep: HashMap<String, HashSet<String>> = HashMap::new();
    let mut scaff = String::new();
    for line in raw_lines(scaf) {
        let fields = split_ws(chop(line));
        let name = fields.first().copied().unwrap_or("");
        match name.find('>').filter(|i| name.len() > i + 1) {
            Some(i) => {
                scaff = name[i + 1..].to_string();
                keep.insert(scaff.clone(), HashSet::new());
            }
            None => {
                keep.entry(scaff.clone())
                    .or_default()
                    .insert(name.to_string());
            }
        }
    }
    let kept = |s: &str, c: &str| keep.get(s).is_some_and(|set| set.contains(c));
    let (mut paf, mut gaps) = (String::new(), String::new());
    let mut previous = String::new();
    let mut prev_scaffold_end = -1.0;
    let (mut scaffold_name, mut scaffold_length) = (String::new(), String::new());
    let (mut identity, mut contig_len) = (0.0, String::new());
    let (mut contig_start, mut contig_end, mut scaffold_start, mut scaffold_end) =
        (0.0, 0.0, 0.0, 0.0);
    let mut orientation = "";
    let mut seen = false;
    for line in sort_coords(coords).into_iter().skip(3) {
        let f: Vec<&str> = line.split('\t').collect();
        let get = |i: usize| f.get(i).copied().unwrap_or("");
        let name = get(10);
        if name != previous {
            if !previous.is_empty() {
                let mapping_length = contig_end - contig_start;
                let fraction = mapping_length / perl_value(&contig_len);
                let in_scaffold = kept(&scaffold_name, &previous);
                if (identity >= 0.95 || in_scaffold)
                    && (fraction >= 0.90
                        || (only_scaffold_contigs && in_scaffold && fraction > 0.55))
                {
                    let _ = writeln!(
                        paf,
                        "{previous}\t{contig_len}\t{}\t{}\t{orientation}\t{scaffold_name}\t{scaffold_length}\t{}\t{}\t{}\t{}\t60\ttp:A:P\tcm:i:175\ts1:i:915\ts2:i:0",
                        perl_format(contig_start),
                        perl_format(contig_end),
                        perl_format(scaffold_start),
                        perl_format(scaffold_end),
                        perl_format(((scaffold_end - scaffold_start) * identity).trunc()),
                        perl_format(scaffold_end - scaffold_start)
                    );
                    if prev_scaffold_end > 0.0 && prev_scaffold_end < scaffold_start {
                        let _ = writeln!(
                            gaps,
                            "{}\t{previous}\t{name}",
                            perl_format(scaffold_start - prev_scaffold_end)
                        );
                    }
                }
                prev_scaffold_end = scaffold_end;
            }
            previous = name.to_string();
            seen = true;
            scaffold_name = get(9).to_string();
            scaffold_length = get(7).to_string();
            identity = perl_value(get(6)) / 100.0;
            contig_len = get(8).to_string();
            contig_start = perl_value(get(2));
            contig_end = perl_value(get(3));
            orientation = "+";
            if contig_end < contig_start {
                orientation = "-";
                std::mem::swap(&mut contig_start, &mut contig_end);
            }
            contig_start -= 1.0;
            scaffold_start = perl_value(get(0)) - 1.0;
            scaffold_end = perl_value(get(1));
        } else if orientation == "+" {
            if contig_end < perl_value(get(3)) {
                contig_end = perl_value(get(3));
                scaffold_end = perl_value(get(1));
            }
        } else if contig_start > perl_value(get(3)) - 1.0 {
            contig_start = perl_value(get(3)) - 1.0;
            scaffold_end = perl_value(get(1));
        }
    }
    // the last contig, always
    if seen {
        let _ = writeln!(
            paf,
            "{previous}\t{contig_len}\t{}\t{}\t{orientation}\t{scaffold_name}\t{scaffold_length}\t{}\t{}\t{}\t{}\t60\ttp:A:P\tcm:i:0\ts1:i:0\ts2:i:0",
            perl_format(contig_start),
            perl_format(contig_end),
            perl_format(scaffold_start),
            perl_format(scaffold_end),
            perl_format(((scaffold_end - scaffold_start) * 0.99).trunc()),
            perl_format(scaffold_end - scaffold_start)
        );
    } else {
        paf.push_str("\t\t\t\t\t\t\t\t\t0\t0\t60\ttp:A:P\tcm:i:0\ts1:i:0\ts2:i:0\n");
    }
    (paf, gaps)
}

/// `get_contig_in_gap`: the contigs of the PAF file, from the contig file,
/// each once with its whole sequence on one line (what `remap` reads).
fn get_contig_in_gap(paf: &str, contig_file: &str) -> Result<String> {
    let mapped: HashSet<&str> = raw_lines(paf)
        .map(|line| line.split('\t').next().unwrap_or(""))
        .collect();
    let mut out = String::new();
    let mut print = false;
    for line in raw_lines(&read(contig_file)?) {
        // not `chop`: the last line may have no newline
        let line = line.strip_suffix('\n').unwrap_or(line);
        if let Some(header) = line.strip_prefix('>') {
            let id = split_ws(header).first().copied().unwrap_or("");
            if print {
                out.push('\n');
            }
            print = mapped.contains(id);
            if print {
                out.push('>');
                out.push_str(id);
                out.push('\n');
            }
        } else if print {
            out.push_str(line);
        }
    }
    if print {
        out.push('\n');
    }
    Ok(out)
}

// ---- remap ----------------------------------------------------------------

struct Placement {
    contig_start: i64,
    contig_end: i64,
    scaffold_start: i64,
    scaffold_end: i64,
}

fn py_int(text: &str) -> Result<i64> {
    text.trim()
        .parse()
        .with_context(|| format!("Remap: invalid literal for int(): {text:?}"))
}

/// The contig sequences put back in place of the consensus: writes `<assembly
/// base>_remapped.fasta` and `_scaff_to_rescue.dat`; returns its stdout.
fn remap(paf: &str, assembly_path: &str, contig_path: &str) -> Result<(String, String, String)> {
    let mut stdout = String::from("Reading mapping file\n\n");
    // valid_contig_map: scaffold -> placements (insertion order)
    let mut valid: HashMap<String, Vec<(String, Placement)>> = HashMap::new();
    let mut mapping_ids: HashMap<String, Vec<String>> = HashMap::new();
    for (i, line) in raw_lines(paf).enumerate() {
        let f: Vec<&str> = line.split_whitespace().collect();
        let get = |i: usize| f.get(i).copied().unwrap_or("");
        let name = get(0);
        let (len, start, end) = (py_int(get(1))?, py_int(get(2))?, py_int(get(3))?);
        if len == 0 {
            bail!("Remap: division by zero (contig {name} of length 0)");
        }
        if (start - end).abs() as f64 / len as f64 <= 0.0 {
            continue;
        }
        // every placement is kept: a repeated contig fills each of its
        // copies
        let id = format!("{name}_XXX_{}", i + 1);
        mapping_ids
            .entry(name.to_string())
            .or_default()
            .push(format!("{}{id}", get(4)));
        valid.entry(get(5).to_string()).or_default().push((
            id,
            Placement {
                contig_start: start,
                contig_end: end,
                scaffold_start: py_int(get(7))?,
                scaffold_end: py_int(get(8))?,
            },
        ));
    }

    stdout.push_str("Reading contig sequence files\n\n");
    let mut contig_seq: HashMap<String, String> = HashMap::new();
    let mut extract: Option<String> = None;
    for line in raw_lines(&read(contig_path)?) {
        if let Some(name) = extract.take() {
            let seq = line.trim_end();
            let mut reverse: Option<String> = None;
            for id in &mapping_ids[&name] {
                let (ori, rest) = id.split_at(1);
                let s = if ori == "-" {
                    // IUPAC complement
                    if reverse.is_none() {
                        reverse = Some(reverse_complement(seq));
                    }
                    reverse.clone().unwrap_or_default()
                } else {
                    seq.to_string()
                };
                contig_seq.insert(rest.to_string(), s);
            }
        } else {
            let key = line.get(1..).unwrap_or("").trim_end();
            if mapping_ids.contains_key(key) {
                extract = Some(key.to_string());
            }
        }
    }

    stdout.push_str("Processing scaffold file\n\n");
    let (mut remapped, mut rescue) = (String::new(), String::new());
    let mut process = |scaff: &str, assembly: &str| -> Result<()> {
        let Some(placements) = valid.get(scaff) else {
            rescue.push_str(scaff);
            rescue.push('\n');
            return Ok(());
        };
        // start sites: one contig per start (the last one), sorted
        let mut starts: Vec<(i64, usize)> = Vec::new();
        for (i, (_, p)) in placements.iter().enumerate() {
            match starts.iter_mut().find(|(s, _)| *s == p.scaffold_start) {
                Some(entry) => entry.1 = i,
                None => starts.push((p.scaffold_start, i)),
            }
        }
        starts.sort();
        let initial = assembly.as_bytes();
        let mut out: Vec<u8> = Vec::new();
        let mut pos: i64 = 0;
        let index_error = || anyhow::anyhow!("Remap: IndexError in scaffold {scaff}");
        for (start, i) in starts {
            let (id, p) = &placements[i];
            for k in pos..start {
                out.push(
                    *initial
                        .get(usize::try_from(k).map_err(|_| index_error())?)
                        .ok_or_else(index_error)?,
                );
            }
            pos = pos.max(start);
            let cut = pos - start + p.contig_start;
            let insert = contig_seq
                .get(id)
                .ok_or_else(|| anyhow::anyhow!("Remap: KeyError {id}"))?
                .as_bytes();
            for k in cut..p.contig_end {
                out.push(
                    *insert
                        .get(usize::try_from(k).map_err(|_| index_error())?)
                        .ok_or_else(index_error)?,
                );
            }
            pos = pos.max(p.scaffold_end);
        }
        if let Ok(p) = usize::try_from(pos) {
            out.extend_from_slice(initial.get(p..).unwrap_or(&[]));
        }
        let _ = writeln!(remapped, ">{scaff}\n{}", String::from_utf8_lossy(&out));
        Ok(())
    };
    let (mut scaff, mut assembly) = (String::new(), String::new());
    for line in raw_lines(&read(assembly_path)?) {
        if let Some(rest) = line.strip_prefix('>') {
            if !scaff.is_empty() {
                let _ = write!(stdout, " *** Process {scaff}\n\n");
                process(&scaff, &assembly)?;
            }
            assembly.clear();
            scaff = rest.trim_end().to_string();
        } else {
            assembly.push_str(line.trim_end());
        }
    }
    let _ = write!(stdout, " *** Process last {scaff}\n\n");
    process(&scaff, &assembly)?;
    Ok((remapped, rescue, stdout))
}

/// `".".join(path.split(".")[:-1])`.
fn python_stem(path: &str) -> &str {
    path.rsplit_once('.').map_or("", |(stem, _)| stem)
}

// ---- gap filling ------------------------------------------------------------

/// `merge_contig_from_scaf`: the scaffolds in their order, filled or not.
fn merge_contig_from_scaf(filled: &str, scaffold_seq: &str, scaffolds: &Scaffolds) -> String {
    let mut seqs: HashMap<String, Option<String>> = HashMap::new();
    let mut lines = raw_lines(scaffold_seq);
    while let Some(line) = lines.next() {
        let Some(rest) = line
            .strip_prefix('>')
            .map(|r| r.strip_suffix('\n').unwrap_or(r))
            .filter(|r| !r.is_empty())
        else {
            continue;
        };
        let name = split_ws(rest).first().copied().unwrap_or("").to_string();
        let next = lines.next();
        if !scaffolds.multi.contains(&name) {
            seqs.insert(name, next.map(str::to_string));
        }
    }
    let mut lines = raw_lines(filled);
    while let Some(line) = lines.next() {
        let Some(rest) = line
            .strip_prefix('>')
            .map(|r| r.strip_suffix('\n').unwrap_or(r))
            .filter(|r| !r.is_empty())
        else {
            continue;
        };
        seqs.insert(rest.to_string(), lines.next().map(str::to_string));
    }
    let mut out = String::new();
    for name in &scaffolds.order {
        let seq = seqs.get(name).cloned().flatten();
        let length = seq.as_ref().map_or(-1, |s| s.len() as i64 - 1);
        let _ = write!(
            out,
            ">{name} length: {length} cov: 0.0\n{}",
            seq.unwrap_or_default()
        );
    }
    out
}

/// The read pools' upper size, from the long-read file: a read is pooled
/// at most once per scaffold, and in practice once (it bridges contigs of
/// one scaffold); gzip level 1 compresses long-read FASTQ to about 60 %
/// (qualities compress badly), a gzipped input is taken at its size plus
/// 10 % (level 1 compresses a little less than the usual level 6).
fn pool_estimate(read_file: &str) -> Option<u64> {
    let size = fs::metadata(read_file).ok()?.len();
    let gzipped = opera_ms_io::seq::is_gzip(Path::new(read_file)).ok()?;
    Some(if gzipped {
        size / 10 * 11
    } else {
        size / 10 * 6
    })
}

/// A warning when the gap-filling directory's file system has less free
/// space than the read pools may need (`None`: enough, or unknown).
fn space_warning(wd: &str, read_file: &str) -> Option<String> {
    let needed = pool_estimate(read_file)?;
    let free = opera_ms_sys::fs::free_space(Path::new(wd)).ok()?;
    (free < needed).then(|| {
        format!(
            " *** warning: {} GB free in {wd}; gap filling's read pools may need up to {} GB",
            free / 1_000_000_000,
            needed.div_ceil(1_000_000_000)
        )
    })
}

/// The racon jobs: the scaffolds by decreasing work (the size of
/// their read pool), largest first, and each with
/// ceil(work / (total work / threads)) threads, at most `nproc`: a scaffold
/// bigger than a thread's fair share gets several threads instead of
/// running alone at the end. minimap2 and racon results don't depend on
/// their thread count.
fn racon_plan(work: Vec<(String, u64)>, nproc: usize) -> Vec<(String, usize)> {
    let mut jobs = work;
    jobs.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    let total: u64 = jobs.iter().map(|j| j.1).sum();
    let share = (total / nproc as u64).max(1);
    jobs.into_iter()
        .map(|(name, w)| {
            let threads = usize::try_from(w.div_ceil(share))
                .unwrap_or(nproc)
                .clamp(1, nproc);
            (name, threads)
        })
        .collect()
}

/// The racon job of a scaffold: its pre-consensus and read pool are put in
/// `<wd>/<scaff>/` (`scaffold_racon`'s inputs), and the directory removed
/// once `racon.fa` is read; it is left for inspection on errors.
#[allow(clippy::too_many_arguments)]
fn racon_job(
    rec: &StepRecorder<'_>,
    wd: &str,
    scaff: &str,
    pools: &Pools,
    pre: &PreConsensus,
    a: &Args<'_>,
    threads: usize,
) -> Result<(String, String)> {
    let dir = format!("{wd}/{scaff}");
    remove_if_exists(&dir)?;
    fs::create_dir_all(&dir).with_context(|| format!("cannot create {dir}"))?;
    if let Some(record) = pre.get(scaff)? {
        let path = format!("{dir}/pre_consensus.fa");
        fs::write(&path, record).with_context(|| format!("cannot write {path}"))?;
    }
    pools.write_pool(scaff, &format!("{dir}/{POOL}"))?;
    scaffold_racon(rec, &dir, a.racon_dir, a.minimap2_dir, threads)?;
    let racon = read(&format!("{dir}/racon.fa"))?;
    remove_if_exists(&dir)?;
    Ok((scaff.to_string(), racon))
}

/// `racon.fa` of a scaffold.
fn scaffold_racon(
    rec: &StepRecorder<'_>,
    dir: &str,
    racon_dir: &str,
    minimap2_dir: &str,
    threads: usize,
) -> Result<()> {
    let threads = threads.max(1).to_string();
    let pre = format!("{dir}/pre_consensus.fa");
    let pool = format!("{dir}/{POOL}");
    let racon_fa = format!("{dir}/racon.fa");
    if !Path::new(&pool).exists() {
        fs::copy(&pre, &racon_fa).with_context(|| format!("cannot copy {pre}"))?;
        return Ok(());
    }
    let paf = format!("{pre}.filtered.paf");
    let report = Pipeline::new(
        Cmd::new(format!("{minimap2_dir}/minimap2"))
            .args(["--cs=short", "-w5", "-m0", "-t", &threads])
            .args([&pre, &pool]),
    )
    .stdout(Output::File(paf.clone().into()))
    .stderr(Output::File(format!("{paf}.err").into()))
    .run()
    .with_context(|| format!("Error during minimap2. Please see {paf}.err for details"))?;
    rec.record(&report);
    // racon 1.x with the parameters of the racon 0.x OPERA-MS 0.9 used
    // (window 500, match 5, mismatch -4) and a single linear gap penalty
    // (-8, its gap opening); `-u` keeps the targets without coverage; the
    // header tags racon 1.x adds (LN, RC, XC) are dropped. Threads: see
    // `racon_threads`.
    let raw = format!("{racon_fa}.raw");
    let args = [
        "-u", "-t", &threads, "-w", "500", "-m", "5", "-x", "-4", "-g", "-8", &pool, &paf, &pre,
    ];
    let report = Pipeline::new(Cmd::new(format!("{racon_dir}/racon")).args(args))
        .stdout(Output::File(raw.clone().into()))
        .stderr(Output::File(format!("{racon_fa}.err").into()))
        .run()
        .with_context(|| format!("Error during racon. Please see {racon_fa}.err for details"))?;
    rec.record(&report);
    let text = read(&raw)?;
    let mut out = String::with_capacity(text.len());
    for line in text.split_inclusive('\n') {
        match line.strip_prefix('>') {
            Some(header) => {
                let name = header.split_whitespace().next().unwrap_or("");
                out.push('>');
                out.push_str(name);
                out.push('\n');
            }
            None => out.push_str(line),
        }
    }
    write(&racon_fa, &out)?;
    fs::remove_file(&raw).with_context(|| format!("cannot remove {raw}"))?;
    if fs::metadata(&racon_fa).is_ok_and(|m| m.len() == 0) {
        fs::copy(&pre, &racon_fa).with_context(|| format!("cannot copy {pre}"))?;
    }
    Ok(())
}

/// One tiling: nucmer, PAF, contigs in the gaps, remap, rescue of the
/// scaffolds without tiling. Returns the remapped consensus path.
#[allow(clippy::too_many_arguments)]
fn tiling(
    rec: &StepRecorder<'_>,
    a: &Args<'_>,
    scaf: &str,
    consensus: &str,
    dir_res: &str,
    name: &str,
    only_scaffold_contigs: bool,
    pre: &PreConsensus,
    log: &Mutex<String>,
) -> Result<String> {
    let tilling_dir = format!("{}/TILLING", a.working_dir);
    let base = format!("{tilling_dir}/{name}");
    map_contigs_minimap2(
        rec,
        consensus,
        a.contig_file,
        dir_res,
        &format!("{base}.coords"),
        a.nproc,
        a.minimap2_dir,
    )
    .with_context(|| "Error during tilling generation")?;
    let (paf, gaps) = get_paf_file(
        scaf,
        &read(&format!("{base}.coords"))?,
        only_scaffold_contigs,
    );
    write(&format!("{base}.paf"), &paf)?;
    write(&format!("{}/S1_gap_size.dat", a.working_dir), &gaps)?;
    let contigs = get_contig_in_gap(&paf, a.contig_file)?;
    let contig_fa = format!("{base}_contig.fa");
    write(&contig_fa, &contigs)?;
    let (remapped, rescue, stdout) =
        remap(&paf, consensus, &contig_fa).with_context(|| "Error during contig remapping")?;
    let stem = python_stem(consensus);
    let remapped_path = format!("{stem}_remapped.fasta");
    write(&remapped_path, &remapped)?;
    write(&format!("{stem}_scaff_to_rescue.dat"), &rescue)?;
    let out_name = if only_scaffold_contigs {
        "tilling_1_remap.out"
    } else {
        "tilling_2_remap.out"
    };
    write(&format!("{tilling_dir}/{out_name}"), &stdout)?;
    // rescue the pre-consensus of the scaffolds without tiling
    for line in raw_lines(&rescue) {
        let scaff = chop(line);
        match pre.get(scaff)? {
            Some(data) => append(&remapped_path, &data)?,
            None => {
                if let Ok(mut log) = log.lock() {
                    let _ = writeln!(log, "no pre-consensus of {scaff} to rescue");
                }
            }
        }
    }
    Ok(remapped_path)
}

/// The gap filling of the scaffolds.
pub fn gap_filling(rec: &StepRecorder<'_>, a: &Args<'_>) -> Result<()> {
    let wd = a.working_dir;
    fs::create_dir_all(wd).with_context(|| format!("cannot create {wd}"))?;
    let log = Mutex::new(String::new());
    let result = gap_filling_steps(rec, a, &log);
    let text = log.into_inner().unwrap_or_default();
    append(&format!("{wd}/gap_filling.err"), text.as_bytes())?;
    result
}

fn gap_filling_steps(rec: &StepRecorder<'_>, a: &Args<'_>, log: &Mutex<String>) -> Result<()> {
    let wd = a.working_dir;
    let time_path = format!("{wd}/../../time.log");
    let opera_lr_dir = format!("{wd}/../");
    let mut start = Instant::now();
    let scaf = read(&format!("{opera_lr_dir}/scaffolds.scaf"))?;
    let scaffolds = split_scaffold_file(&scaf, perl_value(a.scaffold_length_threshold));
    let to_fill = &scaffolds.to_fill;
    time_log(&time_path, "g_dir_construction", start)?;
    if let Some(warning) = space_warning(wd, a.read_file) {
        opera_ms_core::note!("{warning}");
    }

    start = Instant::now();
    let (pools, pre) = rec.substep("extract_read_sequence", |_| {
        extract_read_sequence(
            &format!("{}/edge_read_info.dat", a.read_mapping_dir),
            a.contig_file,
            &opera_lr_dir,
            a.read_file,
            wd,
            to_fill,
            a.nproc.max(1),
        )
        .context("Error in during read sequence extraction")
    })?;
    time_log(&time_path, "g_extract_read", start)?;
    let bytes: u64 = to_fill.iter().map(|n| pools.pool_size(n)).sum();
    opera_ms_core::log::file_only(&format!(
        " *** gap filling: {} read pools, {} MB (gzipped)",
        pools.members.len(),
        bytes / 1_000_000
    ));

    // consensus of each scaffold
    start = Instant::now();
    let work = to_fill
        .iter()
        .map(|n| (n.clone(), pools.pool_size(n)))
        .collect();
    let plan = racon_plan(work, a.nproc.max(1));
    let racon: HashMap<String, String> = rec.substep("racon", |sub| {
        let (pools, pre) = (&pools, &pre);
        let jobs: Vec<_> = plan
            .iter()
            .map(|(name, threads)| move || racon_job(sub, wd, name, pools, pre, a, *threads))
            .collect();
        run_parallel(jobs, a.nproc.max(1))
            .map(|results| results.into_iter().collect())
            .map_err(|e| anyhow::anyhow!("Error during consensus sequence generation: {e}"))
    })?;
    pools.remove()?;
    let tilling_dir = format!("{wd}/TILLING");
    fs::create_dir_all(&tilling_dir).with_context(|| format!("cannot create {tilling_dir}"))?;
    let consensus = format!("{tilling_dir}/consensus.fa");
    let mut text = String::new();
    let mut names: Vec<&String> = racon
        .keys()
        .filter(|n| !n.starts_with('.') && n.contains("opera_scaffold_"))
        .collect();
    names.sort();
    for n in names {
        for line in raw_lines(&racon[n]) {
            text.push_str(&line.replacen("Consensus_", "", 1));
        }
    }
    drop(racon);
    write(&consensus, &text)?;
    drop(text);
    time_log(&time_path, "g_racon", start)?;

    start = Instant::now();
    let first = rec.substep("tiling_1", |sub| {
        tiling(
            sub,
            a,
            &scaf,
            &consensus,
            wd,
            "first_tilling",
            true,
            &pre,
            log,
        )
    })?;
    time_log(&time_path, "g_tilling_1", start)?;
    start = Instant::now();
    let second = rec.substep("tiling_2", |sub| {
        tiling(
            sub,
            a,
            &scaf,
            &first,
            &format!("{tilling_dir}/"),
            "second_tilling",
            false,
            &pre,
            log,
        )
    })?;
    time_log(&time_path, "g_tilling_2", start)?;
    remove_if_exists(&pre.path)?;

    let filled = merge_contig_from_scaf(
        &read(&second)?,
        &read(&format!("{opera_lr_dir}/scaffoldSeq.fasta"))?,
        &scaffolds,
    );
    write(&format!("{wd}/../scaffoldSeq.fasta.filled"), &filled)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn substr_like_perl() {
        let cases = [
            (0, 2, Some("ab")),
            (1, 5, Some("bc")),
            (3, 1, Some("")),
            (4, 1, None),
            (-2, 1, Some("b")),
            (-5, 3, Some("a")),
            (-5, 1, None),
            (-4, 2, Some("a")),
            (-4, 1, Some("")),
            (-6, 3, Some("")),
            (-3, 5, Some("abc")),
        ];
        for (o, l, expected) in cases {
            assert_eq!(perl_substr("abc", o, l), expected, "({o}, {l})");
        }
    }

    /// A scaffold's pool, rebuilt from the buckets, has the bytes a pool
    /// file of its own had (one member per flush), and the reads in order.
    #[test]
    fn read_pools_in_buckets() {
        let dir = tempfile::tempdir().unwrap();
        let wd = dir.path().to_str().unwrap();
        let names: Vec<String> = (0..5).map(|i| format!("opera_scaffold_{i}")).collect();
        let mut pools = Pools::new(wd, &names, 1).unwrap();
        assert_eq!(pools.ends.len(), 4);
        let mut expected: HashMap<&str, Vec<u8>> = HashMap::new();
        for round in 0..3 {
            for (i, name) in names.iter().enumerate().take(4) {
                let record = format!("@r{round}_{i}\nACGT\n+\nIIII\n");
                pools.add(name, &record).unwrap();
                expected
                    .entry(name)
                    .or_default()
                    .extend(gzip_member(record.as_bytes()).unwrap());
            }
            pools.flush().unwrap();
        }
        let out = format!("{wd}/pool.fastq.gz");
        for name in &names[..4] {
            assert!(pools.write_pool(name, &out).unwrap());
            assert_eq!(fs::read(&out).unwrap(), expected[name.as_str()]);
            assert_eq!(pools.pool_size(name), expected[name.as_str()].len() as u64);
        }
        let mut text = String::new();
        flate2::read::MultiGzDecoder::new(fs::File::open(&out).unwrap())
            .read_to_string(&mut text)
            .unwrap();
        assert_eq!(text.lines().filter(|l| l.starts_with('@')).count(), 3);
        assert!(text.starts_with("@r0_3\n") && text.contains("@r2_3\n"));
        // no reads: no pool
        assert!(!pools.write_pool(&names[4], &format!("{wd}/none")).unwrap());
        assert_eq!(pools.pool_size(&names[4]), 0);
        // 5 scaffolds, 4 bucket files
        assert_eq!(fs::read_dir(format!("{wd}/{POOL_DIR}")).unwrap().count(), 4);
        pools.remove().unwrap();
        assert!(!Path::new(&format!("{wd}/{POOL_DIR}")).exists());
    }

    #[test]
    fn scaffolds_to_fill() {
        let scaf = ">s1\tlength: 5000\tcov: 1\nc1\tBE\t100\t10\nc2\tEB\t100\t0\n>s2\tlength: 900\tcov: 1\nc3\tBE\t900\t0\n";
        let s = split_scaffold_file(scaf, 500.0);
        assert_eq!(s.order, ["s1", "s2"]);
        assert_eq!(s.to_fill.len(), 1);
        assert!(s.multi.contains("s1"));
        let edges = scaffold_edges(scaf);
        assert_eq!(edges["c1-vs-c2"].scaff, "s1");
        assert_eq!(edges["c1-vs-c2"].ori, "+");
    }

    #[test]
    fn coords_sorted_like_gnu_sort() {
        let coords = "/a/ref /a/query\nNUCMER\n\n[S1]\t[E1]\n10\t20\t1\t11\t11\t11\t100.00\t50\t11\ts2\tc2\n5\t9\t1\t5\t5\t5\t100.00\t50\t5\ts1\tc1\n1\t4\t1\t4\t4\t4\t100.00\t50\t4\ts2\tc3\n";
        let sorted = sort_coords(coords);
        assert_eq!(sorted[0], "");
        assert_eq!(sorted[1], "/a/ref /a/query");
        assert_eq!(sorted[2], "NUCMER");
        assert!(sorted[3].ends_with("s1\tc1"));
        assert!(sorted[4].ends_with("s2\tc3"));
    }

    #[test]
    fn reverse_complement_iupac() {
        assert_eq!(reverse_complement("ACGTNacgtn"), "nacgtNACGT");
        // K <-> M
        assert_eq!(reverse_complement("RYKM"), "KMRY");
    }

    #[test]
    fn tiling_coords_from_minimap2() {
        // c1 on s1 (+), c2 on s1 (-) next to it with an overlap, c3 inside
        // c1.s interval with a lower identity (dropped), c4 below 95 %
        let paf = "c1\t1000\t0\t1000\t+\ts1\t5000\t100\t1100\t990\t1000\t60\n\
                   c2\t800\t0\t800\t-\ts1\t5000\t1050\t1850\t792\t800\t60\n\
                   c3\t500\t0\t500\t+\ts1\t5000\t300\t800\t480\t500\t60\n\
                   c4\t900\t0\t900\t+\ts1\t5000\t2000\t2900\t800\t900\t60\n";
        let coords = coords_from_paf(paf, "ref.fa", "query.fa");
        let lines: Vec<&str> = coords.lines().collect();
        assert_eq!(lines.len(), 6);
        assert!(lines[3].starts_with("[S1]"));
        assert_eq!(
            lines[4],
            "101\t1100\t1\t1000\t1000\t1000\t99.00\t5000\t1000\ts1\tc1"
        );
        assert_eq!(
            lines[5],
            "1051\t1850\t800\t1\t800\t800\t99.00\t5000\t800\ts1\tc2"
        );
        // read back as get_paf_file does: c1 then c2
        let (paf, _) = get_paf_file(">s1\nc1\t+\nc2\t-\n", &coords, true);
        assert!(
            paf.starts_with("c1\t1000\t0\t1000\t+\ts1\t5000\t100\t1100"),
            "{paf}"
        );
        assert!(
            paf.contains("c2\t800\t0\t800\t-\ts1\t5000\t1050\t1850"),
            "{paf}"
        );
    }

    #[test]
    fn contigs_in_gap_whole_and_once() {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("c.fa");
        // (no newline at the end: the last base must stay)
        std::fs::write(&fa, ">a x\nACGT\nTT\n>b\nGG\n>c\nCCC\nA").unwrap();
        // c placed twice, a once (multi-line), b not placed
        let paf = "c\t4\nc\t4\na\t6\n";
        let out = get_contig_in_gap(paf, fa.to_str().unwrap()).unwrap();
        assert_eq!(out, ">a\nACGTTT\n>c\nCCCA\n");
    }
}
