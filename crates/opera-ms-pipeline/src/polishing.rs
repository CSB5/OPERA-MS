//! Short-read polishing: POLCA's method as pypolca runs it,
//! with the variant calling made parallel, optionally after Polypolish.
//!
//! pypolca maps the short reads (`bwa mem`, `samtools sort`; here minibwa,
//! whose index is multithreaded), calls
//! variants with one `freebayes -m 0 --min-coverage 3 -R 0 -p 1 -F 0.2 -E 0`
//! over the whole assembly (single-threaded: most of its time), and applies
//! the well-supported ones (`fix_consensus_from_vcf`). Here the contigs are
//! split into groups of similar total length and one freebayes runs per
//! group (`--targets`), in parallel; the corrections are applied in Rust
//! with pypolca's rules in `--careful` mode (alt count >= 4 and >= 3 x the
//! ref count; multi-allelic sites skipped; per contig, from the last
//! variant to the first, each only if the reference still matches). The
//! output doesn't depend on the number of threads: fixed `minibwa map -K`,
//! and freebayes calls each contig on its own whatever its group.
//!
//! Polypolish (option `--polisher polypolish-pypolca`): every placement of
//! each read file mapped separately (`bwa mem -a` in Polypolish's
//! instructions; here `minibwa map -N 1000 --outn=1000 --outs=0 -p 0`, close
//! to it), `polypolish filter`, `polypolish polish`; its
//! output is then polished as above.

use std::collections::HashMap;
use std::fmt::Write as _;
use std::fs;

use anyhow::{Context, Result, bail};
use opera_ms_core::bench::StepRecorder;
use opera_ms_core::parallel::run_parallel;
use opera_ms_core::runner::{Cmd, Output, Pipeline};

use crate::coverage::BWA_BATCH_BASES;

/// The polishing method (`--polisher`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Polisher {
    /// POLCA's method (pypolca), parallel; the default.
    Pypolca,
    /// Polypolish, then POLCA's method.
    PolypolishPypolca,
}

impl Polisher {
    pub fn name(self) -> &'static str {
        match self {
            Polisher::Pypolca => "pypolca",
            Polisher::PolypolishPypolca => "polypolish-pypolca",
        }
    }
}

/// freebayes groups per thread (load balancing).
const GROUPS_PER_THREAD: usize = 4;

/// pypolca `--careful`: minimum alt allele count and alt / ref ratio.
const MIN_ALT: u64 = 4;
const MIN_RATIO: u64 = 3;

pub struct Args<'a> {
    /// Working directory (created, emptied of its temporary files at the end).
    pub dir: &'a str,
    /// The contigs to polish (FASTA).
    pub contigs: &'a str,
    pub read1: &'a str,
    pub read2: &'a str,
    pub nproc: usize,
    pub polisher: Polisher,
    /// Tool directories, with a trailing slash.
    pub minibwa_dir: &'a str,
    pub samtools_dir: &'a str,
    pub freebayes_dir: &'a str,
    pub polypolish_dir: &'a str,
}

/// minibwa's index files of `asm` (`<asm>.l2b`, `<asm>.mbw`).
fn index_files(asm: &str) -> [String; 2] {
    [format!("{asm}.l2b"), format!("{asm}.mbw")]
}

fn minibwa_index(rec: &StepRecorder<'_>, a: &Args<'_>, asm: &str, err: Output) -> Result<()> {
    let n = a.nproc.max(1).to_string();
    run(
        rec,
        Pipeline::new(
            Cmd::new(format!("{}minibwa", a.minibwa_dir)).args(["index", "-t", &n, asm, asm]),
        )
        .stdout(Output::Null)
        .stderr(err),
        "minibwa index (polishing)",
    )
}

/// What the polishing did, for the log.
#[derive(Debug, Default)]
pub struct Summary {
    pub variants: usize,
    pub applied: usize,
}

fn read(path: &str) -> Result<String> {
    fs::read_to_string(path).with_context(|| format!("cannot read {path}"))
}

fn write(path: &str, text: &str) -> Result<()> {
    fs::write(path, text).with_context(|| format!("cannot write {path}"))
}

fn remove(path: &str) -> Result<()> {
    match fs::remove_file(path) {
        Err(e) if e.kind() != std::io::ErrorKind::NotFound => {
            Err(e).with_context(|| format!("cannot remove {path}"))
        }
        _ => Ok(()),
    }
}

fn run(rec: &StepRecorder<'_>, pipeline: Pipeline, what: &str) -> Result<()> {
    let report = pipeline
        .run()
        .with_context(|| format!("Error during {what}"))?;
    rec.record(&report);
    Ok(())
}

/// Polishes `a.contigs`; returns the polished FASTA text.
pub fn polish(rec: &StepRecorder<'_>, a: &Args<'_>) -> Result<(String, Summary)> {
    fs::create_dir_all(a.dir).with_context(|| format!("cannot create {}", a.dir))?;
    let input = match a.polisher {
        Polisher::Pypolca => a.contigs.to_string(),
        Polisher::PolypolishPypolca => polypolish(rec, a)?,
    };
    pypolca(rec, a, &input)
}

/// Polypolish: returns the path of its polished assembly.
fn polypolish(rec: &StepRecorder<'_>, a: &Args<'_>) -> Result<String> {
    let d = a.dir;
    let asm = format!("{d}/polypolish_input.fa");
    fs::copy(a.contigs, &asm).with_context(|| format!("cannot copy {}", a.contigs))?;
    let n = a.nproc.max(1).to_string();
    let err = |name: &str| Output::File(format!("{d}/{name}.err").into());
    minibwa_index(rec, a, &asm, err("polypolish_minibwa_index"))?;
    // every placement, each read file on its own, as Polypolish needs
    for (reads, sam) in [(a.read1, "r1.sam"), (a.read2, "r2.sam")] {
        run(
            rec,
            Pipeline::new(Cmd::new(format!("{}minibwa", a.minibwa_dir)).args([
                "map",
                "-x",
                "sr",
                "-P",
                "-N",
                "1000",
                "--outn=1000",
                "--outs=0",
                "-p",
                "0",
                "-t",
                &n,
                "-K",
                BWA_BATCH_BASES,
                &asm,
                reads,
            ]))
            .stdout(Output::File(format!("{d}/{sam}").into()))
            .stderr(err("polypolish_minibwa")),
            "minibwa map (Polypolish)",
        )?;
    }
    let polypolish = format!("{}polypolish", a.polypolish_dir);
    run(
        rec,
        Pipeline::new(Cmd::new(&polypolish).args([
            "filter",
            "--in1",
            &format!("{d}/r1.sam"),
            "--in2",
            &format!("{d}/r2.sam"),
            "--out1",
            &format!("{d}/f1.sam"),
            "--out2",
            &format!("{d}/f2.sam"),
        ]))
        .stdout(Output::Null)
        .stderr(err("polypolish_filter")),
        "polypolish filter",
    )?;
    let out = format!("{d}/polypolish.fasta");
    run(
        rec,
        Pipeline::new(Cmd::new(&polypolish).args([
            "polish",
            &asm,
            &format!("{d}/f1.sam"),
            &format!("{d}/f2.sam"),
        ]))
        .stdout(Output::File(out.clone().into()))
        .stderr(err("polypolish_polish")),
        "polypolish polish",
    )?;
    for name in ["r1.sam", "r2.sam", "f1.sam", "f2.sam"] {
        remove(&format!("{d}/{name}"))?;
    }
    remove(&asm)?;
    for path in index_files(&asm) {
        remove(&path)?;
    }
    Ok(out)
}

/// The contigs of a FASTA text, in order: (name, sequence).
fn fasta_records(text: &str) -> Vec<(String, String)> {
    let mut records: Vec<(String, String)> = Vec::new();
    for line in text.lines() {
        if let Some(header) = line.strip_prefix('>') {
            let name = header.split_whitespace().next().unwrap_or("").to_string();
            records.push((name, String::new()));
        } else if let Some((_, seq)) = records.last_mut() {
            seq.push_str(line.trim_end());
        }
    }
    records
}

/// Contig indices in `groups` groups of similar total length (longest
/// first, each to the currently lightest group; ties by index):
/// deterministic.
fn group_contigs(lengths: &[usize], groups: usize) -> Vec<Vec<usize>> {
    let groups = groups.max(1).min(lengths.len().max(1));
    let mut order: Vec<usize> = (0..lengths.len()).collect();
    order.sort_by(|&x, &y| lengths[y].cmp(&lengths[x]).then(x.cmp(&y)));
    let mut out: Vec<(usize, Vec<usize>)> = vec![(0, Vec::new()); groups];
    for i in order {
        let lightest = (0..groups).min_by_key(|&g| (out[g].0, g)).unwrap_or(0);
        out[lightest].0 += lengths[i];
        out[lightest].1.push(i);
    }
    out.into_iter()
        .map(|(_, mut g)| {
            g.sort_unstable();
            g
        })
        .filter(|g| !g.is_empty())
        .collect()
}

/// POLCA's method on `input`.
fn pypolca(rec: &StepRecorder<'_>, a: &Args<'_>, input: &str) -> Result<(String, Summary)> {
    let d = a.dir;
    let asm = format!("{d}/assembly.fa");
    // one sequence line per contig, names = first word (as pypolca reads them)
    let records = fasta_records(&read(input)?);
    if records.is_empty() {
        bail!("{input}: no contig to polish");
    }
    let mut text = String::new();
    for (name, seq) in &records {
        let _ = writeln!(text, ">{name}\n{seq}");
    }
    write(&asm, &text)?;
    drop(text);
    let n = a.nproc.max(1).to_string();
    let err = |name: &str| Output::File(format!("{d}/{name}.err").into());
    minibwa_index(rec, a, &asm, err("minibwa_index"))?;
    let bam = format!("{d}/assembly.bam");
    run(
        rec,
        // paired-end mapping of the two read files
        Pipeline::new(Cmd::new(format!("{}minibwa", a.minibwa_dir)).args([
            "map",
            "-x",
            "sr",
            "-t",
            &n,
            "-K",
            BWA_BATCH_BASES,
            &asm,
            a.read1,
            a.read2,
        ]))
        .pipe(
            Cmd::new(format!("{}samtools", a.samtools_dir))
                .args(["sort", "-@", &n, "-o", &bam, "-"]),
        )
        .stdout(Output::Null)
        .stderr(err("minibwa")),
        "minibwa map | samtools sort (polishing)",
    )?;
    for (args, what) in [
        (vec!["index", "-@", &n, &bam], "samtools index"),
        (vec!["faidx", &asm], "samtools faidx"),
    ] {
        run(
            rec,
            Pipeline::new(Cmd::new(format!("{}samtools", a.samtools_dir)).args(args))
                .stdout(Output::Null)
                .stderr(err("samtools")),
            what,
        )?;
    }

    // freebayes, one per group of contigs, in parallel
    let lengths: Vec<usize> = records.iter().map(|(_, s)| s.len()).collect();
    // more groups than threads, the largest started first: the time of a
    // group depends on its coverage and variants, not only its length, so
    // finer groups balance the threads better (the calls don't change:
    // freebayes calls each contig on its own)
    let groups = group_contigs(&lengths, GROUPS_PER_THREAD * a.nproc.max(1));
    let freebayes = format!("{}freebayes", a.freebayes_dir);
    let mut jobs = Vec::new();
    let mut order: Vec<usize> = (0..groups.len()).collect();
    let group_length = |g: usize| groups[g].iter().map(|&i| lengths[i]).sum::<usize>();
    order.sort_by_key(|&g| (std::cmp::Reverse(group_length(g)), g));
    for g in order {
        let members = &groups[g];
        let bed = format!("{d}/group_{g}.bed");
        let mut text = String::new();
        for &i in members {
            let _ = writeln!(text, "{}\t0\t{}", records[i].0, lengths[i]);
        }
        write(&bed, &text)?;
        let vcf = format!("{d}/group_{g}.vcf");
        let pipeline = Pipeline::new(Cmd::new(&freebayes).args([
            "-m",
            "0",
            "--min-coverage",
            "3",
            "-R",
            "0",
            "-p",
            "1",
            "-F",
            "0.2",
            "-E",
            "0",
            "-b",
            &bam,
            "-f",
            &asm,
            "-t",
            &bed,
            "-v",
            &vcf,
        ]))
        .stdout(Output::Null)
        .stderr(Output::File(format!("{d}/group_{g}.freebayes.err").into()));
        jobs.push(move || -> Result<()> {
            let report = pipeline.run()?;
            rec.record(&report);
            Ok(())
        });
    }
    run_parallel(jobs, a.nproc.max(1))
        .map_err(|e| anyhow::anyhow!("Error during freebayes: {e}"))?;

    // the variant lines of every group, by contig (each contig is in one
    // group, its variants in position order)
    let mut variants: HashMap<String, Vec<String>> = HashMap::new();
    for g in 0..groups.len() {
        let vcf = read(&format!("{d}/group_{g}.vcf"))?;
        for line in vcf.lines().filter(|l| !l.starts_with('#') && !l.is_empty()) {
            let contig = line.split('\t').next().unwrap_or("").to_string();
            variants.entry(contig).or_default().push(line.to_string());
        }
    }
    let mut summary = Summary::default();
    let mut vcf_all = String::new();
    let mut out = String::new();
    for (name, seq) in &records {
        let lines = variants.remove(name).unwrap_or_default();
        for l in &lines {
            vcf_all.push_str(l);
            vcf_all.push('\n');
        }
        summary.variants += lines.len();
        let (fixed, applied) = apply_variants(seq, &lines);
        summary.applied += applied;
        let _ = writeln!(out, ">{name}\n{fixed}");
    }
    // keep the variants (small), remove the rest
    write(&format!("{d}/polishing.vcf"), &vcf_all)?;
    for g in 0..groups.len() {
        remove(&format!("{d}/group_{g}.bed"))?;
        remove(&format!("{d}/group_{g}.vcf"))?;
    }
    for path in [bam.clone(), format!("{bam}.bai"), format!("{asm}.fai")] {
        remove(&path)?;
    }
    remove(&asm)?;
    for path in index_files(&asm) {
        remove(&path)?;
    }
    if input != a.contigs {
        remove(input)?;
    }
    Ok((out, summary))
}

/// A sample field of a VCF line, by its FORMAT name.
fn sample_field<'a>(fields: &[&'a str], key: &str) -> Option<&'a str> {
    let format = fields.get(8)?.split(':');
    let values: Vec<&str> = fields.get(9)?.split(':').collect();
    format
        .enumerate()
        .find(|(_, k)| *k == key)
        .and_then(|(i, _)| values.get(i).copied())
}

/// pypolca's `fix_consensus_from_vcf` for one contig (`--careful`): the
/// supported single-allele variants, applied from the last to the first,
/// each only if the sequence at its position still matches its reference
/// allele (pypolca's test: a sequence with a base in `acgtnACGTN` that
/// differs, case-insensitively). Returns the sequence and the number of
/// variants applied.
fn apply_variants(seq: &str, lines: &[String]) -> (String, usize) {
    let mut fixes: Vec<(usize, &str, &str)> = Vec::new();
    for line in lines {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 10 || f[4].contains(',') {
            continue;
        }
        let count = |key: &str| {
            sample_field(&f, key)
                .and_then(|v| v.parse::<u64>().ok())
                .unwrap_or(0)
        };
        let (ref_count, alt_count) = (count("RO"), count("AO"));
        if alt_count < MIN_ALT || alt_count < MIN_RATIO * ref_count {
            continue;
        }
        let Ok(pos) = f[1].parse::<usize>() else {
            continue;
        };
        if pos == 0 {
            continue;
        }
        fixes.push((pos, f[3], f[4]));
    }
    let mut out = seq.to_string();
    let mut applied = 0;
    for &(pos, reference, alt) in fixes.iter().rev() {
        let start = (pos - 1).min(out.len());
        let end = (start + reference.len()).min(out.len());
        let original = &out[start..end];
        let mismatch = original.bytes().any(|c| b"acgtnACGTN".contains(&c))
            && !original.eq_ignore_ascii_case(reference);
        if mismatch {
            continue;
        }
        out.replace_range(start..end, alt);
        applied += 1;
    }
    (out, applied)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn vcf(contig: &str, pos: usize, r: &str, a: &str, ro: u64, ao: u64) -> String {
        format!(
            "{contig}\t{pos}\t.\t{r}\t{a}\t1\t.\t.\tGT:DP:AD:RO:QR:AO:QA:GL\t1:9:{ro},{ao}:{ro}:0:{ao}:0:0"
        )
    }

    #[test]
    fn variants_applied_with_careful_thresholds() {
        let seq = "ACGTACGTAC";
        let lines = vec![
            vcf("c", 2, "C", "T", 0, 5),   // applied
            vcf("c", 4, "T", "TT", 1, 3),  // alt < 4: skipped
            vcf("c", 6, "C", "G", 2, 5),   // 5 < 3 x 2: skipped
            vcf("c", 8, "T", "G,A", 0, 9), // multi-allelic: skipped
            vcf("c", 9, "AC", "A", 0, 10), // deletion, applied
            vcf("c", 3, "T", "A", 0, 10),  // reference doesn't match: skipped
        ];
        let (fixed, applied) = apply_variants(seq, &lines);
        assert_eq!(fixed, "ATGTACGTA");
        assert_eq!(applied, 2);
    }

    #[test]
    fn groups_are_balanced_and_deterministic() {
        let g = group_contigs(&[100, 10, 60, 50, 5], 2);
        // 100 -> g0, 60 -> g1, 50 -> g1 (110), 10 -> g0 (110), 5 -> g0 (tie)
        assert_eq!(g, vec![vec![0, 1, 4], vec![2, 3]]);
        assert_eq!(group_contigs(&[7], 4), vec![vec![0]]);
    }

    #[test]
    fn fasta_names_are_the_first_word() {
        let r = fasta_records(">a x\nAC\nGT\n>b\nTT\n");
        assert_eq!(
            r,
            vec![("a".into(), "ACGT".into()), ("b".into(), "TT".into())]
        );
    }
}
