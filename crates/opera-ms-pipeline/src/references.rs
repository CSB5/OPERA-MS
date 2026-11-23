//! The reference genomes of a run: `sylph profile` of the
//! short reads against the user's sylph database (`--sylph-db`, GTDB
//! species representatives) gives the genomes present in the sample (all
//! the genomes sylph reports); each one is taken from the genome cache
//! (`--genome-cache`) or downloaded from NCBI by its accession, and linked
//! into a per-run folder: `<accession>.fna.gz` files and
//! `genomes_length.txt`. The GTDB accession is the species label (each
//! representative is one species).
//!
//! Downloads: the NCBI datasets report of the accession gives its assembly
//! name, hence the FTP path of `<accession>_<assembly>_genomic.fna.gz`
//! (the reports of the genomes not in the store yet are looked up 100 per
//! request, an accession missing from them alone; `NCBI_API_KEY`, when
//! set, is sent with these requests and never written in the log);
//! a few downloads at a time; a file is written under a temporary name,
//! checked (complete gzip FASTA) and renamed, so an interrupted download
//! never looks complete; runs sharing a cache take a lock per genome
//! (`<accession>.lock`, created atomically) so that a genome is downloaded
//! once. Without `--genome-cache` the genomes go to a temporary folder of
//! the run, deleted at the end. A genome that can't be obtained is a
//! warning (skipped) or, with `--strict-genomes`, an error.

use std::collections::HashMap;
use std::fmt::Write as _;
use std::fs::{self, OpenOptions};
use std::io::{BufRead, Read};
use std::path::Path;
use std::time::{Duration, Instant};

use anyhow::{Context, Result, bail};
use opera_ms_core::bench::StepRecorder;
use opera_ms_core::parallel::run_parallel;
use opera_ms_core::runner::{Cmd, Output, Pipeline};

/// Parallel downloads (polite to NCBI's servers).
const DOWNLOADS: usize = 4;
/// How long a run waits for another run's download of the same genome.
const LOCK_WAIT: Duration = Duration::from_secs(1800);
const NCBI_REPORT: &str = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession";
const NCBI_FTP: &str = "https://ftp.ncbi.nlm.nih.gov/genomes/all";
/// Accessions per NCBI datasets report request.
const REPORT_BATCH: usize = 100;

pub struct Args<'a> {
    pub sylph_db: &'a str,
    pub read1: &'a str,
    pub read2: &'a str,
    /// The run's genome folder (`<inter>/reference_genomes`).
    pub out_dir: &'a str,
    /// `--genome-cache`, or `None`: a temporary folder of the run.
    pub cache: Option<&'a str>,
    pub strict: bool,
    /// `--offline`: a genome not in the store is missing (no NCBI request).
    pub offline: bool,
    pub nproc: usize,
    pub sylph_dir: &'a str,
    pub curl_dir: &'a str,
}

/// What `prepare` found.
#[derive(Debug, Default)]
pub struct Summary {
    pub candidates: usize,
    pub from_cache: usize,
    pub downloaded: usize,
    pub missing: Vec<String>,
}

/// Where downloaded genomes go when there is no `--genome-cache`.
pub fn temporary_store(out_dir: &str) -> String {
    format!("{out_dir}/downloaded")
}

/// `GCF_000188175.1` from a sylph genome file name
/// (`.../GCF_000188175.1_genomic.fna.gz`).
fn accession(genome_file: &str) -> Option<String> {
    let name = genome_file.rsplit('/').next()?;
    let mut parts = name.split('_');
    let (prefix, number) = (parts.next()?, parts.next()?);
    matches!(prefix, "GCF" | "GCA").then(|| format!("{prefix}_{number}"))
}

/// The accessions of a `sylph profile` TSV, sorted, each once (sylph's row
/// order varies between runs; the genome files and lists made from them
/// are then the same from run to run).
fn profile_accessions(profile: &str) -> Vec<String> {
    let mut lines = profile.lines();
    let header: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
    let Some(col) = header.iter().position(|h| *h == "Genome_file") else {
        return Vec::new();
    };
    let mut out: Vec<String> = lines
        .filter_map(|line| line.split('\t').nth(col).and_then(accession))
        .collect();
    out.sort();
    out.dedup();
    out
}

/// A genome file's content summary: sequences, total length, plasmid
/// (its first header names a plasmid). `Err`: not a complete gzip FASTA.
fn inspect_genome(path: &Path) -> Result<(usize, u64, bool)> {
    let mut reader = opera_ms_io::seq::open_maybe_gzip(path)?;
    let (mut records, mut length, mut plasmid) = (0, 0u64, false);
    let mut line = Vec::new();
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        if line.first() == Some(&b'>') {
            if records == 0 {
                plasmid = line.windows(7).any(|w| w == b"plasmid");
            }
            records += 1;
        } else {
            length += line.iter().filter(|b| !b.is_ascii_whitespace()).count() as u64;
        }
    }
    if records == 0 || length == 0 {
        bail!("{} has no sequence", path.display());
    }
    Ok((records, length, plasmid))
}

/// Accession -> assembly name, from an NCBI datasets report
/// (`reports[].accession`, `reports[].assembly_info.assembly_name`).
fn assembly_names(json: &str) -> Result<HashMap<String, String>> {
    let value: serde_json::Value =
        serde_json::from_str(json).context("the NCBI report is not JSON")?;
    let mut names = HashMap::new();
    for report in value["reports"].as_array().into_iter().flatten() {
        let accession = report["accession"].as_str();
        let name = report["assembly_info"]["assembly_name"].as_str();
        if let (Some(accession), Some(name)) = (accession, name) {
            names.insert(accession.to_string(), name.to_string());
        }
    }
    Ok(names)
}

/// `NCBI_API_KEY` (optional; raises the datasets API's limit from 5 to 10
/// requests per second) as a curl header file (`-H @file`, mode 0600), so
/// that the key never appears in the commands of the main log; the file
/// is removed when this is dropped.
struct ApiKey {
    header: Option<String>,
}

impl ApiKey {
    fn new(dir: &str) -> Result<ApiKey> {
        let key = std::env::var("NCBI_API_KEY").unwrap_or_default();
        let key = key.trim();
        if key.is_empty() {
            return Ok(ApiKey { header: None });
        }
        use std::io::Write as _;
        use std::os::unix::fs::OpenOptionsExt;
        let path = format!("{dir}/ncbi_api_key.header");
        let _ = fs::remove_file(&path);
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .mode(0o600)
            .open(&path)
            .and_then(|mut f| writeln!(f, "api-key: {key}"))
            .with_context(|| format!("cannot write {path}"))?;
        Ok(ApiKey { header: Some(path) })
    }

    /// curl's arguments for the key (none without it).
    fn curl_args(&self) -> Vec<String> {
        self.header
            .iter()
            .flat_map(|h| ["-H".to_string(), format!("@{h}")])
            .collect()
    }
}

impl Drop for ApiKey {
    fn drop(&mut self) {
        if let Some(path) = &self.header {
            let _ = fs::remove_file(path);
        }
    }
}

/// `curl` of `url` into `out` (retries, errors in `err`).
fn curl_get(
    rec: &StepRecorder<'_>,
    curl: &str,
    url: &str,
    out: &str,
    err: &str,
    extra: &[String],
) -> Result<()> {
    let pipeline = Pipeline::new(
        Cmd::new(curl)
            .args(["-sSfL", "--retry", "3", "--retry-delay", "5"])
            .args(extra)
            .args(["-o", out, url]),
    )
    .stdout(Output::Null)
    .stderr(Output::File(err.into()));
    let r = pipeline
        .run()
        .with_context(|| format!("cannot download {url}"))?;
    rec.record(&r);
    Ok(())
}

/// The assembly names of `accessions`, `REPORT_BATCH` per datasets request;
/// a failed request only leaves its accessions out (they are looked up
/// alone before their download).
fn lookup_names(
    rec: &StepRecorder<'_>,
    a: &Args<'_>,
    dir: &str,
    accessions: &[&String],
    key: &ApiKey,
) -> HashMap<String, String> {
    let curl = format!("{}curl", a.curl_dir);
    let mut names = HashMap::new();
    for (i, batch) in accessions.chunks(REPORT_BATCH).enumerate() {
        let list: Vec<&str> = batch.iter().map(|s| s.as_str()).collect();
        let url = format!(
            "{NCBI_REPORT}/{}/dataset_report?page_size={}",
            list.join(","),
            REPORT_BATCH
        );
        let out = format!("{dir}/ncbi_reports_{i}.json");
        let err = format!("{dir}/ncbi_reports_{i}.curl.err");
        let found = curl_get(rec, &curl, &url, &out, &err, &key.curl_args())
            .and_then(|()| fs::read_to_string(&out).with_context(|| format!("cannot read {out}")))
            .and_then(|json| assembly_names(&json));
        match found {
            Ok(found) => {
                names.extend(found);
                let _ = fs::remove_file(&err);
            }
            Err(e) => opera_ms_core::log::file_only(&format!(
                " *** NCBI reports {i}: {e:#} (looked up one by one)"
            )),
        }
        let _ = fs::remove_file(&out);
    }
    names
}

/// The NCBI FTP URL of an accession's genomic FASTA.
fn ftp_url(accession: &str, assembly_name: &str) -> Option<String> {
    let (prefix, number) = accession.split_once('_')?;
    let digits = number.split('.').next()?;
    if digits.len() != 9 {
        return None;
    }
    // NCBI's FTP replaces spaces and some punctuation by `_`
    let name: String = assembly_name
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || matches!(c, '.' | '-' | '_') {
                c
            } else {
                '_'
            }
        })
        .collect();
    let full = format!("{accession}_{name}");
    Some(format!(
        "{NCBI_FTP}/{prefix}/{}/{}/{}/{full}/{full}_genomic.fna.gz",
        &digits[0..3],
        &digits[3..6],
        &digits[6..9]
    ))
}

/// Downloads one genome into `store` (atomic), its assembly name from
/// `name` or else from its own NCBI report. `Err`: why it failed.
fn download(
    rec: &StepRecorder<'_>,
    a: &Args<'_>,
    store: &str,
    acc: &str,
    name: Option<&String>,
    key: &ApiKey,
) -> Result<()> {
    let curl = format!("{}curl", a.curl_dir);
    let err = format!("{store}/{acc}.curl.err");
    let name = match name {
        Some(name) => name.clone(),
        None => {
            let report = format!("{store}/{acc}.report.json.tmp");
            let url = format!("{NCBI_REPORT}/{acc}/dataset_report");
            curl_get(rec, &curl, &url, &report, &err, &key.curl_args())?;
            let json =
                fs::read_to_string(&report).with_context(|| format!("cannot read {report}"))?;
            let _ = fs::remove_file(&report);
            assembly_names(&json)?
                .remove(acc)
                .with_context(|| format!("{acc}: no assembly name in the NCBI report"))?
        }
    };
    let url = ftp_url(acc, &name).with_context(|| format!("{acc}: unexpected accession"))?;
    let tmp = format!("{store}/{acc}.fna.gz.tmp");
    curl_get(rec, &curl, &url, &tmp, &err, &[])?;
    inspect_genome(Path::new(&tmp)).with_context(|| format!("{acc}: incomplete download"))?;
    fs::rename(&tmp, format!("{store}/{acc}.fna.gz"))
        .with_context(|| format!("cannot rename {tmp}"))?;
    let _ = fs::remove_file(format!("{store}/{acc}.curl.err"));
    Ok(())
}

/// Makes sure `store/<acc>.fna.gz` exists: already there (true), or
/// downloaded under the genome's lock (false).
fn fetch(
    rec: &StepRecorder<'_>,
    a: &Args<'_>,
    store: &str,
    acc: &str,
    name: Option<&String>,
    key: &ApiKey,
) -> Result<bool> {
    let file = format!("{store}/{acc}.fna.gz");
    let lock = format!("{store}/{acc}.lock");
    let start = Instant::now();
    loop {
        if Path::new(&file).is_file() {
            return Ok(true);
        }
        match OpenOptions::new().write(true).create_new(true).open(&lock) {
            Ok(_) => break,
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => {
                if start.elapsed() > LOCK_WAIT {
                    bail!(
                        "{acc}: {lock} held for more than {} min",
                        LOCK_WAIT.as_secs() / 60
                    );
                }
                std::thread::sleep(Duration::from_secs(5));
            }
            Err(e) => return Err(e).with_context(|| format!("cannot create {lock}")),
        }
    }
    let result = if Path::new(&file).is_file() {
        Ok(true)
    } else if a.offline {
        Err(anyhow::anyhow!("{acc}: not in {store} (--offline)"))
    } else {
        download(rec, a, store, acc, name, key).map(|()| false)
    };
    let _ = fs::remove_file(&lock);
    result
}

/// Profiles the sample, fetches its genomes and writes the run's genome
/// folder. Returns what happened, for the log.
pub fn prepare(rec: &StepRecorder<'_>, a: &Args<'_>) -> Result<Summary> {
    let out = a.out_dir;
    fs::create_dir_all(out).with_context(|| format!("cannot create {out}"))?;
    let profile_path = format!("{out}/sylph_profile.tsv");
    let pipeline = Pipeline::new(Cmd::new(format!("{}sylph", a.sylph_dir)).args([
        "profile",
        "-t",
        &a.nproc.max(1).to_string(),
        a.sylph_db,
        "-1",
        a.read1,
        "-2",
        a.read2,
    ]))
    .stdout(Output::File(profile_path.clone().into()))
    .stderr(Output::File(format!("{out}/sylph_profile.err").into()));
    let report = pipeline.run().context("Error during sylph profile")?;
    rec.record(&report);
    let accessions = profile_accessions(
        &fs::read_to_string(&profile_path).context("cannot read the sylph profile")?,
    );

    let store = match a.cache {
        Some(dir) => dir.to_string(),
        None => temporary_store(out),
    };
    fs::create_dir_all(&store).with_context(|| format!("cannot create {store}"))?;
    let key = ApiKey::new(out)?;
    let to_download: Vec<&String> = accessions
        .iter()
        .filter(|acc| !Path::new(&format!("{store}/{acc}.fna.gz")).is_file())
        .collect();
    let names = if a.offline {
        HashMap::new()
    } else {
        lookup_names(rec, a, out, &to_download, &key)
    };
    if !to_download.is_empty() && !a.offline {
        opera_ms_core::log::file_only(&format!(
            " *** NCBI reports: {} of {} genomes to download found in {} requests",
            names.len(),
            to_download.len(),
            to_download.len().div_ceil(REPORT_BATCH)
        ));
    }
    let results: Vec<std::sync::Mutex<Option<Result<bool>>>> = accessions
        .iter()
        .map(|_| std::sync::Mutex::new(None))
        .collect();
    let jobs: Vec<_> = accessions
        .iter()
        .zip(&results)
        .map(|(acc, slot)| {
            let (store, names, key) = (&store, &names, &key);
            move || -> Result<()> {
                let r = fetch(rec, a, store, acc, names.get(acc), key);
                if let Ok(mut s) = slot.lock() {
                    *s = Some(r);
                }
                Ok(())
            }
        })
        .collect();
    run_parallel(jobs, DOWNLOADS).map_err(|e| anyhow::anyhow!("{e}"))?;

    let genome_dir = format!("{out}/genomes");
    if Path::new(&genome_dir).exists() {
        fs::remove_dir_all(&genome_dir).with_context(|| format!("cannot remove {genome_dir}"))?;
    }
    fs::create_dir_all(&genome_dir).with_context(|| format!("cannot create {genome_dir}"))?;
    let mut summary = Summary {
        candidates: accessions.len(),
        ..Summary::default()
    };
    let mut lengths = String::new();
    let mut plasmids = String::new();
    for (acc, slot) in accessions.iter().zip(results) {
        let result = slot.into_inner().ok().flatten();
        match result {
            Some(Ok(cached)) => {
                let src = format!("{store}/{acc}.fna.gz");
                let (records, length, plasmid) = inspect_genome(Path::new(&src))?;
                let link = format!("{genome_dir}/{acc}.fna.gz");
                let target =
                    fs::canonicalize(&src).with_context(|| format!("cannot find {src}"))?;
                std::os::unix::fs::symlink(&target, &link)
                    .with_context(|| format!("cannot link {link}"))?;
                // genomes_length.txt: name, sequences, length
                let _ = writeln!(lengths, "{acc}.fna.gz\t{records}\t{length}");
                if plasmid {
                    let _ = writeln!(plasmids, "{acc}.fna.gz");
                }
                if cached {
                    summary.from_cache += 1;
                } else {
                    summary.downloaded += 1;
                }
            }
            Some(Err(e)) => {
                if a.strict {
                    return Err(e.context(format!("genome {acc} (--strict-genomes)")));
                }
                opera_ms_core::note!(" *** warning: genome {acc} skipped: {e:#}");
                if key.header.is_some() && !a.offline && summary.missing.is_empty() {
                    opera_ms_core::note!(
                        " *** (NCBI_API_KEY is set: NCBI refuses every request with a wrong key)"
                    );
                }
                summary.missing.push(acc.clone());
            }
            None => bail!("genome {acc}: no result"),
        }
    }
    fs::write(format!("{genome_dir}/genomes_length.txt"), lengths)
        .context("cannot write genomes_length.txt")?;
    fs::write(format!("{genome_dir}/plasmids.txt"), plasmids)
        .context("cannot write plasmids.txt")?;
    Ok(summary)
}

/// The genome files of the run's folder, sorted (skani references).
pub fn genome_files(out_dir: &str) -> Result<Vec<String>> {
    let dir = format!("{out_dir}/genomes");
    let mut files: Vec<String> = fs::read_dir(&dir)
        .with_context(|| format!("cannot list {dir}"))?
        .filter_map(|e| e.ok().and_then(|e| e.file_name().into_string().ok()))
        .filter(|n| n.ends_with(".fna.gz"))
        .map(|n| format!("{dir}/{n}"))
        .collect();
    files.sort();
    Ok(files)
}

/// `skani dist` of `queries` against the run's genomes, written in the
/// `mash dist` layout the next steps read (reference, query, distance =
/// 1 - ANI/100, p-value, shared hashes), at most `top` hits per query
/// (the closest; ties to the first genomes by name).
pub fn skani_as_mash(
    rec: &StepRecorder<'_>,
    skani_dir: &str,
    genomes: &[String],
    queries: &[String],
    top: usize,
    nproc: usize,
    work: &str,
) -> Result<String> {
    if genomes.is_empty() || queries.is_empty() {
        return Ok(String::new());
    }
    let ref_list = format!("{work}.refs.txt");
    let query_list = format!("{work}.queries.txt");
    fs::write(&ref_list, genomes.join("\n") + "\n")
        .with_context(|| format!("cannot write {ref_list}"))?;
    fs::write(&query_list, queries.join("\n") + "\n")
        .with_context(|| format!("cannot write {query_list}"))?;
    let tsv = format!("{work}.skani.tsv");
    let pipeline = Pipeline::new(Cmd::new(format!("{skani_dir}skani")).args([
        "dist",
        "-t",
        &nproc.max(1).to_string(),
        "--rl",
        &ref_list,
        "--ql",
        &query_list,
        "-o",
        &tsv,
    ]))
    .stdout(Output::Null)
    .stderr(Output::File(format!("{work}.skani.err").into()));
    let report = pipeline.run().context("Error during skani dist")?;
    rec.record(&report);
    let mut text = String::new();
    fs::File::open(&tsv)
        .and_then(|mut f| f.read_to_string(&mut text))
        .with_context(|| format!("cannot read {tsv}"))?;
    // a fixed order (query, distance, genome), then the first `top` hits of
    // each query: skani's threads write tied hits in any order, and its own
    // `-n` keeps any `top` of them, so ties are cut here, always to the
    // first genomes by name
    let mut hits: Vec<(String, String, f64)> = Vec::new();
    for line in text.lines().skip(1) {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 3 {
            continue;
        }
        let ani: f64 = f[2].parse().unwrap_or(0.0);
        hits.push((f[1].to_string(), f[0].to_string(), 1.0 - ani / 100.0));
    }
    let out = mash_lines(hits, top);
    for path in [&ref_list, &query_list] {
        let _ = fs::remove_file(path);
    }
    Ok(out)
}

/// skani hits (query, genome, distance) as `mash dist` lines, sorted by
/// query, distance, genome, at most `top` per query.
fn mash_lines(mut hits: Vec<(String, String, f64)>, top: usize) -> String {
    hits.sort_by(|a, b| a.0.cmp(&b.0).then(a.2.total_cmp(&b.2)).then(a.1.cmp(&b.1)));
    let mut out = String::new();
    let mut kept = 0;
    for (i, (query, genome, dist)) in hits.iter().enumerate() {
        if i == 0 || hits[i - 1].0 != *query {
            kept = 0;
        }
        if kept < top {
            kept += 1;
            let _ = writeln!(out, "{genome}\t{query}\t{dist:.6}\t0\t0/0");
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn accessions_from_sylph() {
        let profile = "Sample_file\tGenome_file\tTaxonomic_abundance\n\
                       s\tgtdb_genomes_reps_r232/database/GCF/000/188/175/GCF_000188175.1_genomic.fna.gz\t14\n\
                       s\tdb/GCA/900/755/095/GCA_900755095.1_genomic.fna.gz\t11\n\
                       s\tdb/GCF/000/188/175/GCF_000188175.1_genomic.fna.gz\t1\n";
        assert_eq!(
            profile_accessions(profile),
            vec!["GCA_900755095.1".to_string(), "GCF_000188175.1".to_string()]
        );
    }

    #[test]
    fn top_hits_cut_ties_by_genome_name() {
        // skani's order, ties included: 3 genomes at distance 0 for q1
        let hit = |q: &str, g: &str, d: f64| (q.to_string(), g.to_string(), d);
        let hits = vec![
            hit("q2", "g9", 0.05),
            hit("q1", "gC", 0.0),
            hit("q1", "gD", 0.01),
            hit("q1", "gA", 0.0),
            hit("q2", "g1", 0.05),
            hit("q1", "gB", 0.0),
            hit("q2", "g5", 0.02),
        ];
        let expected = "gA\tq1\t0.000000\t0\t0/0\ngB\tq1\t0.000000\t0\t0/0\n\
                        g5\tq2\t0.020000\t0\t0/0\ng1\tq2\t0.050000\t0\t0/0\n";
        assert_eq!(mash_lines(hits.clone(), 2), expected);
        let mut reversed = hits;
        reversed.reverse();
        assert_eq!(mash_lines(reversed, 2), expected);
        assert_eq!(mash_lines(Vec::new(), 2), "");
    }

    #[test]
    fn ncbi_paths() {
        assert_eq!(
            ftp_url("GCF_000188175.1", "ASM18817v1").as_deref(),
            Some(
                "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/175/GCF_000188175.1_ASM18817v1/GCF_000188175.1_ASM18817v1_genomic.fna.gz"
            )
        );
        assert_eq!(
            ftp_url("GCA_964372435.1", "ERP125453 co1/307")
                .as_deref()
                .map(|u| u.ends_with("GCA_964372435.1_ERP125453_co1_307_genomic.fna.gz")),
            Some(true)
        );
    }

    #[test]
    fn names_from_a_report_of_several_accessions() {
        // nested "accession" fields (BioProject, BioSample) are not reports
        let json = r#"{"reports":[
            {"accession":"GCF_000005845.2","paired_accession":"GCA_000005845.2",
             "assembly_info":{"assembly_name":"ASM584v2",
                              "bioproject_lineage":[{"bioprojects":[{"accession":"PRJNA225"}]}],
                              "biosample":{"accession":"SAMN02604091"}}},
            {"accession":"GCA_900755095.1","assembly_info":{"assembly_name":"ERS473029_19"}},
            {"accession":"GCA_1","assembly_info":{}}],"total_count":3}"#;
        let names = assembly_names(json).unwrap();
        assert_eq!(names.len(), 2);
        assert_eq!(names["GCF_000005845.2"], "ASM584v2");
        assert_eq!(names["GCA_900755095.1"], "ERS473029_19");
        assert!(assembly_names("not json").is_err());
        assert!(assembly_names("{}").unwrap().is_empty());
    }
}
