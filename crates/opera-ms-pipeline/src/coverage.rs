//! Short-read mapping and coverage estimation.
//!
//! From `cov_estimate.config`: writes the bundling config
//! (`lib_1_bundles/lib_1.ebconfig`) and the coverage sigma config, the read
//! size, builds the minibwa index of the contigs, and runs
//! `minibwa map -u | short_read_analysis` on the interleaved, renamed read
//! pairs (`@N.1` / `@N.2`).
//!
//! Errors stop the run (missing files, colour-space reads, a read file
//! shorter than the other, truncated input, any process of the pipeline
//! failing). Read files are decompressed according to their content (gzip
//! magic bytes), not their `.gz` name.

use std::collections::HashMap;
use std::fs;
use std::io::{self, BufRead, BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result, bail, ensure};
use opera_ms_core::bench::StepRecorder;
use opera_ms_core::runner::{Cmd, Output, Pipeline};
use opera_ms_io::seq::open_maybe_gzip;

/// Number of reads looked at for the read size (`get_read_size`).
const READ_SIZE_SAMPLE: usize = 1000;

/// `minibwa map -K`: bases per batch, fixed so that the mapping doesn't
/// depend on the number of threads.
pub const BWA_BATCH_BASES: &str = "200000000";

/// The settings of `cov_estimate.config`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CovConfig {
    pub contigs_file: String,
    pub lib: String,
    pub mapping_files: String,
    pub output_dir: String,
    pub samtools_dir: String,
    pub tools_dir: String,
    pub kmer_size: String,
    pub contig_len_thr: Option<String>,
    pub contig_edge_len: Option<String>,
    pub contig_window_len: Option<String>,
    pub illumina_read_1: String,
    pub illumina_read_2: String,
    /// `NUM_PROCESSOR`: threads for the mapping.
    pub num_processor: String,
}

/// The absolute path (Perl's `File::Spec->rel2abs` on Unix), then `canonpath`
/// (repeated `/` collapsed, `/./` removed, trailing `/` removed).
pub fn rel2abs(path: &str, cwd: &str) -> String {
    let full = if path.starts_with('/') {
        path.to_string()
    } else {
        format!("{cwd}/{path}")
    };
    let parts: Vec<&str> = full
        .split('/')
        .filter(|p| !p.is_empty() && *p != ".")
        .collect();
    format!("/{}", parts.join("/"))
}

/// Parses `cov_estimate.config`: first two
/// whitespace-separated fields of each line, the last occurrence of a key
/// wins, `#` lines are comments.
pub fn parse_cov_config(text: &str, cwd: &str) -> Result<CovConfig> {
    let mut values: HashMap<&str, &str> = HashMap::new();
    for line in text.lines() {
        if line.starts_with('#') {
            continue;
        }
        let mut fields = line.split_whitespace();
        if let (Some(key), Some(value)) = (fields.next(), fields.next()) {
            values.insert(key, value);
        }
    }
    let get = |key: &str| -> Result<String> {
        values
            .get(key)
            .map(|v| v.to_string())
            .with_context(|| format!("cov_estimate.config: {key} missing"))
    };
    Ok(CovConfig {
        contigs_file: rel2abs(&get("CONTIGS_FILE")?, cwd),
        lib: get("LIB")?,
        mapping_files: get("MAPPING_FILES")?,
        output_dir: values
            .get("OUTPUT_DIR")
            .map_or("OperaMS".to_string(), |v| v.to_string()),
        samtools_dir: get("SAMTOOLS_DIR")?,
        tools_dir: get("TOOLS_DIR")?,
        kmer_size: get("KMER_SIZE")?,
        contig_len_thr: values.get("CONTIG_LEN_THR").map(|v| v.to_string()),
        contig_edge_len: values.get("CONTIG_EDGE_LEN").map(|v| v.to_string()),
        contig_window_len: values.get("CONTIG_WINDOW_LEN").map(|v| v.to_string()),
        illumina_read_1: get("ILLUMINA_READ_1")?,
        illumina_read_2: get("ILLUMINA_READ_2")?,
        num_processor: get("NUM_PROCESSOR")?,
    })
}

/// Edge bundle size threshold (`EDGE_BUNDLESIZE_THRESHOLD`, default 2).
const BUNDLE_SIZE_THRESHOLD: &str = "2";

/// `generate__bundler_config_file`: `lib_1.ebconfig`.
pub fn bundler_config(c: &CovConfig, lib_bundles_dir: &str) -> String {
    format!(
        "\nsamtools_dir={}\n\noutput_folder={lib_bundles_dir}\n\ncontig_file={}\n\n[LIB]\nmap_file={}\
         \nfilter_repeat=no\n\ncluster_threshold={BUNDLE_SIZE_THRESHOLD}\n\ncluster_increased_step=5\n\
         \nkmer={}\n",
        c.samtools_dir, c.contigs_file, c.lib, c.kmer_size
    )
}

/// The coverage sigma config (no precomputed
/// coverage: the coverage is written to `contigs_<window>_<edge>`).
pub fn coverage_sigma_config(c: &CovConfig, sigma_dir: &str, edges_files: &str) -> String {
    let window = c.contig_window_len.as_deref().unwrap_or("0");
    let edge = c.contig_edge_len.as_deref().unwrap_or("0");
    format!(
        "\nmapping_files={}\n\nsigma_contigs_file={sigma_dir}/contigs_{window}_{edge}\n\
         \ncontigs_file_type=SOAP\n\ncontigs_file={}\n\nedges_files={edges_files}\n\
         \noutput_dir={sigma_dir}\n\nedge_bundle_size_thr={BUNDLE_SIZE_THRESHOLD}\n\nkmer_size={}\n\
         \ncontig_len_thr={}\n\ncontig_edge_len={edge}\n\ncontig_window_len={window}\n\
         \npdist_type=NegativeBinomial\n\nsamtools_dir={}\n",
        c.mapping_files,
        c.contigs_file,
        c.kmer_size,
        c.contig_len_thr.as_deref().unwrap_or("500"),
        c.samtools_dir
    )
}

fn read_line(input: &mut dyn BufRead, line: &mut Vec<u8>) -> io::Result<bool> {
    line.clear();
    Ok(input.read_until(b'\n', line)? > 0)
}

/// `get_read_size`: the longest sequence line (newline included) of the
/// first 1000 records, read as 4-line records.
pub fn read_size(reads: &mut dyn BufRead) -> io::Result<usize> {
    let mut size = 0;
    let mut line = Vec::new();
    for _ in 0..READ_SIZE_SAMPLE {
        if !read_line(reads, &mut line)? {
            break;
        }
        if !read_line(reads, &mut line)? {
            break;
        }
        size = size.max(line.len());
        read_line(reads, &mut line)?;
        read_line(reads, &mut line)?;
    }
    Ok(size)
}

/// FASTA or FASTQ, from the first line that isn't a `#` comment
/// (`checkReadFormat`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadFormat {
    Fasta,
    Fastq,
}

/// `checkReadFormat`; colour-space reads (a digit 0-3 in the line after
/// the header) aren't supported.
pub fn read_format(reads: &mut dyn BufRead, name: &str) -> Result<ReadFormat> {
    let mut line = Vec::new();
    let format = loop {
        if !read_line(reads, &mut line)? {
            bail!("{name}: reads format is not correct (empty file)");
        }
        match line.first() {
            Some(b'#') => continue,
            Some(b'>') => break ReadFormat::Fasta,
            Some(b'@') => break ReadFormat::Fastq,
            _ => bail!("{name}: reads format is not correct"),
        }
    };
    read_line(reads, &mut line)?;
    if line.iter().any(|b| (b'0'..=b'3').contains(b)) {
        bail!("{name}: colour-space reads are not supported (bwa 0.7 has no colour-space index)");
    }
    Ok(format)
}

/// Reads the next line that doesn't start with `#` (the record header).
fn next_header(input: &mut dyn BufRead, line: &mut Vec<u8>) -> io::Result<bool> {
    loop {
        if !read_line(input, line)? {
            return Ok(false);
        }
        if line.first() != Some(&b'#') {
            return Ok(true);
        }
    }
}

fn truncated(what: &str) -> io::Error {
    io::Error::new(
        io::ErrorKind::UnexpectedEof,
        format!("{what}: truncated record"),
    )
}

/// Copies one record (its header already read) as `N.<mate>`.
fn copy_record(
    input: &mut dyn BufRead,
    line: &mut Vec<u8>,
    num: u64,
    mate: u8,
    format: ReadFormat,
    out: &mut dyn Write,
) -> io::Result<()> {
    let name = if mate == 1 { "read 1" } else { "read 2" };
    match format {
        ReadFormat::Fasta => {
            out.write_all(format!(">{num}.{mate}\n").as_bytes())?;
            if !read_line(input, line)? {
                return Err(truncated(name));
            }
            out.write_all(line)
        }
        ReadFormat::Fastq => {
            out.write_all(format!("@{num}.{mate}\n").as_bytes())?;
            if !read_line(input, line)? {
                return Err(truncated(name));
            }
            out.write_all(line)?;
            if !read_line(input, line)? {
                return Err(truncated(name));
            }
            out.write_all(format!("+{num}.{mate}\n").as_bytes())?;
            if !read_line(input, line)? {
                return Err(truncated(name));
            }
            out.write_all(line)
        }
    }
}

/// `readAndMap`: read pairs interleaved and renamed `N.1` / `N.2`, the
/// sequence and quality lines copied as they are.
pub fn interleave(
    first: &mut dyn BufRead,
    second: &mut dyn BufRead,
    format: ReadFormat,
    out: &mut dyn Write,
) -> io::Result<()> {
    let mut line = Vec::new();
    let mut num: u64 = 1;
    while next_header(first, &mut line)? {
        copy_record(first, &mut line, num, 1, format, out)?;
        if !next_header(second, &mut line)? {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!("the read 2 file has fewer reads than the read 1 file (read {num})"),
            ));
        }
        copy_record(second, &mut line, num, 2, format, out)?;
        num += 1;
    }
    out.flush()
}

/// The commands the coverage estimation needs from the installation.
pub struct Programs<'a> {
    /// `opera-ms __worker short_read_analysis` as a command.
    pub short_read_analysis: Cmd,
    pub recorder: &'a StepRecorder<'a>,
}

fn log_and_run(rec: &StepRecorder<'_>, pipeline: Pipeline) -> Result<()> {
    match pipeline.run() {
        Ok(report) => {
            rec.record(&report);
            Ok(())
        }
        Err(err) => {
            if let Some(report) = err.report() {
                rec.record(report);
            }
            Err(err.into())
        }
    }
}

/// The coverage estimation of `config_file`.
pub fn estimate(config_file: &str, cwd: &str, programs: Programs<'_>) -> Result<()> {
    let text =
        fs::read_to_string(config_file).with_context(|| format!("cannot read {config_file}"))?;
    let c = parse_cov_config(&text, cwd)?;
    ensure!(
        Path::new(&c.contigs_file).exists(),
        "Contigs file: {} not found",
        c.contigs_file
    );
    ensure!(
        Path::new(&c.samtools_dir).exists(),
        "Samtools not found at: {}",
        c.samtools_dir
    );
    let results = format!("{}/intermediate_files", c.output_dir);
    fs::create_dir_all(&results).with_context(|| format!("cannot create {results}"))?;

    // bundling config
    let lib_bundles_dir = format!("{results}/lib_1_bundles");
    let edges_files = format!("{lib_bundles_dir}/clusters_opera");
    fs::create_dir_all(&lib_bundles_dir)
        .with_context(|| format!("cannot create {lib_bundles_dir}"))?;
    let bundler_conf = format!("{lib_bundles_dir}/lib_1.ebconfig");
    fs::write(&bundler_conf, bundler_config(&c, &lib_bundles_dir))
        .with_context(|| format!("cannot write {bundler_conf}"))?;

    // sigma config and read size
    let sigma_dir = format!("{results}/coverage_estimation");
    fs::create_dir_all(&sigma_dir).with_context(|| format!("cannot create {sigma_dir}"))?;
    let sigma_conf = format!("{sigma_dir}/sigma.config");
    fs::write(
        &sigma_conf,
        coverage_sigma_config(&c, &sigma_dir, &edges_files),
    )
    .with_context(|| format!("cannot write {sigma_conf}"))?;
    let read1 = &c.illumina_read_1;
    let size = read_size(&mut *open_maybe_gzip(Path::new(read1))?)
        .with_context(|| format!("cannot read {read1}"))?;
    let read_size_file = format!("{sigma_dir}/read_size.dat");
    fs::write(&read_size_file, format!("{size}\n"))
        .with_context(|| format!("cannot write {read_size_file}"))?;

    // the index, then the mapping pipeline
    let folder = format!("{sigma_dir}/");
    let preprocess_err = format!("{sigma_dir}/preprocess_reads.err");
    fs::File::create(&preprocess_err).with_context(|| format!("cannot write {preprocess_err}"))?;
    let format = read_format(&mut *open_maybe_gzip(Path::new(read1))?, read1)?;
    // minibwa from TOOLS_DIR (samtools is no longer part of the mapping)
    let in_dir = |dir: &str, name: &str| {
        if dir.ends_with('/') || dir.is_empty() {
            format!("{dir}{name}")
        } else {
            format!("{dir}/{name}")
        }
    };
    let contig_name = c.contigs_file.rsplit('/').next().unwrap_or(&c.contigs_file);
    let index = format!("{folder}{contig_name}");
    // minibwa: multithreaded index
    let mapper = in_dir(&c.tools_dir, "minibwa");
    let index_files = [format!("{index}.l2b"), format!("{index}.mbw")];
    if !Path::new(&index_files[0]).is_file() {
        log_and_run(
            programs.recorder,
            Pipeline::new(Cmd::new(&mapper).args([
                "index",
                "-t",
                &c.num_processor,
                &c.contigs_file,
                &index,
            ]))
            .stdout(Output::Append(preprocess_err.clone().into()))
            .stderr(Output::Append(preprocess_err.clone().into())),
        )
        .context("Error during minibwa index")?;
    }
    // single-end mapping of the interleaved pairs (the bundler pairs them
    // by name); a fixed batch size
    let map = Cmd::new(&mapper).args([
        "map",
        "-x",
        "sr",
        "-P",
        // only the mapped reads (was `| samtools view -S -h -F 0x4 -`)
        "-u",
        "-t",
        &c.num_processor,
        "-K",
        BWA_BATCH_BASES,
        &index,
        "-",
    ]);
    let pipeline = Pipeline::new(map.stderr(Output::Append(preprocess_err.clone().into())))
        .pipe(
            programs
                .short_read_analysis
                .args([&bundler_conf, &sigma_conf]),
        )
        .stdout(Output::File(
            format!("{folder}short_read_analysis.out").into(),
        ))
        .stderr(Output::File(
            format!("{folder}short_read_analysis.err").into(),
        ));
    let (r1, r2) = (c.illumina_read_1.clone(), c.illumina_read_2.clone());
    let mut first = open_maybe_gzip(Path::new(&r1))?;
    let mut second = open_maybe_gzip(Path::new(&r2))?;
    let result = pipeline.run_fed(move |pipe| {
        let mut out = BufWriter::with_capacity(1 << 20, pipe);
        interleave(&mut *first, &mut *second, format, &mut out)
    });
    match result {
        Ok(report) => programs.recorder.record(&report),
        Err(err) => {
            if let Some(report) = err.report() {
                programs.recorder.record(report);
            }
            return Err(err).context(format!(
                "Error during short read pre-processing. Please see {preprocess_err} \
                 {folder}short_read_analysis.out {folder}short_read_analysis.err for details"
            ));
        }
    }
    // the index is only needed for this mapping
    for path in &index_files {
        fs::remove_file(path).with_context(|| format!("cannot remove {path}"))?;
    }
    for required in [
        format!("{sigma_dir}/assembly_size.dat"),
        format!("{lib_bundles_dir}/lib.txt"),
    ] {
        ensure!(
            Path::new(&required).exists(),
            "Error during short read pre-processing: {required} was not written. Please see \
             {folder}short_read_analysis.out and {folder}short_read_analysis.err for details"
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    const RUN: &str = "/data/ecoli_run";

    fn fixture(name: &str) -> String {
        let root = env!("CARGO_MANIFEST_DIR");
        let path = format!("{root}/../opera-ms-io/tests/fixtures/cpp_config/{name}");
        std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"))
    }

    /// The `cov_estimate.config` of the E. coli reference run (fixture
    /// paths; the tools directory as in the reference).
    fn reference_config() -> CovConfig {
        let text = format!(
            "ILLUMINA_READ_1 /data/reads/ecoliK12.mixed.short.R1.fastq\n\
             ILLUMINA_READ_2 /data/reads/ecoliK12.mixed.short.R2.fastq\n\
             OUTPUT_DIR {RUN}\nCONTIG_EDGE_LEN 80\nCONTIG_WINDOW_LEN 340\nCONTIG_LEN_THR 500\n\
             KMER_SIZE 60\nCONTIGS_FILE /data/ecoli_megahit_contigs.fa\n\
             NUM_PROCESSOR 4\n\nLIB {RUN}/intermediate_files/read_mapping/opera.bam\n\n\
             MAPPING_FILES {RUN}/intermediate_files/read_mapping/opera.bam\n\
             SAMTOOLS_DIR /data/tools/\n\
             TOOLS_DIR /data/tools/\n"
        );
        parse_cov_config(&text, "/tmp").unwrap()
    }

    #[test]
    fn configs_match_the_reference_run() {
        let c = reference_config();
        let inter = format!("{RUN}/intermediate_files");
        assert_eq!(
            bundler_config(&c, &format!("{inter}/lib_1_bundles")),
            fixture("lib_1.ebconfig")
        );
        assert_eq!(
            coverage_sigma_config(
                &c,
                &format!("{inter}/coverage_estimation"),
                &format!("{inter}/lib_1_bundles/clusters_opera")
            ),
            fixture("coverage_estimation.sigma.config")
        );
    }

    #[test]
    fn last_contigs_file_wins_and_rel2abs() {
        let c = parse_cov_config(
            "CONTIGS_FILE a.fa\nLIB l\nMAPPING_FILES m\nSAMTOOLS_DIR s\nTOOLS_DIR t\nKMER_SIZE 60\n\
             ILLUMINA_READ_1 r1\nILLUMINA_READ_2 r2\nNUM_PROCESSOR 8\n\nCONTIGS_FILE /out//inter/./mega/final.fa\n",
            "/cwd",
        )
        .unwrap();
        assert_eq!(c.contigs_file, "/out/inter/mega/final.fa");
        assert_eq!(rel2abs("x/y.fa", "/cwd/"), "/cwd/x/y.fa");
        assert_eq!(c.output_dir, "OperaMS");
    }

    #[test]
    fn read_size_counts_the_newline() {
        let reads = "@a\nACGT\n+\nIIII\n@b\nACGTACG\n+\nIIIIIII\n";
        assert_eq!(read_size(&mut reads.as_bytes()).unwrap(), 8);
    }

    #[test]
    fn interleaving_renames_reads() {
        let r1 = "@x/1\nAC\n+x\nII\n#c\n@y/1\nGG\n+\nJJ\n";
        let r2 = "@x/2\nTT\n+\nKK\n@y/2\nCC\n+\nLL\n";
        let mut out = Vec::new();
        interleave(
            &mut r1.as_bytes(),
            &mut r2.as_bytes(),
            ReadFormat::Fastq,
            &mut out,
        )
        .unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "@1.1\nAC\n+1.1\nII\n@1.2\nTT\n+1.2\nKK\n@2.1\nGG\n+2.1\nJJ\n@2.2\nCC\n+2.2\nLL\n"
        );
        let mut out = Vec::new();
        interleave(
            &mut ">a\nAC\n".as_bytes(),
            &mut ">b\nGT\n".as_bytes(),
            ReadFormat::Fasta,
            &mut out,
        )
        .unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), ">1.1\nAC\n>1.2\nGT\n");
    }

    #[test]
    fn shorter_second_file_is_an_error() {
        let mut out = Vec::new();
        let err = interleave(
            &mut "@a\nA\n+\nI\n@b\nC\n+\nI\n".as_bytes(),
            &mut "@a\nA\n+\nI\n".as_bytes(),
            ReadFormat::Fastq,
            &mut out,
        )
        .unwrap_err();
        assert!(err.to_string().contains("fewer reads"), "{err}");
    }

    #[test]
    fn read_formats() {
        assert_eq!(
            read_format(&mut "#x\n@r\nACGT\n".as_bytes(), "r").unwrap(),
            ReadFormat::Fastq
        );
        assert_eq!(
            read_format(&mut ">r\nACGT\n".as_bytes(), "r").unwrap(),
            ReadFormat::Fasta
        );
        assert!(read_format(&mut "@r\nT0123\n".as_bytes(), "r").is_err());
        assert!(read_format(&mut "hello\n".as_bytes(), "r").is_err());
    }
}
