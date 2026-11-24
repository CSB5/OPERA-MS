//! The pipeline: its stages, directories and resume logic.
//!
//! - every command is checked, and any failure stops the run;
//! - sigma and OPERA-LG run in `opera-ms __worker` processes;
//! - external tools come from `PATH`;
//! - resume uses done markers (see `opera_ms_core::stage`);
//! - `benchmark.tsv` in the output directory: time and memory per step.
//!
//! Paths are strings built as OPERA-MS 0.9 built them (double slashes
//! included), because they end up in the files the steps write.

use std::collections::HashMap;
use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use anyhow::{Context, Result, bail, ensure};
use opera_ms_core::bench::{Benchmark, Status, StepRecorder};
use opera_ms_core::parallel::run_parallel;
use opera_ms_core::runner::{Cmd, Output, Pipeline, RunError};
use opera_ms_core::stage::{Decision, Stages, perl_localtime};
use opera_ms_io::cpp_config;
use opera_ms_io::seq;

use crate::files;
use crate::options::{Options, VERSION};

/// Where the OPERA-MS installation's programs are (`check_dependency`).
#[derive(Debug, Clone)]
pub struct Install {
    /// `dirname(rel2abs($0)) . "/"`: the installation directory, with a
    /// trailing slash.
    pub opera_ms_dir: String,
    /// The `opera-ms` executable, for worker processes.
    pub worker: PathBuf,
}

impl Install {
    /// The installation of the running executable.
    pub fn of_current_exe() -> Result<Self> {
        let worker = std::env::current_exe().context("cannot find the opera-ms executable")?;
        let dir = worker
            .parent()
            .context("the opera-ms executable has no directory")?;
        Ok(Self {
            opera_ms_dir: format!("{}/", dir.display()),
            worker,
        })
    }
}

/// The tool directories (with a trailing slash), resolved from `PATH` by
/// `check_dependency` before the first step; empty for tools the run
/// doesn't need.
#[derive(Debug, Clone, Default)]
struct Tools {
    samtools: String,
    minibwa: String,
    minimap2: String,
    sylph: String,
    skani: String,
    curl: String,
    megahit: String,
    racon: String,
    freebayes: String,
    polypolish: String,
}

/// The settings after `setup_directory_path` (absolute paths as strings).
#[derive(Debug, Clone)]
struct Settings {
    output_dir: String,
    inter_dir: String,
    illumina_read_1: String,
    illumina_read_2: String,
    /// Plain or gzip-compressed.
    long_read: String,
    /// Changes during the run, as `$opera_ms_option{"CONTIGS_FILE"}`.
    contigs_file: Option<String>,
}

/// A run of the pipeline.
pub struct Run {
    opt: Options,
    install: Install,
    tools: Tools,
    s: Settings,
    stages: Stages,
    bench: Arc<Benchmark>,
    /// `ASSEMBLED_CONTIG_FILE`: the contigs come from MEGAHIT.
    assembled_contig_file: Option<bool>,
}

fn absolute(main_dir: &str, path: &str) -> String {
    if path.starts_with('/') {
        path.to_string()
    } else {
        format!("{main_dir}{path}")
    }
}

/// The first lines of a run's log: version, command line, machine.
fn log_header() {
    let args: Vec<String> = std::env::args().collect();
    let cpus = std::thread::available_parallelism().map_or(0, |n| n.get());
    let ram_gb = fs::read_to_string("/proc/meminfo")
        .ok()
        .and_then(|m| {
            m.lines()
                .find(|l| l.starts_with("MemTotal:"))
                .and_then(|l| l.split_whitespace().nth(1)?.parse::<f64>().ok())
        })
        .map_or(0.0, |kb| kb / 1024.0 / 1024.0);
    opera_ms_core::note!("==== OPERA-MS {VERSION}");
    opera_ms_core::note!("command: {}", args.join(" "));
    opera_ms_core::note!("machine: {cpus} CPUs, {ram_gb:.1} GB RAM");
}

/// `init_dir`: `rm -r $dir` if it exists, then `mkdir -p $dir`.
fn init_dir(dir: &str) -> Result<()> {
    let path = Path::new(dir);
    if path.is_dir() {
        fs::remove_dir_all(path).with_context(|| format!("cannot remove {dir}"))?;
    }
    fs::create_dir_all(path).with_context(|| format!("cannot create {dir}"))
}

fn read(path: &str) -> Result<String> {
    let bytes = fs::read(path).with_context(|| format!("cannot read {path}"))?;
    String::from_utf8(bytes).with_context(|| format!("{path} is not UTF-8 text"))
}

fn write(path: &str, content: &str) -> Result<()> {
    fs::write(path, content).with_context(|| format!("cannot write {path}"))
}

fn append(path: &str, content: &str) -> Result<()> {
    fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)
        .and_then(|mut f| f.write_all(content.as_bytes()))
        .with_context(|| format!("cannot write {path}"))
}

fn copy(from: &str, to: &str) -> Result<()> {
    fs::copy(from, to)
        .map(|_| ())
        .with_context(|| format!("cannot copy {from} to {to}"))
}

fn rename(from: &str, to: &str) -> Result<()> {
    fs::rename(from, to).with_context(|| format!("cannot move {from} to {to}"))
}

fn symlink(target: &str, link: &str) -> Result<()> {
    std::os::unix::fs::symlink(target, link)
        .with_context(|| format!("cannot link {link} to {target}"))
}

fn exists(path: &str) -> bool {
    Path::new(path).exists()
}

/// Files of `dir` whose name matches `prefix*suffix`, sorted bytewise (as
/// `dash` sorts a glob).
fn glob(dir: &str, matches: impl Fn(&str) -> bool) -> Result<Vec<String>> {
    let mut names: Vec<String> = fs::read_dir(dir)
        .with_context(|| format!("cannot list {dir}"))?
        .filter_map(|e| e.ok().and_then(|e| e.file_name().into_string().ok()))
        .filter(|n| matches(n))
        .collect();
    names.sort();
    Ok(names.into_iter().map(|n| format!("{dir}/{n}")).collect())
}

/// Runs a pipeline, logging it, and adds its usage to the step.
fn run(rec: &StepRecorder<'_>, pipeline: Pipeline) -> Result<(), RunError> {
    match pipeline.run() {
        Ok(report) => {
            rec.record(&report);
            Ok(())
        }
        Err(err) => {
            if let Some(report) = err.report() {
                rec.record(report);
            }
            Err(err)
        }
    }
}

fn seconds_since(start: Instant) -> u64 {
    start.elapsed().as_secs()
}

impl Run {
    /// `setup_directory_path`, `check_dependency`, `check_data`,
    /// `check_database`.
    pub fn prepare(opt: Options, install: Install) -> Result<Self> {
        let main_dir = format!(
            "{}/",
            std::env::current_dir()
                .context("cannot get the current directory")?
                .display()
        );
        let need = |value: &Option<String>, option: &str| -> Result<String> {
            value
                .clone()
                .with_context(|| format!("{option} is required"))
        };
        let output_dir = need(
            &opt.output_dir.as_ref().map(|p| p.display().to_string()),
            "--out-dir (OUTPUT_DIR)",
        )?;
        let output_dir = if output_dir.starts_with('/') {
            output_dir
        } else {
            format!("{main_dir}{output_dir}/")
        };
        let inter_dir = format!("{output_dir}/intermediate_files");
        fs::create_dir_all(&inter_dir).with_context(|| format!("cannot create {inter_dir}"))?;
        // the main log, from here on
        let log_path = format!("{output_dir}/opera-ms.log");
        opera_ms_core::log::open(Path::new(&log_path))
            .with_context(|| format!("cannot write {log_path}"))?;
        log_header();

        // --prepare-genomes needs only the short reads (the long
        // reads, if given, are ignored)
        let long_read = match &opt.long_read {
            _ if opt.prepare_genomes => String::new(),
            long_read => absolute(
                &main_dir,
                &need(
                    &long_read.as_ref().map(|p| p.display().to_string()),
                    "--long-read (LONG_READ)",
                )?,
            ),
        };
        let s = Settings {
            illumina_read_1: absolute(
                &main_dir,
                &need(&opt.illumina_read_1, "--short-read1 (ILLUMINA_READ_1)")?,
            ),
            illumina_read_2: absolute(
                &main_dir,
                &need(&opt.illumina_read_2, "--short-read2 (ILLUMINA_READ_2)")?,
            ),
            long_read,
            contigs_file: opt
                .contigs_file
                .as_ref()
                .map(|p| absolute(&main_dir, &p.display().to_string())),
            output_dir,
            inter_dir,
        };

        let stages = Stages::new(Path::new(&s.inter_dir), opt.nb_total_stage(), false);
        let bench = Arc::new(
            Benchmark::create(format!("{}/benchmark.tsv", s.output_dir))
                .context("cannot create benchmark.tsv")?,
        );
        let mut run = Self {
            opt,
            install,
            tools: Tools::default(),
            s,
            stages,
            bench,
            assembled_contig_file: None,
        };
        run.check_dependency()?;
        run.check_data()?;
        if run.opt.ref_clustering {
            match &run.opt.sylph_db {
                None => bail!(
                    "reference clustering needs a sylph database of GTDB genomes: give it with \
                     --sylph-db <file.syl2db>, or run without reference clustering \
                     (--no-ref-clustering)"
                ),
                Some(db) if !exists(db) => bail!("--sylph-db: {db} not found"),
                Some(_) => {}
            }
        }
        Ok(run)
    }

    /// The tools the run needs, from `PATH`; the version and directory of
    /// each go to the log; any tool not found stops the run before its
    /// first step.
    fn check_dependency(&mut self) -> Result<()> {
        use crate::deps::Tool;
        let needed_base = [Tool::Samtools, Tool::Minibwa, Tool::Minimap2];
        let mut needed = if self.opt.prepare_genomes {
            vec![Tool::Sylph, Tool::Curl]
        } else {
            needed_base.to_vec()
        };
        if self.s.contigs_file.is_none() && !self.opt.prepare_genomes {
            needed.push(Tool::Megahit);
        }
        if self.opt.ref_clustering && !self.opt.prepare_genomes {
            needed.extend([Tool::Sylph, Tool::Skani]);
            if !self.opt.offline {
                needed.push(Tool::Curl);
            }
        }
        if self.opt.gap_filling && !self.opt.prepare_genomes {
            needed.extend([Tool::Racon, Tool::Minimap2]);
        }
        if !self.opt.no_polishing && !self.opt.prepare_genomes {
            needed.push(Tool::Freebayes);
            if self.opt.polisher == crate::polishing::Polisher::PolypolishPypolca {
                needed.push(Tool::Polypolish);
            }
        }
        needed.sort_by_key(|t| t.name());
        needed.dedup();
        let mut missing = Vec::new();
        for tool in needed {
            let found = match crate::deps::resolve(tool) {
                Ok(found) => found,
                Err(problem) => {
                    missing.push(format!("{} ({problem})", tool.name()));
                    continue;
                }
            };
            opera_ms_core::note!(" *** {} {} ({})", tool.name(), found.version, found.dir);
            let dir = found.dir;
            let t = &mut self.tools;
            match tool {
                Tool::Samtools => t.samtools = dir,
                Tool::Minibwa => t.minibwa = dir,
                Tool::Minimap2 => t.minimap2 = dir,
                Tool::Sylph => t.sylph = dir,
                Tool::Skani => t.skani = dir,
                Tool::Curl => t.curl = dir,
                Tool::Megahit => t.megahit = dir,
                Tool::Racon => t.racon = dir,
                Tool::Freebayes => t.freebayes = dir,
                Tool::Polypolish => t.polypolish = dir,
            }
        }
        ensure!(
            missing.is_empty(),
            "required programs not functional in PATH: {}. See opera-ms check-dependency",
            missing.join("; ")
        );
        Ok(())
    }

    /// `check_data`: the reads exist. Gzipped long reads are read as they
    /// are.
    fn check_data(&mut self) -> Result<()> {
        for path in [
            &self.s.illumina_read_1,
            &self.s.illumina_read_2,
            &self.s.long_read,
        ] {
            if path.is_empty() && self.opt.prepare_genomes {
                continue;
            }
            ensure!(exists(path), "Error : {path} does not exist");
        }
        Ok(())
    }

    /// The settings part of `cov_estimate.config`.
    fn cov_estimate_config_base(&self) -> String {
        let mut values: HashMap<&str, String> = HashMap::new();
        values.insert("ILLUMINA_READ_1", self.s.illumina_read_1.clone());
        values.insert("ILLUMINA_READ_2", self.s.illumina_read_2.clone());
        values.insert("LONG_READ", self.s.long_read.clone());
        values.insert("OUTPUT_DIR", self.s.output_dir.clone());
        values.insert("CONTIG_EDGE_LEN", self.opt.contig_edge_len.to_string());
        values.insert("CONTIG_WINDOW_LEN", self.opt.contig_window_len.to_string());
        values.insert("CONTIG_LEN_THR", self.opt.contig_len_thr.to_string());
        values.insert("KMER_SIZE", self.opt.kmer_size.to_string());
        if let Some(contigs) = &self.s.contigs_file {
            values.insert("CONTIGS_FILE", contigs.clone());
        }
        values.insert("NUM_PROCESSOR", self.opt.num_processor.to_string());
        files::utils_config(&values)
    }

    fn dir(&self, name: &str) -> String {
        format!("{}/{name}", self.s.inter_dir)
    }

    fn mapping_dir(&self) -> String {
        self.dir("read_mapping")
    }
    fn cov_dir(&self) -> String {
        self.dir("coverage_estimation")
    }
    fn sigma_dir(&self) -> String {
        self.dir("hierarchical_clustering")
    }
    fn ref_clustering_dir(&self) -> String {
        self.dir("reference_clustering")
    }
    fn strain_dir(&self) -> String {
        self.dir("strain_analysis")
    }
    fn opera_lr_dir(&self) -> String {
        self.dir("opera_long_read")
    }
    /// The run's reference genomes: sylph profile, genome links,
    /// `genomes_length.txt`.
    fn reference_genomes_dir(&self) -> String {
        self.dir("reference_genomes")
    }
    fn cov_file(&self) -> String {
        format!("{}/{}", self.cov_dir(), self.opt.contigs_cov_file())
    }
    fn contigs_file(&self) -> Result<String> {
        self.s
            .contigs_file
            .clone()
            .context("no contig file (the short-read assembly didn't run)")
    }

    fn time_log(&self, step: &str, seconds: u64) -> Result<()> {
        append(
            &format!("{}/time.log", self.s.inter_dir),
            &files::time_log_line(step, seconds),
        )
    }

    fn worker(&self, program: &str) -> Cmd {
        Cmd::new(&self.install.worker).args(["__worker", program])
    }

    /// Runs one stage, with its benchmark row.
    fn step(
        &mut self,
        name: &str,
        run: impl FnOnce(&mut Self, &StepRecorder<'_>) -> Result<()>,
    ) -> Result<()> {
        let bench = Arc::clone(&self.bench);
        bench.step(name, |rec| run(self, rec))
    }

    /// The whole pipeline.
    pub fn execute(mut self) -> Result<()> {
        if self.opt.prepare_genomes {
            return self.prepare_genomes_only();
        }
        let result = self.execute_stages();
        // downloaded genomes without --genome-cache: deleted
        if self.opt.genome_cache.is_none() {
            let store = crate::references::temporary_store(&self.reference_genomes_dir());
            if Path::new(&store).is_dir() {
                if let Err(e) = fs::remove_dir_all(&store) {
                    opera_ms_core::note!(" *** warning: cannot remove {store}: {e}");
                }
            }
        }
        self.bench.finish(if result.is_ok() {
            Status::Ok
        } else {
            Status::Failed
        });
        result?;
        opera_ms_core::say!("\n[{}]\tOPERA-MS completed", perl_localtime());
        Ok(())
    }

    /// `--prepare-genomes`: the reference genomes of the sample
    /// in the genome cache, nothing else.
    fn prepare_genomes_only(mut self) -> Result<()> {
        opera_ms_core::say!(
            "\n[{}]\tReference genomes (--prepare-genomes)",
            perl_localtime()
        );
        let result = self.step("Reference genomes", |r, rec| {
            r.reference_genomes(rec).map(|_| ())
        });
        self.bench.finish(if result.is_ok() {
            Status::Ok
        } else {
            Status::Failed
        });
        result?;
        opera_ms_core::say!(
            "\n[{}]\tReference genomes ready in {}: run OPERA-MS with the same --genome-cache (and --offline on a machine without internet)",
            perl_localtime(),
            self.opt.genome_cache.as_deref().unwrap_or("")
        );
        Ok(())
    }

    /// The sample's reference genomes: sylph, then the cache or NCBI,
    /// linked into `reference_genomes/genomes`.
    fn reference_genomes(&self, rec: &StepRecorder<'_>) -> Result<crate::references::Summary> {
        let sylph_db = self.opt.sylph_db.clone().unwrap_or_default();
        let summary = rec
            .substep("genomes", |sub| {
                crate::references::prepare(
                    sub,
                    &crate::references::Args {
                        sylph_db: &sylph_db,
                        read1: &self.s.illumina_read_1,
                        read2: &self.s.illumina_read_2,
                        out_dir: &self.reference_genomes_dir(),
                        cache: self.opt.genome_cache.as_deref(),
                        strict: self.opt.strict_genomes,
                        offline: self.opt.offline,
                        nproc: usize::try_from(self.opt.num_processor).unwrap_or(1),
                        sylph_dir: &self.tools.sylph,
                        curl_dir: &self.tools.curl,
                    },
                )
            })
            .context("Error in during the reference genome selection")?;
        opera_ms_core::note!(
            " *** reference genomes: {} found by sylph, {} from the cache, {} downloaded, {} missing{}",
            summary.candidates,
            summary.from_cache,
            summary.downloaded,
            summary.missing.len(),
            if summary.downloaded > 0 && self.opt.genome_cache.is_none() {
                " (use --genome-cache <dir> to keep downloaded genomes for other runs)"
            } else {
                ""
            }
        );
        Ok(summary)
    }

    fn execute_stages(&mut self) -> Result<()> {
        self.short_read_assembly()?;
        self.assembly_graph_creation()?;
        if self.opt.single_genome {
            self.single_genome_assembly()?;
        } else {
            self.hierarchical_clustering()?;
            self.reference_clustering()?;
            self.strain_clustering_and_assembly()?;
        }
        self.gap_filling()?;
        let contig_info = self.generate_assembly_stats()?;
        self.write_final_assembly(&contig_info)?;
        self.polishing(&contig_info)
    }
}

/// Last line of a text (`tail -n1`), without its newline.
fn last_line(text: &str) -> &str {
    text.trim_end_matches('\n')
        .rsplit('\n')
        .next()
        .unwrap_or("")
}

/// Numeric value of a string as Perl reads it (leading number, else 0),
/// for the values OPERA-MS 0.9 wrote that way.
fn perl_float(text: &str) -> f64 {
    let text = text.trim_start();
    let mut end = 0;
    for (i, c) in text.char_indices() {
        let ok = c.is_ascii_digit() || c == '.' || (i == 0 && (c == '-' || c == '+'));
        if !ok {
            break;
        }
        end = i + c.len_utf8();
    }
    text[..end].parse().unwrap_or(0.0)
}

impl Run {
    /// One stage of `check_completed`: progress line, resume decision,
    /// benchmark row, and a check that its output exists before it's marked
    /// done.
    fn stage(
        &mut self,
        id: &str,
        message: &str,
        number: u32,
        skip: bool,
        done_file: &str,
        body: impl FnOnce(&mut Self, &StepRecorder<'_>) -> Result<()>,
    ) -> Result<()> {
        match self.stages.begin(id, message, number, skip) {
            Decision::Run => {
                self.stages.invalidate(id)?;
                let done_file = done_file.to_string();
                self.step(message, |run, rec| {
                    body(run, rec)?;
                    ensure!(
                        exists(&done_file),
                        "{message}: the stage finished but {done_file} is missing"
                    );
                    Ok(())
                })?;
                self.stages.complete(id)?;
            }
            Decision::Skip | Decision::AlreadyDone => self.bench.skipped(message),
        }
        Ok(())
    }

    // ---- stage 1 -------------------------------------------------------

    fn short_read_assembly(&mut self) -> Result<()> {
        if let Some(user_contigs) = self.s.contigs_file.clone() {
            self.assembled_contig_file = Some(false);
            let time = perl_localtime();
            opera_ms_core::say!(
                "\n[{time}]\tShort read assembly [1/{}]",
                self.opt.nb_total_stage()
            );
            opera_ms_core::say!("[{time}]\tSkip [contig file provided as input]");
            // process_user_assembly
            let dir = self.dir("user_assembly");
            init_dir(&dir)?;
            let out = format!("{dir}/contigs.fasta");
            let mut text = String::new();
            seq::open_maybe_gzip(Path::new(&user_contigs))?
                .read_to_string(&mut text)
                .with_context(|| format!("cannot read {user_contigs}"))?;
            write(&out, &files::multi_to_linear_fa(&text))?;
            self.s.contigs_file = Some(out);
            return Ok(());
        }
        self.assembled_contig_file = Some(true);
        let start = Instant::now();
        self.run_megahit()?;
        self.time_log("short_read_assembly", seconds_since(start))
    }

    fn run_megahit(&mut self) -> Result<()> {
        let inter = self.s.inter_dir.clone();
        let dir = self.dir("megahit_assembly");
        let contigs = format!("{dir}/final.contigs.fa");
        self.s.contigs_file = Some(contigs.clone());
        self.stage(
            "short_read_assembly",
            "Short read assembly [MEGAHIT]",
            1,
            false,
            &contigs,
            |r, rec| {
                if Path::new(&dir).is_dir() {
                    fs::remove_dir_all(&dir).with_context(|| format!("cannot remove {dir}"))?;
                }
                let cmd = Cmd::new(format!("{}/megahit", r.tools.megahit)).args([
                    "--num-cpu-threads",
                    &r.opt.num_processor.to_string(),
                    "-1",
                    &r.s.illumina_read_1,
                    "-2",
                    &r.s.illumina_read_2,
                    "-o",
                    &dir,
                ]);
                run(
                    rec,
                    Pipeline::new(cmd)
                        .stdout(Output::File(format!("{inter}/megahit.out").into()))
                        .stderr(Output::File(format!("{inter}/megahit.err").into())),
                )
                .context("Error during MEGAHIT assembly")?;
                let intermediate = format!("{dir}/intermediate_contigs");
                if Path::new(&intermediate).is_dir() {
                    fs::remove_dir_all(&intermediate)
                        .with_context(|| format!("cannot remove {intermediate}"))?;
                }
                Ok(())
            },
        )
    }

    // ---- stages 2 and 3 --------------------------------------------------

    fn assembly_graph_creation(&mut self) -> Result<()> {
        let mapping = self.mapping_dir();
        let done = format!("{mapping}/config");
        self.stage(
            "assembly_graph",
            "Long-read mapping and assembly graph generation",
            2,
            false,
            &done,
            |r, rec| {
                init_dir(&mapping)?;
                let start = Instant::now();
                let contigs = r.contigs_file()?;
                let (kmer, nproc) = (r.opt.kmer_size.to_string(), r.opt.num_processor.to_string());
                crate::long_read::opera_long_read(
                    rec,
                    &crate::long_read::Args {
                        contig_file: &contigs,
                        kmer: &kmer,
                        long_read_file: &r.s.long_read,
                        output_dir: &mapping,
                        nproc: &nproc,
                        minimap2_dir: &r.tools.minimap2,
                    },
                )
                .context("Error in the long read processing")?;
                r.time_log("graph_construction", seconds_since(start))
            },
        )?;

        let cov = self.cov_dir();
        let cov_file = self.cov_file();
        self.stage(
            "coverage_estimation",
            "Short-read mapping and coverage estimation",
            3,
            false,
            &cov_file.clone(),
            |r, rec| {
                fs::create_dir_all(&cov).with_context(|| format!("cannot create {cov}"))?;
                let config = format!("{cov}/cov_estimate.config");
                write(&config, &r.cov_estimate_config_base())?;
                let contigs = r.contigs_file()?;
                let assembled = (r.assembled_contig_file == Some(true)).then_some(contigs.as_str());
                append(
                    &config,
                    &files::cov_estimate_config_additions(
                        assembled,
                        &r.s.illumina_read_1,
                        &mapping,
                        &r.tools.samtools,
                        // TOOLS_DIR: where cov_estimate finds minibwa
                        &r.tools.minibwa,
                    ),
                )?;
                let start = Instant::now();
                let cwd = std::env::current_dir()?.display().to_string();
                crate::coverage::estimate(
                    &config,
                    &cwd,
                    crate::coverage::Programs {
                        short_read_analysis: r.worker("short_read_analysis"),
                        recorder: rec,
                    },
                )
                .context("Error in during coverage estimation (cov_estimate)")?;
                r.time_log("cov_estimation", seconds_since(start))?;
                let start = Instant::now();
                crate::edge_filter::filter_edge_coverage_ratio(&mapping, &cov_file)
                    .context("Error in during filter_edge_coverage_ratio")?;
                r.time_log("cov_filtering", seconds_since(start))?;
                Ok(())
            },
        )
    }

    // ---- stage 4 (metagenome) --------------------------------------------

    fn hierarchical_clustering(&mut self) -> Result<()> {
        let sigma = self.sigma_dir();
        let done = format!("{sigma}/r_value_evaluation_summary.dat");
        self.stage(
            "hierarchical_clustering",
            "Hierarchical clustering",
            4,
            false,
            &done,
            |r, rec| {
                init_dir(&sigma)?;
                let start = Instant::now();
                // generate_sigma_config_file
                let mapping = r.mapping_dir();
                let mapping_config = read(&format!("{mapping}/config"))?;
                for i in 0..cpp_config::EDGE_SETS {
                    let edge_dir = format!("{sigma}/pairedEdges_i{i}");
                    fs::create_dir_all(&edge_dir)
                        .with_context(|| format!("cannot create {edge_dir}"))?;
                    symlink(
                        &format!("{mapping}/pairedEdges_i{i}"),
                        &format!("{edge_dir}/pairedEdges_i{i}"),
                    )?;
                    write(
                        &format!("{edge_dir}/lib.txt"),
                        &cpp_config::sigma_lib_txt(&mapping_config, i)?,
                    )?;
                }
                let config = format!("{sigma}/sigma.config");
                write(
                    &config,
                    &cpp_config::hierarchical_sigma_config(
                        &read(&format!("{}/sigma.config", r.cov_dir()))?,
                        &sigma,
                    ),
                )?;
                run(
                    rec,
                    Pipeline::new(r.worker("sigma").arg(&config))
                        .stdout(Output::File(format!("{sigma}/sigma.out").into()))
                        .stderr(Output::File(format!("{sigma}/sigma.err").into())),
                )
                .context("Error during sigma execution")?;
                let dir = r.install.opera_ms_dir.clone();
                // refine_r.err gets the estimate's notes
                crate::cluster_eval::refine_r_estimate(
                    rec,
                    &r.cov_file(),
                    &sigma,
                    usize::try_from(r.opt.num_processor).unwrap_or(1),
                    &r.opt.cluster_len_threshold.to_string(),
                    &format!("{dir}/bin/sigma"),
                    &r.worker("sigma"),
                )
                .context("Error during the R refinement (refine_r_estimate)")?;
                r.time_log("sigma", seconds_since(start))
            },
        )
    }

    // ---- stage 5 ------------------------------------------------------------

    fn reference_clustering(&mut self) -> Result<()> {
        let ref_dir = self.ref_clustering_dir();
        let done = format!("{ref_dir}/cluster_species.dat");
        let skip = !self.opt.ref_clustering;
        self.stage(
            "reference_clustering",
            "Reference based clustering",
            5,
            skip,
            &done,
            |r, rec| {
                init_dir(&ref_dir)?;
                let start = Instant::now();
                // the reference genomes: sylph, then the cache or NCBI
                let genomes = r.reference_genomes_dir();
                r.reference_genomes(rec)?;
                crate::ref_clustering::sequence_similarity_clustering(
                    rec,
                    &crate::ref_clustering::Args {
                        inter_dir: &r.s.inter_dir,
                        ref_dir: &ref_dir,
                        mapping_dir: &r.mapping_dir(),
                        sigma_dir: &r.sigma_dir(),
                        contig_file: &r.contigs_file()?,
                        nproc: usize::try_from(r.opt.num_processor).unwrap_or(1),
                        genome_dir: &format!("{genomes}/genomes"),
                        skani_dir: &r.tools.skani,
                        minimap2_dir: &r.tools.minimap2,
                    },
                )
                .context("Error in during reference based clustering")?;
                r.time_log("ref_clustering", seconds_since(start))
            },
        )
    }

    // ---- stages 6 and 7 -----------------------------------------------------

    fn strain_clustering_and_assembly(&mut self) -> Result<()> {
        let strain = self.strain_dir();
        let both = self.opt.ref_clustering && self.opt.strain_clustering;
        let done = format!("{strain}/done");
        self.stage("strain_clustering", "Strain clustering and assembly", 6, !both, &done, |r, rec| {
            init_dir(&strain)?;
            let start = Instant::now();
            crate::coverage_clustering::coverage_clustering(&crate::coverage_clustering::Args {
                inter_dir: &r.s.inter_dir,
                out_dir: &strain,
                contig_file: &r.contigs_file()?,
                ref_dir: &r.ref_clustering_dir(),
                opera_ms_db: &format!("{}/genomes", r.reference_genomes_dir()),
            })
            .context("Error in during strain clustering")?;

            // Assembly of species with multiple strains
            let sigma = r.sigma_dir();
            let summary = read(&format!("{sigma}/r_value_evaluation_summary.dat"))?;
            let global_r = last_line(&summary).split('\t').next().unwrap_or("").to_string();
            let interval = read(&format!("{sigma}/r_estimate_value.dat"))?;
            let first = interval.split('\n').next().unwrap_or("");
            let values: Vec<f64> = first.split_whitespace().map(perl_float).collect();
            ensure!(values.len() >= 2, "{sigma}/r_estimate_value.dat: expected two values");
            let mut jobs = Vec::new();
            for line in read(&format!("{strain}/reference_length.dat"))?.split_inclusive('\n') {
                let fields: Vec<&str> = line.split('\t').collect();
                let species = fields.first().copied().unwrap_or("").to_string();
                let species_dir = format!("{strain}/{species}");
                let log = format!("{strain}/{species}.log");
                jobs.push((species, species_dir, log));
            }
            let parallel = usize::try_from(r.opt.num_processor).unwrap_or(1).max(1);
            let run_ref: &Self = r;
            let global_r = &global_r;
            let jobs: Vec<_> = jobs
                .into_iter()
                .map(|(species, species_dir, log)| {
                    move || {
                        rec.substep(&species, |sub| -> Result<()> {
                            let strains = crate::cluster_strain::cluster_strain(&species_dir, global_r, &log)
                                .with_context(|| format!("strain clustering of {species}"))?;
                            for id in strains {
                                let dir = format!("{species_dir}/STRAIN_{id}");
                                run_ref.run_opera_lg_logs(
                                    sub,
                                    &dir,
                                    &format!("{dir}/opera.config"),
                                    ["log.txt", "opera_lg.err", "log_fast.txt", "opera_lg_fast.err"],
                                )?;
                            }
                            Ok(())
                        })
                    }
                })
                .collect();
            run_parallel(jobs, parallel)
                .map_err(|e| anyhow::anyhow!("{e}"))
                .with_context(|| format!("Error in during strain clustering assembly. Please see {strain}/*.log for details"))?;
            fs::File::create(&done).with_context(|| format!("cannot create {done}"))?;
            r.time_log("strain_clustering", seconds_since(start))
        })?;

        let opera_lr = self.opera_lr_dir();
        let done = format!("{opera_lr}/statistics");
        self.stage(
            "other_clusters",
            "Assembly of other clusters",
            7,
            false,
            &done,
            |r, rec| {
                init_dir(&opera_lr)?;
                let start = Instant::now();
                let ref_dir = r.ref_clustering_dir();
                if both {
                    let mut text = read(&format!("{ref_dir}/single_strain_species.fa"))?;
                    let strain_dirs = glob(&strain, |_| true)?;
                    for d in strain_dirs {
                        let excluded = format!("{d}/excluded_contigs.fa");
                        if exists(&excluded) {
                            text.push_str(&read(&excluded)?);
                        }
                    }
                    write(
                        &format!("{ref_dir}/single_strain_species_with_excluded.fa"),
                        &text,
                    )?;
                }
                let mut edge_dir = r.sigma_dir();
                let mut cluster_file = format!("{edge_dir}/clusters");
                let mut contig_file = r.contigs_file()?;
                if r.opt.ref_clustering {
                    edge_dir = ref_dir.clone();
                    cluster_file = format!("{ref_dir}/clusters_seq_similarity");
                    if r.opt.strain_clustering {
                        contig_file = format!("{ref_dir}/single_strain_species_with_excluded.fa");
                        cluster_file = format!("{ref_dir}/clusters_single_strain");
                    }
                }
                crate::cluster_filter::filter_cluster_coverage(
                    &r.cov_file(),
                    &cluster_file,
                    &edge_dir,
                    "1.5",
                    &format!("{edge_dir}/NO_REPEAT"),
                )
                .context("Error in during filter_cluster_coverage")?;
                let mapping = r.mapping_dir();
                let config = cpp_config::long_read_opera_config(
                    &read(&format!("{opera_lr}/../lib_1_bundles/lib_1.ebconfig"))?,
                    &read(&format!("{mapping}/config"))?,
                    &edge_dir,
                    &opera_lr,
                    &contig_file,
                )?;
                write(&format!("{opera_lr}/opera.config"), &config)?;
                r.run_opera_lg(rec, &opera_lr, &format!("{opera_lr}/opera.config"))?;
                r.time_log("other_assembly", seconds_since(start))
            },
        )
    }

    /// `run_opera_lg`: OPERA-LG with the timeout; OPERA-LG-fast if
    /// OPERA-LG didn't write `scaffoldSeq.fasta` (timeout or failure).
    fn run_opera_lg(&self, rec: &StepRecorder<'_>, dir: &str, config: &str) -> Result<()> {
        self.run_opera_lg_logs(
            rec,
            dir,
            config,
            [
                "opera_lg.out",
                "opera_lg.err",
                "opera_lg_fast.out",
                "opera_lg_fast.err",
            ],
        )
    }

    /// OPERA-LG with the timeout, then OPERA-LG-fast if there is no
    /// `scaffoldSeq.fasta`; `logs`: stdout and stderr of each, in `dir`.
    fn run_opera_lg_logs(
        &self,
        rec: &StepRecorder<'_>,
        dir: &str,
        config: &str,
        logs: [&str; 4],
    ) -> Result<()> {
        let [out, err, fast_out, fast_err] = logs;
        let first = run(
            rec,
            Pipeline::new(self.worker("OPERA-LG").arg(config))
                .stdout(Output::File(format!("{dir}/{out}").into()))
                .stderr(Output::File(format!("{dir}/{err}").into()))
                .timeout(self.opt.opera_lg_timeout),
        );
        if exists(&format!("{dir}/scaffoldSeq.fasta")) {
            first.with_context(|| {
                format!("Error in during OPERA-LG assembly. Please see {dir}/{out} and {dir}/{err} for details")
            })?;
            opera_ms_core::note!(" *** OPERA-LG produced {dir}/scaffoldSeq.fasta");
            return Ok(());
        }
        match &first {
            Err(err) if err.is_timeout() => {
                opera_ms_core::note!(" *** OPERA-LG timed out: running OPERA-LG-fast")
            }
            Err(err) => opera_ms_core::note!(" *** OPERA-LG failed ({err}): running OPERA-LG-fast"),
            Ok(()) => opera_ms_core::note!(
                " *** OPERA-LG wrote no scaffoldSeq.fasta: running OPERA-LG-fast"
            ),
        }
        run(
            rec,
            Pipeline::new(self.worker("OPERA-LG-fast").arg(config))
                .stdout(Output::File(format!("{dir}/{fast_out}").into()))
                .stderr(Output::File(format!("{dir}/{fast_err}").into())),
        )
        .with_context(|| {
            format!("Error in during OPERA-LG assembly. Please see {dir}/{out}, {dir}/{err}, {dir}/{fast_out} and {dir}/{fast_err} for details")
        })?;
        opera_ms_core::note!(" *** OPERA-LG-fast produced {dir}/scaffoldSeq.fasta");
        Ok(())
    }

    // ---- stage 4 (single genome) --------------------------------------------

    fn single_genome_assembly(&mut self) -> Result<()> {
        let opera_lr = self.opera_lr_dir();
        let done = format!("{opera_lr}/scaffoldSeq.fasta");
        let skip = !self.opt.single_genome;
        self.stage(
            "single_genome",
            "Single genome assembly",
            4,
            skip,
            &done,
            |r, rec| {
                let start = Instant::now();
                let path = Path::new(&opera_lr);
                if let Ok(meta) = fs::symlink_metadata(path) {
                    if meta.is_dir() {
                        fs::remove_dir_all(path)
                    } else {
                        fs::remove_file(path)
                    }
                    .with_context(|| format!("cannot remove {opera_lr}"))?;
                }
                symlink(&r.mapping_dir(), &opera_lr)?;
                let results = format!("{opera_lr}/results");
                init_dir(&results)?;
                let cov = r.cov_file();
                opera_ms_core::note!(
                    " *** convert_sigma_to_opera_lg_coverage_file \t{cov}\t{results}/contigs"
                );
                write(
                    &format!("{results}/contigs"),
                    &files::sigma_to_opera_lg_coverage(&read(&cov)?),
                )?;
                r.run_opera_lg(rec, &opera_lr, &format!("{opera_lr}/config"))?;
                for name in ["scaffolds.scaf", "scaffoldSeq.fasta"] {
                    rename(&format!("{results}/{name}"), &format!("{opera_lr}/{name}"))?;
                }
                r.time_log("assembly", seconds_since(start))
            },
        )
    }

    // ---- stage 8 ------------------------------------------------------------

    fn gap_filling(&mut self) -> Result<()> {
        let opera_lr = self.opera_lr_dir();
        let done = format!("{opera_lr}/scaffoldSeq.fasta.filled");
        self.stage("gap_filling", "Gap filling", 8, false, &done, |r, rec| {
            let gap_dir = format!("{opera_lr}/GAPFILLING");
            init_dir(&gap_dir)?;
            let start = Instant::now();
            let both = r.opt.ref_clustering && r.opt.strain_clustering;
            let scaff = format!("{opera_lr}/scaffolds.scaf");
            let scaff_seq = format!("{opera_lr}/scaffoldSeq.fasta");
            let contig_sizes = format!("{opera_lr}/contigs");
            if both {
                for (file, single) in [
                    (&scaff, format!("{opera_lr}/scaffolds_single.scaf")),
                    (&scaff_seq, format!("{opera_lr}/scaffoldSeq_single.fasta")),
                    (&contig_sizes, format!("{opera_lr}/contigs_single")),
                ] {
                    if !exists(&single) {
                        rename(file, &single)?;
                    }
                    copy(&single, file)?;
                }
                // The multiple-strain scaffolds, in directory order
                // (file-system order, not sorted, as in OPERA-MS 0.9).
                let strain = r.strain_dir();
                let mut scaff_strain_id = 1;
                for entry in
                    fs::read_dir(&strain).with_context(|| format!("cannot list {strain}"))?
                {
                    let name = entry?.file_name().to_string_lossy().into_owned();
                    let mut strain_id = 1;
                    loop {
                        let dir = format!("{strain}/{name}/STRAIN_{strain_id}/");
                        let strain_scaff = format!("{dir}/scaffolds.scaf");
                        if !exists(&strain_scaff) {
                            break;
                        }
                        let rename_line = |line: &str| {
                            line.replacen(">opera", &format!(">strain{scaff_strain_id}_opera"), 1)
                        };
                        let text: String = read(&strain_scaff)?
                            .split_inclusive('\n')
                            .map(rename_line)
                            .collect();
                        append(&scaff, &text)?;
                        let text: String = read(&format!("{dir}/scaffoldSeq.fasta"))?
                            .split_inclusive('\n')
                            .map(rename_line)
                            .collect();
                        append(&scaff_seq, &text)?;
                        let text: String = read(&format!("{dir}/contigs"))?
                            .split_inclusive('\n')
                            .filter(|l| !l.contains("Length"))
                            .collect();
                        append(&contig_sizes, &text)?;
                        strain_id += 1;
                        scaff_strain_id += 1;
                    }
                }
            }
            if r.opt.gap_filling {
                let t = &r.tools;
                crate::gap_filling::gap_filling(
                    rec,
                    &crate::gap_filling::Args {
                        working_dir: &gap_dir,
                        contig_file: &r.contigs_file()?,
                        read_file: &r.s.long_read,
                        read_mapping_dir: &r.mapping_dir(),
                        nproc: usize::try_from(r.opt.num_processor).unwrap_or(1),
                        scaffold_length_threshold: &r.opt.filled_scaff_length.to_string(),
                        racon_dir: &t.racon,
                        minimap2_dir: &t.minimap2,
                    },
                )
                .context("Error in during gap_filling")?;
            } else {
                symlink(&scaff_seq, &done)?;
            }
            r.time_log("gap_filling", seconds_since(start))
        })
    }

    // ---- final outputs --------------------------------------------------------

    /// `generate_assembly_stats`; returns the path of `contig_info.txt`.
    fn generate_assembly_stats(&mut self) -> Result<String> {
        let out = self.s.output_dir.clone();
        let contig_info = format!("{out}/contig_info.txt");
        let bench = Arc::clone(&self.bench);
        bench.step("Assembly stats", |_rec| -> Result<()> {
            let start = Instant::now();
            let opera_lr = self.opera_lr_dir();
            let mut cluster_file = if self.opt.single_genome {
                "NULL".to_string()
            } else {
                format!("{}/clusters", self.sigma_dir())
            };
            let mut species_file = "NULL".to_string();
            let mut strain_dir = "NULL".to_string();
            if self.opt.ref_clustering {
                let ref_dir = self.ref_clustering_dir();
                cluster_file = format!("{ref_dir}/clusters_seq_similarity");
                species_file = format!("{ref_dir}/cluster_species.dat");
                strain_dir = self.strain_dir();
            }
            // `cat read_size.dat`, chop
            let mut read_size = read(&format!("{}/read_size.dat", self.cov_dir()))?;
            read_size.pop();
            crate::final_stats::generate_assembly_stats(&crate::final_stats::AssemblyStatsArgs {
                scaffold_file: &format!("{opera_lr}/scaffolds.scaf"),
                cluster_file: &cluster_file,
                read_size: &read_size,
                coverage_file: &self.cov_file(),
                long_read_coverage_file: &format!("{}/opera.map.cov", self.mapping_dir()),
                species_file: &species_file,
                strain_dir: &strain_dir,
                scaffold_seq_file: &format!("{opera_lr}/scaffoldSeq.fasta.filled"),
                outfile: &contig_info,
            })
            .context("Error in during assembly statistics generation")?;
            let stats = files::assembly_stats(&read(&contig_info)?, &perl_localtime());
            opera_ms_core::say!("{stats}");
            write(&format!("{out}/assembly.stats"), &format!("{stats}\n"))?;
            self.time_log("stats_generation", seconds_since(start))
        })?;
        Ok(contig_info)
    }

    fn write_final_assembly(&mut self, contig_info: &str) -> Result<()> {
        let opera_lr = self.opera_lr_dir();
        let text = files::final_contigs(
            &read(contig_info)?,
            &read(&format!("{opera_lr}/scaffolds.scaf.cname"))?,
            &read(&format!("{opera_lr}/scaffoldSeq.fasta.filled"))?,
        );
        write(&format!("{}/contigs.fasta", self.s.output_dir), &text)?;
        self.write_cluster_assembly(contig_info)
    }

    fn write_cluster_assembly(&mut self, contig_info: &str) -> Result<()> {
        if !self.opt.ref_clustering {
            return Ok(());
        }
        let out = self.s.output_dir.clone();
        // the cluster FASTAs: kept with --cluster-fasta, otherwise
        // only written for the cluster statistics, then removed (one
        // file per cluster; contig_info.txt gives them)
        let user_dir = format!("{out}/opera_ms_clusters/all");
        let cluster_dir = if self.opt.cluster_fasta {
            user_dir
        } else {
            if exists(&user_dir) {
                fs::remove_dir_all(&user_dir)
                    .with_context(|| format!("cannot remove {user_dir}"))?;
                let _ = fs::remove_dir(format!("{out}/opera_ms_clusters"));
            }
            self.dir("cluster_sequences")
        };
        init_dir(&cluster_dir)?;
        let polished = format!("{out}/contigs.polished.fasta");
        let final_contig = if exists(&polished) {
            polished
        } else {
            format!("{out}/contigs.fasta")
        };
        for file in files::cluster_sequences(&read(contig_info)?, &read(&final_contig)?) {
            write(&format!("{cluster_dir}/{}", file.file_name), &file.content)?;
        }
        // generate_cluster_stats
        let bench = Arc::clone(&self.bench);
        bench.step("Cluster statistics", |rec| -> Result<()> {
            let mash_dist = format!("{cluster_dir}/mash.dist");
            if !exists(&mash_dist) {
                // skani of the clusters against the run's reference
                // genomes, in the `mash.dist` layout of OPERA-MS 0.9
                let inputs = glob(&cluster_dir, |n| n.ends_with(".fasta"))?;
                let genomes = crate::references::genome_files(&self.reference_genomes_dir())?;
                let text = crate::references::skani_as_mash(
                    rec,
                    &self.tools.skani,
                    &genomes,
                    &inputs,
                    1,
                    usize::try_from(self.opt.num_processor).unwrap_or(1),
                    &format!("{cluster_dir}/clusters"),
                )
                .context("Error during skani on the cluster sequences")?;
                write(&mash_dist, &text)?;
            }
            let text = crate::final_stats::cluster_info(&out, &cluster_dir)
                .context("Error during cluster_info")?;
            write(&format!("{out}/cluster_info.txt"), &text)?;
            if !self.opt.cluster_fasta {
                for file in glob(&cluster_dir, |n| n.ends_with(".fasta"))? {
                    fs::remove_file(&file).with_context(|| format!("cannot remove {file}"))?;
                }
            }
            Ok(())
        })
    }

    // ---- stage 9 ------------------------------------------------------------

    fn polishing(&mut self, contig_info: &str) -> Result<()> {
        if self.opt.no_polishing {
            return Ok(());
        }
        let out = self.s.output_dir.clone();
        let polished = format!("{out}/contigs.polished.fasta");
        self.stage(
            "polishing",
            "Polishing",
            9,
            false,
            &polished.clone(),
            |r, rec| {
                let dir = r.dir("polished_assembly");
                init_dir(&dir)?;
                // POLCA's method, parallel; Polypolish first
                // with --polisher polypolish-pypolca
                let t = &r.tools;
                let (fasta, summary) = crate::polishing::polish(
                    rec,
                    &crate::polishing::Args {
                        dir: &dir,
                        contigs: &format!("{out}/contigs.fasta"),
                        read1: &r.s.illumina_read_1,
                        read2: &r.s.illumina_read_2,
                        nproc: usize::try_from(r.opt.num_processor).unwrap_or(1),
                        polisher: r.opt.polisher,
                        minibwa_dir: &t.minibwa,
                        samtools_dir: &t.samtools,
                        freebayes_dir: &t.freebayes,
                        polypolish_dir: &t.polypolish,
                    },
                )
                .context("Error in during polishing")?;
                opera_ms_core::note!(
                    " *** polishing ({}): {} variants called, {} corrections applied",
                    r.opt.polisher.name(),
                    summary.variants,
                    summary.applied
                );
                write(&polished, &fasta)
            },
        )?;
        self.write_cluster_assembly(contig_info)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn perl_string_helpers() {
        assert_eq!(last_line("a\tb\nc\td\n"), "c\td");
        assert_eq!(last_line("x"), "x");
        assert_eq!(perl_float("79.657812"), 79.657812);
        assert_eq!(perl_float(" 238.285988\n"), 238.285988);
        assert_eq!(perl_float("abc"), 0.0);
        assert_eq!(absolute("/cwd/", "out"), "/cwd/out");
        assert_eq!(absolute("/cwd/", "/abs"), "/abs");
    }

    #[test]
    fn glob_is_sorted_bytewise() {
        let dir = tempfile::tempdir().unwrap();
        for name in [
            "cluster_10.fasta",
            "cluster_1.fasta",
            "unclustered.fasta",
            "mash.dist",
        ] {
            fs::write(dir.path().join(name), "").unwrap();
        }
        let d = dir.path().display().to_string();
        let names: Vec<String> = glob(&d, |n| n.contains("fa"))
            .unwrap()
            .into_iter()
            .map(|p| p.rsplit('/').next().unwrap().to_string())
            .collect();
        assert_eq!(
            names,
            ["cluster_1.fasta", "cluster_10.fasta", "unclustered.fasta"]
        );
    }
}
