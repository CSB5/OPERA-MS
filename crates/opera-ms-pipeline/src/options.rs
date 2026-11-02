//! Command line.
//!
//! Invocation:
//! - `opera-ms` with no argument prints the help;
//! - `opera-ms check-dependency`;
//! - `opera-ms [options]`: a run.
//!
//! Options follow Perl's Getopt::Long, as OPERA-MS 0.9 did: options
//! declared with `!` also accept `--no-<name>` and `--no<name>` (the last
//! one given wins; `--no-polishing` is itself such an option, so
//! `--no-no-polishing` turns polishing back on), unique prefixes of option
//! names are accepted, and so is a single leading dash (`-out-dir`).

use std::ffi::OsString;
use std::path::PathBuf;
use std::time::Duration;

use clap::{ArgAction, Parser};
use opera_ms_core::duration::parse_timeout;

/// Software version printed in the logs.
pub const VERSION: &str = "v2.0.0";

/// All settings of a run: the defaults, then the options, then `update`.
#[derive(Debug, Clone, PartialEq)]
pub struct Options {
    pub contigs_file: Option<PathBuf>,
    pub illumina_read_1: Option<String>,
    pub illumina_read_2: Option<String>,
    pub long_read: Option<PathBuf>,
    pub output_dir: Option<PathBuf>,
    /// `--sylph-db`: the sylph database of GTDB genomes.
    pub sylph_db: Option<String>,
    /// `--genome-cache`: where reference genomes are kept between runs.
    pub genome_cache: Option<String>,
    /// `--strict-genomes`: a reference genome that can't be obtained is an
    /// error instead of a warning.
    pub strict_genomes: bool,
    /// `--prepare-genomes`: only the sylph profile and the genome downloads
    /// into `--genome-cache`, then stop (for runs on machines without
    /// internet).
    pub prepare_genomes: bool,
    /// `--offline`: never contact NCBI; a genome not in `--genome-cache` is
    /// missing.
    pub offline: bool,
    pub single_genome: bool,
    pub ref_clustering: bool,
    pub strain_clustering: bool,
    pub gap_filling: bool,
    pub no_polishing: bool,
    pub contig_len_thr: i64,
    pub contig_edge_len: i64,
    pub contig_window_len: i64,
    pub kmer_size: i64,
    pub num_processor: i64,
    /// Maximal cluster size during hierarchical clustering (-1: no limit).
    pub cluster_len_threshold: i64,
    /// Gap filling: shortest scaffold to fill.
    pub filled_scaff_length: i64,
    /// Polishing method.
    pub polisher: crate::polishing::Polisher,
    /// OPERA-LG timeout before falling back to OPERA-LG-fast;
    /// `None`: no timeout.
    pub opera_lg_timeout: Option<Duration>,
}

impl Default for Options {
    /// The default settings.
    fn default() -> Self {
        Self {
            contigs_file: None,
            illumina_read_1: None,
            illumina_read_2: None,
            long_read: None,
            output_dir: None,
            sylph_db: None,
            genome_cache: None,
            strict_genomes: false,
            prepare_genomes: false,
            offline: false,
            single_genome: false,
            ref_clustering: true,
            strain_clustering: true,
            gap_filling: true,
            no_polishing: false,
            contig_len_thr: 500,
            contig_edge_len: 80,
            contig_window_len: 340,
            kmer_size: 60,
            num_processor: 2,
            cluster_len_threshold: 6_000_000,
            filled_scaff_length: 499,
            polisher: crate::polishing::Polisher::Pypolca,
            opera_lg_timeout: Some(Duration::from_secs(300)),
        }
    }
}

impl Options {
    /// Single-genome mode turns off both
    /// clusterings.
    fn update(mut self) -> Self {
        if self.single_genome {
            self.ref_clustering = false;
            self.strain_clustering = false;
        }
        self
    }

    /// Number of stages shown in the progress lines (`[3/9]`).
    pub fn nb_total_stage(&self) -> u32 {
        let base = if self.single_genome { 6 } else { 8 };
        base + u32::from(!self.no_polishing)
    }

    /// Name of the contig coverage file (`CONTIGS_COV_FILE`).
    pub fn contigs_cov_file(&self) -> String {
        format!(
            "contigs_{}_{}",
            self.contig_window_len, self.contig_edge_len
        )
    }
}

/// What the command line asks for.
#[derive(Debug, Clone, PartialEq)]
pub enum Invocation {
    Help,
    Version,
    CheckDependency,
    Run(Box<Options>),
}

/// A command line error.
#[derive(Debug, thiserror::Error)]
pub enum OptionsError {
    #[error("{0}")]
    Clap(String),
    #[error("{0}")]
    Invalid(String),
}

/// Parses the command line (without the program name).
pub fn parse_args(args: &[OsString]) -> Result<Invocation, OptionsError> {
    let first = args.first().and_then(|a| a.to_str());
    match first {
        None => return Ok(Invocation::Help),
        Some("check-dependency") => return Ok(Invocation::CheckDependency),
        Some("--version" | "-version" | "-V") => return Ok(Invocation::Version),
        _ => {}
    }
    parse_options(args)
}

/// The command-line options.
#[derive(Debug, Parser)]
#[command(
    name = "opera-ms",
    disable_help_flag = true,
    disable_version_flag = true,
    infer_long_args = true,
    no_binary_name = true
)]
struct Cli {
    #[arg(long = "contig-file")]
    contig_file: Option<PathBuf>,
    #[arg(long = "short-read1")]
    short_read1: Option<String>,
    #[arg(long = "short-read2")]
    short_read2: Option<String>,
    #[arg(long = "long-read")]
    long_read: Option<PathBuf>,
    #[arg(long = "out-dir")]
    out_dir: Option<PathBuf>,
    #[arg(long = "sylph-db")]
    sylph_db: Option<String>,
    #[arg(long = "genome-cache")]
    genome_cache: Option<String>,
    #[arg(long = "strict-genomes", action = ArgAction::SetTrue)]
    strict_genomes: bool,
    #[arg(long = "prepare-genomes", action = ArgAction::SetTrue)]
    prepare_genomes: bool,
    #[arg(long = "offline", action = ArgAction::SetTrue)]
    offline: bool,

    #[arg(long = "single-genome", action = ArgAction::SetTrue, overrides_with = "no_single_genome")]
    single_genome: bool,
    #[arg(long = "no-single-genome", alias = "nosingle-genome", action = ArgAction::SetTrue, overrides_with = "single_genome")]
    no_single_genome: bool,
    #[arg(long = "strain-clustering", action = ArgAction::SetTrue, overrides_with = "no_strain_clustering")]
    strain_clustering: bool,
    #[arg(long = "no-strain-clustering", alias = "nostrain-clustering", action = ArgAction::SetTrue, overrides_with = "strain_clustering")]
    no_strain_clustering: bool,
    #[arg(long = "ref-clustering", action = ArgAction::SetTrue, overrides_with = "no_ref_clustering")]
    ref_clustering: bool,
    #[arg(long = "no-ref-clustering", alias = "noref-clustering", action = ArgAction::SetTrue, overrides_with = "ref_clustering")]
    no_ref_clustering: bool,
    #[arg(long = "gap-filling", action = ArgAction::SetTrue, overrides_with = "no_gap_filling")]
    gap_filling: bool,
    #[arg(long = "no-gap-filling", alias = "nogap-filling", action = ArgAction::SetTrue, overrides_with = "gap_filling")]
    no_gap_filling: bool,
    #[arg(long = "no-polishing", action = ArgAction::SetTrue, overrides_with = "no_no_polishing")]
    no_polishing: bool,
    #[arg(long = "no-no-polishing", alias = "nono-polishing", hide = true, action = ArgAction::SetTrue, overrides_with = "no_polishing")]
    no_no_polishing: bool,

    #[arg(long = "contig-len-thr")]
    contig_len_thr: Option<i64>,
    #[arg(long = "contig-edge-len")]
    contig_edge_len: Option<i64>,
    #[arg(long = "contig-window-len")]
    contig_window_len: Option<i64>,
    #[arg(long = "kmer-size")]
    kmer_size: Option<i64>,
    #[arg(long = "num-processors")]
    num_processors: Option<i64>,
    #[arg(long = "cluster-len", hide = true)]
    cluster_len: Option<i64>,
    #[arg(long = "opera-lg-timeout")]
    opera_lg_timeout: Option<String>,
    #[arg(long = "polisher")]
    polisher: Option<String>,

    #[arg(long = "help", action = ArgAction::SetTrue)]
    help: bool,
}

/// Turns Perl-style `-name` / `-name=value` into `--name...`, for the
/// long options only (OPERA-MS has no short options).
fn double_dash(args: &[OsString]) -> Vec<OsString> {
    args.iter()
        .map(|arg| match arg.to_str() {
            Some(a) if a.len() > 2 && a.starts_with('-') && !a.starts_with("--") => {
                let name = a[1..].split('=').next().unwrap_or("");
                if name.chars().next().is_some_and(|c| c.is_ascii_alphabetic()) {
                    OsString::from(format!("-{a}"))
                } else {
                    arg.clone()
                }
            }
            _ => arg.clone(),
        })
        .collect()
}

/// Resolves a `x!` option pair: the last one given wins, else `default`.
fn negatable(yes: bool, no: bool, default: bool) -> bool {
    if no {
        false
    } else if yes {
        true
    } else {
        default
    }
}

fn parse_options(args: &[OsString]) -> Result<Invocation, OptionsError> {
    let cli = Cli::try_parse_from(double_dash(args))
        .map_err(|e| OptionsError::Clap(format!("Error in command line arguments.\n{e}")))?;
    if cli.help {
        return Ok(Invocation::Help);
    }
    let defaults = Options::default();
    let options = Options {
        contigs_file: cli.contig_file,
        illumina_read_1: cli.short_read1,
        illumina_read_2: cli.short_read2,
        long_read: cli.long_read,
        output_dir: cli.out_dir,
        sylph_db: cli.sylph_db,
        genome_cache: cli.genome_cache,
        strict_genomes: cli.strict_genomes,
        prepare_genomes: cli.prepare_genomes,
        offline: cli.offline,
        single_genome: negatable(
            cli.single_genome,
            cli.no_single_genome,
            defaults.single_genome,
        ),
        ref_clustering: negatable(
            cli.ref_clustering,
            cli.no_ref_clustering,
            defaults.ref_clustering,
        ),
        strain_clustering: negatable(
            cli.strain_clustering,
            cli.no_strain_clustering,
            defaults.strain_clustering,
        ),
        gap_filling: negatable(cli.gap_filling, cli.no_gap_filling, defaults.gap_filling),
        no_polishing: negatable(cli.no_polishing, cli.no_no_polishing, defaults.no_polishing),
        contig_len_thr: cli.contig_len_thr.unwrap_or(defaults.contig_len_thr),
        contig_edge_len: cli.contig_edge_len.unwrap_or(defaults.contig_edge_len),
        contig_window_len: cli.contig_window_len.unwrap_or(defaults.contig_window_len),
        kmer_size: cli.kmer_size.unwrap_or(defaults.kmer_size),
        num_processor: cli.num_processors.unwrap_or(defaults.num_processor),
        cluster_len_threshold: cli.cluster_len.unwrap_or(defaults.cluster_len_threshold),
        polisher: match cli.polisher.as_deref() {
            None | Some("pypolca") => crate::polishing::Polisher::Pypolca,
            Some("polypolish-pypolca") => crate::polishing::Polisher::PolypolishPypolca,
            Some(other) => {
                return Err(OptionsError::Invalid(format!(
                    "--polisher: {other}: use pypolca or polypolish-pypolca"
                )));
            }
        },
        opera_lg_timeout: match cli.opera_lg_timeout {
            Some(t) => parse_timeout(&t)
                .map_err(|e| OptionsError::Invalid(format!("--opera-lg-timeout: {e}")))?,
            None => defaults.opera_lg_timeout,
        },
        ..defaults
    };
    let options = options.update();
    check_genome_options(&options)?;
    Ok(Invocation::Run(Box::new(options)))
}

/// `--prepare-genomes` and `--offline` need a genome cache (the genomes
/// prepared for, or used by, another run), and exclude each other.
fn check_genome_options(o: &Options) -> Result<(), OptionsError> {
    let invalid = |m: &str| Err(OptionsError::Invalid(m.to_string()));
    if o.prepare_genomes && o.offline {
        return invalid("--prepare-genomes downloads genomes: it can't be --offline");
    }
    if o.prepare_genomes && !o.ref_clustering {
        return invalid(
            "--prepare-genomes prepares the reference clustering: \
             not with --no-ref-clustering or --single-genome",
        );
    }
    if (o.prepare_genomes || o.offline) && o.genome_cache.is_none() {
        let option = if o.offline {
            "--offline"
        } else {
            "--prepare-genomes"
        };
        return Err(OptionsError::Invalid(format!(
            "{option} needs --genome-cache <dir> (where the genomes are kept for the runs)"
        )));
    }
    Ok(())
}

/// The help text.
pub fn help() -> String {
    format!(
        "opera-ms: OPERA-MS {VERSION} (Rust)
contacts: Jean-Sebastien Gounot <jean_sebastien@a-star.edu.sg>

Usage:
  opera-ms [options] --short-read1 <pe1> --short-read2 <pe2> --long-read <lr> --out-dir <out_dir>
  opera-ms --prepare-genomes --short-read1 <pe1> --short-read2 <pe2> --sylph-db <db> --genome-cache <dir> --out-dir <out_dir>
  opera-ms check-dependency

Required arguments:

      --short-read1          STR   fasta file of illumina read1 <pe1>
      --short-read2          STR   fasta file of illumina read2 <pe2>
      --long-read            STR   fasta/fastq file of long reads <lr>, plain or gzip-compressed
      --out-dir              STR   output directory for scaffolding results <out_dir>

Optional arguments:

    Algorithm options:
      --no-ref-clustering          disable reference level clustering
      --no-strain-clustering       disable strain level clustering
      --no-gap-filling             disable gap-filling stage
      --no-polishing               disable assembly polishing
      --polisher             STR   short-read polishing: pypolca (POLCA's method) or polypolish-pypolca [pypolca]
      --kmer-size            INT   kmer value used to assemble contigs [60]
      --contig-len-thr       INT   contig length threshold for clustering; contigs smaller than the threshold will be filtered out [500]
      --contig-edge-len      INT   during contig coverage calculation, number of bases filtered out from each contig end, to avoid biases due to lower mapping efficiency [80]
      --contig-window-len    INT   window length in which the coverage estimation is performed. We recommend using contig-len-thr - 2 * contig-edge-len as the value [340]
      --opera-lg-timeout     DUR   time given to OPERA-LG before falling back to OPERA-LG-fast, e.g. 300, 5m, 1h; 0 = no timeout [5m]

   Other arguments:
      --sylph-db             STR   sylph database of GTDB genomes (.syl2db), needed for reference level clustering
      --genome-cache         STR   folder keeping the reference genomes between runs (looked up first, downloads added)
      --strict-genomes             stop when a reference genome can't be obtained (default: warning, genome skipped)
      --prepare-genomes            only find the sample's reference genomes (sylph) and download them into --genome-cache, then stop
                                   (prepares a run on a machine without internet: run it there with the same --genome-cache and --offline)
      --offline                    never download genomes: a genome missing from --genome-cache is skipped (or an error with --strict-genomes)
      --contig-file          STR   path to the contig file, if the short-reads have been assembled previously [default assembly using MEGAHIT]
      --num-processors       INT   number of processors to use (note that 2 is the minimum) [2]
      --version                    print the version

"
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    fn args(list: &[&str]) -> Vec<OsString> {
        list.iter().map(OsString::from).collect()
    }

    fn run(list: &[&str]) -> Options {
        match parse_args(&args(list)).unwrap() {
            Invocation::Run(options) => *options,
            other => panic!("unexpected {other:?}"),
        }
    }

    #[test]
    fn no_argument_is_help() {
        assert_eq!(parse_args(&[]).unwrap(), Invocation::Help);
        assert_eq!(parse_args(&args(&["--help"])).unwrap(), Invocation::Help);
    }

    #[test]
    fn genome_preparation_options() {
        let base = [
            "--short-read1",
            "R1.fq.gz",
            "--short-read2",
            "R2.fq.gz",
            "--sylph-db",
            "db.syl2db",
            "--out-dir",
            "out",
        ];
        let with = |extra: &[&str]| {
            let mut list = base.to_vec();
            list.extend_from_slice(extra);
            parse_args(&args(&list))
        };
        let o = run(&[&base[..], &["--prepare-genomes", "--genome-cache", "c"]].concat());
        assert!(o.prepare_genomes && !o.offline);
        let o = run(&[&base[..], &["--offline", "--genome-cache", "c"]].concat());
        assert!(o.offline && !o.prepare_genomes);
        for bad in [
            &["--prepare-genomes"][..],
            &["--offline"],
            &["--prepare-genomes", "--offline", "--genome-cache", "c"],
            &[
                "--prepare-genomes",
                "--genome-cache",
                "c",
                "--no-ref-clustering",
            ],
            &[
                "--prepare-genomes",
                "--genome-cache",
                "c",
                "--single-genome",
            ],
        ] {
            assert!(
                matches!(with(bad), Err(OptionsError::Invalid(_))),
                "{bad:?}"
            );
        }
    }

    #[test]
    fn subcommands() {
        assert_eq!(
            parse_args(&args(&["check-dependency"])).unwrap(),
            Invocation::CheckDependency
        );
    }

    #[test]
    fn readme_command_line() {
        let o = run(&[
            "--contig-file",
            "c.fa",
            "--short-read1",
            "R1.fq.gz",
            "--short-read2",
            "R2.fq.gz",
            "--long-read",
            "lr.fq",
            "--no-ref-clustering",
            "--out-dir",
            "out",
        ]);
        assert_eq!(o.contigs_file.as_deref(), Some(Path::new("c.fa")));
        assert_eq!(o.illumina_read_1.as_deref(), Some("R1.fq.gz"));
        assert!(!o.ref_clustering);
        assert!(o.strain_clustering);
        assert_eq!(o.output_dir.as_deref(), Some(Path::new("out")));
        assert_eq!(o.nb_total_stage(), 9);
        assert_eq!(o.opera_lg_timeout, Some(Duration::from_secs(300)));
    }

    #[test]
    fn defaults_match_set_default_value() {
        let o = run(&["--out-dir", "o"]);
        assert_eq!((o.kmer_size, o.num_processor), (60, 2));
        assert_eq!(
            (o.contig_edge_len, o.contig_window_len, o.contig_len_thr),
            (80, 340, 500)
        );
        assert_eq!(o.contigs_cov_file(), "contigs_340_80");
        assert_eq!(o.cluster_len_threshold, 6_000_000);
        assert!(o.ref_clustering && o.strain_clustering && o.gap_filling && !o.no_polishing);
        assert_eq!(o.sylph_db, None);
        assert!(!o.strict_genomes);
    }

    #[test]
    fn perl_negation_semantics() {
        assert!(!run(&["--no-gap-filling"]).gap_filling);
        assert!(!run(&["--nogap-filling"]).gap_filling);
        // last one wins
        assert!(run(&["--no-gap-filling", "--gap-filling"]).gap_filling);
        assert!(!run(&["--gap-filling", "--no-gap-filling"]).gap_filling);
        assert!(run(&["--no-polishing"]).no_polishing);
        assert!(!run(&["--no-polishing", "--no-no-polishing"]).no_polishing);
        assert!(!run(&["--nono-polishing"]).no_polishing);
        assert_eq!(run(&["--no-polishing"]).nb_total_stage(), 8);
    }

    #[test]
    fn single_genome_turns_off_clustering() {
        let o = run(&["--single-genome", "--ref-clustering"]);
        assert!(o.single_genome && !o.ref_clustering && !o.strain_clustering);
        assert_eq!(o.nb_total_stage(), 7);
    }

    #[test]
    fn abbreviations_and_single_dash() {
        let o = run(&["-out-dir", "o", "--num-proc", "8", "-kmer-size=31"]);
        assert_eq!(o.output_dir.as_deref(), Some(Path::new("o")));
        assert_eq!((o.num_processor, o.kmer_size), (8, 31));
        // ambiguous prefix
        assert!(parse_args(&args(&["--contig", "5"])).is_err());
    }

    #[test]
    fn invalid_values_are_errors() {
        assert!(parse_args(&args(&["--kmer-size", "abc"])).is_err());
        // options of OPERA-MS 0.9 that no longer exist
        assert!(parse_args(&args(&["--long-read-mapper", "blasr"])).is_err());
        assert!(parse_args(&args(&["--short-read-assembler", "spades"])).is_err());
        assert!(parse_args(&args(&["some.config"])).is_err());
        assert!(parse_args(&args(&["--opera-lg-timeout", "5x"])).is_err());
        assert!(parse_args(&args(&["--unknown-option"])).is_err());
    }

    #[test]
    fn opera_lg_timeout() {
        assert_eq!(
            run(&["--opera-lg-timeout", "1h"]).opera_lg_timeout,
            Some(Duration::from_secs(3600))
        );
        assert_eq!(run(&["--opera-lg-timeout", "0"]).opera_lg_timeout, None);
    }
}
