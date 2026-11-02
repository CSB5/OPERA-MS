//! `opera-ms`: the OPERA-MS command line.

#![forbid(unsafe_code)]

use std::ffi::OsString;
use std::process::ExitCode;

use anyhow::{Result, bail};
use opera_ms_pipeline::options::{Invocation, Options, OptionsError, VERSION, help, parse_args};
use opera_ms_pipeline::run::{Install, Run};

/// Exit code after Ctrl-C / SIGTERM (128 + SIGINT, as shells do).
const INTERRUPTED: u8 = 130;

/// Hidden subcommand running one C++ program in this process:
/// `opera-ms __worker <program> <args...>`, started by the pipeline.
const WORKER: &str = "__worker";

fn main() -> ExitCode {
    let args: Vec<OsString> = std::env::args_os().skip(1).collect();
    if args.first().is_some_and(|a| a == WORKER) {
        return worker(&args[1..]);
    }

    tracing_subscriber::fmt()
        .with_writer(std::io::stderr)
        .with_max_level(tracing::Level::WARN)
        .with_target(false)
        .init();

    // Child processes run in their own process groups (see the runner), so
    // the terminal's Ctrl-C doesn't reach them: kill them ourselves.
    if let Err(err) = ctrlc::set_handler(|| {
        opera_ms_core::runner::kill_all_running();
        eprintln!("\nopera-ms: interrupted, child processes killed");
        std::process::exit(i32::from(INTERRUPTED));
    }) {
        eprintln!("opera-ms: warning: cannot install the Ctrl-C handler: {err}");
    }

    match run(&args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(err) => {
            // also in the main log, if the run got that far
            opera_ms_core::note!("opera-ms: error: {err:#}");
            ExitCode::FAILURE
        }
    }
}

/// Runs sigma, OPERA-LG(-fast) or short_read_analysis in this process and
/// exits with its exit code. The C++ may also call `exit()` itself.
fn worker(args: &[OsString]) -> ExitCode {
    use opera_ms_sys::cpp::{self, Program};

    let mut args = args.iter().map(|a| a.to_string_lossy().into_owned());
    let name = args.next().unwrap_or_default();
    let Some(program) = Program::from_name(&name) else {
        eprintln!("opera-ms {WORKER}: unknown program {name:?}");
        return ExitCode::FAILURE;
    };
    let args: Vec<String> = args.collect();
    match cpp::run(program, &args) {
        // exit codes are 0..=255 on Unix; C++ `return -1` means 255
        Ok(code) => ExitCode::from((code & 0xff) as u8),
        Err(err) => {
            eprintln!("opera-ms {WORKER}: {name} threw an exception: {err}");
            ExitCode::FAILURE
        }
    }
}

fn run(args: &[OsString]) -> Result<()> {
    let invocation = match parse_args(args) {
        Ok(invocation) => invocation,
        Err(OptionsError::Clap(message)) => bail!("{message}\n{}", help()),
        Err(err) => return Err(err.into()),
    };
    match invocation {
        Invocation::Help => {
            print!("{}", help());
            Ok(())
        }
        Invocation::Version => {
            println!("opera-ms {VERSION}");
            Ok(())
        }
        Invocation::CheckDependency => match opera_ms_pipeline::deps::check_dependency() {
            Ok(report) => {
                print!("{report}");
                Ok(())
            }
            Err(report) => bail!("{report}"),
        },
        Invocation::Run(options) => run_pipeline(&options),
    }
}

fn run_pipeline(options: &Options) -> Result<()> {
    let install = Install::of_current_exe()?;
    Run::prepare(options.clone(), install)?.execute()
}
