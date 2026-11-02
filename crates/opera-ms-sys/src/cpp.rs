//! The C++ programs: sigma, OPERA-LG (and its fast
//! variant) and short_read_analysis, called through their `main()`.
//!
//! Each function runs the program **in the calling process** and returns
//! its exit code. The C++ keeps global state (`Configure::*`, `Sigma::*`),
//! may call `exit()`, reads stdin and writes to stdout/stderr, so call
//! these only from a dedicated worker process, once per process:
//! `opera-ms __worker <program> ...`, started by the runner.
//! A C++ exception becomes an `Err`.

#[cxx::bridge(namespace = "opera_ms")]
mod ffi {
    unsafe extern "C++" {
        include!("opera-ms-sys/cpp/bridge.h");

        fn run_sigma(args: &Vec<String>) -> Result<i32>;
        fn run_opera_lg(args: &Vec<String>, fast: bool) -> Result<i32>;
        fn run_short_read_analysis(args: &Vec<String>) -> Result<i32>;
    }
}

/// A C++ program run by a worker.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Program {
    /// `sigma <config>`
    Sigma,
    /// `OPERA-LG <config>`
    OperaLg,
    /// `OPERA-LG-fast <config>`
    OperaLgFast,
    /// `short_read_analysis <opera-config> <sigma-config>`, SAM on stdin
    ShortReadAnalysis,
}

impl Program {
    /// The name used on the worker command line and in logs.
    pub fn name(self) -> &'static str {
        match self {
            Program::Sigma => "sigma",
            Program::OperaLg => "OPERA-LG",
            Program::OperaLgFast => "OPERA-LG-fast",
            Program::ShortReadAnalysis => "short_read_analysis",
        }
    }

    pub fn from_name(name: &str) -> Option<Self> {
        [
            Program::Sigma,
            Program::OperaLg,
            Program::OperaLgFast,
            Program::ShortReadAnalysis,
        ]
        .into_iter()
        .find(|p| p.name() == name)
    }
}

/// Runs `program` with `args` in this process and returns its exit code.
/// See the module documentation: only from a worker process, once.
pub fn run(program: Program, args: &[String]) -> Result<i32, cxx::Exception> {
    let args: Vec<String> = args.to_vec();
    match program {
        Program::Sigma => ffi::run_sigma(&args),
        Program::OperaLg => ffi::run_opera_lg(&args, false),
        Program::OperaLgFast => ffi::run_opera_lg(&args, true),
        Program::ShortReadAnalysis => ffi::run_short_read_analysis(&args),
    }
}
