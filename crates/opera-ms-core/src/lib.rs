//! Pipeline machinery shared by the OPERA-MS steps.
//!
//! - [`runner`]: external commands and pipelines, with every failure
//!   reported.
//! - [`bench`]: per-step benchmark, `benchmark.tsv`.
//! - [`duration`]: the `--opera-lg-timeout` format.
//! - [`parallel`]: independent jobs, at most N at a time (`xargs -P`).
//! - [`stage`]: progress lines and resume.
//! - [`log`]: the main log, `say!` / `note!`.

#![forbid(unsafe_code)]

pub mod bench;
pub mod duration;
pub mod log;
pub mod parallel;
pub mod runner;
pub mod stage;
