//! The only OPERA-MS crate allowed to contain `unsafe` code. Every other
//! crate has `#![forbid(unsafe_code)]` and uses the safe wrappers defined
//! here.
//!
//! - [`process`]: libc calls the standard library doesn't expose
//!   (`waitid` without reaping, `wait4` with resource usage, signals to
//!   process groups, `getrusage`).
//! - [`fs`]: `statvfs` (free disk space).
//! - [`cpp`]: the C++ bridge to sigma, OPERA-LG and short_read_analysis
//!   (cxx), to be called from worker processes only.

#![deny(unsafe_op_in_unsafe_fn)]

pub mod cpp;
pub mod fs;
pub mod process;
