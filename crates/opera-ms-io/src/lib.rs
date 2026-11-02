//! File formats used between OPERA-MS steps.
//!
//! - [`seq`]: FASTA / FASTQ records, plain or gzip-compressed.
//! - [`cpp_config`]: sigma / OPERA-LG configuration files.

#![forbid(unsafe_code)]

pub mod cpp_config;
pub mod seq;
