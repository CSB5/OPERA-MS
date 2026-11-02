//! The OPERA-MS pipeline.
//!
//! - [`run`]: the stages.
//! - [`options`]: command line.
//! - [`deps`]: the external tools, from `PATH`.
//! - [`long_read`]: long-read mapping and edges.
//! - [`coverage`]: short-read mapping and coverage estimation.
//! - [`edge_filter`]: long-read edges filtered by coverage ratio.
//! - [`cluster_eval`]: coverage clustering parameter (`r`) estimate.
//! - [`references`]: the run's reference genomes (sylph, cache, NCBI).
//! - [`ref_clustering`]: reference-based clustering.
//! - [`coverage_clustering`]: species and strain sets.
//! - [`cluster_strain`]: strain clustering.
//! - [`cluster_filter`]: repeat contigs of the clusters.
//! - [`gap_filling`]: gap filling.
//! - [`polishing`]: short-read polishing.
//! - [`final_stats`]: assembly and cluster statistics.
//! - [`files`]: the small files the pipeline writes itself.
//! - [`rstats`]: density modes and negative binomial quantiles.
//! - [`stats`]: statistics, with Perl's number formatting.

#![forbid(unsafe_code)]

pub mod cluster_eval;
pub mod cluster_filter;
pub mod cluster_strain;
pub mod coverage;
pub mod coverage_clustering;
pub mod deps;
pub mod edge_filter;
pub mod files;
pub mod final_stats;
pub mod gap_filling;
pub mod long_read;
pub mod options;
pub mod polishing;
pub mod ref_clustering;
pub mod references;
pub mod rstats;
pub mod run;
pub mod stats;
