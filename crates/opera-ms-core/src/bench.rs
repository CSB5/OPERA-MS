//! Per-step benchmark of a run: `<out-dir>/benchmark.tsv`.
//!
//! For every step: wall time, CPU time (user and system, plus CPU-hours)
//! and peak memory, measured by `opera-ms` itself:
//! - CPU time of external processes is exact: `wait4` usage of every
//!   process the step ran (see [`crate::runner`]), plus the CPU time of
//!   `opera-ms` itself during the step (`getrusage`);
//! - `peak_rss_single_mb` is the largest single process (exact, `wait4`);
//! - `peak_rss_total_mb` is the peak of the memory of `opera-ms` and all
//!   its descendants together, sampled from `/proc` (every 0.5 s by
//!   default, so a very short peak can be missed).
//!
//! A row is written (and flushed) as soon as its step ends, so a failed
//! or killed run still shows how far it got and what it used. Top-level
//! steps run one at a time; sub-steps may run in parallel and get their
//! own rows, with `NA` as total memory (the processes running in parallel
//! can't be told apart), and aren't counted again in the `total` row.

use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::thread::{self, JoinHandle};
use std::time::{Duration, Instant, SystemTime};

use chrono::{DateTime, Local};
use opera_ms_sys::process::{self as sys, Usage};

use crate::runner::RunReport;

const HEADER: &str = "step\tstatus\tstart\tend\twall_s\tuser_s\tsys_s\tcpu_hours\t\
                      peak_rss_total_mb\tpeak_rss_single_mb\tn_processes";

/// Outcome of a step, as written in the `status` column.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Status {
    Ok,
    Failed,
    Skipped,
    Timeout,
}

impl Status {
    fn as_str(self) -> &'static str {
        match self {
            Status::Ok => "ok",
            Status::Failed => "failed",
            Status::Skipped => "skipped",
            Status::Timeout => "timeout",
        }
    }
}

#[derive(Debug, Default, Clone, Copy)]
struct Totals {
    user: Duration,
    system: Duration,
    max_rss_kb: u64,
    processes: usize,
}

impl Totals {
    fn add_usage(&mut self, usage: &Usage) {
        self.user += usage.user;
        self.system += usage.system;
        self.max_rss_kb = self.max_rss_kb.max(usage.max_rss_kb);
    }
}

/// Collects the processes run by one step (or sub-step).
pub struct StepRecorder<'a> {
    bench: &'a Benchmark,
    name: String,
    totals: Mutex<Totals>,
    parent: Option<&'a StepRecorder<'a>>,
}

impl<'a> StepRecorder<'a> {
    /// Adds the processes of a run (successful or not) to this step.
    pub fn record(&self, report: &RunReport) {
        self.add(report);
        if let Some(parent) = self.parent {
            parent.add(report);
        }
    }

    fn add(&self, report: &RunReport) {
        if let Ok(mut totals) = self.totals.lock() {
            for process in &report.processes {
                totals.add_usage(&process.usage);
                totals.processes += 1;
            }
        }
    }

    /// Runs a sub-step (e.g. one strain cluster's OPERA-LG run) with its own
    /// row named `<step>/<name>`. Sub-steps may run in parallel.
    pub fn substep<T, E>(
        &self,
        name: &str,
        run: impl FnOnce(&StepRecorder<'_>) -> Result<T, E>,
    ) -> Result<T, E> {
        let sub = StepRecorder {
            bench: self.bench,
            name: format!("{}/{}", self.name, name),
            totals: Mutex::new(Totals::default()),
            parent: Some(self),
        };
        let started = Started::now(false);
        let result = run(&sub);
        let status = if result.is_ok() {
            Status::Ok
        } else {
            Status::Failed
        };
        self.bench
            .write_row(&sub.name, status, &started, sub.totals(), None);
        result
    }

    fn totals(&self) -> Totals {
        self.totals.lock().map(|t| *t).unwrap_or_default()
    }
}

struct Started {
    wall_clock: SystemTime,
    instant: Instant,
    /// CPU of opera-ms itself at the start (top-level steps only).
    self_usage: Option<Usage>,
}

impl Started {
    fn now(with_self_usage: bool) -> Self {
        Self {
            wall_clock: SystemTime::now(),
            instant: Instant::now(),
            self_usage: if with_self_usage {
                sys::self_usage().ok()
            } else {
                None
            },
        }
    }
}

/// Peak memory of opera-ms and its descendants, sampled from /proc.
struct Sampler {
    stop: Arc<AtomicBool>,
    step_peak_kb: Arc<AtomicU64>,
    run_peak_kb: Arc<AtomicU64>,
    thread: Option<JoinHandle<()>>,
}

impl Sampler {
    fn start(interval: Duration) -> Self {
        let stop = Arc::new(AtomicBool::new(false));
        let step_peak_kb = Arc::new(AtomicU64::new(0));
        let run_peak_kb = Arc::new(AtomicU64::new(0));
        let thread = {
            let (stop, step_peak_kb, run_peak_kb) =
                (stop.clone(), step_peak_kb.clone(), run_peak_kb.clone());
            thread::Builder::new()
                .name("memory-sampler".into())
                .spawn(move || {
                    while !stop.load(Ordering::Relaxed) {
                        let total = process_tree_rss_kb(std::process::id());
                        step_peak_kb.fetch_max(total, Ordering::Relaxed);
                        run_peak_kb.fetch_max(total, Ordering::Relaxed);
                        thread::sleep(interval);
                    }
                })
                .ok()
        };
        Self {
            stop,
            step_peak_kb,
            run_peak_kb,
            thread,
        }
    }
}

impl Drop for Sampler {
    fn drop(&mut self) {
        self.stop.store(true, Ordering::Relaxed);
        if let Some(thread) = self.thread.take() {
            // The sampler only reads /proc; a panic there loses no data.
            let _ignored = thread.join();
        }
    }
}

/// Sum of the resident memory (kB) of `root` and all its descendants.
fn process_tree_rss_kb(root: u32) -> u64 {
    let Ok(entries) = fs::read_dir("/proc") else {
        return 0;
    };
    let mut parent_of: Vec<(u32, u32)> = Vec::new();
    for entry in entries.flatten() {
        let Some(pid) = entry
            .file_name()
            .to_str()
            .and_then(|n| n.parse::<u32>().ok())
        else {
            continue;
        };
        // "pid (comm) state ppid ..."; comm may contain spaces or ')'.
        let Ok(stat) = fs::read_to_string(entry.path().join("stat")) else {
            continue;
        };
        let Some(after_comm) = stat.rfind(')').map(|i| &stat[i + 1..]) else {
            continue;
        };
        if let Some(ppid) = after_comm
            .split_whitespace()
            .nth(1)
            .and_then(|p| p.parse().ok())
        {
            parent_of.push((pid, ppid));
        }
    }
    let mut tree = vec![root];
    let mut index = 0;
    while index < tree.len() {
        let parent = tree[index];
        tree.extend(
            parent_of
                .iter()
                .filter(|(_, ppid)| *ppid == parent)
                .map(|(pid, _)| *pid),
        );
        index += 1;
    }
    tree.iter().map(|pid| vm_rss_kb(*pid)).sum()
}

fn vm_rss_kb(pid: u32) -> u64 {
    fs::read_to_string(format!("/proc/{pid}/status"))
        .ok()
        .and_then(|status| {
            status
                .lines()
                .find_map(|l| l.strip_prefix("VmRSS:"))
                .and_then(|v| v.split_whitespace().next().and_then(|n| n.parse().ok()))
        })
        .unwrap_or(0)
}

/// The benchmark file of a run.
pub struct Benchmark {
    path: PathBuf,
    out: Mutex<BufWriter<File>>,
    sampler: Sampler,
    started: Started,
    run_totals: Mutex<Totals>,
}

impl Benchmark {
    /// Creates (truncates) the benchmark file and starts the memory sampler.
    pub fn create(path: impl AsRef<Path>) -> io::Result<Self> {
        Self::with_interval(path, Duration::from_millis(500))
    }

    /// As [`Benchmark::create`], with another memory sampling interval.
    pub fn with_interval(path: impl AsRef<Path>, interval: Duration) -> io::Result<Self> {
        let path = path.as_ref().to_path_buf();
        let mut out = BufWriter::new(File::create(&path)?);
        writeln!(out, "{HEADER}")?;
        out.flush()?;
        Ok(Self {
            path,
            out: Mutex::new(out),
            sampler: Sampler::start(interval),
            started: Started::now(true),
            run_totals: Mutex::new(Totals::default()),
        })
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    /// Runs a top-level step and writes its row: `ok` if `run` returns
    /// `Ok`, `failed` otherwise. Steps must not run in parallel.
    pub fn step<T, E>(
        &self,
        name: &str,
        run: impl FnOnce(&StepRecorder<'_>) -> Result<T, E>,
    ) -> Result<T, E> {
        self.step_with_status(name, |rec| {
            let result = run(rec);
            let status = if result.is_ok() {
                Status::Ok
            } else {
                Status::Failed
            };
            (result, status)
        })
    }

    /// As [`Benchmark::step`], with the status chosen by `run`.
    pub fn step_with_status<R>(
        &self,
        name: &str,
        run: impl FnOnce(&StepRecorder<'_>) -> (R, Status),
    ) -> R {
        let recorder = StepRecorder {
            bench: self,
            name: name.to_string(),
            totals: Mutex::new(Totals::default()),
            parent: None,
        };
        self.sampler.step_peak_kb.store(0, Ordering::Relaxed);
        let started = Started::now(true);
        let (result, status) = run(&recorder);
        let mut totals = recorder.totals();
        if let (Some(before), Ok(after)) = (started.self_usage, sys::self_usage()) {
            totals.user += after.user.saturating_sub(before.user);
            totals.system += after.system.saturating_sub(before.system);
        }
        // Sample once more, so that a step shorter than the interval has a value.
        let now = process_tree_rss_kb(std::process::id());
        let peak = self.sampler.step_peak_kb.load(Ordering::Relaxed).max(now);
        self.write_row(name, status, &started, totals, Some(peak));
        if let Ok(mut run_totals) = self.run_totals.lock() {
            run_totals.user += totals.user;
            run_totals.system += totals.system;
            run_totals.max_rss_kb = run_totals.max_rss_kb.max(totals.max_rss_kb);
            run_totals.processes += totals.processes;
        }
        result
    }

    /// Writes the row of a step skipped because a previous run completed it.
    pub fn skipped(&self, name: &str) {
        let now = Started::now(false);
        self.write_line(&format!(
            "{name}\t{}\t{}\t{}\tNA\tNA\tNA\tNA\tNA\tNA\t0",
            Status::Skipped.as_str(),
            timestamp(now.wall_clock),
            timestamp(now.wall_clock)
        ));
    }

    /// Writes the `total` row (call once, at the end of the run).
    pub fn finish(&self, status: Status) {
        let totals = self.run_totals.lock().map(|t| *t).unwrap_or_default();
        let peak = self.sampler.run_peak_kb.load(Ordering::Relaxed);
        self.write_row("total", status, &self.started, totals, Some(peak));
    }

    fn write_row(
        &self,
        name: &str,
        status: Status,
        started: &Started,
        totals: Totals,
        peak_total_kb: Option<u64>,
    ) {
        let wall = started.instant.elapsed();
        let end = started.wall_clock + wall;
        let cpu = totals.user + totals.system;
        let peak_total = peak_total_kb.map_or_else(|| "NA".to_string(), mb);
        self.write_line(&format!(
            "{name}\t{}\t{}\t{}\t{:.1}\t{:.1}\t{:.1}\t{:.4}\t{peak_total}\t{}\t{}",
            status.as_str(),
            timestamp(started.wall_clock),
            timestamp(end),
            wall.as_secs_f64(),
            totals.user.as_secs_f64(),
            totals.system.as_secs_f64(),
            cpu.as_secs_f64() / 3600.0,
            mb(totals.max_rss_kb),
            totals.processes,
        ));
    }

    fn write_line(&self, line: &str) {
        let result = self
            .out
            .lock()
            .map(|mut out| writeln!(out, "{line}").and_then(|()| out.flush()));
        // The benchmark must never stop a run: report and go on.
        match result {
            Ok(Ok(())) => {}
            Ok(Err(err)) => tracing::warn!("cannot write {}: {err}", self.path.display()),
            Err(_) => tracing::warn!("cannot write {}: lock poisoned", self.path.display()),
        }
    }
}

fn mb(kb: u64) -> String {
    format!("{}", (kb + 512) / 1024)
}

fn timestamp(time: SystemTime) -> String {
    DateTime::<Local>::from(time)
        .format("%Y-%m-%dT%H:%M:%S%.3f")
        .to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::runner::{Cmd, Pipeline};

    fn rows(path: &Path) -> Vec<Vec<String>> {
        fs::read_to_string(path)
            .unwrap()
            .lines()
            .map(|l| l.split('\t').map(String::from).collect())
            .collect()
    }

    #[test]
    fn steps_record_cpu_memory_and_status() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("benchmark.tsv");
        let bench = Benchmark::with_interval(&path, Duration::from_millis(50)).unwrap();

        // ~200 MB held by `tail` (no newline in /dev/zero), plus some CPU
        let result: Result<(), String> = bench.step("memory", |rec| {
            let report = Pipeline::new(Cmd::new("head").args(["-c", "200000000", "/dev/zero"]))
                .pipe(Cmd::new("tail").arg("-n1"))
                .pipe(Cmd::new("wc").arg("-c"))
                .run()
                .map_err(|e| e.to_string())?;
            rec.record(&report);
            Ok(())
        });
        result.unwrap();
        let failed: Result<(), String> = bench.step("broken", |rec| {
            let err = Pipeline::new(Cmd::new("false")).run().unwrap_err();
            rec.record(err.report().unwrap());
            Err(err.to_string())
        });
        assert!(failed.is_err());
        bench.skipped("done before");
        bench.finish(Status::Failed);

        let rows = rows(&path);
        assert_eq!(rows[0].join("\t"), HEADER);
        let (memory, broken, skipped, total) = (&rows[1], &rows[2], &rows[3], &rows[4]);
        assert_eq!((memory[0].as_str(), memory[1].as_str()), ("memory", "ok"));
        assert_eq!(memory[10], "3");
        let single: u64 = memory[9].parse().unwrap();
        assert!(single >= 150, "peak single {single} MB");
        let total_mb: u64 = memory[8].parse().unwrap();
        assert!(
            total_mb >= single,
            "peak total {total_mb} MB < single {single} MB"
        );
        assert_eq!(
            (broken[0].as_str(), broken[1].as_str(), broken[10].as_str()),
            ("broken", "failed", "1")
        );
        assert_eq!(
            (skipped[1].as_str(), skipped[4].as_str()),
            ("skipped", "NA")
        );
        assert_eq!(
            (total[0].as_str(), total[1].as_str(), total[10].as_str()),
            ("total", "failed", "4")
        );
    }

    #[test]
    fn substeps_get_their_own_rows_and_count_once() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("benchmark.tsv");
        let bench = Benchmark::create(&path).unwrap();
        let result: Result<(), String> = bench.step("strains", |rec| {
            thread::scope(|scope| {
                for name in ["s1", "s2"] {
                    scope.spawn(move || {
                        rec.substep(name, |sub| {
                            let report = Pipeline::new(Cmd::new("true"))
                                .run()
                                .map_err(|e| e.to_string())?;
                            sub.record(&report);
                            Ok::<(), String>(())
                        })
                    });
                }
            });
            Ok(())
        });
        result.unwrap();
        bench.finish(Status::Ok);
        let rows = rows(&path);
        let names: Vec<&str> = rows.iter().map(|r| r[0].as_str()).collect();
        assert!(
            names.contains(&"strains/s1") && names.contains(&"strains/s2"),
            "{names:?}"
        );
        let strains = rows.iter().find(|r| r[0] == "strains").unwrap();
        assert_eq!(strains[10], "2");
        let sub = rows.iter().find(|r| r[0] == "strains/s1").unwrap();
        assert_eq!((sub[8].as_str(), sub[10].as_str()), ("NA", "1"));
        assert_eq!(rows.last().unwrap()[10], "2");
    }

    #[test]
    fn process_tree_includes_self() {
        assert!(process_tree_rss_kb(std::process::id()) > 0);
    }
}
