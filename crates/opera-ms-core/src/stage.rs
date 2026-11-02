//! Pipeline stages: progress lines and resume.
//!
//! Progress lines are printed on stdout in the format of OPERA-MS 0.9:
//!
//! ```text
//!
//! [Sun Oct 25 06:48:20 2026]<TAB>Short-read mapping and coverage estimation [3/9]
//! ```
//!
//! followed, for a stage that doesn't run, by `[...]\tSkip ` or
//! `[...]\tStage previously completed`.
//!
//! Resume: once one stage runs, every later stage
//! runs too, even if a previous run completed it. A stage counts as
//! completed only if its marker exists, and a marker is written (as a
//! temporary file renamed into place) only after the stage's outputs have
//! been checked, so a stage interrupted half-way is never
//! skipped. Markers live in `<intermediate dir>/.opera-ms/`, so the output
//! layout doesn't change.

use std::fs;
use std::io;
use std::path::{Path, PathBuf};

use chrono::Local;

/// What to do with a stage.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Decision {
    Run,
    /// Turned off by the options (`Skip`).
    Skip,
    /// Completed by a previous run (`Stage previously completed`).
    AlreadyDone,
}

/// Progress and resume state of a run.
#[derive(Debug)]
pub struct Stages {
    marker_dir: PathBuf,
    total: u32,
    /// Set once a stage runs: every later stage runs too.
    restart: bool,
}

/// The time format of Perl's scalar `localtime` ("Sun Oct 25 06:48:20 2026").
pub fn perl_localtime() -> String {
    Local::now().format("%a %b %e %H:%M:%S %Y").to_string()
}

impl Stages {
    /// `restart`: run every stage even if completed before (a stage was
    /// given in config file mode).
    pub fn new(intermediate_dir: &Path, total: u32, restart: bool) -> Self {
        Self {
            marker_dir: intermediate_dir.join(".opera-ms"),
            total,
            restart,
        }
    }

    fn marker(&self, id: &str) -> PathBuf {
        self.marker_dir.join(format!("{id}.done"))
    }

    /// Prints the progress line of stage `number` and decides whether it
    /// runs. `id` names its marker; `skip`: the stage is turned off.
    pub fn begin(&mut self, id: &str, message: &str, number: u32, skip: bool) -> Decision {
        let time = perl_localtime();
        crate::say!("\n[{time}]\t{message} [{number}/{}]", self.total);
        if skip {
            crate::say!("[{time}]\tSkip ");
            Decision::Skip
        } else if !self.restart && self.marker(id).exists() {
            crate::say!("[{time}]\tStage previously completed");
            Decision::AlreadyDone
        } else {
            self.restart = true;
            Decision::Run
        }
    }

    /// Records that stage `id` completed (call after checking its outputs).
    pub fn complete(&self, id: &str) -> io::Result<()> {
        fs::create_dir_all(&self.marker_dir)?;
        let marker = self.marker(id);
        let tmp = marker.with_extension("done.tmp");
        fs::write(&tmp, format!("{}\n", perl_localtime()))?;
        fs::rename(&tmp, &marker)
    }

    /// Forgets that stage `id` completed (before running it again).
    pub fn invalidate(&self, id: &str) -> io::Result<()> {
        match fs::remove_file(self.marker(id)) {
            Err(err) if err.kind() != io::ErrorKind::NotFound => Err(err),
            _ => Ok(()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn restart_rule_and_markers() {
        let dir = tempfile::tempdir().unwrap();
        let mut stages = Stages::new(dir.path(), 9, false);
        assert_eq!(stages.begin("a", "Stage A", 1, false), Decision::Run);
        stages.complete("a").unwrap();
        assert_eq!(stages.begin("b", "Stage B", 2, true), Decision::Skip);
        stages.complete("c").unwrap();

        // second run: a was completed, so it is skipped; c ran before too,
        // but runs again because b (not completed) runs before it
        let mut stages = Stages::new(dir.path(), 9, false);
        assert_eq!(
            stages.begin("a", "Stage A", 1, false),
            Decision::AlreadyDone
        );
        assert_eq!(stages.begin("b", "Stage B", 2, false), Decision::Run);
        assert_eq!(stages.begin("c", "Stage C", 3, false), Decision::Run);

        // restart: everything runs
        let mut stages = Stages::new(dir.path(), 9, true);
        assert_eq!(stages.begin("a", "Stage A", 1, false), Decision::Run);
        stages.invalidate("a").unwrap();
        stages.invalidate("missing").unwrap();
        let mut stages = Stages::new(dir.path(), 9, false);
        assert_eq!(stages.begin("a", "Stage A", 1, false), Decision::Run);
    }

    #[test]
    fn perl_time_format() {
        let time = perl_localtime();
        // "Sun Oct 25 06:48:20 2026" / "Sun Oct  5 06:48:20 2026"
        assert_eq!(time.len(), 24, "{time}");
        assert_eq!(&time[3..4], " ");
    }
}
