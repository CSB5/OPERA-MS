//! The main log: everything `opera-ms` prints, progress lines
//! (stdout) and commands, notes and errors (stderr), also goes to one file,
//! `<out-dir>/opera-ms.log`, each line prefixed with the time in the file
//! only (the terminal output is unchanged). An error message carries the
//! failing command and the end of its error output (see `runner`), so the
//! main log alone is enough to debug a failed run; each tool's full output
//! is kept in its own file.
//!
//! Use the [`say!`](crate::say) (stdout) and [`note!`](crate::note)
//! (stderr) macros instead of `println!` / `eprintln!`. Every command the
//! runner starts is written to the file only, with its wall time, CPU time
//! and peak memory when it ends (see `runner::Pipeline::run`).

use std::fs::{File, OpenOptions};
use std::io::{self, Write};
use std::path::Path;
use std::sync::Mutex;

static LOG: Mutex<Option<File>> = Mutex::new(None);

/// Starts writing the log to `path` (appended: a resumed run continues the
/// same log).
pub fn open(path: &Path) -> io::Result<()> {
    let file = OpenOptions::new().create(true).append(true).open(path)?;
    *LOG.lock().unwrap_or_else(|e| e.into_inner()) = Some(file);
    Ok(())
}

/// Stops writing the log file.
pub fn close() {
    *LOG.lock().unwrap_or_else(|e| e.into_inner()) = None;
}

fn to_file(text: &str) {
    let mut guard = LOG.lock().unwrap_or_else(|e| e.into_inner());
    if let Some(file) = guard.as_mut() {
        let time = chrono::Local::now().format("%Y-%m-%d %H:%M:%S");
        let mut out = String::with_capacity(text.len() + 32);
        for line in text.split('\n') {
            out.push_str(&format!("{time} | {line}\n"));
        }
        // a log write error mustn't stop the run
        let _ignored = file.write_all(out.as_bytes());
    }
}

/// A line in the log file only (the commands the runner starts and how
/// they ended; too long for the terminal).
pub fn file_only(text: &str) {
    to_file(text);
}

/// A line on stdout and in the log.
pub fn out(text: &str) {
    to_file(text);
    let mut stdout = io::stdout().lock();
    // stdout write errors (closed pipe) mustn't stop the run
    let _ignored = writeln!(stdout, "{text}");
    let _ignored = stdout.flush();
}

/// A line on stderr and in the log.
pub fn err(text: &str) {
    to_file(text);
    let _ignored = writeln!(io::stderr().lock(), "{text}");
}

/// `println!` that also writes the main log.
#[macro_export]
macro_rules! say {
    ($($arg:tt)*) => { $crate::log::out(&format!($($arg)*)) };
}

/// `eprintln!` that also writes the main log.
#[macro_export]
macro_rules! note {
    ($($arg:tt)*) => { $crate::log::err(&format!($($arg)*)) };
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lines_go_to_the_file_with_a_time() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("opera-ms.log");
        open(&path).unwrap();
        crate::note!("a\nb");
        crate::say!("c");
        close();
        crate::note!("not logged");
        let text = std::fs::read_to_string(&path).unwrap();
        // other tests running at the same time may log commands into the
        // shared log: look at this test's lines only
        let lines: Vec<&str> = text
            .lines()
            .filter(|l| l.ends_with(" | a") || l.ends_with(" | b") || l.ends_with(" | c"))
            .collect();
        assert_eq!(lines.len(), 3, "{text}");
        assert!(
            lines[0].ends_with(" | a") && lines[1].ends_with(" | b") && lines[2].ends_with(" | c")
        );
        assert_eq!(lines[0].len(), "2026-10-28 13:05:01 | a".len());
        assert!(!text.contains("not logged"));
    }
}
