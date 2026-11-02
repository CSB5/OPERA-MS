//! The external tools: every tool is the first working copy in `PATH`
//! (normally a conda environment made from `environment.yml`).
//! `opera-ms check-dependency` tests and reports every tool with its
//! version and directory; a run resolves the tools it needs the same way,
//! before its first step.
//!
//! A tool "works" when its test passes: the command's exit status, or, for
//! the tools whose usage message exits with an error (`bwa`), a word in its
//! output.

use std::path::{Path, PathBuf};
use std::time::Duration;

use opera_ms_core::runner::{Cmd, Output, Pipeline};

/// An external tool and its test.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Tool {
    Samtools,
    Minibwa,
    Minimap2,
    Sylph,
    Skani,
    Curl,
    Megahit,
    Racon,
    Freebayes,
    Polypolish,
}

impl Tool {
    pub const ALL: [Tool; 10] = [
        Tool::Samtools,
        Tool::Minibwa,
        Tool::Minimap2,
        Tool::Sylph,
        Tool::Skani,
        Tool::Curl,
        Tool::Megahit,
        Tool::Racon,
        Tool::Freebayes,
        Tool::Polypolish,
    ];

    pub fn name(self) -> &'static str {
        match self {
            Tool::Samtools => "samtools",
            Tool::Minibwa => "minibwa",
            Tool::Minimap2 => "minimap2",
            Tool::Sylph => "sylph",
            Tool::Skani => "skani",
            Tool::Curl => "curl",
            Tool::Megahit => "megahit",
            Tool::Racon => "racon",
            Tool::Freebayes => "freebayes",
            Tool::Polypolish => "polypolish",
        }
    }

    /// The file looked for in `PATH`.
    fn program(self) -> &'static str {
        self.name()
    }

    /// Other files the tool's directory must have.
    fn companions(self) -> &'static [&'static str] {
        &[]
    }

    /// The test: arguments, and a word the output must contain (`None`:
    /// the exit status must be 0). The output also gives the version.
    fn test(self) -> (&'static [&'static str], Option<&'static str>) {
        match self {
            Tool::Minibwa => (&["version"], None),
            _ => (&["--version"], None),
        }
    }
}

/// A working tool.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Found {
    /// The tool's directory, with a trailing slash.
    pub dir: String,
    /// The version line the tool printed (empty if none).
    pub version: String,
}

const TEST_TIMEOUT: Duration = Duration::from_secs(60);

/// Runs a tool's test in `dir`: its version line, or why it failed.
fn works(tool: Tool, dir: &str) -> Result<String, String> {
    let program = format!("{dir}{}", tool.program());
    if !Path::new(&program).is_file() {
        return Err(format!("{program} not found"));
    }
    for companion in tool.companions() {
        let path = format!("{dir}{companion}");
        if !Path::new(&path).is_file() {
            return Err(format!("{path} not found"));
        }
    }
    let (args, word) = tool.test();
    let cmd = Cmd::new(&program).args(args);
    let output = capture(Pipeline::new(cmd).timeout(Some(TEST_TIMEOUT)));
    match (output, word) {
        (Ok((true, text)), None) => Ok(version_line(&text)),
        (Ok((_, text)), Some(w)) if text.contains(w) => Ok(version_line(&text)),
        (Ok((_, text)), _) => {
            let first = text.lines().next().unwrap_or("").trim();
            Err(format!("{program} is not functional ({first})"))
        }
        (Err(e), _) => Err(format!("{program} could not run: {e}")),
    }
}

/// The first line of a tool's output that holds a version number.
fn version_line(text: &str) -> String {
    text.lines()
        .map(str::trim)
        .find(|l| {
            l.as_bytes()
                .windows(3)
                .any(|w| w[0].is_ascii_digit() && w[1] == b'.' && w[2].is_ascii_digit())
        })
        .unwrap_or("")
        .to_string()
}

/// Runs a pipeline with stdout and stderr in a temporary file: whether it
/// succeeded, and its output.
fn capture(pipeline: Pipeline) -> Result<(bool, String), String> {
    use std::sync::atomic::{AtomicUsize, Ordering};
    static COUNTER: AtomicUsize = AtomicUsize::new(0);
    let path = std::env::temp_dir().join(format!(
        "opera-ms-check-{}-{}",
        std::process::id(),
        COUNTER.fetch_add(1, Ordering::SeqCst)
    ));
    std::fs::write(&path, "").map_err(|e| format!("cannot write {}: {e}", path.display()))?;
    let result = pipeline
        .stdout(Output::Append(path.clone()))
        .stderr(Output::Append(path.clone()))
        .run();
    let text = std::fs::read(&path)
        .map(|b| String::from_utf8_lossy(&b).into_owned())
        .unwrap_or_default();
    let _ = std::fs::remove_file(&path);
    match result {
        Ok(_) => Ok((true, text)),
        Err(e) if e.report().is_some() && !e.is_timeout() => Ok((false, text)),
        Err(e) => Err(e.to_string()),
    }
}

/// The directories of `PATH`.
fn path_dirs() -> Vec<PathBuf> {
    std::env::var_os("PATH")
        .map(|p| std::env::split_paths(&p).collect())
        .unwrap_or_default()
}

/// Candidate directories in `PATH`: where the program is.
fn path_candidates(tool: Tool) -> Vec<String> {
    path_dirs()
        .into_iter()
        .filter(|dir| dir.join(tool.program()).is_file())
        .map(|dir| format!("{}/", dir.display()))
        .collect()
}

/// Resolves one tool: the first working copy in `PATH`, or why none was
/// found (the problem of the first candidate, if any).
pub fn resolve(tool: Tool) -> Result<Found, String> {
    let mut first_problem = None;
    for dir in path_candidates(tool) {
        match works(tool, &dir) {
            Ok(version) => return Ok(Found { dir, version }),
            Err(e) => {
                first_problem.get_or_insert(e);
            }
        }
    }
    Err(first_problem.unwrap_or_else(|| format!("{} not found in PATH", tool.program())))
}

/// `opera-ms check-dependency`: one line per tool; `Err` lists the tools
/// found nowhere.
pub fn check_dependency() -> Result<String, String> {
    let mut report = String::new();
    let mut missing = Vec::new();
    for tool in Tool::ALL {
        match resolve(tool) {
            Ok(found) => report.push_str(&format!(
                " *** {:<9} {:<45} {}\n",
                tool.name(),
                found.version,
                found.dir
            )),
            Err(problem) => {
                report.push_str(&format!(" *** {:<9} NOT FOUND: {problem}\n", tool.name()));
                missing.push(tool.name());
            }
        }
    }
    if missing.is_empty() {
        report.push_str("\n *** All software are functional.\n");
        Ok(report)
    } else {
        Err(format!("{report}\nNot functional: {}", missing.join(", ")))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn missing_program_is_reported() {
        let dir = tempfile::tempdir().unwrap();
        let d = format!("{}/", dir.path().display());
        let err = works(Tool::Skani, &d).unwrap_err();
        assert!(err.contains("not found"));
    }

    #[test]
    fn exit_status_and_version() {
        let dir = tempfile::tempdir().unwrap();
        let d = format!("{}/", dir.path().display());
        // a fake skani that works, a fake minibwa that fails
        for (name, body) in [
            ("skani", "echo 'skani 0.4.0'; exit 0"),
            ("minibwa", "echo 'ERROR'; exit 1"),
        ] {
            let path = dir.path().join(name);
            std::fs::write(&path, format!("#!/bin/sh\n{body}\n")).unwrap();
            use std::os::unix::fs::PermissionsExt;
            std::fs::set_permissions(&path, std::fs::Permissions::from_mode(0o755)).unwrap();
        }
        assert_eq!(works(Tool::Skani, &d).unwrap(), "skani 0.4.0");
        assert!(
            works(Tool::Minibwa, &d)
                .unwrap_err()
                .contains("not functional")
        );
    }

    #[test]
    fn version_lines() {
        assert_eq!(
            version_line("\nProgram: bwa\nVersion: 0.7.19-r1273\n"),
            "Version: 0.7.19-r1273"
        );
        assert_eq!(version_line("2.31-r1302\n"), "2.31-r1302");
        assert_eq!(version_line("no number\n"), "");
    }
}
