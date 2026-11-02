//! Running external commands and pipelines, with no silent failures.
//!
//! - No shell: every command is started directly, and a pipeline
//!   (`a | b | c`) is built here by connecting the processes' stdin and
//!   stdout, so that each process is a child we wait for and check.
//! - Every process of a pipeline is checked. If one fails, the others are
//!   killed and the run fails, naming the process that failed first.
//! - Death by a signal is reported with the signal name; `SIGKILL` is
//!   reported as probably the out-of-memory killer, with the peak memory.
//! - Timeouts kill the whole process group of the pipeline.
//! - Each run returns the exit status and resource usage (CPU time, peak
//!   memory) of every process, for the per-step benchmark.
//!
//! Children run in their own process group, so a Ctrl-C on the terminal
//! doesn't reach them: the binary must call [`kill_all_running`] from its
//! signal handler.

use std::collections::HashSet;
use std::ffi::{OsStr, OsString};
use std::fs::{File, OpenOptions};
use std::io::{self, Read, Seek, SeekFrom};
use std::os::fd::AsFd;
use std::os::unix::process::CommandExt;
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::sync::mpsc;
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::{Duration, Instant};

use opera_ms_sys::process::{self as sys, Exit, Usage, signal};

/// Lines of a failed command's stderr quoted in the error message.
const STDERR_TAIL_LINES: usize = 20;

/// Process groups of the pipelines currently running.
static RUNNING_GROUPS: Mutex<Option<HashSet<u32>>> = Mutex::new(None);

fn register_group(pgid: u32) {
    if let Ok(mut groups) = RUNNING_GROUPS.lock() {
        groups.get_or_insert_with(HashSet::new).insert(pgid);
    }
}

fn unregister_group(pgid: u32) {
    if let Ok(mut groups) = RUNNING_GROUPS.lock() {
        if let Some(groups) = groups.as_mut() {
            groups.remove(&pgid);
        }
    }
}

/// Kills every pipeline that is still running (for the Ctrl-C / SIGTERM
/// handler of the binary). The runs then fail as killed.
pub fn kill_all_running() {
    if let Ok(groups) = RUNNING_GROUPS.lock() {
        for pgid in groups.iter().flatten() {
            // Best effort while the program is being interrupted.
            let _ignored = sys::kill_process_group(*pgid, signal::SIGKILL);
        }
    }
}

/// One command: a program and its arguments.
#[derive(Debug, Clone)]
pub struct Cmd {
    program: PathBuf,
    args: Vec<OsString>,
    env: Vec<(OsString, OsString)>,
    current_dir: Option<PathBuf>,
    /// Error output of this command (default: the pipeline's).
    stderr: Option<Output>,
}

impl Cmd {
    pub fn new(program: impl Into<PathBuf>) -> Self {
        Self {
            program: program.into(),
            args: Vec::new(),
            env: Vec::new(),
            current_dir: None,
            stderr: None,
        }
    }

    pub fn arg(mut self, arg: impl AsRef<OsStr>) -> Self {
        self.args.push(arg.as_ref().to_os_string());
        self
    }

    pub fn args<I, S>(mut self, args: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: AsRef<OsStr>,
    {
        self.args
            .extend(args.into_iter().map(|a| a.as_ref().to_os_string()));
        self
    }

    pub fn env(mut self, key: impl AsRef<OsStr>, value: impl AsRef<OsStr>) -> Self {
        self.env
            .push((key.as_ref().to_os_string(), value.as_ref().to_os_string()));
        self
    }

    pub fn current_dir(mut self, dir: impl Into<PathBuf>) -> Self {
        self.current_dir = Some(dir.into());
        self
    }

    /// Error output of this command only (the pipeline's `stderr` is used
    /// otherwise), like a `2>` on one command of a shell pipeline.
    pub fn stderr(mut self, output: Output) -> Self {
        self.stderr = Some(output);
        self
    }

    /// The program name, for messages.
    pub fn program_name(&self) -> String {
        self.program
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_else(|| self.program.display().to_string())
    }

    /// The command line, shell-quoted, for logs and messages.
    pub fn display(&self) -> String {
        std::iter::once(self.program.as_os_str())
            .chain(self.args.iter().map(OsString::as_os_str))
            .map(|a| shell_quote(&a.to_string_lossy()))
            .collect::<Vec<_>>()
            .join(" ")
    }

    fn command(&self) -> Command {
        let mut command = Command::new(&self.program);
        command.args(&self.args);
        for (key, value) in &self.env {
            command.env(key, value);
        }
        if let Some(dir) = &self.current_dir {
            command.current_dir(dir);
        }
        command
    }
}

fn shell_quote(arg: &str) -> String {
    let safe = !arg.is_empty()
        && arg
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || "-_./=:,+@%".contains(c));
    if safe {
        arg.to_string()
    } else {
        format!("'{}'", arg.replace('\'', r"'\''"))
    }
}

/// Where a pipeline reads its input from.
#[derive(Debug, Clone, Default)]
pub enum Input {
    #[default]
    Null,
    File(PathBuf),
    /// The read end of an OS pipe (see [`Pipeline::run_fed`]).
    Pipe(Arc<Mutex<Option<std::io::PipeReader>>>),
}

/// Where a pipeline's output (or error output) goes.
#[derive(Debug, Clone, Default)]
pub enum Output {
    #[default]
    Null,
    /// Truncate and write.
    File(PathBuf),
    /// Append (several commands sharing one log, like `2>>`).
    Append(PathBuf),
    /// This program's stderr (what `run_exe` did with a command's stdout:
    /// it printed it into the log).
    Stderr,
}

impl Output {
    fn open(&self) -> io::Result<Stdio> {
        Ok(match self {
            Output::Null => Stdio::null(),
            Output::File(path) => File::create(path)?.into(),
            Output::Append(path) => OpenOptions::new()
                .create(true)
                .append(true)
                .open(path)?
                .into(),
            Output::Stderr => Stdio::from(std::io::stderr().as_fd().try_clone_to_owned()?),
        })
    }

    fn path(&self) -> Option<&Path> {
        match self {
            Output::Null | Output::Stderr => None,
            Output::File(path) | Output::Append(path) => Some(path),
        }
    }
}

/// One or more commands connected stdout to stdin.
#[derive(Debug, Clone)]
pub struct Pipeline {
    cmds: Vec<Cmd>,
    stdin: Input,
    stdout: Output,
    stderr: Output,
    timeout: Option<Duration>,
}

impl Pipeline {
    /// A pipeline of a single command.
    pub fn new(cmd: Cmd) -> Self {
        Self {
            cmds: vec![cmd],
            stdin: Input::Null,
            stdout: Output::Null,
            stderr: Output::Null,
            timeout: None,
        }
    }

    /// Adds a command reading the previous command's output (`| cmd`).
    pub fn pipe(mut self, cmd: Cmd) -> Self {
        self.cmds.push(cmd);
        self
    }

    pub fn stdin(mut self, input: Input) -> Self {
        self.stdin = input;
        self
    }

    pub fn stdout(mut self, output: Output) -> Self {
        self.stdout = output;
        self
    }

    /// Error output of every command of the pipeline.
    pub fn stderr(mut self, output: Output) -> Self {
        self.stderr = output;
        self
    }

    /// Kills the pipeline after `timeout` (`None`: no timeout).
    pub fn timeout(mut self, timeout: Option<Duration>) -> Self {
        self.timeout = timeout;
        self
    }

    /// The whole pipeline, shell-quoted, for logs and messages.
    pub fn display(&self) -> String {
        let mut text = self
            .cmds
            .iter()
            .map(Cmd::display)
            .collect::<Vec<_>>()
            .join(" | ");
        if let Input::File(path) = &self.stdin {
            text.push_str(&format!(" < {}", shell_quote(&path.to_string_lossy())));
        }

        match &self.stdout {
            Output::File(path) => {
                text.push_str(&format!(" > {}", shell_quote(&path.to_string_lossy())))
            }
            Output::Append(path) => {
                text.push_str(&format!(" >> {}", shell_quote(&path.to_string_lossy())))
            }
            Output::Stderr => text.push_str(" >&2"),
            Output::Null => {}
        }
        match &self.stderr {
            Output::File(path) => {
                text.push_str(&format!(" 2> {}", shell_quote(&path.to_string_lossy())))
            }
            Output::Append(path) => {
                text.push_str(&format!(" 2>> {}", shell_quote(&path.to_string_lossy())))
            }
            Output::Stderr | Output::Null => {}
        }
        text
    }

    /// Runs the pipeline and waits for every process.
    pub fn run(&self) -> Result<RunReport, RunError> {
        // numbered, so that parallel commands and their ends can be paired
        static NEXT: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(1);
        let id = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let command = self.display();
        tracing::info!("running: {command}");
        crate::log::file_only(&format!("$ #{id} {command}"));
        let result = run_pipeline(self);
        match &result {
            Ok(report) => crate::log::file_only(&format!(
                "  #{id} done: {:.1} s wall, {:.1} s CPU, {} MB peak",
                report.wall.as_secs_f64(),
                report.cpu().as_secs_f64(),
                report.max_rss_kb() / 1024
            )),
            Err(err) => crate::log::file_only(&format!("  #{id} failed: {err}")),
        }
        result
    }

    /// Runs the pipeline with the first command's stdin written by `feed`
    /// (in a thread, through an OS pipe). A failure of the pipeline is
    /// reported first (the feeder then sees a broken pipe); a failure of
    /// `feed` alone (e.g. a truncated input file) fails the run too, even
    /// if the commands exited normally on the shortened input.
    pub fn run_fed<F>(mut self, feed: F) -> Result<RunReport, RunError>
    where
        F: FnOnce(&mut std::io::PipeWriter) -> io::Result<()> + Send + 'static,
    {
        let command = self.display();
        let (reader, mut writer) = std::io::pipe().map_err(|source| RunError::Input {
            command: command.clone(),
            source,
            report: RunReport::default(),
        })?;
        self.stdin = Input::Pipe(Arc::new(Mutex::new(Some(reader))));
        let feeder = thread::spawn(move || feed(&mut writer));
        let result = self.run();
        let fed = feeder
            .join()
            .unwrap_or_else(|_| Err(io::Error::other("the input feeder panicked")));
        match (result, fed) {
            (Err(err), _) => Err(err),
            (Ok(report), Err(source)) => Err(RunError::Input {
                command,
                source,
                report,
            }),
            (Ok(report), Ok(())) => Ok(report),
        }
    }
}

/// How one process of a run ended.
#[derive(Debug, Clone)]
pub struct ProcessReport {
    pub program: String,
    pub exit: Exit,
    pub usage: Usage,
    /// Killed by the runner (timeout, or another process of the pipeline
    /// failed), as opposed to failing on its own.
    pub killed_by_runner: bool,
}

/// Exit status and resource usage of every process of a run.
#[derive(Debug, Clone, Default)]
pub struct RunReport {
    pub processes: Vec<ProcessReport>,
    pub wall: Duration,
}

impl RunReport {
    /// Total CPU time (user + system) of all processes.
    pub fn cpu(&self) -> Duration {
        self.processes
            .iter()
            .map(|p| p.usage.user + p.usage.system)
            .sum()
    }

    /// Peak memory of the largest single process, in kilobytes.
    pub fn max_rss_kb(&self) -> u64 {
        self.processes
            .iter()
            .map(|p| p.usage.max_rss_kb)
            .max()
            .unwrap_or(0)
    }
}

/// A command or pipeline that didn't complete successfully.
#[derive(Debug, thiserror::Error)]
pub enum RunError {
    #[error("cannot start {program}: {source}\n  command: {command}")]
    Spawn {
        program: String,
        command: String,
        source: io::Error,
    },
    #[error("cannot open {path} for {command}: {source}")]
    Redirect {
        path: PathBuf,
        command: String,
        source: io::Error,
    },
    #[error("waiting for {program}: {source}")]
    Wait { program: String, source: io::Error },
    #[error("{0}")]
    Failed(Box<Failure>),
    #[error("cannot feed the input of {command}: {source}")]
    Input {
        command: String,
        source: io::Error,
        report: RunReport,
    },
    #[error("timed out after {}, killed\n  command: {command}", humanize(*.after))]
    TimedOut {
        command: String,
        after: Duration,
        report: RunReport,
    },
}

/// Details of a process that failed (see [`RunError::Failed`]).
#[derive(Debug)]
pub struct Failure {
    /// The process that failed first.
    pub program: String,
    /// How it failed ("exited with code 2", "was killed by SIGKILL, ...").
    pub failure: String,
    /// The whole command line.
    pub command: String,
    /// The last lines of the error output.
    pub stderr_tail: String,
    pub stderr_path: Option<PathBuf>,
    pub report: RunReport,
}

impl std::fmt::Display for Failure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{} {}\n  command: {}{}",
            self.program,
            self.failure,
            self.command,
            stderr_context(&self.stderr_tail, &self.stderr_path)
        )
    }
}

impl RunError {
    /// Resource usage of the processes that ran, if any.
    pub fn report(&self) -> Option<&RunReport> {
        match self {
            RunError::Failed(failure) => Some(&failure.report),
            RunError::TimedOut { report, .. } | RunError::Input { report, .. } => Some(report),
            _ => None,
        }
    }

    pub fn is_timeout(&self) -> bool {
        matches!(self, RunError::TimedOut { .. })
    }
}

fn stderr_context(tail: &str, path: &Option<PathBuf>) -> String {
    let mut text = String::new();
    if !tail.is_empty() {
        text.push_str("\n  last lines of its error output:");
        for line in tail.lines() {
            text.push_str("\n    ");
            text.push_str(line);
        }
    }
    if let Some(path) = path {
        text.push_str(&format!("\n  full error output: {}", path.display()));
    }
    text
}

/// "5m", "1h 2m", "42s".
pub fn humanize(d: Duration) -> String {
    let secs = d.as_secs();
    let (h, m, s) = (secs / 3600, (secs % 3600) / 60, secs % 60);
    match (h, m, s) {
        (0, 0, s) => format!("{s}s"),
        (0, m, 0) => format!("{m}m"),
        (0, m, s) => format!("{m}m {s}s"),
        (h, 0, 0) => format!("{h}h"),
        (h, m, _) => format!("{h}h {m}m"),
    }
}

fn describe_failure(exit: Exit, usage: &Usage) -> String {
    match exit {
        Exit::Code(code) => format!("exited with code {code}"),
        Exit::Signal(signal::SIGKILL) => format!(
            "was killed by SIGKILL, probably by the out-of-memory killer (its peak memory was \
             {} MB; check `dmesg`, or lower --num-processors)",
            usage.max_rss_kb / 1024
        ),
        Exit::Signal(sig @ (signal::SIGSEGV | signal::SIGABRT)) => {
            format!("crashed ({})", signal::name(sig))
        }
        Exit::Signal(signal::SIGPIPE) => {
            "was killed by SIGPIPE: the next command of the pipeline stopped reading its output"
                .to_string()
        }
        Exit::Signal(sig) => format!("was killed by {}", signal::name(sig)),
    }
}

/// The last `n` lines of a file (best effort: only the end is read).
fn tail(path: &Path, n: usize) -> String {
    let Ok(mut file) = File::open(path) else {
        return String::new();
    };
    let len = file.metadata().map(|m| m.len()).unwrap_or(0);
    let start = len.saturating_sub(16 * 1024);
    if file.seek(SeekFrom::Start(start)).is_err() {
        return String::new();
    }
    let mut buf = Vec::new();
    if file.read_to_end(&mut buf).is_err() {
        return String::new();
    }
    let text = String::from_utf8_lossy(&buf);
    let lines: Vec<&str> = text.lines().collect();
    lines[lines.len().saturating_sub(n)..].join("\n")
}

struct Started {
    children: Vec<Child>,
    pgid: u32,
}

fn spawn_all(pipeline: &Pipeline, command_text: &str) -> Result<Started, RunError> {
    let redirect = |path: Option<&Path>, source| RunError::Redirect {
        path: path.map(Path::to_path_buf).unwrap_or_default(),
        command: command_text.to_string(),
        source,
    };
    // Every process gets its own append handle on the stderr file (after
    // one truncation for `Output::File`), so their lines can't overwrite
    // each other.
    if let Output::File(path) = &pipeline.stderr {
        File::create(path).map_err(|e| redirect(Some(path), e))?;
    }
    let open_stderr = || -> Result<Stdio, RunError> {
        if matches!(pipeline.stderr, Output::Stderr) {
            return Ok(Stdio::inherit());
        }
        match pipeline.stderr.path() {
            None => Ok(Stdio::null()),
            Some(path) => OpenOptions::new()
                .create(true)
                .append(true)
                .open(path)
                .map(Stdio::from)
                .map_err(|e| redirect(Some(path), e)),
        }
    };

    let mut started = Started {
        children: Vec::new(),
        pgid: 0,
    };
    let mut previous_stdout = None;
    let last = pipeline.cmds.len() - 1;
    for (index, cmd) in pipeline.cmds.iter().enumerate() {
        let mut command = cmd.command();
        let stdin = match (index, &pipeline.stdin, previous_stdout.take()) {
            (_, _, Some(stdout)) => Stdio::from(stdout),
            (0, Input::File(path), None) => File::open(path)
                .map(Stdio::from)
                .map_err(|e| redirect(Some(path), e))?,
            (0, Input::Pipe(reader), None) => match reader.lock().ok().and_then(|mut r| r.take()) {
                Some(reader) => Stdio::from(reader),
                None => Stdio::null(),
            },
            _ => Stdio::null(),
        };
        command.stdin(stdin);
        command.stdout(if index == last {
            pipeline
                .stdout
                .open()
                .map_err(|e| redirect(pipeline.stdout.path(), e))?
        } else {
            Stdio::piped()
        });
        command.stderr(match &cmd.stderr {
            Some(Output::File(path)) => File::create(path)
                .map(Stdio::from)
                .map_err(|e| redirect(Some(path), e))?,
            Some(output) => output.open().map_err(|e| redirect(output.path(), e))?,
            None => open_stderr()?,
        });
        command.process_group(if index == 0 { 0 } else { started.pgid as i32 });

        match command.spawn() {
            Ok(mut child) => {
                if index == 0 {
                    started.pgid = child.id();
                    register_group(started.pgid);
                }
                previous_stdout = child.stdout.take();
                started.children.push(child);
            }
            Err(source) => {
                if !started.children.is_empty() {
                    // Best effort: the spawn error is what gets reported.
                    let _ignored = sys::kill_process_group(started.pgid, signal::SIGKILL);
                    for child in &started.children {
                        let _ignored = sys::reap(child.id());
                    }
                    unregister_group(started.pgid);
                }
                return Err(RunError::Spawn {
                    program: cmd.program_name(),
                    command: command_text.to_string(),
                    source,
                });
            }
        }
    }
    Ok(started)
}

fn run_pipeline(pipeline: &Pipeline) -> Result<RunReport, RunError> {
    let command_text = pipeline.display();
    let start = Instant::now();
    let started = spawn_all(pipeline, &command_text)?;
    let pids: Vec<u32> = started.children.iter().map(Child::id).collect();

    // One thread per process waits for it to exit without reaping it, so
    // its PID and the process group stay reserved until we reap below.
    let (sender, receiver) = mpsc::channel();
    for (index, pid) in pids.iter().copied().enumerate() {
        let sender = sender.clone();
        thread::spawn(move || {
            // The receiver only goes away once every process has exited.
            let _ignored = sender.send((index, sys::wait_exit_nowait(pid)));
        });
    }
    drop(sender);

    let deadline = pipeline.timeout.map(|t| start + t);
    let mut exited = vec![false; pids.len()];
    let mut killed_by_runner = vec![false; pids.len()];
    let mut first_failure: Option<usize> = None;
    let mut timed_out = false;
    let mut wait_error = None;
    let mut remaining = pids.len();
    let kill_rest = |exited: &[bool], killed: &mut [bool]| {
        // ESRCH (all gone) is already Ok; another error leaves the
        // processes to finish on their own, and we still wait for them.
        let _ignored = sys::kill_process_group(started.pgid, signal::SIGKILL);
        for (index, done) in exited.iter().enumerate() {
            if !done {
                killed[index] = true;
            }
        }
    };
    while remaining > 0 {
        let event = match deadline {
            Some(deadline) if !timed_out => {
                match receiver.recv_timeout(deadline.saturating_duration_since(Instant::now())) {
                    Ok(event) => event,
                    Err(mpsc::RecvTimeoutError::Timeout) => {
                        timed_out = true;
                        kill_rest(&exited, &mut killed_by_runner);
                        continue;
                    }
                    Err(mpsc::RecvTimeoutError::Disconnected) => break,
                }
            }
            _ => match receiver.recv() {
                Ok(event) => event,
                Err(_) => break,
            },
        };
        let (index, result) = event;
        exited[index] = true;
        remaining -= 1;
        match result {
            Ok(exit) if !exit.success() && !killed_by_runner[index] && first_failure.is_none() => {
                first_failure = Some(index);
                kill_rest(&exited, &mut killed_by_runner);
            }
            Ok(_) => {}
            Err(source) => {
                wait_error.get_or_insert((index, source));
                kill_rest(&exited, &mut killed_by_runner);
            }
        }
    }

    let mut report = RunReport {
        processes: Vec::new(),
        wall: Duration::ZERO,
    };
    for (index, (pid, cmd)) in pids.iter().zip(&pipeline.cmds).enumerate() {
        match sys::reap(*pid) {
            Ok((exit, usage)) => report.processes.push(ProcessReport {
                program: cmd.program_name(),
                exit,
                usage,
                // Only if it really died of our SIGKILL: a process can exit
                // on its own just before the kill reaches the group.
                killed_by_runner: killed_by_runner[index] && exit == Exit::Signal(signal::SIGKILL),
            }),
            Err(source) => {
                wait_error.get_or_insert((index, source));
            }
        }
    }
    unregister_group(started.pgid);
    report.wall = start.elapsed();

    if let Some((index, source)) = wait_error {
        return Err(RunError::Wait {
            program: pipeline.cmds[index].program_name(),
            source,
        });
    }
    // The process to blame: the first one (in pipeline order) that failed
    // on its own (including a kill from outside, e.g. by the OOM killer).
    // A SIGPIPE death is the consequence of the next process stopping, so
    // it's blamed only if nothing else failed.
    let failed_on_its_own = |p: &ProcessReport| !p.exit.success() && !p.killed_by_runner;
    let failed = report
        .processes
        .iter()
        .position(|p| failed_on_its_own(p) && p.exit != Exit::Signal(signal::SIGPIPE))
        .or_else(|| report.processes.iter().position(failed_on_its_own));
    if timed_out && failed.is_none() {
        return Err(RunError::TimedOut {
            command: command_text,
            after: pipeline.timeout.unwrap_or_default(),
            report,
        });
    }
    if let Some(index) = failed {
        let process = &report.processes[index];
        let stderr_path = pipeline.stderr.path().map(Path::to_path_buf);
        return Err(RunError::Failed(Box::new(Failure {
            program: process.program.clone(),
            failure: describe_failure(process.exit, &process.usage),
            command: command_text,
            stderr_tail: stderr_path
                .as_deref()
                .map(|p| tail(p, STDERR_TAIL_LINES))
                .unwrap_or_default(),
            stderr_path,
            report,
        })));
    }
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sh(script: &str) -> Cmd {
        Cmd::new("sh").args(["-c", script])
    }

    #[test]
    fn successful_pipeline_reports_every_process() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out");
        let report = Pipeline::new(sh("printf 'b\\na\\n'"))
            .pipe(Cmd::new("sort"))
            .stdout(Output::File(out.clone()))
            .run()
            .unwrap();
        assert_eq!(std::fs::read_to_string(out).unwrap(), "a\nb\n");
        assert_eq!(report.processes.len(), 2);
        assert!(report.processes.iter().all(|p| p.exit.success()));
    }

    #[test]
    fn failure_in_the_middle_of_a_pipeline_is_reported() {
        let dir = tempfile::tempdir().unwrap();
        let err_log = dir.path().join("err");
        // "yes" would run forever: it must be killed when the middle fails
        let err = Pipeline::new(Cmd::new("yes"))
            .pipe(sh("head -n 5 >/dev/null; echo 'middle broke' >&2; exit 4"))
            .pipe(Cmd::new("cat"))
            .stderr(Output::File(err_log.clone()))
            .run()
            .unwrap_err();
        match &err {
            RunError::Failed(f) => {
                assert_eq!(f.program, "sh");
                assert_eq!(f.failure, "exited with code 4");
                assert!(f.stderr_tail.contains("middle broke"));
            }
            other => panic!("unexpected {other:?}"),
        }
        assert!(err.to_string().contains("full error output"));
    }

    #[test]
    fn failure_of_the_first_command_is_reported_even_if_the_last_succeeds() {
        let err = Pipeline::new(sh("echo partial; exit 2"))
            .pipe(Cmd::new("cat"))
            .run()
            .unwrap_err();
        assert!(matches!(&err, RunError::Failed(f) if f.failure == "exited with code 2"));
    }

    #[test]
    fn early_exit_of_a_reader_is_reported_as_sigpipe_of_the_writer() {
        // head stops reading after one line: yes gets SIGPIPE
        let err = Pipeline::new(Cmd::new("yes"))
            .pipe(Cmd::new("head").arg("-n1"))
            .run()
            .unwrap_err();
        match err {
            RunError::Failed(f) => {
                assert_eq!(f.program, "yes");
                assert!(f.failure.contains("SIGPIPE"), "{}", f.failure);
            }
            other => panic!("unexpected {other:?}"),
        }
    }

    #[test]
    fn sigkill_from_outside_is_reported_as_probable_oom() {
        let err = Pipeline::new(sh("kill -9 $$")).run().unwrap_err();
        match err {
            RunError::Failed(f) => {
                assert!(f.failure.contains("out-of-memory"), "{}", f.failure)
            }
            other => panic!("unexpected {other:?}"),
        }
    }

    #[test]
    fn crash_is_reported() {
        let err = Pipeline::new(sh("kill -SEGV $$")).run().unwrap_err();
        assert!(matches!(&err, RunError::Failed(f) if f.failure.contains("SIGSEGV")));
    }

    #[test]
    fn timeout_kills_the_whole_process_group() {
        let started = Instant::now();
        // the inner sleep is a grandchild: it must die with the group
        let err = Pipeline::new(sh("sleep 30; echo done"))
            .timeout(Some(Duration::from_millis(300)))
            .run()
            .unwrap_err();
        assert!(err.is_timeout(), "{err}");
        assert!(started.elapsed() < Duration::from_secs(10));
        assert!(err.to_string().contains("timed out after 0s"));
    }

    #[test]
    fn missing_program_is_a_spawn_error() {
        let err = Pipeline::new(Cmd::new("/nonexistent/tool"))
            .run()
            .unwrap_err();
        assert!(matches!(err, RunError::Spawn { .. }), "{err}");
    }

    #[test]
    fn missing_input_file_is_an_error() {
        let err = Pipeline::new(Cmd::new("cat"))
            .stdin(Input::File("/nonexistent/input".into()))
            .run()
            .unwrap_err();
        assert!(matches!(err, RunError::Redirect { .. }), "{err}");
    }

    #[test]
    fn usage_is_measured() {
        let report = Pipeline::new(sh("i=0; while [ $i -lt 300000 ]; do i=$((i+1)); done"))
            .run()
            .unwrap();
        assert!(report.cpu() > Duration::from_millis(50));
        assert!(report.max_rss_kb() > 0);
    }

    #[test]
    fn stderr_is_shared_by_all_processes() {
        let dir = tempfile::tempdir().unwrap();
        let log = dir.path().join("log");
        Pipeline::new(sh("echo one >&2; echo x"))
            .pipe(sh("cat >/dev/null; echo two >&2"))
            .stderr(Output::Append(log.clone()))
            .run()
            .unwrap();
        let text = std::fs::read_to_string(log).unwrap();
        assert!(text.contains("one") && text.contains("two"), "{text}");
    }

    #[test]
    fn fed_input_and_per_command_stderr() {
        let dir = tempfile::tempdir().unwrap();
        let (out, first_err, last_err) = (
            dir.path().join("out"),
            dir.path().join("e1"),
            dir.path().join("e2"),
        );
        Pipeline::new(sh("cat; echo first >&2").stderr(Output::File(first_err.clone())))
            .pipe(sh("tr a-z A-Z; echo last >&2"))
            .stdout(Output::File(out.clone()))
            .stderr(Output::File(last_err.clone()))
            .run_fed(|w| {
                use std::io::Write;
                w.write_all(b"hello\n")
            })
            .unwrap();
        assert_eq!(std::fs::read_to_string(out).unwrap(), "HELLO\n");
        assert_eq!(std::fs::read_to_string(first_err).unwrap(), "first\n");
        assert_eq!(std::fs::read_to_string(last_err).unwrap(), "last\n");
    }

    #[test]
    fn failing_feeder_fails_the_run() {
        let err = Pipeline::new(Cmd::new("cat"))
            .run_fed(|w| {
                use std::io::Write;
                w.write_all(b"partial")?;
                Err(io::Error::other("truncated gzip"))
            })
            .unwrap_err();
        assert!(matches!(err, RunError::Input { .. }), "{err}");
        assert!(err.to_string().contains("truncated gzip"));
    }

    #[test]
    fn display_quotes_arguments() {
        let p = Pipeline::new(Cmd::new("bwa").args(["mem", "a b.fa"]))
            .pipe(Cmd::new("samtools").args(["view", "-"]))
            .stdout(Output::File("o.sam".into()));
        assert_eq!(p.display(), "bwa mem 'a b.fa' | samtools view - > o.sam");
    }

    #[test]
    fn humanize_durations() {
        assert_eq!(humanize(Duration::from_secs(300)), "5m");
        assert_eq!(humanize(Duration::from_secs(3700)), "1h 1m");
        assert_eq!(humanize(Duration::from_secs(42)), "42s");
    }
}
