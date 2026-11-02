//! Safe wrappers around libc process functions.
//!
//! The intended use (by `opera-ms-core`'s runner) avoids PID-reuse races:
//! wait for a child to exit with [`wait_exit_nowait`] (the child stays a
//! zombie, so its PID and process group can't be reused), send signals
//! with [`kill_process_group`] only before reaping, then reap it with
//! [`reap`] to get its exit status and resource usage.

use std::io;
use std::time::Duration;

/// How a process ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Exit {
    /// Normal exit with this exit code.
    Code(i32),
    /// Killed by this signal number.
    Signal(i32),
}

impl Exit {
    /// True for a normal exit with code 0.
    pub fn success(self) -> bool {
        self == Exit::Code(0)
    }
}

/// Resource usage of a reaped child (and of its own children that it
/// waited for), as returned by `wait4`.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Usage {
    /// CPU time in user mode.
    pub user: Duration,
    /// CPU time in kernel mode.
    pub system: Duration,
    /// Peak resident set size in kilobytes (largest single process).
    pub max_rss_kb: u64,
}

fn timeval(tv: libc::timeval) -> Duration {
    Duration::from_secs(tv.tv_sec.max(0) as u64) + Duration::from_micros(tv.tv_usec.max(0) as u64)
}

fn usage(ru: &libc::rusage) -> Usage {
    Usage {
        user: timeval(ru.ru_utime),
        system: timeval(ru.ru_stime),
        max_rss_kb: ru.ru_maxrss.max(0) as u64,
    }
}

fn to_pid(pid: u32) -> io::Result<libc::pid_t> {
    libc::pid_t::try_from(pid)
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, format!("invalid pid {pid}")))
}

/// Blocks until the child `pid` has exited, without reaping it, and
/// returns how it ended.
///
/// After this returns `Ok`, the child is a zombie: its PID (and process
/// group, if it leads one) stays reserved until [`reap`] is called.
pub fn wait_exit_nowait(pid: u32) -> io::Result<Exit> {
    let pid = to_pid(pid)?;
    loop {
        // SAFETY: `info` is a valid, zero-initialised siginfo_t that waitid
        // fills in; P_PID with a positive pid only waits for that child.
        let (ret, info) = unsafe {
            let mut info: libc::siginfo_t = std::mem::zeroed();
            let ret = libc::waitid(
                libc::P_PID,
                pid as libc::id_t,
                &mut info,
                libc::WEXITED | libc::WNOWAIT,
            );
            (ret, info)
        };
        if ret == 0 {
            // SAFETY: waitid succeeded for an exited child, so si_status is
            // the field set by the kernel for SIGCHLD-style info.
            let status = unsafe { info.si_status() };
            return Ok(match info.si_code {
                libc::CLD_EXITED => Exit::Code(status),
                _ => Exit::Signal(status),
            });
        }
        let err = io::Error::last_os_error();
        if err.kind() != io::ErrorKind::Interrupted {
            return Err(err);
        }
    }
}

/// Reaps the child `pid` (blocking until it exits) and returns how it
/// ended and its resource usage.
pub fn reap(pid: u32) -> io::Result<(Exit, Usage)> {
    let pid = to_pid(pid)?;
    loop {
        let mut status: libc::c_int = 0;
        // SAFETY: `ru` is zero-initialised and filled in by wait4; `status`
        // is a valid pointer; a positive pid only waits for that child.
        let (ret, ru) = unsafe {
            let mut ru: libc::rusage = std::mem::zeroed();
            let ret = libc::wait4(pid, &mut status, 0, &mut ru);
            (ret, ru)
        };
        if ret == pid {
            let exit = if libc::WIFSIGNALED(status) {
                Exit::Signal(libc::WTERMSIG(status))
            } else {
                Exit::Code(libc::WEXITSTATUS(status))
            };
            return Ok((exit, usage(&ru)));
        }
        let err = io::Error::last_os_error();
        if err.kind() != io::ErrorKind::Interrupted {
            return Err(err);
        }
    }
}

/// Sends `signal` to every process in the process group `pgid`.
///
/// Only call this for a group whose leader hasn't been reaped yet (see the
/// module documentation), so that the group ID can't have been reused.
pub fn kill_process_group(pgid: u32, signal: i32) -> io::Result<()> {
    let pgid = to_pid(pgid)?;
    if pgid <= 1 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "refusing to signal pgid <= 1",
        ));
    }
    // SAFETY: kill() with a negative pid signals that process group; it
    // has no memory-safety requirements.
    let ret = unsafe { libc::kill(-pgid, signal) };
    if ret == 0 {
        return Ok(());
    }
    let err = io::Error::last_os_error();
    // ESRCH: every process of the group has already exited.
    if err.raw_os_error() == Some(libc::ESRCH) {
        Ok(())
    } else {
        Err(err)
    }
}

/// CPU time used so far by the calling process itself (all threads).
pub fn self_usage() -> io::Result<Usage> {
    // SAFETY: `ru` is zero-initialised and filled in by getrusage.
    let (ret, ru) = unsafe {
        let mut ru: libc::rusage = std::mem::zeroed();
        let ret = libc::getrusage(libc::RUSAGE_SELF, &mut ru);
        (ret, ru)
    };
    if ret == 0 {
        Ok(usage(&ru))
    } else {
        Err(io::Error::last_os_error())
    }
}

/// Signal numbers used by the runner.
pub mod signal {
    pub const SIGKILL: i32 = libc::SIGKILL;
    pub const SIGTERM: i32 = libc::SIGTERM;
    pub const SIGPIPE: i32 = libc::SIGPIPE;
    pub const SIGSEGV: i32 = libc::SIGSEGV;
    pub const SIGABRT: i32 = libc::SIGABRT;

    /// Short name of a signal ("SIGKILL"), or its number.
    pub fn name(signal: i32) -> String {
        match signal {
            libc::SIGHUP => "SIGHUP".into(),
            libc::SIGINT => "SIGINT".into(),
            libc::SIGQUIT => "SIGQUIT".into(),
            libc::SIGILL => "SIGILL".into(),
            libc::SIGABRT => "SIGABRT".into(),
            libc::SIGBUS => "SIGBUS".into(),
            libc::SIGFPE => "SIGFPE".into(),
            libc::SIGKILL => "SIGKILL".into(),
            libc::SIGSEGV => "SIGSEGV".into(),
            libc::SIGPIPE => "SIGPIPE".into(),
            libc::SIGALRM => "SIGALRM".into(),
            libc::SIGTERM => "SIGTERM".into(),
            libc::SIGXCPU => "SIGXCPU".into(),
            libc::SIGXFSZ => "SIGXFSZ".into(),
            other => format!("signal {other}"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::os::unix::process::CommandExt;
    use std::process::Command;

    #[test]
    #[allow(clippy::zombie_processes)] // reaped with reap(), not Child::wait()
    fn reap_reports_exit_code_and_cpu_time() {
        let child = Command::new("sh")
            .args([
                "-c",
                "i=0; while [ $i -lt 200000 ]; do i=$((i+1)); done; exit 3",
            ])
            .spawn()
            .unwrap();
        assert_eq!(wait_exit_nowait(child.id()).unwrap(), Exit::Code(3));
        let (exit, usage) = reap(child.id()).unwrap();
        assert_eq!(exit, Exit::Code(3));
        assert!(usage.user + usage.system > Duration::from_millis(10));
        assert!(usage.max_rss_kb > 0);
    }

    #[test]
    #[allow(clippy::zombie_processes)] // reaped with reap(), not Child::wait()
    fn killed_process_group_reports_signal() {
        let child = Command::new("sleep")
            .arg("30")
            .process_group(0)
            .spawn()
            .unwrap();
        kill_process_group(child.id(), signal::SIGKILL).unwrap();
        assert_eq!(
            wait_exit_nowait(child.id()).unwrap(),
            Exit::Signal(signal::SIGKILL)
        );
        let (exit, _) = reap(child.id()).unwrap();
        assert_eq!(exit, Exit::Signal(signal::SIGKILL));
        assert_eq!(signal::name(signal::SIGKILL), "SIGKILL");
    }

    #[test]
    fn self_usage_is_available() {
        assert!(self_usage().is_ok());
    }
}
