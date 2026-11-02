//! Running independent jobs in parallel, at most N at a time.
//!
//! Unlike `xargs`, a failure stops the step: no new job is started
//! after the first failure, and the error of the first failed job (in job
//! order) is returned, with the number of other failures.

use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::thread;

/// Errors of a parallel run.
#[derive(Debug, thiserror::Error)]
pub enum ParallelError<E> {
    #[error("job {index} failed: {error}{}", others_text(*.other_failures))]
    Job {
        index: usize,
        error: E,
        other_failures: usize,
    },
    #[error("a job panicked")]
    Panicked,
}

fn others_text(n: usize) -> String {
    match n {
        0 => String::new(),
        1 => " (1 other job failed too)".into(),
        n => format!(" ({n} other jobs failed too)"),
    }
}

/// Runs `jobs` with at most `max_parallel` of them at the same time and
/// returns their results in job order.
pub fn run_parallel<T, E, F>(jobs: Vec<F>, max_parallel: usize) -> Result<Vec<T>, ParallelError<E>>
where
    T: Send,
    E: Send,
    F: FnOnce() -> Result<T, E> + Send,
{
    let count = jobs.len();
    let jobs: Vec<Mutex<Option<F>>> = jobs.into_iter().map(|j| Mutex::new(Some(j))).collect();
    let results: Vec<Mutex<Option<Result<T, E>>>> = (0..count).map(|_| Mutex::new(None)).collect();
    let next = AtomicUsize::new(0);
    let failed = AtomicBool::new(false);
    let workers = max_parallel.clamp(1, count.max(1));

    let all_joined = thread::scope(|scope| {
        let handles: Vec<_> = (0..workers)
            .map(|_| {
                scope.spawn(|| {
                    while !failed.load(Ordering::SeqCst) {
                        let index = next.fetch_add(1, Ordering::SeqCst);
                        if index >= count {
                            break;
                        }
                        let job = jobs[index].lock().ok().and_then(|mut j| j.take());
                        let Some(job) = job else { continue };
                        let result = job();
                        if result.is_err() {
                            failed.store(true, Ordering::SeqCst);
                        }
                        if let Ok(mut slot) = results[index].lock() {
                            *slot = Some(result);
                        }
                    }
                })
            })
            .collect();
        handles.into_iter().all(|h| h.join().is_ok())
    });
    if !all_joined {
        return Err(ParallelError::Panicked);
    }

    let mut values = Vec::with_capacity(count);
    let mut first_error = None;
    let mut other_failures = 0;
    for (index, slot) in results.into_iter().enumerate() {
        match slot.into_inner().map_err(|_| ParallelError::Panicked)? {
            Some(Ok(value)) => values.push(value),
            Some(Err(error)) => {
                if first_error.is_none() {
                    first_error = Some((index, error));
                } else {
                    other_failures += 1;
                }
            }
            // not started because of an earlier failure
            None => {}
        }
    }
    match first_error {
        Some((index, error)) => Err(ParallelError::Job {
            index,
            error,
            other_failures,
        }),
        None => Ok(values),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{Duration, Instant};

    #[test]
    fn results_in_job_order_and_bounded_parallelism() {
        let running = AtomicUsize::new(0);
        let peak = AtomicUsize::new(0);
        let jobs: Vec<_> = (0..8)
            .map(|i| {
                let (running, peak) = (&running, &peak);
                move || {
                    let now = running.fetch_add(1, Ordering::SeqCst) + 1;
                    peak.fetch_max(now, Ordering::SeqCst);
                    thread::sleep(Duration::from_millis(30));
                    running.fetch_sub(1, Ordering::SeqCst);
                    Ok::<usize, String>(i * 10)
                }
            })
            .collect();
        let started = Instant::now();
        let results = run_parallel(jobs, 3).unwrap();
        assert_eq!(results, (0..8).map(|i| i * 10).collect::<Vec<_>>());
        assert!(peak.load(Ordering::SeqCst) <= 3);
        assert!(started.elapsed() < Duration::from_millis(8 * 30));
    }

    #[test]
    fn first_failure_stops_new_jobs() {
        let started = AtomicUsize::new(0);
        let jobs: Vec<_> = (0..20)
            .map(|i| {
                let started = &started;
                move || {
                    started.fetch_add(1, Ordering::SeqCst);
                    if i == 1 {
                        Err(format!("job {i} broke"))
                    } else {
                        Ok(i)
                    }
                }
            })
            .collect();
        let err = run_parallel(jobs, 1).unwrap_err();
        match err {
            ParallelError::Job {
                index,
                error,
                other_failures,
            } => {
                assert_eq!(
                    (index, error.as_str(), other_failures),
                    (1, "job 1 broke", 0)
                );
            }
            other => panic!("unexpected {other:?}"),
        }
        assert_eq!(started.load(Ordering::SeqCst), 2);
    }

    #[test]
    fn no_jobs() {
        let jobs: Vec<fn() -> Result<(), String>> = Vec::new();
        assert!(run_parallel(jobs, 4).unwrap().is_empty());
    }
}
