//! The large intermediate text files of a run (the long-read mapping and
//! its edges: tens of GB on large samples) are kept gzipped, as
//! `<name>.gz`. Writing compresses chunks in parallel, each one a gzip
//! member, written in order: a standard multi-member gzip file (`zcat`
//! reads it). Reading takes `<name>.gz`, or the plain `<name>` written by
//! OPERA-MS 2.0.0 (a run resumed with a newer version).

use std::fs::{self, File};
use std::io::{BufWriter, Read, Write};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::Compression;
use flate2::write::GzEncoder;

/// Uncompressed bytes per gzip member.
const CHUNK: usize = 32 << 20;

/// The path of the gzipped file of `path`.
pub fn gz_path(path: &str) -> String {
    format!("{path}.gz")
}

/// Writes `text` to `<path>.gz` with `threads` threads (gzip level 6:
/// about 3 times smaller on these files, against 2 at the fastest level),
/// and removes a plain `path` left by an earlier version.
pub fn write(path: &str, text: &str, threads: usize) -> Result<()> {
    write_chunks(path, text, threads, CHUNK)
}

fn write_chunks(path: &str, text: &str, threads: usize, chunk: usize) -> Result<()> {
    let out = gz_path(path);
    let file = File::create(&out).with_context(|| format!("cannot write {out}"))?;
    let mut file = BufWriter::with_capacity(1 << 20, file);
    let bytes = text.as_bytes();
    let chunks: Vec<&[u8]> = if bytes.is_empty() {
        // one empty member: a valid gzip file
        vec![bytes]
    } else {
        bytes.chunks(chunk.max(1)).collect()
    };
    // a wave of `threads` chunks at a time, written in order
    for wave in chunks.chunks(threads.max(1)) {
        let compressed: Vec<std::io::Result<Vec<u8>>> = std::thread::scope(|scope| {
            let handles: Vec<_> = wave
                .iter()
                .map(|chunk| scope.spawn(move || compress(chunk)))
                .collect();
            handles
                .into_iter()
                .map(|h| {
                    h.join()
                        .unwrap_or_else(|_| Err(std::io::Error::other("gzip thread panicked")))
                })
                .collect()
        });
        for member in compressed {
            let member = member.with_context(|| format!("cannot compress {out}"))?;
            file.write_all(&member)
                .with_context(|| format!("cannot write {out}"))?;
        }
    }
    file.flush()
        .with_context(|| format!("cannot write {out}"))?;
    if Path::new(path).exists() {
        fs::remove_file(path).with_context(|| format!("cannot remove {path}"))?;
    }
    Ok(())
}

fn compress(chunk: &[u8]) -> std::io::Result<Vec<u8>> {
    let mut gz = GzEncoder::new(Vec::with_capacity(chunk.len() / 4), Compression::new(6));
    gz.write_all(chunk)?;
    gz.finish()
}

/// Whether `<path>.gz` or the plain `path` exists.
pub fn exists(path: &str) -> bool {
    Path::new(&gz_path(path)).exists() || Path::new(path).exists()
}

/// The text of `<path>.gz`, else of the plain `path`.
pub fn read(path: &str) -> Result<String> {
    let gz = gz_path(path);
    let (source, mut input): (&str, Box<dyn Read>) = if Path::new(&gz).exists() {
        let file = File::open(&gz).with_context(|| format!("cannot read {gz}"))?;
        (
            &gz,
            Box::new(flate2::read::MultiGzDecoder::new(
                std::io::BufReader::with_capacity(1 << 20, file),
            )),
        )
    } else {
        let file = File::open(path).with_context(|| format!("cannot read {path}"))?;
        (path, Box::new(file))
    };
    let mut text = String::new();
    input
        .read_to_string(&mut text)
        .with_context(|| format!("cannot read {source}"))?;
    Ok(text)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip_and_plain_fallback() {
        let dir = tempfile::tempdir().unwrap();
        let path = format!("{}/edges", dir.path().display());
        // several members, in order
        let chunk = 1000;
        let text: String = (0..2000).map(|i| format!("line {i:>12}\n")).collect();
        assert!(text.len() > 30 * chunk);
        // a plain file of an earlier version: read, then replaced
        fs::write(&path, "old\n").unwrap();
        assert!(exists(&path));
        assert_eq!(read(&path).unwrap(), "old\n");
        for threads in [1, 2, 8] {
            write_chunks(&path, &text, threads, chunk).unwrap();
            assert!(!Path::new(&path).exists());
            assert_eq!(read(&path).unwrap(), text);
        }
        // a standard gzip file (when gzip is installed)
        if let Ok(out) = std::process::Command::new("gzip")
            .args(["-dc", &gz_path(&path)])
            .output()
        {
            assert!(out.status.success());
            assert!(out.stdout == text.as_bytes());
        }
        write(&path, "", 4).unwrap();
        assert_eq!(read(&path).unwrap(), "");
        assert!(read(&format!("{path}.missing")).is_err());
    }
}
