//! FASTA / FASTQ reading, plain or gzip-compressed.
//!
//! Compression is detected from the gzip magic bytes, not the file name,
//! and FASTA vs FASTQ from the first character of the decompressed
//! stream. Multi-member gzip (e.g. bgzip, or concatenated `.gz` files) is
//! read in full. A truncated or corrupt file is an error, never a shorter
//! set of records.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

use flate2::read::MultiGzDecoder;

const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// Errors while reading a sequence file.
#[derive(Debug, thiserror::Error)]
pub enum SeqError {
    #[error("cannot read {path}: {source}")]
    Io { path: PathBuf, source: io::Error },
    #[error("{path}: not a FASTA or FASTQ file (first character {found:?})")]
    UnknownFormat { path: PathBuf, found: char },
    #[error("{path}, line {line}: {message}")]
    Format {
        path: PathBuf,
        line: u64,
        message: String,
    },
}

/// FASTA or FASTQ.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Fasta,
    Fastq,
}

/// One sequence record.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Record {
    /// The header line without its leading `>` / `@`.
    pub header: String,
    pub seq: Vec<u8>,
    /// Quality string (FASTQ only).
    pub qual: Option<Vec<u8>>,
}

impl Record {
    /// The record name: the header up to the first whitespace.
    pub fn name(&self) -> &str {
        self.header.split_whitespace().next().unwrap_or("")
    }
}

/// Returns true if the file starts with the gzip magic bytes.
pub fn is_gzip(path: &Path) -> Result<bool, SeqError> {
    let io_err = |source| SeqError::Io {
        path: path.to_path_buf(),
        source,
    };
    let mut magic = [0u8; 2];
    let mut file = File::open(path).map_err(io_err)?;
    let mut read = 0;
    while read < 2 {
        let n = file.read(&mut magic[read..]).map_err(io_err)?;
        if n == 0 {
            return Ok(false);
        }
        read += n;
    }
    Ok(magic == GZIP_MAGIC)
}

/// Opens a file for reading, decompressing it if it is gzip.
pub fn open_maybe_gzip(path: &Path) -> Result<Box<dyn BufRead + Send>, SeqError> {
    let io_err = |source| SeqError::Io {
        path: path.to_path_buf(),
        source,
    };
    let file = File::open(path).map_err(io_err)?;
    if is_gzip(path)? {
        Ok(Box::new(BufReader::with_capacity(
            1 << 16,
            MultiGzDecoder::new(file),
        )))
    } else {
        Ok(Box::new(BufReader::with_capacity(1 << 16, file)))
    }
}

/// Reads FASTA or FASTQ records from a file (plain or gzip).
pub struct SeqReader {
    path: PathBuf,
    input: Box<dyn BufRead + Send>,
    format: Format,
    line_no: u64,
    /// A FASTA header already read while finishing the previous record.
    pending_header: Option<String>,
    done: bool,
}

impl SeqReader {
    /// Opens `path` and detects its compression and format.
    ///
    /// An empty file is accepted and yields no records.
    pub fn open(path: impl AsRef<Path>) -> Result<Self, SeqError> {
        let path = path.as_ref().to_path_buf();
        let mut input = open_maybe_gzip(&path)?;
        let first = input
            .fill_buf()
            .map_err(|source| SeqError::Io {
                path: path.clone(),
                source,
            })?
            .first()
            .copied();
        let (format, done) = match first {
            None => (Format::Fasta, true),
            Some(b'>') => (Format::Fasta, false),
            Some(b'@') => (Format::Fastq, false),
            Some(other) => {
                return Err(SeqError::UnknownFormat {
                    path,
                    found: other as char,
                });
            }
        };
        Ok(Self {
            path,
            input,
            format,
            line_no: 0,
            pending_header: None,
            done,
        })
    }

    /// The detected format.
    pub fn format(&self) -> Format {
        self.format
    }

    fn read_line(&mut self, buf: &mut String) -> Result<bool, SeqError> {
        buf.clear();
        let n = self.input.read_line(buf).map_err(|source| SeqError::Io {
            path: self.path.clone(),
            source,
        })?;
        if n == 0 {
            return Ok(false);
        }
        self.line_no += 1;
        while buf.ends_with('\n') || buf.ends_with('\r') {
            buf.pop();
        }
        Ok(true)
    }

    fn format_error(&self, message: impl Into<String>) -> SeqError {
        SeqError::Format {
            path: self.path.clone(),
            line: self.line_no,
            message: message.into(),
        }
    }

    fn next_fasta(&mut self) -> Result<Option<Record>, SeqError> {
        let mut line = String::new();
        let header = match self.pending_header.take() {
            Some(header) => header,
            None => loop {
                if !self.read_line(&mut line)? {
                    return Ok(None);
                }
                if line.is_empty() {
                    continue;
                }
                match line.strip_prefix('>') {
                    Some(header) => break header.to_string(),
                    None => return Err(self.format_error("expected a FASTA header ('>')")),
                }
            },
        };
        let mut seq = Vec::new();
        while self.read_line(&mut line)? {
            if let Some(next) = line.strip_prefix('>') {
                self.pending_header = Some(next.to_string());
                break;
            }
            seq.extend_from_slice(line.trim_end().as_bytes());
        }
        Ok(Some(Record {
            header,
            seq,
            qual: None,
        }))
    }

    fn next_fastq(&mut self) -> Result<Option<Record>, SeqError> {
        let mut line = String::new();
        loop {
            if !self.read_line(&mut line)? {
                return Ok(None);
            }
            if !line.is_empty() {
                break;
            }
        }
        let header = match line.strip_prefix('@') {
            Some(header) => header.to_string(),
            None => return Err(self.format_error("expected a FASTQ header ('@')")),
        };
        if !self.read_line(&mut line)? {
            return Err(self.format_error("truncated FASTQ record (no sequence line)"));
        }
        let seq = line.as_bytes().to_vec();
        if !self.read_line(&mut line)? || !line.starts_with('+') {
            return Err(self.format_error("truncated FASTQ record (no '+' line)"));
        }
        if !self.read_line(&mut line)? {
            return Err(self.format_error("truncated FASTQ record (no quality line)"));
        }
        if line.len() != seq.len() {
            return Err(self.format_error(format!(
                "quality length {} differs from sequence length {}",
                line.len(),
                seq.len()
            )));
        }
        Ok(Some(Record {
            header,
            seq,
            qual: Some(line.as_bytes().to_vec()),
        }))
    }
}

impl Iterator for SeqReader {
    type Item = Result<Record, SeqError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        let record = match self.format {
            Format::Fasta => self.next_fasta(),
            Format::Fastq => self.next_fastq(),
        };
        match record {
            Ok(Some(record)) => Some(Ok(record)),
            Ok(None) => {
                self.done = true;
                None
            }
            Err(err) => {
                self.done = true;
                Some(Err(err))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::io::Write;

    const FASTQ: &str = "@r1 desc\nACGT\n+\nIIII\n@r2\nGG\n+r2\nII\n";
    const FASTA: &str = ">c1 len=6\nACG\nTAA\n\n>c2\nGGG\n";

    fn write(dir: &tempfile::TempDir, name: &str, data: &[u8]) -> PathBuf {
        let path = dir.path().join(name);
        std::fs::write(&path, data).unwrap();
        path
    }

    fn gzip(data: &[u8]) -> Vec<u8> {
        let mut enc = GzEncoder::new(Vec::new(), Compression::default());
        enc.write_all(data).unwrap();
        enc.finish().unwrap()
    }

    fn read_all(path: &Path) -> Result<Vec<Record>, SeqError> {
        SeqReader::open(path)?.collect()
    }

    #[test]
    fn plain_and_gzip_fastq_give_the_same_records() {
        let dir = tempfile::tempdir().unwrap();
        let plain = read_all(&write(&dir, "r.fq", FASTQ.as_bytes())).unwrap();
        // the name doesn't matter, only the content
        let gz = read_all(&write(&dir, "r.fq", &gzip(FASTQ.as_bytes()))).unwrap();
        assert_eq!(plain, gz);
        assert_eq!(plain.len(), 2);
        assert_eq!(plain[0].name(), "r1");
        assert_eq!(plain[0].seq, b"ACGT");
        assert_eq!(plain[1].qual.as_deref(), Some(&b"II"[..]));
    }

    #[test]
    fn fasta_multiline_records() {
        let dir = tempfile::tempdir().unwrap();
        let path = write(&dir, "c.fa.gz", &gzip(FASTA.as_bytes()));
        let reader = SeqReader::open(&path).unwrap();
        assert_eq!(reader.format(), Format::Fasta);
        let records: Vec<_> = reader.collect::<Result<_, _>>().unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].seq, b"ACGTAA");
        assert_eq!(records[1].header, "c2");
    }

    #[test]
    fn multi_member_gzip_is_read_in_full() {
        let dir = tempfile::tempdir().unwrap();
        let mut data = gzip(b"@a\nAC\n+\nII\n");
        data.extend(gzip(b"@b\nGT\n+\nII\n"));
        let records = read_all(&write(&dir, "m.fq.gz", &data)).unwrap();
        assert_eq!(
            records.iter().map(|r| r.name()).collect::<Vec<_>>(),
            ["a", "b"]
        );
    }

    #[test]
    fn truncated_gzip_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let data = gzip(FASTQ.repeat(1000).as_bytes());
        let path = write(&dir, "t.fq.gz", &data[..data.len() / 2]);
        assert!(read_all(&path).is_err());
    }

    #[test]
    fn truncated_fastq_record_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let path = write(&dir, "t.fq", b"@a\nACGT\n+\n");
        let err = read_all(&path).unwrap_err();
        assert!(err.to_string().contains("truncated"), "{err}");
    }

    #[test]
    fn unknown_format_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let path = write(&dir, "x.txt", b"hello\n");
        assert!(matches!(
            SeqReader::open(&path),
            Err(SeqError::UnknownFormat { .. })
        ));
    }

    #[test]
    fn empty_file_has_no_records() {
        let dir = tempfile::tempdir().unwrap();
        assert!(read_all(&write(&dir, "e.fa", b"")).unwrap().is_empty());
    }
}
