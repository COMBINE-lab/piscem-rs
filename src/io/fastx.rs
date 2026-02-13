//! FASTX reader — wraps `paraseq` for batched, parallel-aware reading.
//!
//! Provides types for reading single-end and paired-end FASTA/FASTQ files
//! in chunks suitable for the mapping pipeline.
//!
//! Uses `paraseq`'s `RecordSet` for efficient buffered reading, and
//! re-exports `paraseq` types needed by the parallel pipeline.

use anyhow::{Context, Result};

// Re-export paraseq types used by the mapping pipeline.
pub use paraseq::fastq;
pub use paraseq::parallel::{ParallelProcessor, PairedParallelProcessor};
pub use paraseq::Record;

// ---------------------------------------------------------------------------
// ReadPair
// ---------------------------------------------------------------------------

/// A read pair (or single read) with owned byte buffers.
///
/// For single-end reads, `seq2` and `qual2` are `None`.
#[derive(Debug, Clone)]
pub struct ReadPair {
    pub name: Vec<u8>,
    pub seq1: Vec<u8>,
    pub qual1: Vec<u8>,
    pub seq2: Option<Vec<u8>>,
    pub qual2: Option<Vec<u8>>,
}

/// A chunk of reads for batch processing.
pub type ReadChunk = Vec<ReadPair>;

// ---------------------------------------------------------------------------
// FastxConfig
// ---------------------------------------------------------------------------

/// Configuration for FASTX input sources.
#[derive(Debug, Clone)]
pub struct FastxConfig {
    pub read1_paths: Vec<String>,
    pub read2_paths: Vec<String>,
    pub chunk_size: usize,
}

impl Default for FastxConfig {
    fn default() -> Self {
        Self {
            read1_paths: Vec::new(),
            read2_paths: Vec::new(),
            chunk_size: 1000,
        }
    }
}

impl FastxConfig {
    /// Whether this is a paired-end configuration.
    pub fn is_paired(&self) -> bool {
        !self.read2_paths.is_empty()
    }
}

// ---------------------------------------------------------------------------
// FastxSource
// ---------------------------------------------------------------------------

/// Sequential FASTX reader that produces chunks of `ReadPair`.
///
/// Wraps `paraseq::fastq::Reader` for single-end and paired-end reading.
/// For parallel dispatch, use the `paraseq` `ParallelProcessor` /
/// `PairedParallelProcessor` traits directly instead.
pub struct FastxSource {
    config: FastxConfig,
    reader1: fastq::Reader<Box<dyn std::io::Read + Send>>,
    reader2: Option<fastq::Reader<Box<dyn std::io::Read + Send>>>,
    record_set1: fastq::RecordSet,
    record_set2: Option<fastq::RecordSet>,
}

impl FastxSource {
    /// Open FASTX files from the configuration.
    pub fn new(config: FastxConfig) -> Result<Self> {
        let r1 = open_concatenated_readers(&config.read1_paths)?;
        let reader1 = fastq::Reader::new(r1);
        let record_set1 = reader1.new_record_set();

        let (reader2, record_set2) = if config.is_paired() {
            let r2 = open_concatenated_readers(&config.read2_paths)?;
            let rdr2 = fastq::Reader::new(r2);
            let rs2 = rdr2.new_record_set();
            (Some(rdr2), Some(rs2))
        } else {
            (None, None)
        };

        Ok(Self {
            config,
            reader1,
            reader2,
            record_set1,
            record_set2,
        })
    }

    /// Whether this source provides paired-end reads.
    pub fn is_paired(&self) -> bool {
        self.config.is_paired()
    }

    /// Read the next chunk of reads into the provided buffer.
    ///
    /// Returns `Ok(true)` if reads were produced, `Ok(false)` at EOF.
    pub fn next_chunk(&mut self, chunk: &mut ReadChunk) -> Result<bool> {
        chunk.clear();

        if self.is_paired() {
            let rs2 = self.record_set2.as_mut().unwrap();
            let rdr2 = self.reader2.as_mut().unwrap();

            let has1 = self.record_set1.fill(&mut self.reader1)?;
            let has2 = rs2.fill(rdr2)?;
            if !has1 || !has2 {
                return Ok(false);
            }

            let mut iter1 = self.record_set1.iter();
            let mut iter2 = rs2.iter();

            loop {
                let r1 = iter1.next();
                let r2 = iter2.next();
                match (r1, r2) {
                    (Some(Ok(rec1)), Some(Ok(rec2))) => {
                        chunk.push(ReadPair {
                            name: rec1.id().to_vec(),
                            seq1: rec1.seq().into_owned(),
                            qual1: rec1.qual().map(|q| q.to_vec()).unwrap_or_default(),
                            seq2: Some(rec2.seq().into_owned()),
                            qual2: Some(rec2.qual().map(|q| q.to_vec()).unwrap_or_default()),
                        });
                    }
                    (None, None) => break,
                    _ => break,
                }
                if chunk.len() >= self.config.chunk_size {
                    break;
                }
            }
        } else {
            let has = self.record_set1.fill(&mut self.reader1)?;
            if !has {
                return Ok(false);
            }

            for rec in self.record_set1.iter() {
                let rec = rec?;
                chunk.push(ReadPair {
                    name: rec.id().to_vec(),
                    seq1: rec.seq().into_owned(),
                    qual1: rec.qual().map(|q| q.to_vec()).unwrap_or_default(),
                    seq2: None,
                    qual2: None,
                });
                if chunk.len() >= self.config.chunk_size {
                    break;
                }
            }
        }

        Ok(!chunk.is_empty())
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Open a single file with automatic decompression (gzip, zstd, etc.).
fn open_with_decompression(path: &str) -> Result<Box<dyn std::io::Read + Send>> {
    let (reader, _format) = niffler::send::from_path(path)
        .with_context(|| format!("failed to open {}", path))?;
    Ok(reader)
}

/// Open multiple files and concatenate them into a single reader.
///
/// Automatically detects and decompresses gzip, zstd, and other formats
/// via niffler.
fn open_concatenated_readers(
    paths: &[String],
) -> Result<Box<dyn std::io::Read + Send>> {
    use std::io::Read;

    if paths.is_empty() {
        anyhow::bail!("no input files specified");
    }
    if paths.len() == 1 {
        return open_with_decompression(&paths[0]);
    }
    // Collect all file readers, then chain them via a wrapper.
    let mut readers: Vec<Box<dyn Read + Send>> = Vec::with_capacity(paths.len());
    for path in paths {
        readers.push(open_with_decompression(path)?);
    }
    Ok(Box::new(MultiReader { readers, current: 0 }))
}

/// Concatenating reader over multiple boxed readers.
struct MultiReader {
    readers: Vec<Box<dyn std::io::Read + Send>>,
    current: usize,
}

impl std::io::Read for MultiReader {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        while self.current < self.readers.len() {
            let n = self.readers[self.current].read(buf)?;
            if n > 0 {
                return Ok(n);
            }
            self.current += 1;
        }
        Ok(0)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fastx_config_defaults() {
        let config = FastxConfig::default();
        assert!(!config.is_paired());
        assert_eq!(config.chunk_size, 1000);
    }

    #[test]
    fn test_fastx_config_paired() {
        let config = FastxConfig {
            read1_paths: vec!["r1.fq".to_string()],
            read2_paths: vec!["r2.fq".to_string()],
            chunk_size: 500,
        };
        assert!(config.is_paired());
    }

    #[test]
    fn test_read_pair_single_end() {
        let rp = ReadPair {
            name: b"read1".to_vec(),
            seq1: b"ACGT".to_vec(),
            qual1: b"IIII".to_vec(),
            seq2: None,
            qual2: None,
        };
        assert!(rp.seq2.is_none());
    }
}
