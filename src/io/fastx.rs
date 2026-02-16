//! FASTX reader helpers — open and concatenate FASTQ files with decompression.
//!
//! Provides `open_concatenated_readers()` for creating a single `Read` stream
//! from multiple (possibly compressed) FASTQ files, suitable for feeding to
//! paraseq's `Reader`.

use anyhow::{Context, Result};

// Re-export paraseq types used by the mapping pipeline.
pub use paraseq::fastq;
pub use paraseq::parallel::{
    MultiParallelProcessor, PairedParallelProcessor, ParallelProcessor, ParallelReader,
};
pub use paraseq::Record;

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
pub fn open_concatenated_readers(
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
    #[test]
    fn test_open_concatenated_readers_empty() {
        let result = super::open_concatenated_readers(&[]);
        assert!(result.is_err());
    }
}
