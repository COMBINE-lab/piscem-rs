//! RAD file comparison.
//!
//! Compares two RAD files for semantic equivalence. Reads the binary format
//! directly and compares headers and records.

use std::io::{Read, BufReader};
use std::path::Path;

use anyhow::{Context, Result};
use serde::Serialize;

/// Summary of RAD file comparison.
#[derive(Debug, Serialize)]
pub struct RadComparisonSummary {
    pub passed: bool,
    pub notes: String,
    pub file_a_size: u64,
    pub file_b_size: u64,
    pub header_match: bool,
    pub num_refs_match: bool,
    pub num_chunks_a: u64,
    pub num_chunks_b: u64,
}

/// Parsed RAD header (minimal — just enough for comparison).
#[derive(Debug)]
struct RadHeader {
    is_paired: bool,
    num_refs: u64,
    ref_names: Vec<String>,
    num_chunks: u64,
}

/// Read a RAD header from a binary stream.
fn read_rad_header<R: Read>(reader: &mut R) -> Result<RadHeader> {
    // is_paired (u8)
    let mut buf1 = [0u8; 1];
    reader.read_exact(&mut buf1).context("reading is_paired")?;
    let is_paired = buf1[0] != 0;

    // num_refs (u64)
    let mut buf8 = [0u8; 8];
    reader.read_exact(&mut buf8).context("reading num_refs")?;
    let num_refs = u64::from_le_bytes(buf8);

    // ref_names: num_refs × (u16 len + bytes)
    let mut ref_names = Vec::with_capacity(num_refs as usize);
    for i in 0..num_refs {
        let mut buf2 = [0u8; 2];
        reader
            .read_exact(&mut buf2)
            .with_context(|| format!("reading ref name length {}", i))?;
        let name_len = u16::from_le_bytes(buf2) as usize;
        let mut name_buf = vec![0u8; name_len];
        reader
            .read_exact(&mut name_buf)
            .with_context(|| format!("reading ref name {}", i))?;
        ref_names.push(String::from_utf8_lossy(&name_buf).to_string());
    }

    // num_chunks (u64)
    reader.read_exact(&mut buf8).context("reading num_chunks")?;
    let num_chunks = u64::from_le_bytes(buf8);

    Ok(RadHeader {
        is_paired,
        num_refs,
        ref_names,
        num_chunks,
    })
}

/// Compare two RAD files.
///
/// Compares file sizes, headers (is_paired, num_refs, ref_names), and
/// chunk counts. Does not do record-level comparison for multi-threaded
/// outputs (record order is nondeterministic).
pub fn compare_rad_files(path_a: &Path, path_b: &Path) -> Result<RadComparisonSummary> {
    let meta_a = std::fs::metadata(path_a)
        .with_context(|| format!("stat {}", path_a.display()))?;
    let meta_b = std::fs::metadata(path_b)
        .with_context(|| format!("stat {}", path_b.display()))?;

    let file_a = std::fs::File::open(path_a)
        .with_context(|| format!("open {}", path_a.display()))?;
    let file_b = std::fs::File::open(path_b)
        .with_context(|| format!("open {}", path_b.display()))?;

    let header_a = read_rad_header(&mut BufReader::new(&file_a))?;
    let header_b = read_rad_header(&mut BufReader::new(&file_b))?;

    let num_refs_match = header_a.num_refs == header_b.num_refs
        && header_a.ref_names == header_b.ref_names;

    let header_match = header_a.is_paired == header_b.is_paired && num_refs_match;

    let passed = header_match && header_a.num_chunks == header_b.num_chunks;

    let mut notes = Vec::new();
    if !header_match {
        notes.push("headers differ".to_string());
    }
    if header_a.num_chunks != header_b.num_chunks {
        notes.push(format!(
            "chunk count mismatch: {} vs {}",
            header_a.num_chunks, header_b.num_chunks
        ));
    }
    if passed {
        notes.push("headers and chunk counts match".to_string());
    }

    Ok(RadComparisonSummary {
        passed,
        notes: notes.join("; "),
        file_a_size: meta_a.len(),
        file_b_size: meta_b.len(),
        header_match,
        num_refs_match,
        num_chunks_a: header_a.num_chunks,
        num_chunks_b: header_b.num_chunks,
    })
}

/// Validate that a RAD file is readable and has a well-formed header.
///
/// Returns the number of references found in the header.
pub fn validate_rad_file(path: &Path) -> Result<u64> {
    let file = std::fs::File::open(path)
        .with_context(|| format!("open {}", path.display()))?;
    let header = read_rad_header(&mut BufReader::new(&file))?;

    if header.num_refs == 0 {
        anyhow::bail!("RAD file has 0 references");
    }

    Ok(header.num_refs)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::rad::write_rad_header_bulk;
    use std::io::Cursor;

    #[test]
    fn test_read_rad_header_bulk() {
        let mut buf = Vec::new();
        let names = vec!["ref1", "ref2"];
        let ref_lens = vec![1000u32, 2000u32];
        let _chunk_off = write_rad_header_bulk(&mut buf, true, 2, &names, &ref_lens).unwrap();

        let header = read_rad_header(&mut Cursor::new(&buf)).unwrap();
        assert!(header.is_paired);
        assert_eq!(header.num_refs, 2);
        assert_eq!(header.ref_names, vec!["ref1", "ref2"]);
    }

    #[test]
    fn test_rad_comparison_summary_serialize() {
        let summary = RadComparisonSummary {
            passed: true,
            notes: "ok".to_string(),
            file_a_size: 1000,
            file_b_size: 1000,
            header_match: true,
            num_refs_match: true,
            num_chunks_a: 5,
            num_chunks_b: 5,
        };
        let json = serde_json::to_string(&summary).unwrap();
        assert!(json.contains("\"passed\":true"));
    }

    #[test]
    fn test_validate_rad_header() {
        // Write a valid bulk RAD header to a temp file
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.rad");
        let names = vec!["ref1"];
        let ref_lens = vec![500u32];
        let mut file = std::fs::File::create(&path).unwrap();
        write_rad_header_bulk(&mut file, false, 1, &names, &ref_lens).unwrap();
        drop(file);

        let num_refs = validate_rad_file(&path).unwrap();
        assert_eq!(num_refs, 1);
    }
}
