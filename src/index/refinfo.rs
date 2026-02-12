//! Reference sequence metadata — names and lengths.
//!
//! Stores the name and length of every reference sequence (transcript) in the
//! index. This corresponds to the C++ `m_ref_names` / `m_ref_lens` vectors
//! stored in the `.refinfo` file.

use anyhow::{Context, Result, bail};
use std::io::{Read, Write, BufWriter, BufReader};

/// Magic bytes for the piscem-rs refinfo file format.
const REFINFO_MAGIC: &[u8; 8] = b"PRFINF01";

/// Reference sequence metadata.
///
/// Contains the name and length of each reference sequence in the index,
/// stored in reference-ID order (i.e. `ref_names[i]` and `ref_lens[i]`
/// describe reference `i`).
#[derive(Debug, Clone)]
pub struct RefInfo {
    /// Reference sequence names, one per reference.
    ref_names: Vec<String>,
    /// Reference sequence lengths (in bases), one per reference.
    ref_lens: Vec<u64>,
}

impl RefInfo {
    /// Create a new `RefInfo` from parallel name and length vectors.
    ///
    /// # Panics
    /// Panics if `names.len() != lens.len()`.
    pub fn new(names: Vec<String>, lens: Vec<u64>) -> Self {
        assert_eq!(
            names.len(),
            lens.len(),
            "RefInfo: names ({}) and lens ({}) must have equal length",
            names.len(),
            lens.len()
        );
        Self {
            ref_names: names,
            ref_lens: lens,
        }
    }

    /// Number of reference sequences.
    #[inline]
    pub fn num_refs(&self) -> usize {
        self.ref_names.len()
    }

    /// Get the name of reference `i`.
    #[inline]
    pub fn name(&self, i: usize) -> &str {
        &self.ref_names[i]
    }

    /// Get the length of reference `i`.
    #[inline]
    pub fn len(&self, i: usize) -> u64 {
        self.ref_lens[i]
    }

    /// Slice of all reference names.
    #[inline]
    pub fn names(&self) -> &[String] {
        &self.ref_names
    }

    /// Slice of all reference lengths.
    #[inline]
    pub fn lens(&self) -> &[u64] {
        &self.ref_lens
    }

    /// Maximum reference length, or 0 if empty.
    pub fn max_ref_len(&self) -> u64 {
        self.ref_lens.iter().copied().max().unwrap_or(0)
    }

    /// Total bases across all references.
    pub fn total_len(&self) -> u64 {
        self.ref_lens.iter().sum()
    }

    // -----------------------------------------------------------------------
    // Serialization
    // -----------------------------------------------------------------------

    /// Serialize to a writer.
    ///
    /// Format:
    /// ```text
    /// [magic: 8 bytes "PRFINF01"]
    /// [num_refs: u64 LE]
    /// for each reference:
    ///   [name_len: u32 LE] [name_bytes: UTF-8]
    /// [ref_lens: num_refs × u64 LE]
    /// ```
    pub fn save<W: Write>(&self, writer: &mut W) -> Result<()> {
        let mut w = BufWriter::new(writer);
        w.write_all(REFINFO_MAGIC)?;

        let n = self.ref_names.len() as u64;
        w.write_all(&n.to_le_bytes())?;

        // Write names: length-prefixed UTF-8 strings
        for name in &self.ref_names {
            let name_bytes = name.as_bytes();
            let name_len = name_bytes.len() as u32;
            w.write_all(&name_len.to_le_bytes())?;
            w.write_all(name_bytes)?;
        }

        // Write lengths as contiguous u64 LE array
        for &l in &self.ref_lens {
            w.write_all(&l.to_le_bytes())?;
        }

        w.flush()?;
        Ok(())
    }

    /// Deserialize from a reader.
    pub fn load<R: Read>(reader: &mut R) -> Result<Self> {
        let mut r = BufReader::new(reader);

        let mut magic = [0u8; 8];
        r.read_exact(&mut magic)
            .context("failed to read refinfo magic")?;
        if magic != *REFINFO_MAGIC {
            bail!(
                "invalid refinfo magic: expected {:?}, got {:?}",
                REFINFO_MAGIC,
                magic
            );
        }

        let num_refs = read_u64_le(&mut r).context("failed to read num_refs")? as usize;

        // Read names
        let mut ref_names = Vec::with_capacity(num_refs);
        for i in 0..num_refs {
            let name_len = read_u32_le(&mut r)
                .with_context(|| format!("failed to read name length for ref {i}"))?
                as usize;
            let mut name_buf = vec![0u8; name_len];
            r.read_exact(&mut name_buf)
                .with_context(|| format!("failed to read name bytes for ref {i}"))?;
            let name = String::from_utf8(name_buf)
                .with_context(|| format!("invalid UTF-8 in name for ref {i}"))?;
            ref_names.push(name);
        }

        // Read lengths
        let mut ref_lens = Vec::with_capacity(num_refs);
        for i in 0..num_refs {
            let l = read_u64_le(&mut r)
                .with_context(|| format!("failed to read length for ref {i}"))?;
            ref_lens.push(l);
        }

        Ok(Self {
            ref_names,
            ref_lens,
        })
    }
}

fn read_u64_le<R: Read>(reader: &mut R) -> std::io::Result<u64> {
    let mut buf = [0u8; 8];
    reader.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

fn read_u32_le<R: Read>(reader: &mut R) -> std::io::Result<u32> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(u32::from_le_bytes(buf))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_construction() {
        let ri = RefInfo::new(
            vec!["tx1".into(), "tx2".into(), "tx3".into()],
            vec![1000, 2000, 3000],
        );
        assert_eq!(ri.num_refs(), 3);
        assert_eq!(ri.name(0), "tx1");
        assert_eq!(ri.name(2), "tx3");
        assert_eq!(ri.len(0), 1000);
        assert_eq!(ri.len(2), 3000);
        assert_eq!(ri.max_ref_len(), 3000);
        assert_eq!(ri.total_len(), 6000);
    }

    #[test]
    fn test_empty() {
        let ri = RefInfo::new(vec![], vec![]);
        assert_eq!(ri.num_refs(), 0);
        assert_eq!(ri.max_ref_len(), 0);
        assert_eq!(ri.total_len(), 0);
    }

    #[test]
    #[should_panic(expected = "names (2) and lens (3) must have equal length")]
    fn test_mismatched_lengths_panics() {
        RefInfo::new(vec!["a".into(), "b".into()], vec![1, 2, 3]);
    }

    #[test]
    fn test_serialization_roundtrip() {
        let ri = RefInfo::new(
            vec![
                "ENST00000456328.2|ENSG00000290825.1|OTTHUMG00000002860.4|DDX11L2-202|DDX11L2|1657|processed_transcript|".into(),
                "short".into(),
                "".into(),
            ],
            vec![1657, 42, 0],
        );

        let mut buf = Vec::new();
        ri.save(&mut buf).unwrap();

        let ri2 = RefInfo::load(&mut &buf[..]).unwrap();
        assert_eq!(ri2.num_refs(), 3);
        assert_eq!(ri2.name(0), ri.name(0));
        assert_eq!(ri2.name(1), "short");
        assert_eq!(ri2.name(2), "");
        assert_eq!(ri2.len(0), 1657);
        assert_eq!(ri2.len(1), 42);
        assert_eq!(ri2.len(2), 0);
    }

    #[test]
    fn test_roundtrip_many_refs() {
        let n = 10_000;
        let names: Vec<String> = (0..n).map(|i| format!("ref_{i}")).collect();
        let lens: Vec<u64> = (0..n).map(|i| (i as u64) * 1000 + 1).collect();

        let ri = RefInfo::new(names.clone(), lens.clone());

        let mut buf = Vec::new();
        ri.save(&mut buf).unwrap();

        let ri2 = RefInfo::load(&mut &buf[..]).unwrap();
        assert_eq!(ri2.num_refs(), n);
        for i in 0..n {
            assert_eq!(ri2.name(i), names[i]);
            assert_eq!(ri2.len(i), lens[i]);
        }
    }

    #[test]
    fn test_invalid_magic() {
        let data = b"BADMAGIC\x00\x00\x00\x00\x00\x00\x00\x00";
        let result = RefInfo::load(&mut &data[..]);
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(err.contains("invalid refinfo magic"));
    }
}
