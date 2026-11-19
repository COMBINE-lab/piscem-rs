//! Reference sequence metadata — names and lengths.
//!
//! Stores the name and length of every reference sequence (transcript) in the
//! index. This corresponds to the C++ `m_ref_names` / `m_ref_lens` vectors
//! stored in the `.refinfo` file.

use anyhow::{Context, Result, bail};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

/// Magic bytes for the piscem-rs refinfo file format.
const REFINFO_MAGIC: &[u8; 8] = b"PRFINF01";

/// Reference sequence metadata.
///
/// Contains the name and length of each reference sequence in the index,
/// stored in reference-ID order (i.e. `ref_names[i]` and `ref_lens[i]`
/// describe reference `i`).
#[derive(Debug, Clone)]
pub struct RefInfo {
    /// Every reference name, back to back: one allocation instead of one
    /// `String` (24 bytes plus a heap block) per reference.
    names: String,
    /// End of each name in `names`; name `i` starts where name `i - 1` ends.
    name_ends: Vec<usize>,
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
        let mut info = Self {
            names: String::with_capacity(names.iter().map(String::len).sum()),
            name_ends: Vec::with_capacity(names.len()),
            ref_lens: lens,
        };
        for name in &names {
            info.push_name(name);
        }
        info
    }

    fn push_name(&mut self, name: &str) {
        self.names.push_str(name);
        self.name_ends.push(self.names.len());
    }

    /// Number of reference sequences.
    #[inline]
    pub fn num_refs(&self) -> usize {
        self.name_ends.len()
    }

    /// Get the name of reference `i`.
    #[inline]
    pub fn name(&self, i: usize) -> &str {
        let start = if i == 0 { 0 } else { self.name_ends[i - 1] };
        &self.names[start..self.name_ends[i]]
    }

    /// Get the length of reference `i`.
    #[inline]
    pub fn len(&self, i: usize) -> u64 {
        self.ref_lens[i]
    }

    /// All reference names, in reference-ID order.
    pub fn names(&self) -> impl ExactSizeIterator<Item = &str> + '_ {
        (0..self.num_refs()).map(|i| self.name(i))
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

        let n = self.num_refs() as u64;
        w.write_all(&n.to_le_bytes())?;

        // Write names: length-prefixed UTF-8 strings
        for name in self.names() {
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
        Self::load_with_file_size(reader, None)
    }

    /// Deserialize from a file, reserving space for the names before reading them.
    ///
    /// The file size, less the header, name-length prefixes, and reference
    /// lengths, gives the total name bytes without a preliminary scan.
    pub fn load_file(path: impl AsRef<Path>) -> Result<Self> {
        let mut file = File::open(path.as_ref())
            .with_context(|| format!("failed to open {}", path.as_ref().display()))?;
        let metadata = file
            .metadata()
            .context("failed to read refinfo file metadata")?;
        let file_size = metadata.is_file().then_some(metadata.len());
        Self::load_with_file_size(&mut file, file_size)
    }

    fn load_with_file_size<R: Read>(reader: &mut R, file_size: Option<u64>) -> Result<Self> {
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

        let num_refs = read_u64_le(&mut r).context("failed to read num_refs")?;
        let mut remaining_name_bytes = if let Some(file_size) = file_size {
            let fixed_bytes = num_refs
                .checked_mul(12)
                .and_then(|n| n.checked_add(16))
                .context("refinfo fixed-size fields overflow")?;
            let name_bytes = file_size
                .checked_sub(fixed_bytes)
                .context("refinfo file is too short for its reference count")?;
            Some(usize::try_from(name_bytes).context("refinfo names exceed addressable memory")?)
        } else {
            None
        };
        let num_refs =
            usize::try_from(num_refs).context("refinfo reference count exceeds usize")?;
        let mut names = String::new();
        if let Some(name_bytes) = remaining_name_bytes {
            names
                .try_reserve_exact(name_bytes)
                .context("failed to reserve refinfo name storage")?;
        }

        // Read names
        let mut info = Self {
            names,
            name_ends: Vec::with_capacity(num_refs),
            ref_lens: Vec::with_capacity(num_refs),
        };
        let mut name_buf = Vec::new();
        for i in 0..num_refs {
            let name_len = read_u32_le(&mut r)
                .with_context(|| format!("failed to read name length for ref {i}"))?
                as usize;
            if let Some(remaining) = &mut remaining_name_bytes {
                *remaining = remaining.checked_sub(name_len).with_context(|| {
                    format!("name length for ref {i} exceeds refinfo file size")
                })?;
            }
            name_buf.resize(name_len, 0);
            r.read_exact(&mut name_buf)
                .with_context(|| format!("failed to read name bytes for ref {i}"))?;
            let name = std::str::from_utf8(&name_buf)
                .with_context(|| format!("invalid UTF-8 in name for ref {i}"))?;
            info.push_name(name);
        }

        // Read lengths
        for i in 0..num_refs {
            let l = read_u64_le(&mut r)
                .with_context(|| format!("failed to read length for ref {i}"))?;
            info.ref_lens.push(l);
        }
        info.names.shrink_to_fit();

        Ok(info)
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

    #[test]
    fn test_names_empty_and_multibyte() {
        let names = ["", "é-tx", "tx3", ""];
        let ri = RefInfo::new(
            names.iter().map(|s| s.to_string()).collect(),
            vec![1, 2, 3, 4],
        );
        assert!(ri.names().eq(names));
        let mut buf = Vec::new();
        ri.save(&mut buf).unwrap();
        let back = RefInfo::load(&mut &buf[..]).unwrap();
        assert!(back.names().eq(names));
        assert_eq!(back.names().len(), 4);
        assert_eq!(back.name(1), "é-tx");
        assert_eq!(back.lens(), &[1, 2, 3, 4]);
    }

    #[test]
    fn test_file_load_matches_generic_reader() {
        for names in [
            vec![],
            vec!["".to_owned(), "é-tx".to_owned(), "long".repeat(500)],
        ] {
            let original = RefInfo::new(names.clone(), (0..names.len() as u64).collect());
            let mut bytes = Vec::new();
            original.save(&mut bytes).unwrap();
            // The generic reader accepts trailing bytes; keep that behavior.
            for trailer in [&b""[..], &b"trailing bytes"[..]] {
                let mut file = tempfile::NamedTempFile::new().unwrap();
                file.write_all(&bytes).unwrap();
                file.write_all(trailer).unwrap();
                let loaded = RefInfo::load_file(file.path()).unwrap();
                let generic = RefInfo::load(&mut &bytes[..]).unwrap();
                assert!(loaded.names().eq(generic.names()));
                assert_eq!(loaded.lens(), generic.lens());
                let mut saved = Vec::new();
                loaded.save(&mut saved).unwrap();
                assert_eq!(saved, bytes);
            }
        }
    }

    #[test]
    fn test_file_load_rejects_invalid_sizes_and_names() {
        let mut valid = Vec::new();
        RefInfo::new(vec!["x".to_owned()], vec![42])
            .save(&mut valid)
            .unwrap();
        let mut overflow = valid.clone();
        overflow[8..16].copy_from_slice(&u64::MAX.to_le_bytes());
        let mut bad_length = valid.clone();
        bad_length[16..20].copy_from_slice(&u32::MAX.to_le_bytes());
        let mut invalid_utf8 = valid.clone();
        invalid_utf8[20] = 0xff;
        for bytes in [
            valid[..16].to_vec(),
            valid[..valid.len() - 1].to_vec(),
            overflow,
            bad_length,
            invalid_utf8,
        ] {
            let mut file = tempfile::NamedTempFile::new().unwrap();
            file.write_all(&bytes).unwrap();
            assert!(RefInfo::load_file(file.path()).is_err());
        }
    }
}
