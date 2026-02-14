//! Contig (unitig) occurrence table — the inverted tiling index.
//!
//! Maps each unitig in the compacted de Bruijn graph to its packed list of
//! reference occurrences. Each occurrence encodes:
//!
//! ```text
//! [  reference_id  |  position  | orientation ]
//!    (high bits)     (mid bits)    (LSB, 1 bit)
//! ```
//!
//! - **orientation** (bit 0): 1 = forward, 0 = reverse complement
//! - **position** (bits 1..ref_len_bits): position on reference where the unitig occurs
//! - **reference_id** (bits ref_len_bits+1..): which reference transcript
//!
//! The posting lists are stored in a compact bit-packed vector (`BitFieldVec`)
//! with boundaries encoded as a monotone Elias-Fano sequence (sux-rs `EfSeq`).
//!
//! This corresponds to the C++ `basic_contig_table` class.

use anyhow::{Context, Result, bail};
use epserde::deser::Deserialize;
use epserde::ser::Serialize;
use mem_dbg::{MemSize, SizeFlags};
use sux::bits::bit_field_vec::BitFieldVec;
use sux::dict::elias_fano::{EliasFanoBuilder, EfSeq};
use sux::traits::IndexedSeq;
use value_traits::slices::SliceByValue;
use value_traits::slices::SliceByValueMut;
use std::io::{Read, Write};

// ---------------------------------------------------------------------------
// Entry encoding / decoding helpers
// ---------------------------------------------------------------------------

/// Parameters derived from the contig table that govern entry encoding.
///
/// These replace the C++ global-state pattern (`PiscemIndexUtils::ref_shift()`,
/// `PiscemIndexUtils::pos_mask()`).
#[derive(Debug, Clone, Copy)]
pub struct EntryEncoding {
    /// Number of bits used to encode reference positions.
    pub ref_len_bits: u64,
    /// Number of bits used to encode reference IDs.
    pub num_ref_bits: u64,
    /// Total bits per entry: `ref_len_bits + num_ref_bits + 1` (orientation).
    pub entry_width: u64,
    /// Right-shift amount to extract the reference ID: `ref_len_bits + 1`.
    pub ref_shift: u64,
    /// Mask to extract position after right-shifting by 1: `(1 << ref_len_bits) - 1`.
    pub pos_mask: u64,
}

impl EntryEncoding {
    /// Create encoding parameters from the bit widths.
    pub fn new(ref_len_bits: u64, num_ref_bits: u64) -> Self {
        Self {
            ref_len_bits,
            num_ref_bits,
            entry_width: ref_len_bits + num_ref_bits + 1,
            ref_shift: ref_len_bits + 1,
            pos_mask: if ref_len_bits == 0 { 0 } else { (1u64 << ref_len_bits) - 1 },
        }
    }

    /// Compute encoding parameters from data characteristics.
    ///
    /// - `max_ref_len`: length of the longest reference sequence
    /// - `num_refs`: total number of reference sequences
    pub fn from_data(max_ref_len: u64, num_refs: u64) -> Self {
        let ref_len_bits = ceil_log2(max_ref_len + 1);
        let num_ref_bits = ceil_log2(num_refs + 1);
        Self::new(ref_len_bits, num_ref_bits)
    }

    /// Encode a single contig table entry.
    #[inline]
    pub fn encode(&self, ref_id: u32, position: u32, is_fw: bool) -> u64 {
        let mut e = ref_id as u64;
        e <<= self.ref_shift;
        e |= (position as u64) << 1;
        e |= is_fw as u64;
        e
    }

    /// Extract the reference (transcript) ID from a packed entry.
    #[inline]
    pub fn transcript_id(&self, entry: u64) -> u32 {
        (entry >> self.ref_shift) as u32
    }

    /// Extract the reference position from a packed entry.
    #[inline]
    pub fn pos(&self, entry: u64) -> u32 {
        ((entry >> 1) & self.pos_mask) as u32
    }

    /// Extract the orientation from a packed entry (true = forward).
    #[inline]
    pub fn orientation(&self, entry: u64) -> bool {
        (entry & 1) != 0
    }
}

/// Minimum number of bits needed to represent values in `[0, n)`.
/// Returns 0 when `n <= 1`.
pub(crate) fn ceil_log2(n: u64) -> u64 {
    if n <= 1 {
        return 0;
    }
    64 - (n - 1).leading_zeros() as u64
}

// ---------------------------------------------------------------------------
// ContigSpan — view into a contig's posting list
// ---------------------------------------------------------------------------

/// A view over the reference occurrences of a single contig.
///
/// Provides iteration over the packed entries. Equivalent to the C++
/// `sshash::util::contig_span`.
#[derive(Clone)]
pub struct ContigSpan<'a> {
    entries: &'a BitFieldVec<usize>,
    start: usize,
    len: usize,
}

impl<'a> ContigSpan<'a> {
    /// Number of reference occurrences for this contig.
    #[inline]
    pub fn len(&self) -> usize {
        self.len
    }

    /// Whether this contig has no reference occurrences.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Get the raw packed entry at position `i` within this span.
    #[inline]
    pub fn get(&self, i: usize) -> u64 {
        debug_assert!(i < self.len);
        self.entries.index_value(self.start + i) as u64
    }

    /// Iterate over the raw packed entries.
    pub fn iter(&self) -> ContigSpanIter<'a> {
        ContigSpanIter {
            entries: self.entries,
            pos: self.start,
            end: self.start + self.len,
        }
    }
}

impl<'a> IntoIterator for &'a ContigSpan<'a> {
    type Item = u64;
    type IntoIter = ContigSpanIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

/// Iterator over packed entries in a contig span.
pub struct ContigSpanIter<'a> {
    entries: &'a BitFieldVec<usize>,
    pos: usize,
    end: usize,
}

impl Iterator for ContigSpanIter<'_> {
    type Item = u64;

    #[inline]
    fn next(&mut self) -> Option<u64> {
        if self.pos < self.end {
            let v = self.entries.index_value(self.pos) as u64;
            self.pos += 1;
            Some(v)
        } else {
            None
        }
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.end - self.pos;
        (remaining, Some(remaining))
    }
}

impl ExactSizeIterator for ContigSpanIter<'_> {}

// ---------------------------------------------------------------------------
// ContigTable
// ---------------------------------------------------------------------------

/// The inverted tiling index: unitig → packed reference occurrences.
///
/// Corresponds to the C++ `basic_contig_table`.
pub struct ContigTable {
    /// Number of bits used to encode reference positions.
    ref_len_bits: u64,
    /// Number of bits used to encode reference IDs.
    num_ref_bits: u64,
    /// Elias-Fano encoded cumulative offsets into `ctg_entries`.
    /// Length = num_contigs + 1. The entries for contig `i` span
    /// `[ctg_offsets[i], ctg_offsets[i+1])`.
    ctg_offsets: EfSeq,
    /// Bit-packed entry vector. Each entry is `ref_len_bits + num_ref_bits + 1` bits.
    ctg_entries: BitFieldVec<usize>,
}

impl ContigTable {
    /// Number of contigs (unitigs) in the table.
    pub fn num_contigs(&self) -> usize {
        let n = self.ctg_offsets.len();
        if n == 0 { 0 } else { n - 1 }
    }

    /// Total number of occurrence entries across all contigs.
    pub fn num_entries(&self) -> usize {
        self.ctg_entries.len()
    }

    /// The `EntryEncoding` parameters for this table.
    pub fn encoding(&self) -> EntryEncoding {
        EntryEncoding::new(self.ref_len_bits, self.num_ref_bits)
    }

    /// Number of bits used to encode reference positions.
    pub fn ref_len_bits(&self) -> u64 {
        self.ref_len_bits
    }

    /// Number of bits used to encode reference IDs.
    pub fn num_ref_bits(&self) -> u64 {
        self.num_ref_bits
    }

    /// Get the posting list (reference occurrences) for a given contig.
    ///
    /// Returns a `ContigSpan` that can be iterated to yield packed entries.
    /// Use `EntryEncoding` to decode each entry.
    #[inline]
    pub fn contig_entries(&self, contig_id: u64) -> ContigSpan<'_> {
        let start = unsafe { self.ctg_offsets.get_unchecked(contig_id as usize) };
        let end = unsafe { self.ctg_offsets.get_unchecked(contig_id as usize + 1) };
        ContigSpan {
            entries: &self.ctg_entries,
            start,
            len: end - start,
        }
    }

    /// Approximate size in bytes of the in-memory representation.
    pub fn size_bytes(&self) -> usize {
        let ef_bytes = self.ctg_offsets.mem_size(SizeFlags::default());
        let entries_bits = self.ctg_entries.len() * self.ctg_entries.bit_width();
        let entries_bytes = entries_bits.div_ceil(8);
        ef_bytes + entries_bytes + std::mem::size_of::<Self>()
    }

    // -----------------------------------------------------------------------
    // Serialization
    // -----------------------------------------------------------------------

    /// Serialize the contig table to a writer.
    ///
    /// Format:
    /// ```text
    /// [magic: 8 bytes "PCTAB02\0"]
    /// [ref_len_bits: u64 LE]
    /// [num_ref_bits: u64 LE]
    /// [ctg_offsets: epserde EfSeq binary]
    /// [ctg_entries: epserde BitFieldVec binary]
    /// ```
    pub fn save<W: Write>(&self, writer: &mut W) -> Result<()> {
        writer.write_all(CONTIG_TABLE_MAGIC)?;
        writer.write_all(&self.ref_len_bits.to_le_bytes())?;
        writer.write_all(&self.num_ref_bits.to_le_bytes())?;
        unsafe {
            self.ctg_offsets
                .serialize(writer)
                .map_err(std::io::Error::other)
                .context("failed to serialize contig offsets")?;
        }
        unsafe {
            self.ctg_entries
                .serialize(writer)
                .map_err(std::io::Error::other)
                .context("failed to serialize contig entries")?;
        }
        Ok(())
    }

    /// Deserialize a contig table from a reader.
    pub fn load<R: Read>(reader: &mut R) -> Result<Self> {
        let mut magic = [0u8; 8];
        reader.read_exact(&mut magic)
            .context("failed to read contig table magic")?;
        if magic != *CONTIG_TABLE_MAGIC {
            bail!("invalid contig table magic: expected {:?}, got {:?}",
                  CONTIG_TABLE_MAGIC, magic);
        }

        let ref_len_bits = read_u64_le(reader).context("failed to read ref_len_bits")?;
        let num_ref_bits = read_u64_le(reader).context("failed to read num_ref_bits")?;

        let ctg_offsets = unsafe {
            EfSeq::deserialize_full(reader)
                .map_err(std::io::Error::other)
                .context("failed to deserialize contig offsets")?
        };

        let ctg_entries = unsafe {
            BitFieldVec::deserialize_full(reader)
                .map_err(std::io::Error::other)
                .context("failed to deserialize contig entries")?
        };

        Ok(Self {
            ref_len_bits,
            num_ref_bits,
            ctg_offsets,
            ctg_entries,
        })
    }
}

const CONTIG_TABLE_MAGIC: &[u8; 8] = b"PCTAB02\0";

fn read_u64_le<R: Read>(reader: &mut R) -> std::io::Result<u64> {
    let mut buf = [0u8; 8];
    reader.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

impl std::fmt::Debug for ContigTable {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ContigTable")
            .field("num_contigs", &self.num_contigs())
            .field("num_entries", &self.num_entries())
            .field("ref_len_bits", &self.ref_len_bits)
            .field("num_ref_bits", &self.num_ref_bits)
            .field("entry_width", &(self.ref_len_bits + self.num_ref_bits + 1))
            .finish()
    }
}

// ---------------------------------------------------------------------------
// Builder
// ---------------------------------------------------------------------------

/// A decoded occurrence record used during construction.
#[derive(Debug, Clone, Copy)]
pub struct OccurrenceRecord {
    /// Internal unitig rank (0-based, determined by segment file order).
    pub unitig_rank: u32,
    /// Reference sequence ID.
    pub ref_id: u32,
    /// Position on the reference where this unitig occurs.
    pub position: u32,
    /// Orientation on the reference (true = forward).
    pub is_fw: bool,
}

/// Builder for constructing a `ContigTable` from occurrence records.
///
/// Usage:
/// 1. Create with `new(num_unitigs, max_ref_len, num_refs)`.
/// 2. Call `add_occurrence()` for each (unitig, ref, pos, orientation) tuple.
/// 3. Call `build()` to produce the final `ContigTable`.
pub struct ContigTableBuilder {
    num_unitigs: usize,
    encoding: EntryEncoding,
    /// Per-unitig occurrence lists, accumulated during construction.
    occurrences: Vec<Vec<u64>>,
}

impl ContigTableBuilder {
    /// Create a new builder.
    ///
    /// - `num_unitigs`: total number of unitigs (determines offset vector size)
    /// - `max_ref_len`: length of the longest reference (determines position bit width)
    /// - `num_refs`: total number of reference sequences (determines ref_id bit width)
    pub fn new(num_unitigs: usize, max_ref_len: u64, num_refs: u64) -> Self {
        let encoding = EntryEncoding::from_data(max_ref_len, num_refs);
        Self {
            num_unitigs,
            encoding,
            occurrences: vec![Vec::new(); num_unitigs],
        }
    }

    /// Add a single occurrence record.
    pub fn add_occurrence(&mut self, unitig_rank: u32, ref_id: u32, position: u32, is_fw: bool) {
        let encoded = self.encoding.encode(ref_id, position, is_fw);
        self.occurrences[unitig_rank as usize].push(encoded);
    }

    /// Build the final `ContigTable`.
    pub fn build(self) -> ContigTable {
        // Build cumulative offset vector
        let mut offsets = Vec::with_capacity(self.num_unitigs + 1);
        offsets.push(0u64);
        let mut total: u64 = 0;
        for occ_list in &self.occurrences {
            total += occ_list.len() as u64;
            offsets.push(total);
        }

        // Build Elias-Fano from offsets
        let n = offsets.len();
        let u = if n > 0 { offsets[n - 1] as usize + 1 } else { 1 };
        let mut ef_builder = EliasFanoBuilder::new(n, u);
        for &v in &offsets {
            ef_builder.push(v as usize);
        }
        let ctg_offsets = ef_builder.build_with_seq();

        // Build BitFieldVec with packed entries
        let entry_width = self.encoding.entry_width as usize;
        let mut ctg_entries = BitFieldVec::new(entry_width, total as usize);
        let mut idx = 0usize;
        for occ_list in &self.occurrences {
            for &encoded in occ_list {
                ctg_entries.set_value(idx, encoded as usize);
                idx += 1;
            }
        }

        ContigTable {
            ref_len_bits: self.encoding.ref_len_bits,
            num_ref_bits: self.encoding.num_ref_bits,
            ctg_offsets,
            ctg_entries,
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ceil_log2() {
        assert_eq!(ceil_log2(0), 0);
        assert_eq!(ceil_log2(1), 0);
        assert_eq!(ceil_log2(2), 1);
        assert_eq!(ceil_log2(3), 2);
        assert_eq!(ceil_log2(4), 2);
        assert_eq!(ceil_log2(5), 3);
        assert_eq!(ceil_log2(8), 3);
        assert_eq!(ceil_log2(9), 4);
        assert_eq!(ceil_log2(256), 8);
        assert_eq!(ceil_log2(257), 9);
        assert_eq!(ceil_log2(1 << 20), 20);
    }

    #[test]
    fn test_entry_encoding_roundtrip() {
        let enc = EntryEncoding::from_data(1_000_000, 100_000);

        // Verify bit widths
        assert_eq!(enc.ref_len_bits, 20); // ceil_log2(1_000_001) = 20
        assert_eq!(enc.num_ref_bits, 17); // ceil_log2(100_001) = 17
        assert_eq!(enc.entry_width, 38);  // 20 + 17 + 1

        // Roundtrip test
        let test_cases = [
            (0u32, 0u32, true),
            (0, 0, false),
            (99_999, 999_999, true),
            (50_000, 500_000, false),
            (1, 1, true),
        ];

        for (ref_id, position, is_fw) in test_cases {
            let encoded = enc.encode(ref_id, position, is_fw);
            assert_eq!(enc.transcript_id(encoded), ref_id,
                       "ref_id mismatch for ({ref_id}, {position}, {is_fw})");
            assert_eq!(enc.pos(encoded), position,
                       "position mismatch for ({ref_id}, {position}, {is_fw})");
            assert_eq!(enc.orientation(encoded), is_fw,
                       "orientation mismatch for ({ref_id}, {position}, {is_fw})");
        }
    }

    #[test]
    fn test_build_and_query() {
        // 3 unitigs, 2 references, max ref len 1000
        let mut builder = ContigTableBuilder::new(3, 1000, 2);

        // Unitig 0: occurs at ref 0 pos 100 FW, and ref 1 pos 200 RC
        builder.add_occurrence(0, 0, 100, true);
        builder.add_occurrence(0, 1, 200, false);

        // Unitig 1: occurs at ref 0 pos 500 FW
        builder.add_occurrence(1, 0, 500, true);

        // Unitig 2: occurs at ref 0 pos 0 FW, ref 0 pos 300 RC, ref 1 pos 0 FW
        builder.add_occurrence(2, 0, 0, true);
        builder.add_occurrence(2, 0, 300, false);
        builder.add_occurrence(2, 1, 0, true);

        let table = builder.build();
        let enc = table.encoding();

        assert_eq!(table.num_contigs(), 3);
        assert_eq!(table.num_entries(), 6);

        // Query unitig 0
        let span0 = table.contig_entries(0);
        assert_eq!(span0.len(), 2);
        let e0 = span0.get(0);
        assert_eq!(enc.transcript_id(e0), 0);
        assert_eq!(enc.pos(e0), 100);
        assert!(enc.orientation(e0));
        let e1 = span0.get(1);
        assert_eq!(enc.transcript_id(e1), 1);
        assert_eq!(enc.pos(e1), 200);
        assert!(!enc.orientation(e1));

        // Query unitig 1
        let span1 = table.contig_entries(1);
        assert_eq!(span1.len(), 1);

        // Query unitig 2
        let span2 = table.contig_entries(2);
        assert_eq!(span2.len(), 3);

        // Test iteration
        let entries: Vec<u64> = span2.iter().collect();
        assert_eq!(entries.len(), 3);
        assert_eq!(enc.transcript_id(entries[0]), 0);
        assert_eq!(enc.pos(entries[0]), 0);
        assert!(enc.orientation(entries[0]));
    }

    #[test]
    fn test_empty_contigs() {
        // Some unitigs may have no occurrences
        let mut builder = ContigTableBuilder::new(3, 100, 1);
        // Only unitig 1 has an occurrence
        builder.add_occurrence(1, 0, 50, true);
        let table = builder.build();

        assert_eq!(table.num_contigs(), 3);
        assert_eq!(table.num_entries(), 1);

        assert!(table.contig_entries(0).is_empty());
        assert_eq!(table.contig_entries(1).len(), 1);
        assert!(table.contig_entries(2).is_empty());
    }

    #[test]
    fn test_serialization_roundtrip() {
        let mut builder = ContigTableBuilder::new(3, 10_000, 50);
        builder.add_occurrence(0, 0, 100, true);
        builder.add_occurrence(0, 1, 200, false);
        builder.add_occurrence(1, 0, 500, true);
        builder.add_occurrence(2, 0, 0, true);
        builder.add_occurrence(2, 49, 9999, false);
        let table = builder.build();

        // Serialize
        let mut buf = Vec::new();
        table.save(&mut buf).unwrap();

        // Deserialize
        let table2 = ContigTable::load(&mut &buf[..]).unwrap();

        // Verify
        assert_eq!(table2.num_contigs(), table.num_contigs());
        assert_eq!(table2.num_entries(), table.num_entries());
        assert_eq!(table2.ref_len_bits(), table.ref_len_bits());
        assert_eq!(table2.num_ref_bits(), table.num_ref_bits());

        let enc = table2.encoding();
        for contig_id in 0..3u64 {
            let span_orig = table.contig_entries(contig_id);
            let span_loaded = table2.contig_entries(contig_id);
            assert_eq!(span_orig.len(), span_loaded.len());
            for i in 0..span_orig.len() {
                let a = span_orig.get(i);
                let b = span_loaded.get(i);
                assert_eq!(enc.transcript_id(a), enc.transcript_id(b));
                assert_eq!(enc.pos(a), enc.pos(b));
                assert_eq!(enc.orientation(a), enc.orientation(b));
            }
        }
    }

    #[test]
    fn test_large_values() {
        // Test with large reference lengths and many references
        let max_ref_len = 3_000_000_000u64; // 3 billion (human genome scale)
        let num_refs = 200_000u64;

        let mut builder = ContigTableBuilder::new(2, max_ref_len, num_refs);
        builder.add_occurrence(0, 199_999, 2_999_999_999, true);
        builder.add_occurrence(1, 0, 0, false);
        let table = builder.build();
        let enc = table.encoding();

        let span = table.contig_entries(0);
        let e = span.get(0);
        assert_eq!(enc.transcript_id(e), 199_999);
        assert_eq!(enc.pos(e), 2_999_999_999);
        assert!(enc.orientation(e));
    }
}
