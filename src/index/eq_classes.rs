//! Equivalence class map — maps tiles (contigs) to equivalence classes of
//! reference targets.
//!
//! Each tile (unitig) maps to an equivalence class (EC), which is a set of
//! `(transcript_id, orientation)` pairs. Tiles that map to exactly the same set
//! of references with the same orientations share an EC, saving space.
//!
//! Entry packing format:
//! ```text
//! [  transcript_id (high bits)  |  orientation (2 LSBs)  ]
//! ```
//!
//! Orientation values:
//! - `0` (FW): forward only
//! - `1` (RC): reverse complement only
//! - `2` (BOTH): both orientations
//!
//! This corresponds to the C++ `equivalence_class_map` class.

use anyhow::{Context, Result, bail};
use cseq::elias_fano::Sequence as CseqSequence;
use dyn_size_of::GetSize;
use epserde::deser::Deserialize;
use epserde::ser::Serialize;
use sux::bits::bit_field_vec::BitFieldVec;
use value_traits::slices::{SliceByValue, SliceByValueMut};
use std::collections::HashMap;
use std::io::{Read, Write};

use super::contig_table::ceil_log2;

// ---------------------------------------------------------------------------
// Orientation
// ---------------------------------------------------------------------------

/// Orientation of a transcript in an equivalence class entry.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum Orientation {
    /// Forward only.
    Forward = 0,
    /// Reverse complement only.
    ReverseComplement = 1,
    /// Both orientations.
    Both = 2,
}

impl Orientation {
    /// Decode from the 2-bit representation.
    #[inline]
    pub fn from_bits(bits: u64) -> Self {
        match bits & 0x3 {
            0 => Self::Forward,
            1 => Self::ReverseComplement,
            2 => Self::Both,
            _ => Self::Both, // 3 is unused; treat conservatively as Both
        }
    }

    /// Encode to 2-bit representation.
    #[inline]
    pub fn to_bits(self) -> u64 {
        self as u64
    }
}

// ---------------------------------------------------------------------------
// EcSpan — view into an equivalence class's label entries
// ---------------------------------------------------------------------------

/// A view over the label entries of a single equivalence class.
///
/// Each entry encodes a `(transcript_id, orientation)` pair packed as:
/// `(transcript_id << 2) | orientation_bits`.
///
/// Corresponds to the C++ `sshash::util::ec_span`.
#[derive(Clone)]
pub struct EcSpan<'a> {
    entries: &'a BitFieldVec<usize>,
    start: usize,
    len: usize,
}

impl<'a> EcSpan<'a> {
    /// Number of label entries in this equivalence class.
    #[inline]
    pub fn len(&self) -> usize {
        self.len
    }

    /// Whether this equivalence class has no entries.
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
    pub fn iter(&self) -> EcSpanIter<'a> {
        EcSpanIter {
            entries: self.entries,
            pos: self.start,
            end: self.start + self.len,
        }
    }
}

/// Extract the transcript ID from a raw packed EC entry.
#[inline]
pub fn ec_entry_transcript_id(entry: u64) -> u32 {
    (entry >> 2) as u32
}

/// Extract the orientation from a raw packed EC entry.
#[inline]
pub fn ec_entry_orientation(entry: u64) -> Orientation {
    Orientation::from_bits(entry)
}

impl<'a> IntoIterator for &'a EcSpan<'a> {
    type Item = u64;
    type IntoIter = EcSpanIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

/// Iterator over packed entries in an EC span.
pub struct EcSpanIter<'a> {
    entries: &'a BitFieldVec<usize>,
    pos: usize,
    end: usize,
}

impl Iterator for EcSpanIter<'_> {
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

impl ExactSizeIterator for EcSpanIter<'_> {}

// ---------------------------------------------------------------------------
// EqClassMap
// ---------------------------------------------------------------------------

/// Equivalence class map: tile → EC → label entries.
///
/// Maps each tile (contig/unitig) to an equivalence class, and each
/// equivalence class to its set of `(transcript_id, orientation)` entries.
///
/// Corresponds to the C++ `equivalence_class_map`.
pub struct EqClassMap {
    /// Maps tile (contig) ID → equivalence class ID.
    tile_ec_ids: BitFieldVec<usize>,
    /// Elias-Fano encoded cumulative offsets into `label_entries`.
    /// Length = num_ecs + 1.
    label_list_offsets: CseqSequence,
    /// Packed label entries: `(transcript_id << 2) | orientation_bits`.
    label_entries: BitFieldVec<usize>,
}

impl EqClassMap {
    /// Number of tiles (contigs) in the map.
    pub fn num_tiles(&self) -> usize {
        self.tile_ec_ids.len()
    }

    /// Number of distinct equivalence classes.
    pub fn num_ecs(&self) -> usize {
        let n = self.label_list_offsets.len();
        if n == 0 { 0 } else { n - 1 }
    }

    /// Total number of label entries across all ECs.
    pub fn num_label_entries(&self) -> usize {
        self.label_entries.len()
    }

    /// Get the equivalence class ID for a given tile.
    #[inline]
    pub fn ec_for_tile(&self, tile_id: u64) -> u64 {
        self.tile_ec_ids.index_value(tile_id as usize) as u64
    }

    /// Get the label entries for a given equivalence class.
    #[inline]
    pub fn entries_for_ec(&self, ec_id: u64) -> EcSpan<'_> {
        let start =
            unsafe { self.label_list_offsets.get_unchecked(ec_id as usize) } as usize;
        let end =
            unsafe { self.label_list_offsets.get_unchecked(ec_id as usize + 1) } as usize;
        EcSpan {
            entries: &self.label_entries,
            start,
            len: end - start,
        }
    }

    /// Get the label entries for the equivalence class of a given tile.
    ///
    /// Convenience method combining `ec_for_tile` and `entries_for_ec`.
    #[inline]
    pub fn entries_for_tile(&self, tile_id: u64) -> EcSpan<'_> {
        self.entries_for_ec(self.ec_for_tile(tile_id))
    }

    /// Approximate size in bytes of the in-memory representation.
    pub fn size_bytes(&self) -> usize {
        let tile_bits = self.tile_ec_ids.len() * self.tile_ec_ids.bit_width();
        let tile_bytes = (tile_bits + 7) / 8;
        let ef_bytes = self.label_list_offsets.size_bytes();
        let label_bits = self.label_entries.len() * self.label_entries.bit_width();
        let label_bytes = (label_bits + 7) / 8;
        tile_bytes + ef_bytes + label_bytes + std::mem::size_of::<Self>()
    }

    // -----------------------------------------------------------------------
    // Serialization
    // -----------------------------------------------------------------------

    /// Serialize the equivalence class map to a writer.
    ///
    /// Format:
    /// ```text
    /// [magic: 8 bytes "PECTB01\0"]
    /// [tile_ec_ids: epserde BitFieldVec binary]
    /// [label_list_offsets: cseq EF binary]
    /// [label_entries: epserde BitFieldVec binary]
    /// ```
    pub fn save<W: Write>(&self, writer: &mut W) -> Result<()> {
        writer.write_all(EC_TABLE_MAGIC)?;
        unsafe {
            self.tile_ec_ids
                .serialize(writer)
                .map_err(std::io::Error::other)
                .context("failed to serialize tile_ec_ids")?;
        }
        self.label_list_offsets
            .write(writer)
            .context("failed to serialize label_list_offsets")?;
        unsafe {
            self.label_entries
                .serialize(writer)
                .map_err(std::io::Error::other)
                .context("failed to serialize label_entries")?;
        }
        Ok(())
    }

    /// Deserialize an equivalence class map from a reader.
    pub fn load<R: Read>(reader: &mut R) -> Result<Self> {
        let mut magic = [0u8; 8];
        reader
            .read_exact(&mut magic)
            .context("failed to read EC table magic")?;
        if magic != *EC_TABLE_MAGIC {
            bail!(
                "invalid EC table magic: expected {:?}, got {:?}",
                EC_TABLE_MAGIC,
                magic
            );
        }

        let tile_ec_ids = unsafe {
            BitFieldVec::deserialize_full(reader)
                .map_err(std::io::Error::other)
                .context("failed to deserialize tile_ec_ids")?
        };

        let label_list_offsets =
            CseqSequence::read(reader).context("failed to deserialize label_list_offsets")?;

        let label_entries = unsafe {
            BitFieldVec::deserialize_full(reader)
                .map_err(std::io::Error::other)
                .context("failed to deserialize label_entries")?
        };

        Ok(Self {
            tile_ec_ids,
            label_list_offsets,
            label_entries,
        })
    }
}

const EC_TABLE_MAGIC: &[u8; 8] = b"PECTB01\0";

impl std::fmt::Debug for EqClassMap {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("EqClassMap")
            .field("num_tiles", &self.num_tiles())
            .field("num_ecs", &self.num_ecs())
            .field("num_label_entries", &self.num_label_entries())
            .finish()
    }
}

// ---------------------------------------------------------------------------
// Builder
// ---------------------------------------------------------------------------

/// Builder for constructing an `EqClassMap`.
///
/// Usage:
/// 1. Create with `new(num_tiles)`.
/// 2. Call `add_tile(tile_id, labels)` for each tile.
/// 3. Call `build()` to produce the final `EqClassMap`.
pub struct EqClassMapBuilder {
    num_tiles: usize,
    /// EC ID assigned to each tile (indexed by tile_id).
    tile_to_ec: Vec<u32>,
    /// Unique labels in order of first occurrence.
    labels: Vec<Vec<(u32, Orientation)>>,
    /// Deduplication map: sorted label → EC ID.
    label_to_ec: HashMap<Vec<(u32, Orientation)>, u32>,
    /// Largest transcript ID seen.
    max_transcript_id: u32,
}

impl EqClassMapBuilder {
    /// Create a new builder for `num_tiles` tiles.
    pub fn new(num_tiles: usize) -> Self {
        Self {
            num_tiles,
            tile_to_ec: vec![0; num_tiles],
            labels: Vec::new(),
            label_to_ec: HashMap::new(),
            max_transcript_id: 0,
        }
    }

    /// Add a tile with its label (set of `(transcript_id, orientation)` pairs).
    ///
    /// The labels will be sorted internally for consistent EC assignment.
    pub fn add_tile(&mut self, tile_id: usize, mut label: Vec<(u32, Orientation)>) {
        // Sort for consistent hashing/comparison
        label.sort();

        // Track max transcript ID
        for &(tid, _) in &label {
            self.max_transcript_id = self.max_transcript_id.max(tid);
        }

        // Look up or create EC
        let next_ec = self.labels.len() as u32;
        let ec_id = *self.label_to_ec.entry(label.clone()).or_insert_with(|| {
            self.labels.push(label);
            next_ec
        });
        self.tile_to_ec[tile_id] = ec_id;
    }

    /// Build the final `EqClassMap`.
    pub fn build(self) -> EqClassMap {
        let num_ecs = self.labels.len();

        // Bit width for tile_ec_ids: need to represent EC IDs [0, num_ecs)
        let ec_bits = ceil_log2(num_ecs as u64 + 1).max(1) as usize;
        let mut tile_ec_ids = BitFieldVec::new(ec_bits, self.num_tiles);
        for (i, &ec_id) in self.tile_to_ec.iter().enumerate() {
            tile_ec_ids.set_value(i, ec_id as usize);
        }

        // Entry bit width: ceil(log2(max_tid + 1)) + 2 for orientation bits
        let entry_bits = (ceil_log2(self.max_transcript_id as u64 + 1) + 2).max(2) as usize;
        let total_entries: usize = self.labels.iter().map(|l| l.len()).sum();

        // Build cumulative offsets and packed entries
        let mut offsets = Vec::with_capacity(num_ecs + 1);
        offsets.push(0u64);
        let mut cumulative: u64 = 0;

        let mut label_entries = BitFieldVec::new(entry_bits, total_entries);
        let mut idx = 0usize;

        for label in &self.labels {
            for &(tid, ori) in label {
                let packed = ((tid as u64) << 2) | ori.to_bits();
                label_entries.set_value(idx, packed as usize);
                idx += 1;
            }
            cumulative += label.len() as u64;
            offsets.push(cumulative);
        }

        let label_list_offsets = CseqSequence::with_items_from_slice(&offsets);

        EqClassMap {
            tile_ec_ids,
            label_list_offsets,
            label_entries,
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
    fn test_orientation_roundtrip() {
        for ori in [Orientation::Forward, Orientation::ReverseComplement, Orientation::Both] {
            assert_eq!(Orientation::from_bits(ori.to_bits()), ori);
        }
    }

    #[test]
    fn test_entry_packing() {
        // transcript_id=42, orientation=FW → packed = (42 << 2) | 0 = 168
        let packed = (42u64 << 2) | Orientation::Forward.to_bits();
        assert_eq!(ec_entry_transcript_id(packed), 42);
        assert_eq!(ec_entry_orientation(packed), Orientation::Forward);

        // transcript_id=99, orientation=RC → packed = (99 << 2) | 1 = 397
        let packed = (99u64 << 2) | Orientation::ReverseComplement.to_bits();
        assert_eq!(ec_entry_transcript_id(packed), 99);
        assert_eq!(ec_entry_orientation(packed), Orientation::ReverseComplement);

        // transcript_id=0, orientation=BOTH → packed = (0 << 2) | 2 = 2
        let packed = (0u64 << 2) | Orientation::Both.to_bits();
        assert_eq!(ec_entry_transcript_id(packed), 0);
        assert_eq!(ec_entry_orientation(packed), Orientation::Both);
    }

    #[test]
    fn test_build_and_query() {
        // 5 tiles, 3 references
        let mut builder = EqClassMapBuilder::new(5);

        // Tiles 0, 2 share the same label: {(0, FW), (1, RC)}
        builder.add_tile(0, vec![(0, Orientation::Forward), (1, Orientation::ReverseComplement)]);
        builder.add_tile(2, vec![(1, Orientation::ReverseComplement), (0, Orientation::Forward)]); // same after sort

        // Tile 1: {(2, Both)}
        builder.add_tile(1, vec![(2, Orientation::Both)]);

        // Tiles 3, 4 share: {(0, FW)}
        builder.add_tile(3, vec![(0, Orientation::Forward)]);
        builder.add_tile(4, vec![(0, Orientation::Forward)]);

        let ec_map = builder.build();

        // Should have 3 distinct ECs
        assert_eq!(ec_map.num_tiles(), 5);
        assert_eq!(ec_map.num_ecs(), 3);
        assert_eq!(ec_map.num_label_entries(), 4); // 2 + 1 + 1

        // Tiles 0 and 2 should share the same EC
        assert_eq!(ec_map.ec_for_tile(0), ec_map.ec_for_tile(2));

        // Tiles 3 and 4 should share the same EC
        assert_eq!(ec_map.ec_for_tile(3), ec_map.ec_for_tile(4));

        // Tile 1 should have a different EC
        assert_ne!(ec_map.ec_for_tile(0), ec_map.ec_for_tile(1));
        assert_ne!(ec_map.ec_for_tile(1), ec_map.ec_for_tile(3));

        // Check entries for tile 0's EC
        let span0 = ec_map.entries_for_tile(0);
        assert_eq!(span0.len(), 2);
        let entries: Vec<u64> = span0.iter().collect();
        // After sorting: (0, FW), (1, RC)
        assert_eq!(ec_entry_transcript_id(entries[0]), 0);
        assert_eq!(ec_entry_orientation(entries[0]), Orientation::Forward);
        assert_eq!(ec_entry_transcript_id(entries[1]), 1);
        assert_eq!(ec_entry_orientation(entries[1]), Orientation::ReverseComplement);

        // Check entries for tile 1's EC
        let span1 = ec_map.entries_for_tile(1);
        assert_eq!(span1.len(), 1);
        let e = span1.get(0);
        assert_eq!(ec_entry_transcript_id(e), 2);
        assert_eq!(ec_entry_orientation(e), Orientation::Both);

        // Check entries for tile 3's EC
        let span3 = ec_map.entries_for_tile(3);
        assert_eq!(span3.len(), 1);
        let e = span3.get(0);
        assert_eq!(ec_entry_transcript_id(e), 0);
        assert_eq!(ec_entry_orientation(e), Orientation::Forward);
    }

    #[test]
    fn test_single_tile_single_ec() {
        let mut builder = EqClassMapBuilder::new(1);
        builder.add_tile(0, vec![(5, Orientation::Both)]);
        let ec_map = builder.build();

        assert_eq!(ec_map.num_tiles(), 1);
        assert_eq!(ec_map.num_ecs(), 1);
        assert_eq!(ec_map.ec_for_tile(0), 0);

        let span = ec_map.entries_for_ec(0);
        assert_eq!(span.len(), 1);
        assert_eq!(ec_entry_transcript_id(span.get(0)), 5);
        assert_eq!(ec_entry_orientation(span.get(0)), Orientation::Both);
    }

    #[test]
    fn test_many_tiles_many_ecs() {
        // Each tile gets its own unique EC
        let n = 100;
        let mut builder = EqClassMapBuilder::new(n);
        for i in 0..n {
            builder.add_tile(i, vec![(i as u32, Orientation::Forward)]);
        }
        let ec_map = builder.build();

        assert_eq!(ec_map.num_tiles(), n);
        assert_eq!(ec_map.num_ecs(), n);

        for i in 0..n {
            let span = ec_map.entries_for_tile(i as u64);
            assert_eq!(span.len(), 1);
            assert_eq!(ec_entry_transcript_id(span.get(0)), i as u32);
        }
    }

    #[test]
    fn test_serialization_roundtrip() {
        let mut builder = EqClassMapBuilder::new(5);
        builder.add_tile(0, vec![(0, Orientation::Forward), (1, Orientation::ReverseComplement)]);
        builder.add_tile(1, vec![(2, Orientation::Both)]);
        builder.add_tile(2, vec![(0, Orientation::Forward), (1, Orientation::ReverseComplement)]);
        builder.add_tile(3, vec![(0, Orientation::Forward)]);
        builder.add_tile(4, vec![(0, Orientation::Forward)]);
        let ec_map = builder.build();

        // Serialize
        let mut buf = Vec::new();
        ec_map.save(&mut buf).unwrap();

        // Deserialize
        let ec_map2 = EqClassMap::load(&mut &buf[..]).unwrap();

        // Verify structure
        assert_eq!(ec_map2.num_tiles(), ec_map.num_tiles());
        assert_eq!(ec_map2.num_ecs(), ec_map.num_ecs());
        assert_eq!(ec_map2.num_label_entries(), ec_map.num_label_entries());

        // Verify all tile → EC mappings
        for tile_id in 0..5u64 {
            assert_eq!(
                ec_map2.ec_for_tile(tile_id),
                ec_map.ec_for_tile(tile_id),
                "tile {tile_id} EC mismatch"
            );
        }

        // Verify all EC entries
        for ec_id in 0..ec_map.num_ecs() as u64 {
            let span1 = ec_map.entries_for_ec(ec_id);
            let span2 = ec_map2.entries_for_ec(ec_id);
            assert_eq!(span1.len(), span2.len(), "EC {ec_id} length mismatch");
            for i in 0..span1.len() {
                assert_eq!(span1.get(i), span2.get(i), "EC {ec_id} entry {i} mismatch");
            }
        }
    }

    #[test]
    fn test_invalid_magic() {
        let data = b"BADMAGIC";
        let result = EqClassMap::load(&mut &data[..]);
        assert!(result.is_err());
    }

    #[test]
    fn test_large_transcript_ids() {
        let mut builder = EqClassMapBuilder::new(2);
        builder.add_tile(0, vec![(200_000, Orientation::Forward), (199_999, Orientation::ReverseComplement)]);
        builder.add_tile(1, vec![(0, Orientation::Both)]);
        let ec_map = builder.build();

        let span = ec_map.entries_for_tile(0);
        assert_eq!(span.len(), 2);
        // After sorting by (tid, ori): (199_999, RC), (200_000, FW)
        let entries: Vec<u64> = span.iter().collect();
        assert_eq!(ec_entry_transcript_id(entries[0]), 199_999);
        assert_eq!(ec_entry_orientation(entries[0]), Orientation::ReverseComplement);
        assert_eq!(ec_entry_transcript_id(entries[1]), 200_000);
        assert_eq!(ec_entry_orientation(entries[1]), Orientation::Forward);
    }
}
