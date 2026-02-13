//! RAD format binary writer.
//!
//! Provides a buffer-based writer for RAD (Reduced Alignment Data) format
//! output, used by salmon/alevin-fry. Includes functions for writing headers
//! and records in both single-cell and bulk modes.
//!
//! Port of C++ `rad_writer` + header/record functions from
//! `piscem-cpp/include/rad/util.hpp`.

use std::io::Write;

use anyhow::Result;

use crate::mapping::hits::{MappingType, SimpleHit};

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// RAD tag type IDs (matching C++ `rad::TagType`).
const TAG_U8: u8 = 1;
const TAG_U16: u8 = 2;
const TAG_U32: u8 = 3;
const TAG_U64: u8 = 4;
const TAG_ARRAY: u8 = 7;
const TAG_STRING: u8 = 8;

// ---------------------------------------------------------------------------
// RadWriter — binary buffer
// ---------------------------------------------------------------------------

/// Binary buffer for RAD format output.
///
/// Accumulates bytes in an internal buffer, then flushes to a writer.
pub struct RadWriter {
    buf: Vec<u8>,
}

impl RadWriter {
    /// Create a new empty writer.
    pub fn new() -> Self {
        Self {
            buf: Vec::with_capacity(4096),
        }
    }

    /// Number of bytes in the buffer.
    #[inline]
    pub fn len(&self) -> usize {
        self.buf.len()
    }

    /// Whether the buffer is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.buf.is_empty()
    }

    /// Write a u8 value.
    #[inline]
    pub fn write_u8(&mut self, v: u8) {
        self.buf.push(v);
    }

    /// Write a u16 value (little-endian).
    #[inline]
    pub fn write_u16(&mut self, v: u16) {
        self.buf.extend_from_slice(&v.to_le_bytes());
    }

    /// Write a u32 value (little-endian).
    #[inline]
    pub fn write_u32(&mut self, v: u32) {
        self.buf.extend_from_slice(&v.to_le_bytes());
    }

    /// Write a u64 value (little-endian).
    #[inline]
    pub fn write_u64(&mut self, v: u64) {
        self.buf.extend_from_slice(&v.to_le_bytes());
    }

    /// Write a length-prefixed string (u16 length + bytes).
    pub fn write_string(&mut self, s: &str) {
        let len = s.len() as u16;
        self.write_u16(len);
        self.buf.extend_from_slice(s.as_bytes());
    }

    /// Write a tag description: length-prefixed name + type byte.
    fn write_tag_desc(&mut self, name: &str, type_id: u8) {
        self.write_string(name);
        self.write_u8(type_id);
    }

    /// Write an array tag description: name + TAG_ARRAY + length_type + elem_type.
    fn write_array_tag_desc(&mut self, name: &str, length_type: u8, elem_type: u8) {
        self.write_string(name);
        self.write_u8(TAG_ARRAY);
        self.write_u8(length_type);
        self.write_u8(elem_type);
    }

    /// Overwrite a u32 at a specific byte offset (for backpatching counts).
    pub fn write_u32_at_offset(&mut self, offset: usize, v: u32) {
        self.buf[offset..offset + 4].copy_from_slice(&v.to_le_bytes());
    }

    /// Overwrite a u64 at a specific byte offset.
    pub fn write_u64_at_offset(&mut self, offset: usize, v: u64) {
        self.buf[offset..offset + 8].copy_from_slice(&v.to_le_bytes());
    }

    /// Flush the buffer contents to a writer and clear the buffer.
    pub fn flush_to<W: Write>(&mut self, writer: &mut W) -> Result<()> {
        writer.write_all(&self.buf)?;
        self.buf.clear();
        Ok(())
    }

    /// Clear the buffer without writing.
    pub fn clear(&mut self) {
        self.buf.clear();
    }

    /// Access the underlying buffer.
    pub fn as_bytes(&self) -> &[u8] {
        &self.buf
    }
}

impl Default for RadWriter {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// 2-bit base packing (for barcode/UMI)
// ---------------------------------------------------------------------------

/// Pack a nucleotide sequence into a 2-bit encoded u64 (MSB-first).
///
/// Uses encoding A=0, C=1, G=2, T=3 (matching C++ `Kmer::fromChars()`).
/// This is **different** from the sshash-rs encoding (which uses A=0, C=1,
/// G=3, T=2).
///
/// Sequences longer than 32 bases cannot be packed into u64.
pub fn pack_bases_2bit(seq: &[u8]) -> u64 {
    let mut packed: u64 = 0;
    let k = seq.len();
    for (i, &base) in seq.iter().enumerate() {
        let code: u64 = match base {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => 0, // N handled by barcode recovery
        };
        let shift = 2 * (k - 1 - i); // MSB-first
        packed |= code << shift;
    }
    packed
}

// ---------------------------------------------------------------------------
// RAD header writers
// ---------------------------------------------------------------------------

/// Write a RAD header for single-cell mode.
///
/// Returns `(chunk_count_offset, Option<read_length_offset>)` for backpatching.
///
/// The chunk_count_offset is where the u64 num_chunks placeholder is.
/// The read_length_offset (only present when with_position=true) is where
/// the u32 read_length placeholder is.
pub fn write_rad_header_sc<W: Write>(
    writer: &mut W,
    num_refs: u64,
    ref_names: &[&str],
    bc_len: u16,
    umi_len: u16,
    with_position: bool,
) -> Result<(u64, Option<u64>)> {
    let mut buf = RadWriter::new();

    // is_paired flag (always 0 for SC)
    buf.write_u8(0);

    // Reference count and names
    buf.write_u64(num_refs);
    for name in ref_names {
        buf.write_string(name);
    }

    // num_chunks placeholder
    let chunk_count_offset = buf.len() as u64;
    buf.write_u64(0);

    // --- Tag descriptions ---

    // File-level tags: 3 (no position) or 4 (with position)
    let num_file_tags: u16 = if with_position { 4 } else { 3 };
    buf.write_u16(num_file_tags);
    buf.write_tag_desc("cblen", TAG_U16);
    buf.write_tag_desc("ulen", TAG_U16);
    buf.write_tag_desc("known_rad_type", TAG_STRING);
    if with_position {
        buf.write_tag_desc("rlen", TAG_U32);
    }

    // Read-level tags: 2 (barcode + UMI)
    buf.write_u16(2);
    // Barcode: u32 if bc_len <= 16, u64 if bc_len <= 32
    let bc_tag_type = if bc_len <= 16 { TAG_U32 } else { TAG_U64 };
    buf.write_tag_desc("b", bc_tag_type);
    // UMI: u32 if umi_len <= 16, u64 if umi_len <= 32
    let umi_tag_type = if umi_len <= 16 { TAG_U32 } else { TAG_U64 };
    buf.write_tag_desc("u", umi_tag_type);

    // Alignment-level tags: 1 (no position) or 2 (with position)
    let num_aln_tags: u16 = if with_position { 2 } else { 1 };
    buf.write_u16(num_aln_tags);
    buf.write_tag_desc("compressed_ori_refid", TAG_U32);
    if with_position {
        buf.write_tag_desc("pos", TAG_U32);
    }

    // --- File-level tag values ---
    buf.write_u16(bc_len);
    buf.write_u16(umi_len);
    // known_rad_type string
    let rad_type_str = if with_position {
        "sc_rna_pos"
    } else {
        "sc_rna_basic"
    };
    buf.write_string(rad_type_str);

    // read_length placeholder (only with_position)
    let read_length_offset = if with_position {
        let off = buf.len() as u64;
        buf.write_u32(0);
        Some(off)
    } else {
        None
    };

    buf.flush_to(writer)?;
    Ok((chunk_count_offset, read_length_offset))
}

/// Write a RAD header for bulk mode.
///
/// Returns `chunk_count_offset` for backpatching.
pub fn write_rad_header_bulk<W: Write>(
    writer: &mut W,
    is_paired: bool,
    num_refs: u64,
    ref_names: &[&str],
    ref_lengths: &[u32],
) -> Result<u64> {
    let mut buf = RadWriter::new();

    // is_paired flag
    buf.write_u8(if is_paired { 1 } else { 0 });

    // Reference count and names
    buf.write_u64(num_refs);
    for name in ref_names {
        buf.write_string(name);
    }

    // num_chunks placeholder
    let chunk_count_offset = buf.len() as u64;
    buf.write_u64(0);

    // --- Tag descriptions ---

    // File-level tags: 2
    buf.write_u16(2);
    buf.write_tag_desc("known_rad_type", TAG_STRING);
    buf.write_array_tag_desc("ref_lengths", TAG_U32, TAG_U32);

    // Read-level tags: 1
    buf.write_u16(1);
    buf.write_tag_desc("frag_map_type", TAG_U8);

    // Alignment-level tags: 3
    buf.write_u16(3);
    buf.write_tag_desc("compressed_ori_ref", TAG_U32);
    buf.write_tag_desc("pos", TAG_U32);
    buf.write_tag_desc("frag_len", TAG_U16);

    // --- File-level tag values ---
    buf.write_string("bulk_with_pos");
    // ref_lengths array: u32 count + u32 per ref
    buf.write_u32(ref_lengths.len() as u32);
    for &len in ref_lengths {
        buf.write_u32(len);
    }

    buf.flush_to(writer)?;
    Ok(chunk_count_offset)
}

// ---------------------------------------------------------------------------
// RAD record writers
// ---------------------------------------------------------------------------

/// Write a single-cell RAD record.
///
/// Barcode and UMI are pre-packed as 2-bit encoded values.
/// `bc_len` and `umi_len` are in bases (for choosing u32 vs u64 encoding).
///
/// Orientation encoding for `compressed_ori_refid`:
/// - `SingleMapped` / `MappedFirstOrphan` / `MappedPair`: `tid | (is_fw ? 0x80000000 : 0)`
/// - `MappedSecondOrphan`: `tid | (is_fw ? 0 : 0x80000000)` (inverted!)
#[allow(clippy::too_many_arguments)]
pub fn write_sc_record(
    bc_packed: u64,
    umi_packed: u64,
    bc_len: u16,
    umi_len: u16,
    map_type: MappingType,
    accepted_hits: &[SimpleHit],
    with_position: bool,
    writer: &mut RadWriter,
) {
    // Number of mappings
    writer.write_u32(accepted_hits.len() as u32);

    // Barcode (u32 or u64 depending on length)
    if bc_len <= 16 {
        writer.write_u32(bc_packed as u32);
    } else {
        writer.write_u64(bc_packed);
    }

    // UMI (u32 or u64 depending on length)
    if umi_len <= 16 {
        writer.write_u32(umi_packed as u32);
    } else {
        writer.write_u64(umi_packed);
    }

    // Per-alignment data
    let invert_ori = map_type == MappingType::MappedSecondOrphan;
    for hit in accepted_hits {
        let fw_mask: u32 = if invert_ori {
            // MappedSecondOrphan: orientation is inverted
            if hit.is_fw { 0 } else { 0x80000000 }
        } else if hit.is_fw {
            0x80000000
        } else {
            0
        };
        let compressed = hit.tid | fw_mask;
        writer.write_u32(compressed);

        if with_position {
            writer.write_u32(hit.pos as u32);
        }
    }
}

/// Write a bulk RAD record.
///
/// Orientation encoding:
/// `compressed_ori_ref = (tid & 0x3FFFFFFF) | fw_mask | mate_fw_mask`
///
/// Position computation:
/// - Single/orphan: `max(0, pos)`
/// - Pair: `min(pos, mate_pos)`, clamped to 0 (adjusting frag_len)
pub fn write_bulk_record(
    map_type: MappingType,
    accepted_hits: &[SimpleHit],
    writer: &mut RadWriter,
) {
    // Number of mappings
    writer.write_u32(accepted_hits.len() as u32);

    // Fragment mapping type
    writer.write_u8(map_type as u8);

    for hit in accepted_hits {
        let fw_mask: u32 = if hit.is_fw { 0x80000000 } else { 0 };
        let mate_fw_mask: u32 = if hit.mate_is_fw { 0x40000000 } else { 0 };
        let compressed = (hit.tid & 0x3FFFFFFF) | fw_mask | mate_fw_mask;
        writer.write_u32(compressed);

        // Position and fragment length
        let is_pair = map_type == MappingType::MappedPair;
        if is_pair {
            let mut leftmost = hit.pos.min(hit.mate_pos);
            let mut frag_len = hit.fragment_length;
            if leftmost < 0 {
                frag_len += leftmost;
                leftmost = 0;
            }
            writer.write_u32(leftmost as u32);
            writer.write_u16(frag_len as u16);
        } else {
            let leftmost = hit.pos.max(0);
            writer.write_u32(leftmost as u32);
            writer.write_u16(u16::MAX); // no fragment length for non-pairs
        }
    }
}

// ---------------------------------------------------------------------------
// ATAC RAD header + record writers
// ---------------------------------------------------------------------------

/// Write a RAD header for scATAC mode.
///
/// Returns `chunk_count_offset` for backpatching.
///
/// ATAC header format:
/// - is_paired (u8, always 1 for ATAC)
/// - num_refs (u64) + ref_names
/// - num_chunks placeholder (u64)
/// - File-level tags: cblen (u16), known_rad_type (string), ref_lengths (array u32)
/// - Read-level tags: barcode (u32 or u64)
/// - Alignment-level tags: ref (u32), type (u8), start_pos (u32), frag_len (u16)
pub fn write_rad_header_atac<W: Write>(
    writer: &mut W,
    num_refs: u64,
    ref_names: &[&str],
    ref_lengths: &[u32],
    bc_len: u16,
) -> Result<u64> {
    let mut buf = RadWriter::new();

    // is_paired=1 (ATAC is always paired)
    buf.write_u8(1);

    // Reference count and names
    buf.write_u64(num_refs);
    for name in ref_names {
        buf.write_string(name);
    }

    // num_chunks placeholder
    let chunk_count_offset = buf.len() as u64;
    buf.write_u64(0);

    // --- Tag descriptions ---

    // File-level tags: 3 (cblen + known_rad_type + ref_lengths)
    buf.write_u16(3);
    buf.write_tag_desc("cblen", TAG_U16);
    buf.write_tag_desc("known_rad_type", TAG_STRING);
    buf.write_array_tag_desc("ref_lengths", TAG_U32, TAG_U32);

    // Read-level tags: 1 (barcode)
    buf.write_u16(1);
    let bc_tag_type = if bc_len <= 16 { TAG_U32 } else { TAG_U64 };
    buf.write_tag_desc("b", bc_tag_type);

    // Alignment-level tags: 4 (ref, type, start_pos, frag_len)
    buf.write_u16(4);
    buf.write_tag_desc("ref", TAG_U32);
    buf.write_tag_desc("type", TAG_U8);
    buf.write_tag_desc("start_pos", TAG_U32);
    buf.write_tag_desc("frag_len", TAG_U16);

    // --- File-level tag values ---
    buf.write_u16(bc_len);
    buf.write_string("sc_atac");
    // ref_lengths array
    buf.write_u32(ref_lengths.len() as u32);
    for &len in ref_lengths {
        buf.write_u32(len);
    }

    buf.flush_to(writer)?;
    Ok(chunk_count_offset)
}

/// ATAC mapping type codes (written as u8 in the RAD record).
#[repr(u8)]
#[derive(Debug, Clone, Copy)]
pub enum AtacMappingCode {
    Single = 1,
    FirstOrphan = 2,
    SecondOrphan = 3,
    Pair = 4,
    Unmapped = 8,
}

impl From<MappingType> for AtacMappingCode {
    fn from(mt: MappingType) -> Self {
        match mt {
            MappingType::SingleMapped => AtacMappingCode::Single,
            MappingType::MappedFirstOrphan => AtacMappingCode::FirstOrphan,
            MappingType::MappedSecondOrphan => AtacMappingCode::SecondOrphan,
            MappingType::MappedPair => AtacMappingCode::Pair,
            MappingType::Unmapped => AtacMappingCode::Unmapped,
        }
    }
}

/// Write a scATAC RAD record.
///
/// Includes Tn5 shift (+4 to position, -9 from fragment length) when enabled.
pub fn write_atac_record(
    bc_packed: u64,
    bc_len: u16,
    map_type: MappingType,
    accepted_hits: &[SimpleHit],
    tn5_shift: bool,
    writer: &mut RadWriter,
) {
    // Number of mappings
    writer.write_u32(accepted_hits.len() as u32);

    // Barcode (u32 or u64)
    if bc_len <= 16 {
        writer.write_u32(bc_packed as u32);
    } else {
        writer.write_u64(bc_packed);
    }

    let atac_type = AtacMappingCode::from(map_type) as u8;

    for hit in accepted_hits {
        // ref id
        writer.write_u32(hit.tid);
        // type
        writer.write_u8(atac_type);

        let mut leftmost_pos: i32;
        let mut frag_len: u16;

        match map_type {
            MappingType::MappedPair => {
                leftmost_pos = hit.pos.min(hit.mate_pos);
                frag_len = hit.frag_len() as u16;
                if leftmost_pos < 0 {
                    frag_len = (hit.frag_len() + leftmost_pos) as u16;
                    leftmost_pos = 0;
                }
            }
            _ => {
                leftmost_pos = hit.pos.max(0);
                frag_len = u16::MAX; // placeholder for non-pairs
            }
        }

        if tn5_shift {
            leftmost_pos += 4;
            frag_len = frag_len.saturating_sub(9);
        }

        writer.write_u32(leftmost_pos as u32);
        writer.write_u16(frag_len);
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rad_writer_primitives() {
        let mut w = RadWriter::new();
        w.write_u8(0x42);
        w.write_u16(0x1234);
        w.write_u32(0xDEADBEEF);
        w.write_u64(0x0102030405060708);

        assert_eq!(w.len(), 1 + 2 + 4 + 8);

        let b = w.as_bytes();
        assert_eq!(b[0], 0x42);
        assert_eq!(u16::from_le_bytes([b[1], b[2]]), 0x1234);
        assert_eq!(
            u32::from_le_bytes([b[3], b[4], b[5], b[6]]),
            0xDEADBEEF
        );
    }

    #[test]
    fn test_write_at_offset() {
        let mut w = RadWriter::new();
        w.write_u32(0); // placeholder at offset 0
        w.write_u32(0xAAAAAAAA);

        w.write_u32_at_offset(0, 0xBBBBBBBB);

        let b = w.as_bytes();
        assert_eq!(
            u32::from_le_bytes([b[0], b[1], b[2], b[3]]),
            0xBBBBBBBB
        );
        assert_eq!(
            u32::from_le_bytes([b[4], b[5], b[6], b[7]]),
            0xAAAAAAAA
        );
    }

    #[test]
    fn test_write_string() {
        let mut w = RadWriter::new();
        w.write_string("hello");
        let b = w.as_bytes();
        assert_eq!(u16::from_le_bytes([b[0], b[1]]), 5);
        assert_eq!(&b[2..7], b"hello");
    }

    #[test]
    fn test_pack_bases_2bit_encoding() {
        // A=0, C=1, G=2, T=3, MSB-first
        assert_eq!(pack_bases_2bit(b"A"), 0b00);
        assert_eq!(pack_bases_2bit(b"C"), 0b01);
        assert_eq!(pack_bases_2bit(b"G"), 0b10);
        assert_eq!(pack_bases_2bit(b"T"), 0b11);

        // "ACGT" → 00_01_10_11 = 0x1B
        assert_eq!(pack_bases_2bit(b"ACGT"), 0b00011011);

        // "TTTT" → 11_11_11_11 = 0xFF
        assert_eq!(pack_bases_2bit(b"TTTT"), 0xFF);

        // "AAAA" → 0
        assert_eq!(pack_bases_2bit(b"AAAA"), 0);

        // Case insensitive
        assert_eq!(pack_bases_2bit(b"acgt"), pack_bases_2bit(b"ACGT"));
    }

    #[test]
    fn test_pack_bases_2bit_16mer() {
        // 16-base barcode fits in u32 (32 bits)
        let bc = b"ACGTACGTACGTACGT";
        let packed = pack_bases_2bit(bc);
        // Should be non-zero and fit in 32 bits
        assert!(packed > 0);
        assert!(packed <= u32::MAX as u64);
    }

    #[test]
    fn test_sc_header_no_magic() {
        let mut buf = Vec::new();
        let names = vec!["ref1", "ref2"];
        let (chunk_off, rlen_off) =
            write_rad_header_sc(&mut buf, 2, &names, 16, 12, false).unwrap();

        // First byte should be is_paired=0, not 'R'
        assert_eq!(buf[0], 0);
        // No magic bytes
        assert_ne!(&buf[0..4], b"RAD\x01");
        // chunk_count_offset should point to a valid u64
        assert!(chunk_off > 0);
        // No read_length_offset when with_position=false
        assert!(rlen_off.is_none());
    }

    #[test]
    fn test_sc_header_with_position() {
        let mut buf = Vec::new();
        let names = vec!["ref1"];
        let (_, rlen_off) =
            write_rad_header_sc(&mut buf, 1, &names, 16, 12, true).unwrap();

        // Should have read_length_offset when with_position=true
        assert!(rlen_off.is_some());
    }

    #[test]
    fn test_sc_record_packed_bc_umi() {
        let mut w = RadWriter::new();
        let hits = vec![
            SimpleHit {
                is_fw: true,
                tid: 42,
                pos: 100,
                ..SimpleHit::default()
            },
        ];
        let bc_packed = pack_bases_2bit(b"ACGTACGTACGTACGT"); // 16-base BC
        let umi_packed = pack_bases_2bit(b"ACGTACGTACGT"); // 12-base UMI

        write_sc_record(
            bc_packed, umi_packed, 16, 12,
            MappingType::SingleMapped, &hits, false, &mut w,
        );

        let b = w.as_bytes();
        // num_mappings (4) + bc_packed u32 (4) + umi_packed u32 (4) + compressed (4) = 16
        assert_eq!(b.len(), 16);

        // num_mappings = 1
        assert_eq!(u32::from_le_bytes([b[0], b[1], b[2], b[3]]), 1);

        // bc_packed as u32 (bc_len=16 ≤ 16)
        let bc = u32::from_le_bytes([b[4], b[5], b[6], b[7]]);
        assert_eq!(bc, bc_packed as u32);

        // umi_packed as u32 (umi_len=12 ≤ 16)
        let umi = u32::from_le_bytes([b[8], b[9], b[10], b[11]]);
        assert_eq!(umi, umi_packed as u32);

        // compressed: tid=42, is_fw=true → 42 | 0x80000000
        let comp = u32::from_le_bytes([b[12], b[13], b[14], b[15]]);
        assert_eq!(comp, 42 | 0x80000000);
    }

    #[test]
    fn test_sc_second_orphan_inversion() {
        let mut w = RadWriter::new();
        let hits = vec![
            SimpleHit {
                is_fw: true,
                tid: 10,
                pos: 50,
                ..SimpleHit::default()
            },
        ];

        write_sc_record(
            0, 0, 16, 12,
            MappingType::MappedSecondOrphan, &hits, false, &mut w,
        );

        let b = w.as_bytes();
        // For MappedSecondOrphan, orientation is inverted:
        // is_fw=true → fw_mask=0 (not 0x80000000)
        let comp = u32::from_le_bytes([b[12], b[13], b[14], b[15]]);
        assert_eq!(comp, 10); // tid=10, no fw_mask

        // Now test with is_fw=false → should get 0x80000000
        let mut w2 = RadWriter::new();
        let hits2 = vec![
            SimpleHit {
                is_fw: false,
                tid: 10,
                pos: 50,
                ..SimpleHit::default()
            },
        ];
        write_sc_record(
            0, 0, 16, 12,
            MappingType::MappedSecondOrphan, &hits2, false, &mut w2,
        );
        let b2 = w2.as_bytes();
        let comp2 = u32::from_le_bytes([b2[12], b2[13], b2[14], b2[15]]);
        assert_eq!(comp2, 10 | 0x80000000);
    }

    #[test]
    fn test_bulk_record_leftmost_pos() {
        let mut w = RadWriter::new();
        // Single-mapped with negative pos → clamped to 0
        let hits = vec![SimpleHit {
            is_fw: true,
            tid: 5,
            pos: -10,
            ..SimpleHit::default()
        }];
        write_bulk_record(MappingType::SingleMapped, &hits, &mut w);

        let b = w.as_bytes();
        // num_mappings (4) + frag_map_type (1) + compressed (4) + pos (4) + frag_len (2) = 15
        assert_eq!(b.len(), 15);

        // pos should be clamped to 0
        let pos = u32::from_le_bytes([b[9], b[10], b[11], b[12]]);
        assert_eq!(pos, 0);

        // frag_len should be u16::MAX for non-pair
        let frag = u16::from_le_bytes([b[13], b[14]]);
        assert_eq!(frag, u16::MAX);
    }

    #[test]
    fn test_bulk_record_pair_frag_len() {
        let mut w = RadWriter::new();
        let hits = vec![SimpleHit {
            is_fw: true,
            mate_is_fw: false,
            tid: 3,
            pos: 100,
            mate_pos: 200,
            fragment_length: 250,
            ..SimpleHit::default()
        }];
        write_bulk_record(MappingType::MappedPair, &hits, &mut w);

        let b = w.as_bytes();
        // num_mappings (4) + frag_map_type (1) + compressed (4) + pos (4) + frag_len (2) = 15
        assert_eq!(b.len(), 15);

        // frag_map_type should be MappedPair = 4
        assert_eq!(b[4], MappingType::MappedPair as u8);

        // leftmost_pos = min(100, 200) = 100
        let pos = u32::from_le_bytes([b[9], b[10], b[11], b[12]]);
        assert_eq!(pos, 100);

        // frag_len = 250
        let frag = u16::from_le_bytes([b[13], b[14]]);
        assert_eq!(frag, 250);
    }

    #[test]
    fn test_bulk_header_format() {
        let mut buf = Vec::new();
        let names = vec!["ref1", "ref2"];
        let ref_lens = vec![1000u32, 2000u32];
        let chunk_off = write_rad_header_bulk(&mut buf, true, 2, &names, &ref_lens).unwrap();

        // First byte should be is_paired=1
        assert_eq!(buf[0], 1);
        // No RAD magic
        assert_ne!(&buf[0..4], b"RAD\x01");
        assert!(chunk_off > 0);
    }

    #[test]
    fn test_atac_header() {
        let mut buf = Vec::new();
        let names = vec!["chr1", "chr2"];
        let ref_lens = vec![248956422u32, 242193529u32];
        let chunk_off =
            write_rad_header_atac(&mut buf, 2, &names, &ref_lens, 16).unwrap();

        // First byte: is_paired=1
        assert_eq!(buf[0], 1);
        assert!(chunk_off > 0);
    }

    #[test]
    fn test_atac_record_basic() {
        let mut w = RadWriter::new();
        let hits = vec![SimpleHit {
            is_fw: true,
            tid: 0,
            pos: 1000,
            ..SimpleHit::default()
        }];
        write_atac_record(0xABCD, 16, MappingType::SingleMapped, &hits, false, &mut w);

        let b = w.as_bytes();
        // num_mappings (4) + bc u32 (4) + ref (4) + type (1) + pos (4) + frag_len (2) = 19
        assert_eq!(b.len(), 19);

        // type = 1 (Single)
        assert_eq!(b[12], 1);

        // pos = 1000
        let pos = u32::from_le_bytes([b[13], b[14], b[15], b[16]]);
        assert_eq!(pos, 1000);
    }

    #[test]
    fn test_atac_record_tn5_shift() {
        let mut w = RadWriter::new();
        let hits = vec![SimpleHit {
            is_fw: true,
            tid: 0,
            pos: 100,
            mate_pos: 300,
            fragment_length: 250,
            ..SimpleHit::default()
        }];
        write_atac_record(0, 16, MappingType::MappedPair, &hits, true, &mut w);

        let b = w.as_bytes();
        // leftmost = min(100, 300) = 100, + Tn5 shift 4 = 104
        let pos = u32::from_le_bytes([b[13], b[14], b[15], b[16]]);
        assert_eq!(pos, 104);

        // frag_len = 250 - 9 = 241
        let frag = u16::from_le_bytes([b[17], b[18]]);
        assert_eq!(frag, 241);
    }

    #[test]
    fn test_sc_record_with_position() {
        let mut w = RadWriter::new();
        let hits = vec![
            SimpleHit {
                is_fw: true,
                tid: 7,
                pos: 42,
                ..SimpleHit::default()
            },
        ];

        write_sc_record(
            0, 0, 16, 12,
            MappingType::SingleMapped, &hits, true, &mut w,
        );

        let b = w.as_bytes();
        // num_mappings (4) + bc u32 (4) + umi u32 (4) + compressed (4) + pos (4) = 20
        assert_eq!(b.len(), 20);

        // pos at the end
        let pos = u32::from_le_bytes([b[16], b[17], b[18], b[19]]);
        assert_eq!(pos, 42);
    }
}
