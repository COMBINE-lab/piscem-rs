//! RAD format binary writer.
//!
//! Provides a buffer-based writer for RAD (Reduced Alignment Data) format
//! output, used by salmon/alevin-fry. Includes functions for writing headers
//! and records in both single-cell and bulk modes.
//!
//! Port of C++ `rad_writer` + header/record functions.

use std::io::Write;

use anyhow::Result;

use crate::mapping::hits::{MappingType, SimpleHit};

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

    /// Write an i32 value (little-endian).
    #[inline]
    pub fn write_i32(&mut self, v: i32) {
        self.buf.extend_from_slice(&v.to_le_bytes());
    }

    /// Write a length-prefixed string (u16 length + bytes).
    pub fn write_string(&mut self, s: &str) {
        let len = s.len() as u16;
        self.write_u16(len);
        self.buf.extend_from_slice(s.as_bytes());
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
// RAD header writers
// ---------------------------------------------------------------------------

/// Write a RAD header for single-cell mode.
///
/// Returns `(num_chunks_offset, barcode_count_offset)` for backpatching.
pub fn write_rad_header_sc<W: Write>(
    writer: &mut W,
    num_refs: u32,
    ref_names: &[&str],
    ref_lengths: &[u64],
    barcode_len: u16,
    umi_len: u16,
) -> Result<(u64, Option<u64>)> {
    let mut buf = RadWriter::new();

    // RAD magic "RAD\x01"
    buf.write_u8(b'R');
    buf.write_u8(b'A');
    buf.write_u8(b'D');
    buf.write_u8(1); // version

    // Number of references
    buf.write_u32(num_refs);

    // Reference names (length-prefixed)
    for name in ref_names {
        buf.write_string(name);
    }

    // Reference lengths
    for &len in ref_lengths {
        buf.write_u64(len);
    }

    // Barcode and UMI lengths
    buf.write_u16(barcode_len);
    buf.write_u16(umi_len);

    // Number of chunks placeholder (will be backpatched)
    let num_chunks_offset = buf.len() as u64;
    buf.write_u64(0);

    buf.flush_to(writer)?;
    Ok((num_chunks_offset, None))
}

/// Write a RAD header for bulk mode.
///
/// Returns `num_chunks_offset` for backpatching.
pub fn write_rad_header_bulk<W: Write>(
    writer: &mut W,
    num_refs: u32,
    ref_names: &[&str],
    ref_lengths: &[u64],
) -> Result<u64> {
    let mut buf = RadWriter::new();

    buf.write_u8(b'R');
    buf.write_u8(b'A');
    buf.write_u8(b'D');
    buf.write_u8(1); // version

    buf.write_u32(num_refs);

    for name in ref_names {
        buf.write_string(name);
    }
    for &len in ref_lengths {
        buf.write_u64(len);
    }

    // Number of chunks placeholder
    let num_chunks_offset = buf.len() as u64;
    buf.write_u64(0);

    buf.flush_to(writer)?;
    Ok(num_chunks_offset)
}

// ---------------------------------------------------------------------------
// RAD record writers
// ---------------------------------------------------------------------------

/// Write a single-cell RAD record.
///
/// Orientation encoding: `compressed_ori_refid = tid | (is_fw ? 0x80000000 : 0)`
pub fn write_sc_record(
    barcode: &[u8],
    umi: &[u8],
    map_type: MappingType,
    accepted_hits: &[SimpleHit],
    writer: &mut RadWriter,
) {
    // Barcode bytes
    writer.buf.extend_from_slice(barcode);
    // UMI bytes
    writer.buf.extend_from_slice(umi);

    // Mapping type
    writer.write_u8(map_type as u8);

    // Number of hits
    let num_hits = accepted_hits.len() as u32;
    writer.write_u32(num_hits);

    // Hits: compressed (tid | orientation)
    for hit in accepted_hits {
        let compressed = hit.tid | if hit.is_fw { 0x80000000 } else { 0 };
        writer.write_u32(compressed);
    }
}

/// Write a bulk RAD record.
///
/// Orientation encoding:
/// `compressed_ori_ref = (tid & 0x3FFFFFFF) | fw_mask | mate_fw_mask`
pub fn write_bulk_record(
    map_type: MappingType,
    accepted_hits: &[SimpleHit],
    writer: &mut RadWriter,
) {
    // Mapping type
    writer.write_u8(map_type as u8);

    // Number of hits
    let num_hits = accepted_hits.len() as u32;
    writer.write_u32(num_hits);

    for hit in accepted_hits {
        let fw_mask: u32 = if hit.is_fw { 0x80000000 } else { 0 };
        let mate_fw_mask: u32 = if hit.mate_is_fw { 0x40000000 } else { 0 };
        let compressed = (hit.tid & 0x3FFFFFFF) | fw_mask | mate_fw_mask;
        writer.write_u32(compressed);
        writer.write_i32(hit.pos);
        if hit.has_mate() {
            writer.write_i32(hit.mate_pos);
            writer.write_i32(hit.fragment_length);
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
    fn test_write_sc_record() {
        let mut w = RadWriter::new();
        let barcode = b"\x01\x02\x03\x04";
        let umi = b"\x05\x06";
        let hits = vec![
            SimpleHit {
                is_fw: true,
                tid: 42,
                ..SimpleHit::default()
            },
            SimpleHit {
                is_fw: false,
                tid: 99,
                ..SimpleHit::default()
            },
        ];
        write_sc_record(barcode, umi, MappingType::SingleMapped, &hits, &mut w);

        let b = w.as_bytes();
        // barcode (4) + umi (2) + map_type (1) + num_hits (4) + 2 * compressed (4) = 19
        assert_eq!(b.len(), 19);
        // map_type at offset 6
        assert_eq!(b[6], MappingType::SingleMapped as u8);
        // num_hits at offset 7
        assert_eq!(u32::from_le_bytes([b[7], b[8], b[9], b[10]]), 2);
        // first hit: tid=42, is_fw=true → 42 | 0x80000000
        let h0 = u32::from_le_bytes([b[11], b[12], b[13], b[14]]);
        assert_eq!(h0, 42 | 0x80000000);
        // second hit: tid=99, is_fw=false → 99
        let h1 = u32::from_le_bytes([b[15], b[16], b[17], b[18]]);
        assert_eq!(h1, 99);
    }

    #[test]
    fn test_write_bulk_record() {
        let mut w = RadWriter::new();
        let hits = vec![SimpleHit {
            is_fw: true,
            mate_is_fw: false,
            tid: 10,
            pos: 500,
            ..SimpleHit::default()
        }];
        write_bulk_record(MappingType::SingleMapped, &hits, &mut w);

        let b = w.as_bytes();
        // map_type (1) + num_hits (4) + compressed (4) + pos (4) = 13
        // (no mate — has_mate returns false because mate_pos == INVALID_MATE_POS)
        assert_eq!(b.len(), 13);
        assert_eq!(b[0], MappingType::SingleMapped as u8);
    }

    #[test]
    fn test_write_string() {
        let mut w = RadWriter::new();
        w.write_string("hello");
        let b = w.as_bytes();
        assert_eq!(u16::from_le_bytes([b[0], b[1]]), 5);
        assert_eq!(&b[2..7], b"hello");
    }
}
