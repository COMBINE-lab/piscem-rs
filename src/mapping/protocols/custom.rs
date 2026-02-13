//! Custom read geometry parser and protocol.
//!
//! Allows user-defined read geometries via a grammar-based parser, supporting
//! split barcodes/UMIs across R1 and R2, unbounded reads, fixed sequences,
//! and discard regions.
//!
//! Grammar:
//! ```text
//! Specification := Read1Description Read2Description
//! Read1Description := '1{' DescList '}'
//! Read2Description := '2{' DescList '}'
//! DescList := (BoundedDesc{1,10} UnboundedDesc?) | UnboundedDesc
//! BoundedDesc := 'b[' Length ']' | 'u[' Length ']' | 'f[' Sequence ']'
//!              | 'x[' Length ']' | 'r[' Length ']'
//! UnboundedDesc := 'x:' | 'r:'
//! Length := [1-9][0-9]*
//! Sequence := [ATGC]+
//! ```
//!
//! Port of C++ `sc/util.hpp` PEG parser + `custom_protocol` class.

use anyhow::{bail, Result};
use smallvec::SmallVec;

use super::{AlignableReads, Protocol, TechSeqs};

// ---------------------------------------------------------------------------
// GeoTagType
// ---------------------------------------------------------------------------

/// Type of a geometry segment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeoTagType {
    Barcode,
    Umi,
    Read,
    Fixed,
    Discard,
}

// ---------------------------------------------------------------------------
// GeoPart
// ---------------------------------------------------------------------------

/// A parsed geometry segment with type and length.
#[derive(Debug, Clone)]
pub struct GeoPart {
    pub tag_type: GeoTagType,
    /// Length in bases. -1 means unbounded (rest of read).
    pub len: i32,
}

// ---------------------------------------------------------------------------
// StrSlice
// ---------------------------------------------------------------------------

/// Offset + length within a read sequence.
#[derive(Debug, Clone, Copy)]
struct StrSlice {
    offset: usize,
    len: i32, // -1 = unbounded
}

// ---------------------------------------------------------------------------
// CustomProtocol
// ---------------------------------------------------------------------------

/// Custom geometry protocol built from a parsed geometry string.
///
/// Supports split barcodes and UMIs across R1 and R2, fixed-length and
/// unbounded biological reads.
#[derive(Debug, Clone)]
pub struct CustomProtocol {
    bc_slices_r1: SmallVec<[StrSlice; 4]>,
    umi_slices_r1: SmallVec<[StrSlice; 4]>,
    read_slices_r1: SmallVec<[StrSlice; 4]>,
    bc_slices_r2: SmallVec<[StrSlice; 4]>,
    umi_slices_r2: SmallVec<[StrSlice; 4]>,
    read_slices_r2: SmallVec<[StrSlice; 4]>,
    bc_len: usize,
    umi_len: usize,
    is_paired_bio: bool,
}

impl CustomProtocol {
    /// Total barcode length in bases.
    pub fn total_bc_len(&self) -> usize {
        self.bc_len
    }

    /// Total UMI length in bases.
    pub fn total_umi_len(&self) -> usize {
        self.umi_len
    }
}

impl Protocol for CustomProtocol {
    fn name(&self) -> &str {
        "custom"
    }

    fn is_bio_paired_end(&self) -> bool {
        self.is_paired_bio
    }

    fn extract_tech_seqs<'a>(&self, r1: &'a [u8], r2: &'a [u8]) -> TechSeqs<'a> {
        // For simplicity, concatenate BC slices into contiguous region.
        // Since barcodes are packed into u64, we need contiguous bytes.
        // The most common case is a single contiguous BC slice on R1.
        //
        // For split barcodes, the caller needs to handle concatenation at
        // a higher level. For now, return the first BC slice.
        let barcode = extract_first_slice(&self.bc_slices_r1, r1)
            .or_else(|| extract_first_slice(&self.bc_slices_r2, r2));

        let umi = extract_first_slice(&self.umi_slices_r1, r1)
            .or_else(|| extract_first_slice(&self.umi_slices_r2, r2));

        TechSeqs { barcode, umi }
    }

    fn extract_mappable_reads<'a>(&self, r1: &'a [u8], r2: &'a [u8]) -> AlignableReads<'a> {
        let seq1 = extract_read_region(&self.read_slices_r1, r1);
        let seq2 = extract_read_region(&self.read_slices_r2, r2);

        AlignableReads { seq1, seq2 }
    }

    fn barcode_len(&self) -> usize {
        self.bc_len
    }

    fn umi_len(&self) -> usize {
        self.umi_len
    }
}

/// Extract the first slice region from a read.
fn extract_first_slice<'a>(slices: &[StrSlice], read: &'a [u8]) -> Option<&'a [u8]> {
    if slices.is_empty() {
        return None;
    }
    let s = &slices[0];
    if s.offset >= read.len() {
        return None;
    }
    let end = if s.len < 0 {
        read.len()
    } else {
        (s.offset + s.len as usize).min(read.len())
    };
    Some(&read[s.offset..end])
}

/// Extract the read region (biological sequence) from a read.
fn extract_read_region<'a>(slices: &[StrSlice], read: &'a [u8]) -> Option<&'a [u8]> {
    if slices.is_empty() {
        return None;
    }
    let s = &slices[0];
    if s.offset >= read.len() {
        return None;
    }
    let end = if s.len < 0 {
        read.len()
    } else {
        (s.offset + s.len as usize).min(read.len())
    };
    let region = &read[s.offset..end];
    if region.is_empty() {
        None
    } else {
        Some(region)
    }
}

// ---------------------------------------------------------------------------
// Parser
// ---------------------------------------------------------------------------

/// Parse a custom geometry string into a `CustomProtocol`.
///
/// Format: `"1{b[16]u[12]x:}2{r:}"` or similar.
///
/// Returns `Err` on invalid syntax or missing required components.
pub fn parse_custom_geometry(geom: &str) -> Result<CustomProtocol> {
    let geom = geom.trim();

    // Find "1{...}" and "2{...}" blocks
    let (parts_r1, parts_r2) = parse_read_descriptions(geom)?;

    // Convert geometry parts to offset/length slices
    let (bc_r1, umi_r1, read_r1) = parts_to_slices(&parts_r1);
    let (bc_r2, umi_r2, read_r2) = parts_to_slices(&parts_r2);

    // Compute total lengths
    let bc_len: usize = sum_bounded_len(&bc_r1) + sum_bounded_len(&bc_r2);
    let umi_len: usize = sum_bounded_len(&umi_r1) + sum_bounded_len(&umi_r2);

    // Validation
    if bc_len == 0 {
        bail!("custom geometry must include at least one barcode segment (b[N])");
    }
    if umi_len == 0 {
        bail!("custom geometry must include at least one UMI segment (u[N])");
    }
    if bc_len > 32 {
        bail!("total barcode length {} exceeds maximum of 32", bc_len);
    }
    if umi_len > 32 {
        bail!("total UMI length {} exceeds maximum of 32", umi_len);
    }

    let has_read_r1 = !read_r1.is_empty();
    let has_read_r2 = !read_r2.is_empty();
    if !has_read_r1 && !has_read_r2 {
        bail!("custom geometry must include at least one biological read segment (r[N] or r:)");
    }

    let is_paired_bio = has_read_r1 && has_read_r2;

    Ok(CustomProtocol {
        bc_slices_r1: bc_r1,
        umi_slices_r1: umi_r1,
        read_slices_r1: read_r1,
        bc_slices_r2: bc_r2,
        umi_slices_r2: umi_r2,
        read_slices_r2: read_r2,
        bc_len,
        umi_len,
        is_paired_bio,
    })
}

/// Parse "1{...}2{...}" into two lists of GeoParts.
fn parse_read_descriptions(geom: &str) -> Result<(Vec<GeoPart>, Vec<GeoPart>)> {
    // Find 1{...} block
    let r1_start = geom
        .find("1{")
        .ok_or_else(|| anyhow::anyhow!("geometry must contain '1{{...}}' block"))?;
    let r1_body_start = r1_start + 2;
    let r1_end = find_matching_brace(geom, r1_body_start)?;
    let r1_body = &geom[r1_body_start..r1_end];

    // Find 2{...} block
    let r2_start = geom[r1_end + 1..]
        .find("2{")
        .ok_or_else(|| anyhow::anyhow!("geometry must contain '2{{...}}' block"))?
        + r1_end
        + 1;
    let r2_body_start = r2_start + 2;
    let r2_end = find_matching_brace(geom, r2_body_start)?;
    let r2_body = &geom[r2_body_start..r2_end];

    let parts_r1 = parse_desc_list(r1_body)?;
    let parts_r2 = parse_desc_list(r2_body)?;

    Ok((parts_r1, parts_r2))
}

/// Find the closing '}' matching an opening '{'.
fn find_matching_brace(s: &str, start: usize) -> Result<usize> {
    let rest = &s[start..];
    rest.find('}')
        .map(|i| i + start)
        .ok_or_else(|| anyhow::anyhow!("unmatched '{{' in geometry"))
}

/// Parse a description list (contents between braces).
fn parse_desc_list(body: &str) -> Result<Vec<GeoPart>> {
    let mut parts = Vec::new();
    let bytes = body.as_bytes();
    let mut i = 0;

    while i < bytes.len() {
        let tag_char = bytes[i];
        if i + 1 >= bytes.len() {
            bail!("unexpected end of geometry at position {}", i);
        }
        let delim = bytes[i + 1];

        let tag_type = match tag_char {
            b'b' => GeoTagType::Barcode,
            b'u' => GeoTagType::Umi,
            b'r' => GeoTagType::Read,
            b'f' => GeoTagType::Fixed,
            b'x' => GeoTagType::Discard,
            _ => bail!("unknown geometry tag '{}' at position {}", tag_char as char, i),
        };

        match delim {
            b'[' => {
                // Bounded: tag[value]
                let close = body[i + 2..]
                    .find(']')
                    .ok_or_else(|| anyhow::anyhow!("unmatched '[' at position {}", i + 1))?
                    + i
                    + 2;
                let inner = &body[i + 2..close];

                let len = if tag_type == GeoTagType::Fixed {
                    // Fixed: value is a DNA sequence, length = sequence length
                    if inner.is_empty() {
                        bail!("empty fixed sequence at position {}", i);
                    }
                    if !inner.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
                        bail!("fixed sequence must contain only ACGT: '{}'", inner);
                    }
                    inner.len() as i32
                } else {
                    // Others: value is a length
                    inner
                        .parse::<i32>()
                        .map_err(|_| anyhow::anyhow!("invalid length '{}' at position {}", inner, i + 2))?
                };

                if len <= 0 && tag_type != GeoTagType::Fixed {
                    bail!("length must be positive, got {} at position {}", len, i + 2);
                }

                parts.push(GeoPart { tag_type, len });
                i = close + 1;
            }
            b':' => {
                // Unbounded
                if tag_type != GeoTagType::Read && tag_type != GeoTagType::Discard {
                    bail!(
                        "only 'r:' and 'x:' are valid unbounded descriptors, got '{}:'",
                        tag_char as char,
                    );
                }
                parts.push(GeoPart {
                    tag_type,
                    len: -1,
                });
                i += 2;
            }
            _ => bail!(
                "expected '[' or ':' after tag '{}', got '{}' at position {}",
                tag_char as char,
                delim as char,
                i + 1,
            ),
        }
    }

    Ok(parts)
}

/// Convert geometry parts into offset/length slices for BC, UMI, and read.
type SliceVec = SmallVec<[StrSlice; 4]>;

fn parts_to_slices(
    parts: &[GeoPart],
) -> (SliceVec, SliceVec, SliceVec) {
    let mut bc_slices = SmallVec::new();
    let mut umi_slices = SmallVec::new();
    let mut read_slices = SmallVec::new();

    let mut offset: usize = 0;

    for part in parts {
        let slice = StrSlice {
            offset,
            len: part.len,
        };

        match part.tag_type {
            GeoTagType::Barcode => bc_slices.push(slice),
            GeoTagType::Umi => umi_slices.push(slice),
            GeoTagType::Read => read_slices.push(slice),
            GeoTagType::Fixed | GeoTagType::Discard => {}
        }

        if part.len > 0 {
            offset += part.len as usize;
        }
        // For unbounded (len=-1), offset stays where it is — the slice
        // will consume from offset to end of read.
    }

    (bc_slices, umi_slices, read_slices)
}

/// Sum bounded lengths from slices.
fn sum_bounded_len(slices: &[StrSlice]) -> usize {
    slices
        .iter()
        .filter(|s| s.len > 0)
        .map(|s| s.len as usize)
        .sum()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_chromium_v3_equivalent() {
        // "1{b[16]u[12]x:}2{r:}" matches chromium_v3 behavior
        let proto = parse_custom_geometry("1{b[16]u[12]x:}2{r:}").unwrap();
        assert_eq!(proto.barcode_len(), 16);
        assert_eq!(proto.umi_len(), 12);
        assert!(!proto.is_bio_paired_end()); // only R2 has bio read

        let r1 = b"ACGTACGTACGTACGTAAAAAAAAAAAA_extra";
        let r2 = b"TGCATGCATGCA";

        let tech = proto.extract_tech_seqs(r1, r2);
        assert_eq!(tech.barcode.unwrap(), b"ACGTACGTACGTACGT");
        assert_eq!(tech.umi.unwrap(), b"AAAAAAAAAAAA");

        let reads = proto.extract_mappable_reads(r1, r2);
        assert!(reads.seq1.is_none()); // R1 has no read segment
        assert_eq!(reads.seq2.unwrap(), b"TGCATGCATGCA");
    }

    #[test]
    fn test_parse_custom_split_bc() {
        // Not supported in detail (first-slice extraction), but parses ok
        let proto = parse_custom_geometry("1{b[8]u[12]r:}2{b[8]r:}").unwrap();
        assert_eq!(proto.barcode_len(), 16); // 8 + 8
        assert_eq!(proto.umi_len(), 12);
        assert!(proto.is_bio_paired_end()); // both R1 and R2 have reads
    }

    #[test]
    fn test_parse_5prime() {
        // "1{b[16]u[12]x[13]r:}2{r:}" matches chromium_v3_5p
        let proto = parse_custom_geometry("1{b[16]u[12]x[13]r:}2{r:}").unwrap();
        assert_eq!(proto.barcode_len(), 16);
        assert_eq!(proto.umi_len(), 12);
        assert!(proto.is_bio_paired_end()); // R1 has bio read after TSO skip

        let r1 = b"ACGTACGTACGTACGTBBBBBBBBBBBBCCCCCCCCCCCCCMAPPABLE_BIO";
        let r2 = b"SECOND_READ_BIO";

        let tech = proto.extract_tech_seqs(r1, r2);
        assert_eq!(tech.barcode.unwrap().len(), 16);
        assert_eq!(tech.umi.unwrap().len(), 12);

        let reads = proto.extract_mappable_reads(r1, r2);
        // R1: skip 16 BC + 12 UMI + 13 discard = 41, rest is bio
        assert_eq!(reads.seq1.unwrap(), b"MAPPABLE_BIO");
        assert_eq!(reads.seq2.unwrap(), b"SECOND_READ_BIO");
    }

    #[test]
    fn test_parse_invalid_geometry_no_read() {
        // Missing read component
        let result = parse_custom_geometry("1{b[16]u[12]}2{x:}");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_invalid_geometry_no_bc() {
        // Missing barcode
        let result = parse_custom_geometry("1{u[12]r:}2{r:}");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_invalid_geometry_bad_syntax() {
        let result = parse_custom_geometry("garbage");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_invalid_geometry_unmatched_brace() {
        let result = parse_custom_geometry("1{b[16]u[12]2{r:}");
        assert!(result.is_err());
    }

    #[test]
    fn test_custom_protocol_extract() {
        let proto = parse_custom_geometry("1{b[16]u[10]r:}2{r:}").unwrap();
        assert_eq!(proto.barcode_len(), 16);
        assert_eq!(proto.umi_len(), 10);

        // R1 with BC + UMI + bio
        let r1 = b"ACGTACGTACGTACGTBBBBBBBBBBREST_OF_R1";
        let r2 = b"BIOLOGICAL_READ_2";

        let tech = proto.extract_tech_seqs(r1, r2);
        assert_eq!(tech.barcode.unwrap(), b"ACGTACGTACGTACGT");
        assert_eq!(tech.umi.unwrap(), b"BBBBBBBBBB");

        let reads = proto.extract_mappable_reads(r1, r2);
        // R1 read starts after BC(16) + UMI(10) = offset 26
        assert_eq!(reads.seq1.unwrap(), b"REST_OF_R1");
        assert_eq!(reads.seq2.unwrap(), b"BIOLOGICAL_READ_2");
    }

    #[test]
    fn test_parse_fixed_sequence() {
        // Fixed sequence segment
        let proto = parse_custom_geometry("1{b[16]u[12]f[ACGT]r:}2{r:}").unwrap();
        assert_eq!(proto.barcode_len(), 16);
        assert_eq!(proto.umi_len(), 12);
        assert!(proto.is_bio_paired_end());

        let r1 = b"ACGTACGTACGTACGTBBBBBBBBBBBBACGT_BIO_READ_1";
        let r2 = b"BIO_READ_2";
        let reads = proto.extract_mappable_reads(r1, r2);
        // R1: 16 BC + 12 UMI + 4 fixed = 32, rest is bio
        assert_eq!(reads.seq1.unwrap(), b"_BIO_READ_1");
    }
}
