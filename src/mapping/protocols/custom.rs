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

use anyhow::{Result, bail};
use smallvec::SmallVec;

use smallvec::smallvec;

use super::{AlignableReads, BarcodeDesc, BarcodeRole, Protocol, TechSeqs};

// ---------------------------------------------------------------------------
// GeoTagType
// ---------------------------------------------------------------------------

/// Type of a geometry segment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeoTagType {
    /// Cell barcode (the `b` tag, or `b` when no `s` is present)
    Barcode,
    /// Sample/library barcode (the `s` tag — syntactic sugar for b0 with Sample role)
    SampleBarcode,
    /// Numbered barcode at a specific level (the `bN` tags, e.g., b0, b1, b2)
    NumberedBarcode(u8),
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
/// unbounded biological reads. Also supports multi-barcode protocols
/// via `s[N]` (sample barcode) or numbered `bN[L]` tags.
#[derive(Debug, Clone)]
pub struct CustomProtocol {
    bc_slices_r1: SmallVec<[StrSlice; 4]>,
    umi_slices_r1: SmallVec<[StrSlice; 4]>,
    read_slices_r1: SmallVec<[StrSlice; 4]>,
    bc_slices_r2: SmallVec<[StrSlice; 4]>,
    umi_slices_r2: SmallVec<[StrSlice; 4]>,
    read_slices_r2: SmallVec<[StrSlice; 4]>,
    /// Per-level sample barcode slices (for multi-barcode protocols).
    /// sample_bc_slices_r1[level] and sample_bc_slices_r2[level] give the
    /// extraction slices for barcode at that level.
    sample_bc_slices_r1: SmallVec<[SmallVec<[StrSlice; 4]>; 2]>,
    sample_bc_slices_r2: SmallVec<[SmallVec<[StrSlice; 4]>; 2]>,
    /// Lengths of each barcode level (index = level, value = length in bases).
    /// For single-barcode protocols, this is empty (bc_len is used instead).
    multi_bc_lens: SmallVec<[usize; 4]>,
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
        let umi = extract_first_slice(&self.umi_slices_r1, r1)
            .or_else(|| extract_first_slice(&self.umi_slices_r2, r2));

        if self.is_multi_barcode() {
            // Multi-barcode: extract each level's barcode
            let mut barcodes = SmallVec::new();
            for (slices_r1, slices_r2) in self
                .sample_bc_slices_r1
                .iter()
                .zip(self.sample_bc_slices_r2.iter())
            {
                let bc = extract_first_slice(slices_r1, r1)
                    .or_else(|| extract_first_slice(slices_r2, r2));
                barcodes.push(bc);
            }
            TechSeqs { barcodes, umi }
        } else {
            // Single-barcode: use the standard bc_slices
            let barcode = extract_first_slice(&self.bc_slices_r1, r1)
                .or_else(|| extract_first_slice(&self.bc_slices_r2, r2));
            TechSeqs {
                barcodes: smallvec![barcode],
                umi,
            }
        }
    }

    fn extract_mappable_reads<'a>(&self, r1: &'a [u8], r2: &'a [u8]) -> AlignableReads<'a> {
        if self.is_paired_bio {
            // Both R1 and R2 have biological read segments
            AlignableReads::Paired {
                read1: extract_read_region(&self.read_slices_r1, r1).unwrap_or(&[]),
                read2: extract_read_region(&self.read_slices_r2, r2).unwrap_or(&[]),
            }
        } else if !self.read_slices_r1.is_empty() {
            // Bio read is on R1 — the geometry specified r[N] or r: in the 1{...} block
            AlignableReads::Single(extract_read_region(&self.read_slices_r1, r1).unwrap_or(&[]))
        } else {
            // Bio read is on R2 — the geometry specified r[N] or r: in the 2{...} block
            AlignableReads::Single(extract_read_region(&self.read_slices_r2, r2).unwrap_or(&[]))
        }
    }

    fn barcode_len(&self) -> usize {
        self.bc_len
    }

    fn umi_len(&self) -> usize {
        self.umi_len
    }

    fn num_barcodes(&self) -> usize {
        if self.multi_bc_lens.is_empty() {
            1
        } else {
            self.multi_bc_lens.len()
        }
    }

    fn barcode_descs(&self) -> Vec<BarcodeDesc> {
        if self.multi_bc_lens.is_empty() {
            // Single-barcode: default behavior
            vec![BarcodeDesc {
                tag_name: "b".to_string(),
                role: BarcodeRole::Cell,
                len: self.bc_len as u16,
            }]
        } else {
            self.multi_bc_lens
                .iter()
                .enumerate()
                .map(|(i, &len)| {
                    let role = if i == 0 && self.multi_bc_lens.len() > 1 {
                        BarcodeRole::Sample
                    } else {
                        BarcodeRole::Cell
                    };
                    BarcodeDesc {
                        tag_name: format!("b{}", i),
                        role,
                        len: len as u16,
                    }
                })
                .collect()
        }
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
    let sr1 = parts_to_slices(&parts_r1);
    let sr2 = parts_to_slices(&parts_r2);

    // Determine if this is a multi-barcode geometry
    let has_multi_bc = !sr1.multi_bc.is_empty() || !sr2.multi_bc.is_empty();

    // When `s` is present alongside `b`, renumber: s -> b0 (sample), b -> b1 (cell)
    // When only numbered barcodes are used, take them as-is.
    let has_sample_tag = parts_r1
        .iter()
        .chain(parts_r2.iter())
        .any(|p| p.tag_type == GeoTagType::SampleBarcode);
    let has_plain_bc = !sr1.bc.is_empty() || !sr2.bc.is_empty();

    // Build multi-barcode extraction plan
    let (multi_bc_lens, sample_bc_slices_r1, sample_bc_slices_r2) = if has_multi_bc {
        // Collect all levels from both reads
        let mut all_levels: SmallVec<[u8; 4]> = SmallVec::new();
        for (level, _) in sr1.multi_bc.iter().chain(sr2.multi_bc.iter()) {
            if !all_levels.contains(level) {
                all_levels.push(*level);
            }
        }

        // If `s` + `b` pattern: add the plain `b` as the next level
        if has_sample_tag && has_plain_bc {
            let next_level = all_levels.iter().max().map(|m| m + 1).unwrap_or(1);
            all_levels.push(next_level);
        }

        all_levels.sort();
        let num_levels = all_levels.len();

        let mut slices_r1: SmallVec<[SmallVec<[StrSlice; 4]>; 2]> = SmallVec::new();
        let mut slices_r2: SmallVec<[SmallVec<[StrSlice; 4]>; 2]> = SmallVec::new();
        let mut lens: SmallVec<[usize; 4]> = SmallVec::new();

        for &level in &all_levels {
            // Get slices for this level from R1
            let r1_slices: SmallVec<[StrSlice; 4]> = sr1
                .multi_bc
                .iter()
                .filter(|(l, _)| *l == level)
                .flat_map(|(_, s)| s.iter().copied())
                .collect();
            // Get slices for this level from R2
            let r2_slices: SmallVec<[StrSlice; 4]> = sr2
                .multi_bc
                .iter()
                .filter(|(l, _)| *l == level)
                .flat_map(|(_, s)| s.iter().copied())
                .collect();

            // If this is the plain-bc-as-next-level case
            if has_sample_tag && has_plain_bc && level == all_levels[num_levels - 1] {
                // The plain `b` tag becomes this level
                let combined_r1: SmallVec<[StrSlice; 4]> =
                    r1_slices.iter().chain(sr1.bc.iter()).copied().collect();
                let combined_r2: SmallVec<[StrSlice; 4]> =
                    r2_slices.iter().chain(sr2.bc.iter()).copied().collect();
                let len = sum_bounded_len(&combined_r1) + sum_bounded_len(&combined_r2);
                slices_r1.push(combined_r1);
                slices_r2.push(combined_r2);
                lens.push(len);
            } else {
                let len = sum_bounded_len(&r1_slices) + sum_bounded_len(&r2_slices);
                slices_r1.push(r1_slices);
                slices_r2.push(r2_slices);
                lens.push(len);
            }
        }

        (lens, slices_r1, slices_r2)
    } else {
        (SmallVec::new(), SmallVec::new(), SmallVec::new())
    };

    // Compute total lengths for standard barcode
    let bc_len: usize = if has_multi_bc && has_sample_tag && has_plain_bc {
        // The plain `b` length is captured in multi_bc_lens (last level)
        *multi_bc_lens.last().unwrap_or(&0)
    } else {
        sum_bounded_len(&sr1.bc) + sum_bounded_len(&sr2.bc)
    };
    let umi_len: usize = sum_bounded_len(&sr1.umi) + sum_bounded_len(&sr2.umi);

    // Validation
    if bc_len == 0 && multi_bc_lens.is_empty() {
        bail!("custom geometry must include at least one barcode segment (b[N], s[N], or bN[L])");
    }
    if umi_len == 0 {
        bail!("custom geometry must include at least one UMI segment (u[N])");
    }
    for (i, &len) in multi_bc_lens.iter().enumerate() {
        if len > 32 {
            bail!("barcode level {} length {} exceeds maximum of 32", i, len);
        }
    }
    if bc_len > 32 {
        bail!("total barcode length {} exceeds maximum of 32", bc_len);
    }
    if umi_len > 32 {
        bail!("total UMI length {} exceeds maximum of 32", umi_len);
    }

    let has_read_r1 = !sr1.read.is_empty();
    let has_read_r2 = !sr2.read.is_empty();
    if !has_read_r1 && !has_read_r2 {
        bail!("custom geometry must include at least one biological read segment (r[N] or r:)");
    }

    let is_paired_bio = has_read_r1 && has_read_r2;

    Ok(CustomProtocol {
        bc_slices_r1: sr1.bc,
        umi_slices_r1: sr1.umi,
        read_slices_r1: sr1.read,
        bc_slices_r2: sr2.bc,
        umi_slices_r2: sr2.umi,
        read_slices_r2: sr2.read,
        sample_bc_slices_r1,
        sample_bc_slices_r2,
        multi_bc_lens,
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

        let (tag_type, extra_consumed) = match tag_char {
            b's' => (GeoTagType::SampleBarcode, 0),
            b'b' => {
                // Check for numbered barcode: b0[N], b1[N], etc.
                if i + 2 < bytes.len() && bytes[i + 1].is_ascii_digit() {
                    let level = bytes[i + 1] - b'0';
                    // delim is now the digit; the actual delimiter is the next char
                    (GeoTagType::NumberedBarcode(level), 1)
                } else {
                    (GeoTagType::Barcode, 0)
                }
            }
            b'u' => (GeoTagType::Umi, 0),
            b'r' => (GeoTagType::Read, 0),
            b'f' => (GeoTagType::Fixed, 0),
            b'x' => (GeoTagType::Discard, 0),
            _ => bail!(
                "unknown geometry tag '{}' at position {}",
                tag_char as char,
                i
            ),
        };
        // Advance past any extra consumed characters (e.g., the digit in b0)
        let i_delim = i + 1 + extra_consumed;
        if i_delim >= bytes.len() {
            bail!("unexpected end of geometry at position {}", i_delim);
        }
        let delim = bytes[i_delim];

        match delim {
            b'[' => {
                // Bounded: tag[value]
                let value_start = i_delim + 1;
                let close = body[value_start..]
                    .find(']')
                    .ok_or_else(|| anyhow::anyhow!("unmatched '[' at position {}", i_delim))?
                    + value_start;
                let inner = &body[value_start..close];

                let len = if tag_type == GeoTagType::Fixed {
                    // Fixed: value is a DNA sequence, length = sequence length
                    if inner.is_empty() {
                        bail!("empty fixed sequence at position {}", i);
                    }
                    if !inner
                        .bytes()
                        .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
                    {
                        bail!("fixed sequence must contain only ACGT: '{}'", inner);
                    }
                    inner.len() as i32
                } else {
                    // Others: value is a length
                    inner.parse::<i32>().map_err(|_| {
                        anyhow::anyhow!("invalid length '{}' at position {}", inner, value_start)
                    })?
                };

                if len <= 0 && tag_type != GeoTagType::Fixed {
                    bail!(
                        "length must be positive, got {} at position {}",
                        len,
                        value_start
                    );
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
                parts.push(GeoPart { tag_type, len: -1 });
                i = i_delim + 1;
            }
            _ => bail!(
                "expected '[' or ':' after tag '{}', got '{}' at position {}",
                tag_char as char,
                delim as char,
                i_delim,
            ),
        }
    }

    Ok(parts)
}

/// Result of converting geometry parts to slices.
/// Contains standard barcode slices, UMI, read, plus per-level multi-barcode slices.
type SliceVec = SmallVec<[StrSlice; 4]>;

struct SliceResult {
    bc: SliceVec,
    umi: SliceVec,
    read: SliceVec,
    /// Per-level multi-barcode slices (for `s` and `bN` tags).
    /// Maps level -> slices for that level on this read.
    multi_bc: SmallVec<[(u8, SliceVec); 4]>,
}

fn parts_to_slices(parts: &[GeoPart]) -> SliceResult {
    let mut bc_slices = SmallVec::new();
    let mut umi_slices = SmallVec::new();
    let mut read_slices = SmallVec::new();
    let mut multi_bc: SmallVec<[(u8, SliceVec); 4]> = SmallVec::new();

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
            GeoTagType::SampleBarcode => {
                // `s` is syntactic sugar for numbered barcode at level 0
                let entry = multi_bc.iter_mut().find(|(l, _)| *l == 0);
                if let Some((_, slices)) = entry {
                    slices.push(slice);
                } else {
                    multi_bc.push((0, smallvec![slice]));
                }
            }
            GeoTagType::NumberedBarcode(level) => {
                let entry = multi_bc.iter_mut().find(|(l, _)| *l == level);
                if let Some((_, slices)) = entry {
                    slices.push(slice);
                } else {
                    multi_bc.push((level, smallvec![slice]));
                }
            }
            GeoTagType::Fixed | GeoTagType::Discard => {}
        }

        if part.len > 0 {
            offset += part.len as usize;
        }
        // For unbounded (len=-1), offset stays where it is — the slice
        // will consume from offset to end of read.
    }

    SliceResult {
        bc: bc_slices,
        umi: umi_slices,
        read: read_slices,
        multi_bc,
    }
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
        assert_eq!(tech.barcode().unwrap(), b"ACGTACGTACGTACGT");
        assert_eq!(tech.umi.unwrap(), b"AAAAAAAAAAAA");

        // 3' protocol: single bio read from R2
        let reads = proto.extract_mappable_reads(r1, r2);
        match reads {
            AlignableReads::Single(bio) => assert_eq!(bio, b"TGCATGCATGCA"),
            _ => panic!("SE custom protocol should return Single"),
        }
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
    fn test_custom_split_bc_extract() {
        // Geometry: R1 has 8bp BC + 12bp UMI + bio read,
        //           R2 has 8bp BC + bio read.
        let proto = parse_custom_geometry("1{b[8]u[12]r:}2{b[8]r:}").unwrap();

        let r1 = b"AAAACCCCGGGGTTTTAAAAMAPPABLE_R1";
        //         |--BC1--|---UMI----|--bio--->
        //         8       12              ...
        let r2 = b"TTTTGGGGBIOLOGICAL_R2";
        //         |--BC2--|--bio-------->
        //         8            ...

        // extract_tech_seqs returns the first BC slice only (R1's 8bp),
        // not the concatenated 16bp split barcode.
        let tech = proto.extract_tech_seqs(r1, r2);
        assert_eq!(tech.barcode().unwrap(), b"AAAACCCC");
        assert_eq!(tech.umi.unwrap(), b"GGGGTTTTAAAA");

        // extract_mappable_reads returns paired bio segments from both reads
        let reads = proto.extract_mappable_reads(r1, r2);
        match reads {
            AlignableReads::Paired { read1, read2 } => {
                assert_eq!(read1, b"MAPPABLE_R1");
                assert_eq!(read2, b"BIOLOGICAL_R2");
            }
            _ => panic!("split BC protocol should return Paired"),
        }
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
        assert_eq!(tech.barcode().unwrap().len(), 16);
        assert_eq!(tech.umi.unwrap().len(), 12);

        let reads = proto.extract_mappable_reads(r1, r2);
        match reads {
            AlignableReads::Paired { read1, read2 } => {
                assert_eq!(read1, b"MAPPABLE_BIO");
                assert_eq!(read2, b"SECOND_READ_BIO");
            }
            _ => panic!("5' custom protocol should return Paired"),
        }
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
        assert_eq!(tech.barcode().unwrap(), b"ACGTACGTACGTACGT");
        assert_eq!(tech.umi.unwrap(), b"BBBBBBBBBB");

        let reads = proto.extract_mappable_reads(r1, r2);
        match reads {
            AlignableReads::Paired { read1, read2 } => {
                assert_eq!(read1, b"REST_OF_R1");
                assert_eq!(read2, b"BIOLOGICAL_READ_2");
            }
            _ => panic!("PE custom protocol should return Paired"),
        }
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
        match reads {
            AlignableReads::Paired { read1, read2 } => {
                // R1: 16 BC + 12 UMI + 4 fixed = 32, rest is bio
                assert_eq!(read1, b"_BIO_READ_1");
                assert_eq!(read2, b"BIO_READ_2");
            }
            _ => panic!("PE fixed sequence protocol should return Paired"),
        }
    }

    #[test]
    fn test_parse_sample_barcode_sugar() {
        // s[8] = sample barcode on R1, b[16] = cell barcode
        let proto = parse_custom_geometry("1{s[8]b[16]u[12]x:}2{r:}").unwrap();
        assert!(proto.is_multi_barcode());
        assert_eq!(proto.num_barcodes(), 2);

        let descs = proto.barcode_descs();
        assert_eq!(descs[0].tag_name, "b0");
        assert_eq!(descs[0].role, BarcodeRole::Sample);
        assert_eq!(descs[0].len, 8);
        assert_eq!(descs[1].tag_name, "b1");
        assert_eq!(descs[1].role, BarcodeRole::Cell);
        assert_eq!(descs[1].len, 16);

        // R1 = 8bp sample + 16bp cell + 12bp UMI + extra
        let r1 = b"SSSSSSSSACGTACGTACGTACGTAAAAAAAAAAAA_extra";
        let r2 = b"BIO_READ";

        let tech = proto.extract_tech_seqs(r1, r2);
        assert_eq!(tech.barcodes.len(), 2);
        assert_eq!(tech.barcodes[0].unwrap(), b"SSSSSSSS");
        assert_eq!(tech.barcodes[1].unwrap(), b"ACGTACGTACGTACGT");
        assert_eq!(tech.umi.unwrap(), b"AAAAAAAAAAAA");
    }

    #[test]
    fn test_parse_sample_barcode_on_r2() {
        // Flex-like: cell BC on R1, sample BC on R2 after probe sequence
        let proto = parse_custom_geometry("1{b[16]u[12]x:}2{r[25]s[8]}").unwrap();
        assert!(proto.is_multi_barcode());
        assert_eq!(proto.num_barcodes(), 2);

        let descs = proto.barcode_descs();
        assert_eq!(descs[0].role, BarcodeRole::Sample);
        assert_eq!(descs[0].len, 8);
        assert_eq!(descs[1].role, BarcodeRole::Cell);
        assert_eq!(descs[1].len, 16);
    }

    #[test]
    fn test_parse_numbered_barcodes() {
        // Numbered barcodes: b0[8]b1[16]
        let proto = parse_custom_geometry("1{b0[8]b1[16]u[12]x:}2{r:}").unwrap();
        assert!(proto.is_multi_barcode());
        assert_eq!(proto.num_barcodes(), 2);

        let descs = proto.barcode_descs();
        assert_eq!(descs[0].tag_name, "b0");
        assert_eq!(descs[0].len, 8);
        assert_eq!(descs[1].tag_name, "b1");
        assert_eq!(descs[1].len, 16);

        let r1 = b"SSSSSSSSACGTACGTACGTACGTAAAAAAAAAAAA_extra";
        let r2 = b"BIO_READ";

        let tech = proto.extract_tech_seqs(r1, r2);
        assert_eq!(tech.barcodes.len(), 2);
        assert_eq!(tech.barcodes[0].unwrap(), b"SSSSSSSS");
        assert_eq!(tech.barcodes[1].unwrap(), b"ACGTACGTACGTACGT");
    }

    #[test]
    fn test_single_bc_unchanged() {
        // Standard single-barcode: should NOT be multi-barcode
        let proto = parse_custom_geometry("1{b[16]u[12]x:}2{r:}").unwrap();
        assert!(!proto.is_multi_barcode());
        assert_eq!(proto.num_barcodes(), 1);
    }
}
