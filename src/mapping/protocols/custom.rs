//! Custom read geometry protocol.
//!
//! Wraps the `seq_geom_parser` crate to provide geometry parsing and
//! extraction behind the piscem `Protocol` trait. This module is a thin
//! adapter — all parsing, validation, and extraction logic lives in
//! `seq_geom_parser`.

use anyhow::Result;

use seq_geom_parser::{self as sgp, CompiledGeom};

use super::{BarcodeDesc, BarcodeRole, ExtractedSeqs, Protocol};

// ---------------------------------------------------------------------------
// CustomProtocol
// ---------------------------------------------------------------------------

/// Custom geometry protocol built from a parsed geometry string.
///
/// Wraps a [`CompiledGeom`] from `seq_geom_parser` and adapts it to the
/// piscem [`Protocol`] trait.
#[derive(Debug, Clone)]
pub struct CustomProtocol {
    compiled: CompiledGeom,
}

impl CustomProtocol {
    /// Total primary barcode length in bases.
    pub fn total_bc_len(&self) -> usize {
        let meta = self.compiled.meta();
        if meta.num_bc_levels > 0 {
            // Primary barcode is the last level (cell BC)
            *meta.bc_lens.last().unwrap_or(&0)
        } else {
            0
        }
    }

    /// Total UMI length in bases.
    pub fn total_umi_len(&self) -> usize {
        self.compiled.meta().umi_len
    }
}

impl Protocol for CustomProtocol {
    fn name(&self) -> &str {
        "custom"
    }

    fn is_bio_paired_end(&self) -> bool {
        self.compiled.meta().is_paired_read
    }

    fn extract<'a>(&self, r1: &'a [u8], r2: &'a [u8]) -> ExtractedSeqs<'a> {
        self.compiled.extract(r1, r2)
    }

    fn barcode_len(&self) -> usize {
        self.total_bc_len()
    }

    fn umi_len(&self) -> usize {
        self.total_umi_len()
    }

    fn num_barcodes(&self) -> usize {
        self.compiled.meta().num_bc_levels
    }

    fn barcode_descs(&self) -> Vec<BarcodeDesc> {
        let meta = self.compiled.meta();
        if meta.num_bc_levels <= 1 {
            vec![BarcodeDesc {
                tag_name: "b".to_string(),
                role: BarcodeRole::Cell,
                len: meta.bc_lens.first().copied().unwrap_or(0) as u16,
            }]
        } else {
            meta.bc_lens
                .iter()
                .enumerate()
                .map(|(i, &len)| {
                    let role = if i == 0 && meta.num_bc_levels > 1 {
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

    fn normalization_meta(&self) -> Option<&sgp::NormalizationMeta> {
        Some(&self.compiled.meta().normalization)
    }
}

// ---------------------------------------------------------------------------
// Parsing entry point
// ---------------------------------------------------------------------------

/// Parse a custom geometry string into a [`CustomProtocol`].
///
/// Returns an error with descriptive messages if the geometry is invalid.
///
/// # Examples
///
/// ```ignore
/// let proto = parse_custom_geometry("1{b[16]u[12]x:}2{r:}").unwrap();
/// assert_eq!(proto.barcode_len(), 16);
/// assert_eq!(proto.umi_len(), 12);
/// ```
pub fn parse_custom_geometry(geom: &str) -> Result<CustomProtocol> {
    // Parse with seq_geom_parser
    let fragment_geom = sgp::parse_geometry(geom).map_err(|errors| {
        let msg = sgp::format_errors(geom, &errors);
        anyhow::anyhow!(
            "Could not parse geometry '{}'. Parse errors:\n{}",
            geom,
            msg
        )
    })?;

    // Validate
    sgp::validate_geometry(&fragment_geom)
        .map_err(|e| anyhow::anyhow!("Geometry validation failed for '{}': {}", geom, e))?;

    // Compile extraction plan
    let compiled = CompiledGeom::from_fragment_geom(&fragment_geom)
        .map_err(|e| anyhow::anyhow!("Failed to compile geometry '{}': {}", geom, e))?;

    Ok(CustomProtocol { compiled })
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_chromium_v3() {
        let proto = parse_custom_geometry("1{b[16]u[12]x:}2{r:}").unwrap();
        assert_eq!(proto.barcode_len(), 16);
        assert_eq!(proto.umi_len(), 12);
        assert!(!proto.is_bio_paired_end());
        assert_eq!(proto.num_barcodes(), 1);
    }

    #[test]
    fn parse_chromium_v2() {
        let proto = parse_custom_geometry("1{b[16]u[10]x:}2{r:}").unwrap();
        assert_eq!(proto.barcode_len(), 16);
        assert_eq!(proto.umi_len(), 10);
    }

    #[test]
    fn parse_5prime_with_tso() {
        let proto = parse_custom_geometry("1{b[16]u[12]x[13]r:}2{r:}").unwrap();
        assert_eq!(proto.barcode_len(), 16);
        assert_eq!(proto.umi_len(), 12);
        assert!(proto.is_bio_paired_end());
    }

    #[test]
    fn parse_flex_v1() {
        let proto = parse_custom_geometry("1{b[16]u[12]x:}2{r[50]x[18]s[8]x:}").unwrap();
        assert_eq!(proto.num_barcodes(), 2);
        let descs = proto.barcode_descs();
        assert_eq!(descs[0].role, BarcodeRole::Sample);
        assert_eq!(descs[0].len, 8);
        assert_eq!(descs[1].role, BarcodeRole::Cell);
        assert_eq!(descs[1].len, 16);
    }

    #[test]
    fn parse_flex_v2() {
        let proto =
            parse_custom_geometry("1{b[16]u[12]x[0-3]f[TTGCTAGGACCG]s[10]x:}2{r:}").unwrap();
        assert_eq!(proto.num_barcodes(), 2);
        let descs = proto.barcode_descs();
        assert_eq!(descs[0].len, 10); // sample BC
        assert_eq!(descs[1].len, 16); // cell BC
    }

    #[test]
    fn extract_v3() {
        let proto = parse_custom_geometry("1{b[16]u[12]x:}2{r:}").unwrap();
        let r1 = b"ACGTACGTACGTACGTAAAAAAAAAAAA_extra";
        let r2 = b"BIOLOGICAL_READ";

        let seqs = proto.extract(r1, r2);
        assert_eq!(seqs.barcodes[0], Some(&r1[..16]));
        assert_eq!(seqs.umi, Some(&r1[16..28]));
        assert_eq!(seqs.reads.len(), 1);
        assert_eq!(seqs.reads[0], r2.as_slice());
    }

    #[test]
    fn extract_flex_v1() {
        let proto = parse_custom_geometry("1{b[16]u[12]x:}2{r[50]x[18]s[8]x:}").unwrap();

        let r1 = b"CELLBARCODEACGTUMI_AAAAAA_extra";
        let mut r2 = vec![b'N'; 80];
        r2[68..76].copy_from_slice(b"SAMPLEBC");

        let seqs = proto.extract(r1, &r2);
        assert_eq!(seqs.barcodes.len(), 2);
        assert_eq!(seqs.barcodes[0], Some(&r2[68..76])); // sample (b0)
        assert_eq!(seqs.barcodes[1], Some(&r1[..16])); // cell (b1)
    }

    #[test]
    fn extract_flex_v2_anchor() {
        let proto =
            parse_custom_geometry("1{b[16]u[12]x[0-3]f[TTGCTAGGACCG]s[10]x:}2{r:}").unwrap();

        let bc = b"ACGTACGTACGTACGT";
        let umi = b"AAAAAAAAAAAA";
        let gap = b"NN";
        let anchor = b"TTGCTAGGACCG";
        let sample = b"SAMPLEBC10";
        let mut r1 = Vec::new();
        r1.extend_from_slice(bc);
        r1.extend_from_slice(umi);
        r1.extend_from_slice(gap);
        r1.extend_from_slice(anchor);
        r1.extend_from_slice(sample);
        r1.extend_from_slice(b"extra");

        let r2 = b"BIO_READ";

        let seqs = proto.extract(&r1, r2);
        assert_eq!(seqs.barcodes[0], Some(sample.as_slice())); // sample
        assert_eq!(seqs.barcodes[1], Some(bc.as_slice())); // cell
        assert_eq!(seqs.umi, Some(umi.as_slice()));
    }

    #[test]
    fn extract_boundary_resolved_geometry() {
        let proto = parse_custom_geometry("1{r:f[ACAGT]b[9-11]}2{u[12]x:}").unwrap();

        let read_prefix = b"BIOREAD";
        let anchor = b"ACAGT";
        let barcode = b"BARCODE09";
        let mut r1 = Vec::new();
        r1.extend_from_slice(read_prefix);
        r1.extend_from_slice(anchor);
        r1.extend_from_slice(barcode);

        let r2 = b"TTTTTTTTTTTTtail";
        let seqs = proto.extract(&r1, r2);
        assert_eq!(seqs.reads[0], read_prefix.as_slice());
        assert_eq!(seqs.barcodes[0], Some(barcode.as_slice()));
        assert_eq!(seqs.umi, Some(&r2[..12]));
    }

    #[test]
    fn normalization_meta_for_variable_barcode_geometry() {
        let proto = parse_custom_geometry("1{s[7-8]b[16]u[12]f[ACGT]x:}2{r:}").unwrap();
        let meta = proto.normalization_meta().unwrap();
        assert!(meta.any_normalization);
        assert_eq!(meta.bc_needs_padding.as_slice(), &[true, false]);
        assert!(!meta.umi_needs_padding);
    }

    #[test]
    fn normalization_meta_for_variable_umi_geometry() {
        let proto = parse_custom_geometry("1{b[16]u[10-12]f[ACGT]x:}2{r:}").unwrap();
        let meta = proto.normalization_meta().unwrap();
        assert!(meta.any_normalization);
        assert_eq!(meta.bc_needs_padding.as_slice(), &[false]);
        assert!(meta.umi_needs_padding);
    }

    #[test]
    fn invalid_geometry_error() {
        let result = parse_custom_geometry("1{b[16]u[12]z:}2{r:}");
        assert!(result.is_err());
        let msg = result.unwrap_err().to_string();
        assert!(msg.contains("Could not parse geometry"));
    }

    #[test]
    fn missing_umi_error() {
        let result = parse_custom_geometry("1{b[16]x:}2{r:}");
        assert!(result.is_err());
        let msg = result.unwrap_err().to_string();
        assert!(msg.contains("UMI"));
    }
}
