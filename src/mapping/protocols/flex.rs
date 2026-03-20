//! 10x Chromium Flex (Fixed RNA Profiling) protocol.
//!
//! Flex uses probe-based capture with a multiplexing/sample barcode.
//! Read structure (native 10x Flex):
//!   R1: [16bp cell BC][12bp UMI][probe capture sequence...]
//!   R2: [probe sequence (biological read)][sample/mux BC at fixed offset]
//!
//! The sample barcode position on R2 is configurable.

use smallvec::smallvec;

use super::{BarcodeDesc, BarcodeRole, ExtractedSeqs, Protocol};

// ---------------------------------------------------------------------------
// ChromiumFlexProtocol
// ---------------------------------------------------------------------------

/// 10x Chromium Flex protocol with bespoke optimized extraction.
///
/// This protocol extracts:
/// - Cell barcode + UMI from R1 (same as standard Chromium)
/// - Sample/multiplexing barcode from R2 at a configurable offset
/// - Biological read: R2 probe sequence (up to the sample BC)
#[derive(Debug, Clone)]
pub struct ChromiumFlexProtocol {
    cell_bc_len: usize,
    umi_len: usize,
    sample_bc_len: usize,
    /// Byte offset on R2 where the sample barcode starts
    sample_bc_offset: usize,
}

impl ChromiumFlexProtocol {
    /// Create a Flex protocol with default Chromium cell barcode (16bp) and UMI (12bp).
    ///
    /// `sample_bc_len`: length of the sample/multiplexing barcode (e.g., 8 for 16-plex)
    /// `sample_bc_offset`: position on R2 where the sample BC starts (after probe sequence)
    pub fn new(sample_bc_len: usize, sample_bc_offset: usize) -> Self {
        Self {
            cell_bc_len: 16,
            umi_len: 12,
            sample_bc_len,
            sample_bc_offset,
        }
    }

    /// Create from a named variant.
    ///
    /// Recognized names: "chromium_flex" (defaults to 8bp sample BC at offset 25).
    pub fn from_name(name: &str) -> Option<Self> {
        match name.to_lowercase().as_str() {
            // Default 16-plex: 8bp sample BC. Offset is protocol-dependent;
            // 25 is a placeholder — users should specify via geometry grammar
            // for exact positioning.
            "chromium_flex" => Some(Self::new(8, 25)),
            _ => None,
        }
    }

    /// Sample barcode length in bases.
    pub fn sample_bc_len(&self) -> usize {
        self.sample_bc_len
    }
}

impl Protocol for ChromiumFlexProtocol {
    fn name(&self) -> &str {
        "chromium_flex"
    }

    fn is_bio_paired_end(&self) -> bool {
        false // Only the probe sequence on R2 is mapped
    }

    fn extract<'a>(&self, r1: &'a [u8], r2: &'a [u8]) -> ExtractedSeqs<'a> {
        let cell_bc = if r1.len() >= self.cell_bc_len {
            Some(&r1[..self.cell_bc_len])
        } else {
            None
        };
        let umi = if r1.len() >= self.cell_bc_len + self.umi_len {
            Some(&r1[self.cell_bc_len..self.cell_bc_len + self.umi_len])
        } else {
            None
        };
        let sample_bc = if r2.len() >= self.sample_bc_offset + self.sample_bc_len {
            Some(&r2[self.sample_bc_offset..self.sample_bc_offset + self.sample_bc_len])
        } else {
            None
        };
        let end = self.sample_bc_offset.min(r2.len());
        ExtractedSeqs {
            barcodes: smallvec![sample_bc, cell_bc],
            umi,
            reads: smallvec![&r2[..end]],
        }
    }

    fn barcode_len(&self) -> usize {
        self.cell_bc_len
    }

    fn umi_len(&self) -> usize {
        self.umi_len
    }

    fn num_barcodes(&self) -> usize {
        2
    }

    fn barcode_descs(&self) -> Vec<BarcodeDesc> {
        vec![
            BarcodeDesc {
                tag_name: "b0".to_string(),
                role: BarcodeRole::Sample,
                len: self.sample_bc_len as u16,
            },
            BarcodeDesc {
                tag_name: "b1".to_string(),
                role: BarcodeRole::Cell,
                len: self.cell_bc_len as u16,
            },
        ]
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_flex_basic() {
        let proto = ChromiumFlexProtocol::new(8, 25);
        assert_eq!(proto.name(), "chromium_flex");
        assert_eq!(proto.barcode_len(), 16);
        assert_eq!(proto.umi_len(), 12);
        assert_eq!(proto.num_barcodes(), 2);
        assert!(proto.is_multi_barcode());
        assert!(!proto.is_bio_paired_end());
    }

    #[test]
    fn test_flex_extract() {
        let proto = ChromiumFlexProtocol::new(8, 25);

        let r1 = b"ACGTACGTACGTACGTAAAAAAAAAAAA_extra";
        let r2 = b"PROBE_SEQUENCE_25_BASES__SSSSSSSS";

        let seqs = proto.extract(r1, r2);
        assert_eq!(seqs.barcodes[0].unwrap(), b"SSSSSSSS"); // sample BC
        assert_eq!(seqs.barcodes[1].unwrap(), b"ACGTACGTACGTACGT"); // cell BC
        assert_eq!(seqs.umi.unwrap(), b"AAAAAAAAAAAA");
        assert_eq!(seqs.reads.len(), 1);
        assert_eq!(seqs.reads[0], b"PROBE_SEQUENCE_25_BASES__");
    }

    #[test]
    fn test_flex_barcode_descs() {
        let proto = ChromiumFlexProtocol::new(8, 25);
        let descs = proto.barcode_descs();
        assert_eq!(descs.len(), 2);
        assert_eq!(descs[0].tag_name, "b0");
        assert_eq!(descs[0].role, BarcodeRole::Sample);
        assert_eq!(descs[0].len, 8);
        assert_eq!(descs[1].tag_name, "b1");
        assert_eq!(descs[1].role, BarcodeRole::Cell);
        assert_eq!(descs[1].len, 16);
    }

    #[test]
    fn test_flex_from_name() {
        assert!(ChromiumFlexProtocol::from_name("chromium_flex").is_some());
        assert!(ChromiumFlexProtocol::from_name("CHROMIUM_FLEX").is_some());
        assert!(ChromiumFlexProtocol::from_name("unknown").is_none());
    }

    #[test]
    fn test_flex_short_reads() {
        let proto = ChromiumFlexProtocol::new(8, 25);

        // R2 too short for sample BC
        let r1 = b"ACGTACGTACGTACGTAAAAAAAAAAAA";
        let r2 = b"SHORT";

        let seqs = proto.extract(r1, r2);
        assert!(seqs.barcodes[0].is_none()); // sample BC not extractable
        assert!(seqs.barcodes[1].is_some()); // cell BC still ok
    }
}
