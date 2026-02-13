//! scATAC-seq protocol.
//!
//! scATAC reads are always paired-end with a separate barcode file.
//! R1 and R2 are biological reads; the barcode is extracted from a
//! separate FASTQ file.
//!
//! Port of C++ `pesc_sc_atac.cpp`.

use super::{AlignableReads, Protocol, TechSeqs};

// ---------------------------------------------------------------------------
// ScatacProtocol
// ---------------------------------------------------------------------------

/// scATAC-seq protocol definition.
#[derive(Debug, Clone)]
pub struct ScatacProtocol {
    /// Barcode length in bases.
    pub bc_len: usize,
    /// Whether to apply Tn5 transposase shift (+4/-9).
    pub tn5_shift: bool,
}

impl ScatacProtocol {
    /// Create with default settings.
    pub fn new(bc_len: usize) -> Self {
        Self {
            bc_len,
            tn5_shift: true,
        }
    }
}

impl Protocol for ScatacProtocol {
    fn name(&self) -> &str {
        "scatac"
    }

    fn is_bio_paired_end(&self) -> bool {
        true // ATAC always maps both R1 and R2 as biological reads
    }

    fn extract_tech_seqs<'a>(&self, _r1: &'a [u8], _r2: &'a [u8]) -> TechSeqs<'a> {
        // Barcodes come from a separate file in scATAC, not from R1/R2.
        // The CLI handles barcode extraction from the third file.
        TechSeqs {
            barcode: None,
            umi: None,
        }
    }

    fn extract_mappable_reads<'a>(&self, r1: &'a [u8], r2: &'a [u8]) -> AlignableReads<'a> {
        // Both R1 and R2 are biological reads
        AlignableReads {
            seq1: if r1.is_empty() { None } else { Some(r1) },
            seq2: if r2.is_empty() { None } else { Some(r2) },
        }
    }

    fn barcode_len(&self) -> usize {
        self.bc_len
    }

    fn umi_len(&self) -> usize {
        0 // scATAC has no UMI
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scatac_protocol_basic() {
        let proto = ScatacProtocol::new(16);
        assert_eq!(proto.name(), "scatac");
        assert!(proto.is_bio_paired_end());
        assert_eq!(proto.barcode_len(), 16);
        assert_eq!(proto.umi_len(), 0);
        assert!(proto.tn5_shift);
    }

    #[test]
    fn test_scatac_extract_reads() {
        let proto = ScatacProtocol::new(16);

        let r1 = b"ACGTACGTACGT";
        let r2 = b"TGCATGCATGCA";

        // Both reads are biological
        let reads = proto.extract_mappable_reads(r1, r2);
        assert_eq!(reads.seq1.unwrap(), b"ACGTACGTACGT");
        assert_eq!(reads.seq2.unwrap(), b"TGCATGCATGCA");

        // Tech seqs come from separate barcode file, not R1/R2
        let tech = proto.extract_tech_seqs(r1, r2);
        assert!(tech.barcode.is_none());
        assert!(tech.umi.is_none());
    }
}
