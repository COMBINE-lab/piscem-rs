//! scATAC-seq protocol.
//!
//! scATAC reads are always paired-end with a separate barcode file.
//! R1 and R2 are biological reads; the barcode is extracted from a
//! separate FASTQ file.

use super::{ExtractedSeqs, Protocol};

/// scATAC-seq protocol definition.
#[derive(Debug, Clone)]
pub struct ScatacProtocol {
    pub bc_len: usize,
    pub tn5_shift: bool,
}

impl ScatacProtocol {
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
        true
    }

    fn extract<'a>(&self, r1: &'a [u8], r2: &'a [u8]) -> ExtractedSeqs<'a> {
        // Barcodes come from a separate file in scATAC, not from R1/R2.
        ExtractedSeqs {
            barcodes: smallvec::smallvec![],
            umi: None,
            reads: smallvec::smallvec![r1, r2],
        }
    }

    fn barcode_len(&self) -> usize {
        self.bc_len
    }

    fn umi_len(&self) -> usize {
        0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scatac_extract() {
        let proto = ScatacProtocol::new(16);
        let seqs = proto.extract(b"ACGTACGTACGT", b"TGCATGCATGCA");
        assert!(seqs.barcodes.is_empty());
        assert!(seqs.umi.is_none());
        assert_eq!(seqs.reads.len(), 2);
        assert_eq!(seqs.reads[0], b"ACGTACGTACGT");
        assert_eq!(seqs.reads[1], b"TGCATGCATGCA");
    }
}
