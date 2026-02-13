//! Bulk RNA-seq protocol.
//!
//! Simple protocol that maps all reads directly (no barcode/UMI extraction).

use super::{AlignableReads, Protocol, TechSeqs};

/// Bulk sequencing protocol.
#[derive(Debug, Clone)]
pub struct BulkProtocol {
    pub is_paired: bool,
}

impl BulkProtocol {
    pub fn new(is_paired: bool) -> Self {
        Self { is_paired }
    }
}

impl Protocol for BulkProtocol {
    fn name(&self) -> &str {
        "bulk"
    }

    fn is_bio_paired_end(&self) -> bool {
        self.is_paired
    }

    fn extract_tech_seqs<'a>(&self, _r1: &'a [u8], _r2: &'a [u8]) -> TechSeqs<'a> {
        TechSeqs {
            barcode: None,
            umi: None,
        }
    }

    fn extract_mappable_reads<'a>(&self, r1: &'a [u8], r2: &'a [u8]) -> AlignableReads<'a> {
        AlignableReads {
            seq1: Some(r1),
            seq2: if self.is_paired { Some(r2) } else { None },
        }
    }

    fn barcode_len(&self) -> usize {
        0
    }

    fn umi_len(&self) -> usize {
        0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bulk_single_end() {
        let proto = BulkProtocol::new(false);
        assert_eq!(proto.name(), "bulk");
        assert!(!proto.is_bio_paired_end());
        assert_eq!(proto.barcode_len(), 0);
        assert_eq!(proto.umi_len(), 0);

        let reads = proto.extract_mappable_reads(b"ACGT", b"");
        assert_eq!(reads.seq1.unwrap(), b"ACGT");
        assert!(reads.seq2.is_none());

        let tech = proto.extract_tech_seqs(b"ACGT", b"");
        assert!(tech.barcode.is_none());
        assert!(tech.umi.is_none());
    }

    #[test]
    fn test_bulk_paired_end() {
        let proto = BulkProtocol::new(true);
        assert!(proto.is_bio_paired_end());

        let reads = proto.extract_mappable_reads(b"ACGT", b"TGCA");
        assert_eq!(reads.seq1.unwrap(), b"ACGT");
        assert_eq!(reads.seq2.unwrap(), b"TGCA");
    }
}
