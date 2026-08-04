//! Bulk RNA-seq protocol.
//!
//! Simple protocol that maps all reads directly (no barcode/UMI extraction).

use super::{ExtractedSeqs, Protocol};

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

    fn bio_read_file(&self) -> usize {
        // Bulk reads are biological from the first file onward.
        0
    }

    fn is_bio_paired_end(&self) -> bool {
        self.is_paired
    }

    fn extract<'a>(&self, r1: &'a [u8], r2: &'a [u8]) -> ExtractedSeqs<'a> {
        let mut reads = smallvec::SmallVec::new();
        if self.is_paired {
            reads.push(r1);
            reads.push(r2);
        } else {
            reads.push(r1);
        }
        ExtractedSeqs {
            barcodes: smallvec::smallvec![],
            umi: None,
            reads,
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
        let seqs = proto.extract(b"ACGT", b"");
        assert!(seqs.barcodes.is_empty());
        assert!(seqs.umi.is_none());
        assert_eq!(seqs.reads.len(), 1);
        assert_eq!(seqs.reads[0], b"ACGT");
    }

    #[test]
    fn test_bulk_paired_end() {
        let proto = BulkProtocol::new(true);
        let seqs = proto.extract(b"ACGT", b"TGCA");
        assert_eq!(seqs.reads.len(), 2);
        assert_eq!(seqs.reads[0], b"ACGT");
        assert_eq!(seqs.reads[1], b"TGCA");
    }
}
