//! Protocol abstractions for read mapping.
//!
//! Different sequencing protocols (scRNA-seq, bulk RNA-seq, scATAC-seq) share
//! the same mapping kernel but differ in how they extract barcodes, UMIs, and
//! biological sequences from raw reads. The `Protocol` trait captures this
//! variation so protocol-specific logic can be plugged into the generic
//! mapping pipeline.

pub use seq_geom_parser::ExtractedSeqs;

pub mod bulk;
pub mod custom;
pub mod flex;
pub mod scatac;
pub mod scrna;

// ---------------------------------------------------------------------------
// BarcodeRole / BarcodeDesc
// ---------------------------------------------------------------------------

/// The semantic role of a barcode in the experiment hierarchy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BarcodeRole {
    /// Sample/library-level barcode (outermost collation level)
    Sample,
    /// Cell-level barcode
    Cell,
    /// Feature barcode (e.g., CITE-seq antibody tag)
    Feature,
}

/// Descriptor for a single barcode component in a multi-barcode protocol.
#[derive(Debug, Clone)]
pub struct BarcodeDesc {
    /// Short tag name used in RAD header (e.g., "b0", "b1")
    pub tag_name: String,
    /// Semantic role
    pub role: BarcodeRole,
    /// Length in bases
    pub len: u16,
}

// ---------------------------------------------------------------------------
// Protocol trait
// ---------------------------------------------------------------------------

/// Trait for sequencing protocol definitions.
///
/// Each protocol knows how to extract technical (barcode/UMI) and biological
/// (mappable) sequences from raw read pairs into a unified [`ExtractedSeqs`].
///
/// Implementations are expected to be lightweight and `Send + Sync` for use
/// across worker threads.
pub trait Protocol: Send + Sync {
    /// Human-readable protocol name (e.g., "chromium_v3", "bulk").
    fn name(&self) -> &str;

    /// Whether the biological read is paired-end for mapping purposes.
    fn is_bio_paired_end(&self) -> bool;

    /// Extract all sequences (barcodes, UMI, biological reads) from a read pair.
    ///
    /// Returns an [`ExtractedSeqs`] containing:
    /// - `barcodes`: barcode at each level (single-barcode: 1 entry, multi: N entries)
    /// - `umi`: the UMI sequence
    /// - `reads`: biological read(s) for mapping (1 for SE, 2 for PE)
    fn extract<'a>(&self, r1: &'a [u8], r2: &'a [u8]) -> ExtractedSeqs<'a>;

    /// Expected primary (cell) barcode length in bases (0 for bulk).
    fn barcode_len(&self) -> usize;

    /// Expected UMI length in bases (0 for bulk).
    fn umi_len(&self) -> usize;

    /// Number of distinct barcode components. Default: 1.
    fn num_barcodes(&self) -> usize {
        1
    }

    /// Descriptors for each barcode component. Default: single cell barcode.
    fn barcode_descs(&self) -> Vec<BarcodeDesc> {
        vec![BarcodeDesc {
            tag_name: "b".to_string(),
            role: BarcodeRole::Cell,
            len: self.barcode_len() as u16,
        }]
    }

    /// Whether this protocol has multiple barcode components.
    fn is_multi_barcode(&self) -> bool {
        self.num_barcodes() > 1
    }

    /// Convenience: get the primary (first) barcode from an extraction result.
    fn primary_barcode<'a>(&self, seqs: &ExtractedSeqs<'a>) -> Option<&'a [u8]> {
        seqs.barcodes.first().copied().flatten()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extracted_seqs_fields() {
        let seqs = ExtractedSeqs {
            barcodes: smallvec::smallvec![Some(b"ACGT".as_slice())],
            umi: Some(b"TTTT".as_slice()),
            reads: smallvec::smallvec![b"BIO".as_slice()],
        };
        assert_eq!(seqs.barcodes[0], Some(b"ACGT".as_slice()));
        assert_eq!(seqs.umi, Some(b"TTTT".as_slice()));
        assert_eq!(seqs.reads.len(), 1);
    }

    #[test]
    fn test_extracted_seqs_empty() {
        let seqs = ExtractedSeqs {
            barcodes: smallvec::smallvec![],
            umi: None,
            reads: smallvec::smallvec![],
        };
        assert!(seqs.barcodes.is_empty());
        assert!(seqs.umi.is_none());
        assert!(seqs.reads.is_empty());
    }
}
