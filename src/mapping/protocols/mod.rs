//! Protocol abstractions for read mapping.
//!
//! Different sequencing protocols (scRNA-seq, bulk RNA-seq, scATAC-seq) share
//! the same mapping kernel but differ in how they extract barcodes, UMIs, and
//! biological sequences from raw reads. The `Protocol` trait captures this
//! variation so protocol-specific logic can be plugged into the generic
//! mapping pipeline.

use smallvec::SmallVec;

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
// AlignableReads
// ---------------------------------------------------------------------------

/// Biological sequences extracted from a read (pair) for mapping.
///
/// For single-end protocols, only `seq1` is `Some`.
/// For paired-end protocols, both `seq1` and `seq2` are `Some`.
pub struct AlignableReads<'a> {
    pub seq1: Option<&'a [u8]>,
    pub seq2: Option<&'a [u8]>,
}

// ---------------------------------------------------------------------------
// TechSeqs
// ---------------------------------------------------------------------------

/// Technical sequences (barcodes, UMI) extracted from a read pair.
///
/// For single-barcode protocols, `barcodes` has exactly one entry.
/// For multi-barcode protocols (e.g., 10x Flex), `barcodes` has N entries
/// matching the order of `Protocol::barcode_descs()`.
/// For bulk protocols, `barcodes` is empty and `umi` is `None`.
pub struct TechSeqs<'a> {
    pub barcodes: SmallVec<[Option<&'a [u8]>; 2]>,
    pub umi: Option<&'a [u8]>,
}

impl<'a> TechSeqs<'a> {
    /// Convenience: get the primary (first) barcode. For single-barcode
    /// protocols this is the cell barcode. For backward compatibility.
    pub fn barcode(&self) -> Option<&'a [u8]> {
        self.barcodes.first().copied().flatten()
    }
}

// ---------------------------------------------------------------------------
// Protocol trait
// ---------------------------------------------------------------------------

/// Trait for sequencing protocol definitions.
///
/// Each protocol knows how to extract technical (barcode/UMI) and biological
/// (mappable) sequences from raw read pairs. Implementations are expected to
/// be lightweight and `Send + Sync` for use across worker threads.
pub trait Protocol: Send + Sync {
    /// Human-readable protocol name (e.g., "chromium_v3", "bulk").
    fn name(&self) -> &str;

    /// Whether the biological read is paired-end for mapping purposes.
    fn is_bio_paired_end(&self) -> bool;

    /// Extract technical sequences (barcodes, UMI) from raw reads.
    ///
    /// `r1` is always present. `r2` may be empty for single-end data.
    fn extract_tech_seqs<'a>(&self, r1: &'a [u8], r2: &'a [u8]) -> TechSeqs<'a>;

    /// Extract the biological sequence(s) to be mapped.
    ///
    /// `r1` is always present. `r2` may be empty for single-end data.
    fn extract_mappable_reads<'a>(&self, r1: &'a [u8], r2: &'a [u8]) -> AlignableReads<'a>;

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
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_alignable_reads_single_end() {
        let seq = b"ACGTACGT";
        let ar = AlignableReads {
            seq1: Some(seq),
            seq2: None,
        };
        assert!(ar.seq1.is_some());
        assert!(ar.seq2.is_none());
    }

    #[test]
    fn test_tech_seqs_none() {
        let ts = TechSeqs {
            barcodes: smallvec::smallvec![],
            umi: None,
        };
        assert!(ts.barcode().is_none());
        assert!(ts.umi.is_none());
    }
}
