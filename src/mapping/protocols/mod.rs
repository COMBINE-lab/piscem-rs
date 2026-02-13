//! Protocol abstractions for read mapping.
//!
//! Different sequencing protocols (scRNA-seq, bulk RNA-seq, scATAC-seq) share
//! the same mapping kernel but differ in how they extract barcodes, UMIs, and
//! biological sequences from raw reads. The `Protocol` trait captures this
//! variation so protocol-specific logic can be plugged into the generic
//! mapping pipeline.

pub mod bulk;
pub mod scatac;
pub mod scrna;

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

/// Technical sequences (barcode, UMI) extracted from a read pair.
///
/// For bulk protocols, both fields are `None`.
pub struct TechSeqs<'a> {
    pub barcode: Option<&'a [u8]>,
    pub umi: Option<&'a [u8]>,
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

    /// Extract technical sequences (barcode, UMI) from raw reads.
    ///
    /// `r1` is always present. `r2` may be empty for single-end data.
    fn extract_tech_seqs<'a>(&self, r1: &'a [u8], r2: &'a [u8]) -> TechSeqs<'a>;

    /// Extract the biological sequence(s) to be mapped.
    ///
    /// `r1` is always present. `r2` may be empty for single-end data.
    fn extract_mappable_reads<'a>(&self, r1: &'a [u8], r2: &'a [u8]) -> AlignableReads<'a>;

    /// Expected barcode length in bases (0 for bulk).
    fn barcode_len(&self) -> usize;

    /// Expected UMI length in bases (0 for bulk).
    fn umi_len(&self) -> usize;
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
            barcode: None,
            umi: None,
        };
        assert!(ts.barcode.is_none());
        assert!(ts.umi.is_none());
    }
}
