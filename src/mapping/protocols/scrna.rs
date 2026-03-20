//! 10x Chromium scRNA-seq protocols.
//!
//! Port of C++ `piscem-cpp/include/sc/util.hpp` protocol classes.

use smallvec::smallvec;

use super::{ExtractedSeqs, Protocol};

// ---------------------------------------------------------------------------
// ChromiumVersion
// ---------------------------------------------------------------------------

/// Supported 10x Chromium chemistry versions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChromiumVersion {
    V2,
    V2_5p,
    V3,
    V3_5p,
    V4_3p,
}

// ---------------------------------------------------------------------------
// ChromiumProtocol
// ---------------------------------------------------------------------------

/// 10x Chromium scRNA-seq protocol definition.
#[derive(Debug, Clone)]
pub struct ChromiumProtocol {
    version: ChromiumVersion,
    bc_len: usize,
    umi_len: usize,
    tso_len: usize,
}

impl ChromiumProtocol {
    /// Create from a geometry name string (case-insensitive).
    ///
    /// Recognized names: "chromium_v2", "chromium_v2_5p", "chromium_v3",
    /// "chromium_v3_5p", "chromium_v4_3p".
    pub fn from_name(name: &str) -> Option<Self> {
        let name_lower = name.to_lowercase();
        let version = match name_lower.as_str() {
            "chromium_v2" => ChromiumVersion::V2,
            "chromium_v2_5p" => ChromiumVersion::V2_5p,
            "chromium_v3" => ChromiumVersion::V3,
            "chromium_v3_5p" => ChromiumVersion::V3_5p,
            "chromium_v4_3p" => ChromiumVersion::V4_3p,
            _ => return None,
        };
        Some(Self::new(version))
    }

    /// Create from an explicit version.
    pub fn new(version: ChromiumVersion) -> Self {
        let (bc_len, umi_len, tso_len) = match version {
            ChromiumVersion::V2 => (16, 10, 0),
            ChromiumVersion::V2_5p => (16, 10, 13),
            ChromiumVersion::V3 => (16, 12, 0),
            ChromiumVersion::V3_5p => (16, 12, 13),
            ChromiumVersion::V4_3p => (16, 12, 0),
        };
        Self {
            version,
            bc_len,
            umi_len,
            tso_len,
        }
    }

    /// The chromium chemistry version.
    pub fn version(&self) -> ChromiumVersion {
        self.version
    }

    /// Whether this is a 5' protocol (has TSO sequence).
    pub fn is_5prime(&self) -> bool {
        self.tso_len > 0
    }

    /// Total technical prefix length in R1 (bc + umi + tso).
    pub fn tech_prefix_len(&self) -> usize {
        self.bc_len + self.umi_len + self.tso_len
    }
}

impl Protocol for ChromiumProtocol {
    fn name(&self) -> &str {
        match self.version {
            ChromiumVersion::V2 => "chromium_v2",
            ChromiumVersion::V2_5p => "chromium_v2_5p",
            ChromiumVersion::V3 => "chromium_v3",
            ChromiumVersion::V3_5p => "chromium_v3_5p",
            ChromiumVersion::V4_3p => "chromium_v4_3p",
        }
    }

    fn is_bio_paired_end(&self) -> bool {
        // 5' protocols produce paired biological reads
        self.is_5prime()
    }

    fn extract<'a>(&self, r1: &'a [u8], r2: &'a [u8]) -> ExtractedSeqs<'a> {
        let barcode = if r1.len() >= self.bc_len {
            Some(&r1[..self.bc_len])
        } else {
            None
        };
        let umi = if r1.len() >= self.bc_len + self.umi_len {
            Some(&r1[self.bc_len..self.bc_len + self.umi_len])
        } else {
            None
        };
        let mut reads = smallvec::SmallVec::new();
        if self.is_5prime() {
            let start = self.tech_prefix_len().min(r1.len());
            reads.push(&r1[start..]);
            reads.push(r2);
        } else {
            reads.push(r2);
        }
        ExtractedSeqs {
            barcodes: smallvec![barcode],
            umi,
            reads,
        }
    }

    fn barcode_len(&self) -> usize {
        self.bc_len
    }

    fn umi_len(&self) -> usize {
        self.umi_len
    }
}

// ---------------------------------------------------------------------------
// Barcode recovery
// ---------------------------------------------------------------------------

/// Attempt to recover a barcode with a single N base.
///
/// If the barcode has exactly 1 N, replace it with A and return
/// `Some(recovered)`. Returns `None` if there are 0 or >1 N bases.
pub fn recover_barcode(bc: &[u8]) -> Option<Vec<u8>> {
    let n_count = bc.iter().filter(|&&b| b == b'N' || b == b'n').count();
    if n_count != 1 {
        return None;
    }
    let mut recovered = bc.to_vec();
    for b in &mut recovered {
        if *b == b'N' || *b == b'n' {
            *b = b'A';
        }
    }
    Some(recovered)
}

/// Check if a barcode contains any N bases.
pub fn barcode_has_n(bc: &[u8]) -> bool {
    bc.iter().any(|&b| b == b'N' || b == b'n')
}

/// Count N bases in a barcode.
pub fn count_ns(bc: &[u8]) -> usize {
    bc.iter().filter(|&&b| b == b'N' || b == b'n').count()
}

/// Check whether a sequence contains only valid ACGT bases (upper or lowercase).
///
/// Matches C++ `fromChars()` validation: returns `false` if any character
/// is not in {A, C, G, T, a, c, g, t}.
#[inline]
pub fn is_all_acgt(seq: &[u8]) -> bool {
    seq.iter()
        .all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_v3_extract() {
        let proto = ChromiumProtocol::new(ChromiumVersion::V3);
        assert_eq!(proto.name(), "chromium_v3");
        assert_eq!(proto.barcode_len(), 16);
        assert_eq!(proto.umi_len(), 12);
        assert!(!proto.is_bio_paired_end());
        assert!(!proto.is_5prime());

        // R1 = 16 BC + 12 UMI + bio, R2 = bio read
        let r1 = b"ACGTACGTACGTACGTAAAAAAAAAAAA_extra_bio";
        let r2 = b"TGCATGCATGCA";

        let seqs = proto.extract(r1, r2);
        assert_eq!(seqs.barcodes[0].unwrap(), b"ACGTACGTACGTACGT");
        assert_eq!(seqs.umi.unwrap(), b"AAAAAAAAAAAA");
        assert_eq!(seqs.reads.len(), 1);
        assert_eq!(seqs.reads[0], b"TGCATGCATGCA");
    }

    #[test]
    fn test_v2_5p_paired() {
        let proto = ChromiumProtocol::new(ChromiumVersion::V2_5p);
        assert_eq!(proto.name(), "chromium_v2_5p");
        assert_eq!(proto.barcode_len(), 16);
        assert_eq!(proto.umi_len(), 10);
        assert!(proto.is_bio_paired_end());
        assert!(proto.is_5prime());

        let r1 = b"ACGTACGTACGTACGTBBBBBBBBBBCCCCCCCCCCCCCMAPPABLE_BIO";
        let r2 = b"SECOND_READ_BIO";

        let seqs = proto.extract(r1, r2);
        assert_eq!(seqs.barcodes[0].unwrap().len(), 16);
        assert_eq!(seqs.umi.unwrap().len(), 10);
        assert_eq!(seqs.reads.len(), 2);
        assert_eq!(seqs.reads[0], b"MAPPABLE_BIO");
        assert_eq!(seqs.reads[1], b"SECOND_READ_BIO");
    }

    #[test]
    fn test_recover_barcode() {
        // No N → None (no recovery needed)
        assert!(recover_barcode(b"ACGTACGTACGTACGT").is_none());

        // One N → recovery
        let result = recover_barcode(b"ACNTACGTACGTACGT").unwrap();
        assert_eq!(result, b"ACATACGTACGTACGT");

        // Two N's → None (can't recover)
        assert!(recover_barcode(b"NCNTACGTACGTACGT").is_none());
    }

    #[test]
    fn test_from_name() {
        assert!(ChromiumProtocol::from_name("chromium_v3").is_some());
        assert!(ChromiumProtocol::from_name("CHROMIUM_V3").is_some());
        assert!(ChromiumProtocol::from_name("chromium_v4_3p").is_some());
        assert!(ChromiumProtocol::from_name("unknown").is_none());

        let p = ChromiumProtocol::from_name("chromium_v2").unwrap();
        assert_eq!(p.version(), ChromiumVersion::V2);
        assert_eq!(p.barcode_len(), 16);
        assert_eq!(p.umi_len(), 10);
    }

    #[test]
    fn test_barcode_has_n() {
        assert!(!barcode_has_n(b"ACGTACGT"));
        assert!(barcode_has_n(b"ACNTACGT"));
        assert!(barcode_has_n(b"ACnTACGT"));
    }

    #[test]
    fn test_v4_3p() {
        let proto = ChromiumProtocol::new(ChromiumVersion::V4_3p);
        assert_eq!(proto.name(), "chromium_v4_3p");
        assert_eq!(proto.barcode_len(), 16);
        assert_eq!(proto.umi_len(), 12);
        assert!(!proto.is_5prime());
    }
}
