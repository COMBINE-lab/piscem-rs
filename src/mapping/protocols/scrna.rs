//! 10x Chromium scRNA-seq protocols.
//!
//! Port of C++ `piscem-cpp/include/sc/util.hpp` protocol classes.

use super::{AlignableReads, Protocol, TechSeqs};

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

    fn extract_tech_seqs<'a>(&self, r1: &'a [u8], _r2: &'a [u8]) -> TechSeqs<'a> {
        let bc_end = self.bc_len.min(r1.len());
        let umi_end = (self.bc_len + self.umi_len).min(r1.len());
        TechSeqs {
            barcode: Some(&r1[..bc_end]),
            umi: Some(&r1[bc_end..umi_end]),
        }
    }

    fn extract_mappable_reads<'a>(&self, r1: &'a [u8], r2: &'a [u8]) -> AlignableReads<'a> {
        if self.is_5prime() {
            // 5' protocols: map the remainder of R1 (after BC+UMI+TSO) AND R2
            let start = self.tech_prefix_len().min(r1.len());
            let bio_r1 = &r1[start..];
            AlignableReads {
                seq1: if bio_r1.is_empty() { None } else { Some(bio_r1) },
                seq2: Some(r2),
            }
        } else {
            // 3' protocols: map R2 only
            AlignableReads {
                seq1: Some(r2),
                seq2: None,
            }
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

        let tech = proto.extract_tech_seqs(r1, r2);
        assert_eq!(tech.barcode.unwrap(), b"ACGTACGTACGTACGT");
        assert_eq!(tech.umi.unwrap(), b"AAAAAAAAAAAA");

        // 3' protocol: map R2 only
        let reads = proto.extract_mappable_reads(r1, r2);
        assert_eq!(reads.seq1.unwrap(), b"TGCATGCATGCA");
        assert!(reads.seq2.is_none());
    }

    #[test]
    fn test_v2_5p_paired() {
        let proto = ChromiumProtocol::new(ChromiumVersion::V2_5p);
        assert_eq!(proto.name(), "chromium_v2_5p");
        assert_eq!(proto.barcode_len(), 16);
        assert_eq!(proto.umi_len(), 10);
        assert!(proto.is_bio_paired_end());
        assert!(proto.is_5prime());

        // TSO = 13 bases. Tech prefix = 16 + 10 + 13 = 39 bases.
        let r1 = b"ACGTACGTACGTACGTBBBBBBBBBBCCCCCCCCCCCCCMAPPABLE_BIO";
        let r2 = b"SECOND_READ_BIO";

        let tech = proto.extract_tech_seqs(r1, r2);
        assert_eq!(tech.barcode.unwrap().len(), 16);
        assert_eq!(tech.umi.unwrap().len(), 10);

        let reads = proto.extract_mappable_reads(r1, r2);
        // 5' protocol: R1 after tech prefix + R2
        assert_eq!(reads.seq1.unwrap(), b"MAPPABLE_BIO");
        assert_eq!(reads.seq2.unwrap(), b"SECOND_READ_BIO");
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
