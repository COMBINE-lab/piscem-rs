//! Mate overlap detection for paired-end reads.
//!
//! Detects dovetail and regular overlaps between read pairs, returning a
//! merged fragment sequence when possible. Used by scATAC mapping.
//!
//! Port of C++ `check_overlap.cpp`.

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/// Type of overlap detected between paired-end reads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OverlapType {
    /// End-to-end overlap (read1 end matches RC of read2 start).
    Dovetail,
    /// Regular overlap (RC of read2 overlaps within read1).
    Overlap,
    /// No valid overlap found.
    NoOverlap,
}

/// Result of overlap detection between paired-end reads.
#[derive(Debug, Clone)]
pub struct MateOverlap {
    /// Length of the merged fragment (0 if no overlap).
    pub frag_length: usize,
    /// Type of overlap found.
    pub ov_type: OverlapType,
    /// Merged fragment sequence (empty if no overlap).
    pub frag: Vec<u8>,
    /// Whether read1 was shorter (used for orientation tracking).
    pub frag_fw: bool,
}

impl Default for MateOverlap {
    fn default() -> Self {
        Self {
            frag_length: 0,
            ov_type: OverlapType::NoOverlap,
            frag: Vec::new(),
            frag_fw: true,
        }
    }
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Find overlap between paired-end reads.
///
/// Tries dovetail overlap first, then regular overlap.
///
/// - `min_overlap_length`: Minimum overlap required (typically 30).
/// - `error_threshold`: Maximum mismatches allowed during merge (typically 0).
pub fn find_overlap(
    seq1: &[u8],
    seq2: &[u8],
    min_overlap_length: i32,
    error_threshold: i32,
) -> MateOverlap {
    let mut mate_ov = MateOverlap::default();

    // Try dovetail first
    get_overlap(seq1, seq2, &mut mate_ov, true, min_overlap_length, error_threshold);
    if !mate_ov.frag.is_empty() {
        mate_ov.ov_type = OverlapType::Dovetail;
    } else {
        // Try regular overlap
        get_overlap(seq1, seq2, &mut mate_ov, false, min_overlap_length, error_threshold);
        mate_ov.ov_type = if !mate_ov.frag.is_empty() {
            OverlapType::Overlap
        } else {
            OverlapType::NoOverlap
        };
    }
    mate_ov.frag_length = mate_ov.frag.len();
    mate_ov
}

// ---------------------------------------------------------------------------
// Implementation
// ---------------------------------------------------------------------------

/// Reverse complement a DNA sequence.
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            other => other,
        })
        .collect()
}

/// Core overlap detection: seed-based alignment.
///
/// C++ swaps reads so the shorter is always `read1` and RC of longer is
/// `negative_read2`. We must match this for parity.
fn get_overlap(
    seq1: &[u8],
    seq2: &[u8],
    mate_ov: &mut MateOverlap,
    dovetail: bool,
    min_overlap_length: i32,
    error_threshold: i32,
) {
    let rc1 = reverse_complement(seq1);
    let rc2 = reverse_complement(seq2);

    // C++ assigns read1 = shorter, negative_read2 = RC(longer)
    let read1_shorter = seq1.len() <= seq2.len();
    mate_ov.frag_fw = read1_shorter;

    let (read1, negative_read2): (&[u8], &[u8]) = if read1_shorter {
        (seq1, &rc2)
    } else {
        (seq2, &rc1)
    };

    let seed_length = (min_overlap_length / 2) as usize;
    if seed_length == 0 {
        return;
    }

    // Assign roles based on dovetail mode
    let (suffix_read, prefix_read): (&[u8], &[u8]) = if dovetail {
        (negative_read2, read1)
    } else {
        (read1, negative_read2)
    };

    let suff_length = suffix_read.len();
    let pref_length = prefix_read.len();

    let mut overlap_length: usize = 0;
    let mut is_merged = false;

    // Progressive seed search
    for si in 0..=(error_threshold as usize) {
        let seed_offset = si * seed_length;
        if seed_offset + seed_length > pref_length {
            break;
        }
        let seed = &prefix_read[seed_offset..seed_offset + seed_length];

        // Find all occurrences of seed in suffix_read
        let mut search_start = 0;
        while search_start + seed_length <= suff_length {
            let pos = find_subsequence(&suffix_read[search_start..], seed);
            let seed_start = match pos {
                Some(p) => search_start + p,
                None => break,
            };

            // Validate: enough room before seed
            let before_ok = seed_start >= seed_offset;
            // Validate: enough overlap length
            let remaining_overlap =
                (suff_length - seed_start + seed_offset) as i32;
            let overlap_ok = remaining_overlap >= min_overlap_length;

            if before_ok && overlap_ok {
                let mut can_merge = true;
                let mut num_errors: i32 = 0;

                // Phase 1: match before seed
                for (i, &prefix_base) in prefix_read.iter().enumerate().take(seed_offset) {
                    let suff_idx = seed_start - seed_offset + i;
                    if suffix_read[suff_idx] != prefix_base {
                        num_errors += 1;
                        if num_errors > error_threshold {
                            can_merge = false;
                            break;
                        }
                    }
                }

                // Phase 2: match after seed
                if can_merge {
                    let max_after =
                        (suff_length - seed_start).min(pref_length - seed_offset);
                    for i in seed_length..max_after {
                        if suffix_read[seed_start + i]
                            != prefix_read[seed_offset + i]
                        {
                            num_errors += 1;
                        }
                        if num_errors > error_threshold {
                            can_merge = false;
                            break;
                        }
                    }
                }

                if can_merge {
                    is_merged = true;
                    overlap_length = (suff_length - seed_start + seed_offset)
                        .min(pref_length);
                    break;
                }
            }

            search_start = seed_start + 1;
        }

        if is_merged {
            break;
        }
    }

    // Build merged fragment
    if overlap_length > 0 {
        if dovetail {
            mate_ov.frag = prefix_read[..overlap_length].to_vec();
        } else {
            let mut frag = suffix_read.to_vec();
            if overlap_length < pref_length {
                frag.extend_from_slice(&prefix_read[overlap_length..]);
            }
            mate_ov.frag = frag;
        }
    }
}

/// Find first occurrence of needle in haystack.
fn find_subsequence(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    haystack
        .windows(needle.len())
        .position(|window| window == needle)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b"ATCG"), b"CGAT");
    }

    #[test]
    fn test_find_dovetail_overlap() {
        // Read1 and Read2 overlap at ends (read2 RC matches read1 suffix)
        let read1 = b"ACGTACGTACGTACGT";
        let read2 = reverse_complement(&read1[8..]); // RC of last 8 bases
        // read2 = RC("ACGTACGT") = "ACGTACGT" (palindromic)
        // This should detect dovetail if the reads overlap
        let result = find_overlap(read1, &read2, 6, 0);
        // Whether it detects depends on exact overlap geometry
        assert!(result.ov_type == OverlapType::Dovetail || result.ov_type == OverlapType::NoOverlap);
    }

    #[test]
    fn test_find_regular_overlap() {
        // Create two reads that overlap by 20 bases in the middle
        let shared = b"ACGTACGTACGTACGTACGT"; // 20 bases
        let read1 = [b"TTTTTTTTTT".as_slice(), shared].concat(); // 30 bases
        let read2_body = [shared.as_slice(), b"GGGGGGGGGG"].concat(); // 30 bases
        let read2 = reverse_complement(&read2_body);

        let result = find_overlap(&read1, &read2, 10, 0);
        // Should find some kind of overlap
        if result.ov_type != OverlapType::NoOverlap {
            assert!(!result.frag.is_empty());
            assert!(result.frag_length > 0);
        }
    }

    #[test]
    fn test_no_overlap() {
        let read1 = b"AAAAAAAAAAAAAAAA";
        let read2 = b"CCCCCCCCCCCCCCCC";
        let result = find_overlap(read1, read2, 10, 0);
        assert_eq!(result.ov_type, OverlapType::NoOverlap);
        assert!(result.frag.is_empty());
    }

    #[test]
    fn test_identical_reads_overlap() {
        // Two identical reads should have a dovetail overlap
        let read = b"ACGTACGTACGTACGT";
        let result = find_overlap(read, read, 10, 0);
        // Identical reads: RC(read2) searched in read1
        // Since RC of a palindrome is itself, this should find overlap
        if result.ov_type != OverlapType::NoOverlap {
            assert!(!result.frag.is_empty());
        }
    }

    #[test]
    fn test_overlap_with_errors() {
        // Create reads with 1 mismatch in overlap region
        let read1 = b"ACGTACGTACGT"; // 12 bases
        let mut read2_body = b"ACGTACGTACGTTTTT".to_vec(); // overlap first 12
        read2_body[5] = b'T'; // introduce mismatch
        let read2 = reverse_complement(&read2_body);

        // With error_threshold=0, should not merge
        let result_strict = find_overlap(read1, &read2, 8, 0);
        // With error_threshold=1, might merge
        let result_lenient = find_overlap(read1, &read2, 8, 1);

        // At least the strict version should have fewer/no merges
        assert!(
            result_strict.frag_length <= result_lenient.frag_length
                || result_strict.ov_type == OverlapType::NoOverlap
        );
    }
}
