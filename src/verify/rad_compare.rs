//! RAD file comparison.
//!
//! Compares two RAD files for semantic equivalence. Includes both a fast
//! header-only comparison and a full record-level multiset comparison using
//! libradicl to parse the RAD format.

use std::io::{Read, BufReader};
use std::path::Path;

use anyhow::{Context, Result};
use serde::Serialize;

/// Summary of RAD file header comparison (fast check).
#[derive(Debug, Serialize)]
pub struct RadComparisonSummary {
    pub passed: bool,
    pub notes: String,
    pub file_a_size: u64,
    pub file_b_size: u64,
    pub header_match: bool,
    pub num_refs_match: bool,
    pub num_chunks_a: u64,
    pub num_chunks_b: u64,
}

/// Parsed RAD header (minimal — just enough for comparison).
#[derive(Debug)]
struct RadHeader {
    is_paired: bool,
    num_refs: u64,
    ref_names: Vec<String>,
    num_chunks: u64,
}

/// Read a RAD header from a binary stream.
fn read_rad_header<R: Read>(reader: &mut R) -> Result<RadHeader> {
    // is_paired (u8)
    let mut buf1 = [0u8; 1];
    reader.read_exact(&mut buf1).context("reading is_paired")?;
    let is_paired = buf1[0] != 0;

    // num_refs (u64)
    let mut buf8 = [0u8; 8];
    reader.read_exact(&mut buf8).context("reading num_refs")?;
    let num_refs = u64::from_le_bytes(buf8);

    // ref_names: num_refs × (u16 len + bytes)
    let mut ref_names = Vec::with_capacity(num_refs as usize);
    for i in 0..num_refs {
        let mut buf2 = [0u8; 2];
        reader
            .read_exact(&mut buf2)
            .with_context(|| format!("reading ref name length {}", i))?;
        let name_len = u16::from_le_bytes(buf2) as usize;
        let mut name_buf = vec![0u8; name_len];
        reader
            .read_exact(&mut name_buf)
            .with_context(|| format!("reading ref name {}", i))?;
        ref_names.push(String::from_utf8_lossy(&name_buf).to_string());
    }

    // num_chunks (u64)
    reader.read_exact(&mut buf8).context("reading num_chunks")?;
    let num_chunks = u64::from_le_bytes(buf8);

    Ok(RadHeader {
        is_paired,
        num_refs,
        ref_names,
        num_chunks,
    })
}

/// Compare two RAD files (header-only fast check).
///
/// Compares file sizes, headers (is_paired, num_refs, ref_names), and
/// chunk counts. Does not do record-level comparison.
pub fn compare_rad_files(path_a: &Path, path_b: &Path) -> Result<RadComparisonSummary> {
    let meta_a = std::fs::metadata(path_a)
        .with_context(|| format!("stat {}", path_a.display()))?;
    let meta_b = std::fs::metadata(path_b)
        .with_context(|| format!("stat {}", path_b.display()))?;

    let file_a = std::fs::File::open(path_a)
        .with_context(|| format!("open {}", path_a.display()))?;
    let file_b = std::fs::File::open(path_b)
        .with_context(|| format!("open {}", path_b.display()))?;

    let header_a = read_rad_header(&mut BufReader::new(&file_a))?;
    let header_b = read_rad_header(&mut BufReader::new(&file_b))?;

    let num_refs_match = header_a.num_refs == header_b.num_refs
        && header_a.ref_names == header_b.ref_names;

    let header_match = header_a.is_paired == header_b.is_paired && num_refs_match;

    let passed = header_match && header_a.num_chunks == header_b.num_chunks;

    let mut notes = Vec::new();
    if !header_match {
        notes.push("headers differ".to_string());
    }
    if header_a.num_chunks != header_b.num_chunks {
        notes.push(format!(
            "chunk count mismatch: {} vs {}",
            header_a.num_chunks, header_b.num_chunks
        ));
    }
    if passed {
        notes.push("headers and chunk counts match".to_string());
    }

    Ok(RadComparisonSummary {
        passed,
        notes: notes.join("; "),
        file_a_size: meta_a.len(),
        file_b_size: meta_b.len(),
        header_match,
        num_refs_match,
        num_chunks_a: header_a.num_chunks,
        num_chunks_b: header_b.num_chunks,
    })
}

/// Validate that a RAD file is readable and has a well-formed header.
///
/// Returns the number of references found in the header.
pub fn validate_rad_file(path: &Path) -> Result<u64> {
    let file = std::fs::File::open(path)
        .with_context(|| format!("open {}", path.display()))?;
    let header = read_rad_header(&mut BufReader::new(&file))?;

    if header.num_refs == 0 {
        anyhow::bail!("RAD file has 0 references");
    }

    Ok(header.num_refs)
}

// ---------------------------------------------------------------------------
// Record-level multiset comparison (cfg(test) only — uses libradicl dev-dep)
// ---------------------------------------------------------------------------

/// Summary of full record-level RAD comparison.
#[cfg(any(test, feature = "parity-test"))]
#[derive(Debug, Serialize)]
pub struct RecordComparisonSummary {
    pub header_match: bool,
    pub total_records_a: usize,
    pub total_records_b: usize,
    pub matching_records: usize,
    pub missing_in_a: usize,
    pub missing_in_b: usize,
    pub passed: bool,
    pub notes: String,
    pub first_mismatches: Vec<String>,
    /// Records with same target set but different positions/flen/ori
    pub same_targets_diff_detail: usize,
    /// Records with completely different target sets
    pub different_targets: usize,
    /// Records that only differ in fragment length
    pub diff_frag_len_only: usize,
    /// Records that only differ in position
    pub diff_pos_only: usize,
    /// Records that differ in number of alignments
    pub diff_num_alns: usize,
}

/// A canonicalized bulk RAD record for multiset comparison.
///
/// Each alignment is represented as a tuple of (ref_id, orientation, pos, frag_len)
/// and sorted for deterministic comparison.
#[cfg(any(test, feature = "parity-test"))]
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct CanonicalBulkRecord {
    pub frag_type: u8,
    pub alignments: Vec<(u32, u8, u32, u16)>, // (ref_id, ori, pos, frag_len)
}

#[cfg(any(test, feature = "parity-test"))]
impl std::fmt::Display for CanonicalBulkRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "frag_type={}, alns=[", self.frag_type)?;
        for (i, (r, o, p, fl)) in self.alignments.iter().enumerate() {
            if i > 0 {
                write!(f, ", ")?;
            }
            write!(f, "(ref={r}, ori={o}, pos={p}, flen={fl})")?;
        }
        write!(f, "]")
    }
}

/// Read all bulk RAD records from a file using libradicl.
///
/// Returns the header info and a Vec of canonicalized records.
#[cfg(any(test, feature = "parity-test"))]
pub fn read_bulk_rad_records(
    path: &Path,
) -> Result<(libradicl::header::RadPrelude, Vec<CanonicalBulkRecord>)> {
    use libradicl::chunk::Chunk;
    use libradicl::header::RadPrelude;
    use libradicl::record::{MappedRecord, PiscemBulkReadRecord, PiscemBulkRecordContext, RecordContext};

    let file = std::fs::File::open(path)
        .with_context(|| format!("open {}", path.display()))?;
    let mut reader = BufReader::new(file);

    let prelude = RadPrelude::from_bytes(&mut reader)
        .context("parsing RAD prelude")?;

    // Parse file-level tags (advances the reader past tag values)
    let _file_tag_map = prelude.file_tags.try_parse_tags_from_bytes(&mut reader)
        .context("parsing file-level tags")?;

    // Get record parsing context
    let ctx = PiscemBulkRecordContext::get_context_from_tag_section(
        &prelude.file_tags,
        &prelude.read_tags,
        &prelude.aln_tags,
    ).context("getting bulk record context")?;

    let num_chunks = prelude.hdr.num_chunks;
    let mut all_records = Vec::new();

    for _ in 0..num_chunks {
        let chunk = Chunk::<PiscemBulkReadRecord>::from_bytes(&mut reader, &ctx);
        for rec in &chunk.reads {
            if rec.is_empty() {
                continue;
            }
            let canonical = canonicalize_bulk_record(rec);
            all_records.push(canonical);
        }
    }

    Ok((prelude, all_records))
}

/// Convert a PiscemBulkReadRecord to a canonical form for comparison.
#[cfg(any(test, feature = "parity-test"))]
fn canonicalize_bulk_record(
    rec: &libradicl::record::PiscemBulkReadRecord,
) -> CanonicalBulkRecord {
    use libradicl::rad_types::MappedFragmentOrientation;

    let mut alignments: Vec<(u32, u8, u32, u16)> = rec
        .refs
        .iter()
        .zip(rec.dirs.iter())
        .zip(rec.positions.iter())
        .zip(rec.frag_lengths.iter())
        .map(|(((&ref_id, &dir), &pos), &flen)| {
            let ori: u8 = match dir {
                MappedFragmentOrientation::Forward => 0,
                MappedFragmentOrientation::Reverse => 1,
                MappedFragmentOrientation::ForwardForward => 2,
                MappedFragmentOrientation::ForwardReverse => 3,
                MappedFragmentOrientation::ReverseForward => 4,
                MappedFragmentOrientation::ReverseReverse => 5,
                MappedFragmentOrientation::Unknown => 6,
            };
            (ref_id, ori, pos, flen)
        })
        .collect();

    // Sort alignments for deterministic comparison
    alignments.sort();

    CanonicalBulkRecord {
        frag_type: rec.frag_type,
        alignments,
    }
}

/// Compare two bulk RAD files at the record level using multiset comparison.
///
/// Since thread interleaving produces non-deterministic chunk ordering, records
/// are compared as multisets (frequency maps). Handles different reference
/// orderings between files by translating ref IDs to a common namespace via
/// ref names.
#[cfg(any(test, feature = "parity-test"))]
pub fn compare_bulk_rad_full(
    path_a: &Path,
    path_b: &Path,
) -> Result<RecordComparisonSummary> {
    let (prelude_a, records_a) = read_bulk_rad_records(path_a)
        .with_context(|| format!("reading {}", path_a.display()))?;
    let (prelude_b, records_b) = read_bulk_rad_records(path_b)
        .with_context(|| format!("reading {}", path_b.display()))?;

    // Compare headers: is_paired, ref_count, and ref_names as a set
    let paired_match = prelude_a.hdr.is_paired == prelude_b.hdr.is_paired;
    let count_match = prelude_a.hdr.ref_count == prelude_b.hdr.ref_count;

    // Check if ref names are the same SET (allow different ordering)
    let names_a: std::collections::HashSet<&str> =
        prelude_a.hdr.ref_names.iter().map(|s| s.as_str()).collect();
    let names_b: std::collections::HashSet<&str> =
        prelude_b.hdr.ref_names.iter().map(|s| s.as_str()).collect();
    let names_match = names_a == names_b;

    let header_match = paired_match && count_match && names_match;

    // Build ref ID translation: B's ref_id → A's ref_id (via ref name)
    let name_to_id_a: std::collections::HashMap<&str, u32> = prelude_a
        .hdr
        .ref_names
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i as u32))
        .collect();

    let b_to_a_id: Vec<u32> = prelude_b
        .hdr
        .ref_names
        .iter()
        .map(|name| {
            name_to_id_a
                .get(name.as_str())
                .copied()
                .unwrap_or(u32::MAX)
        })
        .collect();

    let total_a = records_a.len();
    let total_b = records_b.len();

    // Build frequency map for A (already uses A's ref IDs)
    let mut freq_a: std::collections::HashMap<CanonicalBulkRecord, u64> =
        std::collections::HashMap::with_capacity(total_a);
    for rec in &records_a {
        *freq_a.entry(rec.clone()).or_insert(0) += 1;
    }

    // Build frequency map for B, translating ref IDs to A's namespace
    let mut freq_b: std::collections::HashMap<CanonicalBulkRecord, u64> =
        std::collections::HashMap::with_capacity(total_b);
    for rec in &records_b {
        let translated = translate_record(rec, &b_to_a_id);
        *freq_b.entry(translated).or_insert(0) += 1;
    }

    // Compare frequency maps
    let mut matching: usize = 0;
    let mut missing_in_a: usize = 0;
    let mut missing_in_b: usize = 0;
    let mut first_mismatches: Vec<String> = Vec::new();
    let max_mismatch_report = 10;

    // Check records in A against B
    for (rec, &count_a) in &freq_a {
        let count_b = freq_b.get(rec).copied().unwrap_or(0);
        if count_a == count_b {
            matching += count_a as usize;
        } else if count_a > count_b {
            matching += count_b as usize;
            let diff = (count_a - count_b) as usize;
            missing_in_b += diff;
            if first_mismatches.len() < max_mismatch_report {
                first_mismatches.push(format!(
                    "in A but not B ({diff}x): {rec}"
                ));
            }
        } else {
            matching += count_a as usize;
            let diff = (count_b - count_a) as usize;
            missing_in_a += diff;
            if first_mismatches.len() < max_mismatch_report {
                first_mismatches.push(format!(
                    "in B but not A ({diff}x): {rec}"
                ));
            }
        }
    }

    // Check records only in B (not in A at all)
    for (rec, &count_b) in &freq_b {
        if !freq_a.contains_key(rec) {
            missing_in_a += count_b as usize;
            if first_mismatches.len() < max_mismatch_report {
                first_mismatches.push(format!(
                    "in B but not A ({count_b}x): {rec}"
                ));
            }
        }
    }

    let passed = header_match && missing_in_a == 0 && missing_in_b == 0 && total_a == total_b;

    // --- Detailed mismatch categorization ---
    // For records in A but not in B, check if there's a B record with the same target set.
    // Build a target-set index for B's unmatched records.
    let mut b_by_targets: std::collections::HashMap<Vec<u32>, Vec<&CanonicalBulkRecord>> =
        std::collections::HashMap::new();
    for (rec, &count_b) in &freq_b {
        let count_a = freq_a.get(rec).copied().unwrap_or(0);
        if count_b > count_a {
            // This record has excess in B
            let fp = target_set_fingerprint(rec);
            b_by_targets.entry(fp).or_default().push(rec);
        }
    }

    let mut same_targets_diff_detail: usize = 0;
    let mut different_targets: usize = 0;
    let mut diff_frag_len_only: usize = 0;
    let mut diff_pos_only: usize = 0;
    let mut diff_num_alns: usize = 0;
    let mut detail_examples: Vec<String> = Vec::new();
    let max_detail_examples = 5;

    // For records in A not in B, look for coarse matches in B
    for (rec_a, &count_a) in &freq_a {
        let count_b = freq_b.get(rec_a).copied().unwrap_or(0);
        if count_a > count_b {
            let excess = (count_a - count_b) as usize;
            let fp = target_set_fingerprint(rec_a);
            if let Some(b_recs) = b_by_targets.get(&fp) {
                // Found a B record with the same target set!
                same_targets_diff_detail += excess;
                // Classify the difference
                for b_rec in b_recs.iter().take(1) {
                    let kind = classify_detail_diff(rec_a, b_rec);
                    match kind {
                        "pos_only" => diff_pos_only += excess,
                        "flen_only" => diff_frag_len_only += excess,
                        "num_alns" => diff_num_alns += excess,
                        _ => {}
                    }
                    if detail_examples.len() < max_detail_examples {
                        detail_examples.push(format!(
                            "diff={kind}: A={rec_a}\n           B={}",
                            b_recs[0]
                        ));
                    }
                }
            } else {
                different_targets += excess;
            }
        }
    }

    let mut notes_parts = Vec::new();
    if !header_match {
        if !paired_match {
            notes_parts.push("is_paired differs".to_string());
        }
        if !count_match {
            notes_parts.push(format!(
                "ref_count differs: {} vs {}",
                prelude_a.hdr.ref_count, prelude_b.hdr.ref_count
            ));
        }
        if !names_match {
            notes_parts.push("ref_names differ (as sets)".to_string());
        }
    }
    if total_a != total_b {
        notes_parts.push(format!("record count mismatch: {total_a} vs {total_b}"));
    }
    if missing_in_a > 0 {
        notes_parts.push(format!("{missing_in_a} records in B missing from A"));
    }
    if missing_in_b > 0 {
        notes_parts.push(format!("{missing_in_b} records in A missing from B"));
    }
    if passed {
        notes_parts.push(format!("all {matching} records match"));
    }

    // Add detail examples to first_mismatches
    for ex in &detail_examples {
        if first_mismatches.len() < max_mismatch_report + max_detail_examples {
            first_mismatches.push(ex.clone());
        }
    }

    Ok(RecordComparisonSummary {
        header_match,
        total_records_a: total_a,
        total_records_b: total_b,
        matching_records: matching,
        missing_in_a,
        missing_in_b,
        passed,
        notes: notes_parts.join("; "),
        first_mismatches,
        same_targets_diff_detail,
        different_targets,
        diff_frag_len_only,
        diff_pos_only,
        diff_num_alns,
    })
}

/// Translate ref IDs in a record from one namespace to another.
#[cfg(any(test, feature = "parity-test"))]
fn translate_record(
    rec: &CanonicalBulkRecord,
    id_map: &[u32],
) -> CanonicalBulkRecord {
    let mut translated = rec.clone();
    for aln in &mut translated.alignments {
        if (aln.0 as usize) < id_map.len() {
            aln.0 = id_map[aln.0 as usize];
        }
    }
    // Re-sort after translation (ref IDs changed)
    translated.alignments.sort();
    translated
}

/// Coarse fingerprint: just the sorted target IDs (ignoring pos, flen, ori).
#[cfg(any(test, feature = "parity-test"))]
fn target_set_fingerprint(rec: &CanonicalBulkRecord) -> Vec<u32> {
    let mut tids: Vec<u32> = rec.alignments.iter().map(|a| a.0).collect();
    tids.sort();
    tids
}

/// Classify how two records with the same target set differ.
#[cfg(any(test, feature = "parity-test"))]
fn classify_detail_diff(a: &CanonicalBulkRecord, b: &CanonicalBulkRecord) -> &'static str {
    if a.alignments.len() != b.alignments.len() {
        return "num_alns";
    }
    // Same number of alignments, same target IDs (caller guarantees)
    // Check if only positions differ
    let pos_differ = a.alignments.iter().zip(b.alignments.iter())
        .any(|(aa, bb)| aa.2 != bb.2);
    let flen_differ = a.alignments.iter().zip(b.alignments.iter())
        .any(|(aa, bb)| aa.3 != bb.3);
    let ori_differ = a.alignments.iter().zip(b.alignments.iter())
        .any(|(aa, bb)| aa.1 != bb.1);
    let frag_type_differ = a.frag_type != b.frag_type;

    if frag_type_differ {
        "frag_type"
    } else if pos_differ && !flen_differ && !ori_differ {
        "pos_only"
    } else if flen_differ && !pos_differ && !ori_differ {
        "flen_only"
    } else if pos_differ && flen_differ && !ori_differ {
        "pos_and_flen"
    } else if ori_differ {
        "orientation"
    } else {
        "unknown"
    }
}

// ---------------------------------------------------------------------------
// SC record-level comparison
// ---------------------------------------------------------------------------

/// A canonicalized SC RAD record for multiset comparison.
///
/// Each alignment is (ref_id, direction). Sorted for deterministic comparison.
#[cfg(any(test, feature = "parity-test"))]
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct CanonicalScRecord {
    pub bc: u64,
    pub umi: u64,
    pub alignments: Vec<(u32, bool)>, // (ref_id, forward)
}

#[cfg(any(test, feature = "parity-test"))]
impl std::fmt::Display for CanonicalScRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "bc={}, umi={}, alns=[", self.bc, self.umi)?;
        for (i, (r, d)) in self.alignments.iter().enumerate() {
            if i > 0 {
                write!(f, ", ")?;
            }
            write!(f, "(ref={r}, fw={d})")?;
        }
        write!(f, "]")
    }
}

/// Read all SC RAD records from a file using libradicl.
#[cfg(any(test, feature = "parity-test"))]
pub fn read_sc_rad_records(
    path: &Path,
) -> Result<(libradicl::header::RadPrelude, Vec<CanonicalScRecord>)> {
    use libradicl::chunk::Chunk;
    use libradicl::header::RadPrelude;
    use libradicl::record::{AlevinFryReadRecord, AlevinFryRecordContext, MappedRecord, RecordContext};

    let file = std::fs::File::open(path)
        .with_context(|| format!("open {}", path.display()))?;
    let mut reader = BufReader::new(file);

    let prelude = RadPrelude::from_bytes(&mut reader)
        .context("parsing RAD prelude")?;

    // Parse file-level tags
    let _file_tag_map = prelude.file_tags.try_parse_tags_from_bytes(&mut reader)
        .context("parsing file-level tags")?;

    // Get record parsing context
    let ctx = AlevinFryRecordContext::get_context_from_tag_section(
        &prelude.file_tags,
        &prelude.read_tags,
        &prelude.aln_tags,
    ).context("getting SC record context")?;

    let num_chunks = prelude.hdr.num_chunks;
    let mut all_records = Vec::new();

    for _ in 0..num_chunks {
        let chunk = Chunk::<AlevinFryReadRecord>::from_bytes(&mut reader, &ctx);
        for rec in &chunk.reads {
            if rec.is_empty() {
                continue;
            }
            let mut alignments: Vec<(u32, bool)> = rec.refs.iter()
                .zip(rec.dirs.iter())
                .map(|(&ref_id, &dir)| (ref_id, dir))
                .collect();
            alignments.sort();

            all_records.push(CanonicalScRecord {
                bc: rec.bc,
                umi: rec.umi,
                alignments,
            });
        }
    }

    Ok((prelude, all_records))
}

/// Compare two SC RAD files at the record level using multiset comparison.
#[cfg(any(test, feature = "parity-test"))]
pub fn compare_sc_rad_full(
    path_a: &Path,
    path_b: &Path,
) -> Result<RecordComparisonSummary> {
    let (prelude_a, records_a) = read_sc_rad_records(path_a)
        .with_context(|| format!("reading {}", path_a.display()))?;
    let (prelude_b, records_b) = read_sc_rad_records(path_b)
        .with_context(|| format!("reading {}", path_b.display()))?;

    // Compare headers
    let paired_match = prelude_a.hdr.is_paired == prelude_b.hdr.is_paired;
    let count_match = prelude_a.hdr.ref_count == prelude_b.hdr.ref_count;
    let names_a: std::collections::HashSet<&str> =
        prelude_a.hdr.ref_names.iter().map(|s| s.as_str()).collect();
    let names_b: std::collections::HashSet<&str> =
        prelude_b.hdr.ref_names.iter().map(|s| s.as_str()).collect();
    let names_match = names_a == names_b;
    let header_match = paired_match && count_match && names_match;

    // Build ref ID translation: B's ref_id → A's ref_id (via ref name)
    let name_to_id_a: std::collections::HashMap<&str, u32> = prelude_a
        .hdr.ref_names.iter().enumerate()
        .map(|(i, n)| (n.as_str(), i as u32))
        .collect();
    let b_to_a_id: Vec<u32> = prelude_b
        .hdr.ref_names.iter()
        .map(|name| name_to_id_a.get(name.as_str()).copied().unwrap_or(u32::MAX))
        .collect();

    let total_a = records_a.len();
    let total_b = records_b.len();

    // Build frequency map for A
    let mut freq_a: std::collections::HashMap<CanonicalScRecord, u64> =
        std::collections::HashMap::with_capacity(total_a);
    for rec in &records_a {
        *freq_a.entry(rec.clone()).or_insert(0) += 1;
    }

    // Build frequency map for B, translating ref IDs
    let mut freq_b: std::collections::HashMap<CanonicalScRecord, u64> =
        std::collections::HashMap::with_capacity(total_b);
    for rec in &records_b {
        let translated = translate_sc_record(rec, &b_to_a_id);
        *freq_b.entry(translated).or_insert(0) += 1;
    }

    // Compare frequency maps
    let mut matching: usize = 0;
    let mut missing_in_a: usize = 0;
    let mut missing_in_b: usize = 0;
    let mut first_mismatches: Vec<String> = Vec::new();
    let max_mismatch_report = 10;

    for (rec, &count_a) in &freq_a {
        let count_b = freq_b.get(rec).copied().unwrap_or(0);
        if count_a == count_b {
            matching += count_a as usize;
        } else if count_a > count_b {
            matching += count_b as usize;
            let diff = (count_a - count_b) as usize;
            missing_in_b += diff;
            if first_mismatches.len() < max_mismatch_report {
                first_mismatches.push(format!("in A but not B ({diff}x): {rec}"));
            }
        } else {
            matching += count_a as usize;
            let diff = (count_b - count_a) as usize;
            missing_in_a += diff;
            if first_mismatches.len() < max_mismatch_report {
                first_mismatches.push(format!("in B but not A ({diff}x): {rec}"));
            }
        }
    }
    for (rec, &count_b) in &freq_b {
        if !freq_a.contains_key(rec) {
            missing_in_a += count_b as usize;
            if first_mismatches.len() < max_mismatch_report {
                first_mismatches.push(format!("in B but not A ({count_b}x): {rec}"));
            }
        }
    }

    let passed = header_match && missing_in_a == 0 && missing_in_b == 0 && total_a == total_b;

    let mut notes_parts = Vec::new();
    if !header_match {
        if !paired_match { notes_parts.push("is_paired differs".to_string()); }
        if !count_match {
            notes_parts.push(format!(
                "ref_count differs: {} vs {}", prelude_a.hdr.ref_count, prelude_b.hdr.ref_count
            ));
        }
        if !names_match { notes_parts.push("ref_names differ".to_string()); }
    }
    if total_a != total_b {
        notes_parts.push(format!("record count mismatch: {total_a} vs {total_b}"));
    }
    if missing_in_a > 0 { notes_parts.push(format!("{missing_in_a} records in B missing from A")); }
    if missing_in_b > 0 { notes_parts.push(format!("{missing_in_b} records in A missing from B")); }
    if passed { notes_parts.push(format!("all {matching} records match")); }

    Ok(RecordComparisonSummary {
        header_match,
        total_records_a: total_a,
        total_records_b: total_b,
        matching_records: matching,
        missing_in_a,
        missing_in_b,
        passed,
        notes: notes_parts.join("; "),
        first_mismatches,
        // SC records don't have position/frag_len detail categories
        same_targets_diff_detail: 0,
        different_targets: 0,
        diff_frag_len_only: 0,
        diff_pos_only: 0,
        diff_num_alns: 0,
    })
}

/// Translate ref IDs in an SC record from one namespace to another.
#[cfg(any(test, feature = "parity-test"))]
fn translate_sc_record(
    rec: &CanonicalScRecord,
    id_map: &[u32],
) -> CanonicalScRecord {
    let mut translated = rec.clone();
    for aln in &mut translated.alignments {
        if (aln.0 as usize) < id_map.len() {
            aln.0 = id_map[aln.0 as usize];
        }
    }
    translated.alignments.sort();
    translated
}

// ---------------------------------------------------------------------------
// ATAC record-level comparison (custom parser — libradicl lacks ATAC support)
// ---------------------------------------------------------------------------

/// A canonicalized ATAC RAD record for multiset comparison.
///
/// Each alignment is (ref_id, type_code, start_pos, frag_len).
/// Sorted for deterministic comparison.
#[cfg(any(test, feature = "parity-test"))]
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct CanonicalAtacRecord {
    pub bc: u64,
    pub alignments: Vec<(u32, u8, u32, u16)>, // (ref_id, type, start_pos, frag_len)
}

#[cfg(any(test, feature = "parity-test"))]
impl std::fmt::Display for CanonicalAtacRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "bc={}, alns=[", self.bc)?;
        for (i, (r, t, p, fl)) in self.alignments.iter().enumerate() {
            if i > 0 {
                write!(f, ", ")?;
            }
            write!(f, "(ref={r}, type={t}, pos={p}, flen={fl})")?;
        }
        write!(f, "]")
    }
}

/// Parsed ATAC RAD header with tag type metadata needed for record parsing.
#[cfg(any(test, feature = "parity-test"))]
#[derive(Debug)]
pub struct AtacRadHeader {
    pub is_paired: bool,
    pub num_refs: u64,
    pub ref_names: Vec<String>,
    pub num_chunks: u64,
    pub bc_tag_type: u8, // TAG_U32=3 or TAG_U64=4
}

/// Read an ATAC RAD header, including tag descriptions and file-level tag values.
///
/// This is a custom parser because libradicl doesn't support ATAC records.
/// Handles both C++ (v_u64 ref_lengths) and Rust (array_u32 ref_lengths) formats.
#[cfg(any(test, feature = "parity-test"))]
fn read_atac_rad_header<R: Read>(reader: &mut R) -> Result<AtacRadHeader> {
    // is_paired (u8)
    let mut buf1 = [0u8; 1];
    reader.read_exact(&mut buf1).context("reading is_paired")?;
    let is_paired = buf1[0] != 0;

    // num_refs (u64)
    let mut buf8 = [0u8; 8];
    reader.read_exact(&mut buf8).context("reading num_refs")?;
    let num_refs = u64::from_le_bytes(buf8);

    // ref_names
    let mut ref_names = Vec::with_capacity(num_refs as usize);
    let mut buf2 = [0u8; 2];
    for _ in 0..num_refs {
        reader.read_exact(&mut buf2).context("reading ref name length")?;
        let name_len = u16::from_le_bytes(buf2) as usize;
        let mut name_buf = vec![0u8; name_len];
        reader.read_exact(&mut name_buf).context("reading ref name")?;
        ref_names.push(String::from_utf8_lossy(&name_buf).to_string());
    }

    // num_chunks (u64)
    reader.read_exact(&mut buf8).context("reading num_chunks")?;
    let num_chunks = u64::from_le_bytes(buf8);

    // --- Tag descriptions ---
    // We need to parse these to know the barcode type and to skip past them.

    // File-level tags
    reader.read_exact(&mut buf2).context("reading num_file_tags")?;
    let num_file_tags = u16::from_le_bytes(buf2);
    let mut file_tag_descs = Vec::new();
    for _ in 0..num_file_tags {
        let (name, type_id) = read_tag_desc(reader)?;
        file_tag_descs.push((name, type_id));
    }

    // Read-level tags
    reader.read_exact(&mut buf2).context("reading num_read_tags")?;
    let num_read_tags = u16::from_le_bytes(buf2);
    let mut bc_tag_type = 3u8; // default TAG_U32
    for _ in 0..num_read_tags {
        let (name, type_id) = read_tag_desc(reader)?;
        if name == "b" {
            bc_tag_type = type_id;
        }
    }

    // Alignment-level tags
    reader.read_exact(&mut buf2).context("reading num_aln_tags")?;
    let num_aln_tags = u16::from_le_bytes(buf2);
    for _ in 0..num_aln_tags {
        let (_name, _type_id) = read_tag_desc(reader)?;
    }

    // --- File-level tag values ---
    // We need to skip these. Parse based on the tag types we collected.
    for (_name, type_id) in &file_tag_descs {
        skip_tag_value(reader, *type_id)?;
    }

    Ok(AtacRadHeader {
        is_paired,
        num_refs,
        ref_names,
        num_chunks,
        bc_tag_type,
    })
}

/// Read a tag description: name (u16 len + bytes) + type_id (u8).
/// For array tags (type=7), also reads len_type (u8) + elem_type (u8).
#[cfg(any(test, feature = "parity-test"))]
fn read_tag_desc<R: Read>(reader: &mut R) -> Result<(String, u8)> {
    let mut buf2 = [0u8; 2];
    reader.read_exact(&mut buf2).context("reading tag name length")?;
    let name_len = u16::from_le_bytes(buf2) as usize;
    let mut name_buf = vec![0u8; name_len];
    reader.read_exact(&mut name_buf).context("reading tag name")?;
    let name = String::from_utf8_lossy(&name_buf).to_string();

    let mut buf1 = [0u8; 1];
    reader.read_exact(&mut buf1).context("reading tag type")?;
    let type_id = buf1[0];

    // Array tags (TAG_ARRAY=7) have two additional bytes: len_type + elem_type
    if type_id == 7 {
        let mut extra = [0u8; 2];
        reader.read_exact(&mut extra).context("reading array len/elem types")?;
        let _len_type = extra[0]; // type of count field (unused, always u32 in practice)
        let elem_type = extra[1]; // type of each element
        // We mark this as type 100+elem_type for internal tracking
        return Ok((name, 100 + elem_type));
    }

    Ok((name, type_id))
}

/// Skip a tag value based on its type.
#[cfg(any(test, feature = "parity-test"))]
fn skip_tag_value<R: Read>(reader: &mut R, type_id: u8) -> Result<()> {
    let mut buf = [0u8; 8];
    match type_id {
        1 => { reader.read_exact(&mut buf[..1])?; } // u8
        2 => { reader.read_exact(&mut buf[..2])?; } // u16
        3 => { reader.read_exact(&mut buf[..4])?; } // u32
        4 => { reader.read_exact(&mut buf[..8])?; } // u64
        8 => {
            // string: u16 len + bytes
            let mut buf2 = [0u8; 2];
            reader.read_exact(&mut buf2)?;
            let len = u16::from_le_bytes(buf2) as usize;
            let mut skip = vec![0u8; len];
            reader.read_exact(&mut skip)?;
        }
        t if t >= 100 => {
            // Array: the count field size depends on the len_type from the tag desc.
            // C++ libradicl uses u64 for v_u64 count: `add(Type::u64(val.val().size()))`
            // Rust uses u32 count: `write_u32(ref_lengths.len() as u32)`
            // We read the len_type to determine count size, but since we encoded
            // len_type into type_id too, we can't distinguish here. Instead, try
            // the element type as a proxy: if elements are u64, count is u64 (C++ style).
            let elem_type = t - 100;
            let elem_size = match elem_type {
                1 => 1usize, // u8
                2 => 2,      // u16
                3 => 4,      // u32
                4 => 8,      // u64
                _ => anyhow::bail!("unknown array element type {}", elem_type),
            };
            // C++ v_u64 uses u64 count; Rust uses u32 count.
            // Determine based on elem_type: u64 elements => u64 count (C++ pattern)
            let count = if elem_type == 4 {
                let mut buf8 = [0u8; 8];
                reader.read_exact(&mut buf8)?;
                u64::from_le_bytes(buf8) as usize
            } else {
                let mut buf4 = [0u8; 4];
                reader.read_exact(&mut buf4)?;
                u32::from_le_bytes(buf4) as usize
            };
            let mut skip = vec![0u8; count * elem_size];
            reader.read_exact(&mut skip)?;
        }
        _ => anyhow::bail!("unknown tag type {}", type_id),
    }
    Ok(())
}

/// Read all ATAC RAD records from a file using a custom parser.
///
/// Returns the header info and a Vec of canonicalized records.
#[cfg(any(test, feature = "parity-test"))]
pub fn read_atac_rad_records(
    path: &Path,
) -> Result<(AtacRadHeader, Vec<CanonicalAtacRecord>)> {
    let file = std::fs::File::open(path)
        .with_context(|| format!("open {}", path.display()))?;
    let mut reader = BufReader::new(file);

    let header = read_atac_rad_header(&mut reader)?;

    let bc_is_u64 = header.bc_tag_type == 4; // TAG_U64

    let mut all_records = Vec::new();

    // num_chunks=0 means "unknown" (C++ libradicl never backpatches it).
    // In that case, read until EOF.
    let known_chunks = header.num_chunks > 0;
    let mut chunks_read: u64 = 0;

    loop {
        if known_chunks && chunks_read >= header.num_chunks {
            break;
        }

        // Chunk header: num_bytes (u32) + num_reads (u32)
        // C++ chunk_sz includes the 8-byte header itself.
        let mut buf4 = [0u8; 4];
        match reader.read_exact(&mut buf4) {
            Ok(()) => {}
            Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break,
            Err(e) => return Err(e).context("reading chunk num_bytes"),
        }
        let num_bytes = u32::from_le_bytes(buf4) as usize;

        reader.read_exact(&mut buf4)
            .with_context(|| format!("reading chunk {} num_reads", chunks_read))?;
        let num_reads = u32::from_le_bytes(buf4);

        // Read the entire chunk data as a byte buffer (num_bytes - 8 for header)
        let data_len = num_bytes.saturating_sub(8);
        let mut chunk_data = vec![0u8; data_len];
        reader.read_exact(&mut chunk_data)
            .with_context(|| format!("reading chunk {} data ({} bytes)", chunks_read, data_len))?;

        // Parse records from the chunk buffer
        let mut cursor = std::io::Cursor::new(&chunk_data);
        for _ in 0..num_reads {
            use std::io::Read as _;
            let mut b4 = [0u8; 4];

            // num_mappings (u32)
            cursor.read_exact(&mut b4)?;
            let num_mappings = u32::from_le_bytes(b4);

            // barcode
            let bc: u64 = if bc_is_u64 {
                let mut b8 = [0u8; 8];
                cursor.read_exact(&mut b8)?;
                u64::from_le_bytes(b8)
            } else {
                cursor.read_exact(&mut b4)?;
                u32::from_le_bytes(b4) as u64
            };

            // Per-alignment: ref (u32) + type (u8) + start_pos (u32) + frag_len (u16)
            let mut alignments = Vec::with_capacity(num_mappings as usize);
            for _ in 0..num_mappings {
                cursor.read_exact(&mut b4)?;
                let ref_id = u32::from_le_bytes(b4);

                let mut b1 = [0u8; 1];
                cursor.read_exact(&mut b1)?;
                let atype = b1[0];

                cursor.read_exact(&mut b4)?;
                let start_pos = u32::from_le_bytes(b4);

                let mut b2 = [0u8; 2];
                cursor.read_exact(&mut b2)?;
                let frag_len = u16::from_le_bytes(b2);

                alignments.push((ref_id, atype, start_pos, frag_len));
            }

            if num_mappings > 0 {
                alignments.sort();
                all_records.push(CanonicalAtacRecord { bc, alignments });
            }
        }

        chunks_read += 1;
    }

    Ok((header, all_records))
}

/// Compare two ATAC RAD files at the record level using multiset comparison.
///
/// Since thread interleaving produces non-deterministic chunk ordering, records
/// are compared as multisets. Handles different reference orderings.
#[cfg(any(test, feature = "parity-test"))]
pub fn compare_atac_rad_full(
    path_a: &Path,
    path_b: &Path,
) -> Result<RecordComparisonSummary> {
    let (header_a, records_a) = read_atac_rad_records(path_a)
        .with_context(|| format!("reading {}", path_a.display()))?;
    let (header_b, records_b) = read_atac_rad_records(path_b)
        .with_context(|| format!("reading {}", path_b.display()))?;

    // Compare headers
    let paired_match = header_a.is_paired == header_b.is_paired;
    let count_match = header_a.num_refs == header_b.num_refs;
    let names_a: std::collections::HashSet<&str> =
        header_a.ref_names.iter().map(|s| s.as_str()).collect();
    let names_b: std::collections::HashSet<&str> =
        header_b.ref_names.iter().map(|s| s.as_str()).collect();
    let names_match = names_a == names_b;
    let header_match = paired_match && count_match && names_match;

    // Build ref ID translation: B's ref_id → A's ref_id
    let name_to_id_a: std::collections::HashMap<&str, u32> = header_a
        .ref_names.iter().enumerate()
        .map(|(i, n)| (n.as_str(), i as u32))
        .collect();
    let b_to_a_id: Vec<u32> = header_b
        .ref_names.iter()
        .map(|name| name_to_id_a.get(name.as_str()).copied().unwrap_or(u32::MAX))
        .collect();

    let total_a = records_a.len();
    let total_b = records_b.len();

    // Build frequency map for A
    let mut freq_a: std::collections::HashMap<CanonicalAtacRecord, u64> =
        std::collections::HashMap::with_capacity(total_a);
    for rec in &records_a {
        *freq_a.entry(rec.clone()).or_insert(0) += 1;
    }

    // Build frequency map for B, translating ref IDs
    let mut freq_b: std::collections::HashMap<CanonicalAtacRecord, u64> =
        std::collections::HashMap::with_capacity(total_b);
    for rec in &records_b {
        let translated = translate_atac_record(rec, &b_to_a_id);
        *freq_b.entry(translated).or_insert(0) += 1;
    }

    // Compare frequency maps
    let mut matching: usize = 0;
    let mut missing_in_a: usize = 0;
    let mut missing_in_b: usize = 0;
    let mut first_mismatches: Vec<String> = Vec::new();
    let max_mismatch_report = 10;

    for (rec, &count_a) in &freq_a {
        let count_b = freq_b.get(rec).copied().unwrap_or(0);
        if count_a == count_b {
            matching += count_a as usize;
        } else if count_a > count_b {
            matching += count_b as usize;
            let diff = (count_a - count_b) as usize;
            missing_in_b += diff;
            if first_mismatches.len() < max_mismatch_report {
                first_mismatches.push(format!("in A but not B ({diff}x): {rec}"));
            }
        } else {
            matching += count_a as usize;
            let diff = (count_b - count_a) as usize;
            missing_in_a += diff;
            if first_mismatches.len() < max_mismatch_report {
                first_mismatches.push(format!("in B but not A ({diff}x): {rec}"));
            }
        }
    }
    for (rec, &count_b) in &freq_b {
        if !freq_a.contains_key(rec) {
            missing_in_a += count_b as usize;
            if first_mismatches.len() < max_mismatch_report {
                first_mismatches.push(format!("in B but not A ({count_b}x): {rec}"));
            }
        }
    }

    let passed = header_match && missing_in_a == 0 && missing_in_b == 0 && total_a == total_b;

    // Detailed mismatch categorization
    let mut b_by_targets: std::collections::HashMap<Vec<u32>, Vec<&CanonicalAtacRecord>> =
        std::collections::HashMap::new();
    for (rec, &count_b) in &freq_b {
        let count_a = freq_a.get(rec).copied().unwrap_or(0);
        if count_b > count_a {
            let fp: Vec<u32> = rec.alignments.iter().map(|a| a.0).collect();
            b_by_targets.entry(fp).or_default().push(rec);
        }
    }

    let mut same_targets_diff_detail: usize = 0;
    let mut different_targets: usize = 0;
    let mut diff_frag_len_only: usize = 0;
    let mut diff_pos_only: usize = 0;
    let mut diff_num_alns: usize = 0;

    for (rec_a, &count_a) in &freq_a {
        let count_b = freq_b.get(rec_a).copied().unwrap_or(0);
        if count_a > count_b {
            let excess = (count_a - count_b) as usize;
            let fp: Vec<u32> = rec_a.alignments.iter().map(|a| a.0).collect();
            if let Some(b_recs) = b_by_targets.get(&fp) {
                same_targets_diff_detail += excess;
                for b_rec in b_recs.iter().take(1) {
                    if rec_a.alignments.len() != b_rec.alignments.len() {
                        diff_num_alns += excess;
                    } else {
                        let pos_differ = rec_a.alignments.iter().zip(b_rec.alignments.iter())
                            .any(|(a, b)| a.2 != b.2);
                        let flen_differ = rec_a.alignments.iter().zip(b_rec.alignments.iter())
                            .any(|(a, b)| a.3 != b.3);
                        if pos_differ && !flen_differ {
                            diff_pos_only += excess;
                        } else if flen_differ && !pos_differ {
                            diff_frag_len_only += excess;
                        }
                    }
                }
            } else {
                different_targets += excess;
            }
        }
    }

    let mut notes_parts = Vec::new();
    if !header_match {
        if !paired_match { notes_parts.push("is_paired differs".to_string()); }
        if !count_match {
            notes_parts.push(format!("ref_count differs: {} vs {}", header_a.num_refs, header_b.num_refs));
        }
        if !names_match { notes_parts.push("ref_names differ".to_string()); }
    }
    if total_a != total_b {
        notes_parts.push(format!("record count mismatch: {total_a} vs {total_b}"));
    }
    if missing_in_a > 0 { notes_parts.push(format!("{missing_in_a} records in B missing from A")); }
    if missing_in_b > 0 { notes_parts.push(format!("{missing_in_b} records in A missing from B")); }
    if passed { notes_parts.push(format!("all {matching} records match")); }

    Ok(RecordComparisonSummary {
        header_match,
        total_records_a: total_a,
        total_records_b: total_b,
        matching_records: matching,
        missing_in_a,
        missing_in_b,
        passed,
        notes: notes_parts.join("; "),
        first_mismatches,
        same_targets_diff_detail,
        different_targets,
        diff_frag_len_only,
        diff_pos_only,
        diff_num_alns,
    })
}

/// Translate ref IDs in an ATAC record from one namespace to another.
#[cfg(any(test, feature = "parity-test"))]
fn translate_atac_record(
    rec: &CanonicalAtacRecord,
    id_map: &[u32],
) -> CanonicalAtacRecord {
    let mut translated = rec.clone();
    for aln in &mut translated.alignments {
        if (aln.0 as usize) < id_map.len() {
            aln.0 = id_map[aln.0 as usize];
        }
    }
    translated.alignments.sort();
    translated
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::rad::{
        RadWriter, write_bulk_record, write_rad_header_bulk,
    };
    use crate::mapping::hits::{MappingType, SimpleHit};
    use std::io::{Cursor, Write};

    #[test]
    fn test_read_rad_header_bulk() {
        let mut buf = Vec::new();
        let names = vec!["ref1", "ref2"];
        let ref_lens = vec![1000u32, 2000u32];
        let _chunk_off = write_rad_header_bulk(&mut buf, true, 2, &names, &ref_lens).unwrap();

        let header = read_rad_header(&mut Cursor::new(&buf)).unwrap();
        assert!(header.is_paired);
        assert_eq!(header.num_refs, 2);
        assert_eq!(header.ref_names, vec!["ref1", "ref2"]);
    }

    #[test]
    fn test_rad_comparison_summary_serialize() {
        let summary = RadComparisonSummary {
            passed: true,
            notes: "ok".to_string(),
            file_a_size: 1000,
            file_b_size: 1000,
            header_match: true,
            num_refs_match: true,
            num_chunks_a: 5,
            num_chunks_b: 5,
        };
        let json = serde_json::to_string(&summary).unwrap();
        assert!(json.contains("\"passed\":true"));
    }

    #[test]
    fn test_validate_rad_header() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.rad");
        let names = vec!["ref1"];
        let ref_lens = vec![500u32];
        let mut file = std::fs::File::create(&path).unwrap();
        write_rad_header_bulk(&mut file, false, 1, &names, &ref_lens).unwrap();
        drop(file);

        let num_refs = validate_rad_file(&path).unwrap();
        assert_eq!(num_refs, 1);
    }

    // --- Helper to write a complete bulk RAD file for testing ---

    fn write_test_bulk_rad_file(
        path: &Path,
        is_paired: bool,
        ref_names: &[&str],
        ref_lens: &[u32],
        chunks: &[Vec<(MappingType, Vec<SimpleHit>)>],
    ) {
        let mut file = std::fs::File::create(path).unwrap();
        let chunk_count_offset = write_rad_header_bulk(
            &mut file,
            is_paired,
            ref_names.len() as u64,
            ref_names,
            ref_lens,
        )
        .unwrap();

        for chunk_records in chunks {
            let mut rad = RadWriter::new();
            // Chunk header: num_bytes (u32) + num_reads (u32)
            rad.write_u32(0); // placeholder
            rad.write_u32(0); // placeholder
            let mut count: u32 = 0;
            for (mt, hits) in chunk_records {
                write_bulk_record(*mt, hits, &mut rad);
                count += 1;
            }
            let total_bytes = rad.len() as u32;
            rad.write_u32_at_offset(0, total_bytes);
            rad.write_u32_at_offset(4, count);
            file.write_all(rad.as_bytes()).unwrap();
        }

        // Backpatch num_chunks
        use std::io::{Seek, SeekFrom};
        let num_chunks = chunks.len() as u64;
        file.seek(SeekFrom::Start(chunk_count_offset)).unwrap();
        file.write_all(&num_chunks.to_le_bytes()).unwrap();
    }

    #[test]
    fn test_identical_bulk_files_pass() {
        let dir = tempfile::tempdir().unwrap();
        let path_a = dir.path().join("a.rad");
        let path_b = dir.path().join("b.rad");

        let names = vec!["ref1", "ref2"];
        let lens = vec![1000u32, 2000u32];
        let hits = vec![SimpleHit {
            is_fw: true,
            mate_is_fw: false,
            tid: 0,
            pos: 100,
            mate_pos: 200,
            fragment_length: 150,
            ..SimpleHit::default()
        }];
        let chunks = vec![vec![(MappingType::MappedPair, hits)]];

        write_test_bulk_rad_file(&path_a, true, &names, &lens, &chunks);
        write_test_bulk_rad_file(&path_b, true, &names, &lens, &chunks);

        let result = compare_bulk_rad_full(&path_a, &path_b).unwrap();
        assert!(result.passed, "Expected pass: {}", result.notes);
        assert!(result.header_match);
        assert_eq!(result.total_records_a, 1);
        assert_eq!(result.total_records_b, 1);
        assert_eq!(result.matching_records, 1);
        assert_eq!(result.missing_in_a, 0);
        assert_eq!(result.missing_in_b, 0);
    }

    #[test]
    fn test_different_record_counts_fail() {
        let dir = tempfile::tempdir().unwrap();
        let path_a = dir.path().join("a.rad");
        let path_b = dir.path().join("b.rad");

        let names = vec!["ref1"];
        let lens = vec![1000u32];
        let hit1 = SimpleHit {
            is_fw: true,
            tid: 0,
            pos: 100,
            ..SimpleHit::default()
        };
        let hit2 = SimpleHit {
            is_fw: false,
            tid: 0,
            pos: 200,
            ..SimpleHit::default()
        };

        let chunks_a = vec![vec![
            (MappingType::SingleMapped, vec![hit1]),
            (MappingType::SingleMapped, vec![hit2]),
        ]];
        let chunks_b = vec![vec![(MappingType::SingleMapped, vec![hit1])]];

        write_test_bulk_rad_file(&path_a, false, &names, &lens, &chunks_a);
        write_test_bulk_rad_file(&path_b, false, &names, &lens, &chunks_b);

        let result = compare_bulk_rad_full(&path_a, &path_b).unwrap();
        assert!(!result.passed);
        assert_eq!(result.total_records_a, 2);
        assert_eq!(result.total_records_b, 1);
        assert_eq!(result.missing_in_b, 1);
    }

    #[test]
    fn test_same_count_different_records_detected() {
        let dir = tempfile::tempdir().unwrap();
        let path_a = dir.path().join("a.rad");
        let path_b = dir.path().join("b.rad");

        let names = vec!["ref1"];
        let lens = vec![1000u32];
        let hit_a = SimpleHit {
            is_fw: true,
            tid: 0,
            pos: 100,
            ..SimpleHit::default()
        };
        let hit_b = SimpleHit {
            is_fw: true,
            tid: 0,
            pos: 999, // different position
            ..SimpleHit::default()
        };

        let chunks_a = vec![vec![(MappingType::SingleMapped, vec![hit_a])]];
        let chunks_b = vec![vec![(MappingType::SingleMapped, vec![hit_b])]];

        write_test_bulk_rad_file(&path_a, false, &names, &lens, &chunks_a);
        write_test_bulk_rad_file(&path_b, false, &names, &lens, &chunks_b);

        let result = compare_bulk_rad_full(&path_a, &path_b).unwrap();
        assert!(!result.passed);
        assert!(!result.first_mismatches.is_empty());
    }

    #[test]
    fn test_canonicalization_alignment_order_irrelevant() {
        // Two records with same alignments in different order should match.
        // We test this by creating two files with records in different chunks
        // (simulating thread interleaving).
        let dir = tempfile::tempdir().unwrap();
        let path_a = dir.path().join("a.rad");
        let path_b = dir.path().join("b.rad");

        let names = vec!["ref1", "ref2"];
        let lens = vec![1000u32, 2000u32];
        let hit1 = SimpleHit {
            is_fw: true,
            tid: 0,
            pos: 50,
            ..SimpleHit::default()
        };
        let hit2 = SimpleHit {
            is_fw: false,
            tid: 0,
            pos: 200,
            ..SimpleHit::default()
        };

        // File A: records in chunk 1
        let chunks_a = vec![vec![
            (MappingType::SingleMapped, vec![hit1]),
            (MappingType::SingleMapped, vec![hit2]),
        ]];
        // File B: records split across two chunks (different thread ordering)
        let chunks_b = vec![
            vec![(MappingType::SingleMapped, vec![hit2])],
            vec![(MappingType::SingleMapped, vec![hit1])],
        ];

        write_test_bulk_rad_file(&path_a, false, &names, &lens, &chunks_a);
        write_test_bulk_rad_file(&path_b, false, &names, &lens, &chunks_b);

        let result = compare_bulk_rad_full(&path_a, &path_b).unwrap();
        assert!(result.passed, "Expected pass: {}", result.notes);
        assert_eq!(result.matching_records, 2);
    }
}
