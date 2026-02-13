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
