/// Detailed diagnostic for ATAC RAD mismatches.
///
/// Reads both C++ and Rust RAD files, identifies all mismatched records,
/// and categorizes them to help identify root causes.
use std::collections::HashMap;
use std::path::Path;

use piscem_rs::verify::rad_compare::{read_atac_rad_records, CanonicalAtacRecord};

fn main() {
    let cpp_rad = Path::new("/tmp/cpp_atac_out/map.rad");
    let rust_rad = Path::new("/tmp/rust_atac_out2/map.rad");

    let (header_a, records_a) = read_atac_rad_records(cpp_rad).expect("read C++");
    let (header_b, records_b) = read_atac_rad_records(rust_rad).expect("read Rust");

    // Build ref ID translation: Rust → C++
    let name_to_id_a: HashMap<&str, u32> = header_a
        .ref_names
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i as u32))
        .collect();
    let b_to_a_id: Vec<u32> = header_b
        .ref_names
        .iter()
        .map(|name| name_to_id_a.get(name.as_str()).copied().unwrap_or(u32::MAX))
        .collect();

    // Translate Rust records to C++ ref namespace
    let records_b_translated: Vec<CanonicalAtacRecord> = records_b
        .iter()
        .map(|rec| {
            let mut t = rec.clone();
            for aln in &mut t.alignments {
                if (aln.0 as usize) < b_to_a_id.len() {
                    aln.0 = b_to_a_id[aln.0 as usize];
                }
            }
            t.alignments.sort();
            t
        })
        .collect();

    println!("Records: C++={}, Rust={}", records_a.len(), records_b_translated.len());

    // Build frequency maps
    let mut freq_a: HashMap<CanonicalAtacRecord, u64> = HashMap::with_capacity(records_a.len());
    for rec in &records_a {
        *freq_a.entry(rec.clone()).or_insert(0) += 1;
    }
    let mut freq_b: HashMap<CanonicalAtacRecord, u64> =
        HashMap::with_capacity(records_b_translated.len());
    for rec in &records_b_translated {
        *freq_b.entry(rec.clone()).or_insert(0) += 1;
    }

    // Collect records only in C++ (excess in A)
    let mut only_a: Vec<&CanonicalAtacRecord> = Vec::new();
    for (rec, &count_a) in &freq_a {
        let count_b = freq_b.get(rec).copied().unwrap_or(0);
        if count_a > count_b {
            for _ in 0..(count_a - count_b) {
                only_a.push(rec);
            }
        }
    }
    // Collect records only in Rust (excess in B)
    let mut only_b: Vec<&CanonicalAtacRecord> = Vec::new();
    for (rec, &count_b) in &freq_b {
        let count_a = freq_a.get(rec).copied().unwrap_or(0);
        if count_b > count_a {
            for _ in 0..(count_b - count_a) {
                only_b.push(rec);
            }
        }
    }

    println!("Only in C++: {}", only_a.len());
    println!("Only in Rust: {}", only_b.len());

    // Index by barcode for pairing
    let mut a_by_bc: HashMap<u64, Vec<&CanonicalAtacRecord>> = HashMap::new();
    for rec in &only_a {
        a_by_bc.entry(rec.bc).or_default().push(rec);
    }
    let mut b_by_bc: HashMap<u64, Vec<&CanonicalAtacRecord>> = HashMap::new();
    for rec in &only_b {
        b_by_bc.entry(rec.bc).or_default().push(rec);
    }

    // Count unique barcodes that appear in mismatches
    let all_bc: std::collections::HashSet<u64> = a_by_bc
        .keys()
        .chain(b_by_bc.keys())
        .copied()
        .collect();
    println!("Unique barcodes with mismatches: {}", all_bc.len());

    // Find paired mismatches (same BC appears in both only_a and only_b)
    let mut paired_bc = 0usize;
    let mut unpaired_a_only = 0usize;
    let mut unpaired_b_only = 0usize;

    // Categorize paired mismatches
    let mut cat_same_targets_same_nalns = 0usize;
    let mut cat_same_targets_diff_nalns = 0usize;
    let mut cat_diff_targets = 0usize;
    let mut cat_subset_targets = 0usize; // one is a subset of the other

    // Detail categories for same-target mismatches
    let mut detail_pos_diff_only = 0usize;
    let mut detail_flen_diff_only = 0usize;
    let mut detail_type_diff_only = 0usize;
    let mut detail_pos_and_flen = 0usize;
    let mut detail_other = 0usize;

    // Position difference histogram
    let mut pos_diffs: Vec<i64> = Vec::new();
    let mut flen_diffs: Vec<i64> = Vec::new();

    // Number of alignments histogram
    let mut nalns_a_hist: HashMap<usize, usize> = HashMap::new();
    let mut nalns_b_hist: HashMap<usize, usize> = HashMap::new();

    let mut examples: Vec<String> = Vec::new();
    let max_examples = 30;

    for &bc in &all_bc {
        let a_recs = a_by_bc.get(&bc);
        let b_recs = b_by_bc.get(&bc);

        match (a_recs, b_recs) {
            (Some(a_list), Some(b_list)) => {
                paired_bc += 1;

                // For each record in A, try to find best match in B
                for a_rec in a_list {
                    *nalns_a_hist.entry(a_rec.alignments.len()).or_insert(0) += 1;

                    // Find a B record with closest match
                    let mut best_match: Option<&&CanonicalAtacRecord> = None;
                    let mut best_score = 0usize;

                    let a_tids: std::collections::HashSet<u32> =
                        a_rec.alignments.iter().map(|a| a.0).collect();

                    for b_rec in b_list {
                        let b_tids: std::collections::HashSet<u32> =
                            b_rec.alignments.iter().map(|a| a.0).collect();
                        let overlap = a_tids.intersection(&b_tids).count();
                        if overlap > best_score {
                            best_score = overlap;
                            best_match = Some(b_rec);
                        }
                    }

                    if let Some(b_rec) = best_match {
                        *nalns_b_hist.entry(b_rec.alignments.len()).or_insert(0) += 1;
                        let b_tids: std::collections::HashSet<u32> =
                            b_rec.alignments.iter().map(|a| a.0).collect();

                        if a_tids == b_tids {
                            // Same target set
                            if a_rec.alignments.len() == b_rec.alignments.len() {
                                cat_same_targets_same_nalns += 1;

                                // Classify detail differences
                                let mut pos_differs = false;
                                let mut flen_differs = false;
                                let mut type_differs = false;

                                for (a_aln, b_aln) in
                                    a_rec.alignments.iter().zip(b_rec.alignments.iter())
                                {
                                    if a_aln.2 != b_aln.2 {
                                        pos_differs = true;
                                        pos_diffs
                                            .push(a_aln.2 as i64 - b_aln.2 as i64);
                                    }
                                    if a_aln.3 != b_aln.3 {
                                        flen_differs = true;
                                        flen_diffs
                                            .push(a_aln.3 as i64 - b_aln.3 as i64);
                                    }
                                    if a_aln.1 != b_aln.1 {
                                        type_differs = true;
                                    }
                                }

                                match (pos_differs, flen_differs, type_differs) {
                                    (true, false, false) => detail_pos_diff_only += 1,
                                    (false, true, false) => detail_flen_diff_only += 1,
                                    (false, false, true) => detail_type_diff_only += 1,
                                    (true, true, false) => detail_pos_and_flen += 1,
                                    _ => detail_other += 1,
                                }

                                if examples.len() < max_examples {
                                    examples.push(format!(
                                        "SAME_TARGETS bc={} nalns={}\n  C++:  {}\n  Rust: {}",
                                        bc,
                                        a_rec.alignments.len(),
                                        a_rec,
                                        b_rec
                                    ));
                                }
                            } else {
                                cat_same_targets_diff_nalns += 1;
                                if examples.len() < max_examples {
                                    examples.push(format!(
                                        "SAME_TARGETS_DIFF_NALNS bc={} a_nalns={} b_nalns={}\n  C++:  {}\n  Rust: {}",
                                        bc,
                                        a_rec.alignments.len(),
                                        b_rec.alignments.len(),
                                        a_rec,
                                        b_rec
                                    ));
                                }
                            }
                        } else {
                            // Different target sets
                            let a_is_subset = a_tids.is_subset(&b_tids);
                            let b_is_subset = b_tids.is_subset(&a_tids);
                            if a_is_subset || b_is_subset {
                                cat_subset_targets += 1;
                                if examples.len() < max_examples {
                                    let dir = if a_is_subset {
                                        "C++ ⊂ Rust"
                                    } else {
                                        "Rust ⊂ C++"
                                    };
                                    examples.push(format!(
                                        "SUBSET_TARGETS ({dir}) bc={}\n  C++:  {}\n  Rust: {}",
                                        bc, a_rec, b_rec
                                    ));
                                }
                            } else {
                                cat_diff_targets += 1;
                                if examples.len() < max_examples {
                                    examples.push(format!(
                                        "DIFF_TARGETS bc={}\n  C++:  {}\n  Rust: {}",
                                        bc, a_rec, b_rec
                                    ));
                                }
                            }
                        }
                    }
                }
            }
            (Some(a_list), None) => {
                unpaired_a_only += a_list.len();
                for a_rec in a_list {
                    *nalns_a_hist.entry(a_rec.alignments.len()).or_insert(0) += 1;
                }
                if examples.len() < max_examples && unpaired_a_only <= 3 {
                    for a_rec in a_list {
                        examples.push(format!("ONLY_CPP bc={}\n  C++: {}", bc, a_rec));
                    }
                }
            }
            (None, Some(b_list)) => {
                unpaired_b_only += b_list.len();
                for b_rec in b_list {
                    *nalns_b_hist.entry(b_rec.alignments.len()).or_insert(0) += 1;
                }
                if examples.len() < max_examples && unpaired_b_only <= 3 {
                    for b_rec in b_list {
                        examples.push(format!("ONLY_RUST bc={}\n  Rust: {}", bc, b_rec));
                    }
                }
            }
            _ => {}
        }
    }

    println!("\n=== Mismatch Categories ===");
    println!("Paired (same BC in both):     {}", paired_bc);
    println!("  Same targets, same #alns:   {}", cat_same_targets_same_nalns);
    println!("    pos diff only:            {}", detail_pos_diff_only);
    println!("    flen diff only:           {}", detail_flen_diff_only);
    println!("    type diff only:           {}", detail_type_diff_only);
    println!("    pos+flen diff:            {}", detail_pos_and_flen);
    println!("    other:                    {}", detail_other);
    println!("  Same targets, diff #alns:   {}", cat_same_targets_diff_nalns);
    println!("  Subset targets:             {}", cat_subset_targets);
    println!("  Different targets:          {}", cat_diff_targets);
    println!("Unpaired (C++ only):          {}", unpaired_a_only);
    println!("Unpaired (Rust only):         {}", unpaired_b_only);

    if !pos_diffs.is_empty() {
        pos_diffs.sort();
        let min = pos_diffs[0];
        let max = pos_diffs[pos_diffs.len() - 1];
        let median = pos_diffs[pos_diffs.len() / 2];
        println!("\nPosition diffs (C++ - Rust): min={min}, max={max}, median={median}, n={}", pos_diffs.len());
    }

    if !flen_diffs.is_empty() {
        flen_diffs.sort();
        let min = flen_diffs[0];
        let max = flen_diffs[flen_diffs.len() - 1];
        let median = flen_diffs[flen_diffs.len() / 2];
        println!("Frag len diffs (C++ - Rust): min={min}, max={max}, median={median}, n={}", flen_diffs.len());
    }

    println!("\n#alignments histogram (C++ mismatches):");
    let mut hist_keys: Vec<usize> = nalns_a_hist.keys().copied().collect();
    hist_keys.sort();
    for k in &hist_keys {
        println!("  {} alns: {} records", k, nalns_a_hist[k]);
    }

    println!("#alignments histogram (Rust mismatches):");
    let mut hist_keys: Vec<usize> = nalns_b_hist.keys().copied().collect();
    hist_keys.sort();
    for k in &hist_keys {
        println!("  {} alns: {} records", k, nalns_b_hist[k]);
    }

    println!("\n=== Examples ===");
    for ex in &examples {
        println!("{}\n", ex);
    }
}
