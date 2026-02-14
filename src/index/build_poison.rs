//! Poison table builder — scans decoy sequences to identify edge k-mers.
//!
//! The "edge method" detects k-mers at boundaries between regions present in
//! the reference index and regions absent from it. When a read maps to such a
//! boundary, the mapping is likely spurious and should be discarded.
//!
//! Port of C++ `poison_table_builder.cpp`.

use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Mutex;

use anyhow::{Context, Result};
use crossbeam::channel;
use sshash_lib::{Kmer, KmerBits};
use tracing::info;

use super::poison_table::{LabeledPoisonOcc, PoisonTable};
use super::reference_index::ReferenceIndex;
use crate::mapping::kmer_value::CanonicalKmer;
use crate::mapping::projected_hits::ProjectedHits;
use crate::mapping::streaming_query::PiscemStreamingQuery;

// ---------------------------------------------------------------------------
// FASTA reading
// ---------------------------------------------------------------------------

/// A single FASTA sequence.
struct FastaRecord {
    seq: Vec<u8>,
}

/// Read FASTA records and send them through a channel.
fn read_fasta_records<R: Read>(reader: R, sender: &channel::Sender<FastaRecord>) -> Result<u64> {
    let buf = BufReader::with_capacity(1 << 20, reader);
    let mut current_seq: Vec<u8> = Vec::new();
    let mut num_seqs: u64 = 0;

    for line in buf.lines() {
        let line = line.context("reading FASTA line")?;
        if line.starts_with('>') {
            if !current_seq.is_empty() {
                num_seqs += 1;
                if sender
                    .send(FastaRecord {
                        seq: std::mem::take(&mut current_seq),
                    })
                    .is_err()
                {
                    break;
                }
            }
        } else {
            current_seq.extend(line.bytes().map(|b| b.to_ascii_uppercase()));
        }
    }

    if !current_seq.is_empty() {
        num_seqs += 1;
        let _ = sender.send(FastaRecord { seq: current_seq });
    }

    Ok(num_seqs)
}

// ---------------------------------------------------------------------------
// Edge-method state machine
// ---------------------------------------------------------------------------

/// State for the poison k-mer detection state machine.
struct PoisonKmerState<'a> {
    first: bool,
    predecessor_present: bool,
    predecessor_canonical: CanonicalKmer,
    predecessor_phits: Option<ProjectedHits<'a>>,
}

impl<'a> PoisonKmerState<'a> {
    fn new() -> Self {
        Self {
            first: true,
            predecessor_present: false,
            predecessor_canonical: CanonicalKmer::new(0),
            predecessor_phits: None,
        }
    }

    fn reset(&mut self) {
        self.first = true;
        self.predecessor_present = false;
        self.predecessor_canonical = CanonicalKmer::new(0);
        self.predecessor_phits = None;
    }

    /// Inspect a k-mer and update state, possibly emitting poison occurrences.
    fn inspect_and_update(
        &mut self,
        canonical: CanonicalKmer,
        phits: Option<ProjectedHits<'a>>,
        k: u32,
        occs: &mut Vec<LabeledPoisonOcc>,
    ) {
        let present = phits.as_ref().is_some_and(|ph| !ph.is_empty());

        if !self.first {
            if self.predecessor_present && !present {
                // Case I: predecessor in index, current not → current is poison.
                if let Some(ref pred_ph) = self.predecessor_phits {
                    let offset: i64 = if pred_ph.hit_fw_on_contig() { 1 } else { -1 };
                    let cpos = pred_ph.contig_pos() as i64 + offset;
                    let max_pos = pred_ph.contig_len() as i64 - k as i64;
                    let cpos = cpos.max(0).min(max_pos) as u32;

                    occs.push(LabeledPoisonOcc {
                        canonical_kmer: canonical,
                        unitig_id: pred_ph.contig_id(),
                        unitig_pos: cpos,
                    });
                }
            } else if !self.predecessor_present && present {
                // Case II: predecessor not in index, current is → predecessor is poison.
                if let Some(ref curr_ph) = phits {
                    let offset: i64 = if curr_ph.hit_fw_on_contig() { -1 } else { 1 };
                    let cpos = curr_ph.contig_pos() as i64 + offset;
                    let max_pos = curr_ph.contig_len() as i64 - k as i64;
                    let cpos = cpos.max(0).min(max_pos) as u32;

                    occs.push(LabeledPoisonOcc {
                        canonical_kmer: self.predecessor_canonical,
                        unitig_id: curr_ph.contig_id(),
                        unitig_pos: cpos,
                    });
                }
            }
        }

        self.first = false;
        self.predecessor_present = present;
        self.predecessor_canonical = canonical;
        self.predecessor_phits = phits;
    }
}

// ---------------------------------------------------------------------------
// Per-sequence scanner
// ---------------------------------------------------------------------------

/// Scan a single decoy sequence for poison k-mers.
fn scan_sequence_for_poison<'a, const K: usize>(
    seq: &[u8],
    index: &'a ReferenceIndex,
    query: &mut PiscemStreamingQuery<'a, K>,
    state: &mut PoisonKmerState<'a>,
    occs: &mut Vec<LabeledPoisonOcc>,
) where
    Kmer<K>: KmerBits,
{
    let k = index.k();
    let k32 = k as u32;

    state.reset();
    query.reset();

    let mask: u64 = if k >= 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };

    let mut fw_kmer: u64 = 0;
    let mut rc_kmer: u64 = 0;
    let mut valid_bases: u32 = 0;

    for (i, &b) in seq.iter().enumerate() {
        let bits = match b {
            b'A' | b'a' => 0u64,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => {
                // Invalid base — reset k-mer rolling state and streaming
                // query (positions won't be consecutive after the gap).
                // Do NOT reset the poison state machine: C++ keeps it
                // running across N gaps, only resetting between sequences.
                valid_bases = 0;
                fw_kmer = 0;
                rc_kmer = 0;
                query.reset();
                continue;
            }
        };

        fw_kmer = ((fw_kmer << 2) | bits) & mask;
        rc_kmer = (rc_kmer >> 2) | ((bits ^ 3) << (2 * (k32 - 1)));
        valid_bases += 1;

        if valid_bases >= k32 {
            let read_pos = i + 1 - k;
            let canonical = CanonicalKmer::new(fw_kmer.min(rc_kmer));

            // Get string slice for streaming query lookup
            let kmer_str = std::str::from_utf8(&seq[read_pos..read_pos + k]).unwrap();
            let result = query.lookup_at(kmer_str, read_pos as i32);
            let phits = index.resolve_lookup(&result);

            state.inspect_and_update(canonical, phits, k32, occs);
        }
    }
}

// ---------------------------------------------------------------------------
// Public builder entry point
// ---------------------------------------------------------------------------

/// Build a poison table by scanning decoy FASTA files.
pub fn build_poison_table<const K: usize>(
    index: &ReferenceIndex,
    decoy_paths: &[String],
    threads: usize,
) -> Result<PoisonTable>
where
    Kmer<K>: KmerBits,
{
    let threads = if threads == 0 {
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
    } else {
        threads
    };

    info!(
        "Building poison table with {} threads from {} decoy file(s)",
        threads,
        decoy_paths.len()
    );

    let total_seqs = AtomicU64::new(0);
    let total_occs_count = AtomicU64::new(0);
    let all_occs: Mutex<Vec<LabeledPoisonOcc>> = Mutex::new(Vec::new());

    let (sender, receiver) = channel::bounded::<FastaRecord>(threads * 4);

    let decoy_paths_ref = decoy_paths;
    let total_seqs_ref = &total_seqs;
    let total_occs_ref = &total_occs_count;
    let all_occs_ref = &all_occs;

    crossbeam::scope(|s| {
        // Producer thread: reads FASTA files and sends records
        let sender_for_producer = sender.clone();
        s.spawn(move |_| {
            for path_str in decoy_paths_ref {
                let path = Path::new(path_str);
                info!("Reading decoy sequences from {}", path.display());
                let (reader, _format) = niffler::send::from_path(path)
                    .unwrap_or_else(|e| {
                        panic!("failed to open {} for decompression: {}", path.display(), e)
                    });
                let n = read_fasta_records(reader, &sender_for_producer)
                    .unwrap_or_else(|e| panic!("error reading {}: {}", path.display(), e));
                total_seqs_ref.fetch_add(n, Ordering::Relaxed);
            }
            drop(sender_for_producer);
        });

        // Drop original sender so channel closes when producer finishes
        drop(sender);

        // Worker threads
        let recv_ref = &receiver;
        for _ in 0..threads {
            s.spawn(move |_| {
                let recv = recv_ref.clone();
                let mut query = PiscemStreamingQuery::<K>::new(index.dict());
                let mut state = PoisonKmerState::new();
                let mut local_occs: Vec<LabeledPoisonOcc> = Vec::new();

                for record in recv {
                    scan_sequence_for_poison::<K>(
                        &record.seq,
                        index,
                        &mut query,
                        &mut state,
                        &mut local_occs,
                    );
                }

                total_occs_ref.fetch_add(local_occs.len() as u64, Ordering::Relaxed);
                all_occs_ref.lock().unwrap().append(&mut local_occs);
            });
        }
    })
    .map_err(|e| anyhow::anyhow!("thread panicked: {:?}", e))?;

    let num_seqs = total_seqs.load(Ordering::Relaxed);
    let raw_occs = total_occs_count.load(Ordering::Relaxed);
    info!(
        "Scanned {} decoy sequences, found {} raw poison occurrences",
        num_seqs, raw_occs,
    );

    let occs = all_occs.into_inner().unwrap();
    let table = PoisonTable::build_from_occs(occs)?;

    info!(
        "Poison table: {} distinct k-mers, {} occurrences, max occ = {}",
        table.num_poison_kmers(),
        table.num_poison_occs(),
        table.max_poison_occ(),
    );

    Ok(table)
}

/// Verify that no poison k-mers are present in the reference dictionary.
///
/// Returns the number of poison k-mers found in the dictionary (should be 0).
pub fn verify_poison_table<const K: usize>(
    table: &PoisonTable,
    index: &ReferenceIndex,
) -> usize
where
    Kmer<K>: KmerBits,
{
    let k = index.k();
    let mut found_in_dict = 0usize;
    let mut checked = 0usize;

    for &kmer in table.keys() {
        // Convert packed canonical k-mer to a string
        let raw = kmer.as_u64();
        let mut s = String::with_capacity(k);
        for i in (0..k).rev() {
            let bits = (raw >> (2 * i)) & 3;
            s.push(match bits {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                3 => 'T',
                _ => unreachable!(),
            });
        }

        // Look up in dictionary via streaming query (single lookup, reset each time)
        let mut query = PiscemStreamingQuery::<K>::new(index.dict());
        let result = query.lookup_at(&s, 0);
        let phits = index.resolve_lookup(&result);
        let present = phits.as_ref().is_some_and(|ph| !ph.is_empty());
        if present {
            found_in_dict += 1;
            if found_in_dict <= 10 {
                info!(
                    "WARNING: poison k-mer {} (canonical={}) FOUND in dictionary!",
                    s, raw,
                );
            }
        }
        checked += 1;
        if checked % 500_000 == 0 {
            info!("Verified {}/{} poison k-mers ({} in dict so far)",
                checked, table.num_poison_kmers(), found_in_dict);
        }
    }

    info!(
        "Poison verification: {}/{} k-mers found in dictionary",
        found_in_dict, table.num_poison_kmers()
    );
    found_in_dict
}
