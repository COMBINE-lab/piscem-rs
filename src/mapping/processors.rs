//! Parallel processors for paraseq-based mapping pipelines.
//!
//! Each processor implements the appropriate paraseq parallel processing trait
//! (`ParallelProcessor`, `PairedParallelProcessor`, or `MultiParallelProcessor`)
//! and encapsulates all per-thread mapping state.
//!
//! Per-thread mutable state is stored in `Option<ThreadState>` and initialized
//! lazily on first batch. The custom `Clone` impl sets `state: None` so each
//! cloned processor (one per worker thread) gets fresh state.

use std::sync::atomic::Ordering;
use std::sync::Mutex;

use indicatif::ProgressBar;
use paraseq::parallel::{MultiParallelProcessor, PairedParallelProcessor, ParallelProcessor};
use paraseq::Record;
use smallvec::SmallVec;
use sshash_lib::{Kmer, KmerBits};

use crate::index::reference_index::ReferenceIndex;
use crate::io::rad::{
    RadWriter, pack_bases_2bit, write_atac_record, write_bulk_record, write_sc_record,
};
use crate::io::threads::{MappingStats, OutputInfo};
use crate::mapping::binning::BinPos;
use crate::mapping::cache::MappingCache;
use crate::mapping::filters::PoisonState;
use crate::mapping::hit_searcher::{HitSearcher, SkippingStrategy};
use crate::mapping::hits::MappingType;
use crate::mapping::map_fragment::{
    map_pe_fragment, map_pe_fragment_atac, map_se_fragment, map_se_fragment_atac,
};
use crate::mapping::merge_pairs::{remove_duplicate_hits_pub, simple_hit_cmp_bins};
use crate::mapping::overlap::{find_overlap, OverlapType};
use crate::mapping::protocols::scrna::{barcode_has_n, count_ns, is_all_acgt, recover_barcode};
use crate::mapping::protocols::Protocol;
use crate::mapping::sketch_hit_simple::SketchHitInfoSimple;
use crate::mapping::streaming_query::PiscemStreamingQuery;
use crate::mapping::unitig_end_cache::UnitigEndCache;

// ===========================================================================
// Constants
// ===========================================================================

/// Maximum mapped reads per RAD chunk before flushing (matches C++ max_chunk_reads).
const MAX_CHUNK_READS: u32 = 5000;

// ===========================================================================
// Common thread state helpers
// ===========================================================================

/// Common per-thread mapping state shared across all processor types.
struct CommonThreadState<'a, const K: usize>
where
    Kmer<K>: KmerBits,
{
    hs: HitSearcher<'a>,
    query: PiscemStreamingQuery<'a, K>,
    cache_out: MappingCache<SketchHitInfoSimple>,
    cache_left: MappingCache<SketchHitInfoSimple>,
    cache_right: MappingCache<SketchHitInfoSimple>,
    poison_state: PoisonState<'a>,
    rad_writer: RadWriter,
    local_reads: u64,
    local_mapped: u64,
    local_poisoned: u64,
    num_reads_in_chunk: u32,
    chunk_in_progress: bool,
}

impl<'a, const K: usize> CommonThreadState<'a, K>
where
    Kmer<K>: KmerBits,
{
    fn new(index: &'a ReferenceIndex, end_cache: &'a UnitigEndCache) -> Self {
        Self {
            hs: HitSearcher::new(index),
            query: PiscemStreamingQuery::<K>::with_cache(index.dict(), end_cache),
            cache_out: MappingCache::new(K),
            cache_left: MappingCache::new(K),
            cache_right: MappingCache::new(K),
            poison_state: PoisonState::new(index.poison_table()),
            rad_writer: RadWriter::with_capacity(150_000),
            local_reads: 0,
            local_mapped: 0,
            local_poisoned: 0,
            num_reads_in_chunk: 0,
            chunk_in_progress: false,
        }
    }

    /// Ensure a RAD chunk header has been written. Call before processing records.
    fn ensure_chunk_started(&mut self) {
        if !self.chunk_in_progress {
            self.rad_writer.clear();
            self.rad_writer.write_u32(0); // placeholder num_bytes
            self.rad_writer.write_u32(0); // placeholder num_reads
            self.num_reads_in_chunk = 0;
            self.chunk_in_progress = true;
        }
    }

    /// Finalize and flush the current chunk if it has any reads.
    fn finalize_chunk(&mut self, output: &OutputInfo) {
        if !self.chunk_in_progress {
            return;
        }
        let total_bytes = self.rad_writer.len() as u32;
        self.rad_writer.write_u32_at_offset(0, total_bytes);
        self.rad_writer.write_u32_at_offset(4, self.num_reads_in_chunk);
        if self.num_reads_in_chunk > 0 {
            let mut file = output.rad_file.lock().unwrap();
            self.rad_writer.flush_to(&mut *file).ok();
            output.num_chunks.fetch_add(1, Ordering::Relaxed);
        }
        self.chunk_in_progress = false;
    }

    /// Flush if we've accumulated enough mapped reads.
    fn maybe_flush_chunk(&mut self, output: &OutputInfo) {
        if self.num_reads_in_chunk >= MAX_CHUNK_READS {
            self.finalize_chunk(output);
        }
    }

    fn flush_stats(&self, stats: &MappingStats) {
        stats
            .num_reads
            .fetch_add(self.local_reads, Ordering::Relaxed);
        stats
            .num_mapped
            .fetch_add(self.local_mapped, Ordering::Relaxed);
        stats
            .num_poisoned
            .fetch_add(self.local_poisoned, Ordering::Relaxed);
    }
}

// ===========================================================================
// BulkProcessor
// ===========================================================================

/// Parallel processor for bulk RNA-seq mapping.
///
/// Implements `PairedParallelProcessor` for PE and `ParallelProcessor` for SE.
pub struct BulkProcessor<'a, const K: usize>
where
    Kmer<K>: KmerBits,
{
    index: &'a ReferenceIndex,
    end_cache: &'a UnitigEndCache,
    output: &'a OutputInfo,
    stats: &'a MappingStats,
    progress: &'a ProgressBar,
    strat: SkippingStrategy,
    state: Option<CommonThreadState<'a, K>>,
}

impl<'a, const K: usize> BulkProcessor<'a, K>
where
    Kmer<K>: KmerBits,
{
    pub fn new(
        index: &'a ReferenceIndex,
        end_cache: &'a UnitigEndCache,
        output: &'a OutputInfo,
        stats: &'a MappingStats,
        strat: SkippingStrategy,
        progress: &'a ProgressBar,
    ) -> Self {
        Self {
            index,
            end_cache,
            output,
            stats,
            progress,
            strat,
            state: None,
        }
    }
}

impl<const K: usize> Clone for BulkProcessor<'_, K>
where
    Kmer<K>: KmerBits,
{
    fn clone(&self) -> Self {
        Self {
            index: self.index,
            end_cache: self.end_cache,
            output: self.output,
            stats: self.stats,
            progress: self.progress,
            strat: self.strat,
            state: None,
        }
    }
}

// Safety: all fields are either `Copy` references or `None`.
unsafe impl<const K: usize> Send for BulkProcessor<'_, K> where Kmer<K>: KmerBits {}

// --- Bulk PE ---

impl<'a, 'r, const K: usize> PairedParallelProcessor<paraseq::fastq::RefRecord<'r>>
    for BulkProcessor<'a, K>
where
    Kmer<K>: KmerBits,
{
    fn process_record_pair_batch(
        &mut self,
        record_pairs: impl Iterator<Item = (paraseq::fastq::RefRecord<'r>, paraseq::fastq::RefRecord<'r>)>,
    ) -> paraseq::Result<()> {
        let index = self.index;
        let end_cache = self.end_cache;
        let strat = self.strat;
        let s = self
            .state
            .get_or_insert_with(|| CommonThreadState::new(index, end_cache));
        s.ensure_chunk_started();

        let mut batch_reads: u64 = 0;
        for (rec1, rec2) in record_pairs {
            s.local_reads += 1;
            batch_reads += 1;
            let seq1 = rec1.seq();
            let seq2 = rec2.seq();

            s.poison_state.paired_for_mapping = true;
            map_pe_fragment::<K>(
                &seq1,
                &seq2,
                &mut s.hs,
                &mut s.query,
                &mut s.cache_left,
                &mut s.cache_right,
                &mut s.cache_out,
                index,
                &mut s.poison_state,
                strat,
            );

            if s.poison_state.is_poisoned() {
                s.local_poisoned += 1;
                continue;
            }

            if s.cache_out.map_type != MappingType::Unmapped {
                s.local_mapped += 1;
                write_bulk_record(
                    s.cache_out.map_type,
                    &s.cache_out.accepted_hits,
                    &mut s.rad_writer,
                );
                s.num_reads_in_chunk += 1;
            }
        }
        self.progress.inc(batch_reads);
        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::Result<()> {
        let output = self.output;
        if let Some(s) = &mut self.state {
            s.maybe_flush_chunk(output);
        }
        Ok(())
    }

    fn on_thread_complete(&mut self) -> paraseq::Result<()> {
        let output = self.output;
        let stats = self.stats;
        if let Some(s) = &mut self.state {
            s.finalize_chunk(output);
            s.flush_stats(stats);
        }
        Ok(())
    }
}

// --- Bulk SE ---

impl<'a, 'r, const K: usize> ParallelProcessor<paraseq::fastq::RefRecord<'r>>
    for BulkProcessor<'a, K>
where
    Kmer<K>: KmerBits,
{
    fn process_record_batch(
        &mut self,
        records: impl Iterator<Item = paraseq::fastq::RefRecord<'r>>,
    ) -> paraseq::Result<()> {
        let index = self.index;
        let end_cache = self.end_cache;
        let strat = self.strat;
        let s = self
            .state
            .get_or_insert_with(|| CommonThreadState::new(index, end_cache));
        s.ensure_chunk_started();

        let mut batch_reads: u64 = 0;
        for rec in records {
            s.local_reads += 1;
            batch_reads += 1;
            let seq1 = rec.seq();

            s.poison_state.paired_for_mapping = false;
            map_se_fragment::<K>(
                &seq1,
                &mut s.hs,
                &mut s.query,
                &mut s.cache_out,
                index,
                &mut s.poison_state,
                strat,
            );

            if s.poison_state.is_poisoned() {
                s.local_poisoned += 1;
                continue;
            }

            if s.cache_out.map_type != MappingType::Unmapped {
                s.local_mapped += 1;
                write_bulk_record(
                    s.cache_out.map_type,
                    &s.cache_out.accepted_hits,
                    &mut s.rad_writer,
                );
                s.num_reads_in_chunk += 1;
            }
        }
        self.progress.inc(batch_reads);
        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::Result<()> {
        let output = self.output;
        if let Some(s) = &mut self.state {
            s.maybe_flush_chunk(output);
        }
        Ok(())
    }

    fn on_thread_complete(&mut self) -> paraseq::Result<()> {
        let output = self.output;
        let stats = self.stats;
        if let Some(s) = &mut self.state {
            s.finalize_chunk(output);
            s.flush_stats(stats);
        }
        Ok(())
    }
}

// ===========================================================================
// ScrnaProcessor
// ===========================================================================

/// Per-thread state for scRNA mapping (extends common state).
struct ScrnaThreadState<'a, const K: usize>
where
    Kmer<K>: KmerBits,
{
    common: CommonThreadState<'a, K>,
    local_rlen_samples: Vec<u32>,
}

/// Parallel processor for scRNA-seq mapping.
pub struct ScrnaProcessor<'a, const K: usize>
where
    Kmer<K>: KmerBits,
{
    index: &'a ReferenceIndex,
    end_cache: &'a UnitigEndCache,
    output: &'a OutputInfo,
    stats: &'a MappingStats,
    progress: &'a ProgressBar,
    strat: SkippingStrategy,
    protocol: &'a dyn Protocol,
    bc_len: u16,
    umi_len: u16,
    with_position: bool,
    read_length_samples: &'a Mutex<Vec<u32>>,
    state: Option<ScrnaThreadState<'a, K>>,
}

impl<'a, const K: usize> ScrnaProcessor<'a, K>
where
    Kmer<K>: KmerBits,
{
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        index: &'a ReferenceIndex,
        end_cache: &'a UnitigEndCache,
        output: &'a OutputInfo,
        stats: &'a MappingStats,
        strat: SkippingStrategy,
        protocol: &'a dyn Protocol,
        bc_len: u16,
        umi_len: u16,
        with_position: bool,
        read_length_samples: &'a Mutex<Vec<u32>>,
        progress: &'a ProgressBar,
    ) -> Self {
        Self {
            index,
            end_cache,
            output,
            stats,
            progress,
            strat,
            protocol,
            bc_len,
            umi_len,
            with_position,
            read_length_samples,
            state: None,
        }
    }
}

impl<const K: usize> Clone for ScrnaProcessor<'_, K>
where
    Kmer<K>: KmerBits,
{
    fn clone(&self) -> Self {
        Self {
            index: self.index,
            end_cache: self.end_cache,
            output: self.output,
            stats: self.stats,
            progress: self.progress,
            strat: self.strat,
            protocol: self.protocol,
            bc_len: self.bc_len,
            umi_len: self.umi_len,
            with_position: self.with_position,
            read_length_samples: self.read_length_samples,
            state: None,
        }
    }
}

unsafe impl<const K: usize> Send for ScrnaProcessor<'_, K> where Kmer<K>: KmerBits {}

impl<'a, 'r, const K: usize> PairedParallelProcessor<paraseq::fastq::RefRecord<'r>>
    for ScrnaProcessor<'a, K>
where
    Kmer<K>: KmerBits,
{
    fn process_record_pair_batch(
        &mut self,
        record_pairs: impl Iterator<Item = (paraseq::fastq::RefRecord<'r>, paraseq::fastq::RefRecord<'r>)>,
    ) -> paraseq::Result<()> {
        let index = self.index;
        let end_cache = self.end_cache;
        let strat = self.strat;
        let protocol = self.protocol;
        let bc_len = self.bc_len;
        let umi_len = self.umi_len;
        let with_position = self.with_position;
        let is_bio_paired = protocol.is_bio_paired_end();
        let max_rlen_samples: usize = 10;

        let st = self.state.get_or_insert_with(|| ScrnaThreadState {
            common: CommonThreadState::new(index, end_cache),
            local_rlen_samples: Vec::new(),
        });
        let s = &mut st.common;
        s.ensure_chunk_started();

        let mut batch_reads: u64 = 0;
        for (rec1, rec2) in record_pairs {
            s.local_reads += 1;
            batch_reads += 1;

            let r1 = rec1.seq();
            let r2 = rec2.seq();

            // Extract technical sequences
            // C++ returns nullptr if R1 too short for BC or UMI → skip read
            let tech = protocol.extract_tech_seqs(&r1, &r2);
            let bc_raw = match tech.barcode {
                Some(bc) if !bc.is_empty() => bc,
                _ => continue,
            };
            let umi_raw = match tech.umi {
                Some(umi) if !umi.is_empty() => umi,
                _ => continue,
            };

            // Barcode validation + recovery (matching C++ recover_barcode + fromChars)
            // C++ uses find_first_not_of("ACTGactg") which catches any non-ACGT char
            let n_count = count_ns(bc_raw);
            if n_count > 1 {
                continue;
            }
            let recovered_bc = if n_count == 1 {
                recover_barcode(bc_raw)
            } else {
                None
            };
            let bc_to_pack = match &recovered_bc {
                Some(bc) => bc.as_slice(),
                None => {
                    // No recovery needed/attempted — validate original BC
                    if !is_all_acgt(bc_raw) {
                        continue;
                    }
                    bc_raw
                }
            };
            // Validate recovered BC (C++ fromChars check)
            if !is_all_acgt(bc_to_pack) {
                continue;
            }

            // UMI validation (matching C++ umi_kmer.fromChars check)
            if !is_all_acgt(umi_raw) {
                continue;
            }

            let bc_packed = pack_bases_2bit(bc_to_pack);
            let umi_packed = pack_bases_2bit(umi_raw);

            // Extract mappable reads
            let alignable = protocol.extract_mappable_reads(&r1, &r2);

            if is_bio_paired {
                let seq1 = alignable.seq1.unwrap_or(&[]);
                let seq2 = alignable.seq2.unwrap_or(&[]);
                if seq1.is_empty() && seq2.is_empty() {
                    continue;
                }
                s.poison_state.paired_for_mapping = true;
                map_pe_fragment::<K>(
                    seq1,
                    seq2,
                    &mut s.hs,
                    &mut s.query,
                    &mut s.cache_left,
                    &mut s.cache_right,
                    &mut s.cache_out,
                    index,
                    &mut s.poison_state,
                    strat,
                );

                if with_position && st.local_rlen_samples.len() < max_rlen_samples {
                    st.local_rlen_samples.push(seq2.len() as u32);
                }
            } else {
                let seq1 = alignable.seq1.unwrap_or(&[]);
                if seq1.is_empty() {
                    continue;
                }
                s.poison_state.paired_for_mapping = false;
                map_se_fragment::<K>(
                    seq1,
                    &mut s.hs,
                    &mut s.query,
                    &mut s.cache_out,
                    index,
                    &mut s.poison_state,
                    strat,
                );

                if with_position && st.local_rlen_samples.len() < max_rlen_samples {
                    st.local_rlen_samples.push(seq1.len() as u32);
                }
            }

            if s.poison_state.is_poisoned() {
                s.local_poisoned += 1;
                continue;
            }

            if s.cache_out.map_type != MappingType::Unmapped {
                s.local_mapped += 1;
                write_sc_record(
                    bc_packed,
                    umi_packed,
                    bc_len,
                    umi_len,
                    s.cache_out.map_type,
                    &s.cache_out.accepted_hits,
                    with_position,
                    &mut s.rad_writer,
                );
                s.num_reads_in_chunk += 1;
            }
        }
        self.progress.inc(batch_reads);
        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::Result<()> {
        let output = self.output;
        if let Some(st) = &mut self.state {
            st.common.maybe_flush_chunk(output);
        }
        Ok(())
    }

    fn on_thread_complete(&mut self) -> paraseq::Result<()> {
        let output = self.output;
        let stats = self.stats;
        let with_position = self.with_position;
        let read_length_samples = self.read_length_samples;
        if let Some(st) = &mut self.state {
            st.common.finalize_chunk(output);
            st.common.flush_stats(stats);
            if with_position && !st.local_rlen_samples.is_empty() {
                let mut samples = read_length_samples.lock().unwrap();
                samples.extend_from_slice(&st.local_rlen_samples);
            }
        }
        Ok(())
    }
}

// ===========================================================================
// ScatacProcessor
// ===========================================================================

/// Parallel processor for scATAC-seq mapping (triple-file input).
pub struct ScatacProcessor<'a, const K: usize>
where
    Kmer<K>: KmerBits,
{
    index: &'a ReferenceIndex,
    end_cache: &'a UnitigEndCache,
    output: &'a OutputInfo,
    stats: &'a MappingStats,
    progress: &'a ProgressBar,
    binning: &'a BinPos,
    bc_len: u16,
    tn5_shift: bool,
    min_overlap: i32,
    state: Option<CommonThreadState<'a, K>>,
}

impl<'a, const K: usize> ScatacProcessor<'a, K>
where
    Kmer<K>: KmerBits,
{
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        index: &'a ReferenceIndex,
        end_cache: &'a UnitigEndCache,
        output: &'a OutputInfo,
        stats: &'a MappingStats,
        binning: &'a BinPos,
        bc_len: u16,
        tn5_shift: bool,
        min_overlap: i32,
        progress: &'a ProgressBar,
    ) -> Self {
        Self {
            index,
            end_cache,
            output,
            stats,
            progress,
            binning,
            bc_len,
            tn5_shift,
            min_overlap,
            state: None,
        }
    }
}

impl<const K: usize> Clone for ScatacProcessor<'_, K>
where
    Kmer<K>: KmerBits,
{
    fn clone(&self) -> Self {
        Self {
            index: self.index,
            end_cache: self.end_cache,
            output: self.output,
            stats: self.stats,
            progress: self.progress,
            binning: self.binning,
            bc_len: self.bc_len,
            tn5_shift: self.tn5_shift,
            min_overlap: self.min_overlap,
            state: None,
        }
    }
}

unsafe impl<const K: usize> Send for ScatacProcessor<'_, K> where Kmer<K>: KmerBits {}

impl<'a, 'r, const K: usize> MultiParallelProcessor<paraseq::fastq::RefRecord<'r>>
    for ScatacProcessor<'a, K>
where
    Kmer<K>: KmerBits,
{
    fn process_multi_record_batch(
        &mut self,
        multi_records: impl Iterator<
            Item = SmallVec<[paraseq::fastq::RefRecord<'r>; paraseq::MAX_ARITY]>,
        >,
    ) -> paraseq::Result<()> {
        let index = self.index;
        let end_cache = self.end_cache;
        let binning = self.binning;
        let bc_len = self.bc_len;
        let tn5_shift = self.tn5_shift;
        let min_overlap = self.min_overlap;

        let s = self
            .state
            .get_or_insert_with(|| CommonThreadState::new(index, end_cache));
        s.ensure_chunk_started();

        let mut batch_reads: u64 = 0;
        for multi in multi_records {
            s.local_reads += 1;
            batch_reads += 1;

            // multi[0] = R1 (genomic left), multi[1] = barcode, multi[2] = R2 (genomic right)
            if multi.len() < 3 {
                continue;
            }
            let r1 = multi[0].seq();
            let bc_raw = multi[1].seq();
            let r2 = multi[2].seq();

            if r1.is_empty() || r2.is_empty() {
                continue;
            }

            // Check barcode for Ns (truncate to bc_len if longer)
            let bc_end = (bc_len as usize).min(bc_raw.len());
            let bc_slice = &bc_raw[..bc_end];

            let n_count = count_ns(bc_slice);
            if n_count > 1 {
                continue;
            }
            let recovered_bc = if n_count == 1 {
                recover_barcode(bc_slice)
            } else {
                None
            };
            let bc_to_pack = match &recovered_bc {
                Some(bc) => bc.as_slice(),
                None => {
                    if barcode_has_n(bc_slice) {
                        continue;
                    }
                    bc_slice
                }
            };
            let bc_packed = pack_bases_2bit(bc_to_pack);

            s.poison_state.clear();
            s.cache_out.clear();

            // Try mate overlap detection
            let mate_ov = find_overlap(&r1, &r2, min_overlap, 0);

            if mate_ov.ov_type != OverlapType::NoOverlap && !mate_ov.frag.is_empty() {
                // Map merged fragment as single read (bin-based)
                s.poison_state.paired_for_mapping = false;
                map_se_fragment_atac::<K>(
                    &mate_ov.frag,
                    &mut s.hs,
                    &mut s.query,
                    &mut s.cache_out,
                    index,
                    &mut s.poison_state,
                    binning,
                );

                // Sort+dedup hits
                if !s.cache_out.accepted_hits.is_empty() {
                    s.cache_out.accepted_hits.sort_by(simple_hit_cmp_bins);
                    remove_duplicate_hits_pub(&mut s.cache_out.accepted_hits);
                    s.cache_out.map_type = if !s.cache_out.accepted_hits.is_empty() {
                        MappingType::MappedPair
                    } else {
                        MappingType::Unmapped
                    };

                    let (r1_len, r2_len) = if r1.len() <= r2.len() {
                        (r1.len() as i32, r2.len() as i32)
                    } else {
                        (r2.len() as i32, r1.len() as i32)
                    };

                    for hit in &mut s.cache_out.accepted_hits {
                        hit.fragment_length = mate_ov.frag_length as i32;
                        hit.mate_pos = if hit.is_fw {
                            hit.pos + hit.fragment_length - r2_len - 1
                        } else {
                            hit.pos + r1_len - hit.fragment_length - 1
                        };
                        let ref_len = index.ref_len(hit.tid as usize) as i32;
                        if hit.mate_pos > ref_len {
                            hit.mate_pos = hit.pos;
                        }
                        if mate_ov.ov_type == OverlapType::Dovetail {
                            hit.mate_pos = hit.pos;
                        }
                    }
                }
            } else {
                // No overlap: map both ends independently, then merge (bin-aware)
                s.poison_state.paired_for_mapping = true;
                map_pe_fragment_atac::<K>(
                    &r1,
                    &r2,
                    &mut s.hs,
                    &mut s.query,
                    &mut s.cache_left,
                    &mut s.cache_right,
                    &mut s.cache_out,
                    index,
                    &mut s.poison_state,
                    binning,
                );

                if s.cache_out.map_type == MappingType::MappedFirstOrphan {
                    for hit in &mut s.cache_out.accepted_hits {
                        hit.fragment_length = r1.len() as i32;
                    }
                } else if s.cache_out.map_type == MappingType::MappedSecondOrphan {
                    for hit in &mut s.cache_out.accepted_hits {
                        hit.fragment_length = r2.len() as i32;
                    }
                }
            }

            if s.poison_state.is_poisoned() {
                s.local_poisoned += 1;
                continue;
            }

            if s.cache_out.map_type != MappingType::Unmapped {
                s.local_mapped += 1;
                write_atac_record(
                    bc_packed,
                    bc_len,
                    s.cache_out.map_type,
                    &s.cache_out.accepted_hits,
                    tn5_shift,
                    &mut s.rad_writer,
                );
                s.num_reads_in_chunk += 1;
            }
        }
        self.progress.inc(batch_reads);
        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::Result<()> {
        let output = self.output;
        if let Some(s) = &mut self.state {
            s.maybe_flush_chunk(output);
        }
        Ok(())
    }

    fn on_thread_complete(&mut self) -> paraseq::Result<()> {
        let output = self.output;
        let stats = self.stats;
        if let Some(s) = &mut self.state {
            s.finalize_chunk(output);
            s.flush_stats(stats);
        }
        Ok(())
    }
}
