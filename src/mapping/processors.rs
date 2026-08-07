//! Parallel processors for paraseq-based mapping pipelines.
//!
//! Each processor implements the appropriate paraseq parallel processing trait
//! (`ParallelProcessor`, `PairedParallelProcessor`, or `MultiParallelProcessor`)
//! and encapsulates all per-thread mapping state.
//!
//! Per-thread mutable state is stored in `Option<ThreadState>` and initialized
//! lazily on first batch. The custom `Clone` impl sets `state: None` so each
//! cloned processor (one per worker thread) gets fresh state.

use std::sync::Mutex;
use std::sync::atomic::Ordering;

use crate::hash::{FixedMap, fixed_map};
use indicatif::ProgressBar;
use paraseq::Record;
use paraseq::parallel::{MultiParallelProcessor, PairedParallelProcessor, ParallelProcessor};
use seq_geom_parser::normalize::PadBuf;
use smallvec::SmallVec;
use sshash_lib::{Kmer, KmerBits, KmerDictionary};

use crate::index::contig_table::{ContigTable, ContigTableLike};
use crate::index::reference_index::ReferenceIndex;
use crate::io::rad::{
    RadWriter, pack_bases_2bit, write_atac_record, write_bulk_record, write_sc_record,
    write_sc_record_multi_bc,
};
use crate::io::threads::{MappingStats, OutputInfo};
use crate::mapping::binning::BinPos;
use crate::mapping::cache::MappingCache;
use crate::mapping::filters::PoisonState;
use crate::mapping::hit_searcher::{HitSearcher, SkippingStrategy};
use crate::mapping::hits::MappingType;
use crate::mapping::hits::SketchHitInfo;
use crate::mapping::map_fragment::{
    map_pe_fragment, map_pe_fragment_atac, map_se_fragment, map_se_fragment_atac,
};
use crate::mapping::merge_pairs::{remove_duplicate_hits_pub, simple_hit_cmp_bins};
use crate::mapping::overlap::{OverlapType, find_overlap};
use crate::mapping::protocols::Protocol;
use crate::mapping::protocols::scrna::{barcode_has_n, count_ns, is_all_acgt, recover_barcode};
use crate::mapping::sketch_hit_simple::SketchHitInfoSimple;
use crate::mapping::streaming_query::PiscemStreamingQuery;
use crate::mapping::unitig_end_cache::UnitigEndCache;

// ===========================================================================
// Constants
// ===========================================================================

/// Maximum mapped reads per RAD chunk before flushing (matches C++ max_chunk_reads).
const MAX_CHUNK_READS: u32 = 5000;

// ===========================================================================
// MappingOpts
// ===========================================================================

/// Per-cache mapping parameters exposed via CLI flags.
///
/// These are applied to all `MappingCache` instances when a worker thread
/// initializes its state. Defaults match the C++ piscem values.
#[derive(Clone, Copy, Debug)]
pub struct MappingOpts {
    /// Max k-mer occurrence count considered in the first mapping pass.
    pub max_hit_occ: usize,
    /// Max occurrence for the recovery pass (used when all k-mers exceeded
    /// `max_hit_occ` but at least one is below this value).
    pub max_hit_occ_recover: usize,
    /// Reads with more than this many accepted mappings are discarded.
    pub max_read_occ: usize,
    /// Maximum equivalence class cardinality for ambiguous hit filtering.
    pub max_ec_card: u32,
}

impl Default for MappingOpts {
    fn default() -> Self {
        Self {
            max_hit_occ: 256,
            max_hit_occ_recover: 1024,
            max_read_occ: 2500,
            max_ec_card: 4096,
        }
    }
}

impl MappingOpts {
    /// Apply these opts to a `MappingCache`, overriding its defaults.
    pub(crate) fn apply_to<S: crate::mapping::hits::SketchHitInfo>(
        &self,
        cache: &mut crate::mapping::cache::MappingCache<S>,
    ) {
        cache.max_hit_occ = self.max_hit_occ;
        cache.max_hit_occ_recover = self.max_hit_occ_recover;
        cache.attempt_occ_recover = self.max_hit_occ_recover > self.max_hit_occ;
        cache.max_read_occ = self.max_read_occ;
        cache.max_ec_card = self.max_ec_card;
    }
}

// ===========================================================================
// Common thread state helpers
// ===========================================================================

/// Common per-thread mapping state shared across all processor types.
///
/// Generic over `S: SketchHitInfo` to support both the default
/// `SketchHitInfoSimple` (no structural constraints) and
/// `SketchHitInfoChained` (structural constraints enabled via `-c`).
struct CommonThreadState<
    'a,
    const K: usize,
    S: SketchHitInfo + Send + 'static = SketchHitInfoSimple,
    D: KmerDictionary + 'a = sshash_lib::Dictionary,
    C: ContigTableLike + 'a = ContigTable,
> where
    Kmer<K>: KmerBits,
{
    hs: HitSearcher<'a, D, C>,
    query: PiscemStreamingQuery<'a, K, D>,
    cache_out: MappingCache<S>,
    cache_left: MappingCache<S>,
    cache_right: MappingCache<S>,
    poison_state: PoisonState<'a>,
    rad_writer: RadWriter,
    local_reads: u64,
    local_mapped: u64,
    local_poisoned: u64,
    num_reads_in_chunk: u32,
    chunk_in_progress: bool,
}

impl<
    'a,
    const K: usize,
    S: SketchHitInfo + Send + 'static,
    D: KmerDictionary + 'a,
    C: ContigTableLike + 'a,
> CommonThreadState<'a, K, S, D, C>
where
    Kmer<K>: KmerBits,
{
    fn new(
        index: &'a ReferenceIndex<D, C>,
        end_cache: Option<&'a UnitigEndCache>,
        opts: &MappingOpts,
    ) -> Self {
        let query = match end_cache {
            Some(cache) => PiscemStreamingQuery::<K, D>::with_cache(index.dict(), cache),
            None => PiscemStreamingQuery::<K, D>::new(index.dict()),
        };
        let mut cache_out = MappingCache::new(K);
        opts.apply_to(&mut cache_out);
        let mut cache_left = MappingCache::new(K);
        opts.apply_to(&mut cache_left);
        let mut cache_right = MappingCache::new(K);
        opts.apply_to(&mut cache_right);
        Self {
            hs: HitSearcher::new(index),
            query,
            cache_out,
            cache_left,
            cache_right,
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
        self.rad_writer
            .write_u32_at_offset(4, self.num_reads_in_chunk);
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
///
/// Generic over `S: SketchHitInfo` so that structural constraints can be enabled
/// at compile time by instantiating with `SketchHitInfoChained` rather than the
/// default `SketchHitInfoSimple`.
pub struct BulkProcessor<
    'a,
    const K: usize,
    S: SketchHitInfo + Send + 'static = SketchHitInfoSimple,
    D: KmerDictionary + 'a = sshash_lib::Dictionary,
    C: ContigTableLike + 'a = ContigTable,
> where
    Kmer<K>: KmerBits,
{
    index: &'a ReferenceIndex<D, C>,
    end_cache: Option<&'a UnitigEndCache>,
    output: &'a OutputInfo,
    stats: &'a MappingStats,
    progress: &'a ProgressBar,
    strat: SkippingStrategy,
    opts: MappingOpts,
    state: Option<CommonThreadState<'a, K, S, D, C>>,
}

impl<
    'a,
    const K: usize,
    S: SketchHitInfo + Send + 'static,
    D: KmerDictionary + 'a,
    C: ContigTableLike + 'a,
> BulkProcessor<'a, K, S, D, C>
where
    Kmer<K>: KmerBits,
{
    pub fn new(
        index: &'a ReferenceIndex<D, C>,
        end_cache: Option<&'a UnitigEndCache>,
        output: &'a OutputInfo,
        stats: &'a MappingStats,
        strat: SkippingStrategy,
        opts: MappingOpts,
        progress: &'a ProgressBar,
    ) -> Self {
        Self {
            index,
            end_cache,
            output,
            stats,
            progress,
            strat,
            opts,
            state: None,
        }
    }
}

impl<const K: usize, S: SketchHitInfo + Send + 'static, D: KmerDictionary, C: ContigTableLike> Clone
    for BulkProcessor<'_, K, S, D, C>
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
            opts: self.opts,
            state: None,
        }
    }
}

// Safety: all shared fields are `Copy` references; `state` is always `None` at clone time.
unsafe impl<
    const K: usize,
    S: SketchHitInfo + Send + 'static,
    D: KmerDictionary,
    C: ContigTableLike,
> Send for BulkProcessor<'_, K, S, D, C>
where
    Kmer<K>: KmerBits,
{
}

// --- Bulk PE ---

impl<
    'a,
    'r,
    const K: usize,
    S: SketchHitInfo + Send + 'static,
    D: KmerDictionary + 'a,
    C: ContigTableLike + 'a,
> PairedParallelProcessor<paraseq::fastx::RefRecord<'r>> for BulkProcessor<'a, K, S, D, C>
where
    Kmer<K>: KmerBits,
{
    fn process_record_pair_batch(
        &mut self,
        record_pairs: impl Iterator<
            Item = (paraseq::fastx::RefRecord<'r>, paraseq::fastx::RefRecord<'r>),
        >,
    ) -> paraseq::Result<()> {
        let index = self.index;
        let end_cache = self.end_cache;
        let strat = self.strat;
        let s = self
            .state
            .get_or_insert_with(|| CommonThreadState::new(index, end_cache, &self.opts));
        s.ensure_chunk_started();

        // Times the mapping work itself, for the thread broker. Published every
        // few hundred reads rather than once per batch: a 16384-record batch
        // takes far longer than one sampling window, so a per-batch counter
        // would leave windows reading zero work -- reporting maximum starvation
        // for a thread mapping flat out. See `thread_broker::BusyMeter`.
        let mut busy = self.stats.busy.timer();
        let mut batch_reads: u64 = 0;
        for (rec1, rec2) in record_pairs {
            s.local_reads += 1;
            batch_reads += 1;
            busy.tick();
            let seq1 = rec1.seq();
            let seq2 = rec2.seq();

            s.poison_state.paired_for_mapping = true;
            map_pe_fragment(
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

impl<
    'a,
    'r,
    const K: usize,
    S: SketchHitInfo + Send + 'static,
    D: KmerDictionary + 'a,
    C: ContigTableLike + 'a,
> ParallelProcessor<paraseq::fastx::RefRecord<'r>> for BulkProcessor<'a, K, S, D, C>
where
    Kmer<K>: KmerBits,
{
    fn process_record_batch(
        &mut self,
        records: impl Iterator<Item = paraseq::fastx::RefRecord<'r>>,
    ) -> paraseq::Result<()> {
        let index = self.index;
        let end_cache = self.end_cache;
        let strat = self.strat;
        let s = self
            .state
            .get_or_insert_with(|| CommonThreadState::new(index, end_cache, &self.opts));
        s.ensure_chunk_started();

        // Times the mapping work itself, for the thread broker. Published every
        // few hundred reads rather than once per batch: a 16384-record batch
        // takes far longer than one sampling window, so a per-batch counter
        // would leave windows reading zero work -- reporting maximum starvation
        // for a thread mapping flat out. See `thread_broker::BusyMeter`.
        let mut busy = self.stats.busy.timer();
        let mut batch_reads: u64 = 0;
        for rec in records {
            s.local_reads += 1;
            batch_reads += 1;
            busy.tick();
            let seq1 = rec.seq();

            s.poison_state.paired_for_mapping = false;
            map_se_fragment(
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

/// Unmapped barcode counts — two variants for zero overhead in the common case.
///
/// For single-barcode protocols: `Single(AHashMap<u64, u32>)` — same as the
/// pre-multi-barcode code, no SmallVec allocation per lookup.
///
/// For multi-barcode protocols: `Multi(AHashMap<SmallVec<[u64; 2]>, u32>)` —
/// stores per-field barcode values for structured output.
enum UnmappedBcCounts {
    /// Single barcode field — zero overhead compared to pre-multi-barcode code.
    Single {
        counts: FixedMap<u64, u32>,
        bc_len: u16,
    },
    /// Multiple barcode fields — stores per-field values.
    Multi {
        counts: FixedMap<SmallVec<[u64; 2]>, u32>,
        bc_lens: SmallVec<[u16; 2]>,
    },
}

#[derive(Clone, Debug, Default)]
struct TechNormalizationPlan {
    any_normalization: bool,
    bc_needs_padding: SmallVec<[bool; 4]>,
    umi_needs_padding: bool,
}

struct TechNormalizationState {
    bc_pad_buf: PadBuf,
    umi_pad_buf: PadBuf,
}

impl TechNormalizationState {
    fn new() -> Self {
        Self {
            bc_pad_buf: PadBuf::new(),
            umi_pad_buf: PadBuf::new(),
        }
    }
}

impl UnmappedBcCounts {
    fn new_single(bc_len: u16) -> Self {
        UnmappedBcCounts::Single {
            counts: fixed_map(),
            bc_len,
        }
    }

    fn new_multi(bc_lens: SmallVec<[u16; 2]>) -> Self {
        UnmappedBcCounts::Multi {
            counts: fixed_map(),
            bc_lens,
        }
    }

    fn is_empty(&self) -> bool {
        match self {
            UnmappedBcCounts::Single { counts, .. } => counts.is_empty(),
            UnmappedBcCounts::Multi { counts, .. } => counts.is_empty(),
        }
    }

    /// Increment count for a single-barcode key.
    #[inline]
    fn inc_single(&mut self, bc: u64) {
        if let UnmappedBcCounts::Single { counts, .. } = self {
            *counts.entry(bc).or_insert(0) += 1;
        }
    }

    /// Increment count for a multi-barcode key.
    #[inline]
    fn inc_multi(&mut self, bcs: &SmallVec<[u64; 2]>) {
        if let UnmappedBcCounts::Multi { counts, .. } = self {
            *counts.entry(bcs.clone()).or_insert(0) += 1;
        }
    }

    /// Write a single barcode value at its declared width (bases → bytes).
    #[inline]
    fn write_bc_at_width(buf: &mut Vec<u8>, val: u64, bc_len_bases: u16) {
        if bc_len_bases <= 4 {
            buf.extend_from_slice(&(val as u8).to_le_bytes());
        } else if bc_len_bases <= 8 {
            buf.extend_from_slice(&(val as u16).to_le_bytes());
        } else if bc_len_bases <= 16 {
            buf.extend_from_slice(&(val as u32).to_le_bytes());
        } else {
            buf.extend_from_slice(&val.to_le_bytes());
        }
    }

    /// Write all accumulated counts to the given writer.
    /// The format header must have already been written to the file.
    fn flush_to(&self, writer: &mut impl std::io::Write) {
        match self {
            UnmappedBcCounts::Single { counts, bc_len } => {
                let mut buf = Vec::with_capacity(counts.len() * 8);
                for (&bc, &count) in counts {
                    Self::write_bc_at_width(&mut buf, bc, *bc_len);
                    buf.extend_from_slice(&count.to_le_bytes());
                }
                writer.write_all(&buf).ok();
            }
            UnmappedBcCounts::Multi { counts, bc_lens } => {
                let mut buf = Vec::with_capacity(counts.len() * 16);
                for (bc_fields, &count) in counts {
                    for (i, &bc_len) in bc_lens.iter().enumerate() {
                        let val = if i < bc_fields.len() {
                            bc_fields[i]
                        } else {
                            0u64
                        };
                        Self::write_bc_at_width(&mut buf, val, bc_len);
                    }
                    buf.extend_from_slice(&count.to_le_bytes());
                }
                writer.write_all(&buf).ok();
            }
        }
    }
}

#[inline]
fn normalized_barcode_packed(
    raw: &[u8],
    target_len: usize,
    needs_padding: bool,
    pad_buf: &mut PadBuf,
) -> Option<u64> {
    let n_count = count_ns(raw);
    if n_count > 1 {
        return None;
    }

    let recovered = if n_count == 1 {
        recover_barcode(raw)
    } else {
        None
    };
    let seq = match &recovered {
        Some(bc) => bc.as_slice(),
        None => raw,
    };

    if !is_all_acgt(seq) {
        return None;
    }

    if needs_padding && seq.len() < target_len {
        Some(pack_bases_2bit(pad_buf.pad(seq, target_len)))
    } else {
        Some(pack_bases_2bit(seq))
    }
}

#[inline]
fn normalized_umi_packed(
    raw: &[u8],
    target_len: usize,
    needs_padding: bool,
    pad_buf: &mut PadBuf,
) -> Option<u64> {
    if !is_all_acgt(raw) {
        return None;
    }

    if needs_padding && raw.len() < target_len {
        Some(pack_bases_2bit(pad_buf.pad(raw, target_len)))
    } else {
        Some(pack_bases_2bit(raw))
    }
}

/// Per-thread state for scRNA mapping (extends common state).
struct ScrnaThreadState<
    'a,
    const K: usize,
    S: SketchHitInfo + Send + 'static = SketchHitInfoSimple,
    D: KmerDictionary + 'a = sshash_lib::Dictionary,
    C: ContigTableLike + 'a = ContigTable,
> where
    Kmer<K>: KmerBits,
{
    common: CommonThreadState<'a, K, S, D, C>,
    local_rlen_samples: Vec<u32>,
    unmapped_bc_counts: UnmappedBcCounts,
    normalization: Option<TechNormalizationState>,
}

/// Parallel processor for scRNA-seq mapping.
///
/// Generic over `S: SketchHitInfo` so that structural constraints can be enabled
/// by instantiating with `SketchHitInfoChained`.
pub struct ScrnaProcessor<
    'a,
    const K: usize,
    S: SketchHitInfo + Send + 'static = SketchHitInfoSimple,
    D: KmerDictionary + 'a = sshash_lib::Dictionary,
    C: ContigTableLike + 'a = ContigTable,
> where
    Kmer<K>: KmerBits,
{
    index: &'a ReferenceIndex<D, C>,
    end_cache: Option<&'a UnitigEndCache>,
    output: &'a OutputInfo,
    stats: &'a MappingStats,
    progress: &'a ProgressBar,
    strat: SkippingStrategy,
    opts: MappingOpts,
    protocol: &'a dyn Protocol,
    bc_len: u16,
    umi_len: u16,
    normalization_plan: TechNormalizationPlan,
    with_position: bool,
    read_length_samples: &'a Mutex<Vec<u32>>,
    state: Option<ScrnaThreadState<'a, K, S, D, C>>,
}

impl<
    'a,
    const K: usize,
    S: SketchHitInfo + Send + 'static,
    D: KmerDictionary + 'a,
    C: ContigTableLike + 'a,
> ScrnaProcessor<'a, K, S, D, C>
where
    Kmer<K>: KmerBits,
{
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        index: &'a ReferenceIndex<D, C>,
        end_cache: Option<&'a UnitigEndCache>,
        output: &'a OutputInfo,
        stats: &'a MappingStats,
        strat: SkippingStrategy,
        opts: MappingOpts,
        protocol: &'a dyn Protocol,
        bc_len: u16,
        umi_len: u16,
        with_position: bool,
        read_length_samples: &'a Mutex<Vec<u32>>,
        progress: &'a ProgressBar,
    ) -> Self {
        let normalization_plan = protocol
            .normalization_meta()
            .map(|meta| TechNormalizationPlan {
                any_normalization: meta.any_normalization,
                bc_needs_padding: meta.bc_needs_padding.clone(),
                umi_needs_padding: meta.umi_needs_padding,
            })
            .unwrap_or_default();
        Self {
            index,
            end_cache,
            output,
            stats,
            progress,
            strat,
            opts,
            protocol,
            bc_len,
            umi_len,
            normalization_plan,
            with_position,
            read_length_samples,
            state: None,
        }
    }
}

impl<const K: usize, S: SketchHitInfo + Send + 'static, D: KmerDictionary, C: ContigTableLike> Clone
    for ScrnaProcessor<'_, K, S, D, C>
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
            opts: self.opts,
            protocol: self.protocol,
            bc_len: self.bc_len,
            umi_len: self.umi_len,
            normalization_plan: self.normalization_plan.clone(),
            with_position: self.with_position,
            read_length_samples: self.read_length_samples,
            state: None,
        }
    }
}

// Safety: all shared fields are `Copy` references; `state` is always `None` at clone time.
unsafe impl<
    const K: usize,
    S: SketchHitInfo + Send + 'static,
    D: KmerDictionary,
    C: ContigTableLike,
> Send for ScrnaProcessor<'_, K, S, D, C>
where
    Kmer<K>: KmerBits,
{
}

impl<
    'a,
    'r,
    const K: usize,
    S: SketchHitInfo + Send + 'static,
    D: KmerDictionary + 'a,
    C: ContigTableLike + 'a,
> PairedParallelProcessor<paraseq::fastx::RefRecord<'r>> for ScrnaProcessor<'a, K, S, D, C>
where
    Kmer<K>: KmerBits,
{
    fn process_record_pair_batch(
        &mut self,
        record_pairs: impl Iterator<
            Item = (paraseq::fastx::RefRecord<'r>, paraseq::fastx::RefRecord<'r>),
        >,
    ) -> paraseq::Result<()> {
        let index = self.index;
        let end_cache = self.end_cache;
        let strat = self.strat;
        let protocol = self.protocol;
        let bc_len = self.bc_len;
        let umi_len = self.umi_len;
        let normalization_plan = &self.normalization_plan;
        let with_position = self.with_position;
        let max_rlen_samples: usize = 10;

        let is_multi_bc = protocol.is_multi_barcode();
        let bc_descs = protocol.barcode_descs();

        let st = self
            .state
            .get_or_insert_with(|| ScrnaThreadState::<K, S, D, C> {
                common: CommonThreadState::new(index, end_cache, &self.opts),
                local_rlen_samples: Vec::new(),
                unmapped_bc_counts: if is_multi_bc {
                    let bc_lens = protocol.barcode_descs().iter().map(|d| d.len).collect();
                    UnmappedBcCounts::new_multi(bc_lens)
                } else {
                    UnmappedBcCounts::new_single(bc_len)
                },
                normalization: normalization_plan
                    .any_normalization
                    .then(TechNormalizationState::new),
            });
        let s = &mut st.common;
        s.ensure_chunk_started();

        // Pre-compute multi-barcode lens (invariant across reads) and
        // allocate the packed barcode buffer once, cleared each iteration.
        let multi_bc_lens: SmallVec<[u16; 2]> = if is_multi_bc {
            bc_descs.iter().map(|d| d.len).collect()
        } else {
            SmallVec::new()
        };
        let mut multi_bc_packed: SmallVec<[u64; 2]> = SmallVec::new();

        // See the note on the other batch loops: published within the batch,
        // not at its end.
        let mut busy = self.stats.busy.timer();
        let mut batch_reads: u64 = 0;
        for (rec1, rec2) in record_pairs {
            s.local_reads += 1;
            batch_reads += 1;
            busy.tick();

            let r1 = rec1.seq();
            let r2 = rec2.seq();

            // Extract all sequences (barcodes, UMI, bio reads) in one call.
            let seqs = protocol.extract(&r1, &r2);

            let (umi_packed, bc_packed) = if normalization_plan.any_normalization {
                let norm = st
                    .normalization
                    .as_mut()
                    .expect("normalization state missing for variable-length protocol");
                let umi_raw = match seqs.umi {
                    Some(umi) if !umi.is_empty() => umi,
                    _ => continue,
                };
                let umi_packed = match normalized_umi_packed(
                    umi_raw,
                    umi_len as usize,
                    normalization_plan.umi_needs_padding,
                    &mut norm.umi_pad_buf,
                ) {
                    Some(umi) => umi,
                    None => continue,
                };

                let bc_packed = if is_multi_bc {
                    multi_bc_packed.clear();
                    let mut all_valid = true;
                    for (i, _desc) in bc_descs.iter().enumerate() {
                        let bc_raw = match seqs.barcodes.get(i) {
                            Some(Some(bc)) if !bc.is_empty() => *bc,
                            _ => {
                                all_valid = false;
                                break;
                            }
                        };
                        let needs_padding = normalization_plan
                            .bc_needs_padding
                            .get(i)
                            .copied()
                            .unwrap_or(false);
                        match normalized_barcode_packed(
                            bc_raw,
                            multi_bc_lens[i] as usize,
                            needs_padding,
                            &mut norm.bc_pad_buf,
                        ) {
                            Some(bc) => multi_bc_packed.push(bc),
                            None => {
                                all_valid = false;
                                break;
                            }
                        }
                    }
                    if !all_valid {
                        continue;
                    }
                    *multi_bc_packed.last().unwrap_or(&0)
                } else {
                    let bc_raw = match seqs.barcodes.first().copied().flatten() {
                        Some(bc) if !bc.is_empty() => bc,
                        _ => continue,
                    };
                    match normalized_barcode_packed(
                        bc_raw,
                        bc_len as usize,
                        normalization_plan
                            .bc_needs_padding
                            .first()
                            .copied()
                            .unwrap_or(false),
                        &mut norm.bc_pad_buf,
                    ) {
                        Some(bc) => bc,
                        None => continue,
                    }
                };
                (umi_packed, bc_packed)
            } else {
                // UMI validation (matching C++ umi_kmer.fromChars check)
                let umi_raw = match seqs.umi {
                    Some(umi) if !umi.is_empty() => umi,
                    _ => continue,
                };
                if !is_all_acgt(umi_raw) {
                    continue;
                }
                let umi_packed = pack_bases_2bit(umi_raw);

                let bc_packed = if is_multi_bc {
                    multi_bc_packed.clear();
                    let mut all_valid = true;
                    for (i, _desc) in bc_descs.iter().enumerate() {
                        let bc_raw = match seqs.barcodes.get(i) {
                            Some(Some(bc)) if !bc.is_empty() => *bc,
                            _ => {
                                all_valid = false;
                                break;
                            }
                        };
                        let n_count = count_ns(bc_raw);
                        if n_count > 1 {
                            all_valid = false;
                            break;
                        }
                        let recovered = if n_count == 1 {
                            recover_barcode(bc_raw)
                        } else {
                            None
                        };
                        let bc_to_pack = match &recovered {
                            Some(r) => {
                                if !is_all_acgt(r) {
                                    all_valid = false;
                                    break;
                                }
                                r.as_slice()
                            }
                            None => {
                                if !is_all_acgt(bc_raw) {
                                    all_valid = false;
                                    break;
                                }
                                bc_raw
                            }
                        };
                        multi_bc_packed.push(pack_bases_2bit(bc_to_pack));
                    }
                    if !all_valid {
                        continue;
                    }
                    *multi_bc_packed.last().unwrap_or(&0)
                } else {
                    let bc_raw = match seqs.barcodes.first().copied().flatten() {
                        Some(bc) if !bc.is_empty() => bc,
                        _ => continue,
                    };
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
                            if !is_all_acgt(bc_raw) {
                                continue;
                            }
                            bc_raw
                        }
                    };
                    if !is_all_acgt(bc_to_pack) {
                        continue;
                    }
                    pack_bases_2bit(bc_to_pack)
                };
                (umi_packed, bc_packed)
            };

            match seqs.reads.len() {
                2 => {
                    let read1 = seqs.reads[0];
                    let read2 = seqs.reads[1];
                    if read1.is_empty() && read2.is_empty() {
                        continue;
                    }
                    s.poison_state.paired_for_mapping = true;
                    map_pe_fragment(
                        read1,
                        read2,
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
                        st.local_rlen_samples.push(read2.len() as u32);
                    }
                }
                1 => {
                    let bio_seq = seqs.reads[0];
                    if bio_seq.is_empty() {
                        continue;
                    }
                    s.poison_state.paired_for_mapping = false;
                    map_se_fragment(
                        bio_seq,
                        &mut s.hs,
                        &mut s.query,
                        &mut s.cache_out,
                        index,
                        &mut s.poison_state,
                        strat,
                    );

                    if with_position && st.local_rlen_samples.len() < max_rlen_samples {
                        st.local_rlen_samples.push(bio_seq.len() as u32);
                    }
                }
                _ => continue, // no biological reads
            }

            if s.poison_state.is_poisoned() {
                s.local_poisoned += 1;
                continue;
            }

            if s.cache_out.map_type != MappingType::Unmapped {
                s.local_mapped += 1;
                if is_multi_bc {
                    write_sc_record_multi_bc(
                        &multi_bc_packed,
                        &multi_bc_lens,
                        umi_packed,
                        umi_len,
                        s.cache_out.map_type,
                        &s.cache_out.accepted_hits,
                        with_position,
                        &mut s.rad_writer,
                    );
                } else {
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
                }
                s.num_reads_in_chunk += 1;
            } else if is_multi_bc {
                st.unmapped_bc_counts.inc_multi(&multi_bc_packed);
            } else {
                st.unmapped_bc_counts.inc_single(bc_packed);
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
            // Write unmapped barcode counts to shared file using the
            // self-describing format (header already written by map_scrna.rs).
            if let Some(ref unmapped_file) = output.unmapped_bc_file
                && !st.unmapped_bc_counts.is_empty()
            {
                let mut file = unmapped_file.lock().unwrap();
                st.unmapped_bc_counts.flush_to(&mut *file);
            }
        }
        Ok(())
    }
}

// ===========================================================================
// ScatacProcessor
// ===========================================================================

/// Per-thread state for scATAC mapping (extends common state).
struct ScatacThreadState<
    'a,
    const K: usize,
    D: KmerDictionary + 'a = sshash_lib::Dictionary,
    C: ContigTableLike + 'a = ContigTable,
> where
    Kmer<K>: KmerBits,
{
    common: CommonThreadState<'a, K, SketchHitInfoSimple, D, C>,
    unmapped_bc_counts: FixedMap<u64, u32>,
}

/// Parallel processor for scATAC-seq mapping.
///
/// Supports both paired-end (arity=3: R1, barcode, R2) and single-end
/// (arity=2: reads, barcode) via `is_paired`.
pub struct ScatacProcessor<
    'a,
    const K: usize,
    D: KmerDictionary + 'a = sshash_lib::Dictionary,
    C: ContigTableLike + 'a = ContigTable,
> where
    Kmer<K>: KmerBits,
{
    index: &'a ReferenceIndex<D, C>,
    end_cache: Option<&'a UnitigEndCache>,
    output: &'a OutputInfo,
    stats: &'a MappingStats,
    progress: &'a ProgressBar,
    binning: &'a BinPos,
    bc_len: u16,
    tn5_shift: bool,
    min_overlap: i32,
    is_paired: bool,
    opts: MappingOpts,
    state: Option<ScatacThreadState<'a, K, D, C>>,
}

impl<'a, const K: usize, D: KmerDictionary + 'a, C: ContigTableLike + 'a>
    ScatacProcessor<'a, K, D, C>
where
    Kmer<K>: KmerBits,
{
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        index: &'a ReferenceIndex<D, C>,
        end_cache: Option<&'a UnitigEndCache>,
        output: &'a OutputInfo,
        stats: &'a MappingStats,
        binning: &'a BinPos,
        bc_len: u16,
        tn5_shift: bool,
        min_overlap: i32,
        is_paired: bool,
        opts: MappingOpts,
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
            is_paired,
            opts,
            state: None,
        }
    }
}

impl<const K: usize, D: KmerDictionary, C: ContigTableLike> Clone for ScatacProcessor<'_, K, D, C>
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
            is_paired: self.is_paired,
            opts: self.opts,
            state: None,
        }
    }
}

unsafe impl<const K: usize, D: KmerDictionary, C: ContigTableLike> Send
    for ScatacProcessor<'_, K, D, C>
where
    Kmer<K>: KmerBits,
{
}

impl<'a, 'r, const K: usize, D: KmerDictionary + 'a, C: ContigTableLike + 'a>
    MultiParallelProcessor<paraseq::fastx::RefRecord<'r>> for ScatacProcessor<'a, K, D, C>
where
    Kmer<K>: KmerBits,
{
    fn process_multi_record_batch(
        &mut self,
        multi_records: impl Iterator<
            Item = SmallVec<[paraseq::fastx::RefRecord<'r>; paraseq::MAX_ARITY]>,
        >,
    ) -> paraseq::Result<()> {
        let index = self.index;
        let end_cache = self.end_cache;
        let binning = self.binning;
        let bc_len = self.bc_len;
        let tn5_shift = self.tn5_shift;
        let min_overlap = self.min_overlap;
        let is_paired = self.is_paired;

        let st = self.state.get_or_insert_with(|| ScatacThreadState {
            common: CommonThreadState::new(index, end_cache, &self.opts),
            unmapped_bc_counts: fixed_map(),
        });
        let s = &mut st.common;
        s.ensure_chunk_started();

        // Times the mapping work itself, for the thread broker. Published every
        // few hundred reads rather than once per batch: a 16384-record batch
        // takes far longer than one sampling window, so a per-batch counter
        // would leave windows reading zero work -- reporting maximum starvation
        // for a thread mapping flat out. See `thread_broker::BusyMeter`.
        let mut busy = self.stats.busy.timer();
        let mut batch_reads: u64 = 0;
        for multi in multi_records {
            s.local_reads += 1;
            batch_reads += 1;
            busy.tick();

            if is_paired {
                // PE mode: multi[0]=R1, multi[1]=barcode, multi[2]=R2
                if multi.len() < 3 {
                    continue;
                }
            } else {
                // SE mode: multi[0]=bio read, multi[1]=barcode
                if multi.len() < 2 {
                    continue;
                }
            }

            let r1 = multi[0].seq();
            let bc_raw = multi[1].seq();

            if r1.is_empty() {
                continue;
            }

            // For PE, also require r2
            if is_paired && multi[2].seq().is_empty() {
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

            if is_paired {
                let r2 = multi[2].seq();

                // Try mate overlap detection
                let mate_ov = find_overlap(&r1, &r2, min_overlap, 0);

                if mate_ov.ov_type != OverlapType::NoOverlap && !mate_ov.frag.is_empty() {
                    // Map merged fragment as single read (bin-based)
                    s.poison_state.paired_for_mapping = false;
                    map_se_fragment_atac(
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
                    map_pe_fragment_atac(
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
            } else {
                // SE mode: map the single biological read
                s.poison_state.paired_for_mapping = false;
                map_se_fragment_atac(
                    &r1,
                    &mut s.hs,
                    &mut s.query,
                    &mut s.cache_out,
                    index,
                    &mut s.poison_state,
                    binning,
                );

                if !s.cache_out.accepted_hits.is_empty() {
                    s.cache_out.accepted_hits.sort_by(simple_hit_cmp_bins);
                    remove_duplicate_hits_pub(&mut s.cache_out.accepted_hits);
                    s.cache_out.map_type = if !s.cache_out.accepted_hits.is_empty() {
                        MappingType::SingleMapped
                    } else {
                        MappingType::Unmapped
                    };
                    for hit in &mut s.cache_out.accepted_hits {
                        hit.fragment_length = r1.len() as i32;
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
            } else {
                *st.unmapped_bc_counts.entry(bc_packed).or_insert(0) += 1;
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
        if let Some(st) = &mut self.state {
            st.common.finalize_chunk(output);
            st.common.flush_stats(stats);
            // Write unmapped barcode counts to shared file
            if let Some(ref unmapped_file) = output.unmapped_bc_file
                && !st.unmapped_bc_counts.is_empty()
            {
                let mut buf = Vec::with_capacity(st.unmapped_bc_counts.len() * 12);
                for (&bc, &count) in &st.unmapped_bc_counts {
                    buf.extend_from_slice(&bc.to_le_bytes());
                    buf.extend_from_slice(&count.to_le_bytes());
                }
                let mut file = unmapped_file.lock().unwrap();
                use std::io::Write;
                file.write_all(&buf).ok();
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalized_barcode_packing_uses_fixed_width() {
        let mut pad = PadBuf::new();
        let raw = normalized_barcode_packed(b"ACGTACGTA", 10, false, &mut pad).unwrap();
        let padded = normalized_barcode_packed(b"ACGTACGTA", 10, true, &mut pad).unwrap();
        assert_ne!(raw, padded);
    }

    #[test]
    fn normalized_umi_packing_uses_fixed_width() {
        let mut pad = PadBuf::new();
        let raw = normalized_umi_packed(b"ACGTACGTAC", 12, false, &mut pad).unwrap();
        let padded = normalized_umi_packed(b"ACGTACGTAC", 12, true, &mut pad).unwrap();
        assert_ne!(raw, padded);
    }
}
