//! CLI command for scATAC-seq mapping.
//!
//! scATAC reads use triple-file input: R1 (genomic), barcode, R2 (genomic).
//! Supports mate overlap detection (merged mapping), bin-based hit collection,
//! and optional Tn5 shift.
//!
//! Key differences from scRNA mapping:
//! - Uses `get_raw_hits_sketch_everykmer` (not STRICT/PERMISSIVE skipping)
//! - Bin-based hit accumulation with threshold filtering
//! - Bin-aware paired-end merge
//! - Poison filtering disabled by default (`--no-poison` is true)

use std::io::{Seek, SeekFrom, Write};
use std::path::PathBuf;
use std::sync::atomic::Ordering;
#[cfg(feature = "rapidgzip")]
use std::time::Duration;
use std::time::Instant;

use anyhow::{Context, Result};
use clap::Args;
use tracing::info;

use sshash_lib::{Kmer, KmerBits, KmerDictionary, dispatch_on_k};

use crate::index::contig_table::ContigTableLike;

use super::DictKind;
use crate::index::reference_index::{ReferenceIndex, tiny_artifacts_exist};
use crate::io::fastx::{Collection, CollectionType};
use crate::io::map_info::{MapInfoParams, write_map_info};
use crate::io::rad::write_rad_header_atac;
use crate::io::threads::{MappingStats, OutputInfo};
use crate::mapping::binning::BinPos;
use crate::mapping::processors::{MappingOpts, SCATAC_PROGRESS_FLUSH_EVERY, ScatacProcessor};
use crate::mapping::unitig_end_cache::UnitigEndCache;

use super::map_bulk::make_progress_bar;

/// Same-binary validation hook for the scATAC broker progress cadence.
///
/// This is intentionally not a user-facing CLI option. It exists so the
/// 64-record default can be compared with the generic 256-record cadence
/// without changing any executable code between benchmark arms.
const SCATAC_PROGRESS_FLUSH_EVERY_ENV: &str = "PISCEM_SCATAC_PROGRESS_FLUSH_EVERY";
/// Same-binary validation hook for scATAC reader batch geometry.
const SCATAC_READER_BATCH_SIZE_ENV: &str = "PISCEM_SCATAC_READER_BATCH_SIZE";

fn parse_progress_flush_every(value: Option<&str>) -> Result<Option<u64>> {
    let Some(value) = value else {
        return Ok(None);
    };
    let flush_every = value.parse::<u64>().with_context(|| {
        format!("invalid {SCATAC_PROGRESS_FLUSH_EVERY_ENV}={value:?}; expected a positive integer")
    })?;
    anyhow::ensure!(
        flush_every > 0,
        "invalid {SCATAC_PROGRESS_FLUSH_EVERY_ENV}=0; value must be non-zero"
    );
    Ok(Some(flush_every))
}

fn parse_reader_batch_size(value: Option<&str>) -> Result<Option<usize>> {
    let Some(value) = value else {
        return Ok(None);
    };
    let batch_size = value.parse::<usize>().with_context(|| {
        format!("invalid {SCATAC_READER_BATCH_SIZE_ENV}={value:?}; expected a positive integer")
    })?;
    anyhow::ensure!(
        batch_size > 0,
        "invalid {SCATAC_READER_BATCH_SIZE_ENV}=0; value must be non-zero"
    );
    Ok(Some(batch_size))
}

fn scatac_pipeline_geometry(
    allocation: crate::io::fastx::DecodeAllocation,
    _budget: usize,
    progress_override: Option<u64>,
    batch_override: Option<usize>,
) -> crate::io::fastx::PipelineTuning {
    let (default_batch, default_progress) = match allocation {
        crate::io::fastx::DecodeAllocation::Adaptive => (
            crate::io::fastx::SCATAC_READER_BATCH_SIZE,
            SCATAC_PROGRESS_FLUSH_EVERY,
        ),
        crate::io::fastx::DecodeAllocation::Serial => (
            crate::io::fastx::SCATAC_READER_BATCH_SIZE,
            thread_broker::DEFAULT_FLUSH_EVERY,
        ),
        crate::io::fastx::DecodeAllocation::PinnedAggregate { .. }
        | crate::io::fastx::DecodeAllocation::PinnedPerFile { .. } => (
            crate::io::fastx::SCATAC_READER_BATCH_SIZE,
            thread_broker::DEFAULT_FLUSH_EVERY,
        ),
    };
    crate::io::fastx::PipelineTuning {
        reader_batch_size: batch_override.unwrap_or(default_batch),
        progress_flush_every: progress_override.unwrap_or(default_progress),
    }
}

#[cfg(feature = "rapidgzip")]
fn scatac_broker_config(budget: usize) -> thread_broker::BrokerConfig {
    let mut config = crate::io::fastx::broker_config_for_budget(budget);
    // Stable scATAC runs converge to one split and hold it. Retain responsive
    // regime-change handling, but make its ordinary cadence sparse by default;
    // the environment override can restore 25 ms monitoring for applications
    // that need rapid post-convergence reaction.
    config.steady_probe_interval = Some(Duration::from_secs(5));
    // scATAC has measured allocation-dependent consumer scaling. Confirm a
    // disagreement between the opening and cost model under the crate's
    // bounded startup defaults; do not encode a performance belief in either
    // safety floor. The true shared-pool decoder floor remains one.
    config.opening_policy = thread_broker::OpeningPolicy::Bracket(Default::default());
    config
}

fn scatac_initial_decode_slots(budget: usize) -> usize {
    // Provenance: 2M-read 10x scATAC grid, 2026-08-07, whose t8 peak is
    // producer five and whose t32 low end is 1--3. Four is intentionally only
    // an opening hint: the startup bracket can leave it in either direction.
    4.min(budget.saturating_sub(1)).max(1)
}

#[derive(Args, Debug)]
#[command(group(
    clap::ArgGroup::new("bio_reads")
        .required(true)
        .args(["reads", "read1"])
))]
pub struct MapScatacArgs {
    /// Index prefix path
    #[arg(short = 'i', long)]
    pub index: PathBuf,
    /// Single-end biological FASTQ files (comma-separated); mutually exclusive with -1/-2
    #[arg(short = 'r', long = "reads", value_delimiter = ',',
          conflicts_with_all = ["read1", "read2"])]
    pub reads: Vec<PathBuf>,
    /// Read 1 FASTQ files (genomic left, comma-separated); requires -2
    #[arg(
        short = '1',
        long,
        value_delimiter = ',',
        requires = "read2",
        conflicts_with = "reads"
    )]
    pub read1: Vec<PathBuf>,
    /// Read 2 FASTQ files (genomic right, comma-separated); requires -1
    #[arg(
        short = '2',
        long,
        value_delimiter = ',',
        requires = "read1",
        conflicts_with = "reads"
    )]
    pub read2: Vec<PathBuf>,
    /// Barcode FASTQ files (comma-separated)
    #[arg(short = 'b', long, value_delimiter = ',')]
    pub barcode: Vec<PathBuf>,
    /// Output directory
    #[arg(short = 'o', long)]
    pub output: PathBuf,
    /// Total execution-slot budget shared by mapping and gzip decoding
    #[arg(short = 't', long, default_value = "16")]
    pub threads: usize,

    /// Gzip decoder selection: `auto`, `serial`, `parallel`, or `parallel=N`.
    ///
    /// `auto` adapts the aggregate mapping/decode split during the real run.
    /// `serial` gives mapping the full budget. `parallel` forces the parallel
    /// path but still adapts its split; `parallel=N` fixes N slots per
    /// decoder-capable input and disables adaptation. Non-regular inputs remain
    /// serial because the parallel decoder requires positional reads.
    #[arg(long, default_value = "auto", value_name = "MODE")]
    pub decoder: String,
    /// Path to a JSON file overriding thread and decoder policy.
    ///
    /// Every field is optional and defaults to a measured value, so a file need
    /// only name what it changes; an unknown field is an error rather than a
    /// silent no-op. Currently:
    ///
    /// `{"parallel_decode": {"min_threads_per_stream": 8}}`
    ///
    /// That value is how many threads must be available per gzip input before
    /// the parallel decoder is used at all. Below it the serial decoder already
    /// supplies one inflate stream per input for free, on threads that also map.
    #[arg(long, value_name = "FILE")]
    pub thread_policy: Option<PathBuf>,
    /// Barcode length in bases
    #[arg(long, default_value = "16")]
    pub bc_len: usize,
    /// Disable Tn5 transposase shift
    #[arg(long)]
    pub no_tn5_shift: bool,
    /// Disable poison k-mer filtering (default: true, matching C++ behavior)
    #[arg(long, default_value = "true")]
    pub no_poison: bool,
    /// Enable ambiguous hit checking (EC table loading); disabled by default for ATAC
    #[arg(long)]
    pub check_ambig_hits: bool,
    /// Maximum equivalence class cardinality for ambiguous hit resolution
    #[arg(long, default_value = "256")]
    pub max_ec_card: u32,
    /// Maximum k-mer occurrence count considered in the first mapping pass
    #[arg(long, default_value = "256")]
    pub max_hit_occ: usize,
    /// Maximum occurrence for the recovery pass
    #[arg(long, default_value = "1024")]
    pub max_hit_occ_recover: usize,
    /// Reads with more than this many accepted mappings are discarded
    #[arg(long, default_value = "2500")]
    pub max_read_occ: usize,
    /// UnitigEndCache capacity (number of entries)
    #[arg(long, default_value = "5000000")]
    pub end_cache_capacity: usize,
    /// Bin size for genomic binning
    #[arg(long, default_value = "1000")]
    pub bin_size: u64,
    /// Overlap between adjacent bins
    #[arg(long, default_value = "300")]
    pub bin_overlap: u64,
    /// Hit threshold fraction for bin-based filtering
    #[arg(long, default_value = "0.7")]
    pub thr: f32,
    /// Minimum overlap length for mate merging
    #[arg(long, default_value = "30")]
    pub min_overlap: i32,
    /// K-mer skipping strategy (ignored for ATAC — always uses every-kmer)
    #[arg(long)]
    pub skipping_strategy: Option<String>,
    /// Suppress progress output
    #[arg(short = 'q', long)]
    pub quiet: bool,
    /// K-mer dictionary backend: `sshash` (compact, default) or `tiny`
    /// (hashbrown-backed, faster but higher memory).
    #[arg(long, value_enum, default_value_t = DictKind::Auto)]
    pub dict: DictKind,
}

pub fn run(args: MapScatacArgs) -> Result<()> {
    let start = Instant::now();

    if args.skipping_strategy.is_some() {
        info!("Note: --skipping-strategy is ignored for scATAC (always uses every-kmer mode)");
    }

    info!(
        "scATAC mapping (bc_len={}, tn5_shift={}, bin_size={}, overlap={}, thr={:.2})",
        args.bc_len, !args.no_tn5_shift, args.bin_size, args.bin_overlap, args.thr,
    );

    let is_paired = !args.read1.is_empty();

    let effective_dict = args.dict.resolve_for_map(&args.index);
    let load_start = Instant::now();
    if matches!(effective_dict, super::DictKind::Tiny) && tiny_artifacts_exist(&args.index) {
        info!(
            "Loading prebuilt Tiny index artifacts from {}.{{tdct,tct}}",
            args.index.display()
        );
        let index = crate::index::reference_index::ReferenceIndex::<
            tiny_dict::TinyDictionary,
            crate::index::contig_table::TinyContigTable,
        >::load_tiny(&args.index, args.check_ambig_hits, !args.no_poison)?;
        info!(
            "Index loaded: k={}, {} refs ({:.2}s)",
            index.k(),
            index.num_refs(),
            load_start.elapsed().as_secs_f64()
        );
        return run_scatac_with_index(args, &index, start, is_paired);
    }
    info!("Loading index from {}", args.index.display());
    let index = ReferenceIndex::load(&args.index, args.check_ambig_hits, !args.no_poison)?;
    info!(
        "Index loaded: k={}, {} refs ({:.2}s)",
        index.k(),
        index.num_refs(),
        load_start.elapsed().as_secs_f64()
    );

    if matches!(effective_dict, super::DictKind::Tiny) {
        info!("Converting sshash index to Tiny (in-memory)");
        let convert_start = Instant::now();
        let k = index.k();
        let tiny_index = sshash_lib::dispatch_on_k!(k, K => index.into_tiny_full::<K>());
        let inline_frac = tiny_index.contig_table().inline_fraction();
        info!(
            "Tiny index ready ({:.2}s); single-ref inline fraction = {:.2}%",
            convert_start.elapsed().as_secs_f64(),
            inline_frac * 100.0,
        );
        return run_scatac_with_index(args, &tiny_index, start, is_paired);
    }

    run_scatac_with_index(args, &index, start, is_paired)
}

#[allow(clippy::too_many_arguments)]
fn run_scatac_with_index<D, C>(
    args: MapScatacArgs,
    index: &ReferenceIndex<D, C>,
    start: Instant,
    is_paired: bool,
) -> Result<()>
where
    D: KmerDictionary + Sync,
    C: ContigTableLike + Sync,
{
    // Create output directory and RAD file
    let out_dir = args.output.clone();
    std::fs::create_dir_all(&out_dir)
        .with_context(|| format!("failed to create output directory: {}", out_dir.display()))?;
    let rad_path = out_dir.join("map.rad");
    let mut rad_file = std::fs::File::create(&rad_path)
        .with_context(|| format!("failed to create {}", rad_path.display()))?;

    // Write ATAC RAD header
    let ref_names: Vec<&str> = (0..index.num_refs()).map(|i| index.ref_name(i)).collect();
    let ref_lengths: Vec<u32> = (0..index.num_refs())
        .map(|i| index.ref_len(i) as u32)
        .collect();
    let bc_len = args.bc_len as u16;
    let chunk_count_offset = write_rad_header_atac(
        &mut rad_file,
        index.num_refs() as u64,
        &ref_names,
        &ref_lengths,
        bc_len,
    )?;

    // Create binning scheme
    let binning = BinPos::new(index, args.bin_size, args.bin_overlap, args.thr);
    info!("Binning: {} total bins", binning.num_bins());

    // Create unmapped barcode count file
    let unmapped_bc_path = out_dir.join("unmapped_bc_count.bin");
    let unmapped_bc_file = std::fs::File::create(&unmapped_bc_path)
        .with_context(|| format!("failed to create {}", unmapped_bc_path.display()))?;

    // Setup shared state
    let stats = MappingStats::new();
    let output_info = OutputInfo {
        num_chunks: std::sync::atomic::AtomicUsize::new(0),
        rad_file: std::sync::Mutex::new(std::io::BufWriter::new(
            rad_file
                .try_clone()
                .context("failed to clone RAD file handle")?,
        )),
        unmapped_bc_file: Some(std::sync::Mutex::new(std::io::BufWriter::new(
            unmapped_bc_file,
        ))),
    };

    // Setup progress bar
    let progress = make_progress_bar(args.quiet);

    let k = index.k();
    let tn5_shift = !args.no_tn5_shift;
    let min_overlap = args.min_overlap;
    let end_cache = UnitigEndCache::new(args.end_cache_capacity);
    let num_threads = args.threads.max(1);
    let decoder_pref = crate::io::calibrate::DecoderPreference::parse(&args.decoder)
        .map_err(|e| anyhow::anyhow!("invalid --decoder value: {e}"))?;
    let thread_policy = crate::io::policy::ThreadPolicy::load(args.thread_policy.as_deref())?;
    let progress_flush_override = parse_progress_flush_every(
        std::env::var(SCATAC_PROGRESS_FLUSH_EVERY_ENV)
            .ok()
            .as_deref(),
    )?;
    let reader_batch_override =
        parse_reader_batch_size(std::env::var(SCATAC_READER_BATCH_SIZE_ENV).ok().as_deref());
    let reader_batch_override = reader_batch_override?;
    let opts = MappingOpts {
        max_hit_occ: args.max_hit_occ,
        max_hit_occ_recover: args.max_hit_occ_recover,
        max_read_occ: args.max_read_occ,
        max_ec_card: args.max_ec_card,
    };

    let index_k = index.k();
    let index_m = index.m();
    let index_num_refs = index.num_refs();
    let sig_info_owned = index.ref_sig_info().cloned();

    let outcome = dispatch_on_k!(k, K => {
        let (bio_paths, r2_paths) = if is_paired {
            (args.read1.as_slice(), args.read2.as_slice())
        } else {
            (args.reads.as_slice(), [].as_slice())
        };
        run_atac_pipeline::<K, _, _>(
            bio_paths, &args.barcode, r2_paths,
            &output_info, &stats,
            index, &binning, bc_len, tn5_shift, min_overlap, Some(&end_cache),
            is_paired, progress_flush_override, reader_batch_override, opts,
            num_threads, decoder_pref, thread_policy, &progress,
        )?
    });

    progress.finish_and_clear();

    // Backpatch num_chunks
    let num_chunks = output_info.num_chunks.load(Ordering::Relaxed) as u64;
    drop(output_info);
    rad_file.seek(SeekFrom::Start(chunk_count_offset))?;
    rad_file.write_all(&num_chunks.to_le_bytes())?;
    drop(rad_file);

    let elapsed = start.elapsed().as_secs_f64();
    let (num_reads, num_mapped, num_poisoned) = stats.summary();

    info!(
        "Mapped {}/{} reads ({:.1}%), {} poisoned, {:.1}s",
        num_mapped,
        num_reads,
        if num_reads > 0 {
            num_mapped as f64 / num_reads as f64 * 100.0
        } else {
            0.0
        },
        num_poisoned,
        elapsed,
    );

    // Write map_info.json
    write_map_info(&MapInfoParams {
        path: &out_dir.join("map_info.json"),
        mode: "sc-atac",
        num_reads,
        num_mapped,
        num_poisoned,
        elapsed_secs: elapsed,
        mapping_elapsed_secs: outcome.mapping_elapsed_secs,
        sig_info: sig_info_owned.as_ref(),
        piscem_rs_version: crate::VERSION,
        num_threads: outcome.execution_plan.effective_budget,
        execution_plan: Some(&outcome.execution_plan),
        broker_report: outcome.broker_report.as_ref(),
        broker_failure: outcome.broker_failure.as_ref(),
        producer_measurement: outcome.producer_measurement.as_ref(),
        consumer_measurement: Some(&outcome.consumer_measurement),
        pipeline_tuning: outcome.tuning.as_ref(),
        index_path: &args.index,
        k: index_k,
        m: index_m,
        num_refs: index_num_refs,
        skipping_strategy: "every-kmer",
    })?;

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn run_atac_pipeline<const K: usize, D: KmerDictionary + Sync, C: ContigTableLike + Sync>(
    bio_paths: &[PathBuf],
    barcode_paths: &[PathBuf],
    read2_paths: &[PathBuf],
    output: &OutputInfo,
    stats: &MappingStats,
    index: &ReferenceIndex<D, C>,
    binning: &BinPos,
    bc_len: u16,
    tn5_shift: bool,
    min_overlap: i32,
    end_cache: Option<&UnitigEndCache>,
    is_paired: bool,
    progress_flush_override: Option<u64>,
    reader_batch_override: Option<usize>,
    opts: MappingOpts,
    num_threads: usize,
    decoder_pref: crate::io::calibrate::DecoderPreference,
    thread_policy: crate::io::policy::ThreadPolicy,
    progress: &indicatif::ProgressBar,
) -> Result<crate::io::fastx::PipelineOutcome>
where
    Kmer<K>: KmerBits,
{
    let mut processor = ScatacProcessor::<K, D, C>::new(
        index,
        end_cache,
        output,
        stats,
        binning,
        bc_len,
        tn5_shift,
        min_overlap,
        is_paired,
        opts,
        progress,
    );

    let arity = if is_paired { 3 } else { 2 };
    let num_input_files = bio_paths.len() * arity;

    // scATAC opens `arity` files per sample. Preserve that grouping so mixed
    // regular/stream input is handled per file without inspecting a stream.
    let groups: Vec<crate::io::calibrate::ReadGroup> = (0..bio_paths.len())
        .map(|i| {
            let mut g = vec![bio_paths[i].clone(), barcode_paths[i].clone()];
            if is_paired {
                g.push(read2_paths[i].clone());
            }
            g
        })
        .collect();
    let decision = crate::io::calibrate::choose_decoder(
        &groups,
        num_threads,
        decoder_pref,
        thread_policy.parallel_decode,
    );
    #[cfg_attr(not(feature = "rapidgzip"), allow(unused_mut))]
    let mut plan = crate::io::fastx::plan_thread_budget(
        num_threads,
        num_input_files,
        decision.parallel,
        decoder_pref,
    )?;
    if plan.adaptive() {
        // scATAC is the measured allocation-dependent modality. Four is one
        // provenance-backed opening hint at every budget; the bounded startup
        // bracket, rather than a budget-specific safety floor, decides whether
        // to leave it. Fixed and serial requests retain exact semantics.
        plan.decode_slots = scatac_initial_decode_slots(plan.effective_budget);
        plan.map_threads = plan.effective_budget - plan.decode_slots;
    }

    // One pool for the whole run, sized to the entire budget: `workers` is an
    // immutable maximum and a refused `set_worker_limit` would silently desync
    // the broker's accounting from the pool's.
    #[cfg(feature = "rapidgzip")]
    let decode_pool = if plan.parallel_gzip() {
        Some(
            rapidgzip_core::DecoderPool::builder()
                .workers(plan.effective_budget)
                .initial_worker_limit(plan.decode_slots)
                .build()
                .map_err(|e| anyhow::anyhow!("could not create decoder pool: {e}"))?,
        )
    } else {
        None
    };
    #[cfg(feature = "rapidgzip")]
    let mut handles: Vec<rapidgzip_core::DecoderHandle> = Vec::new();

    // Open every input exactly once before choosing reader geometry. Decoder
    // handle reconciliation can turn an initially adaptive plan into serial
    // (for example, regular but non-gzip input), and only the final plan should
    // pay adaptive scATAC's smaller batches/finer progress cadence.
    let open = |path: &std::path::Path| -> Result<_> {
        #[cfg(feature = "rapidgzip")]
        let opened =
            crate::io::fastx::open_input_pooled(path, decode_pool.as_ref(), plan.effective_budget)?;
        #[cfg(not(feature = "rapidgzip"))]
        let opened = crate::io::fastx::open_input(path, 0)?;
        Ok(opened)
    };

    let mut opened_inputs = Vec::with_capacity(num_input_files);
    for i in 0..bio_paths.len() {
        for path in [&bio_paths[i], &barcode_paths[i]] {
            let opened = open(path)?;
            #[cfg(feature = "rapidgzip")]
            handles.extend(opened.handle.clone());
            opened_inputs.push((path.clone(), opened));
        }
        if is_paired {
            let path = &read2_paths[i];
            let opened = open(path)?;
            #[cfg(feature = "rapidgzip")]
            handles.extend(opened.handle.clone());
            opened_inputs.push((path.clone(), opened));
        }
    }

    #[cfg(feature = "rapidgzip")]
    {
        plan.reconcile_parallel_decoders(handles.len());
        if plan.parallel_gzip()
            && let Some(pool) = &decode_pool
        {
            pool.set_worker_limit(plan.decode_slots)
                .map_err(|e| anyhow::anyhow!("could not apply decoder execution plan: {e}"))?;
        }
    }

    let tuning = scatac_pipeline_geometry(
        plan.allocation,
        plan.effective_budget,
        progress_flush_override,
        reader_batch_override,
    );
    stats.set_broker_progress_flush_every(tuning.progress_flush_every);
    tracing::info!(
        reader_batch_size = tuning.reader_batch_size,
        progress_flush_every = tuning.progress_flush_every,
        adaptive = plan.adaptive(),
        "configured scATAC pipeline geometry"
    );

    let readers = opened_inputs
        .into_iter()
        .map(|(path, opened)| {
            crate::io::fastx::scatac_reader(opened.reader, tuning.reader_batch_size)
                .map_err(|e| anyhow::anyhow!("failed to open {}: {}", path.display(), e))
        })
        .collect::<Result<Vec<_>>>()?;

    tracing::info!(
        requested_budget = plan.requested_budget,
        effective_budget = plan.effective_budget,
        mapping_threads = plan.map_threads,
        decode_slots = plan.decode_slots,
        allocation = ?plan.allocation,
        "thread execution plan"
    );

    // The mapping side, resizable for the same reason.
    let map_pool = paraseq::parallel::ThreadPool::with_max(plan.map_threads, plan.effective_budget);
    #[cfg(feature = "rapidgzip")]
    let consumer_floor =
        crate::io::fastx::collection_share_floor(readers.len(), arity, plan.map_threads);

    #[cfg(feature = "rapidgzip")]
    let broker = match (&decode_pool, plan.adaptive()) {
        (Some(pool), true) => {
            let broker_config = scatac_broker_config(plan.effective_budget);
            match crate::io::broker::DecodeProducer::new(pool.clone(), handles.clone()) {
                Ok(producer) => {
                    let built = thread_broker::ThreadBroker::builder_with(
                        crate::io::broker::MappingConsumer::new(map_pool.clone(), stats)
                            .with_progress_cadence(
                                stats,
                                tuning.progress_flush_every,
                                thread_broker::DEFAULT_FLUSH_EVERY,
                            ),
                        producer,
                        crate::io::broker::broker_config_from_environment(broker_config).map_err(
                            |e| anyhow::anyhow!("invalid thread broker probe interval: {e}"),
                        )?,
                    )
                    .budget(plan.effective_budget)
                    .initial_producer_slots(plan.decode_slots)
                    .min_consumer_threads(consumer_floor)
                    .steady_state_policy(
                        crate::io::broker::broker_policy_from_environment()
                            .map_err(|e| anyhow::anyhow!("invalid thread broker policy: {e}"))?,
                    )
                    .build()
                    .map_err(|e| anyhow::anyhow!("invalid thread broker configuration: {e}"))?;
                    crate::io::broker::AdvisoryBroker::start(
                        built,
                        plan.map_threads,
                        plan.decode_slots,
                    )
                }
                Err(error) => crate::io::broker::AdvisoryBroker::failed(
                    crate::io::broker::BrokerFailureStage::ProducerMeasurementStartup,
                    error,
                    plan.map_threads,
                    plan.decode_slots,
                ),
            }
        }
        _ => crate::io::broker::AdvisoryBroker::disabled(),
    };

    #[cfg(feature = "rapidgzip")]
    let fixed_decode_measurement = crate::io::broker::FixedDecodeMeasurement::from_environment(
        handles.clone(),
        plan.adaptive(),
    )
    .map_err(|e| anyhow::anyhow!("invalid fixed decode measurement control: {e}"))?;

    let measurement_at_start = stats.measurement_snapshot();
    let wall_start = std::time::Instant::now();

    let collection = Collection::new(readers, CollectionType::Multi { arity })
        .map_err(|e| anyhow::anyhow!("failed to create collection: {}", e))?;
    let result = collection
        .process_parallel_multi_pool(&mut processor, &map_pool, None)
        .map_err(|e| anyhow::anyhow!("mapping failed: {}", e));
    let mapping_elapsed = wall_start.elapsed();
    let consumer_measurement =
        stats.log_measurement(measurement_at_start, mapping_elapsed, plan.effective_budget);

    #[cfg(feature = "rapidgzip")]
    let fixed_producer_measurement =
        fixed_decode_measurement.map(crate::io::broker::FixedDecodeMeasurement::finish);

    #[cfg(feature = "rapidgzip")]
    let mut broker_diagnostics = broker.finish();
    #[cfg(feature = "rapidgzip")]
    if let Some(r) = &mut broker_diagnostics.report {
        if let Some(measurement) = &mut r.producer_measurement {
            crate::io::broker::refresh_measurement_cpu(measurement, &handles);
        }
        tracing::info!(
            "thread broker: settled at {} mapping / {} decode of {} \
             ({} moves, {} reverts, {} resurveys, converged {}, {:?} to settle)",
            r.final_consumer_threads,
            r.final_producer_limit,
            plan.effective_budget,
            r.moves,
            r.reverts,
            r.resurveys,
            r.converged(),
            r.time_to_converge,
        );
        // The measured decode share of total per-read cost, which is what the
        // split is solved from. Worth logging on its own: if a split looks
        // wrong, this says whether the model was misinformed or merely
        // overruled by the cap.
        if let Some(m) = r.final_model {
            tracing::info!(
                "thread broker: decode is {:.0}% of per-read cost -> wanted {} slots, \
                 usable ceiling {}",
                m.producer_cost_share * 100.0,
                m.ideal_producer_slots,
                if m.useful_cap == usize::MAX {
                    "none".to_string()
                } else {
                    m.useful_cap.to_string()
                },
            );
        }
    }
    #[cfg(feature = "rapidgzip")]
    let broker_report = broker_diagnostics.report;
    #[cfg(feature = "rapidgzip")]
    let broker_failure = broker_diagnostics.failure;
    #[cfg(not(feature = "rapidgzip"))]
    let broker_report = None;
    #[cfg(not(feature = "rapidgzip"))]
    let broker_failure = None;

    #[cfg(feature = "rapidgzip")]
    let producer_measurement = broker_report
        .as_ref()
        .and_then(|report| report.producer_measurement)
        .or(fixed_producer_measurement);
    #[cfg(not(feature = "rapidgzip"))]
    let producer_measurement = None;

    result?;
    Ok(crate::io::fastx::PipelineOutcome {
        execution_plan: plan,
        broker_report,
        broker_failure,
        producer_measurement,
        consumer_measurement,
        tuning: Some(tuning),
        mapping_elapsed_secs: mapping_elapsed.as_secs_f64(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scatac_pipeline_geometry_is_policy_scoped_and_validates_overrides() {
        assert_eq!(parse_progress_flush_every(None).unwrap(), None);
        assert_eq!(parse_progress_flush_every(Some("256")).unwrap(), Some(256));
        assert!(parse_progress_flush_every(Some("0")).is_err());
        assert!(parse_progress_flush_every(Some("not-a-number")).is_err());
        assert_eq!(parse_reader_batch_size(None).unwrap(), None);
        assert_eq!(parse_reader_batch_size(Some("2048")).unwrap(), Some(2048));
        assert!(parse_reader_batch_size(Some("0")).is_err());

        assert_eq!(
            scatac_pipeline_geometry(crate::io::fastx::DecodeAllocation::Adaptive, 8, None, None),
            crate::io::fastx::PipelineTuning {
                reader_batch_size: crate::io::fastx::SCATAC_READER_BATCH_SIZE,
                progress_flush_every: SCATAC_PROGRESS_FLUSH_EVERY,
            }
        );
        assert_eq!(
            scatac_pipeline_geometry(crate::io::fastx::DecodeAllocation::Adaptive, 32, None, None),
            crate::io::fastx::PipelineTuning {
                reader_batch_size: crate::io::fastx::SCATAC_READER_BATCH_SIZE,
                progress_flush_every: SCATAC_PROGRESS_FLUSH_EVERY,
            }
        );
        assert_eq!(
            scatac_pipeline_geometry(crate::io::fastx::DecodeAllocation::Serial, 8, None, None),
            crate::io::fastx::PipelineTuning {
                reader_batch_size: crate::io::fastx::SCATAC_READER_BATCH_SIZE,
                progress_flush_every: thread_broker::DEFAULT_FLUSH_EVERY,
            }
        );
        assert_eq!(
            scatac_pipeline_geometry(
                crate::io::fastx::DecodeAllocation::PinnedAggregate {
                    slots: 2,
                    source: crate::io::fastx::DecodePinSource::AggregateEnvironment,
                },
                32,
                Some(17),
                Some(33),
            ),
            crate::io::fastx::PipelineTuning {
                reader_batch_size: 33,
                progress_flush_every: 17,
            }
        );
        assert_eq!(
            scatac_pipeline_geometry(
                crate::io::fastx::DecodeAllocation::PinnedPerFile {
                    slots_per_file: 2,
                    source: crate::io::fastx::DecodePinSource::CliPerFile,
                },
                32,
                None,
                None,
            ),
            crate::io::fastx::PipelineTuning {
                reader_batch_size: crate::io::fastx::SCATAC_READER_BATCH_SIZE,
                progress_flush_every: thread_broker::DEFAULT_FLUSH_EVERY,
            }
        );
    }

    #[cfg(feature = "rapidgzip")]
    #[test]
    fn scatac_defaults_to_sparse_responsive_monitoring() {
        let config = scatac_broker_config(8);
        assert_eq!(config.steady_probe_interval, Some(Duration::from_secs(5)));
        assert_eq!(
            config.opening_policy,
            thread_broker::OpeningPolicy::Bracket(Default::default())
        );
        assert_eq!(scatac_broker_config(32).min_producer_slots, 1);
        assert_eq!(scatac_broker_config(64).min_producer_slots, 1);
        assert_eq!(scatac_initial_decode_slots(8), 4);
        assert_eq!(scatac_initial_decode_slots(32), 4);
        assert_eq!(scatac_initial_decode_slots(64), 4);
    }
}
