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
use std::time::Instant;

use anyhow::{Context, Result};
use clap::Args;
use tracing::info;

use sshash_lib::{Kmer, KmerBits, KmerDictionary, dispatch_on_k};

use crate::index::contig_table::ContigTableLike;

use super::DictKind;
use crate::index::reference_index::{ReferenceIndex, tiny_artifacts_exist};
use crate::io::fastx::{Collection, CollectionType, reader_with_batch_size};
use crate::io::map_info::{MapInfoParams, write_map_info};
use crate::io::rad::write_rad_header_atac;
use crate::io::threads::{MappingStats, OutputInfo};
use crate::mapping::binning::BinPos;
use crate::mapping::processors::{MappingOpts, ScatacProcessor};
use crate::mapping::unitig_end_cache::UnitigEndCache;

use super::map_bulk::make_progress_bar;

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
    /// Number of mapping threads
    #[arg(short = 't', long, default_value = "16")]
    pub threads: usize,

    /// Gzip decoder selection: `auto`, `serial`, `parallel`, or `parallel=N`.
    ///
    /// `auto` (the default) decides from the input: forced rules first, then a
    /// brief measurement of how starved the mapping threads are. Override it
    /// when you know something the probe cannot — a slow network filesystem, a
    /// shared node where spending cores on decode is antisocial, or to reproduce
    /// a measurement. `parallel` still yields on inputs that are not regular
    /// files, where the parallel decoder degrades to sequential anyway.
    #[arg(long, default_value = "auto", value_name = "MODE")]
    pub decoder: String,
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

    dispatch_on_k!(k, K => {
        let (bio_paths, r2_paths) = if is_paired {
            (args.read1.as_slice(), args.read2.as_slice())
        } else {
            (args.reads.as_slice(), [].as_slice())
        };
        run_atac_pipeline::<K, _, _>(
            bio_paths, &args.barcode, r2_paths,
            &output_info, &stats,
            index, &binning, bc_len, tn5_shift, min_overlap, Some(&end_cache),
            is_paired, opts,
            num_threads, decoder_pref, &progress,
        )?;
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
        sig_info: sig_info_owned.as_ref(),
        piscem_rs_version: crate::VERSION,
        num_threads,
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
    opts: MappingOpts,
    num_threads: usize,
    decoder_pref: crate::io::calibrate::DecoderPreference,
    progress: &indicatif::ProgressBar,
) -> Result<()>
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

    // scATAC opens `arity` files per sample -- biological read, barcode read,
    // and mate when paired -- so a multi-sample run already has many inflate
    // streams and is usually served by serial decode alone. Measured on 10x
    // `atac_pbmc_5k` against a k=23 chr1+chr2 index at `-t 32`, the parallel
    // decoder cost +1.0% wall and +1.1% CPU: the supervisor found no starvation
    // to answer. Calibrating anyway is what makes that a measured outcome for
    // each run rather than an assumption carried from one dataset.
    let groups: Vec<crate::io::calibrate::ReadGroup> = (0..bio_paths.len())
        .map(|i| {
            let mut g = vec![bio_paths[i].clone(), barcode_paths[i].clone()];
            if is_paired {
                g.push(read2_paths[i].clone());
            }
            g
        })
        .collect();
    let decision = crate::io::calibrate::choose_decoder(&groups, num_threads, decoder_pref);
    let mut plan = crate::io::fastx::plan_thread_budget(num_threads, num_input_files);
    if !decision.parallel {
        plan.parallel_gzip = false;
        plan.decode_budget = 0;
    } else if let crate::io::calibrate::DecoderPreference::Parallel {
        workers_per_file: Some(w),
    } = decoder_pref
    {
        // A named number is a request to *use* it, not to converge toward it.
        plan.decode_budget = w
            .saturating_mul(num_input_files)
            .max(1)
            .min(num_threads - 1);
        // The budget is the point. `decode_budget` was overwritten without
        // touching `map_threads`, so a pinned request spent both sides of the
        // split in full: measured at `-t 32`, `--decoder parallel=16` ran an
        // average of 41 concurrent threads and looked 41% faster than `auto`
        // purely by using half as much machine again.
        plan.map_threads = num_threads.saturating_sub(plan.decode_budget).max(1);
    }
    #[cfg(feature = "rapidgzip")]
    let pinned = matches!(
        decoder_pref,
        crate::io::calibrate::DecoderPreference::Parallel {
            workers_per_file: Some(_)
        }
    );

    // One pool for the whole run, sized to the entire budget: `workers` is an
    // immutable maximum and a refused `set_worker_limit` would silently desync
    // the broker's accounting from the pool's.
    #[cfg(feature = "rapidgzip")]
    let decode_pool = if plan.parallel_gzip {
        rapidgzip_core::DecoderPool::builder()
            .workers(num_threads.max(2))
            .initial_worker_limit(plan.decode_budget.clamp(1, num_threads - 1))
            .build()
            .ok()
    } else {
        None
    };
    let num_threads = plan.map_threads;
    #[cfg(feature = "rapidgzip")]
    let mut handles: Vec<rapidgzip_core::DecoderHandle> = Vec::new();

    // `mut` is required only with the `rapidgzip` feature, where the closure
    // pushes into `handles`; without it the body has nothing to mutate.
    #[allow(unused_mut)]
    let mut open = |path: &std::path::Path| -> Result<_> {
        #[cfg(feature = "rapidgzip")]
        let opened = crate::io::fastx::open_input_pooled(path, decode_pool.as_ref(), num_threads)?;
        #[cfg(not(feature = "rapidgzip"))]
        let opened = crate::io::fastx::open_input(path, 0)?;
        #[cfg(feature = "rapidgzip")]
        if let Some(h) = opened.handle {
            handles.push(h);
        }
        reader_with_batch_size(opened.reader)
            .map_err(|e| anyhow::anyhow!("failed to open {}: {}", path.display(), e))
    };

    // Scoped so the closure's mutable borrow of `handles` ends here, leaving
    // the vector free for the supervisor below.
    let mut readers = Vec::with_capacity(num_input_files);
    {
        for i in 0..bio_paths.len() {
            readers.push(open(&bio_paths[i])?);
            readers.push(open(&barcode_paths[i])?);
            if is_paired {
                readers.push(open(&read2_paths[i])?);
            }
        }
    }

    // The mapping side, resizable for the same reason.
    let map_pool = paraseq::parallel::ThreadPool::with_max(plan.map_threads, num_threads);

    #[cfg(feature = "rapidgzip")]
    let broker = match (&decode_pool, pinned) {
        (Some(pool), false) => thread_broker::ThreadBroker::builder(
            crate::io::broker::MappingConsumer::new(map_pool.clone(), stats),
            crate::io::broker::DecodeProducer::new(pool.clone(), handles.clone()),
        )
        .budget(num_threads)
        .initial_producer_slots(plan.decode_budget)
        .build()
        .map(|b| b.start())
        // Loud, not swallowed. A misconfigured broker silently returning `None`
        // means the run quietly keeps whatever split it opened with, which is
        // exactly how a self-inconsistent default pair went unnoticed: the
        // decode limit stayed at its opening value and consumer utilisation
        // fell from 93% to 32% with nothing said.
        .inspect_err(|e| tracing::warn!("thread broker disabled: {e}"))
        .ok(),
        _ => None,
    };

    let collection = Collection::new(readers, CollectionType::Multi { arity })
        .map_err(|e| anyhow::anyhow!("failed to create collection: {}", e))?;
    collection
        .process_parallel_multi_pool(&mut processor, &map_pool, None)
        .map_err(|e| anyhow::anyhow!("mapping failed: {}", e))?;

    #[cfg(feature = "rapidgzip")]
    if let Some(broker) = broker {
        let r = broker.finish();
        tracing::info!(
            "thread broker: settled at {} mapping / {} decode of {num_threads} \
             ({} moves, {} reverts, {} resurveys, converged {}, {:?} to settle)",
            r.final_consumer_threads,
            r.final_producer_limit,
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

    Ok(())
}
