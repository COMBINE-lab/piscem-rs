//! CLI command for bulk RNA-seq mapping.

use std::io::{Seek, SeekFrom, Write};
use std::path::PathBuf;
use std::sync::atomic::Ordering;
use std::time::{Duration, Instant};

use anyhow::{Context, Result};
use clap::Args;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use tracing::info;

use sshash_lib::{Kmer, KmerBits, KmerDictionary, dispatch_on_k};

use crate::index::contig_table::ContigTableLike;

use super::DictKind;
use crate::index::reference_index::{ReferenceIndex, tiny_artifacts_exist};
use crate::io::fastx::{
    Collection, CollectionType, reader_with_batch_size,
};
use crate::io::map_info::{MapInfoParams, write_map_info};
use crate::io::rad::write_rad_header_bulk;
use crate::io::threads::{MappingStats, OutputInfo};
use crate::mapping::chain_state::SketchHitInfoChained;
use crate::mapping::hit_searcher::SkippingStrategy;
use crate::mapping::hits::SketchHitInfo;
use crate::mapping::processors::{BulkProcessor, MappingOpts};
use crate::mapping::sketch_hit_simple::SketchHitInfoSimple;

#[derive(Args, Debug)]
#[command(group(
    clap::ArgGroup::new("input_reads")
        .required(true)
        .args(["reads", "read1"])
))]
pub struct MapBulkArgs {
    /// Index prefix path
    #[arg(short = 'i', long)]
    pub index: PathBuf,
    /// Single-end FASTQ files (comma-separated); mutually exclusive with -1/-2
    #[arg(short = 'r', long = "reads", value_delimiter = ',',
          conflicts_with_all = ["read1", "read2"])]
    pub reads: Vec<PathBuf>,
    /// Read 1 FASTQ files (comma-separated); requires -2
    #[arg(
        short = '1',
        long,
        value_delimiter = ',',
        requires = "read2",
        conflicts_with = "reads"
    )]
    pub read1: Vec<PathBuf>,
    /// Read 2 FASTQ files (comma-separated); requires -1
    #[arg(
        short = '2',
        long,
        value_delimiter = ',',
        requires = "read1",
        conflicts_with = "reads"
    )]
    pub read2: Vec<PathBuf>,
    /// Output file stem (e.g. foo/bar/sample); creates foo/bar/sample.rad and foo/bar/sample.map_info.json
    #[arg(short = 'o', long)]
    pub output: PathBuf,
    /// Number of mapping threads
    #[arg(short = 't', long, default_value = "16")]
    pub threads: usize,
    /// Gzip decoder selection: `auto`, `serial`, `parallel`, or `parallel=N`.
    ///
    /// `auto` (the default) decides from the input: forced rules first, then a
    /// brief measurement of decode and mapping rates. Override it when you know
    /// something the probe cannot — a slow network filesystem, a shared node
    /// where spending cores on decode is antisocial, or to reproduce a
    /// measurement. `parallel` still yields on inputs that are not regular
    /// files, where the parallel decoder degrades to sequential anyway.
    #[arg(long, default_value = "auto", value_name = "MODE")]
    pub decoder: String,
    /// K-mer skipping strategy (permissive or strict)
    #[arg(long, default_value = "permissive")]
    pub skipping_strategy: String,
    /// Disable poison k-mer filtering
    #[arg(long)]
    pub no_poison: bool,
    /// Ignore highly-ambiguous hits rather than using EC-based fallback;
    /// mutually exclusive with --max-ec-card
    #[arg(long, conflicts_with = "max_ec_card")]
    pub ignore_ambig_hits: bool,
    /// Maximum equivalence class cardinality for ambiguous hit resolution
    #[arg(long, default_value = "4096", conflicts_with = "ignore_ambig_hits")]
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
    /// Apply structural constraints: require k-mers to form positionally
    /// coherent chains (max stretch 31 bp, up to 8 chains per orientation).
    /// Equivalent to the C++ `-c`/`--struct-constraints` flag.
    #[arg(short = 'c', long)]
    pub struct_constraints: bool,
    /// Suppress progress output
    #[arg(short = 'q', long)]
    pub quiet: bool,
    /// K-mer dictionary backend: `sshash` (compact, default) or `tiny`
    /// (hashbrown-backed, faster but higher memory).
    #[arg(long, value_enum, default_value_t = DictKind::Auto)]
    pub dict: DictKind,
}

pub fn run(args: MapBulkArgs) -> Result<()> {
    let start = Instant::now();

    match args.skipping_strategy.to_lowercase().as_str() {
        "permissive" | "strict" => {}
        other => anyhow::bail!("unknown skipping strategy: {}", other),
    };

    let is_paired = !args.read1.is_empty();
    info!(
        "Mapping {} reads ({})",
        if is_paired {
            "paired-end"
        } else {
            "single-end"
        },
        args.skipping_strategy,
    );

    // --ignore-ambig-hits disables EC table loading
    let check_ambig = !args.ignore_ambig_hits;

    // Resolve `--dict auto` against on-disk artifacts. When effective dict is
    // tiny and .tdct/.tct exist, load the Tiny index directly. Otherwise load
    // sshash and optionally convert to Tiny in memory.
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
        >::load_tiny(&args.index, check_ambig, !args.no_poison)?;
        let load_secs = load_start.elapsed().as_secs_f64();
        info!(
            "Index loaded: k={}, {} refs ({:.2}s)",
            index.k(),
            index.num_refs(),
            load_secs
        );
        return run_bulk_with_index(args, &index, start, is_paired, load_secs);
    }
    info!("Loading index from {}", args.index.display());
    let index = ReferenceIndex::load(&args.index, check_ambig, !args.no_poison)?;
    let load_secs = load_start.elapsed().as_secs_f64();
    info!(
        "Index loaded: k={}, {} refs ({:.2}s)",
        index.k(),
        index.num_refs(),
        load_secs
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
        return run_bulk_with_index(args, &tiny_index, start, is_paired, load_secs);
    }

    run_bulk_with_index(args, &index, start, is_paired, load_secs)
}

#[allow(clippy::too_many_arguments)]
fn run_bulk_with_index<D, C>(
    args: MapBulkArgs,
    index: &ReferenceIndex<D, C>,
    start: Instant,
    is_paired: bool,
    load_secs: f64,
) -> Result<()>
where
    D: KmerDictionary + Sync,
    C: ContigTableLike + Sync,
{
    let strat = match args.skipping_strategy.to_lowercase().as_str() {
        "permissive" => SkippingStrategy::Permissive,
        "strict" => SkippingStrategy::Strict,
        other => anyhow::bail!("unknown skipping strategy: {}", other),
    };

    // Treat -o as a file stem: create parent dirs, then append extensions
    if let Some(parent) = args.output.parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("failed to create output directory: {}", parent.display()))?;
    }
    let mut rad_path = args.output.clone();
    rad_path.add_extension("rad");
    let mut rad_file = std::fs::File::create(&rad_path)
        .with_context(|| format!("failed to create {}", rad_path.display()))?;

    // Write bulk RAD header
    let ref_names: Vec<&str> = (0..index.num_refs()).map(|i| index.ref_name(i)).collect();
    let ref_lengths: Vec<u32> = (0..index.num_refs())
        .map(|i| index.ref_len(i) as u32)
        .collect();
    let chunk_count_offset = write_rad_header_bulk(
        &mut rad_file,
        is_paired,
        index.num_refs() as u64,
        &ref_names,
        &ref_lengths,
    )?;

    // Setup shared output state
    let stats = MappingStats::new();
    let output_info = OutputInfo {
        num_chunks: std::sync::atomic::AtomicUsize::new(0),
        rad_file: std::sync::Mutex::new(std::io::BufWriter::new(
            rad_file
                .try_clone()
                .context("failed to clone RAD file handle")?,
        )),
        unmapped_bc_file: None,
    };

    // Setup progress bar
    let progress = make_progress_bar(args.quiet);

    let k = index.k();
    let num_threads = args.threads.max(1);
    let decoder_pref = crate::io::calibrate::DecoderPreference::parse(&args.decoder)
        .map_err(|e| anyhow::anyhow!("invalid --decoder value: {e}"))?;
    let opts = MappingOpts {
        max_hit_occ: args.max_hit_occ,
        max_hit_occ_recover: args.max_hit_occ_recover,
        max_read_occ: args.max_read_occ,
        max_ec_card: if args.ignore_ambig_hits {
            0
        } else {
            args.max_ec_card
        },
    };
    let struct_constraints = args.struct_constraints;

    let index_k = index.k();
    let index_m = index.m();
    let index_num_refs = index.num_refs();
    let sig_info_owned = index.ref_sig_info().cloned();

    dispatch_on_k!(k, K => {
        let (r1_paths, r2_paths) = if is_paired {
            (args.read1.as_slice(), args.read2.as_slice())
        } else {
            (args.reads.as_slice(), [].as_slice())
        };
        if struct_constraints {
            run_bulk_pipeline::<K, SketchHitInfoChained, _, _>(
                r1_paths, r2_paths,
                &output_info, &stats,
                index, strat, opts, is_paired,
                num_threads, decoder_pref, &progress,
            )?;
        } else {
            run_bulk_pipeline::<K, SketchHitInfoSimple, _, _>(
                r1_paths, r2_paths,
                &output_info, &stats,
                index, strat, opts, is_paired,
                num_threads, decoder_pref, &progress,
            )?;
        }
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

    let mapping_secs = elapsed - load_secs;
    info!(
        "Mapped {}/{} reads ({:.1}%), {} poisoned, {:.2}s total ({:.2}s index load + {:.2}s mapping)",
        num_mapped,
        num_reads,
        if num_reads > 0 {
            num_mapped as f64 / num_reads as f64 * 100.0
        } else {
            0.0
        },
        num_poisoned,
        elapsed,
        load_secs,
        mapping_secs,
    );

    // Write map_info.json
    let mut map_info_path = args.output.clone();
    map_info_path.add_extension("map_info.json");
    write_map_info(&MapInfoParams {
        path: &map_info_path,
        mode: "bulk",
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
        skipping_strategy: &args.skipping_strategy,
    })?;

    Ok(())
}

/// Run the Tier 1 calibration probe when nothing forces the decoder choice.
///
/// Maps a small prefix of the first input on one thread through the real
/// kernel, so the measured per-thread rate is the rate this index and these
/// reads actually achieve — not a guess from index size, which two data points
/// could not support.
///
/// Returns `None` when the choice was already forced, when there is no input,
/// or when the probe itself fails; every one of those falls back to the
/// ratio rule in `plan_thread_budget`, which is what ran before this existed.
#[allow(clippy::too_many_arguments)]
fn calibrate_decoder<const K: usize, S: SketchHitInfo, D: KmerDictionary, C: ContigTableLike>(
    first_path: Option<&std::path::PathBuf>,
    num_files: usize,
    num_threads: usize,
    index: &ReferenceIndex<D, C>,
    opts: &crate::mapping::processors::MappingOpts,
    strat: SkippingStrategy,
    pref: crate::io::calibrate::DecoderPreference,
) -> Option<crate::io::calibrate::Decision>
where
    Kmer<K>: KmerBits,
{
    use crate::io::calibrate;

    let path = first_path?;
    let kind = calibrate::classify_input(path);
    // Tell the user when their flag cannot be carried out, rather than quietly
    // doing something else.
    let compression = calibrate::detect_compression(path, kind);
    if let Some(conflict) = calibrate::preference_conflict(pref, kind, compression) {
        tracing::warn!("{}", conflict);
    }

    // An explicit request outranks measurement, but not the forcings: a
    // preference cannot make a pipe seekable.
    if let Some(chosen) = calibrate::preference_choice(pref, kind) {
        tracing::info!(
            "decoder selected by request: {} ({:?})",
            if chosen.parallel { "parallel" } else { "serial" },
            chosen.reason
        );
        return Some(chosen);
    }
    if let Some(forced) = calibrate::forced_choice(kind, num_files, num_threads) {
        tracing::debug!("decoder choice forced: {:?}", forced.reason);
        return Some(forced);
    }

    // Ambiguous: measure. One thread's worth of mapping state, matching what a
    // worker would build.
    let mut query =
        crate::mapping::streaming_query::PiscemStreamingQuery::<K, D>::new(index.dict());
    let mut hs = crate::mapping::hit_searcher::HitSearcher::new(index);
    let mut cache = crate::mapping::cache::MappingCache::<S>::new(K);
    opts.apply_to(&mut cache);
    let mut poison = crate::mapping::filters::PoisonState::new(index.poison_table());

    let rates = calibrate::probe_path(path, calibrate::DEFAULT_PROBE_RECORDS, |seq| {
        crate::mapping::engine::map_read::<K, S, D, C>(
            seq, &mut cache, &mut hs, &mut query, index, &mut poison, strat,
        );
    })
    .ok()
    .flatten()?;

    let decision = calibrate::decide(&rates, num_files, num_threads, calibrate::DEFAULT_MARGIN);
    tracing::info!(
        "decoder calibration: {} records, {:.3} GB/s per stream, {:.3} GB/s per mapping thread -> {} ({:?})",
        rates.records_sampled,
        rates.decode_gbps_per_stream,
        rates.map_gbps_per_thread,
        if decision.parallel { "parallel" } else { "serial" },
        decision.reason
    );
    Some(decision)
}


#[allow(clippy::too_many_arguments)]
fn run_bulk_pipeline<
    const K: usize,
    S: SketchHitInfo + Send + 'static,
    D: KmerDictionary + Sync,
    C: ContigTableLike + Sync,
>(
    read1_paths: &[PathBuf],
    read2_paths: &[PathBuf],
    output: &OutputInfo,
    stats: &MappingStats,
    index: &ReferenceIndex<D, C>,
    strat: SkippingStrategy,
    opts: MappingOpts,
    is_paired: bool,
    num_threads: usize,
    decoder_pref: crate::io::calibrate::DecoderPreference,
    progress: &ProgressBar,
) -> Result<()>
where
    Kmer<K>: KmerBits,
{
    let mut processor =
        BulkProcessor::<K, S, D, C>::new(index, None, output, stats, strat, opts, progress);

    // Same decode/map budget split as the scRNA path; see `plan_thread_budget`.
    let num_input_files = if is_paired {
        read1_paths.len() * 2
    } else {
        read1_paths.len()
    };
    let mut plan = crate::io::fastx::plan_thread_budget(num_threads, num_input_files);
    // Tier 1: when nothing forces the choice, measure instead of assuming.
    if let Some(decision) = calibrate_decoder::<K, S, D, C>(
        read1_paths.first(),
        num_input_files,
        num_threads,
        index,
        &opts,
        strat,
        decoder_pref,
    ) {
        // Only an explicit user request may force serial. The measured
        // supply/demand comparison is logged but does NOT override the
        // threads-per-file rule: validated against the crossover surface it got
        // 4 of 8 points wrong, every one of them predicting serial where the
        // parallel decoder measurably won (by up to 1.92x), while the ratio
        // rule got 7 of 8. See `calibrate` for why margin tuning cannot fix it.
        if !decision.parallel
            && decision.reason != crate::io::calibrate::Reason::MeasuredConsumerBound
        {
            plan.parallel_gzip = false;
            plan.decode_budget = 0;
            plan.per_file_ceiling = 0;
            plan.initial_per_file = 0;
        } else if let crate::io::calibrate::DecoderPreference::Parallel {
            workers_per_file: Some(w),
        } = decoder_pref
        {
            // An explicit worker count is a ceiling *and* a starting point:
            // someone naming a number wants it used, not ratcheted up to.
            plan.parallel_gzip = true;
            plan.per_file_ceiling = w;
            plan.initial_per_file = w;
            plan.decode_budget = w.saturating_mul(num_input_files).max(1);
        }
    }
    // Only exists on the rapidgzip path; without it there is nothing to
    // supervise and the type would be uninferable.
    #[cfg(feature = "rapidgzip")]
    let mut handles: Vec<rapidgzip_core::DecoderHandle> = Vec::new();

    let mut readers = Vec::with_capacity(num_input_files);
    if is_paired {
        for (r1_path, r2_path) in read1_paths.iter().zip(read2_paths.iter()) {
            for path in [r1_path, r2_path] {
                let o = crate::io::fastx::open_input(
                    path,
                    plan.per_file_ceiling,
                    plan.initial_per_file,
                )?;
                #[cfg(feature = "rapidgzip")]
                handles.extend(o.handle.clone());
                readers.push(
                    reader_with_batch_size(o.reader).map_err(|e| {
                        anyhow::anyhow!("failed to open {}: {}", path.display(), e)
                    })?,
                );
            }
        }
    } else {
        for r1_path in read1_paths {
            let o = crate::io::fastx::open_input(
                r1_path,
                plan.per_file_ceiling,
                plan.initial_per_file,
            )?;
            #[cfg(feature = "rapidgzip")]
            handles.extend(o.handle.clone());
            readers.push(
                reader_with_batch_size(o.reader)
                    .map_err(|e| anyhow::anyhow!("failed to open {}: {}", r1_path.display(), e))?,
            );
        }
    }

    #[cfg(feature = "rapidgzip")]
    let budget = crate::io::decode_budget::DecodeBudget::spawn(handles, plan.decode_budget);

    let result = if is_paired {
        Collection::new(readers, CollectionType::Paired)
            .map_err(|e| anyhow::anyhow!("failed to create collection: {}", e))
            .and_then(|c| {
                c.process_parallel_paired(&mut processor, num_threads, None)
                    .map_err(|e| anyhow::anyhow!("mapping failed: {}", e))
            })
    } else {
        Collection::new(readers, CollectionType::Single)
            .map_err(|e| anyhow::anyhow!("failed to create collection: {}", e))
            .and_then(|c| {
                c.process_parallel(&mut processor, num_threads, None)
                    .map_err(|e| anyhow::anyhow!("mapping failed: {}", e))
            })
    };

    #[cfg(feature = "rapidgzip")]
    if let Some(budget) = budget {
        let report = budget.finish();
        tracing::info!(
            "decoder threads: peak {} worker + {} auxiliary (budget {}); peak busy {}; \
             calibration settled at {:?}",
            report.peak_worker_threads,
            report.peak_auxiliary_threads,
            plan.decode_budget,
            report.peak_busy_workers,
            report.converged_workers,
        );
    }

    result
}

/// Create a progress bar for mapping (shared across all CLI commands).
pub(crate) fn make_progress_bar(quiet: bool) -> ProgressBar {
    if quiet {
        ProgressBar::hidden()
    } else {
        let pb = ProgressBar::new_spinner();
        pb.set_draw_target(ProgressDrawTarget::stderr_with_hz(1));
        pb.set_style(
            ProgressStyle::with_template(
                "{spinner:.green} [{elapsed_precise}] {human_pos} reads processed ({per_sec})",
            )
            .unwrap(),
        );
        pb.enable_steady_tick(Duration::from_millis(1_000));
        pb
    }
}
