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
use crate::io::fastx::{Collection, CollectionType, reader_with_batch_size};
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

    // The run's input as `paraseq` will see it: one group per logical record.
    // `choose_decoder` needs groups rather than a flat list so it can pick a
    // group that is safe to probe and scale the probe's thread share itself.
    let groups: Vec<crate::io::calibrate::ReadGroup> = if is_paired {
        read1_paths
            .iter()
            .zip(read2_paths.iter())
            .map(|(a, b)| vec![a.clone(), b.clone()])
            .collect()
    } else {
        read1_paths.iter().map(|a| vec![a.clone()]).collect()
    };

    let decision = crate::io::calibrate::choose_decoder(&groups, num_threads, decoder_pref);
    let mut plan = crate::io::fastx::plan_thread_budget(num_threads, num_input_files);
    if !decision.parallel {
        plan.parallel_gzip = false;
        plan.decode_budget = 0;
    } else if let crate::io::calibrate::DecoderPreference::Parallel {
        workers_per_file: Some(w),
    } = decoder_pref
    {
        // A named number is a request to *use* it, not to converge toward it,
        // so it pins the pool and the broker stays out of the way.
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

    // One pool for the whole run, sized to the *entire* budget.
    //
    // `workers` is an immutable maximum and `set_worker_limit` is refused above
    // it, so sizing it to the expected split would let the broker grant threads
    // the pool then refuses -- and since the broker tracks what it asked for
    // rather than reading it back, its accounting would silently diverge.
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

    #[cfg(feature = "rapidgzip")]
    let mut handles: Vec<rapidgzip_core::DecoderHandle> = Vec::new();

    let mut readers = Vec::with_capacity(num_input_files);
    {
        let mut paths: Vec<&PathBuf> = Vec::with_capacity(num_input_files);
        if is_paired {
            for (a, b) in read1_paths.iter().zip(read2_paths.iter()) {
                paths.push(a);
                paths.push(b);
            }
        } else {
            paths.extend(read1_paths.iter());
        }
        for path in paths {
            #[cfg(feature = "rapidgzip")]
            let o = crate::io::fastx::open_input_pooled(path, decode_pool.as_ref(), num_threads)?;
            #[cfg(not(feature = "rapidgzip"))]
            let o = crate::io::fastx::open_input(path, 0)?;
            #[cfg(feature = "rapidgzip")]
            handles.extend(o.handle.clone());
            readers.push(
                reader_with_batch_size(o.reader)
                    .map_err(|e| anyhow::anyhow!("failed to open {}: {}", path.display(), e))?,
            );
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

    let busy_at_start = stats.busy.nanos();
    let wall_start = std::time::Instant::now();

    let result = if is_paired {
        Collection::new(readers, CollectionType::Paired)
            .map_err(|e| anyhow::anyhow!("failed to create collection: {}", e))
            .and_then(|c| {
                c.process_parallel_paired_pool(&mut processor, &map_pool, None)
                    .map_err(|e| anyhow::anyhow!("mapping failed: {}", e))
            })
    } else {
        Collection::new(readers, CollectionType::Single)
            .map_err(|e| anyhow::anyhow!("failed to create collection: {}", e))
            .and_then(|c| {
                c.process_parallel_pool(&mut processor, &map_pool, None)
                    .map_err(|e| anyhow::anyhow!("mapping failed: {}", e))
            })
    };

    {
        let busy = stats.busy.nanos().saturating_sub(busy_at_start);
        let wall = wall_start.elapsed().as_nanos() as u64;
        let cap = wall.saturating_mul(num_threads.max(1) as u64);
        tracing::info!(
            "consumer meter: {:.2} thread-seconds mapping over {:.2} s wall x {} threads              -> {:.1}% utilisation ({:.1}% idle)",
            busy as f64 / 1e9,
            wall as f64 / 1e9,
            num_threads,
            if cap > 0 {
                busy as f64 * 100.0 / cap as f64
            } else {
                0.0
            },
            if cap > 0 {
                100.0 - busy as f64 * 100.0 / cap as f64
            } else {
                0.0
            },
        );
    }

    #[cfg(feature = "rapidgzip")]
    if let Some(broker) = broker {
        let r = broker.finish();
        tracing::info!(
            "thread broker: settled at {} mapping / {} decode of {num_threads} \
             ({} moves, {} reverts, {} resurveys, converged {}, {:?} to settle{}{})",
            r.final_consumer_threads,
            r.final_producer_limit,
            r.moves,
            r.reverts,
            r.resurveys,
            r.converged(),
            r.time_to_converge,
            if r.source_bound_samples > 0 {
                format!(", {} source-bound samples", r.source_bound_samples)
            } else {
                String::new()
            },
            if r.inelastic_samples > 0 {
                format!(", {} inelastic samples", r.inelastic_samples)
            } else {
                String::new()
            },
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
