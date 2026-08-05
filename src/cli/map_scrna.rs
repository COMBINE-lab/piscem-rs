//! CLI command for single-cell RNA-seq mapping.

use std::io::{Seek, SeekFrom, Write};
use std::path::PathBuf;
use std::sync::Mutex;
use std::sync::atomic::Ordering;
use std::time::Instant;

use anyhow::{Context, Result};
use clap::Args;
use tracing::info;

use sshash_lib::{Kmer, KmerBits, KmerDictionary, dispatch_on_k};

use crate::index::contig_table::ContigTableLike;

use super::DictKind;
use crate::index::reference_index::{ReferenceIndex, tiny_artifacts_exist};
use crate::io::fastx::{
    Collection, CollectionType, reader_with_batch_size,
};
use crate::io::map_info::{MapInfoParams, write_map_info};
use crate::io::rad::{write_rad_header_sc, write_rad_header_sc_multi_bc};
use crate::io::threads::{MappingStats, OutputInfo};
use crate::mapping::chain_state::SketchHitInfoChained;
use crate::mapping::hit_searcher::SkippingStrategy;
use crate::mapping::hits::SketchHitInfo;
use crate::mapping::processors::{MappingOpts, ScrnaProcessor};
use crate::mapping::protocols::Protocol;
use crate::mapping::protocols::custom::parse_custom_geometry;
use crate::mapping::protocols::flex::ChromiumFlexProtocol;
use crate::mapping::protocols::scrna::ChromiumProtocol;
use crate::mapping::sketch_hit_simple::SketchHitInfoSimple;

use super::map_bulk::make_progress_bar;

#[derive(Args, Debug)]
pub struct MapScrnaArgs {
    /// Index prefix path
    #[arg(short = 'i', long)]
    pub index: PathBuf,
    /// Read 1 FASTQ files (comma-separated)
    #[arg(short = '1', long, value_delimiter = ',')]
    pub read1: Vec<PathBuf>,
    /// Read 2 FASTQ files (comma-separated)
    #[arg(short = '2', long, value_delimiter = ',')]
    pub read2: Vec<PathBuf>,
    /// Output directory
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
    /// Protocol geometry (e.g., chromium_v3, chromium_v2_5p)
    #[arg(short = 'g', long)]
    pub geometry: String,
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
    /// Include mapping positions in RAD output
    #[arg(long)]
    pub with_position: bool,
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

pub fn run(args: MapScrnaArgs) -> Result<()> {
    let start = Instant::now();

    // Validate skipping strategy early; full parse happens inside the generic runner.
    match args.skipping_strategy.to_lowercase().as_str() {
        "permissive" | "strict" => {}
        other => anyhow::bail!("unknown skipping strategy: {}", other),
    };

    // Parse protocol geometry: try built-in names first, then custom geometry
    let protocol: Box<dyn Protocol> =
        if let Some(chromium) = ChromiumProtocol::from_name(&args.geometry) {
            Box::new(chromium)
        } else if let Some(flex) = ChromiumFlexProtocol::from_name(&args.geometry) {
            Box::new(flex)
        } else {
            match parse_custom_geometry(&args.geometry) {
                Ok(custom) => Box::new(custom),
                Err(e) => anyhow::bail!(
                    "unknown geometry '{}' (not a built-in name, and custom parse failed: {})",
                    args.geometry,
                    e,
                ),
            }
        };
    if protocol.is_multi_barcode() {
        let descs = protocol.barcode_descs();
        let bc_desc_str: Vec<String> = descs
            .iter()
            .map(|d| format!("{}={}", d.tag_name, d.len))
            .collect();
        info!(
            "Protocol: {} (multi-barcode: [{}], umi_len={})",
            protocol.name(),
            bc_desc_str.join(", "),
            protocol.umi_len(),
        );
    } else {
        info!(
            "Protocol: {} (bc_len={}, umi_len={})",
            protocol.name(),
            protocol.barcode_len(),
            protocol.umi_len(),
        );
    }

    // --ignore-ambig-hits disables EC table loading
    let check_ambig = !args.ignore_ambig_hits;

    // Load index. Resolve `--dict auto` against on-disk artifacts. When the
    // effective dict is tiny and .tdct/.tct exist, load the Tiny-backed index
    // directly. Otherwise load the sshash index and, if tiny was explicitly
    // requested without prebuilt artifacts, convert in memory.
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
        info!(
            "Index loaded: k={}, {} refs ({:.2}s)",
            index.k(),
            index.num_refs(),
            load_start.elapsed().as_secs_f64()
        );
        return run_scrna_with_index(args, &index, start, protocol.as_ref());
    }
    info!("Loading index from {}", args.index.display());
    let index = ReferenceIndex::load(&args.index, check_ambig, !args.no_poison)?;
    info!(
        "Index loaded: k={}, {} refs ({:.2}s)",
        index.k(),
        index.num_refs(),
        load_start.elapsed().as_secs_f64()
    );

    // If --dict tiny was requested but no prebuilt artifacts exist, convert
    // in memory now and hand off to the same generic runner the on-disk path
    // uses. This keeps the conversion code path alive as a fallback.
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
        return run_scrna_with_index(args, &tiny_index, start, protocol.as_ref());
    }

    run_scrna_with_index(args, &index, start, protocol.as_ref())
}

#[allow(clippy::too_many_arguments)]
fn run_scrna_with_index<D, C>(
    args: MapScrnaArgs,
    index: &ReferenceIndex<D, C>,
    start: Instant,
    protocol: &dyn Protocol,
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

    // Create output directory and RAD file
    let out_dir = args.output.clone();
    std::fs::create_dir_all(&out_dir)
        .with_context(|| format!("failed to create output directory: {}", out_dir.display()))?;
    let rad_path = out_dir.join("map.rad");
    let mut rad_file = std::fs::File::create(&rad_path)
        .with_context(|| format!("failed to create {}", rad_path.display()))?;

    // Write SC RAD header (single or multi-barcode)
    let ref_names: Vec<&str> = (0..index.num_refs()).map(|i| index.ref_name(i)).collect();
    let bc_len = protocol.barcode_len() as u16;
    let umi_len = protocol.umi_len() as u16;
    let (chunk_count_offset, read_length_offset) = if protocol.is_multi_barcode() {
        let descs = protocol.barcode_descs();
        write_rad_header_sc_multi_bc(
            &mut rad_file,
            index.num_refs() as u64,
            &ref_names,
            &descs,
            umi_len,
            args.with_position,
        )?
    } else {
        write_rad_header_sc(
            &mut rad_file,
            index.num_refs() as u64,
            &ref_names,
            bc_len,
            umi_len,
            args.with_position,
        )?
    };

    // Create unmapped barcode count file with self-describing header.
    // The header encodes the number and widths of barcode fields so the
    // reader can parse records without hardcoding barcode widths.
    let unmapped_bc_path = out_dir.join("unmapped_bc_count.bin");
    let mut unmapped_bc_file = std::fs::File::create(&unmapped_bc_path)
        .with_context(|| format!("failed to create {}", unmapped_bc_path.display()))?;

    // Write the unmapped BC format header
    {
        use std::io::Write;
        let descs = protocol.barcode_descs();
        // version
        unmapped_bc_file.write_all(&[1u8])?;
        // num_fields
        unmapped_bc_file.write_all(&[descs.len() as u8])?;
        // field types: use the RAD tag type IDs based on barcode length in bases.
        // 2 bits per base: <=4bp fits u8, <=8bp fits u16, <=16bp fits u32, <=32bp fits u64
        for desc in &descs {
            let type_id: u8 = if desc.len <= 4 {
                1
            }
            // U8  (≤8 bits)
            else if desc.len <= 8 {
                2
            }
            // U16 (≤16 bits)
            else if desc.len <= 16 {
                3
            }
            // U32 (≤32 bits)
            else {
                4
            }; // U64 (≤64 bits)
            unmapped_bc_file.write_all(&[type_id])?;
        }
    }

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

    let read_length_samples: Mutex<Vec<u32>> = Mutex::new(Vec::new());

    // Setup progress bar
    let progress = make_progress_bar(args.quiet);

    let k = index.k();
    let with_position = args.with_position;
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

    // Capture values needed post-dispatch.
    let index_k = index.k();
    let index_m = index.m();
    let index_num_refs = index.num_refs();
    let sig_info_owned = index.ref_sig_info().cloned();

    // Dispatch on K and hit-info type, then run the pipeline via paraseq
    dispatch_on_k!(k, K => {
        if struct_constraints {
            run_scrna_pipeline::<K, SketchHitInfoChained, _, _>(
                &args.read1, &args.read2,
                &output_info, &stats,
                index, strat, opts, protocol, bc_len, umi_len,
                with_position, &read_length_samples,
                num_threads, decoder_pref, &progress,
            )?;
        } else {
            run_scrna_pipeline::<K, SketchHitInfoSimple, _, _>(
                &args.read1, &args.read2,
                &output_info, &stats,
                index, strat, opts, protocol, bc_len, umi_len,
                with_position, &read_length_samples,
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

    // Backpatch read length if --with-position
    if let Some(rlen_offset) = read_length_offset {
        let samples = read_length_samples.lock().unwrap();
        if !samples.is_empty() {
            let mode = compute_mode(&samples);
            let mode_count = samples.iter().filter(|&&v| v == mode).count();
            let total = samples.len();
            if mode_count as f64 / total as f64 >= 0.9 {
                rad_file.seek(SeekFrom::Start(rlen_offset))?;
                rad_file.write_all(&mode.to_le_bytes())?;
                info!(
                    "Read length mode: {} ({}/{} samples)",
                    mode, mode_count, total
                );
            } else {
                info!(
                    "Warning: read lengths are heterogeneous (mode {} = {}/{})",
                    mode, mode_count, total,
                );
            }
        }
    }
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
        mode: "sc-rna",
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
fn run_scrna_pipeline<
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
    protocol: &dyn Protocol,
    bc_len: u16,
    umi_len: u16,
    with_position: bool,
    read_length_samples: &Mutex<Vec<u32>>,
    num_threads: usize,
    decoder_pref: crate::io::calibrate::DecoderPreference,
    progress: &indicatif::ProgressBar,
) -> Result<()>
where
    Kmer<K>: KmerBits,
{
    let mut processor = ScrnaProcessor::<K, S, D, C>::new(
        index,
        None,
        output,
        stats,
        strat,
        opts,
        protocol,
        bc_len,
        umi_len,
        with_position,
        read_length_samples,
        progress,
    );

    // Split the user's thread budget between mapping and decompression.
    //
    // `num_threads` here is already the mapping allotment; `decode_budget`
    // is what the decoders may collectively spend on top of it. Both derive
    // from the same user-supplied `-t`, see `plan_thread_budget`.
    let mut plan = crate::io::fastx::plan_thread_budget(num_threads, read1_paths.len() * 2);
    // Tier 1: when nothing forces the choice, measure instead of assuming.
    if let Some(decision) = calibrate_decoder::<K, S, D, C>(
        // Ask the geometry which file carries mappable sequence rather than
        // assuming read 2. It is read 2 in every current chemistry, but a
        // custom geometry says where `r` is, and timing the technical read
        // would measure the wrong sequence.
        if protocol.bio_read_file() == 0 {
            read1_paths.first()
        } else {
            read2_paths.first()
        },
        read1_paths.len() * 2,
        num_threads,
        index,
        &opts,
        strat,
        decoder_pref,
    ) {
        if !decision.parallel {
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
            plan.decode_budget = w.saturating_mul(read1_paths.len() * 2).max(1);
        }
    }

    let mut pairs = Vec::with_capacity(read1_paths.len());
    // Only exists on the rapidgzip path; without it there is nothing to
    // supervise and the type would be uninferable.
    #[cfg(feature = "rapidgzip")]
    let mut handles: Vec<rapidgzip_core::DecoderHandle> = Vec::new();
    for (r1_path, r2_path) in read1_paths.iter().zip(read2_paths.iter()) {
        let o1 = crate::io::fastx::open_input(r1_path, plan.per_file_ceiling, plan.initial_per_file)?;
        let o2 = crate::io::fastx::open_input(r2_path, plan.per_file_ceiling, plan.initial_per_file)?;
        #[cfg(feature = "rapidgzip")]
        handles.extend([o1.handle.clone(), o2.handle.clone()].into_iter().flatten());
        let r1 = reader_with_batch_size(o1.reader)
            .map_err(|e| anyhow::anyhow!("failed to open {}: {}", r1_path.display(), e))?;
        let r2 = reader_with_batch_size(o2.reader)
            .map_err(|e| anyhow::anyhow!("failed to open {}: {}", r2_path.display(), e))?;
        pairs.push((r1, r2));
    }

    // Supervise the decoders under one shared budget. Upstream adapts within a
    // decoder; this bounds the sum across them, which is the part it does not do.
    #[cfg(feature = "rapidgzip")]
    let budget = crate::io::decode_budget::DecodeBudget::spawn(handles, plan.decode_budget);

    let mut readers = Vec::with_capacity(pairs.len() * 2);
    for (r1, r2) in pairs {
        readers.push(r1);
        readers.push(r2);
    }
    let collection = Collection::new(readers, CollectionType::Paired)
        .map_err(|e| anyhow::anyhow!("failed to create collection: {}", e))?;
    let result = collection
        .process_parallel_paired(&mut processor, num_threads, None)
        .map_err(|e| anyhow::anyhow!("mapping failed: {}", e));

    #[cfg(feature = "rapidgzip")]
    if let Some(budget) = budget {
        let report = budget.finish();
        tracing::info!(
            "decoder threads: peak {} worker + {} auxiliary (budget {}, unused {}); \
             peak busy {}; calibration settled at {:?}",
            report.peak_worker_threads,
            report.peak_auxiliary_threads,
            plan.decode_budget,
            report.unused_budget,
            report.peak_busy_workers,
            report.converged_workers,
        );
    }

    result
}

/// Compute the mode (most frequent value) of a u32 slice.
fn compute_mode(values: &[u32]) -> u32 {
    let mut counts = std::collections::HashMap::new();
    for &v in values {
        *counts.entry(v).or_insert(0u64) += 1;
    }
    counts
        .into_iter()
        .max_by_key(|&(_, count)| count)
        .map(|(val, _)| val)
        .unwrap_or(0)
}
