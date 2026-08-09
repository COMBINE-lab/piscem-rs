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
use crate::io::fastx::{Collection, CollectionType, reader_with_batch_size};
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
    let thread_policy = crate::io::policy::ThreadPolicy::load(args.thread_policy.as_deref())?;
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
    let outcome = dispatch_on_k!(k, K => {
        if struct_constraints {
            run_scrna_pipeline::<K, SketchHitInfoChained, _, _>(
                &args.read1, &args.read2,
                &output_info, &stats,
                index, strat, opts, protocol, bc_len, umi_len,
                with_position, &read_length_samples,
                num_threads, decoder_pref, thread_policy, &progress,
            )?
        } else {
            run_scrna_pipeline::<K, SketchHitInfoSimple, _, _>(
                &args.read1, &args.read2,
                &output_info, &stats,
                index, strat, opts, protocol, bc_len, umi_len,
                with_position, &read_length_samples,
                num_threads, decoder_pref, thread_policy, &progress,
            )?
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
        skipping_strategy: &args.skipping_strategy,
    })?;

    Ok(())
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
    thread_policy: crate::io::policy::ThreadPolicy,
    progress: &indicatif::ProgressBar,
) -> Result<crate::io::fastx::PipelineOutcome>
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

    // The run's input as `paraseq` will see it: one group per logical record.
    // Both mates are listed because both consume decoder capacity even though
    // only the biological read is mapped.
    let groups: Vec<crate::io::calibrate::ReadGroup> = read1_paths
        .iter()
        .zip(read2_paths.iter())
        .map(|(a, b)| vec![a.clone(), b.clone()])
        .collect();

    // Ask the geometry which file carries mappable sequence rather than
    let decision = crate::io::calibrate::choose_decoder(
        &groups,
        num_threads,
        decoder_pref,
        thread_policy.parallel_decode,
    );
    #[cfg_attr(not(feature = "rapidgzip"), allow(unused_mut))]
    let mut plan = crate::io::fastx::plan_thread_budget(
        num_threads,
        read1_paths.len() * 2,
        decision.parallel,
        decoder_pref,
    )?;
    if plan.adaptive() {
        // Real scRNA and Flex normal-tier runs settle near one quarter of the
        // budget: 2/8, 8--9/32, and 13--16/64. Their generic mapping-biased
        // one-eighth opening therefore spends about two seconds growing into a
        // split already known to be safe. This call-site opening removes that
        // bounded startup tax without changing bulk modalities or fixed pins.
        plan.decode_slots = scrna_initial_decode_slots(plan.effective_budget);
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
    let mut pairs = Vec::with_capacity(read1_paths.len());
    for (r1_path, r2_path) in read1_paths.iter().zip(read2_paths.iter()) {
        #[cfg(feature = "rapidgzip")]
        let (o1, o2) = (
            crate::io::fastx::open_input_pooled(
                r1_path,
                decode_pool.as_ref(),
                plan.effective_budget,
            )?,
            crate::io::fastx::open_input_pooled(
                r2_path,
                decode_pool.as_ref(),
                plan.effective_budget,
            )?,
        );
        #[cfg(not(feature = "rapidgzip"))]
        let (o1, o2) = (
            crate::io::fastx::open_input(r1_path, 0)?,
            crate::io::fastx::open_input(r2_path, 0)?,
        );
        #[cfg(feature = "rapidgzip")]
        handles.extend([o1.handle.clone(), o2.handle.clone()].into_iter().flatten());
        let r1 = reader_with_batch_size(o1.reader)
            .map_err(|e| anyhow::anyhow!("failed to open {}: {}", r1_path.display(), e))?;
        let r2 = reader_with_batch_size(o2.reader)
            .map_err(|e| anyhow::anyhow!("failed to open {}: {}", r2_path.display(), e))?;
        pairs.push((r1, r2));
    }

    let mut readers = Vec::with_capacity(pairs.len() * 2);
    for (r1, r2) in pairs {
        readers.push(r1);
        readers.push(r2);
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
        crate::io::fastx::collection_share_floor(readers.len(), 2, plan.map_threads);

    #[cfg(feature = "rapidgzip")]
    let broker = match (&decode_pool, plan.adaptive()) {
        (Some(pool), true) => {
            match crate::io::broker::DecodeProducer::new(pool.clone(), handles.clone()) {
                Ok(producer) => {
                    let built = thread_broker::ThreadBroker::builder_with(
                        crate::io::broker::MappingConsumer::new(map_pool.clone(), stats),
                        producer,
                        crate::io::broker::broker_config_from_environment(
                            crate::io::fastx::broker_config_for_budget(plan.effective_budget),
                        )
                        .map_err(|e| {
                            anyhow::anyhow!("invalid thread broker probe interval: {e}")
                        })?,
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

    let collection = Collection::new(readers, CollectionType::Paired)
        .map_err(|e| anyhow::anyhow!("failed to create collection: {}", e))?;
    let result = collection
        .process_parallel_paired_pool(&mut processor, &map_pool, None)
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
        tuning: None,
        mapping_elapsed_secs: mapping_elapsed.as_secs_f64(),
    })
}

fn scrna_initial_decode_slots(budget: usize) -> usize {
    (budget / 4).clamp(1, budget.saturating_sub(1).max(1))
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

#[cfg(test)]
mod tests {
    use super::scrna_initial_decode_slots;

    #[test]
    fn adaptive_single_cell_opening_tracks_measured_quarter_budget() {
        assert_eq!(scrna_initial_decode_slots(2), 1);
        assert_eq!(scrna_initial_decode_slots(8), 2);
        assert_eq!(scrna_initial_decode_slots(32), 8);
        assert_eq!(scrna_initial_decode_slots(64), 16);
    }
}
