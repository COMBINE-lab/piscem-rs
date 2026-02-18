//! CLI command for single-cell RNA-seq mapping.

use std::io::{Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::Ordering;
use std::sync::Mutex;
use std::time::Instant;

use anyhow::{Context, Result};
use clap::Args;
use tracing::info;

use sshash_lib::{Kmer, KmerBits, dispatch_on_k};

use crate::index::reference_index::ReferenceIndex;
use crate::io::fastx::{open_with_decompression, Collection, CollectionType};
use crate::io::map_info::write_map_info;
use crate::io::rad::write_rad_header_sc;
use crate::io::threads::{MappingStats, OutputInfo};
use crate::mapping::hit_searcher::SkippingStrategy;
use crate::mapping::processors::{MappingOpts, ScrnaProcessor};
use crate::mapping::protocols::custom::parse_custom_geometry;
use crate::mapping::protocols::scrna::ChromiumProtocol;
use crate::mapping::protocols::Protocol;

use super::map_bulk::make_progress_bar;

#[derive(Args, Debug)]
pub struct MapScrnaArgs {
    /// Index prefix path
    #[arg(short = 'i', long)]
    pub index: String,
    /// Read 1 FASTQ files (comma-separated)
    #[arg(short = '1', long, value_delimiter = ',')]
    pub read1: Vec<String>,
    /// Read 2 FASTQ files (comma-separated)
    #[arg(short = '2', long, value_delimiter = ',')]
    pub read2: Vec<String>,
    /// Output directory
    #[arg(short = 'o', long)]
    pub output: String,
    /// Number of mapping threads
    #[arg(short = 't', long, default_value = "16")]
    pub threads: usize,
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
    /// Suppress progress output
    #[arg(short = 'q', long)]
    pub quiet: bool,
}

pub fn run(args: MapScrnaArgs) -> Result<()> {
    let start = Instant::now();

    let strat = match args.skipping_strategy.to_lowercase().as_str() {
        "permissive" => SkippingStrategy::Permissive,
        "strict" => SkippingStrategy::Strict,
        other => anyhow::bail!("unknown skipping strategy: {}", other),
    };

    // Parse protocol geometry: try built-in names first, then custom geometry
    let protocol: Box<dyn Protocol> =
        if let Some(chromium) = ChromiumProtocol::from_name(&args.geometry) {
            Box::new(chromium)
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
    info!(
        "Protocol: {} (bc_len={}, umi_len={})",
        protocol.name(),
        protocol.barcode_len(),
        protocol.umi_len(),
    );

    // --ignore-ambig-hits disables EC table loading
    let check_ambig = !args.ignore_ambig_hits;

    // Load index
    let index_prefix = Path::new(&args.index);
    info!("Loading index from {}", index_prefix.display());
    let index = ReferenceIndex::load(index_prefix, check_ambig, !args.no_poison)?;
    info!("Index loaded: k={}, {} refs", index.k(), index.num_refs());

    // Create output directory and RAD file
    let out_dir = PathBuf::from(&args.output);
    std::fs::create_dir_all(&out_dir)
        .with_context(|| format!("failed to create output directory: {}", out_dir.display()))?;
    let rad_path = out_dir.join("map.rad");
    let mut rad_file = std::fs::File::create(&rad_path)
        .with_context(|| format!("failed to create {}", rad_path.display()))?;

    // Write SC RAD header
    let ref_names: Vec<&str> = (0..index.num_refs()).map(|i| index.ref_name(i)).collect();
    let bc_len = protocol.barcode_len() as u16;
    let umi_len = protocol.umi_len() as u16;
    let (chunk_count_offset, read_length_offset) = write_rad_header_sc(
        &mut rad_file,
        index.num_refs() as u64,
        &ref_names,
        bc_len,
        umi_len,
        args.with_position,
    )?;

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

    let read_length_samples: Mutex<Vec<u32>> = Mutex::new(Vec::new());

    // Setup progress bar
    let progress = make_progress_bar(args.quiet);

    let k = index.k();
    let with_position = args.with_position;
    let num_threads = args.threads.max(1);
    let opts = MappingOpts {
        max_hit_occ: args.max_hit_occ,
        max_hit_occ_recover: args.max_hit_occ_recover,
        max_read_occ: args.max_read_occ,
        max_ec_card: if args.ignore_ambig_hits { 0 } else { args.max_ec_card },
    };

    // Dispatch on K and run the pipeline via paraseq
    dispatch_on_k!(k, K => {
        run_scrna_pipeline::<K>(
            &args.read1, &args.read2,
            &output_info, &stats,
            &index, strat, opts, protocol.as_ref(), bc_len, umi_len,
            with_position, &read_length_samples,
            num_threads, &progress,
        )?;
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
                info!("Read length mode: {} ({}/{} samples)", mode, mode_count, total);
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
    let cmdline = format!(
        "piscem-rs map-scrna -i {} -g {} -o {}",
        args.index, args.geometry, args.output
    );
    write_map_info(
        &out_dir.join("map_info.json"),
        num_reads,
        num_mapped,
        num_poisoned,
        &cmdline,
        elapsed,
    )?;

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn run_scrna_pipeline<const K: usize>(
    read1_paths: &[String],
    read2_paths: &[String],
    output: &OutputInfo,
    stats: &MappingStats,
    index: &ReferenceIndex,
    strat: SkippingStrategy,
    opts: MappingOpts,
    protocol: &dyn Protocol,
    bc_len: u16,
    umi_len: u16,
    with_position: bool,
    read_length_samples: &Mutex<Vec<u32>>,
    num_threads: usize,
    progress: &indicatif::ProgressBar,
) -> Result<()>
where
    Kmer<K>: KmerBits,
{
    let mut processor = ScrnaProcessor::<K>::new(
        index, None, output, stats, strat, opts, protocol, bc_len, umi_len, with_position,
        read_length_samples, progress,
    );

    let mut readers = Vec::with_capacity(read1_paths.len() * 2);
    for (r1_path, r2_path) in read1_paths.iter().zip(read2_paths.iter()) {
        readers.push(
            paraseq::fastx::Reader::new(open_with_decompression(r1_path)?)
                .map_err(|e| anyhow::anyhow!("failed to open {}: {}", r1_path, e))?,
        );
        readers.push(
            paraseq::fastx::Reader::new(open_with_decompression(r2_path)?)
                .map_err(|e| anyhow::anyhow!("failed to open {}: {}", r2_path, e))?,
        );
    }
    let collection = Collection::new(readers, CollectionType::Paired)
        .map_err(|e| anyhow::anyhow!("failed to create collection: {}", e))?;
    collection
        .process_parallel_paired(&mut processor, num_threads, None)
        .map_err(|e| anyhow::anyhow!("mapping failed: {}", e))?;

    Ok(())
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
