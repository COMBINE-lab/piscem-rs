//! CLI command for bulk RNA-seq mapping.

use std::io::{Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::Ordering;
use std::time::{Duration, Instant};

use anyhow::{Context, Result};
use clap::Args;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use tracing::info;

use sshash_lib::{Kmer, KmerBits, dispatch_on_k};

use crate::index::reference_index::ReferenceIndex;
use crate::io::fastx::{Collection, CollectionType, open_with_decompression};
use crate::io::map_info::write_map_info;
use crate::io::rad::write_rad_header_bulk;
use crate::io::threads::{MappingStats, OutputInfo};
use crate::mapping::hit_searcher::SkippingStrategy;
use crate::mapping::processors::{BulkProcessor, MappingOpts};

#[derive(Args, Debug)]
#[command(group(
    clap::ArgGroup::new("input_reads")
        .required(true)
        .args(["reads", "read1"])
))]
pub struct MapBulkArgs {
    /// Index prefix path
    #[arg(short = 'i', long)]
    pub index: String,
    /// Single-end FASTQ files (comma-separated); mutually exclusive with -1/-2
    #[arg(short = 'r', long = "reads", value_delimiter = ',',
          conflicts_with_all = ["read1", "read2"])]
    pub reads: Vec<String>,
    /// Read 1 FASTQ files (comma-separated); requires -2
    #[arg(short = '1', long, value_delimiter = ',',
          requires = "read2", conflicts_with = "reads")]
    pub read1: Vec<String>,
    /// Read 2 FASTQ files (comma-separated); requires -1
    #[arg(short = '2', long, value_delimiter = ',',
          requires = "read1", conflicts_with = "reads")]
    pub read2: Vec<String>,
    /// Output directory
    #[arg(short = 'o', long)]
    pub output: String,
    /// Number of mapping threads
    #[arg(short = 't', long, default_value = "16")]
    pub threads: usize,
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
    /// Suppress progress output
    #[arg(short = 'q', long)]
    pub quiet: bool,
}

pub fn run(args: MapBulkArgs) -> Result<()> {
    let start = Instant::now();

    let strat = match args.skipping_strategy.to_lowercase().as_str() {
        "permissive" => SkippingStrategy::Permissive,
        "strict" => SkippingStrategy::Strict,
        other => anyhow::bail!("unknown skipping strategy: {}", other),
    };

    let is_paired = !args.read1.is_empty();
    info!(
        "Mapping {} reads ({})",
        if is_paired { "paired-end" } else { "single-end" },
        args.skipping_strategy,
    );

    // --ignore-ambig-hits disables EC table loading
    let check_ambig = !args.ignore_ambig_hits;

    // Load index
    let index_prefix = Path::new(&args.index);
    info!("Loading index from {}", index_prefix.display());
    let load_start = Instant::now();
    let index = ReferenceIndex::load(index_prefix, check_ambig, !args.no_poison)?;
    let load_secs = load_start.elapsed().as_secs_f64();
    info!(
        "Index loaded: k={}, {} refs ({:.2}s)",
        index.k(),
        index.num_refs(),
        load_secs
    );

    // Create output directory and RAD file
    let out_dir = PathBuf::from(&args.output);
    std::fs::create_dir_all(&out_dir)
        .with_context(|| format!("failed to create output directory: {}", out_dir.display()))?;
    let rad_path = out_dir.join("map.rad");
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
    let opts = MappingOpts {
        max_hit_occ: args.max_hit_occ,
        max_hit_occ_recover: args.max_hit_occ_recover,
        max_read_occ: args.max_read_occ,
        max_ec_card: if args.ignore_ambig_hits { 0 } else { args.max_ec_card },
    };

    // Dispatch on K and run the pipeline via paraseq
    dispatch_on_k!(k, K => {
        let (r1_paths, r2_paths) = if is_paired {
            (args.read1.as_slice(), args.read2.as_slice())
        } else {
            (args.reads.as_slice(), [].as_slice())
        };
        run_bulk_pipeline::<K>(
            r1_paths, r2_paths,
            &output_info, &stats,
            &index, strat, opts, is_paired,
            num_threads, &progress,
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
    write_map_info(
        &out_dir.join("map_info.json"),
        num_reads,
        num_mapped,
        num_poisoned,
        elapsed,
    )?;

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn run_bulk_pipeline<const K: usize>(
    read1_paths: &[String],
    read2_paths: &[String],
    output: &OutputInfo,
    stats: &MappingStats,
    index: &ReferenceIndex,
    strat: SkippingStrategy,
    opts: MappingOpts,
    is_paired: bool,
    num_threads: usize,
    progress: &ProgressBar,
) -> Result<()>
where
    Kmer<K>: KmerBits,
{
    let mut processor = BulkProcessor::<K>::new(index, None, output, stats, strat, opts, progress);

    if is_paired {
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
    } else {
        let mut readers = Vec::with_capacity(read1_paths.len());
        for r1_path in read1_paths {
            readers.push(
                paraseq::fastx::Reader::new(open_with_decompression(r1_path)?)
                    .map_err(|e| anyhow::anyhow!("failed to open {}: {}", r1_path, e))?,
            );
        }
        let collection = Collection::new(readers, CollectionType::Single)
            .map_err(|e| anyhow::anyhow!("failed to create collection: {}", e))?;
        collection
            .process_parallel(&mut processor, num_threads, None)
            .map_err(|e| anyhow::anyhow!("mapping failed: {}", e))?;
    }

    Ok(())
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
