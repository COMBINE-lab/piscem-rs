//! CLI command for bulk RNA-seq mapping.

use std::io::{Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::Ordering;
use std::time::Instant;

use anyhow::{Context, Result};
use clap::Args;
use tracing::info;

use sshash_lib::{Kmer, KmerBits, dispatch_on_k};

use crate::index::reference_index::ReferenceIndex;
use crate::io::fastx::open_concatenated_readers;
use crate::io::map_info::write_map_info;
use crate::io::rad::write_rad_header_bulk;
use crate::io::threads::{MappingStats, OutputInfo};
use crate::mapping::hit_searcher::SkippingStrategy;
use crate::mapping::processors::BulkProcessor;
use crate::mapping::unitig_end_cache::UnitigEndCache;

#[derive(Args, Debug)]
pub struct MapBulkArgs {
    /// Index prefix path
    #[arg(short = 'i', long)]
    pub index: String,
    /// Read 1 FASTQ files (comma-separated)
    #[arg(short = '1', long, value_delimiter = ',')]
    pub read1: Vec<String>,
    /// Read 2 FASTQ files (comma-separated, omit for single-end)
    #[arg(short = '2', long, value_delimiter = ',')]
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
}

pub fn run(args: MapBulkArgs) -> Result<()> {
    let start = Instant::now();

    let strat = match args.skipping_strategy.to_lowercase().as_str() {
        "permissive" => SkippingStrategy::Permissive,
        "strict" => SkippingStrategy::Strict,
        other => anyhow::bail!("unknown skipping strategy: {}", other),
    };

    let is_paired = !args.read2.is_empty();
    info!(
        "Mapping {} reads ({})",
        if is_paired { "paired-end" } else { "single-end" },
        args.skipping_strategy,
    );

    // Load index
    let index_prefix = Path::new(&args.index);
    info!("Loading index from {}", index_prefix.display());
    let load_start = Instant::now();
    let index = ReferenceIndex::load(index_prefix, true, !args.no_poison)?;
    let load_secs = load_start.elapsed().as_secs_f64();
    info!("Index loaded: k={}, {} refs ({:.2}s)", index.k(), index.num_refs(), load_secs);

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
    };

    let k = index.k();
    let end_cache = UnitigEndCache::new(5_000_000);
    let num_threads = args.threads.max(1);

    // Dispatch on K and run the pipeline via paraseq
    dispatch_on_k!(k, K => {
        run_bulk_pipeline::<K>(
            &args.read1, &args.read2,
            &output_info, &stats,
            &index, strat, is_paired, &end_cache,
            num_threads,
        )?;
    });

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
    let cmdline = format!("piscem-rs map-bulk -i {} -o {}", args.index, args.output);
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
fn run_bulk_pipeline<const K: usize>(
    read1_paths: &[String],
    read2_paths: &[String],
    output: &OutputInfo,
    stats: &MappingStats,
    index: &ReferenceIndex,
    strat: SkippingStrategy,
    is_paired: bool,
    end_cache: &UnitigEndCache,
    num_threads: usize,
) -> Result<()>
where
    Kmer<K>: KmerBits,
{
    use paraseq::parallel::ParallelReader;

    let mut processor = BulkProcessor::<K>::new(index, end_cache, output, stats, strat);

    if is_paired {
        let r1 = open_concatenated_readers(read1_paths)?;
        let r2 = open_concatenated_readers(read2_paths)?;
        let reader1 = paraseq::fastq::Reader::new(r1);
        let reader2 = paraseq::fastq::Reader::new(r2);
        reader1
            .process_parallel_paired(reader2, &mut processor, num_threads)
            .map_err(|e| anyhow::anyhow!("mapping failed: {}", e))?;
    } else {
        let r1 = open_concatenated_readers(read1_paths)?;
        let reader1 = paraseq::fastq::Reader::new(r1);
        reader1
            .process_parallel(&mut processor, num_threads)
            .map_err(|e| anyhow::anyhow!("mapping failed: {}", e))?;
    }

    Ok(())
}
