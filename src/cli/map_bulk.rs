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
use crate::io::fastx::{FastxConfig, FastxSource, ReadChunk};
use crate::io::map_info::write_map_info;
use crate::io::rad::{RadWriter, write_bulk_record, write_rad_header_bulk};
use crate::io::threads::{MappingStats, OutputInfo, ThreadConfig, run_mapping_pipeline};
use crate::mapping::cache::MappingCache;
use crate::mapping::filters::PoisonState;
use crate::mapping::hit_searcher::{HitSearcher, SkippingStrategy};
use crate::mapping::hits::MappingType;
use crate::mapping::map_fragment::{map_pe_fragment, map_se_fragment};
use crate::mapping::sketch_hit_simple::SketchHitInfoSimple;
use crate::mapping::streaming_query::PiscemStreamingQuery;
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
    let index = ReferenceIndex::load(index_prefix, true, !args.no_poison)?;
    info!("Index loaded: k={}, {} refs", index.k(), index.num_refs());

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

    // Setup pipeline
    let fastx_config = FastxConfig {
        read1_paths: args.read1.clone(),
        read2_paths: args.read2.clone(),
        chunk_size: 5000,
    };
    let fastx = FastxSource::new(fastx_config)?;

    let thread_config = ThreadConfig {
        threads: args.threads,
    };
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

    // Dispatch on K and run the pipeline
    dispatch_on_k!(k, K => {
        run_bulk_pipeline::<K>(
            fastx, thread_config, &output_info, &stats,
            &index, strat, is_paired, &end_cache,
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
    fastx: FastxSource,
    thread_config: ThreadConfig,
    output: &OutputInfo,
    stats: &MappingStats,
    index: &ReferenceIndex,
    strat: SkippingStrategy,
    is_paired: bool,
    end_cache: &UnitigEndCache,
) -> Result<()>
where
    Kmer<K>: KmerBits,
{
    let worker_fn = |chunk: ReadChunk, output: &OutputInfo, stats: &MappingStats| {
        // Per-thread state
        let mut hs = HitSearcher::new(index);
        let mut query = PiscemStreamingQuery::<K>::with_cache(index.dict(), end_cache);
        let mut cache_out = MappingCache::<SketchHitInfoSimple>::new(K);
        let mut cache_left = MappingCache::<SketchHitInfoSimple>::new(K);
        let mut cache_right = MappingCache::<SketchHitInfoSimple>::new(K);
        let mut poison_state = PoisonState::new(index.poison_table());
        let mut rad_writer = RadWriter::new();

        let mut local_reads: u64 = 0;
        let mut local_mapped: u64 = 0;
        let mut local_poisoned: u64 = 0;

        // Reserve space for chunk header: num_bytes (u32) + num_reads (u32)
        rad_writer.write_u32(0); // placeholder num_bytes
        rad_writer.write_u32(0); // placeholder num_reads
        let mut num_reads_in_chunk: u32 = 0;

        for read_pair in &chunk {
            local_reads += 1;

            if is_paired {
                let seq2 = read_pair.seq2.as_deref().unwrap_or(&[]);
                poison_state.paired_for_mapping = true;
                map_pe_fragment::<K>(
                    &read_pair.seq1,
                    seq2,
                    &mut hs,
                    &mut query,
                    &mut cache_left,
                    &mut cache_right,
                    &mut cache_out,
                    index,
                    &mut poison_state,
                    strat,
                );
            } else {
                poison_state.paired_for_mapping = false;
                map_se_fragment::<K>(
                    &read_pair.seq1,
                    &mut hs,
                    &mut query,
                    &mut cache_out,
                    index,
                    &mut poison_state,
                    strat,
                );
            }

            if poison_state.is_poisoned() {
                local_poisoned += 1;
                continue;
            }

            if cache_out.map_type != MappingType::Unmapped {
                local_mapped += 1;
                write_bulk_record(
                    cache_out.map_type,
                    &cache_out.accepted_hits,
                    &mut rad_writer,
                );
                num_reads_in_chunk += 1;
            }
        }

        // Backpatch chunk header
        let total_bytes = rad_writer.len() as u32;
        rad_writer.write_u32_at_offset(0, total_bytes);
        rad_writer.write_u32_at_offset(4, num_reads_in_chunk);

        // Flush to file under mutex
        if num_reads_in_chunk > 0 {
            let mut file = output.rad_file.lock().unwrap();
            rad_writer.flush_to(&mut *file).ok();
            output.num_chunks.fetch_add(1, Ordering::Relaxed);
        }

        stats.num_reads.fetch_add(local_reads, Ordering::Relaxed);
        stats.num_mapped.fetch_add(local_mapped, Ordering::Relaxed);
        stats
            .num_poisoned
            .fetch_add(local_poisoned, Ordering::Relaxed);
    };

    run_mapping_pipeline(fastx, thread_config, output, stats, worker_fn)
}
