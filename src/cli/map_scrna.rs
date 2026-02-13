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
use crate::io::fastx::{FastxConfig, FastxSource, ReadChunk};
use crate::io::map_info::write_map_info;
use crate::io::rad::{RadWriter, pack_bases_2bit, write_rad_header_sc, write_sc_record};
use crate::io::threads::{MappingStats, OutputInfo, ThreadConfig, run_mapping_pipeline};
use crate::mapping::cache::MappingCache;
use crate::mapping::filters::PoisonState;
use crate::mapping::hit_searcher::{HitSearcher, SkippingStrategy};
use crate::mapping::hits::MappingType;
use crate::mapping::map_fragment::{map_pe_fragment, map_se_fragment};
use crate::mapping::protocols::custom::parse_custom_geometry;
use crate::mapping::protocols::scrna::{ChromiumProtocol, count_ns, recover_barcode};
use crate::mapping::protocols::Protocol;
use crate::mapping::sketch_hit_simple::SketchHitInfoSimple;
use crate::mapping::streaming_query::PiscemStreamingQuery;
use crate::mapping::unitig_end_cache::UnitigEndCache;

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
    /// Include mapping positions in RAD output
    #[arg(long)]
    pub with_position: bool,
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

    // Shared read length samples (for --with-position backpatching)
    let read_length_samples: Mutex<Vec<u32>> = Mutex::new(Vec::new());

    let k = index.k();
    let with_position = args.with_position;
    let end_cache = UnitigEndCache::new(5_000_000);

    // Dispatch on K and run the pipeline
    dispatch_on_k!(k, K => {
        run_scrna_pipeline::<K>(
            fastx, thread_config, &output_info, &stats,
            &index, strat, protocol.as_ref(), bc_len, umi_len,
            with_position, &read_length_samples, &end_cache,
        )?;
    });

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
    fastx: FastxSource,
    thread_config: ThreadConfig,
    output: &OutputInfo,
    stats: &MappingStats,
    index: &ReferenceIndex,
    strat: SkippingStrategy,
    protocol: &dyn Protocol,
    bc_len: u16,
    umi_len: u16,
    with_position: bool,
    read_length_samples: &Mutex<Vec<u32>>,
    end_cache: &UnitigEndCache,
) -> Result<()>
where
    Kmer<K>: KmerBits,
{
    let is_bio_paired = protocol.is_bio_paired_end();

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
        let mut local_rlen_samples: Vec<u32> = Vec::new();
        let max_rlen_samples: usize = 10;

        // Reserve space for chunk header
        rad_writer.write_u32(0); // placeholder num_bytes
        rad_writer.write_u32(0); // placeholder num_reads
        let mut num_reads_in_chunk: u32 = 0;

        for read_pair in &chunk {
            local_reads += 1;

            let r1 = &read_pair.seq1;
            let r2 = read_pair.seq2.as_deref().unwrap_or(&[]);

            // Extract technical sequences
            let tech = protocol.extract_tech_seqs(r1, r2);
            let bc_raw = tech.barcode.unwrap_or(&[]);
            let umi_raw = tech.umi.unwrap_or(&[]);

            // Check barcode for N bases
            let n_count = count_ns(bc_raw);
            if n_count > 1 {
                continue; // can't recover >1 N
            }
            let recovered_bc = if n_count == 1 {
                recover_barcode(bc_raw)
            } else {
                None
            };
            let bc_to_pack = match &recovered_bc {
                Some(bc) => bc.as_slice(),
                None => bc_raw,
            };

            let bc_packed = pack_bases_2bit(bc_to_pack);
            let umi_packed = pack_bases_2bit(umi_raw);

            // Extract mappable reads
            let alignable = protocol.extract_mappable_reads(r1, r2);

            if is_bio_paired {
                let seq1 = alignable.seq1.unwrap_or(&[]);
                let seq2 = alignable.seq2.unwrap_or(&[]);
                if seq1.is_empty() && seq2.is_empty() {
                    continue;
                }
                poison_state.paired_for_mapping = true;
                map_pe_fragment::<K>(
                    seq1,
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

                // Sample read length from R2 for position backpatching
                if with_position && local_rlen_samples.len() < max_rlen_samples {
                    local_rlen_samples.push(seq2.len() as u32);
                }
            } else {
                let seq1 = alignable.seq1.unwrap_or(&[]);
                if seq1.is_empty() {
                    continue;
                }
                poison_state.paired_for_mapping = false;
                map_se_fragment::<K>(
                    seq1,
                    &mut hs,
                    &mut query,
                    &mut cache_out,
                    index,
                    &mut poison_state,
                    strat,
                );

                // Sample read length
                if with_position && local_rlen_samples.len() < max_rlen_samples {
                    local_rlen_samples.push(seq1.len() as u32);
                }
            }

            if poison_state.is_poisoned() {
                local_poisoned += 1;
                continue;
            }

            if cache_out.map_type != MappingType::Unmapped {
                local_mapped += 1;
                write_sc_record(
                    bc_packed,
                    umi_packed,
                    bc_len,
                    umi_len,
                    cache_out.map_type,
                    &cache_out.accepted_hits,
                    with_position,
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

        // Collect read length samples
        if with_position && !local_rlen_samples.is_empty() {
            let mut samples = read_length_samples.lock().unwrap();
            samples.extend_from_slice(&local_rlen_samples);
        }
    };

    run_mapping_pipeline(fastx, thread_config, output, stats, worker_fn)
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
