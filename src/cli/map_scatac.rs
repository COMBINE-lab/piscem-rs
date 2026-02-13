//! CLI command for scATAC-seq mapping.
//!
//! scATAC reads are always paired-end with a separate barcode file.
//! Supports mate overlap detection (merged mapping) and optional Tn5 shift.

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
use crate::io::rad::{RadWriter, pack_bases_2bit, write_atac_record, write_rad_header_atac};
use crate::io::threads::{MappingStats, OutputInfo, ThreadConfig, run_mapping_pipeline};
use crate::mapping::cache::MappingCache;
use crate::mapping::filters::PoisonState;
use crate::mapping::hit_searcher::{HitSearcher, SkippingStrategy};
use crate::mapping::hits::MappingType;
use crate::mapping::map_fragment::{map_pe_fragment, map_se_fragment};
use crate::mapping::overlap::{find_overlap, OverlapType};
use crate::mapping::protocols::scrna::{barcode_has_n, count_ns, recover_barcode};
use crate::mapping::sketch_hit_simple::SketchHitInfoSimple;
use crate::mapping::streaming_query::PiscemStreamingQuery;
use crate::mapping::unitig_end_cache::UnitigEndCache;

#[derive(Args, Debug)]
pub struct MapScatacArgs {
    /// Index prefix path
    #[arg(short = 'i', long)]
    pub index: String,
    /// Read 1 FASTQ files (comma-separated)
    #[arg(short = '1', long, value_delimiter = ',')]
    pub read1: Vec<String>,
    /// Read 2 FASTQ files (comma-separated)
    #[arg(short = '2', long, value_delimiter = ',')]
    pub read2: Vec<String>,
    /// Barcode FASTQ files (comma-separated)
    #[arg(short = 'b', long, value_delimiter = ',')]
    pub barcode: Vec<String>,
    /// Output directory
    #[arg(short = 'o', long)]
    pub output: String,
    /// Number of mapping threads
    #[arg(short = 't', long, default_value = "16")]
    pub threads: usize,
    /// Barcode length in bases
    #[arg(long, default_value = "16")]
    pub bc_len: usize,
    /// Disable Tn5 transposase shift
    #[arg(long)]
    pub no_tn5_shift: bool,
    /// Disable poison k-mer filtering
    #[arg(long)]
    pub no_poison: bool,
    /// K-mer skipping strategy (permissive or strict)
    #[arg(long, default_value = "permissive")]
    pub skipping_strategy: String,
    /// Minimum overlap length for mate merging
    #[arg(long, default_value = "30")]
    pub min_overlap: i32,
}

pub fn run(args: MapScatacArgs) -> Result<()> {
    let start = Instant::now();

    let strat = match args.skipping_strategy.to_lowercase().as_str() {
        "permissive" => SkippingStrategy::Permissive,
        "strict" => SkippingStrategy::Strict,
        other => anyhow::bail!("unknown skipping strategy: {}", other),
    };

    info!("scATAC mapping (bc_len={}, tn5_shift={})",
        args.bc_len,
        !args.no_tn5_shift,
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

    // Setup pipeline
    // Note: scATAC needs paired reads + barcode file.
    // We pass read1 and read2 as paired, barcode is extracted inline.
    // For now, barcodes are passed as read2 of a second paired source,
    // but the simplest approach is to pass barcode files as a third
    // channel. Since our FastxSource supports paired mode, we use R1+R2
    // for biological reads and handle barcodes separately.
    //
    // TODO: Full barcode file support requires extending FastxSource for
    // triple-file input. For now, we accept paired R1+R2 only and the
    // barcode is expected as the first bc_len bases of R1 (fallback mode).
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
    let tn5_shift = !args.no_tn5_shift;
    let min_overlap = args.min_overlap;
    let end_cache = UnitigEndCache::new(5_000_000);

    // Dispatch on K and run the pipeline
    dispatch_on_k!(k, K => {
        run_atac_pipeline::<K>(
            fastx, thread_config, &output_info, &stats,
            &index, strat, bc_len, tn5_shift, min_overlap, &end_cache,
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
    let cmdline = format!(
        "piscem-rs map-scatac -i {} -o {}",
        args.index, args.output
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
fn run_atac_pipeline<const K: usize>(
    fastx: FastxSource,
    thread_config: ThreadConfig,
    output: &OutputInfo,
    stats: &MappingStats,
    index: &ReferenceIndex,
    strat: SkippingStrategy,
    bc_len: u16,
    tn5_shift: bool,
    min_overlap: i32,
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

        // Reserve space for chunk header
        rad_writer.write_u32(0); // placeholder num_bytes
        rad_writer.write_u32(0); // placeholder num_reads
        let mut num_reads_in_chunk: u32 = 0;

        for read_pair in &chunk {
            local_reads += 1;

            let r1 = &read_pair.seq1;
            let r2 = read_pair.seq2.as_deref().unwrap_or(&[]);

            if r1.is_empty() || r2.is_empty() {
                continue;
            }

            // For scATAC, the barcode comes from a separate file.
            // In this implementation, we extract it from the first bc_len
            // bases of R1 as a fallback. A full implementation would use
            // a triple-file reader.
            let bc_end = (bc_len as usize).min(r1.len());
            let bc_raw = &r1[..bc_end];

            // Check barcode for Ns
            let n_count = count_ns(bc_raw);
            if n_count > 1 {
                continue;
            }
            let recovered_bc = if n_count == 1 {
                recover_barcode(bc_raw)
            } else {
                None
            };
            let bc_to_pack = match &recovered_bc {
                Some(bc) => bc.as_slice(),
                None => {
                    if barcode_has_n(bc_raw) {
                        continue;
                    }
                    bc_raw
                }
            };
            let bc_packed = pack_bases_2bit(bc_to_pack);

            // Biological reads: all of R1 (after barcode) and R2
            let bio_r1 = &r1[bc_end..];
            let bio_r2 = r2;

            if bio_r1.is_empty() && bio_r2.is_empty() {
                continue;
            }

            poison_state.clear();
            cache_out.clear();

            // Try mate overlap detection
            let mate_ov = find_overlap(bio_r1, bio_r2, min_overlap, 0);

            if mate_ov.ov_type != OverlapType::NoOverlap && !mate_ov.frag.is_empty() {
                // Map merged fragment as single read
                poison_state.paired_for_mapping = false;
                map_se_fragment::<K>(
                    &mate_ov.frag,
                    &mut hs,
                    &mut query,
                    &mut cache_out,
                    index,
                    &mut poison_state,
                    strat,
                );

                // Set fragment length from overlap
                if cache_out.map_type != MappingType::Unmapped {
                    for hit in &mut cache_out.accepted_hits {
                        hit.fragment_length = mate_ov.frag_length as i32;
                    }
                }
            } else {
                // No overlap: map both ends independently, then merge
                poison_state.paired_for_mapping = true;
                map_pe_fragment::<K>(
                    bio_r1,
                    bio_r2,
                    &mut hs,
                    &mut query,
                    &mut cache_left,
                    &mut cache_right,
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
                write_atac_record(
                    bc_packed,
                    bc_len,
                    cache_out.map_type,
                    &cache_out.accepted_hits,
                    tn5_shift,
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
        }

        stats.num_reads.fetch_add(local_reads, Ordering::Relaxed);
        stats.num_mapped.fetch_add(local_mapped, Ordering::Relaxed);
        stats
            .num_poisoned
            .fetch_add(local_poisoned, Ordering::Relaxed);
    };

    run_mapping_pipeline(fastx, thread_config, output, stats, worker_fn)
}
