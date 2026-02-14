//! CLI command for scATAC-seq mapping.
//!
//! scATAC reads use triple-file input: R1 (genomic), barcode, R2 (genomic).
//! Supports mate overlap detection (merged mapping), bin-based hit collection,
//! and optional Tn5 shift.
//!
//! Key differences from scRNA mapping:
//! - Uses `get_raw_hits_sketch_everykmer` (not STRICT/PERMISSIVE skipping)
//! - Bin-based hit accumulation with threshold filtering
//! - Bin-aware paired-end merge
//! - Poison filtering disabled by default (`--no-poison` is true)

use std::io::{Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::Ordering;
use std::time::Instant;

use anyhow::{Context, Result};
use clap::Args;
use tracing::info;

use sshash_lib::{Kmer, KmerBits, dispatch_on_k};

use crate::index::reference_index::ReferenceIndex;
use crate::io::fastx::{FastxTripleSource, ReadTripletChunk};
use crate::io::map_info::write_map_info;
use crate::io::rad::{RadWriter, pack_bases_2bit, write_atac_record, write_rad_header_atac};
use crate::io::threads::{MappingStats, OutputInfo, ThreadConfig, run_mapping_pipeline_triple};
use crate::mapping::binning::BinPos;
use crate::mapping::cache::MappingCache;
use crate::mapping::filters::PoisonState;
use crate::mapping::hit_searcher::HitSearcher;
use crate::mapping::hits::MappingType;
use crate::mapping::map_fragment::{map_pe_fragment_atac, map_se_fragment_atac};
use crate::mapping::merge_pairs::{remove_duplicate_hits_pub, simple_hit_cmp_bins};
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
    /// Read 1 FASTQ files (genomic left, comma-separated)
    #[arg(short = '1', long, value_delimiter = ',')]
    pub read1: Vec<String>,
    /// Read 2 FASTQ files (genomic right, comma-separated)
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
    /// Disable poison k-mer filtering (default: true, matching C++ behavior)
    #[arg(long, default_value = "true")]
    pub no_poison: bool,
    /// Bin size for genomic binning
    #[arg(long, default_value = "1000")]
    pub bin_size: u64,
    /// Overlap between adjacent bins
    #[arg(long, default_value = "300")]
    pub bin_overlap: u64,
    /// Hit threshold fraction for bin-based filtering
    #[arg(long, default_value = "0.7")]
    pub thr: f32,
    /// Minimum overlap length for mate merging
    #[arg(long, default_value = "30")]
    pub min_overlap: i32,
    /// K-mer skipping strategy (ignored for ATAC — always uses every-kmer)
    #[arg(long)]
    pub skipping_strategy: Option<String>,
}

pub fn run(args: MapScatacArgs) -> Result<()> {
    let start = Instant::now();

    if args.skipping_strategy.is_some() {
        info!("Note: --skipping-strategy is ignored for scATAC (always uses every-kmer mode)");
    }

    info!(
        "scATAC mapping (bc_len={}, tn5_shift={}, bin_size={}, overlap={}, thr={:.2})",
        args.bc_len,
        !args.no_tn5_shift,
        args.bin_size,
        args.bin_overlap,
        args.thr,
    );

    // Load index
    let index_prefix = Path::new(&args.index);
    info!("Loading index from {}", index_prefix.display());
    // C++ ATAC: check_ambig_hits defaults to false, so EC table is NOT loaded
    let index = ReferenceIndex::load(index_prefix, false, !args.no_poison)?;
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

    // Create binning scheme
    let binning = BinPos::new(&index, args.bin_size, args.bin_overlap, args.thr);
    info!("Binning: {} total bins", binning.num_bins());

    // Setup triple-file FASTQ source
    let fastx = FastxTripleSource::new(
        &args.read1,
        &args.barcode,
        &args.read2,
        5000,
        false,
    )?;

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
            &index, &binning, bc_len, tn5_shift, min_overlap, &end_cache,
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
    fastx: FastxTripleSource,
    thread_config: ThreadConfig,
    output: &OutputInfo,
    stats: &MappingStats,
    index: &ReferenceIndex,
    binning: &BinPos,
    bc_len: u16,
    tn5_shift: bool,
    min_overlap: i32,
    end_cache: &UnitigEndCache,
) -> Result<()>
where
    Kmer<K>: KmerBits,
{
    let worker_fn = |chunk: ReadTripletChunk, output: &OutputInfo, stats: &MappingStats| {
        // Per-thread state
        let mut hs = HitSearcher::new(index);
        let mut query = PiscemStreamingQuery::<K>::with_cache(index.dict(), end_cache);
        let mut cache_out = MappingCache::<SketchHitInfoSimple>::new(K);
        let mut cache_left = MappingCache::<SketchHitInfoSimple>::new(K);
        let mut cache_right = MappingCache::<SketchHitInfoSimple>::new(K);
        let mut poison_state = PoisonState::new(index.poison_table());
        let mut rad_writer = RadWriter::with_capacity(150_000);

        let mut local_reads: u64 = 0;
        let mut local_mapped: u64 = 0;
        let mut local_poisoned: u64 = 0;

        // Reserve space for chunk header
        rad_writer.write_u32(0); // placeholder num_bytes
        rad_writer.write_u32(0); // placeholder num_reads
        let mut num_reads_in_chunk: u32 = 0;

        for triplet in &chunk {
            local_reads += 1;

            let r1 = &triplet.seq1;
            let r2 = &triplet.seq2;
            let bc_raw = &triplet.barcode;

            if r1.is_empty() || r2.is_empty() {
                continue;
            }

            // Check barcode for Ns (truncate to bc_len if longer)
            let bc_end = (bc_len as usize).min(bc_raw.len());
            let bc_slice = &bc_raw[..bc_end];

            let n_count = count_ns(bc_slice);
            if n_count > 1 {
                continue;
            }
            let recovered_bc = if n_count == 1 {
                recover_barcode(bc_slice)
            } else {
                None
            };
            let bc_to_pack = match &recovered_bc {
                Some(bc) => bc.as_slice(),
                None => {
                    if barcode_has_n(bc_slice) {
                        continue;
                    }
                    bc_slice
                }
            };
            let bc_packed = pack_bases_2bit(bc_to_pack);

            poison_state.clear();
            cache_out.clear();

            // Try mate overlap detection
            let mate_ov = find_overlap(r1, r2, min_overlap, 0);

            if mate_ov.ov_type != OverlapType::NoOverlap && !mate_ov.frag.is_empty() {
                // Map merged fragment as single read (bin-based)
                poison_state.paired_for_mapping = false;
                map_se_fragment_atac::<K>(
                    &mate_ov.frag,
                    &mut hs,
                    &mut query,
                    &mut cache_out,
                    index,
                    &mut poison_state,
                    binning,
                );

                // Sort+dedup hits (matching C++ pesc_sc_atac.cpp lines 101-104)
                if !cache_out.accepted_hits.is_empty() {
                    cache_out.accepted_hits.sort_by(simple_hit_cmp_bins);
                    remove_duplicate_hits_pub(&mut cache_out.accepted_hits);
                    cache_out.map_type = if !cache_out.accepted_hits.is_empty() {
                        MappingType::MappedPair
                    } else {
                        MappingType::Unmapped
                    };

                    // C++: r1_len = shorter read, r2_len = longer read
                    let (r1_len, r2_len) = if r1.len() <= r2.len() {
                        (r1.len() as i32, r2.len() as i32)
                    } else {
                        (r2.len() as i32, r1.len() as i32)
                    };

                    for hit in &mut cache_out.accepted_hits {
                        hit.fragment_length = mate_ov.frag_length as i32;
                        // Compute mate_pos based on orientation
                        hit.mate_pos = if hit.is_fw {
                            hit.pos + hit.fragment_length - r2_len - 1
                        } else {
                            hit.pos + r1_len - hit.fragment_length - 1
                        };
                        // C++: clamp mate_pos > ref_len (not < 0)
                        let ref_len = index.ref_len(hit.tid as usize) as i32;
                        if hit.mate_pos > ref_len {
                            hit.mate_pos = hit.pos;
                        }
                        if mate_ov.ov_type == OverlapType::Dovetail {
                            hit.mate_pos = hit.pos;
                        }
                    }
                }
            } else {
                // No overlap: map both ends independently, then merge (bin-aware)
                poison_state.paired_for_mapping = true;
                map_pe_fragment_atac::<K>(
                    r1,
                    r2,
                    &mut hs,
                    &mut query,
                    &mut cache_left,
                    &mut cache_right,
                    &mut cache_out,
                    index,
                    &mut poison_state,
                    binning,
                );

                // Set fragment_length for orphan reads (matching C++ lines 204-213)
                if cache_out.map_type == MappingType::MappedFirstOrphan {
                    for hit in &mut cache_out.accepted_hits {
                        hit.fragment_length = r1.len() as i32;
                    }
                } else if cache_out.map_type == MappingType::MappedSecondOrphan {
                    for hit in &mut cache_out.accepted_hits {
                        hit.fragment_length = r2.len() as i32;
                    }
                }
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
            output.num_chunks.fetch_add(1, Ordering::Relaxed);
        }

        stats.num_reads.fetch_add(local_reads, Ordering::Relaxed);
        stats.num_mapped.fetch_add(local_mapped, Ordering::Relaxed);
        stats
            .num_poisoned
            .fetch_add(local_poisoned, Ordering::Relaxed);
    };

    run_mapping_pipeline_triple(fastx, thread_config, output, stats, worker_fn)
}
