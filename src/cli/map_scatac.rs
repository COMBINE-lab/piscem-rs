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
use crate::io::fastx::{open_with_decompression, Collection, CollectionType};
use crate::io::map_info::write_map_info;
use crate::io::rad::write_rad_header_atac;
use crate::io::threads::{MappingStats, OutputInfo};
use crate::mapping::binning::BinPos;
use crate::mapping::processors::{MappingOpts, ScatacProcessor};
use crate::mapping::unitig_end_cache::UnitigEndCache;

use super::map_bulk::make_progress_bar;

#[derive(Args, Debug)]
#[command(group(
    clap::ArgGroup::new("bio_reads")
        .required(true)
        .args(["reads", "read1"])
))]
pub struct MapScatacArgs {
    /// Index prefix path
    #[arg(short = 'i', long)]
    pub index: String,
    /// Single-end biological FASTQ files (comma-separated); mutually exclusive with -1/-2
    #[arg(short = 'r', long = "reads", value_delimiter = ',',
          conflicts_with_all = ["read1", "read2"])]
    pub reads: Vec<String>,
    /// Read 1 FASTQ files (genomic left, comma-separated); requires -2
    #[arg(short = '1', long, value_delimiter = ',',
          requires = "read2", conflicts_with = "reads")]
    pub read1: Vec<String>,
    /// Read 2 FASTQ files (genomic right, comma-separated); requires -1
    #[arg(short = '2', long, value_delimiter = ',',
          requires = "read1", conflicts_with = "reads")]
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
    /// Enable ambiguous hit checking (EC table loading); disabled by default for ATAC
    #[arg(long)]
    pub check_ambig_hits: bool,
    /// Maximum equivalence class cardinality for ambiguous hit resolution
    #[arg(long, default_value = "256")]
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
    /// UnitigEndCache capacity (number of entries)
    #[arg(long, default_value = "5000000")]
    pub end_cache_capacity: usize,
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
    /// Suppress progress output
    #[arg(short = 'q', long)]
    pub quiet: bool,
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

    let is_paired = !args.read1.is_empty();

    // Load index
    let index_prefix = Path::new(&args.index);
    info!("Loading index from {}", index_prefix.display());
    // C++ ATAC: check_ambig_hits defaults to false, so EC table is NOT loaded
    let index = ReferenceIndex::load(index_prefix, args.check_ambig_hits, !args.no_poison)?;
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

    // Setup progress bar
    let progress = make_progress_bar(args.quiet);

    let k = index.k();
    let tn5_shift = !args.no_tn5_shift;
    let min_overlap = args.min_overlap;
    let end_cache = UnitigEndCache::new(args.end_cache_capacity);
    let num_threads = args.threads.max(1);
    let opts = MappingOpts {
        max_hit_occ: args.max_hit_occ,
        max_hit_occ_recover: args.max_hit_occ_recover,
        max_read_occ: args.max_read_occ,
        max_ec_card: args.max_ec_card,
    };

    // Dispatch on K and run the pipeline via paraseq multi-reader
    dispatch_on_k!(k, K => {
        let (bio_paths, r2_paths) = if is_paired {
            (args.read1.as_slice(), args.read2.as_slice())
        } else {
            (args.reads.as_slice(), [].as_slice())
        };
        run_atac_pipeline::<K>(
            bio_paths, &args.barcode, r2_paths,
            &output_info, &stats,
            &index, &binning, bc_len, tn5_shift, min_overlap, Some(&end_cache),
            is_paired, opts,
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
fn run_atac_pipeline<const K: usize>(
    bio_paths: &[String],
    barcode_paths: &[String],
    read2_paths: &[String],
    output: &OutputInfo,
    stats: &MappingStats,
    index: &ReferenceIndex,
    binning: &BinPos,
    bc_len: u16,
    tn5_shift: bool,
    min_overlap: i32,
    end_cache: Option<&UnitigEndCache>,
    is_paired: bool,
    opts: MappingOpts,
    num_threads: usize,
    progress: &indicatif::ProgressBar,
) -> Result<()>
where
    Kmer<K>: KmerBits,
{
    let mut processor = ScatacProcessor::<K>::new(
        index, end_cache, output, stats, binning, bc_len, tn5_shift, min_overlap, is_paired, opts,
        progress,
    );

    let mut readers = Vec::with_capacity(bio_paths.len() * if is_paired { 3 } else { 2 });
    for i in 0..bio_paths.len() {
        readers.push(
            paraseq::fastx::Reader::new(open_with_decompression(&bio_paths[i])?)
                .map_err(|e| anyhow::anyhow!("failed to open {}: {}", bio_paths[i], e))?,
        );
        readers.push(
            paraseq::fastx::Reader::new(open_with_decompression(&barcode_paths[i])?)
                .map_err(|e| anyhow::anyhow!("failed to open {}: {}", barcode_paths[i], e))?,
        );
        if is_paired {
            readers.push(
                paraseq::fastx::Reader::new(open_with_decompression(&read2_paths[i])?)
                    .map_err(|e| anyhow::anyhow!("failed to open {}: {}", read2_paths[i], e))?,
            );
        }
    }
    let arity = if is_paired { 3 } else { 2 };
    let collection = Collection::new(readers, CollectionType::Multi { arity })
        .map_err(|e| anyhow::anyhow!("failed to create collection: {}", e))?;
    collection
        .process_parallel_multi(&mut processor, num_threads, None)
        .map_err(|e| anyhow::anyhow!("mapping failed: {}", e))?;

    Ok(())
}
