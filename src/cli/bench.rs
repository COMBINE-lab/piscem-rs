use anyhow::Result;
use clap::Args;
use sshash_lib::{Dictionary, Kmer, KmerBits, dispatch_on_k};
use std::path::PathBuf;
use std::time::Instant;
use tracing::info;

use crate::index::contig_table::ContigTable;
use crate::index::reference_index::ReferenceIndex;

#[derive(Args, Debug)]
pub struct LookupBenchArgs {
    /// Index prefix
    #[arg(short, long)]
    index: String,

    /// Query FASTA/FASTQ file
    #[arg(short, long)]
    query: String,

    /// Use point lookup instead of streaming
    #[arg(long)]
    point: bool,

    /// Skip contig table locate (measure pure sshash time)
    #[arg(long)]
    no_locate: bool,
}

fn load_sequences(path: &str) -> Result<Vec<Vec<u8>>> {
    use std::io::BufRead;

    let file = std::fs::File::open(path)?;
    let reader: Box<dyn std::io::Read> = if path.ends_with(".gz") {
        Box::new(flate2::read::GzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let buf = std::io::BufReader::new(reader);
    let mut sequences = Vec::new();
    let mut lines = buf.lines();

    let is_fastq = path.contains(".fq") || path.contains(".fastq");

    while let Some(Ok(header)) = lines.next() {
        if header.is_empty() {
            break;
        }
        if let Some(Ok(seq)) = lines.next() {
            sequences.push(seq.into_bytes());
        }
        if is_fastq {
            lines.next(); // +
            lines.next(); // quality
        }
    }

    Ok(sequences)
}

fn run_streaming_bench<const K: usize>(
    ri: &ReferenceIndex<Dictionary, ContigTable>,
    sequences: &[Vec<u8>],
    locate: bool,
) where
    Kmer<K>: KmerBits,
{
    let dict = ri.dict();
    let k = dict.k();
    let encoding = ri.encoding();
    let mut engine = dict.create_streaming_query::<K>();
    let mut num_kmers: u64 = 0;
    let mut found: u64 = 0;
    let mut decoded: u64 = 0;

    let t_start = Instant::now();

    let ct = ri.contig_table();
    let mut cached_entries: Vec<u64> = Vec::with_capacity(64);

    for seq in sequences {
        if seq.len() < k {
            continue;
        }
        engine.reset();
        let n = seq.len() - k + 1;
        let mut prev_string_id = u64::MAX;

        for i in 0..n {
            let kmer_bytes = unsafe { seq.get_unchecked(i..i + k) };
            let result = engine.lookup(kmer_bytes);
            num_kmers += 1;
            if result.is_found() {
                found += 1;
                if locate {
                    let string_id = result.string_id;
                    if string_id != prev_string_id {
                        prev_string_id = string_id;
                        cached_entries.clear();
                        for entry in ct.contig_entries(string_id).iter() {
                            cached_entries.push(entry);
                        }
                    }
                    for &entry in &cached_entries {
                        decoded += encoding.pos(entry) as u64;
                    }
                }
            }
        }
    }

    let elapsed = t_start.elapsed();
    let ns_per_kmer = elapsed.as_nanos() as f64 / num_kmers as f64;
    let stats = engine.stats();

    let label = if locate {
        "streaming lookup+locate"
    } else {
        "streaming lookup"
    };
    println!("==== {} report (piscem-rs):", label);
    println!("num_kmers = {}", num_kmers);
    println!(
        "found_kmers = {} ({:.2}%)",
        found,
        if num_kmers > 0 {
            found as f64 / num_kmers as f64 * 100.0
        } else {
            0.0
        }
    );
    if locate {
        println!("decoded_positions = {}", decoded);
    }
    println!("searches = {}", stats.num_searches);
    println!("extensions = {}", stats.num_extensions);
    println!(
        "extension_ratio = {:.2}",
        if stats.num_searches > 0 {
            stats.num_extensions as f64 / stats.num_searches as f64
        } else {
            0.0
        }
    );
    println!("time_per_kmer = {:.3} ns", ns_per_kmer);
    println!("total_time = {:.6} s", elapsed.as_secs_f64());
}

fn run_point_bench<const K: usize>(
    ri: &ReferenceIndex<Dictionary, ContigTable>,
    sequences: &[Vec<u8>],
) where
    Kmer<K>: KmerBits,
{
    let dict = ri.dict();
    let k = dict.k();
    let encoding = ri.encoding();
    let mut engine = dict.create_streaming_query::<K>();
    let mut num_kmers: u64 = 0;
    let mut found: u64 = 0;
    let mut decoded: u64 = 0;

    let t_start = Instant::now();

    for seq in sequences {
        if seq.len() < k {
            continue;
        }
        let n = seq.len() - k + 1;
        for i in 0..n {
            engine.reset();
            let kmer_bytes = &seq[i..i + k];
            let result = engine.lookup(kmer_bytes);
            num_kmers += 1;
            if result.is_found() {
                found += 1;
                if let Some(proj_hits) = ri.resolve_lookup(&result) {
                    for entry in proj_hits.ref_range().iter() {
                        let ref_pos = proj_hits.decode_hit(entry, &encoding);
                        decoded += ref_pos.pos as u64;
                    }
                }
            }
        }
    }

    let elapsed = t_start.elapsed();
    let ns_per_kmer = elapsed.as_nanos() as f64 / num_kmers as f64;

    println!("==== point lookup+locate report (piscem-rs):");
    println!("num_kmers = {}", num_kmers);
    println!(
        "found_kmers = {} ({:.2}%)",
        found,
        if num_kmers > 0 {
            found as f64 / num_kmers as f64 * 100.0
        } else {
            0.0
        }
    );
    println!("decoded_positions = {}", decoded);
    println!("time_per_kmer = {:.3} ns", ns_per_kmer);
    println!("total_time = {:.6} s", elapsed.as_secs_f64());
}

pub fn run(args: LookupBenchArgs) -> Result<()> {
    let prefix = PathBuf::from(&args.index);
    info!("Loading index from {}", prefix.display());
    let ri: ReferenceIndex<Dictionary, ContigTable> = ReferenceIndex::load(&prefix, false, false)?;
    let k = ri.k();
    info!("Index loaded: k={}, {} contigs", k, ri.num_contigs());

    info!("Loading query sequences from {}...", args.query);
    let sequences = load_sequences(&args.query)?;
    let total_kmers: u64 = sequences
        .iter()
        .filter(|s| s.len() >= k)
        .map(|s| (s.len() - k + 1) as u64)
        .sum();
    info!(
        "Loaded {} sequences, {} k-mer positions",
        sequences.len(),
        total_kmers
    );

    let locate = !args.no_locate;
    if args.point {
        dispatch_on_k!(k, K => {
            run_point_bench::<K>(&ri, &sequences);
        });
    } else {
        dispatch_on_k!(k, K => {
            run_streaming_bench::<K>(&ri, &sequences, locate);
        });
    }

    Ok(())
}
