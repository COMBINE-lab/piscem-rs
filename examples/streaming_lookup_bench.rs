use std::io::BufRead;
use std::path::Path;
use std::time::Instant;

use piscem_rs::index::reference_index::ReferenceIndex;
use sshash_lib::{Kmer, KmerBits, KmerDictionary, KmerStreamingQuery, dispatch_on_k};

fn run_bench<const K: usize, D: KmerDictionary>(dict: &D, sequences: &[Vec<u8>], k: usize)
where
    Kmer<K>: KmerBits,
{
    let mut engine = dict.create_streaming_query::<K>();

    let mut found: u64 = 0;
    let mut num_kmers: u64 = 0;

    eprintln!("starting benchmark...");
    let t_start = Instant::now();

    for seq in sequences {
        if seq.len() < k {
            continue;
        }
        engine.reset();
        let n = seq.len() - k + 1;
        for i in 0..n {
            let kmer_bytes = &seq[i..i + k];
            let result = engine.lookup(kmer_bytes);
            num_kmers += 1;
            if result.is_found() {
                found += 1;
            }
        }
    }

    let elapsed = t_start.elapsed();
    let ns_per_kmer = elapsed.as_nanos() as f64 / num_kmers as f64;
    let extensions = engine.num_extensions();
    let searches = engine.num_searches();

    println!("==== streaming lookup report:");
    println!("num_kmers = {}", num_kmers);
    println!(
        "found_kmers = {} ({:.4}%)",
        found,
        if num_kmers > 0 {
            found as f64 / num_kmers as f64 * 100.0
        } else {
            0.0
        }
    );
    println!("searches = {}", searches);
    println!("extensions = {}", extensions);
    println!(
        "extension_ratio = {:.4}",
        if searches > 0 {
            extensions as f64 / searches as f64
        } else {
            0.0
        }
    );
    println!("time_per_kmer = {:.4} ns", ns_per_kmer);
    println!("total_time = {:.3} s", elapsed.as_secs_f64());
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: streaming_lookup_bench <index_prefix> <query.fastq.gz>");
        std::process::exit(1);
    }
    let index_prefix = &args[1];
    let query_file = &args[2];

    eprintln!("loading index from {}", index_prefix);
    let index = ReferenceIndex::load(Path::new(index_prefix), false, false)?;
    let k = index.k();
    eprintln!("index loaded (k={})", k);

    eprintln!("loading queries from {}", query_file);
    let (reader, _format) = niffler::send::from_path(Path::new(query_file))?;
    let buf = std::io::BufReader::new(reader);

    let mut sequences: Vec<Vec<u8>> = Vec::with_capacity(300_000);
    let mut line_num = 0u64;
    for line in buf.lines() {
        let line = line?;
        line_num += 1;
        // FASTQ: sequence is every 4th line starting at line 2
        if line_num % 4 == 2 {
            sequences.push(line.into_bytes());
        }
    }

    let total_kmers: u64 = sequences
        .iter()
        .filter(|s| s.len() >= k)
        .map(|s| (s.len() - k + 1) as u64)
        .sum();
    eprintln!(
        "loaded {} sequences, {} k-mer positions",
        sequences.len(),
        total_kmers
    );

    let dict = index.dict();
    dispatch_on_k!(k, K => {
        run_bench::<K, _>(dict, &sequences, k);
    });

    Ok(())
}
