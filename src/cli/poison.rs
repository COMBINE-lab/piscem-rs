use std::path::Path;
use std::time::Instant;

use anyhow::Result;
use clap::Args;
#[allow(unused_imports)]
use sshash_lib::{Kmer, KmerBits, dispatch_on_k};
use tracing::info;

use crate::index::build_poison::{build_poison_table, verify_poison_table};
use crate::index::reference_index::ReferenceIndex;

#[derive(Args, Debug)]
pub struct BuildPoisonArgs {
    /// Path prefix of the piscem index
    #[arg(short = 'i', long)]
    pub index: String,
    /// Decoy FASTA file(s) to scan for poison k-mers
    #[arg(short = 'd', long, num_args = 1..)]
    pub decoys: Vec<String>,
    /// Number of threads (0 = all cores)
    #[arg(short = 't', long, default_value = "0")]
    pub threads: usize,
}

pub fn run(args: BuildPoisonArgs) -> Result<()> {
    let start = Instant::now();

    let index_prefix = Path::new(&args.index);
    info!("Loading index from {}", index_prefix.display());
    let index = ReferenceIndex::load(index_prefix, false, false)?;
    info!(
        "Index loaded: k={}, {} refs, {} contigs",
        index.k(),
        index.num_refs(),
        index.num_contigs(),
    );

    let k = index.k();

    let table = dispatch_on_k!(k, K => {
        build_poison_table::<K>(&index, &args.decoys, args.threads)?
    });

    // Verify: no poison k-mers should be in the dictionary
    info!("Verifying poison table against dictionary...");
    let found_in_dict = dispatch_on_k!(k, K => {
        verify_poison_table::<K>(&table, &index)
    });
    if found_in_dict > 0 {
        tracing::warn!(
            "BUG: {} poison k-mers are also present in the reference dictionary!",
            found_in_dict
        );
    }

    // Save poison table
    let poison_path = format!("{}.poison", args.index);
    info!("Saving poison table to {}", poison_path);
    let mut poison_file = std::fs::File::create(&poison_path)?;
    table.save(&mut poison_file)?;

    // Save stats JSON
    let json_path = format!("{}.poison.json", args.index);
    info!("Saving poison stats to {}", json_path);
    let mut json_file = std::fs::File::create(&json_path)?;
    table.save_stats_json(&mut json_file)?;

    let elapsed = start.elapsed().as_secs_f64();
    info!(
        "Poison table built in {:.1}s: {} distinct k-mers, {} occurrences",
        elapsed,
        table.num_poison_kmers(),
        table.num_poison_occs(),
    );

    Ok(())
}
