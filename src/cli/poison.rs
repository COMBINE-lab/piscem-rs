use std::path::PathBuf;
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
    pub index: PathBuf,
    /// Decoy FASTA file(s) to scan for poison k-mers
    #[arg(short = 'd', long, num_args = 1..)]
    pub decoys: Vec<PathBuf>,
    /// Number of threads (0 = all cores)
    #[arg(short = 't', long, default_value = "0")]
    pub threads: usize,
}

pub fn run(args: BuildPoisonArgs) -> Result<()> {
    let start = Instant::now();

    info!("Loading index from {}", args.index.display());
    let index = ReferenceIndex::load(&args.index, false, false)?;
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
    let mut poison_path = args.index.clone();
    poison_path.add_extension("poison");
    info!("Saving poison table to {}", poison_path.display());
    let mut poison_file = std::fs::File::create(&poison_path)?;
    table.save(&mut poison_file)?;

    // Save stats JSON
    let mut json_path = poison_path.clone();
    json_path.add_extension("json");
    info!("Saving poison stats to {}", json_path.display());
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
