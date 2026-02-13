use anyhow::Result;
use clap::Args;

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

pub fn run(_args: BuildPoisonArgs) -> Result<()> {
    anyhow::bail!(
        "build-poison command is not yet fully implemented.\n\
         The PoisonTable data structure is complete (Phase 2), but the full builder \
         requires the streaming query engine from Phase 3."
    )
}
