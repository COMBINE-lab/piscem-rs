use anyhow::Result;
use clap::Args;

#[derive(Args, Debug)]
pub struct BuildPoisonArgs {
    #[arg(long)]
    pub index: String,
    #[arg(long)]
    pub decoys: Vec<String>,
}

pub fn run(_args: BuildPoisonArgs) -> Result<()> {
    anyhow::bail!("build-poison command is not implemented yet (Phase 2 target)")
}
