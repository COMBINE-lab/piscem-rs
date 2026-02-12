use anyhow::Result;
use clap::Args;

#[derive(Args, Debug)]
pub struct MapScatacArgs {
    #[arg(long)]
    pub index: String,
    #[arg(long)]
    pub read1: String,
    #[arg(long)]
    pub read2: Option<String>,
}

pub fn run(_args: MapScatacArgs) -> Result<()> {
    anyhow::bail!("map-scatac command is not implemented yet (Phase 4C target)")
}
