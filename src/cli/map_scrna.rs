use anyhow::Result;
use clap::Args;

#[derive(Args, Debug)]
pub struct MapScrnaArgs {
    #[arg(long)]
    pub index: String,
    #[arg(long)]
    pub read1: String,
    #[arg(long)]
    pub read2: Option<String>,
}

pub fn run(_args: MapScrnaArgs) -> Result<()> {
    anyhow::bail!("map-scrna command is not implemented yet (Phase 4A target)")
}
