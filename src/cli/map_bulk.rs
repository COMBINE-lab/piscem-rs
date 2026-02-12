use anyhow::Result;
use clap::Args;

#[derive(Args, Debug)]
pub struct MapBulkArgs {
    #[arg(long)]
    pub index: String,
    #[arg(long)]
    pub read1: String,
    #[arg(long)]
    pub read2: Option<String>,
}

pub fn run(_args: MapBulkArgs) -> Result<()> {
    anyhow::bail!("map-bulk command is not implemented yet (Phase 4B target)")
}
