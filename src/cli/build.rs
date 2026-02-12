use anyhow::Result;
use clap::Args;

#[derive(Args, Debug)]
pub struct BuildArgs {
    #[arg(long)]
    pub input: String,
    #[arg(long)]
    pub output: String,
}

pub fn run(_args: BuildArgs) -> Result<()> {
    anyhow::bail!("build command is not implemented yet (Phase 1 target)")
}
