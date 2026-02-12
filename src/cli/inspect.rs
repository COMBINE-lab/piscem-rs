use anyhow::Result;
use clap::Args;

#[derive(Args, Debug)]
pub struct InspectArgs {
    #[arg(long)]
    pub index: String,
}

pub fn run(_args: InspectArgs) -> Result<()> {
    anyhow::bail!("inspect command is not implemented yet")
}
