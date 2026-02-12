use anyhow::Result;
use clap::Args;

use crate::verify::parity::{run_phase0_parity, Phase0ParityConfig};

#[derive(Args, Debug)]
pub struct ParityArgs {
    #[arg(long)]
    pub dataset: String,
    #[arg(long)]
    pub cpp_bin_dir: Option<String>,
    #[arg(long)]
    pub output_report: Option<String>,
}

pub fn run(args: ParityArgs) -> Result<()> {
    let cfg = Phase0ParityConfig {
        dataset: args.dataset,
        cpp_bin_dir: args.cpp_bin_dir,
        output_report: args.output_report,
    };
    run_phase0_parity(cfg)
}
