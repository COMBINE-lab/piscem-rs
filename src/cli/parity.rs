use anyhow::Result;
use clap::Args;

use crate::verify::parity::{run_parity, ParityConfig};

#[derive(Args, Debug)]
pub struct ParityArgs {
    /// Dataset directory or index prefix
    #[arg(long)]
    pub dataset: String,
    /// C++ piscem output directory for comparison
    #[arg(long)]
    pub cpp_bin_dir: Option<String>,
    /// Output report path
    #[arg(long)]
    pub output_report: Option<String>,
}

pub fn run(args: ParityArgs) -> Result<()> {
    let cfg = ParityConfig {
        dataset: args.dataset,
        cpp_bin_dir: args.cpp_bin_dir,
        output_report: args.output_report,
    };
    run_parity(cfg)
}
