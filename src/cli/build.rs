use anyhow::Result;
use clap::Args;
use std::path::PathBuf;

use crate::index::build::{BuildConfig, build_index};

#[derive(Args, Debug)]
pub struct BuildArgs {
    /// Input prefix (cuttlefish basename, e.g. path/to/index_cfish)
    #[arg(short = 'i', long)]
    pub input: String,
    /// Output prefix for index files
    #[arg(short = 'o', long)]
    pub output: String,
    /// K-mer length
    #[arg(short = 'k', long)]
    pub klen: usize,
    /// Minimizer length
    #[arg(short = 'm', long)]
    pub mlen: usize,
    /// Number of threads (0 = all cores)
    #[arg(short = 't', long, default_value = "0")]
    pub threads: usize,
    /// Build equivalence class table
    #[arg(long)]
    pub build_ec_table: bool,
    /// Use canonical k-mer mode
    #[arg(long)]
    pub canonical: bool,
    /// Hash seed for dictionary construction
    #[arg(short = 's', long, default_value = "1")]
    pub seed: u64,
    /// Use a single monolithic MPHF instead of partitioned (disables parallel MPHF build)
    #[arg(long)]
    pub single_mphf: bool,
}

pub fn run(args: BuildArgs) -> Result<()> {
    let config = BuildConfig {
        input_prefix: PathBuf::from(&args.input),
        output_prefix: PathBuf::from(&args.output),
        k: args.klen,
        m: args.mlen,
        build_ec_table: args.build_ec_table,
        num_threads: args.threads,
        canonical: args.canonical,
        seed: args.seed,
        single_mphf: args.single_mphf,
    };
    build_index(&config)
}
