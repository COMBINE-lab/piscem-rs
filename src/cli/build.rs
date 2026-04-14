use anyhow::Result;
use clap::Args;
use std::path::PathBuf;

use super::DictKind;
use crate::index::build::{BuildConfig, build_index};

#[derive(Args, Debug)]
pub struct BuildArgs {
    /// Input prefix (cuttlefish basename, e.g. path/to/index_cfish)
    #[arg(short = 'i', long)]
    pub input: PathBuf,
    /// Output prefix for index files
    #[arg(short = 'o', long)]
    pub output: PathBuf,
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
    /// Dictionary artifacts to emit. `auto` (default) emits Tiny artifacts when
    /// the canonical-k-mer count is below the auto threshold, else sshash-only.
    /// `sshash` writes only `.ssi`/`.ctab`; `tiny` additionally writes
    /// `.tdct`/`.tct` so `piscem map-* --dict tiny` can skip the runtime
    /// sshash→tiny conversion.
    #[arg(long, value_enum, default_value_t = DictKind::Auto)]
    pub dict: DictKind,
}

pub fn run(args: BuildArgs) -> Result<()> {
    let config = BuildConfig {
        input_prefix: args.input,
        output_prefix: args.output,
        k: args.klen,
        m: args.mlen,
        build_ec_table: args.build_ec_table,
        num_threads: args.threads,
        canonical: args.canonical,
        seed: args.seed,
        single_mphf: args.single_mphf,
        emit_tiny: match args.dict {
            DictKind::Tiny => Some(true),
            DictKind::Sshash => Some(false),
            DictKind::Auto => None,
        },
    };
    build_index(&config)
}
