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
    /// Deprecated, no effect: the index is always canonical (a k-mer and its
    /// reverse complement share one entry, and lookups report orientation).
    /// Accepted for one release; ignored with a warning
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
    /// Directory for sshash's external minimizer-sort scratch files.
    /// Unset keeps sshash's default (`sshash_tmp` in the current directory).
    #[arg(long)]
    pub tmp_dir: Option<PathBuf>,
    /// RAM ceiling, in GiB, for sshash's external minimizer sort. Unset keeps
    /// sshash's default (8 GiB); a smaller value spills to disk sooner.
    #[arg(long)]
    pub ram_limit_gib: Option<usize>,
}

pub fn run(args: BuildArgs) -> Result<()> {
    if args.canonical {
        tracing::warn!(
            "--canonical is deprecated and has no effect: the index is always canonical"
        );
    }
    let config = BuildConfig {
        input_prefix: args.input,
        output_prefix: args.output,
        k: args.klen,
        m: args.mlen,
        build_ec_table: args.build_ec_table,
        num_threads: args.threads,
        seed: args.seed,
        single_mphf: args.single_mphf,
        emit_tiny: match args.dict {
            DictKind::Tiny => Some(true),
            DictKind::Sshash => Some(false),
            DictKind::Auto => None,
        },
        tmp_dir: args.tmp_dir,
        ram_limit_gib: args.ram_limit_gib,
    };
    build_index(&config)
}
