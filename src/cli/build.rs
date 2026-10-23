use anyhow::Result;
use clap::{Args, ValueEnum};
use sshash_lib::{MphfBuild, Offline};
use std::path::PathBuf;

use super::DictKind;
use crate::index::build::{BuildConfig, build_index};

#[derive(Args, Debug)]
pub struct BuildArgs {
    /// Input prefix (cuttlefish basename, e.g. path/to/index_cfish)
    #[arg(short = 'i', long)]
    pub input: PathBuf,
    /// Output prefix for index files, e.g. `idx/gencode.v49` writes
    /// `idx/gencode.v49.ssi`, `idx/gencode.v49.ctab`, ... (not a directory)
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
    /// How the minimizer MPHF is built. `streaming` (default) builds one
    /// function from the sorted minimizers without collecting them, offline
    /// (in `--tmp-dir`) for large references; `in-memory` builds one function
    /// from a copy of the minimizers (no temporary disk, more memory);
    /// `partitioned` builds ~3M-key functions in parallel (no temporary disk,
    /// ~25% slower minimizer lookups). The index format is the same.
    #[arg(long, value_enum, default_value_t = MphfArg::Streaming)]
    pub mphf: MphfArg,
    /// Where `--mphf streaming` keeps its construction records (16 bytes per
    /// minimizer): `auto` (default) in `--tmp-dir` from 2^24 minimizers on,
    /// `always` there (least memory), `never` in memory (no temporary disk).
    #[arg(long, value_enum, default_value_t = OfflineArg::Auto)]
    pub mphf_offline: OfflineArg,
    /// Deprecated: same as `--mphf in-memory`.
    #[arg(long, hide = true, conflicts_with = "mphf")]
    pub single_mphf: bool,
    /// Dictionary artifacts to emit. `auto` (default) emits Tiny artifacts when
    /// the canonical-k-mer count is below the auto threshold, else sshash-only.
    /// `sshash` writes only `.ssi`/`.ctab`; `tiny` additionally writes
    /// `.tdct`/`.tct` so `piscem map-* --dict tiny` can skip the runtime
    /// sshash→tiny conversion.
    #[arg(long, value_enum, default_value_t = DictKind::Auto)]
    pub dict: DictKind,
    /// Directory for sshash's temporary files: the external minimizer sort,
    /// and the offline MPHF construction records (`--mphf-offline`).
    /// Unset keeps sshash's default (`sshash_tmp` in the current directory).
    #[arg(long)]
    pub tmp_dir: Option<PathBuf>,
    /// RAM ceiling, in GiB, for sshash's external minimizer sort; a smaller
    /// value spills to disk sooner, 0 means unlimited (always sort in memory).
    /// Unset: half the memory this process may use (the least of physical RAM,
    /// RLIMIT_AS and the cgroup limit), at least 8 GiB or all of it if less is
    /// usable; 8 GiB if it cannot be detected.
    #[arg(long)]
    pub ram_limit_gib: Option<usize>,
}

/// `--mphf` values.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum MphfArg {
    Streaming,
    InMemory,
    Partitioned,
}

/// `--mphf-offline` values.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum OfflineArg {
    Auto,
    Always,
    Never,
}

impl BuildArgs {
    /// The sshash MPHF strategy these arguments select.
    pub fn mphf_build(&self) -> MphfBuild {
        if self.single_mphf {
            return MphfBuild::InMemory;
        }
        match self.mphf {
            MphfArg::Streaming => MphfBuild::Streaming(match self.mphf_offline {
                OfflineArg::Auto => Offline::default(),
                OfflineArg::Always => Offline::Always,
                OfflineArg::Never => Offline::Never,
            }),
            MphfArg::InMemory => MphfBuild::InMemory,
            MphfArg::Partitioned => MphfBuild::Partitioned,
        }
    }
}

pub fn run(args: BuildArgs) -> Result<()> {
    if args.canonical {
        tracing::warn!(
            "--canonical is deprecated and has no effect: the index is always canonical"
        );
    }
    if args.single_mphf {
        tracing::warn!("--single-mphf is deprecated: use `--mphf in-memory`");
    }
    if args.mphf != MphfArg::Streaming && args.mphf_offline != OfflineArg::Auto {
        tracing::warn!("--mphf-offline only applies to `--mphf streaming`; ignored");
    }
    let mphf_build = args.mphf_build();
    let config = BuildConfig {
        input_prefix: args.input,
        output_prefix: args.output,
        k: args.klen,
        m: args.mlen,
        build_ec_table: args.build_ec_table,
        num_threads: args.threads,
        seed: args.seed,
        mphf_build,
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
