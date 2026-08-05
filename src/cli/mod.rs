pub mod bench;
pub mod build;
pub mod map_bulk;
pub mod map_scatac;
pub mod map_scrna;
mod parity;
pub mod poison;
pub mod probe_bench;
mod stats;

use anyhow::Result;
use clap::{Parser, Subcommand, ValueEnum};

/// K-mer dictionary backend to use for mapping.
///
/// - `Auto` (default): at build time, emit Tiny artifacts iff the
///   canonical-k-mer count is below `AUTO_TINY_KMER_THRESHOLD`
///   (≈100 M, ~2 GB hashbrown). At map time, load Tiny when its
///   artifacts are present on disk, otherwise load sshash.
/// - `Sshash`: compact sshash dictionary only.
/// - `Tiny`: Tiny hashbrown-backed dictionary — faster per lookup,
///   higher RAM. At map time, loads prebuilt Tiny artifacts if present
///   or converts from sshash on the fly.
#[derive(Copy, Clone, Debug, ValueEnum, Default, PartialEq, Eq)]
#[value(rename_all = "lowercase")]
pub enum DictKind {
    #[default]
    Auto,
    Sshash,
    Tiny,
}

/// Canonical-k-mer threshold below which `--dict auto` selects the Tiny
/// backend at build time. At ~20 bytes/entry, 100 M k-mers ≈ 2 GB RAM.
pub const AUTO_TINY_KMER_THRESHOLD: u64 = 100_000_000;

impl DictKind {
    /// Resolve `Auto` for map-time dispatch based on which artifacts exist
    /// alongside the index prefix. Returns `Tiny` if the `.tdct` + `.tct`
    /// files are present, otherwise `Sshash`. Non-`Auto` variants pass
    /// through unchanged.
    pub fn resolve_for_map(self, index_prefix: &std::path::Path) -> Self {
        match self {
            DictKind::Auto => {
                if crate::index::reference_index::tiny_artifacts_exist(index_prefix) {
                    DictKind::Tiny
                } else {
                    DictKind::Sshash
                }
            }
            other => other,
        }
    }
}

#[derive(Parser, Debug)]
#[command(name = "piscem-rs")]
#[command(about = "Rust implementation of piscem")]
#[command(version = crate::VERSION)]
pub struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    Build(build::BuildArgs),
    MapScrna(map_scrna::MapScrnaArgs),
    MapBulk(map_bulk::MapBulkArgs),
    MapScatac(map_scatac::MapScatacArgs),
    BuildPoison(poison::BuildPoisonArgs),
    Stats(stats::StatsArgs),
    Parity(parity::ParityArgs),
    LookupBench(bench::LookupBenchArgs),
    /// Diagnostic: compare the calibration probe's producer estimate against a
    /// whole-file measurement of the same work, in this binary.
    #[command(hide = true)]
    ProbeBench(probe_bench::ProbeBenchArgs),
}

pub fn run() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Commands::Build(args) => build::run(args),
        Commands::MapScrna(args) => map_scrna::run(args),
        Commands::MapBulk(args) => map_bulk::run(args),
        Commands::MapScatac(args) => map_scatac::run(args),
        Commands::BuildPoison(args) => poison::run(args),
        Commands::Stats(args) => stats::run(args),
        Commands::Parity(args) => parity::run(args),
        Commands::LookupBench(args) => bench::run(args),
        Commands::ProbeBench(args) => probe_bench::run(args),
    }
}
