pub mod build;
pub mod map_bulk;
pub mod map_scatac;
pub mod map_scrna;
mod parity;
pub mod poison;
mod stats;

use anyhow::Result;
use clap::{Parser, Subcommand};

#[derive(Parser, Debug)]
#[command(name = "piscem-rs")]
#[command(about = "Rust implementation of piscem")]
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
    }
}
