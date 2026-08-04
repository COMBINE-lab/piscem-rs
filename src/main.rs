use anyhow::Result;

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

fn main() -> Result<()> {
    tracing_subscriber::fmt()
        // Warnings must reach the user without them having set RUST_LOG:
        // messages like "your --decoder request was overridden" are useless if
        // they are only visible to someone who already suspected a problem.
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("warn")),
        )
        .with_target(false)
        .compact()
        .init();

    piscem_rs::cli::run()
}
