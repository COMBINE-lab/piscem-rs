use anyhow::Result;

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

fn main() -> Result<()> {
    tracing_subscriber::fmt()
        // Default to `warn`, and defer to RUST_LOG only when it is actually
        // set. Warnings have to reach the user unprompted -- a message saying
        // "your --decoder request was overridden" is useless if it is visible
        // only to someone who already suspected a problem.
        //
        // Keyed on the variable existing rather than on parsing succeeding:
        // falling back on a parse failure would silently discard a RUST_LOG the
        // user did set, which is the opposite of deferring to it.
        .with_env_filter(match std::env::var_os("RUST_LOG") {
            Some(_) => tracing_subscriber::EnvFilter::from_default_env(),
            None => tracing_subscriber::EnvFilter::new("warn"),
        })
        .with_target(false)
        .compact()
        .init();

    piscem_rs::cli::run()
}
