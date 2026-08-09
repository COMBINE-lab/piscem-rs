//! Stable help-contract checks for thread/decode semantics.

use clap::CommandFactory;

#[test]
fn mapping_help_describes_the_shared_budget_and_fixed_controls() {
    let mut command = piscem_rs::cli::Cli::command();
    for name in ["map-bulk", "map-scrna", "map-scatac"] {
        let subcommand = command
            .find_subcommand_mut(name)
            .unwrap_or_else(|| panic!("missing {name} subcommand"));
        let help = subcommand.render_long_help().to_string();
        let normalized = help.split_whitespace().collect::<Vec<_>>().join(" ");
        for required in [
            "Total execution-slot budget shared by mapping and gzip decoding",
            "adapts the aggregate mapping/decode split during the real run",
            "fixes N slots per decoder-capable input and disables adaptation",
            "Non-regular inputs remain serial",
        ] {
            assert!(
                normalized.contains(required),
                "{name} help lost `{required}`:\n{help}"
            );
        }
    }
}
