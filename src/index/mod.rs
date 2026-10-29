pub mod build;
pub mod build_poison;
pub mod contig_table;
pub mod eq_classes;
pub mod formats;
pub mod mem_budget;
pub mod paths;
pub mod poison_table;
pub mod reference_index;
pub mod refinfo;

/// Checks the 8-byte magic of a versioned table (`<prefix>NN\0`, e.g.
/// `PCTAB04\0`). A magic of the same table at another version gets an error
/// saying the index must be rebuilt, rather than a bare mismatch.
pub(crate) fn check_table_magic(
    what: &str,
    got: &[u8; 8],
    expected: &[u8; 8],
) -> anyhow::Result<()> {
    if got == expected {
        return Ok(());
    }
    let name = |m: &[u8; 8]| {
        String::from_utf8_lossy(m)
            .trim_end_matches('\0')
            .to_string()
    };
    let prefix = expected
        .iter()
        .position(|b| b.is_ascii_digit())
        .unwrap_or(expected.len());
    if got[..prefix] == expected[..prefix] {
        anyhow::bail!(
            "unsupported {what} format {} (this build reads {}): the index was written by a \
             different piscem release and must be rebuilt",
            name(got),
            name(expected)
        );
    }
    anyhow::bail!(
        "invalid {what} magic: expected {:?}, got {:?}",
        expected,
        got
    )
}
