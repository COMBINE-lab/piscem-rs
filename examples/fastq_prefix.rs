//! Materialize a record-aligned gzip FASTQ prefix for validation matrices.
//!
//! This intentionally small utility is for reproducible benchmark fixtures,
//! not the mapping path. It accepts ordinary four-line FASTQ records from any
//! compression format supported by `niffler`. An optional fourth argument
//! starts a new gzip member after that many records, which makes the dense-
//! member rapidgzip path available to validation fixtures.

use std::ffi::OsString;
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;

use flate2::Compression;
use flate2::write::GzEncoder;

fn usage(program: &OsString) -> ! {
    eprintln!(
        "usage: {} INPUT OUTPUT RECORDS [RECORDS_PER_GZIP_MEMBER]",
        PathBuf::from(program).display()
    );
    std::process::exit(2);
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args_os().collect();
    if !(4..=5).contains(&args.len()) {
        usage(&args[0]);
    }
    let input = PathBuf::from(&args[1]);
    let output = PathBuf::from(&args[2]);
    let records: usize = args[3].to_str().ok_or("RECORDS is not UTF-8")?.parse()?;
    let records_per_member = if let Some(value) = args.get(4) {
        value
            .to_str()
            .ok_or("RECORDS_PER_GZIP_MEMBER is not UTF-8")?
            .parse::<usize>()?
    } else {
        records.max(1)
    };
    if records_per_member == 0 {
        return Err("RECORDS_PER_GZIP_MEMBER must be positive".into());
    }

    let (source, _) = niffler::send::from_path(&input)?;
    let mut source = BufReader::new(source);
    let mut output_file = std::fs::File::create(&output)?;
    let mut lines = [Vec::new(), Vec::new(), Vec::new(), Vec::new()];

    let mut written = 0;
    while written < records {
        let member_end = records.min(written + records_per_member);
        let mut sink = GzEncoder::new(&mut output_file, Compression::fast());
        while written < member_end {
            for line in &mut lines {
                line.clear();
                if source.read_until(b'\n', line)? == 0 {
                    return Err(format!(
                        "{} ended after {written} complete records; requested {records}",
                        input.display()
                    )
                    .into());
                }
            }
            if !lines[0].starts_with(b"@") || !lines[2].starts_with(b"+") {
                return Err(
                    format!("record {} is not a four-line FASTQ record", written + 1).into(),
                );
            }
            let sequence_len = lines[1].strip_suffix(b"\n").unwrap_or(&lines[1]).len();
            let quality_len = lines[3].strip_suffix(b"\n").unwrap_or(&lines[3]).len();
            if sequence_len != quality_len {
                return Err(format!(
                    "record {} has sequence length {sequence_len} but quality length {quality_len}",
                    written + 1
                )
                .into());
            }
            for line in &lines {
                sink.write_all(line)?;
            }
            written += 1;
        }
        sink.finish()?;
    }
    output_file.sync_all()?;
    Ok(())
}
