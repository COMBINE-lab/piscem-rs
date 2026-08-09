//! Emit a stable, order-independent binary representation of a RAD file.
//!
//! This validation utility deliberately writes to stdout so an oracle runner
//! can stream it directly into a cryptographic digest without retaining a
//! second full output file. It is not part of the mapping hot path.

#[cfg(feature = "parity-test")]
use std::io::{self, BufWriter, Write};
#[cfg(feature = "parity-test")]
use std::path::Path;

#[cfg(feature = "parity-test")]
use anyhow::{Context, Result, bail};
#[cfg(feature = "parity-test")]
use piscem_rs::verify::rad_compare::{
    read_atac_rad_records, read_bulk_rad_records, read_sc_multi_rad_records, read_sc_rad_records,
};

#[cfg(feature = "parity-test")]
fn write_u64(writer: &mut impl Write, value: u64) -> io::Result<()> {
    writer.write_all(&value.to_le_bytes())
}

#[cfg(feature = "parity-test")]
fn write_u32(writer: &mut impl Write, value: u32) -> io::Result<()> {
    writer.write_all(&value.to_le_bytes())
}

#[cfg(feature = "parity-test")]
fn write_u16(writer: &mut impl Write, value: u16) -> io::Result<()> {
    writer.write_all(&value.to_le_bytes())
}

#[cfg(feature = "parity-test")]
fn write_bulk(path: &Path, writer: &mut impl Write) -> Result<()> {
    let (_, mut records) = read_bulk_rad_records(path)?;
    records.sort_unstable();
    writer.write_all(b"piscem-canonical-bulk-v1\0")?;
    write_u64(writer, records.len() as u64)?;
    for record in records {
        writer.write_all(&[record.frag_type])?;
        write_u64(writer, record.alignments.len() as u64)?;
        for (reference, orientation, position, fragment_length) in record.alignments {
            write_u32(writer, reference)?;
            writer.write_all(&[orientation])?;
            write_u32(writer, position)?;
            write_u16(writer, fragment_length)?;
        }
    }
    Ok(())
}

#[cfg(feature = "parity-test")]
fn write_scrna(path: &Path, writer: &mut impl Write) -> Result<()> {
    let (_, mut records) = read_sc_rad_records(path)?;
    records.sort_unstable();
    writer.write_all(b"piscem-canonical-scrna-v1\0")?;
    write_u64(writer, records.len() as u64)?;
    for record in records {
        write_u64(writer, record.bc)?;
        write_u64(writer, record.umi)?;
        write_u64(writer, record.alignments.len() as u64)?;
        for (reference, forward) in record.alignments {
            write_u32(writer, reference)?;
            writer.write_all(&[u8::from(forward)])?;
        }
    }
    Ok(())
}

#[cfg(feature = "parity-test")]
fn write_atac(path: &Path, writer: &mut impl Write) -> Result<()> {
    let (_, mut records) = read_atac_rad_records(path)?;
    records.sort_unstable();
    writer.write_all(b"piscem-canonical-atac-v1\0")?;
    write_u64(writer, records.len() as u64)?;
    for record in records {
        write_u64(writer, record.bc)?;
        write_u64(writer, record.alignments.len() as u64)?;
        for (reference, alignment_type, position, fragment_length) in record.alignments {
            write_u32(writer, reference)?;
            writer.write_all(&[alignment_type])?;
            write_u32(writer, position)?;
            write_u16(writer, fragment_length)?;
        }
    }
    Ok(())
}

#[cfg(feature = "parity-test")]
fn write_flex(path: &Path, writer: &mut impl Write) -> Result<()> {
    let (_, mut records) = read_sc_multi_rad_records(path)?;
    records.sort_unstable();
    writer.write_all(b"piscem-canonical-flex-v1\0")?;
    write_u64(writer, records.len() as u64)?;
    for record in records {
        write_u64(writer, record.barcodes.len() as u64)?;
        for barcode in record.barcodes {
            write_u64(writer, barcode)?;
        }
        write_u64(writer, record.umi)?;
        write_u64(writer, record.alignments.len() as u64)?;
        for (reference, forward, position) in record.alignments {
            write_u32(writer, reference)?;
            writer.write_all(&[u8::from(forward), u8::from(position.is_some())])?;
            if let Some(position) = position {
                write_u32(writer, position)?;
            }
        }
    }
    Ok(())
}

#[cfg(feature = "parity-test")]
fn main() -> Result<()> {
    let mut args = std::env::args_os().skip(1);
    let kind = args
        .next()
        .context("usage: canonical_rad {bulk|scrna|flex|atac} PATH")?;
    let path = args
        .next()
        .context("usage: canonical_rad {bulk|scrna|flex|atac} PATH")?;
    if args.next().is_some() {
        bail!("usage: canonical_rad {{bulk|scrna|flex|atac}} PATH");
    }

    let stdout = io::stdout();
    let mut writer = BufWriter::new(stdout.lock());
    match kind.to_str() {
        Some("bulk") => write_bulk(Path::new(&path), &mut writer),
        Some("scrna") => write_scrna(Path::new(&path), &mut writer),
        Some("flex") => write_flex(Path::new(&path), &mut writer),
        Some("atac") => write_atac(Path::new(&path), &mut writer),
        _ => bail!("RAD kind must be bulk, scrna, flex, or atac"),
    }?;
    writer.flush()?;
    Ok(())
}

#[cfg(not(feature = "parity-test"))]
fn main() {
    eprintln!("canonical_rad requires --features parity-test");
}
