//! Compare serialization sizes between cseq and sux-rs Elias-Fano

fn main() {
    // Create test sequences of different sizes
    // For gencode_pc_v44:
    // - .ssi contains string offsets (one EfSeqDict) with ~60K strings, universe ~19M bases
    // - .ctab contains contig offsets (one EfSeq) with ~60K contigs
    // - .ectab contains EC offsets (one EfSeq) with ~22K ECs

    let test_cases = [
        ("Small (1K elements)", 1_000, 1_000_000u64),
        ("Medium (10K elements)", 10_000, 5_000_000u64),
        ("Real .ssi (60K strings)", 60_800, 19_880_000u64),
        ("Large (100K elements)", 100_000, 50_000_000u64),
    ];

    for (name, n, max_val) in test_cases {
        println!("======== {} ========", name);
        run_test(n, max_val);
        println!();
    }
}

fn run_test(n: usize, max_val: u64) {
    use epserde::ser::Serialize;
    use mem_dbg::{MemSize, SizeFlags};
    use sux::dict::elias_fano::EliasFanoBuilder;

    // Generate a monotone sequence with semi-random gaps
    let mut values = Vec::new();
    let mut current = 0u64;
    let gap = max_val / n as u64;
    for _ in 0..n {
        values.push(current);
        current += gap;
    }
    values.push(max_val); // Final boundary
    let actual_n = values.len();

    let u = max_val as usize + 1;
    let mut builder = EliasFanoBuilder::new(actual_n, u);
    for &v in &values {
        builder.push(v as usize);
    }
    let ef = builder.build_with_seq_and_dict();

    // Serialize to measure size
    let mut buf = Vec::new();
    unsafe {
        ef.serialize(&mut buf).unwrap();
    }

    let mem_size = ef.mem_size(SizeFlags::default());

    println!("  Sequence: {} values, range [0, {}]", actual_n, max_val);
    println!(
        "  In-memory size: {} bytes ({:.2} KB)",
        mem_size,
        mem_size as f64 / 1_024.0
    );
    println!(
        "  Serialized size: {} bytes ({:.2} KB)",
        buf.len(),
        buf.len() as f64 / 1_024.0
    );
    println!(
        "  Theoretical bits/elem: {:.2}",
        2.0 + (u as f64 / actual_n as f64).log2()
    );
    println!(
        "  Actual bits/elem: {:.2}",
        (mem_size * 8) as f64 / actual_n as f64
    );
    println!(
        "  Serialized bits/elem: {:.2}",
        (buf.len() * 8) as f64 / actual_n as f64
    );

    // Parse epserde header
    let type_name_start = 29;
    let header_size = if buf.len() > type_name_start + 8 {
        let type_name_len = usize::from_le_bytes([
            buf[type_name_start],
            buf[type_name_start + 1],
            buf[type_name_start + 2],
            buf[type_name_start + 3],
            buf[type_name_start + 4],
            buf[type_name_start + 5],
            buf[type_name_start + 6],
            buf[type_name_start + 7],
        ]);
        type_name_start + 8 + type_name_len
    } else {
        29
    };

    println!(
        "  Epserde header: {} bytes ({:.1}%)",
        header_size,
        100.0 * header_size as f64 / buf.len() as f64
    );
    println!(
        "  Epserde alignment padding: {} bytes",
        buf.len() as isize - mem_size as isize - header_size as isize
    );

    // Breakdown of data structures
    let (n_ef, u_ef, l_ef, _low_bits, _high_bits, _first, _last) = ef.clone().into_parts();
    println!("  EliasFano breakdown:");
    println!(
        "    n={}, u={}, l={} (low bits per element)",
        n_ef, u_ef, l_ef
    );

    // Estimate component sizes
    let low_bits_size = (n_ef * l_ef).div_ceil(64) * 8;
    println!("    Low bits: ~{} bytes", low_bits_size);

    // High bits: n + u/2^l bits, plus Select1 and Select0 inventories
    let high_bits_raw = (n_ef + (u_ef >> l_ef)).div_ceil(64) * 8;
    println!("    High bits raw: ~{} bytes", high_bits_raw);

    // SelectAdaptConst<12, 3> inventory: every 2^12 bits stores absolute rank
    let select1_inventory = ((n_ef + (u_ef >> l_ef)) >> 12) * 8; // rough estimate
    let select0_inventory = ((n_ef + (u_ef >> l_ef)) >> 12) * 8;
    println!(
        "    Select1 inventory (~2^12 spacing): ~{} bytes",
        select1_inventory
    );
    println!(
        "    Select0 inventory (~2^12 spacing): ~{} bytes",
        select0_inventory
    );
}
