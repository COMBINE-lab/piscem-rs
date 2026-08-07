//! The assertion whose absence cost three yanked releases.
//!
//! piscem-rs 0.7.0, 0.7.1 and 0.8.0 all shipped a parallel gzip decoder that was
//! *selected* and then never actually parallelized. Unit tests passed, clippy
//! passed, and mapping output was byte-identical throughout — because the defect
//! was never a wrong answer, it was a right answer that changed nothing.
//!
//! So these tests do not check a decision. They check that **the decision
//! produced its effect**: that workers were really spawned, and that the decoder
//! did not silently fall back to its sequential backend. Both are properties of
//! the running decoder, invisible to every other kind of test in this crate.
//!
//! Ignored by default and requires the `rapidgzip` feature plus real fixtures:
//!
//! ```text
//! cargo test --release --features rapidgzip --test decode_parallelism -- --ignored --nocapture
//! ```
//!
//! Fixture paths are overridable so this can run somewhere else:
//! `PISCEM_TEST_DENSE_GZ`, `PISCEM_TEST_MARKER_GZ`.

#![cfg(feature = "rapidgzip")]

use std::io::Read;
use std::path::PathBuf;

use piscem_rs::io::fastx::open_input;
use rapidgzip_core::DecoderPath;

/// Multi-member archive; takes `DecoderPath::DenseMembers`.
const DENSE_DEFAULT: &str = "/scratch1/rob/flex_bench/mm2/n1_R1_0.gz";
/// Single-member gzip; takes `DecoderPath::MarkerWindow`.
///
/// Both paths are covered deliberately. They are selected by entirely separate
/// admission code upstream, and a regression in one is invisible from the other:
/// `pigz`-written single-member input silently took `Sequential` for the whole
/// of rapidgzip-core 0.2.0 while multi-member input was unaffected.
const MARKER_DEFAULT: &str = "/scratch1/rob/rshash_testing/human_reads/SRR21186103_1.fastq.gz";

fn fixture(env: &str, default: &str) -> Option<PathBuf> {
    let p = PathBuf::from(std::env::var(env).unwrap_or_else(|_| default.to_string()));
    if p.is_file() {
        Some(p)
    } else {
        eprintln!("skipping: fixture {} not present (set {env})", p.display());
        None
    }
}

/// Decode a prefix and report what actually happened.
///
/// Returns `(peak_spawned_workers, path, bytes)`.
///
/// **Worker counts are sampled here rather than taken from
/// `DecodeBudgetReport`.** The supervisor polls every 250 ms, but a 256 MiB
/// prefix of a multi-member archive decodes in well under that, so its first
/// sample lands before any worker exists and its last never happens: it reports
/// a peak of 0 for a decode that really ran eight workers. That is a reporting
/// limitation, not a decode failure, and this test is about the decode.
fn decode_prefix(path: &PathBuf, ceiling: usize, mib: usize) -> (usize, DecoderPath, u64) {
    let opened = open_input(path, ceiling).expect("open fixture");
    let handle = opened
        .handle
        .expect("a regular gzip file must yield a decoder handle when the ceiling is nonzero");

    let mut reader = opened.reader;
    let mut buf = vec![0u8; 1 << 20];
    let mut total = 0u64;
    let mut peak = 0usize;
    let mut path_taken = DecoderPath::Starting;
    while total < (mib as u64) << 20 {
        match reader.read(&mut buf) {
            Ok(0) => break,
            Ok(n) => total += n as u64,
            Err(e) => panic!("read failed: {e}"),
        }
        // Sample while decoding is in flight. Afterwards the workers retire and
        // `spawned_workers` decays to zero, so a single check at the end would
        // race the retirement and flake.
        let s = handle.stats();
        peak = peak.max(s.spawned_workers);
        if s.path != DecoderPath::Starting {
            path_taken = s.path;
        }
    }

    drop(reader);
    (peak, path_taken, total)
}

/// A multi-member archive must really run workers, not merely be offered them.
#[test]
#[ignore = "needs the rapidgzip feature and a real fixture"]
fn dense_member_input_actually_spawns_workers() {
    let Some(f) = fixture("PISCEM_TEST_DENSE_GZ", DENSE_DEFAULT) else {
        return;
    };
    let (peak, path, bytes) = decode_prefix(&f, 8, 256);

    assert!(bytes > 0, "decoded nothing");
    assert_ne!(
        path,
        DecoderPath::Sequential,
        "decoder fell back to the sequential backend: the parallel path was \
         selected and did nothing, which is exactly the defect that yanked 0.7.0-0.8.0"
    );
    assert!(
        peak > 1,
        "peak worker threads {peak}: the decoder was granted a ceiling of 8 and \
         never used more than one, so the decision changed nothing"
    );
    eprintln!(
        "dense: path {path:?}, peak {peak} workers, {} MiB",
        bytes >> 20
    );
}

/// The single-member marker path, admitted by entirely separate upstream code.
#[test]
#[ignore = "needs the rapidgzip feature and a real fixture"]
fn marker_window_input_actually_spawns_workers() {
    let Some(f) = fixture("PISCEM_TEST_MARKER_GZ", MARKER_DEFAULT) else {
        return;
    };
    let (peak, path, bytes) = decode_prefix(&f, 8, 256);

    assert!(bytes > 0, "decoded nothing");
    assert_ne!(
        path,
        DecoderPath::Sequential,
        "single-member gzip fell back to sequential. If this fixture was written \
         by `pigz`, check the rapidgzip-core version: every `pigz` stream took \
         this path before 0.2.1 (COMBINE-lab/rapidgzip-rust#29)"
    );
    assert!(
        peak > 1,
        "peak worker threads {peak}: granted a ceiling of 8 and never used more than one"
    );
    eprintln!(
        "marker: path {path:?}, peak {peak} workers, {} MiB",
        bytes >> 20
    );
}

/// A decoder must not spend more live workers than its ceiling.
///
/// The per-run budget is now the shared `DecoderPool`'s aggregate limit rather
/// than a supervisor summing per-file grants, but the property this guards is
/// the same: what a decoder is allowed must bound what it actually spends.
#[test]
#[ignore = "needs the rapidgzip feature and a real fixture"]
fn the_supervisor_holds_the_budget() {
    let Some(f) = fixture("PISCEM_TEST_DENSE_GZ", DENSE_DEFAULT) else {
        return;
    };
    for budget in [2usize, 4, 8] {
        let (peak, path, _) = decode_prefix(&f, budget, 128);
        assert_ne!(path, DecoderPath::Sequential, "budget {budget}");
        assert!(
            peak <= budget,
            "budget {budget} but {peak} live workers -- the ceiling must bound \
             what a decoder actually spends, since the supervisor only ever \
             raises ceilings and cannot claw one back"
        );
        assert!(peak > 1, "budget {budget} bought only {peak} worker(s)");
    }
}

/// `--decoder serial` must produce no decoder at all, not a parallel decoder
/// that happens to run one worker.
#[test]
#[ignore = "needs the rapidgzip feature and a real fixture"]
fn a_serial_choice_opens_no_decoder() {
    let Some(f) = fixture("PISCEM_TEST_DENSE_GZ", DENSE_DEFAULT) else {
        return;
    };
    // Ceiling 0 is how the serial verdict is plumbed through `open_input`.
    let opened = open_input(&f, 0).expect("open fixture");
    assert!(
        opened.handle.is_none(),
        "the serial path must not construct a parallel decoder"
    );
}
