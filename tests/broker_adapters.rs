//! The `Producer` adapter against real decoder states.
//!
//! The control law is tested against fakes in `thread-broker`. What cannot be
//! tested there is whether piscem's translation of `DecoderPoolStats` and
//! `DecoderStats` into `ProducerPressure` matches what those states actually
//! mean — which is the judgement-heavy half, since upstream reports no
//! source-bound state and it has to be inferred.
//!
//! ```text
//! cargo test --release --features rapidgzip --test broker_adapters -- --ignored --nocapture
//! ```

#![cfg(feature = "rapidgzip")]

use std::io::Read;
use std::path::PathBuf;

use piscem_rs::io::broker::DecodeProducer;
use rapidgzip_core::{Decoder, DecoderPool};
use thread_broker::{Producer, ProducerPressure, Work};

/// Multi-member; parallelises well.
const DENSE: &str = "/scratch1/rob/flex_bench/mm2/n1_R1_0.gz";
/// Below `MINIMUM_GRID_TASKS`, so `rapidgzip` commits to its sequential
/// backend: the input that can never use more concurrency.
const TINY: &str = "/scratch3/rob/tmp/claude-13757/-scratch1-rob-alevin-fry-ecosystem/580bacef-3758-4e5d-b3a2-798ffd2fb3de/scratchpad/fx/tiny.fq.gz";

fn fixture(env: &str, default: &str) -> Option<PathBuf> {
    let p = PathBuf::from(std::env::var(env).unwrap_or_else(|_| default.to_string()));
    if p.is_file() {
        Some(p)
    } else {
        eprintln!("skipping: {} absent (set {env})", p.display());
        None
    }
}

/// Open `path` under a pool, read `mib`, and report what the adapter said while
/// decoding was in flight.
fn observe(path: &PathBuf, pool_workers: usize, limit: usize, mib: usize) -> Vec<ProducerPressure> {
    observe_with_work(path, pool_workers, limit, mib).0
}

/// As [`observe`], also returning the work counters the broker would have seen.
fn observe_with_work(
    path: &PathBuf,
    pool_workers: usize,
    limit: usize,
    mib: usize,
) -> (Vec<ProducerPressure>, Vec<Work>) {
    let pool = DecoderPool::builder()
        .workers(pool_workers)
        .initial_worker_limit(limit)
        .build()
        .expect("pool");
    let decoder = Decoder::builder()
        .decoder_threads(pool_workers)
        .decoder_pool(pool.clone())
        .build()
        .expect("decoder");
    let mut reader = decoder.open(path).expect("open");
    let producer = DecodeProducer::new(pool, vec![reader.handle()]).unwrap();

    let mut buf = vec![0u8; 1 << 20];
    let mut total = 0u64;
    let mut seen = Vec::new();
    let mut work = Vec::new();
    while total < (mib as u64) << 20 {
        match reader.read(&mut buf) {
            Ok(0) => break,
            Ok(n) => total += n as u64,
            Err(e) => panic!("read: {e}"),
        }
        seen.push(producer.pressure());
        work.push(producer.work());
    }
    (seen, work)
}

/// A decoder held at one slot against a file with plenty of parallel work must
/// report starvation, or the broker will never grow it.
#[test]
#[ignore = "needs the rapidgzip feature and a real fixture"]
fn a_throttled_decoder_reports_starvation() {
    let Some(f) = fixture("PISCEM_TEST_DENSE_GZ", DENSE) else {
        return;
    };
    let seen = observe(&f, 16, 1, 128);

    let starved = seen
        .iter()
        .filter(|p| matches!(p, ProducerPressure::Starved))
        .count();
    assert!(
        starved > 0,
        "a decoder throttled to one slot never reported starvation: {seen:?}"
    );
}

/// The one measurement the whole control law rests on.
///
/// `busy_nanos` must be *decode time with blocking excluded*, integrated by
/// rapidgzip's native executing-region counter. Two properties make it usable,
/// and neither is guaranteed by the type:
///
/// * it must actually accumulate — a counter stuck at zero silently reads as
///   "decoding is free", and the model would then hand the entire budget to
///   mapping on every workload;
/// * it must never exceed `limit x elapsed`, or the integration is
///   double-counting slots and the model would over-buy decode concurrency in
///   proportion to the error.
#[test]
#[ignore = "needs the rapidgzip feature and a real fixture"]
fn producer_work_is_bounded_busy_time() {
    let Some(f) = fixture("PISCEM_TEST_DENSE_GZ", DENSE) else {
        return;
    };
    let limit = 4;
    let started = std::time::Instant::now();
    let (_, work) = observe_with_work(&f, 16, limit, 128);
    let elapsed = started.elapsed().as_nanos() as u64;
    let last = *work.last().expect("at least one sample");

    assert!(
        last.busy_nanos > 0,
        "decode busy time never accumulated; the model would price decoding at zero",
    );
    assert!(
        last.items > 0,
        "no decompressed bytes reported: {:?}",
        last.items,
    );
    // Generous: `observe` starts its clock before the decoder does, so the true
    // bound is tighter than this. It still catches a slot-counting error, which
    // would overshoot by a whole multiple.
    assert!(
        last.busy_nanos <= elapsed * limit as u64,
        "busy time {} exceeds {} slots x {} ns elapsed",
        last.busy_nanos,
        limit,
        elapsed,
    );
    assert!(
        work.windows(2)
            .all(|w| w[0].busy_nanos <= w[1].busy_nanos && w[0].items <= w[1].items),
        "work counters went backwards; the broker takes deltas and would read \
         a huge negative as zero",
    );
}

/// The state that dominates all others: an input the decoder can only process
/// serially reports total starvation forever and can never be helped.
#[test]
#[ignore = "needs the rapidgzip feature and a real fixture"]
fn a_sequential_input_reports_inelastic() {
    let Some(f) = fixture("PISCEM_TEST_TINY_GZ", TINY) else {
        return;
    };
    let seen = observe(&f, 16, 8, 8);

    assert!(
        seen.contains(&ProducerPressure::Inelastic),
        "a sequential-only input was not reported inelastic: {seen:?}"
    );
    assert!(
        !seen.iter().any(|p| matches!(p, ProducerPressure::Starved)),
        "an unparallelisable input asked for more threads: {seen:?}"
    );
}

/// A decoder that gets ahead of its consumer must report `Satisfied`, or the
/// broker keeps buying decode capacity that will idle.
///
/// The consumer here deliberately dawdles. A bare `read()` loop is an
/// infinitely fast consumer, so decode is *always* the constraint against one
/// and a decoder at its full limit correctly reports starvation — an earlier
/// version of this test asserted the opposite and was simply wrong about what
/// it had built.
#[test]
#[ignore = "needs the rapidgzip feature and a real fixture"]
fn a_decoder_ahead_of_its_consumer_reports_satisfied() {
    let Some(f) = fixture("PISCEM_TEST_DENSE_GZ", DENSE) else {
        return;
    };

    let pool = DecoderPool::builder()
        .workers(16)
        .initial_worker_limit(16)
        .build()
        .expect("pool");
    let decoder = Decoder::builder()
        .decoder_threads(16)
        .decoder_pool(pool.clone())
        .build()
        .expect("decoder");
    let mut reader = decoder.open(&f).expect("open");
    let producer = DecodeProducer::new(pool, vec![reader.handle()]).unwrap();

    let mut buf = vec![0u8; 1 << 20];
    let mut seen = Vec::new();
    for _ in 0..40 {
        // Slow enough that 16 decode workers comfortably outrun it.
        std::thread::sleep(std::time::Duration::from_millis(20));
        match reader.read(&mut buf) {
            Ok(0) => break,
            Ok(_) => {}
            Err(e) => panic!("read: {e}"),
        }
        seen.push(producer.pressure());
    }

    let satisfied = seen
        .iter()
        .filter(|p| **p == ProducerPressure::Satisfied)
        .count();
    assert!(
        satisfied > 0,
        "a decoder outrunning its consumer never reported satisfied: {seen:?}"
    );
    assert!(
        !seen.contains(&ProducerPressure::Inelastic),
        "a parallelisable input was called inelastic: {seen:?}"
    );
}

/// The limit round-trips, so the broker's arithmetic and the pool's agree.
#[test]
fn the_limit_round_trips() {
    let pool = DecoderPool::builder()
        .workers(24)
        .initial_worker_limit(4)
        .build()
        .expect("pool");
    let producer = DecodeProducer::new(pool, Vec::new()).unwrap();

    assert_eq!(producer.limit(), 4);
    producer.set_limit(12).expect("valid limit");
    assert_eq!(producer.limit(), 12);

    // Above the immutable maximum is refused, and must leave the limit intact
    // rather than silently clamping to something the broker did not choose.
    assert!(producer.set_limit(999).is_err());
    assert_eq!(
        producer.limit(),
        12,
        "a refused limit changed the pool anyway"
    );

    // With no decoders attached there is nothing to size.
    assert_eq!(producer.pressure(), ProducerPressure::Satisfied);
    assert_eq!(producer.active_slots(), 0);
}
