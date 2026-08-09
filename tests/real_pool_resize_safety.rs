//! Budget and correctness checks with both real resizable pools.

#![cfg(feature = "rapidgzip")]

use std::io::Write;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::time::{Duration, Instant};

use flate2::Compression;
use flate2::write::GzEncoder;
use paraseq::fastx::{Collection, CollectionType, Reader, RefRecord};
use paraseq::parallel::ThreadPool;
use paraseq::prelude::*;
use piscem_rs::io::broker::DecodeProducer;
use rapidgzip_core::{Decoder, DecoderPool};
use thread_broker::Producer;

const BUDGET: usize = 16;
const RECORDS_PER_GROUP: usize = 3_200;
const RECORDS_PER_MEMBER: usize = 64;

fn write_fixture(path: &std::path::Path, group: usize) {
    let mut file = std::fs::File::create(path).unwrap();
    for first in (0..RECORDS_PER_GROUP).step_by(RECORDS_PER_MEMBER) {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
        for record in first..(first + RECORDS_PER_MEMBER).min(RECORDS_PER_GROUP) {
            let len = 120 + (group + record) % 31;
            writeln!(encoder, "@g{group}-r{record}").unwrap();
            for base in 0..len {
                encoder
                    .write_all(&[b"ACGT"[(group + record + base) % 4]])
                    .unwrap();
            }
            encoder.write_all(b"\n+\n").unwrap();
            encoder.write_all(&vec![b'I'; len]).unwrap();
            encoder.write_all(b"\n").unwrap();
        }
        file.write_all(&encoder.finish().unwrap()).unwrap();
    }
}

fn process_cpu_seconds() -> f64 {
    let mut usage = std::mem::MaybeUninit::<libc::rusage>::uninit();
    // SAFETY: `getrusage` initializes the supplied `rusage` on success and the
    // pointer is valid for the duration of the call.
    let rc = unsafe { libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()) };
    assert_eq!(rc, 0, "getrusage failed");
    // SAFETY: success above guarantees initialization.
    let usage = unsafe { usage.assume_init() };
    let timeval = |time: libc::timeval| time.tv_sec as f64 + time.tv_usec as f64 / 1e6;
    timeval(usage.ru_utime) + timeval(usage.ru_stime)
}

#[derive(Debug)]
struct RunStats {
    records: u64,
    bases: u64,
    peak_slots: usize,
    cpu_concurrency: f64,
    decode_busy_seconds: f64,
    process_cpu_seconds: f64,
}

fn run(paths: &[std::path::PathBuf], churn: bool) -> RunStats {
    let opening_mapping = BUDGET - 1;
    let floor = piscem_rs::io::fastx::collection_share_floor(paths.len(), 1, opening_mapping);
    let map_pool = ThreadPool::with_max(opening_mapping, opening_mapping);
    let decode_pool = DecoderPool::builder()
        .workers(BUDGET)
        .initial_worker_limit(1)
        .build()
        .unwrap();

    let mut handles = Vec::new();
    let mut readers = Vec::new();
    for path in paths {
        let decoder = Decoder::builder()
            .decoder_threads(BUDGET)
            .decoder_pool(decode_pool.clone())
            .build()
            .unwrap();
        let reader = decoder.open(path).unwrap();
        handles.push(reader.handle());
        readers.push(Reader::new_with_batch_size(reader, 32).unwrap());
    }
    let producer = Arc::new(DecodeProducer::new(decode_pool, handles).unwrap());

    let done = Arc::new(AtomicBool::new(false));
    let peak = Arc::new(AtomicUsize::new(0));
    let violation = Arc::new(AtomicBool::new(false));
    let monitor = {
        let map_pool = map_pool.clone();
        let producer = Arc::clone(&producer);
        let done = Arc::clone(&done);
        let peak = Arc::clone(&peak);
        let violation = Arc::clone(&violation);
        std::thread::spawn(move || {
            while !done.load(Ordering::Acquire) {
                let occupied = map_pool
                    .total_live()
                    .saturating_add(producer.active_slots());
                peak.fetch_max(occupied, Ordering::AcqRel);
                if occupied > BUDGET {
                    violation.store(true, Ordering::Release);
                }
                std::thread::yield_now();
            }
        })
    };

    let resizer = churn.then(|| {
        let map_pool = map_pool.clone();
        let producer = Arc::clone(&producer);
        let done = Arc::clone(&done);
        std::thread::spawn(move || {
            while map_pool.total_live() == 0 && !done.load(Ordering::Acquire) {
                std::thread::yield_now();
            }
            while !done.load(Ordering::Acquire) {
                map_pool.set_threads(floor);
                let deadline = Instant::now() + Duration::from_secs(2);
                while map_pool.total_live() > floor && !done.load(Ordering::Acquire) {
                    assert!(Instant::now() < deadline, "mapping shrink timed out");
                    std::thread::yield_now();
                }
                producer.set_limit(BUDGET - floor).unwrap();
                std::thread::sleep(Duration::from_millis(100));

                producer.set_limit(1).unwrap();
                let deadline = Instant::now() + Duration::from_secs(2);
                while producer.active_slots() > 1 && !done.load(Ordering::Acquire) {
                    assert!(Instant::now() < deadline, "decode shrink timed out");
                    std::thread::yield_now();
                }
                map_pool.set_threads(opening_mapping);
                std::thread::sleep(Duration::from_millis(100));
            }
        })
    });

    let count = Arc::new(AtomicU64::new(0));
    let bases = Arc::new(AtomicU64::new(0));
    let mut processor = {
        let count = Arc::clone(&count);
        let bases = Arc::clone(&bases);
        move |batch: &mut dyn Iterator<Item = RefRecord<'_>>| {
            let mut local_count = 0u64;
            let mut local_bases = 0u64;
            for record in batch {
                local_count += 1;
                local_bases += record.seq().len() as u64;
            }
            count.fetch_add(local_count, Ordering::Relaxed);
            bases.fetch_add(local_bases, Ordering::Relaxed);
            std::thread::sleep(Duration::from_millis(2));
            Ok(())
        }
    };

    let cpu_start = process_cpu_seconds();
    let wall_start = Instant::now();
    Collection::new(readers, CollectionType::Single)
        .unwrap()
        .process_parallel_pool(&mut processor, &map_pool, None)
        .unwrap();
    done.store(true, Ordering::Release);
    if let Some(handle) = resizer {
        handle.join().unwrap();
    }
    monitor.join().unwrap();
    let process_cpu = process_cpu_seconds() - cpu_start;
    let cpu_concurrency = process_cpu / wall_start.elapsed().as_secs_f64();
    let decode_busy_seconds = producer.work().busy_nanos as f64 / 1e9;
    assert!(
        !violation.load(Ordering::Acquire),
        "real pools exceeded {BUDGET} controlled slots; peak was {}",
        peak.load(Ordering::Acquire)
    );

    RunStats {
        records: count.load(Ordering::Relaxed),
        bases: bases.load(Ordering::Relaxed),
        peak_slots: peak.load(Ordering::Relaxed),
        cpu_concurrency,
        decode_busy_seconds,
        process_cpu_seconds: process_cpu,
    }
}

#[test]
fn real_mapping_and_decode_pools_preserve_budget_and_output_under_churn() {
    let dir = tempfile::tempdir().unwrap();
    let paths: Vec<_> = (0..24)
        .map(|group| {
            let path = dir.path().join(format!("group-{group}.fq.gz"));
            write_fixture(&path, group);
            path
        })
        .collect();

    for groups in [1, 2, 8, 24] {
        let fixed = run(&paths[..groups], false);
        let churned = run(&paths[..groups], true);
        assert_eq!(churned.records, (groups * RECORDS_PER_GROUP) as u64);
        assert_eq!(
            (churned.records, churned.bases),
            (fixed.records, fixed.bases),
            "real-pool churn changed output for {groups} logical inputs"
        );
        assert!(churned.peak_slots <= BUDGET);
        // One CPU above the controlled budget is the explicit allowance for
        // decoder coordinators, collection driver threads, and this test's
        // monitor/resizer. This is at least the audit's max(1, 3% of N) rule.
        assert!(
            churned.cpu_concurrency <= (BUDGET + 1) as f64,
            "average process CPU concurrency {:.2} exceeded the {}-slot budget plus one auxiliary CPU",
            churned.cpu_concurrency,
            BUDGET,
        );
        assert!(
            churned.decode_busy_seconds <= churned.process_cpu_seconds * 1.05,
            "sampled decode busy time cannot exceed process CPU ground truth: {churned:?}",
        );
    }
}
