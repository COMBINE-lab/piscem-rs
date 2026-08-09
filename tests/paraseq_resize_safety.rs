//! Adversarial resize checks against paraseq's real shared worker pool.

use std::io::Cursor;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::time::{Duration, Instant};

use paraseq::fastx::{Collection, CollectionType, Reader, RefRecord};
use paraseq::parallel::ThreadPool;
use paraseq::prelude::*;

const BUDGET: usize = 32;
const RECORDS_PER_GROUP: usize = 6_400;

fn input(group: usize) -> Cursor<Vec<u8>> {
    let mut bytes = Vec::new();
    for record in 0..RECORDS_PER_GROUP {
        let len = 20 + (group + record) % 17;
        bytes.extend_from_slice(format!("@g{group}-r{record}\n").as_bytes());
        bytes.extend(std::iter::repeat_n(b'A', len));
        bytes.extend_from_slice(b"\n+\n");
        bytes.extend(std::iter::repeat_n(b'I', len));
        bytes.push(b'\n');
    }
    Cursor::new(bytes)
}

fn run(groups: usize, churn: bool) -> (u64, u64, usize) {
    let opening_mapping = BUDGET - 1;
    let floor = piscem_rs::io::fastx::collection_share_floor(groups, 1, opening_mapping);
    let pool = ThreadPool::with_max(opening_mapping, opening_mapping);
    let producer_active = Arc::new(AtomicUsize::new(1));
    let done = Arc::new(AtomicBool::new(false));
    let peak = Arc::new(AtomicUsize::new(0));
    let violation = Arc::new(AtomicBool::new(false));

    let monitor = {
        let pool = pool.clone();
        let producer_active = Arc::clone(&producer_active);
        let done = Arc::clone(&done);
        let peak = Arc::clone(&peak);
        let violation = Arc::clone(&violation);
        std::thread::spawn(move || {
            while !done.load(Ordering::Acquire) {
                let occupied = pool
                    .total_live()
                    .saturating_add(producer_active.load(Ordering::Acquire));
                peak.fetch_max(occupied, Ordering::AcqRel);
                if occupied > BUDGET {
                    violation.store(true, Ordering::Release);
                }
                std::thread::yield_now();
            }
        })
    };

    let resizer = churn.then(|| {
        let pool = pool.clone();
        let producer_active = Arc::clone(&producer_active);
        let done = Arc::clone(&done);
        std::thread::spawn(move || {
            while pool.total_live() == 0 && !done.load(Ordering::Acquire) {
                std::thread::yield_now();
            }
            while !done.load(Ordering::Acquire) {
                // Consumer shrinks first. Producer growth is forbidden until
                // the real aggregate live count acknowledges the target.
                pool.set_threads(floor);
                let deadline = Instant::now() + Duration::from_secs(2);
                while pool.total_live() > floor && !done.load(Ordering::Acquire) {
                    assert!(Instant::now() < deadline, "consumer shrink timed out");
                    std::thread::yield_now();
                }
                producer_active.store(BUDGET - floor, Ordering::Release);
                std::thread::sleep(Duration::from_millis(100));

                // Reverse direction: producer releases its slots before the
                // mapping target is allowed to grow.
                producer_active.store(1, Ordering::Release);
                pool.set_threads(opening_mapping);
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
            // Keep batches in flight long enough for several 100 ms resize
            // cycles and, importantly, for retirement to lag the request.
            std::thread::sleep(Duration::from_millis(2));
            Ok(())
        }
    };

    let readers = (0..groups)
        .map(|group| Reader::new_with_batch_size(input(group), 32).unwrap())
        .collect();
    Collection::new(readers, CollectionType::Single)
        .unwrap()
        .process_parallel_pool(&mut processor, &pool, None)
        .unwrap();

    done.store(true, Ordering::Release);
    if let Some(handle) = resizer {
        handle.join().unwrap();
    }
    monitor.join().unwrap();
    assert!(
        !violation.load(Ordering::Acquire),
        "actual controlled occupancy exceeded {BUDGET}; peak was {}",
        peak.load(Ordering::Acquire)
    );

    (
        count.load(Ordering::Relaxed),
        bases.load(Ordering::Relaxed),
        peak.load(Ordering::Relaxed),
    )
}

#[test]
fn real_multi_share_pool_preserves_budget_and_records_under_churn() {
    for groups in [1, 2, 8, 24] {
        let fixed = run(groups, false);
        let churned = run(groups, true);
        assert_eq!(churned.0, (groups * RECORDS_PER_GROUP) as u64);
        assert_eq!(
            (churned.0, churned.1),
            (fixed.0, fixed.1),
            "resize churn changed output for {groups} logical groups"
        );
        assert!(churned.2 <= BUDGET);
    }
}
