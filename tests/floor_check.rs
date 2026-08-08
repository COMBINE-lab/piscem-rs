#![cfg(feature = "rapidgzip")]
use rapidgzip_core::{Decoder, DecoderPath, DecoderPool};
use std::io::Read;

/// Does a pool limit of 1 latch a decoder into its sequential backend?
///
/// It did before the shared-pool work: admission read the runtime throttle, so
/// a low limit permanently selected `Sequential`. If admission now reads the
/// configured width instead, the aggregate limit can go to 1 without harm --
/// which decides whether the broker's producer floor has to be per-input or
/// can be a single aggregate slot.
#[test]
#[ignore]
fn a_pool_limit_of_one_does_not_latch_sequential() {
    let f = "/scratch1/rob/flex_bench/mm2/n1_R1_0.gz";
    if !std::path::Path::new(f).is_file() {
        eprintln!("skip");
        return;
    }

    for limit in [1usize, 2, 4] {
        let pool = DecoderPool::builder()
            .workers(32)
            .initial_worker_limit(limit)
            .build()
            .unwrap();
        let decoder = Decoder::builder()
            .decoder_threads(32)
            .decoder_pool(pool.clone())
            .build()
            .unwrap();
        let mut r = decoder.open(f).unwrap();
        let h = r.handle();
        let mut buf = vec![0u8; 1 << 20];
        let mut total = 0u64;
        while total < (64u64 << 20) {
            match r.read(&mut buf) {
                Ok(0) => break,
                Ok(n) => total += n as u64,
                Err(e) => panic!("{e}"),
            }
        }
        let s = h.stats();
        eprintln!(
            "  pool limit {limit}: path={:?} spawned={} desired={}",
            s.path, s.spawned_workers, s.desired_workers
        );
        assert_ne!(
            s.path,
            DecoderPath::Sequential,
            "a pool limit of {limit} latched the decoder sequential"
        );
    }
}
