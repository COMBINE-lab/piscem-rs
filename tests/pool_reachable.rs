// Phase 1 gate: the shared-pool API is reachable from piscem-rs.
#![cfg(feature = "rapidgzip")]
#[test]
fn the_shared_pool_api_is_reachable() {
    let pool = rapidgzip_core::DecoderPool::builder()
        .workers(24)
        .initial_worker_limit(8)
        .build()
        .expect("build pool");

    let s = pool.stats();
    assert_eq!(s.configured_workers, 24);
    assert_eq!(s.worker_limit, 8);
    assert_eq!(s.attached_decoders, 0);

    // The mutable aggregate ceiling: this is what the scheduler will drive.
    pool.set_worker_limit(16).expect("raise limit");
    assert_eq!(pool.stats().worker_limit, 16);
    assert!(
        pool.set_worker_limit(25).is_err(),
        "must not exceed configured"
    );

    // A decoder can be attached to it.
    let _decoder = rapidgzip_core::Decoder::builder()
        .decoder_threads(24)
        .decoder_pool(pool.clone())
        .build()
        .expect("build decoder with pool");
}
