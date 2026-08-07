//! The control law, against a pipeline whose optimum is known analytically.
//!
//! # Why the fake is coupled
//!
//! Two earlier versions of this file tested nothing. The first had a background
//! thread sleeping 200 us per item, which Linux will not honour. The second let
//! the consumer's idleness be set independently of the producer's capacity — so
//! the broker could move every thread it liked and the "starvation" it was
//! chasing never responded. Both passed.
//!
//! What makes this one a test is that supply and demand are joined:
//!
//! ```text
//!     supply     = producer_threads / s_p     (items/sec, capped by the source)
//!     demand     = consumer_threads / s_c
//!     throughput = min(supply, demand)
//! ```
//!
//! and both busy-time counters are then derived from the *same* throughput, so
//! the broker cannot improve one reading without paying for it in the other.
//! The optimum is where supply meets demand:
//!
//! ```text
//!     d* = N * s_p / (s_p + s_c)
//! ```
//!
//! which every test below asserts against directly rather than against a
//! remembered number.

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use thread_broker::{
    BrokerConfig, BrokerReport, Consumer, Producer, ProducerPressure, ThreadBroker, Work,
};

const NANOS: f64 = 1e9;

/// A two-stage pipeline with known costs.
struct Pipeline {
    /// Producer thread-nanoseconds per item.
    s_p: f64,
    /// Consumer thread-nanoseconds per item.
    s_c: f64,
    /// Ceiling on throughput imposed by something neither side controls —
    /// storage bandwidth, or a stream that cannot be decoded in parallel.
    source_cap: Option<f64>,
    /// Multiplier applied to the producer's *reported* busy time.
    ///
    /// The model is only as good as its inputs; this is how a test makes the
    /// model wrong on purpose, to check that the guards catch it.
    producer_bias: f64,
    state: Mutex<State>,
}

struct State {
    consumer_threads: usize,
    producer_slots: usize,
    last: Instant,
    consumer_busy: f64,
    producer_busy: f64,
    items: f64,
}

impl Pipeline {
    fn new(s_p: f64, s_c: f64, consumer_threads: usize, producer_slots: usize) -> Arc<Self> {
        Arc::new(Self {
            s_p,
            s_c,
            source_cap: None,
            producer_bias: 1.0,
            state: Mutex::new(State {
                consumer_threads,
                producer_slots,
                last: Instant::now(),
                consumer_busy: 0.0,
                producer_busy: 0.0,
                items: 0.0,
            }),
        })
    }

    /// Advance the simulation to now.
    ///
    /// Driven by polling rather than by a thread, so it stays exact regardless
    /// of how the scheduler treats the test process: every call integrates over
    /// however much real time actually elapsed.
    fn advance(&self) {
        let mut st = self.state.lock().unwrap();
        let now = Instant::now();
        let dt = now.saturating_duration_since(st.last).as_secs_f64();
        st.last = now;
        if dt <= 0.0 {
            return;
        }

        let supply = st.producer_slots as f64 * NANOS / self.s_p;
        let supply = match self.source_cap {
            Some(cap) => supply.min(cap),
            None => supply,
        };
        let demand = st.consumer_threads as f64 * NANOS / self.s_c;
        let throughput = supply.min(demand);

        let items = throughput * dt;
        st.items += items;
        // Both derived from the same throughput: this is the coupling.
        st.consumer_busy += items * self.s_c;
        st.producer_busy += items * self.s_p * self.producer_bias;
    }

    /// The analytic answer, for tests to assert against.
    fn optimum(&self, budget: usize) -> usize {
        (budget as f64 * self.s_p / (self.s_p + self.s_c)).round() as usize
    }

    fn split(&self) -> (usize, usize) {
        let st = self.state.lock().unwrap();
        (st.consumer_threads, st.producer_slots)
    }
}

struct FakeConsumer(Arc<Pipeline>);
struct FakeProducer(Arc<Pipeline>);

impl Consumer for FakeConsumer {
    fn set_threads(&self, n: usize) {
        self.0.advance();
        self.0.state.lock().unwrap().consumer_threads = n;
    }
    fn live_threads(&self) -> usize {
        self.0.state.lock().unwrap().consumer_threads
    }
    fn work(&self) -> Work {
        self.0.advance();
        let st = self.0.state.lock().unwrap();
        Work {
            busy_nanos: st.consumer_busy as u64,
            items: st.items as u64,
        }
    }
}

impl Producer for FakeProducer {
    fn set_limit(&self, n: usize) {
        self.0.advance();
        self.0.state.lock().unwrap().producer_slots = n;
    }
    fn limit(&self) -> usize {
        self.0.state.lock().unwrap().producer_slots
    }
    fn work(&self) -> Work {
        self.0.advance();
        let st = self.0.state.lock().unwrap();
        Work {
            busy_nanos: st.producer_busy as u64,
            items: st.items as u64,
        }
    }
    fn pressure(&self) -> ProducerPressure {
        let st = self.0.state.lock().unwrap();
        let uncapped = st.producer_slots as f64 * NANOS / self.0.s_p;
        let demand = st.consumer_threads as f64 * NANOS / self.0.s_c;
        match self.0.source_cap {
            // The source is the constraint and more slots cannot change that.
            Some(cap) if uncapped > cap => ProducerPressure::SourceBound,
            _ if uncapped < demand => ProducerPressure::Starved,
            _ => ProducerPressure::Satisfied,
        }
    }
}

/// Fast enough to finish a test in about a second, slow enough that the
/// smoothing and blackout windows still mean something.
fn quick() -> BrokerConfig {
    BrokerConfig {
        sample_interval: Duration::from_millis(20),
        warmup: Duration::from_millis(60),
        smoothing_windows: 3,
        blackout_samples: 3,
        ratify_samples: 3,
        ..BrokerConfig::default()
    }
}

fn run_for(
    pipeline: Arc<Pipeline>,
    budget: usize,
    cfg: BrokerConfig,
    dur: Duration,
) -> BrokerReport {
    run_with(pipeline, budget, cfg, dur, |_| {})
}

fn run_with(
    pipeline: Arc<Pipeline>,
    budget: usize,
    cfg: BrokerConfig,
    dur: Duration,
    during: impl FnOnce(&Arc<Pipeline>),
) -> BrokerReport {
    let start = pipeline.split().1;
    let broker = ThreadBroker::builder_with(
        FakeConsumer(Arc::clone(&pipeline)),
        FakeProducer(Arc::clone(&pipeline)),
        cfg,
    )
    .budget(budget)
    .initial_producer_slots(start)
    .build()
    .expect("valid config")
    .start();

    std::thread::sleep(dur / 2);
    during(&pipeline);
    std::thread::sleep(dur / 2);
    broker.finish()
}

/// Within the deadband, plus one for the rounding in the model itself.
fn assert_near(actual: usize, want: usize, cfg: &BrokerConfig, what: &str) {
    let slack = cfg.deadband_threads + 1;
    assert!(
        actual.abs_diff(want) <= slack,
        "{what}: settled at {actual}, expected {want} +/- {slack}",
    );
}

#[test]
fn solves_a_balanced_split() {
    let cfg = quick();
    // Equal per-item cost, so the budget should divide in half.
    let p = Pipeline::new(1_000.0, 1_000.0, 28, 4);
    let want = p.optimum(32);
    assert_eq!(want, 16);
    let report = run_for(Arc::clone(&p), 32, cfg, Duration::from_millis(1200));
    assert_near(report.final_producer_limit, want, &cfg, "balanced");
}

#[test]
fn solves_a_decode_heavy_split() {
    let cfg = quick();
    // Decoding costs 3x what consuming does, so it should get 3/4 of the budget.
    let p = Pipeline::new(3_000.0, 1_000.0, 28, 4);
    let want = p.optimum(32);
    assert_eq!(want, 24);
    let report = run_for(Arc::clone(&p), 32, cfg, Duration::from_millis(1200));
    assert_near(report.final_producer_limit, want, &cfg, "decode-heavy");
}

#[test]
fn solves_a_mapping_heavy_split() {
    let cfg = quick();
    // The common case: mapping dominates and decode needs very little.
    let p = Pipeline::new(1_000.0, 7_000.0, 28, 4);
    let want = p.optimum(32);
    assert_eq!(want, 4);
    let report = run_for(Arc::clone(&p), 32, cfg, Duration::from_millis(1200));
    assert_near(report.final_producer_limit, want, &cfg, "mapping-heavy");
}

/// The regression test for the bug this control law exists to fix.
///
/// Started *at* the optimum, the previous law still walked away from it: it
/// drove consumer starvation to zero, and starvation stays at zero for every
/// split above the balance point, so nothing pulled back. It settled at 47
/// producer slots of 64 where 32 was right, 44% slower, reporting near-zero
/// error the whole way.
#[test]
fn does_not_walk_away_from_the_optimum() {
    let cfg = quick();
    let want = 32;
    let p = Pipeline::new(1_000.0, 1_000.0, 64 - want, want);
    assert_eq!(p.optimum(64), want);
    let report = run_for(Arc::clone(&p), 64, cfg, Duration::from_millis(1200));
    assert_near(report.final_producer_limit, want, &cfg, "held at optimum");
    assert!(
        report.producer_trajectory.iter().all(|&d| d <= want + 4),
        "drifted above the optimum: {:?}",
        report.producer_trajectory,
    );
}

/// Solving means jumping, so the distance to the answer must not set the cost.
#[test]
fn reaches_a_distant_optimum_in_a_handful_of_moves() {
    let cfg = quick();
    // Start at 2 of 64 with the answer at 48: a 46-thread journey, which the
    // old one-step-at-a-time law could not finish inside a short run at all.
    let p = Pipeline::new(3_000.0, 1_000.0, 62, 2);
    let want = p.optimum(64);
    assert_eq!(want, 48);
    let report = run_for(Arc::clone(&p), 64, cfg, Duration::from_millis(1500));
    assert_near(report.final_producer_limit, want, &cfg, "distant optimum");
    assert!(
        report.moves <= 3,
        "took {} moves to cross 46 threads",
        report.moves,
    );
}

/// A producer that cannot use more slots must not be given them, however
/// expensive per item it looks. This is FLINK-31215 in miniature.
#[test]
fn does_not_buy_slots_the_source_cannot_feed() {
    let cfg = quick();
    let mut pipeline = Pipeline::new(4_000.0, 1_000.0, 60, 4);
    // The unclamped model wants 51 of 64. The source can only sustain what 8
    // decode threads would produce, so everything above that is waste.
    let cap = 8.0 * NANOS / 4_000.0;
    Arc::get_mut(&mut pipeline).unwrap().source_cap = Some(cap);
    assert_eq!(pipeline.optimum(64), 51);

    let report = run_for(Arc::clone(&pipeline), 64, cfg, Duration::from_millis(1500));
    assert!(
        report.final_producer_limit <= 12,
        "bought {} slots against a source that can feed 8; cap was {:?}",
        report.final_producer_limit,
        report.final_model.map(|m| m.useful_cap),
    );
    assert!(
        report.source_bound_samples > 0,
        "never noticed the source was the constraint",
    );
}

/// The model can be wrong. When it is, throughput has to be what settles it.
#[test]
fn reverts_a_move_that_makes_throughput_worse() {
    let cfg = quick();
    // Mapping dominates, so almost the whole budget belongs to the consumer --
    // but the producer over-reports its busy time by 20x, so the model asks for
    // most of the budget instead. Applying that starves the consumer and
    // throughput falls, which is the only signal that can catch this.
    let mut pipeline = Pipeline::new(200.0, 8_000.0, 30, 2);
    Arc::get_mut(&mut pipeline).unwrap().producer_bias = 20.0;

    let report = run_for(Arc::clone(&pipeline), 32, cfg, Duration::from_millis(2000));
    assert!(report.reverts > 0, "kept a move that cost throughput");
    assert!(
        report.final_producer_limit <= 8,
        "left the split at {} producer slots",
        report.final_producer_limit,
    );
}

/// Once it has settled it must stay settled, or the moves cost more than the
/// answer is worth: every one rebuilds per-thread state on both sides.
#[test]
fn settles_and_stops_moving() {
    let cfg = quick();
    let p = Pipeline::new(1_000.0, 3_000.0, 24, 8);
    let report = run_for(Arc::clone(&p), 32, cfg, Duration::from_millis(2000));
    assert!(report.converged(), "never settled: {report:?}");
    assert!(
        report.moves <= 3,
        "still moving after settling: {} moves, trajectory {:?}",
        report.moves,
        report.producer_trajectory,
    );

    // The tail of the trajectory is flat.
    let tail = &report.producer_trajectory[report.producer_trajectory.len() / 2..];
    assert!(
        tail.windows(2).all(|w| w[0] == w[1]),
        "kept moving in the second half: {tail:?}",
    );
}

/// Settling must not mean going deaf. A file ending or a different sample
/// starting changes the costs, and the split has to follow.
#[test]
fn re_solves_after_a_regime_change() {
    let cfg = quick();
    // Starts mapping-heavy, so decode gets little.
    let p = Pipeline::new(1_000.0, 7_000.0, 28, 4);
    let report = run_with(Arc::clone(&p), 32, cfg, Duration::from_millis(3000), |p| {
        // Halfway through, decode becomes the expensive side. `s_c` is only
        // read under the lock inside `advance`, so mutating it here is
        // exactly the abrupt change a new input would cause.
        let ptr = Arc::as_ptr(p) as *mut Pipeline;
        // SAFETY: the broker touches `s_c` only through `&Pipeline` reads
        // that are already serialised by `advance`'s lock, and the test
        // owns the only other reference.
        unsafe {
            (*ptr).s_c = 300.0;
        }
    });
    assert!(
        report.resurveys > 0 || report.final_producer_limit > 8,
        "ignored the regime change: settled at {} with {} resurveys",
        report.final_producer_limit,
        report.resurveys,
    );
}

/// A run whose inputs need no decoding at all must not reserve slots for it.
#[test]
fn hands_the_whole_budget_over_when_there_is_nothing_to_decode() {
    let cfg = quick();
    // A producer that costs essentially nothing per item: uncompressed input.
    let p = Pipeline::new(0.01, 5_000.0, 28, 4);
    let report = run_for(Arc::clone(&p), 32, cfg, Duration::from_millis(1200));
    assert_eq!(
        report.final_producer_limit, cfg.min_producer_slots,
        "reserved decode slots for a pipeline with no decode cost",
    );
    assert_eq!(report.final_consumer_threads, 32 - cfg.min_producer_slots);
}

/// Floors are a contract with the caller, not a suggestion: below them a
/// `rapidgzip` decoder can commit irreversibly to its sequential backend.
#[test]
fn never_breaches_either_floor() {
    let cfg = BrokerConfig {
        min_consumer_threads: 4,
        min_producer_slots: 3,
        ..quick()
    };
    // Wants everything on the consumer; the producer floor must survive it.
    let p = Pipeline::new(1.0, 9_000.0, 12, 4);
    let report = run_for(Arc::clone(&p), 16, cfg, Duration::from_millis(1200));
    assert!(
        report
            .producer_trajectory
            .iter()
            .all(|&d| d >= cfg.min_producer_slots),
        "breached the producer floor: {:?}",
        report.producer_trajectory,
    );
    assert!(
        report
            .consumer_trajectory
            .iter()
            .all(|&c| c >= cfg.min_consumer_threads),
        "breached the consumer floor: {:?}",
        report.consumer_trajectory,
    );
}

/// The budget is the point. Two sides summing above it is oversubscription on a
/// full machine, which is what the whole crate exists to prevent.
#[test]
fn never_oversubscribes_the_budget() {
    let cfg = quick();
    let p = Pipeline::new(2_000.0, 1_000.0, 28, 4);
    let report = run_for(Arc::clone(&p), 32, cfg, Duration::from_millis(1500));
    for (c, d) in report
        .consumer_trajectory
        .iter()
        .zip(&report.producer_trajectory)
    {
        assert_eq!(c + d, 32, "split summed to {} rather than 32", c + d);
    }
}

/// Discarding the startup ramp is load-bearing: every workload measured reported
/// ~99% consumer starvation in its first 75 ms, purely because its threads were
/// still starting.
#[test]
fn makes_no_decision_during_warm_up() {
    let cfg = BrokerConfig {
        sample_interval: Duration::from_millis(20),
        warmup: Duration::from_millis(400),
        smoothing_windows: 3,
        blackout_samples: 3,
        ratify_samples: 3,
        ..BrokerConfig::default()
    };
    let p = Pipeline::new(3_000.0, 1_000.0, 30, 2);
    let report = run_for(Arc::clone(&p), 32, cfg, Duration::from_millis(300));
    assert_eq!(report.moves, 0, "moved during warm-up");
    assert!(report.producer_trajectory.is_empty());
    assert_eq!(report.final_producer_limit, 2);
}

/// The defaults were once the one combination no test exercised, and they did
/// not satisfy the crate's own validation rules.
#[test]
fn the_default_configuration_builds() {
    let p = Pipeline::new(1_000.0, 1_000.0, 30, 2);
    ThreadBroker::builder(FakeConsumer(Arc::clone(&p)), FakeProducer(p))
        .budget(32)
        .build()
        .expect("the default configuration must be valid");
}

#[test]
fn rejects_a_blackout_shorter_than_the_smoothing_window() {
    let p = Pipeline::new(1_000.0, 1_000.0, 30, 2);
    let cfg = BrokerConfig {
        smoothing_windows: 4,
        blackout_samples: 2,
        ..quick()
    };
    assert!(
        ThreadBroker::builder_with(FakeConsumer(Arc::clone(&p)), FakeProducer(p), cfg)
            .budget(32)
            .build()
            .is_err(),
        "a blackout shorter than the smoothing window must be rejected",
    );
}

/// `BusyMeter` exists because a counter published once per batch cannot be
/// sampled faster than a batch completes — and a window with no completed batch
/// reads as zero work, which is maximum starvation for a thread running flat out.
#[test]
fn the_meter_publishes_within_a_batch() {
    let meter = thread_broker::BusyMeter::default();
    let observed = Arc::new(AtomicU64::new(0));

    let mut timer = meter.timer().flush_every(4);
    for i in 1..=8u64 {
        timer.tick();
        std::thread::sleep(Duration::from_millis(1));
        if i == 4 {
            observed.store(meter.work().busy_nanos, Ordering::Relaxed);
        }
    }
    drop(timer);

    assert!(
        observed.load(Ordering::Relaxed) > 0,
        "nothing was published until the batch ended",
    );
    let done = meter.work();
    assert_eq!(done.items, 8, "the tail of the batch was dropped");
    assert!(done.busy_nanos >= observed.load(Ordering::Relaxed));
}
