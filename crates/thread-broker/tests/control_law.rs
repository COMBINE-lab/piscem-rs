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

use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use thread_broker::{
    BrokerConfig, BrokerError, BrokerErrorKind, BrokerReport, Consumer, Producer, ProducerPressure,
    ResizeSide, SteadyStatePolicy, ThreadBroker, Work,
};

const NANOS: f64 = 1e9;

/// A two-stage pipeline with known costs.
struct Pipeline {
    state: Mutex<State>,
}

struct State {
    /// Producer and consumer thread-nanoseconds per item.
    s_p: f64,
    s_c: f64,
    /// Throughput scales as allocated slots raised to this exponent. Values
    /// below one model cache, memory-bandwidth, and serialized-work limits.
    producer_scaling: f64,
    consumer_scaling: f64,
    /// Ceiling on throughput imposed by something neither side controls.
    source_cap: Option<f64>,
    /// Multiplier applied to the producer's reported busy time.
    producer_bias: f64,
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
            state: Mutex::new(State {
                s_p,
                s_c,
                producer_scaling: 1.0,
                consumer_scaling: 1.0,
                source_cap: None,
                producer_bias: 1.0,
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

        let supply = (st.producer_slots as f64).powf(st.producer_scaling) * NANOS / st.s_p;
        let supply = match st.source_cap {
            Some(cap) => supply.min(cap),
            None => supply,
        };
        let demand = (st.consumer_threads as f64).powf(st.consumer_scaling) * NANOS / st.s_c;
        let throughput = supply.min(demand);

        let items = throughput * dt;
        st.items += items;
        // Both derived from the same throughput: this is the coupling.
        st.consumer_busy +=
            items * st.s_c * (st.consumer_threads as f64).powf(1.0 - st.consumer_scaling);
        st.producer_busy += items
            * st.s_p
            * (st.producer_slots as f64).powf(1.0 - st.producer_scaling)
            * st.producer_bias;
    }

    /// The oracle answer, found over every legal discrete split.
    fn optimum(&self, budget: usize) -> usize {
        let st = self.state.lock().unwrap();
        (1..budget)
            .max_by(|&left, &right| {
                Self::throughput_at(&st, budget, left)
                    .partial_cmp(&Self::throughput_at(&st, budget, right))
                    .unwrap()
            })
            .unwrap()
    }

    fn throughput_at(st: &State, budget: usize, producer_slots: usize) -> f64 {
        let consumer_threads = budget - producer_slots;
        let supply = (producer_slots as f64).powf(st.producer_scaling) * NANOS / st.s_p;
        let supply = st.source_cap.map_or(supply, |cap| supply.min(cap));
        let demand = (consumer_threads as f64).powf(st.consumer_scaling) * NANOS / st.s_c;
        supply.min(demand)
    }

    fn throughput_ratio_to_oracle(&self, budget: usize, producer_slots: usize) -> f64 {
        let st = self.state.lock().unwrap();
        let oracle = Self::throughput_at(&st, budget, self.optimum_locked(&st, budget));
        Self::throughput_at(&st, budget, producer_slots) / oracle
    }

    fn optimum_locked(&self, st: &State, budget: usize) -> usize {
        (1..budget)
            .max_by(|&left, &right| {
                Self::throughput_at(st, budget, left)
                    .partial_cmp(&Self::throughput_at(st, budget, right))
                    .unwrap()
            })
            .unwrap()
    }

    fn split(&self) -> (usize, usize) {
        let st = self.state.lock().unwrap();
        (st.consumer_threads, st.producer_slots)
    }
}

struct FakeConsumer(Arc<Pipeline>);
struct FakeProducer(Arc<Pipeline>);

struct DropAwareProducer {
    pipeline: Arc<Pipeline>,
    dropped: Arc<AtomicBool>,
}

impl Drop for DropAwareProducer {
    fn drop(&mut self) {
        self.dropped.store(true, Ordering::Release);
    }
}

struct BufferedProducer {
    pipeline: Arc<Pipeline>,
    buffered: Arc<AtomicU64>,
}

struct ResettingConsumer {
    pipeline: Arc<Pipeline>,
    calls: AtomicUsize,
}

struct ResettingProducer {
    pipeline: Arc<Pipeline>,
    calls: AtomicUsize,
}

type RefusalRule = fn(target: usize, call: usize) -> bool;

struct FaultyConsumer {
    pipeline: Arc<Pipeline>,
    calls: AtomicUsize,
    refuse: RefusalRule,
}

struct FaultyProducer {
    pipeline: Arc<Pipeline>,
    calls: AtomicUsize,
    refuse: RefusalRule,
}

impl FaultyConsumer {
    fn new(pipeline: Arc<Pipeline>, refuse: RefusalRule) -> Self {
        Self {
            pipeline,
            calls: AtomicUsize::new(0),
            refuse,
        }
    }
}

impl FaultyProducer {
    fn new(pipeline: Arc<Pipeline>, refuse: RefusalRule) -> Self {
        Self {
            pipeline,
            calls: AtomicUsize::new(0),
            refuse,
        }
    }
}

impl Consumer for FakeConsumer {
    fn set_threads(&self, n: usize) -> Result<(), thread_broker::ResizeError> {
        self.0.advance();
        self.0.state.lock().unwrap().consumer_threads = n;
        Ok(())
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
    fn set_limit(&self, n: usize) -> Result<(), thread_broker::ResizeError> {
        self.0.advance();
        self.0.state.lock().unwrap().producer_slots = n;
        Ok(())
    }
    fn limit(&self) -> usize {
        self.0.state.lock().unwrap().producer_slots
    }
    fn active_slots(&self) -> usize {
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
        let uncapped = (st.producer_slots as f64).powf(st.producer_scaling) * NANOS / st.s_p;
        let demand = (st.consumer_threads as f64).powf(st.consumer_scaling) * NANOS / st.s_c;
        match st.source_cap {
            // The source is the constraint and more slots cannot change that.
            Some(cap) if uncapped > cap => ProducerPressure::SourceBound,
            _ if uncapped < demand => ProducerPressure::Starved,
            _ => ProducerPressure::Satisfied,
        }
    }
}

impl Producer for DropAwareProducer {
    fn set_limit(&self, n: usize) -> Result<(), thread_broker::ResizeError> {
        FakeProducer(Arc::clone(&self.pipeline)).set_limit(n)
    }

    fn limit(&self) -> usize {
        FakeProducer(Arc::clone(&self.pipeline)).limit()
    }

    fn active_slots(&self) -> usize {
        FakeProducer(Arc::clone(&self.pipeline)).active_slots()
    }

    fn pressure(&self) -> ProducerPressure {
        FakeProducer(Arc::clone(&self.pipeline)).pressure()
    }

    fn work(&self) -> Work {
        FakeProducer(Arc::clone(&self.pipeline)).work()
    }
}

impl Producer for BufferedProducer {
    fn set_limit(&self, n: usize) -> Result<(), thread_broker::ResizeError> {
        FakeProducer(Arc::clone(&self.pipeline)).set_limit(n)
    }

    fn limit(&self) -> usize {
        FakeProducer(Arc::clone(&self.pipeline)).limit()
    }

    fn active_slots(&self) -> usize {
        FakeProducer(Arc::clone(&self.pipeline)).active_slots()
    }

    fn buffered_items(&self) -> Option<u64> {
        Some(self.buffered.load(Ordering::Relaxed))
    }

    fn pressure(&self) -> ProducerPressure {
        FakeProducer(Arc::clone(&self.pipeline)).pressure()
    }

    fn work(&self) -> Work {
        FakeProducer(Arc::clone(&self.pipeline)).work()
    }
}

impl Consumer for ResettingConsumer {
    fn set_threads(&self, n: usize) -> Result<(), thread_broker::ResizeError> {
        FakeConsumer(Arc::clone(&self.pipeline)).set_threads(n)
    }

    fn live_threads(&self) -> usize {
        FakeConsumer(Arc::clone(&self.pipeline)).live_threads()
    }

    fn work(&self) -> Work {
        let work = FakeConsumer(Arc::clone(&self.pipeline)).work();
        if self.calls.fetch_add(1, Ordering::Relaxed) == 3 {
            Work::default()
        } else {
            work
        }
    }
}

impl Producer for ResettingProducer {
    fn set_limit(&self, n: usize) -> Result<(), thread_broker::ResizeError> {
        FakeProducer(Arc::clone(&self.pipeline)).set_limit(n)
    }

    fn limit(&self) -> usize {
        FakeProducer(Arc::clone(&self.pipeline)).limit()
    }

    fn active_slots(&self) -> usize {
        FakeProducer(Arc::clone(&self.pipeline)).active_slots()
    }

    fn work(&self) -> Work {
        let work = FakeProducer(Arc::clone(&self.pipeline)).work();
        if self.calls.fetch_add(1, Ordering::Relaxed) == 3 {
            Work::default()
        } else {
            work
        }
    }

    fn pressure(&self) -> ProducerPressure {
        FakeProducer(Arc::clone(&self.pipeline)).pressure()
    }
}

impl Consumer for FaultyConsumer {
    fn set_threads(&self, n: usize) -> Result<(), thread_broker::ResizeError> {
        let call = self.calls.fetch_add(1, Ordering::Relaxed) + 1;
        if (self.refuse)(n, call) {
            return Err(thread_broker::ResizeError::new(format!(
                "injected consumer refusal on call {call}"
            )));
        }
        FakeConsumer(Arc::clone(&self.pipeline)).set_threads(n)
    }

    fn live_threads(&self) -> usize {
        FakeConsumer(Arc::clone(&self.pipeline)).live_threads()
    }

    fn work(&self) -> Work {
        FakeConsumer(Arc::clone(&self.pipeline)).work()
    }
}

impl Producer for FaultyProducer {
    fn set_limit(&self, n: usize) -> Result<(), thread_broker::ResizeError> {
        let call = self.calls.fetch_add(1, Ordering::Relaxed) + 1;
        if (self.refuse)(n, call) {
            return Err(thread_broker::ResizeError::new(format!(
                "injected producer refusal on call {call}"
            )));
        }
        FakeProducer(Arc::clone(&self.pipeline)).set_limit(n)
    }

    fn limit(&self) -> usize {
        FakeProducer(Arc::clone(&self.pipeline)).limit()
    }

    fn active_slots(&self) -> usize {
        FakeProducer(Arc::clone(&self.pipeline)).active_slots()
    }

    fn work(&self) -> Work {
        FakeProducer(Arc::clone(&self.pipeline)).work()
    }

    fn pressure(&self) -> ProducerPressure {
        FakeProducer(Arc::clone(&self.pipeline)).pressure()
    }
}

fn never_refuse(_: usize, _: usize) -> bool {
    false
}

fn run_faulty(
    pipeline: Arc<Pipeline>,
    consumer_rule: RefusalRule,
    producer_rule: RefusalRule,
) -> Result<BrokerReport, BrokerError> {
    let broker = ThreadBroker::builder_with(
        FaultyConsumer::new(Arc::clone(&pipeline), consumer_rule),
        FaultyProducer::new(Arc::clone(&pipeline), producer_rule),
        BrokerConfig {
            resize_timeout: Duration::from_millis(100),
            ..quick()
        },
    )
    .budget(32)
    .initial_producer_slots(2)
    .build()
    .expect("valid config")
    .start()?;
    std::thread::sleep(Duration::from_millis(800));
    broker.finish()
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
    .start()
    .expect("broker starts");

    std::thread::sleep(dur / 2);
    during(&pipeline);
    std::thread::sleep(dur / 2);
    broker.finish().expect("broker finishes")
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

/// The response-curve policy must not turn its opening split into a hard floor.
/// When ordinary consumer scaling makes the model's floor-directed move safe,
/// that measured move is the discriminator and no upward probe is needed. This
/// is the large-budget counterpart to the negative-scaling test below.
#[test]
fn nonlinear_policy_keeps_a_valid_floor_move_instead_of_the_opening() {
    let cfg = BrokerConfig {
        nonlinear_probes: true,
        ..quick()
    };
    let pipeline = Pipeline::new(100.0, 10_000.0, 28, 4);
    assert_eq!(pipeline.optimum(32), 1);

    let report = run_for(Arc::clone(&pipeline), 32, cfg, Duration::from_millis(1800));
    assert_eq!(report.final_producer_limit, 1, "{report:?}");
    assert_eq!(report.nonlinear_probes, 0, "{report:?}");
    assert!(!report.nonlinear_override, "{report:?}");
    assert!(report.producer_trajectory.contains(&1), "{report:?}");
}

/// A one-point cost model cannot see that adding consumers makes this stage
/// slower. Responsive mode must measure the response curve, while freeze is the
/// explicitly cheaper cost-model-only policy.
#[test]
fn responsive_probes_negative_consumer_scaling_but_freeze_skips_it() {
    let cfg = BrokerConfig {
        nonlinear_probes: true,
        ..quick()
    };
    let responsive = Pipeline::new(100.0, 10_000.0, 6, 2);
    responsive.state.lock().unwrap().consumer_scaling = -1.0;
    assert_eq!(responsive.optimum(8), 7);

    let report = run_for(Arc::clone(&responsive), 8, cfg, Duration::from_millis(1800));
    assert_eq!(report.final_producer_limit, 7, "{report:?}");
    assert!(report.nonlinear_override);
    assert_eq!(report.nonlinear_probes, 3);
    assert_eq!(report.nonlinear_probe_improvements, 3);
    assert!(
        report.producer_trajectory.contains(&1),
        "floor-directed validation was skipped: {report:?}"
    );
    assert!(
        report.reverts >= 1,
        "regressed floor was not restored: {report:?}"
    );

    let frozen = Pipeline::new(100.0, 10_000.0, 6, 2);
    frozen.state.lock().unwrap().consumer_scaling = -1.0;
    let broker = ThreadBroker::builder_with(
        FakeConsumer(Arc::clone(&frozen)),
        FakeProducer(Arc::clone(&frozen)),
        cfg,
    )
    .budget(8)
    .initial_producer_slots(2)
    .steady_state_policy(SteadyStatePolicy::FreezeAfterConvergence)
    .build()
    .unwrap()
    .start()
    .unwrap();
    let report = broker.finish().unwrap();
    assert_eq!(report.nonlinear_probes, 0);
    assert!(!report.nonlinear_override);
    assert_ne!(report.final_producer_limit, 7);
}

/// A geometric endpoint can bracket the discrete optimum without landing on
/// it. Producer 6 is no better than the safe opening at 4, producer 7 is worse,
/// and producer 5 is the true peak; completing the curve at 6 would therefore
/// miss the only allocation that matters.
#[test]
fn nonlinear_search_refines_the_unmeasured_interior_peak() {
    let cfg = BrokerConfig {
        nonlinear_probes: true,
        ..quick()
    };
    let pipeline = Pipeline::new(2_887.0, 1_000.0, 4, 4);
    {
        let mut state = pipeline.state.lock().unwrap();
        state.consumer_scaling = 0.5;
        // Force the isolated busy-cost model to its floor so the response
        // fallback, rather than that deliberately biased model, owns the test.
        state.producer_bias = 0.01;
    }
    assert_eq!(pipeline.optimum(8), 5);

    let report = run_for(Arc::clone(&pipeline), 8, cfg, Duration::from_millis(2500));
    assert_eq!(report.final_producer_limit, 5, "{report:?}");
    assert!(report.nonlinear_override, "{report:?}");
    assert!(report.nonlinear_probes <= 4, "{report:?}");
    assert!(report.moves <= 5, "{report:?}");
    assert!(
        report.producer_trajectory.contains(&5),
        "interior producer split was never measured: {report:?}"
    );
    let transitions: Vec<_> = report
        .producer_trajectory
        .windows(2)
        .filter_map(|pair| (pair[0] != pair[1]).then_some((pair[0], pair[1])))
        .collect();
    assert!(
        transitions.contains(&(4, 5)),
        "interior probe did not restore its retained adjacent baseline: {report:?}"
    );
}

/// Full-calibration freeze must retain the response-search result and release
/// its recurring adapter only after the discrete interior optimum was tested.
#[test]
fn full_calibration_freeze_refines_then_releases_the_producer() {
    let cfg = BrokerConfig {
        nonlinear_probes: true,
        ..quick()
    };
    let pipeline = Pipeline::new(2_887.0, 1_000.0, 4, 4);
    {
        let mut state = pipeline.state.lock().unwrap();
        state.consumer_scaling = 0.5;
        state.producer_bias = 0.01;
    }
    assert_eq!(pipeline.optimum(8), 5);

    let dropped = Arc::new(AtomicBool::new(false));
    let broker = ThreadBroker::builder_with(
        FakeConsumer(Arc::clone(&pipeline)),
        DropAwareProducer {
            pipeline,
            dropped: Arc::clone(&dropped),
        },
        cfg,
    )
    .budget(8)
    .initial_producer_slots(4)
    .steady_state_policy(SteadyStatePolicy::FreezeAfterFullCalibration)
    .build()
    .unwrap()
    .start()
    .unwrap();

    let deadline = Instant::now() + Duration::from_secs(3);
    while !dropped.load(Ordering::Acquire) && Instant::now() < deadline {
        std::thread::sleep(Duration::from_millis(10));
    }
    assert!(
        dropped.load(Ordering::Acquire),
        "full-calibration freeze kept its measurement adapter alive"
    );

    let report = broker.finish().unwrap();
    assert_eq!(
        report.steady_state_policy,
        SteadyStatePolicy::FreezeAfterFullCalibration
    );
    assert!(report.monitoring_stopped_after_convergence, "{report:?}");
    assert!(report.converged(), "{report:?}");
    assert_eq!(report.final_producer_limit, 5, "{report:?}");
    assert!(report.nonlinear_override, "{report:?}");
    assert!(report.nonlinear_probes <= 4, "{report:?}");
    assert!(report.moves <= 5, "{report:?}");
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

/// The cost-share fixed point remains valid when thread-time per item changes
/// with allocation: the observed occupancy term is exactly what moves the next
/// solution toward the nonlinear supply/demand balance. This deliberately
/// violates the old fake's perfect linear-scaling assumption on both sides.
#[test]
fn reaches_the_oracle_with_sublinear_stage_scaling() {
    let cfg = quick();
    let pipeline = Pipeline::new(2_000.0, 4_500.0, 60, 4);
    {
        let mut state = pipeline.state.lock().unwrap();
        state.producer_scaling = 0.62;
        state.consumer_scaling = 0.78;
    }
    let oracle = pipeline.optimum(64);
    let report = run_for(Arc::clone(&pipeline), 64, cfg, Duration::from_millis(2500));
    let ratio = pipeline.throughput_ratio_to_oracle(64, report.final_producer_limit);
    assert!(
        ratio >= 0.90,
        "sublinear auto split {} delivered {:.1}% of oracle split {oracle}: {report:?}",
        report.final_producer_limit,
        ratio * 100.0,
    );
    assert!(
        report.moves <= 5,
        "sublinear response did not settle within five moves: {report:?}"
    );
}

#[test]
fn starvation_vetoes_a_shrink_that_an_allocation_dependent_model_underprices() {
    let cfg = quick();
    let pipeline = Pipeline::new(1_100.0, 3_000.0, 6, 2);
    pipeline.state.lock().unwrap().producer_bias = 0.1;
    assert_eq!(pipeline.optimum(8), 2);
    let report = run_for(Arc::clone(&pipeline), 8, cfg, Duration::from_millis(1200));
    assert_eq!(report.final_producer_limit, 2, "{report:?}");
    assert!(report.pressure_vetoed_shrinks > 0, "{report:?}");
}

/// A producer that cannot use more slots must not be given them, however
/// expensive per item it looks. This is FLINK-31215 in miniature.
#[test]
fn does_not_buy_slots_the_source_cannot_feed() {
    let cfg = quick();
    let pipeline = Pipeline::new(4_000.0, 1_000.0, 60, 4);
    // The unclamped model wants 51 of 64. The source can only sustain what 8
    // decode threads would produce, so everything above that is waste.
    assert_eq!(pipeline.optimum(64), 51);
    let cap = 8.0 * NANOS / 4_000.0;
    pipeline.state.lock().unwrap().source_cap = Some(cap);

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
    let pipeline = Pipeline::new(200.0, 8_000.0, 30, 2);
    pipeline.state.lock().unwrap().producer_bias = 20.0;

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
    #[cfg(target_os = "linux")]
    {
        assert!(report.controller_cpu_nanos.is_some());
        assert_eq!(report.controller_cpu_accounting_failures, 0);
    }
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
        // Halfway through, decode becomes the expensive side. Costs are part
        // of the same locked simulation state as progress, so this is a real
        // regime change without racing the broker thread.
        p.state.lock().unwrap().s_c = 300.0;
    });
    assert!(
        report.resurveys > 0 || report.final_producer_limit > 8,
        "ignored the regime change: settled at {} with {} resurveys",
        report.final_producer_limit,
        report.resurveys,
    );
}

/// A uniform slowdown changes throughput but not the correct split. Reopening
/// here can only spend work and was observed on stable finite-input tails.
#[test]
fn common_mode_rate_changes_do_not_open_a_new_epoch() {
    let cfg = quick();
    let p = Pipeline::new(1_000.0, 3_000.0, 24, 8);
    let report = run_with(
        Arc::clone(&p),
        32,
        cfg,
        Duration::from_millis(3000),
        |pipeline| {
            let mut state = pipeline.state.lock().unwrap();
            state.s_p *= 2.0;
            state.s_c *= 2.0;
        },
    );
    assert_eq!(
        report.resurveys, 0,
        "common-mode slowdown reopened: {report:?}"
    );
    assert_near(report.final_producer_limit, 8, &cfg, "common-mode slowdown");
}

#[test]
fn a_rejected_target_is_reconsidered_only_after_a_real_epoch_change() {
    let cfg = quick();
    let p = Pipeline::new(200.0, 8_000.0, 30, 2);
    p.state.lock().unwrap().producer_bias = 20.0;
    let report = run_with(
        Arc::clone(&p),
        32,
        cfg,
        Duration::from_millis(4000),
        |pipeline| {
            let mut state = pipeline.state.lock().unwrap();
            // The biased first epoch proposes 11 and regresses. In the second
            // epoch, 11 really is the analytic optimum.
            state.s_p = 8_000.0 * 11.0 / 21.0;
            state.producer_bias = 1.0;
        },
    );
    assert!(
        report
            .rejections
            .iter()
            .any(|rejection| rejection.epoch == 0 && rejection.producer_target == 11),
        "first epoch did not reject the biased target: {report:?}",
    );
    assert!(
        report.final_epoch > 0,
        "regime change did not open a new epoch"
    );
    assert_near(report.final_producer_limit, 11, &cfg, "revisited target");
    assert!(
        report.moves <= 6,
        "too many moves across two epochs: {report:?}"
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
        assert!(c + d <= 32, "split oversubscribed with {} slots", c + d);
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
fn freeze_policy_releases_the_producer_after_honest_convergence() {
    let pipeline = Pipeline::new(1_000.0, 1_000.0, 16, 16);
    let dropped = Arc::new(AtomicBool::new(false));
    let broker = ThreadBroker::builder_with(
        FakeConsumer(Arc::clone(&pipeline)),
        DropAwareProducer {
            pipeline,
            dropped: Arc::clone(&dropped),
        },
        quick(),
    )
    .budget(32)
    .initial_producer_slots(16)
    .steady_state_policy(SteadyStatePolicy::FreezeAfterConvergence)
    .build()
    .unwrap()
    .start()
    .unwrap();

    let deadline = Instant::now() + Duration::from_secs(2);
    while !dropped.load(Ordering::Acquire) && Instant::now() < deadline {
        std::thread::sleep(Duration::from_millis(10));
    }
    assert!(
        dropped.load(Ordering::Acquire),
        "freeze policy kept the producer measurement adapter alive"
    );

    let report = broker.finish().unwrap();
    assert_eq!(
        report.steady_state_policy,
        SteadyStatePolicy::FreezeAfterConvergence
    );
    assert!(report.monitoring_stopped_after_convergence);
    assert!(report.converged());
    assert_eq!(report.final_producer_limit, 16);
}

#[test]
fn responsive_policy_remains_live_until_finish() {
    let pipeline = Pipeline::new(1_000.0, 1_000.0, 16, 16);
    let dropped = Arc::new(AtomicBool::new(false));
    let broker = ThreadBroker::builder_with(
        FakeConsumer(Arc::clone(&pipeline)),
        DropAwareProducer {
            pipeline,
            dropped: Arc::clone(&dropped),
        },
        quick(),
    )
    .budget(32)
    .initial_producer_slots(16)
    .steady_state_policy(SteadyStatePolicy::Responsive)
    .build()
    .unwrap()
    .start()
    .unwrap();

    std::thread::sleep(Duration::from_millis(300));
    assert!(
        !dropped.load(Ordering::Acquire),
        "responsive policy terminated before the application finished"
    );
    let report = broker.finish().unwrap();
    assert!(dropped.load(Ordering::Acquire));
    assert_eq!(report.steady_state_policy, SteadyStatePolicy::Responsive);
    assert!(!report.monitoring_stopped_after_convergence);
}

#[test]
fn sparse_responsive_probe_is_interruptible_at_finish() {
    let pipeline = Pipeline::new(1_000.0, 1_000.0, 16, 16);
    let broker = ThreadBroker::builder_with(
        FakeConsumer(Arc::clone(&pipeline)),
        FakeProducer(pipeline),
        quick(),
    )
    .budget(32)
    .initial_producer_slots(16)
    .steady_probe_interval(Duration::from_secs(30))
    .build()
    .unwrap()
    .start()
    .unwrap();

    // Warm-up plus the first ordinary-cadence clean steady window completes
    // well before this. The controller should then be waiting on its 30-second
    // responsive probe, but finish must wake it rather than joining for 30 s.
    std::thread::sleep(Duration::from_millis(350));
    let finishing = Instant::now();
    let report = broker.finish().unwrap();
    assert!(
        finishing.elapsed() < Duration::from_millis(250),
        "finish waited for the sparse probe timeout"
    );
    assert!(report.converged());
    assert_eq!(report.steady_probe_interval, Duration::from_secs(30));
    assert!(report.controller_samples < 10, "{report:?}");
}

#[test]
fn rejects_a_zero_steady_probe_interval() {
    let pipeline = Pipeline::new(1_000.0, 1_000.0, 16, 16);
    assert!(
        ThreadBroker::builder(FakeConsumer(Arc::clone(&pipeline)), FakeProducer(pipeline))
            .budget(32)
            .steady_probe_interval(Duration::ZERO)
            .build()
            .is_err()
    );
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

fn assert_resize_failure(
    error: BrokerError,
    expected_side: ResizeSide,
    expected_request: impl Fn(usize) -> bool,
) {
    match error.kind() {
        BrokerErrorKind::ResizeRefused {
            side, requested, ..
        } => {
            assert_eq!(*side, expected_side);
            assert!(
                expected_request(*requested),
                "unexpected target {requested}"
            );
        }
        other => panic!("expected resize refusal, got {other:?}"),
    }
    let report = error.report().expect("runtime failure carries a report");
    assert!(
        report.peak_controlled_slots <= 32,
        "failure path oversubscribed: {report:?}"
    );
}

#[test]
fn surfaces_a_refused_consumer_shrink_without_growing_the_producer() {
    let p = Pipeline::new(1_000.0, 1_000.0, 30, 2);
    let error = run_faulty(
        Arc::clone(&p),
        |target, call| call > 1 && target < 30,
        never_refuse,
    )
    .expect_err("consumer shrink should fail");
    assert_resize_failure(error, ResizeSide::Consumer, |target| target < 30);
    assert_eq!(p.split(), (30, 2));
}

#[test]
fn surfaces_a_refused_producer_growth_after_a_safe_consumer_shrink() {
    let p = Pipeline::new(1_000.0, 1_000.0, 30, 2);
    let error = run_faulty(Arc::clone(&p), never_refuse, |target, _| target > 2)
        .expect_err("producer growth should fail");
    assert_resize_failure(error, ResizeSide::Producer, |target| target > 2);
    let (consumer, producer) = p.split();
    assert_eq!(producer, 2, "refused growth changed producer occupancy");
    assert!(consumer + producer <= 32);
}

#[test]
fn surfaces_a_refused_producer_shrink_during_rollback() {
    let p = Pipeline::new(200.0, 8_000.0, 30, 2);
    p.state.lock().unwrap().producer_bias = 20.0;
    let error = run_faulty(Arc::clone(&p), never_refuse, |target, call| {
        call > 1 && target == 2
    })
    .expect_err("rollback producer shrink should fail");
    assert_resize_failure(error, ResizeSide::Producer, |target| target == 2);
}

#[test]
fn surfaces_a_refused_consumer_growth_during_rollback() {
    let p = Pipeline::new(200.0, 8_000.0, 30, 2);
    p.state.lock().unwrap().producer_bias = 20.0;
    let error = run_faulty(
        Arc::clone(&p),
        |target, call| call > 1 && target == 30,
        never_refuse,
    )
    .expect_err("rollback consumer growth should fail");
    assert_resize_failure(error, ResizeSide::Consumer, |target| target == 30);
}

struct StuckConsumer(Arc<Pipeline>);

impl Consumer for StuckConsumer {
    fn set_threads(&self, n: usize) -> Result<(), thread_broker::ResizeError> {
        let current = self.0.state.lock().unwrap().consumer_threads;
        if n >= current {
            FakeConsumer(Arc::clone(&self.0)).set_threads(n)?;
        }
        Ok(())
    }

    fn live_threads(&self) -> usize {
        FakeConsumer(Arc::clone(&self.0)).live_threads()
    }

    fn work(&self) -> Work {
        FakeConsumer(Arc::clone(&self.0)).work()
    }
}

#[test]
fn a_shrink_that_never_acknowledges_times_out_before_growth() {
    let p = Pipeline::new(1_000.0, 1_000.0, 30, 2);
    let cfg = BrokerConfig {
        resize_timeout: Duration::from_millis(60),
        ..quick()
    };
    let broker = ThreadBroker::builder_with(
        StuckConsumer(Arc::clone(&p)),
        FakeProducer(Arc::clone(&p)),
        cfg,
    )
    .budget(32)
    .initial_producer_slots(2)
    .build()
    .unwrap()
    .start()
    .unwrap();
    std::thread::sleep(Duration::from_millis(500));
    let error = broker.finish().expect_err("stuck shrink must time out");
    assert!(matches!(
        error.kind(),
        BrokerErrorKind::ResizeTimedOut {
            side: ResizeSide::Consumer,
            ..
        }
    ));
    assert_eq!(p.split(), (30, 2), "producer grew before shrink ack");
    assert!(error.report().unwrap().peak_controlled_slots <= 32);
}

struct PanickingConsumer(Arc<Pipeline>);

impl Consumer for PanickingConsumer {
    fn set_threads(&self, n: usize) -> Result<(), thread_broker::ResizeError> {
        FakeConsumer(Arc::clone(&self.0)).set_threads(n)
    }

    fn live_threads(&self) -> usize {
        FakeConsumer(Arc::clone(&self.0)).live_threads()
    }

    fn work(&self) -> Work {
        panic!("injected sampler panic")
    }
}

#[test]
fn sampler_panics_are_returned_instead_of_default_reports() {
    let p = Pipeline::new(1_000.0, 1_000.0, 30, 2);
    let broker =
        ThreadBroker::builder_with(PanickingConsumer(Arc::clone(&p)), FakeProducer(p), quick())
            .budget(32)
            .initial_producer_slots(2)
            .build()
            .unwrap()
            .start()
            .unwrap();
    let error = broker.finish().expect_err("panic must surface");
    assert!(matches!(error.kind(), BrokerErrorKind::ThreadPanicked));
}

#[test]
fn decreasing_consumer_counters_are_a_structured_failure() {
    let p = Pipeline::new(1_000.0, 1_000.0, 30, 2);
    let broker = ThreadBroker::builder_with(
        ResettingConsumer {
            pipeline: Arc::clone(&p),
            calls: AtomicUsize::new(0),
        },
        FakeProducer(p),
        quick(),
    )
    .budget(32)
    .initial_producer_slots(2)
    .build()
    .unwrap()
    .start()
    .unwrap();
    std::thread::sleep(Duration::from_millis(150));
    let error = broker.finish().expect_err("counter reset must surface");
    assert!(matches!(
        error.kind(),
        BrokerErrorKind::WorkCountersRegressed {
            side: ResizeSide::Consumer,
            ..
        }
    ));
    assert!(error.report().is_some());
}

#[test]
fn decreasing_producer_counters_are_a_structured_failure() {
    let p = Pipeline::new(1_000.0, 1_000.0, 30, 2);
    let broker = ThreadBroker::builder_with(
        FakeConsumer(Arc::clone(&p)),
        ResettingProducer {
            pipeline: p,
            calls: AtomicUsize::new(0),
        },
        quick(),
    )
    .budget(32)
    .initial_producer_slots(2)
    .build()
    .unwrap()
    .start()
    .unwrap();
    std::thread::sleep(Duration::from_millis(150));
    let error = broker.finish().expect_err("counter reset must surface");
    assert!(matches!(
        error.kind(),
        BrokerErrorKind::WorkCountersRegressed {
            side: ResizeSide::Producer,
            ..
        }
    ));
    assert!(error.report().is_some());
}

#[test]
fn initial_split_respects_a_multi_thread_consumer_floor() {
    let p = Pipeline::new(1_000.0, 1_000.0, 1, 31);
    let cfg = BrokerConfig {
        min_consumer_threads: 8,
        ..quick()
    };
    let broker = ThreadBroker::builder_with(
        FakeConsumer(Arc::clone(&p)),
        FakeProducer(Arc::clone(&p)),
        cfg,
    )
    .budget(32)
    .initial_producer_slots(31)
    .build()
    .unwrap()
    .start()
    .unwrap();
    std::thread::sleep(Duration::from_millis(10));
    let report = broker.finish().unwrap();
    assert_eq!(p.split(), (8, 24));
    assert_eq!(
        (report.final_consumer_threads, report.final_producer_limit),
        (8, 24)
    );
}

#[test]
fn changing_buffer_occupancy_is_rejected_as_flow_transient_evidence() {
    let p = Pipeline::new(3_000.0, 1_000.0, 30, 2);
    let buffered = Arc::new(AtomicU64::new(0));
    let stop = Arc::new(std::sync::atomic::AtomicBool::new(false));
    let updater = {
        let buffered = Arc::clone(&buffered);
        let stop = Arc::clone(&stop);
        std::thread::spawn(move || {
            while !stop.load(Ordering::Acquire) {
                buffered.fetch_add(1, Ordering::Relaxed);
                std::thread::sleep(Duration::from_millis(2));
            }
        })
    };
    let cfg = BrokerConfig {
        max_buffer_drift_items: 0,
        max_buffer_drift_fraction: 0.0,
        ..quick()
    };
    let broker = ThreadBroker::builder_with(
        FakeConsumer(Arc::clone(&p)),
        BufferedProducer {
            pipeline: p,
            buffered,
        },
        cfg,
    )
    .budget(32)
    .initial_producer_slots(2)
    .build()
    .unwrap()
    .start()
    .unwrap();
    std::thread::sleep(Duration::from_millis(400));
    let report = broker.finish().unwrap();
    stop.store(true, Ordering::Release);
    updater.join().unwrap();
    assert!(report.flow_transient_windows > 0);
    assert_eq!(report.moves, 0, "acted on mismatched-flow windows");
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

#[test]
fn completed_item_does_not_publish_a_started_record_early() {
    let meter = thread_broker::BusyMeter::default();
    let mut timer = meter.timer().flush_every(1);
    {
        let _item = timer.completed_item();
        assert_eq!(meter.work().items, 0);
    }
    assert_eq!(meter.work().items, 1);
}

#[test]
fn fine_progress_does_not_force_fine_busy_time_clock_reads() {
    let meter = thread_broker::BusyMeter::default();
    let mut timer = meter.timer().progress_every(1).time_every(4);

    std::thread::sleep(Duration::from_millis(1));
    timer.tick();
    let progress_only = meter.work();
    assert_eq!(progress_only.items, 1);
    assert_eq!(progress_only.busy_nanos, 0);

    for _ in 0..3 {
        timer.tick();
    }
    let timed = meter.work();
    assert_eq!(timed.items, 4);
    assert!(timed.busy_nanos > 0);
}

#[test]
fn progress_can_publish_to_a_worker_local_counter() {
    let meter = thread_broker::BusyMeter::default();
    let local = AtomicU64::new(0);
    for expected in [4, 7] {
        let mut timer = meter
            .timer()
            .single_writer_progress_counter(&local)
            .progress_every(1)
            .time_every(4);
        while local.load(Ordering::Relaxed) < expected {
            timer.tick();
        }
        drop(timer);
        assert_eq!(local.load(Ordering::Relaxed), expected);
    }
    assert_eq!(meter.work().items, 0);
    assert!(meter.work().busy_nanos > 0);
}

#[test]
fn max_cadences_publish_only_when_the_batch_ends() {
    let meter = thread_broker::BusyMeter::default();
    let mut timer = meter.timer().progress_every(u64::MAX).time_every(u64::MAX);
    timer.tick();
    assert_eq!(meter.work(), thread_broker::Work::default());

    drop(timer);
    let done = meter.work();
    assert_eq!(done.items, 1);
    assert!(done.busy_nanos > 0);
}
