//! The control loop: measure each stage's true rate, solve for the split, jump.
//!
//! The shape of this loop is taken from DS2 (Kalavri et al., OSDI'18): estimate
//! a rate per stage that is independent of the resources currently assigned to
//! it, then compute the assignment directly rather than searching for it. The
//! guards around the jump — ratify against throughput, remember what failed,
//! blank out measurements taken while the system is still moving — come from
//! Flink's autoscaler (FLIP-271) and from FDP (Suleman et al., PACT'10), which
//! demote hill-climbing to a check on a model-based decision rather than the
//! decision itself.
//!
//! See `DESIGN-thread-broker.md` for the derivation and the citations.

use std::collections::HashSet;
use std::sync::Arc;
use std::thread::JoinHandle;
use std::time::{Duration, Instant};

use crate::{Consumer, Producer, ProducerPressure, Stop, ThreadBroker, Work, idle_fraction, now};

/// What the broker did, for logging and for tests.
#[derive(Debug, Clone, Default)]
pub struct BrokerReport {
    /// Consumer threads at each decision, in order.
    pub consumer_trajectory: Vec<usize>,
    /// Producer limit at each decision, in order.
    pub producer_trajectory: Vec<usize>,
    /// Splits applied, including reverts.
    pub moves: usize,
    /// Moves undone because throughput got worse.
    pub reverts: usize,
    /// Times a settled split was re-opened by a detected regime change.
    pub resurveys: usize,
    /// Windows discarded as carrying too little work to estimate a rate from.
    pub rejected_windows: usize,
    /// Samples where the run was source-bound rather than decode-bound.
    pub source_bound_samples: usize,
    /// Samples where the producer could not use more concurrency at all.
    pub inelastic_samples: usize,
    /// Time from the end of warm-up until the broker settled, if it did.
    pub time_to_converge: Option<Duration>,
    /// Final split.
    pub final_consumer_threads: usize,
    pub final_producer_limit: usize,
    /// Last solved model, for diagnosing a split that looks wrong.
    pub final_model: Option<Model>,
}

impl BrokerReport {
    /// Whether the broker settled and stayed settled.
    pub fn converged(&self) -> bool {
        self.time_to_converge.is_some()
    }
}

/// The solved model at one instant.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Model {
    /// Fraction of the pipeline's total per-item cost that falls on the
    /// producer: `s_p / (s_p + s_c)`. The split is this times the budget.
    pub producer_cost_share: f64,
    /// Producer slots the model asks for, before clamping.
    pub ideal_producer_slots: usize,
    /// Producer slots above which more concurrency has been observed not to
    /// help. `usize::MAX` when nothing has been observed to bound it.
    pub useful_cap: usize,
}

/// A broker sampling in the background. Dropping it stops the loop.
pub struct RunningBroker {
    stop: Arc<Stop>,
    join: Option<JoinHandle<BrokerReport>>,
}

impl RunningBroker {
    /// Stop sampling and collect the report.
    pub fn finish(mut self) -> BrokerReport {
        self.stop.set();
        self.join
            .take()
            .map(|h| h.join().unwrap_or_default())
            .unwrap_or_default()
    }
}

impl Drop for RunningBroker {
    fn drop(&mut self) {
        self.stop.set();
        if let Some(h) = self.join.take() {
            let _ = h.join();
        }
    }
}

/// Consecutive samples of a producer that provably cannot use more slots before
/// its current limit is taken as a ceiling.
///
/// More than one, because both saturating pressures are read from instantaneous
/// pool fields: a single sample can catch an ordinary lull between tasks.
const SATURATION_SAMPLES: usize = 3;

/// A split: consumer threads, producer slots.
type Split = (usize, usize);

/// Where the loop is.
#[derive(Debug, Clone)]
enum Phase {
    /// Gathering cost evidence at the current split.
    Survey,
    /// A move was just applied; measurements are contaminated until it settles.
    Blackout { left: usize, next: Box<Phase> },
    /// The move has settled. Is throughput actually better?
    Ratify {
        left: usize,
        from: Split,
        target: usize,
        baseline: f64,
        items: u64,
        elapsed: Duration,
    },
    /// Settled. Sleeping until the model's answer moves a long way, and stays
    /// moved.
    Steady { drifted_for: usize },
}

impl<C: Consumer + 'static, P: Producer + 'static> ThreadBroker<C, P> {
    /// Apply the starting split and begin sampling.
    ///
    /// The start is deliberately coarse. Under a solving controller it is only
    /// the point at which the first rates are measured, not a position that has
    /// to be walked away from one step at a time — which is what made the
    /// starting point dominate the outcome under the previous control law. It
    /// stays biased toward the consumer because a run that never needs decode
    /// concurrency then never pays for any.
    pub fn start(self) -> RunningBroker {
        let floor = self.config.min_producer_slots.max(1);
        let coarse = (self.budget / 8).clamp(floor, self.budget - 1);
        let initial_producer = self
            .initial_producer_slots
            .unwrap_or(coarse)
            .clamp(floor.min(self.budget - 1), self.budget - 1);
        let initial_consumer = self
            .budget
            .saturating_sub(initial_producer)
            .max(self.config.min_consumer_threads);

        self.consumer.set_threads(initial_consumer);
        self.producer.set_limit(initial_producer);
        let initial = (initial_consumer, initial_producer);

        let stop = Stop::new();
        let stop_thread = Arc::clone(&stop);
        let join = std::thread::Builder::new()
            .name("thread-broker".into())
            .spawn(move || self.run(stop_thread, initial))
            .ok();

        RunningBroker { stop, join }
    }

    fn run(self, stop: Arc<Stop>, initial: Split) -> BrokerReport {
        let cfg = self.config;
        let mut report = BrokerReport::default();

        let started = now();
        let mut prev_consumer = self.consumer.work();
        let mut prev_producer = self.producer.work();
        let mut prev_at = started;

        // The broker's own view of the split, rather than reading it back. Both
        // pools converge toward what they are told over some interval, so
        // reading the current value back conflates "what I asked for" with
        // "what has taken effect yet". The budget must be enforced against the
        // former or two consecutive decisions can double-spend the same threads.
        let mut split = initial;

        let mut costs = Costs::new(cfg.smoothing_windows);
        // Throughput over the same span the ratify stage will measure, so the
        // before/after comparison is like for like. A single pre-move window as
        // the baseline made the test a coin flip: one unlucky window sets a bar
        // the move cannot clear, or an unusually good one hides a genuine
        // regression.
        let mut recent = Recent::new(cfg.ratify_samples);
        let mut caps = Caps::default();
        // Targets that were tried and made throughput worse. FDP calls this the
        // tabu set; without it a model that keeps recomputing the same answer
        // re-applies a move that has already been measured as a regression.
        let mut rejected: HashSet<usize> = HashSet::new();
        let mut phase = Phase::Survey;
        let mut warm_ended: Option<Instant> = None;
        let mut surveyed = 0usize;

        while !stop.is_set() {
            std::thread::sleep(cfg.sample_interval);

            let sampled_at = now();
            let consumer = self.consumer.work();
            let producer = self.producer.work();
            let window = sampled_at.saturating_duration_since(prev_at);
            let dc = consumer.delta(prev_consumer);
            let dp = producer.delta(prev_producer);
            prev_consumer = consumer;
            prev_producer = producer;
            prev_at = sampled_at;

            // Discard the startup ramp entirely, re-basing as we go so the first
            // recorded window contains none of it. Without this every workload
            // reports ~99% consumer starvation in its first 75 ms, because the
            // consumer's threads are still starting and nothing has been
            // processed yet.
            if sampled_at.saturating_duration_since(started) < cfg.warmup {
                continue;
            }
            let warm_ended = *warm_ended.get_or_insert(sampled_at);

            let pressure = self.producer.pressure();
            match pressure {
                ProducerPressure::SourceBound => report.source_bound_samples += 1,
                ProducerPressure::Inelastic => report.inelastic_samples += 1,
                _ => {}
            }

            // Throughput of the pipeline as a whole, which is what the
            // accept/revert decision is actually about. Measured at the
            // consumer, because the producer's output can run ahead into buffers
            // and briefly show progress the pipeline did not make.
            let throughput = rate_per_second(dc.items, window);

            recent.push(dc.items, window);

            let usable = usable_window(dc, dp, window);
            if usable {
                costs.push(dc.busy_nanos, dp.busy_nanos);
                caps.observe(dp, window, split.1, pressure);
            } else {
                report.rejected_windows += 1;
            }

            let live = self.consumer.live_threads();
            let idle = idle_fraction(dc.busy_nanos, live, window);

            phase = match phase {
                Phase::Blackout { left, next } => {
                    if left > 1 {
                        Phase::Blackout {
                            left: left - 1,
                            next,
                        }
                    } else {
                        // Everything measured while the pools were converging is
                        // about a split that no longer exists.
                        costs.clear();
                        *next
                    }
                }

                Phase::Ratify {
                    left,
                    from,
                    target,
                    baseline,
                    items,
                    elapsed,
                } => {
                    let items = items + dc.items;
                    let elapsed = elapsed + window;
                    if left > 1 {
                        Phase::Ratify {
                            left: left - 1,
                            from,
                            target,
                            baseline,
                            items,
                            elapsed,
                        }
                    } else {
                        let achieved = rate_per_second(items, elapsed);
                        // Only a regression reverts; absence of improvement does
                        // not. Near the optimum the surface is locally flat --
                        // measured, 22 and 23 producer slots of 64 differ by 3%,
                        // which is inside this comparison's own noise -- so "no
                        // better" is a perfectly ordinary reading for a *correct*
                        // move. Treating it as failure would revert the broker to
                        // wherever it happened to start, which is the one place
                        // with no evidence for it at all.
                        //
                        // Flink reaches the same conclusion from the other end:
                        // its `detectIneffectiveScaleUp` compares against a
                        // model-*predicted* increase rather than any increase,
                        // and still ships disabled by default.
                        if baseline > 0.0 && achieved < baseline * (1.0 - cfg.regression_tolerance)
                        {
                            tracing::debug!(
                                "thread-broker: reverting producer {} -> {} \
                                 ({:.0} vs {:.0} items/s)",
                                target,
                                from.1,
                                achieved,
                                baseline,
                            );
                            rejected.insert(target);
                            let reverting_from = split;
                            split = from;
                            self.apply(reverting_from, split);
                            report.moves += 1;
                            report.reverts += 1;
                            costs.clear();
                            Phase::Blackout {
                                left: cfg.blackout_samples,
                                next: Box::new(Phase::Steady { drifted_for: 0 }),
                            }
                        } else {
                            // Kept. Re-survey once: the rates may have been
                            // measured under a split that distorted them, and
                            // DS2 reports convergence in at most a few such
                            // rounds. The deadband and the tabu set are what
                            // stop this from cycling.
                            surveyed = 0;
                            Phase::Survey
                        }
                    }
                }

                Phase::Steady { drifted_for } => {
                    // Keep solving, but act only on a large drift. The model is
                    // the change detector: it is computed from the cost ratio
                    // rather than from the outcome, so it moves as soon as the
                    // workload does -- including when the change is an
                    // *opportunity* rather than a loss, which a throughput
                    // monitor cannot see at all.
                    let solved = costs.solve(self.budget, &cfg, &caps);
                    // Refresh the diagnostic even when nothing is acted on: a
                    // report read at the end of a long settled run should
                    // describe the workload as it finished, not as it looked at
                    // the one sample that happened to settle it.
                    if let Some(s) = &solved {
                        report.final_model = Some(s.snapshot);
                    }
                    let drifted = solved
                        .filter(|s| s.target.abs_diff(split.1) >= cfg.resurvey_distance)
                        .filter(|s| !rejected.contains(&s.target));
                    match drifted {
                        // Distance alone is not enough: the noise in the solved
                        // target grows with the budget while a fixed distance
                        // does not, so at 64 threads ordinary jitter clears a
                        // band that 32 threads never approaches. Persistence is
                        // what tells the two apart -- noise does not survive a
                        // second of consecutive windows, a regime change does.
                        Some(solved) if drifted_for + 1 >= cfg.resurvey_samples => {
                            tracing::debug!(
                                "thread-broker: the split the model wants moved {} -> {} \
                                 and stayed moved; re-opening",
                                split.1,
                                solved.target,
                            );
                            report.resurveys += 1;
                            report.time_to_converge = None;
                            surveyed = 0;
                            Phase::Survey
                        }
                        Some(_) => Phase::Steady {
                            drifted_for: drifted_for + 1,
                        },
                        None => {
                            if report.time_to_converge.is_none() {
                                report.time_to_converge =
                                    Some(sampled_at.saturating_duration_since(warm_ended));
                            }
                            Phase::Steady { drifted_for: 0 }
                        }
                    }
                }

                Phase::Survey => {
                    surveyed += usable as usize;
                    match costs.solve(self.budget, &cfg, &caps) {
                        // Not enough evidence yet, or one side has produced no
                        // work at all. Keep looking rather than guessing.
                        _ if surveyed < cfg.smoothing_windows => Phase::Survey,
                        None => Phase::Survey,
                        Some(model) => {
                            let target = model.target;
                            let distance = target.abs_diff(split.1);
                            if distance < cfg.deadband_threads || rejected.contains(&target) {
                                report.final_model = Some(model.snapshot);
                                Phase::Steady { drifted_for: 0 }
                            } else {
                                let from = split;
                                tracing::debug!(
                                    "thread-broker: solved producer {} -> {} \
                                     (cost share {:.2}, cap {})",
                                    from.1,
                                    target,
                                    model.snapshot.producer_cost_share,
                                    fmt_cap(model.snapshot.useful_cap),
                                );
                                split = (self.budget - target, target);
                                self.apply(from, split);
                                report.moves += 1;
                                report.final_model = Some(model.snapshot);
                                Phase::Blackout {
                                    left: cfg.blackout_samples,
                                    next: Box::new(Phase::Ratify {
                                        left: cfg.ratify_samples,
                                        from,
                                        target,
                                        baseline: recent.rate(),
                                        items: 0,
                                        elapsed: Duration::ZERO,
                                    }),
                                }
                            }
                        }
                    }
                }
            };

            report.consumer_trajectory.push(split.0);
            report.producer_trajectory.push(split.1);

            tracing::trace!(
                "thread-broker: {:?} idle {:.1}% pressure {:?} {:.0} items/s \
                 consumer {} (live {live}) producer {}",
                PhaseName(&phase),
                idle * 100.0,
                pressure,
                throughput,
                split.0,
                split.1,
            );
        }

        report.final_consumer_threads = split.0;
        report.final_producer_limit = split.1;
        report
    }

    /// Push a split to both sides, always shrinking before growing.
    ///
    /// The order matters: growing first lets the two sides briefly sum above the
    /// budget, which on a full machine is exactly the oversubscription the budget
    /// exists to prevent.
    ///
    /// The direction comes from the caller rather than from
    /// [`Producer::limit`]. Reading it back would reintroduce the very confusion
    /// the broker avoids by tracking its own view of the split: a pool that has
    /// not yet converged to the last request reports the old value, and the
    /// ordering decision would then be made against a split that no longer
    /// exists.
    fn apply(&self, from: Split, to: Split) {
        let (consumer, producer) = to;
        if producer > from.1 {
            self.consumer.set_threads(consumer);
            self.producer.set_limit(producer);
        } else {
            self.producer.set_limit(producer);
            self.consumer.set_threads(consumer);
        }
    }
}

/// Recent pipeline throughput, over a fixed number of windows.
///
/// Exists only so the ratify comparison averages the same span on both sides of
/// a move. Both halves are `items / elapsed` over `ratify_samples` windows, so
/// neither is more exposed to a single noisy window than the other.
struct Recent {
    ring: std::collections::VecDeque<(u64, Duration)>,
    capacity: usize,
}

impl Recent {
    fn new(windows: usize) -> Self {
        let capacity = windows.max(1);
        Self {
            ring: std::collections::VecDeque::with_capacity(capacity),
            capacity,
        }
    }

    fn push(&mut self, items: u64, window: Duration) {
        if self.ring.len() == self.capacity {
            self.ring.pop_front();
        }
        self.ring.push_back((items, window));
    }

    fn rate(&self) -> f64 {
        let (items, elapsed) = self
            .ring
            .iter()
            .fold((0u64, Duration::ZERO), |(i, d), &(items, w)| {
                (i + items, d + w)
            });
        rate_per_second(items, elapsed)
    }
}

/// Smoothed busy time for the two stages over the recent windows.
///
/// Sums rather than an average of per-window ratios, so a window in which the
/// pipeline barely moved contributes proportionally rather than equally.
struct Costs {
    ring: std::collections::VecDeque<(u64, u64)>,
    capacity: usize,
}

struct Solved {
    target: usize,
    snapshot: Model,
}

impl Costs {
    fn new(windows: usize) -> Self {
        let capacity = windows.max(1);
        Self {
            ring: std::collections::VecDeque::with_capacity(capacity),
            capacity,
        }
    }

    fn push(&mut self, consumer_busy: u64, producer_busy: u64) {
        if self.ring.len() == self.capacity {
            self.ring.pop_front();
        }
        self.ring.push_back((consumer_busy, producer_busy));
    }

    fn clear(&mut self) {
        self.ring.clear();
    }

    /// The closed-form split.
    ///
    /// Writing `s_p` and `s_c` for the two stages' per-item costs, a budget `N`
    /// divided into `d` producer threads and `N - d` consumer threads finishes
    /// both sides together when `d / s_p = (N - d) / s_c`, whose solution is
    ///
    /// ```text
    ///     d* = N * s_p / (s_p + s_c)
    /// ```
    ///
    /// Sanity check the direction: as the producer gets more expensive per item,
    /// `d*` tends to `N`, which is right.
    ///
    /// # Why the busy times can be substituted directly
    ///
    /// Over a window in which the pipeline moved `X` items, the producer spent
    /// `X * s_p` thread-nanoseconds and the consumer `X * s_c`. So
    ///
    /// ```text
    ///     busy_p / (busy_p + busy_c) = s_p / (s_p + s_c)
    /// ```
    ///
    /// exactly — the `X` cancels, and with it any need for the two stages to
    /// count work in the same unit. In particular a run mixing compressed and
    /// uncompressed inputs is handled correctly for free: the producer simply
    /// does no work for the uncompressed ones, its busy time is proportionally
    /// smaller, and the split follows.
    ///
    /// # The assumption in that step
    ///
    /// "The same `X`" is exact only in steady state. Over a window in which the
    /// buffer between the stages is filling or draining, the two sides handled
    /// *different* numbers of items and the estimate is biased by the
    /// difference. Three things bound it — the estimate is summed over
    /// `smoothing_windows` rather than taken from one; a move is followed by a
    /// blackout long enough for the buffer to re-level; and a window in which
    /// the pipeline made no progress is discarded outright by
    /// [`usable_window`] — but none of them is a proof, and a workload whose
    /// compressibility varies in long runs is where to look first if the split
    /// comes out wrong.
    ///
    /// This is solvable at all only because busy time excludes blocking. Wall
    /// time would make both sides a function of `d`, and the equation would say
    /// nothing. DS2 (Kalavri et al., OSDI'18) calls this the *true* rate as
    /// against the *observed* rate, and section 2 of that paper is a catalogue
    /// of the controllers that go wrong by using the latter.
    fn solve(&self, budget: usize, cfg: &crate::BrokerConfig, caps: &Caps) -> Option<Solved> {
        let (consumer, producer) = self
            .ring
            .iter()
            .fold((0u64, 0u64), |(c, p), (dc, dp)| (c + dc, p + dp));
        let total = consumer.checked_add(producer)?;
        if total == 0 {
            return None;
        }

        let share = producer as f64 / total as f64;
        let ideal = (budget as f64 * share).round().max(1.0) as usize;

        let lo = cfg.min_producer_slots.max(1);
        let hi = budget.saturating_sub(cfg.min_consumer_threads).max(lo);
        let useful = caps.useful();
        let target = ideal.clamp(lo, hi).min(useful.max(lo));

        Some(Solved {
            target,
            snapshot: Model {
                producer_cost_share: share,
                ideal_producer_slots: ideal,
                useful_cap: useful,
            },
        })
    }
}

/// The empirical ceiling on producer concurrency.
///
/// The model assumes a stage can use every thread it is given. A producer
/// reading from a saturated disk cannot, and the model has no way to know that:
/// the cost share says how much *work* decoding is, not whether there is
/// anything to decode it from. Left uncapped, the solver buys slots that sit
/// idle — which is FLINK-31215, where a non-scalable bottleneck upstream leads
/// the autoscaler to keep scaling a stage that cannot go any faster.
///
/// Two independent observations bound it, and both are one-sided: a cap is only
/// ever recorded when the producer *demonstrably* failed to use what it had.
///
/// # Why it is a recent window and not a running minimum
///
/// A running minimum never recovers. One quiet stretch — a between-files lull,
/// a stall downstream — would pin the ceiling for the rest of the run, and the
/// broker could not buy decode concurrency again however badly the next file
/// needed it. Bounding it to a recent window makes the cap self-healing, at the
/// cost of forgetting a constraint that has not been observed for a few seconds,
/// which is the right trade: a constraint that still holds will be observed
/// again almost immediately.
struct Caps {
    /// `(concurrency actually achieved, slots granted)` per recent window.
    recent: std::collections::VecDeque<(f64, usize)>,
    saturated_for: usize,
    saturated_at: Option<usize>,
}

/// Windows of history behind the ceiling.
///
/// At the default sampling interval this is a few seconds — long enough to span
/// a decision cycle and its blackout, short enough that a constraint which has
/// genuinely lifted is forgotten quickly.
const CAP_HISTORY: usize = 32;

impl Default for Caps {
    fn default() -> Self {
        Self {
            recent: std::collections::VecDeque::with_capacity(CAP_HISTORY),
            saturated_for: 0,
            saturated_at: None,
        }
    }
}

impl Caps {
    fn observe(&mut self, dp: Work, window: Duration, producer_limit: usize, p: ProducerPressure) {
        // Concurrency the producer actually achieved this window, in
        // thread-equivalents. Busy time excludes blocking, so this is work done
        // rather than slots held.
        let achieved = dp.busy_nanos as f64 / window.as_nanos().max(1) as f64;
        if self.recent.len() == CAP_HISTORY {
            self.recent.pop_front();
        }
        self.recent.push_back((achieved, producer_limit));

        // The only thing `ProducerPressure` is used for. It is a saturating
        // signal -- it cannot tell "just enough" from "far too much" -- so it
        // may veto growth but must never size the split.
        match p {
            ProducerPressure::SourceBound | ProducerPressure::Inelastic => {
                self.saturated_for += 1;
                if self.saturated_for >= SATURATION_SAMPLES {
                    self.saturated_at = Some(match self.saturated_at {
                        Some(prev) => prev.min(producer_limit.max(1)),
                        None => producer_limit.max(1),
                    });
                }
            }
            _ => {
                self.saturated_for = 0;
                self.saturated_at = None;
            }
        }
    }

    /// The largest producer limit worth asking for, or `usize::MAX` if nothing
    /// has been observed to bound it.
    fn useful(&self) -> usize {
        let mut cap = usize::MAX;

        // Slack: if the producer was handed more than a thread's worth of slack
        // and declined to use it, more slots are not what it is short of. The
        // ceiling is the best it managed while it had that slack, plus one so a
        // marginal case is not pinched.
        let slack_seen = self.recent.iter().any(|&(a, limit)| limit as f64 > a + 1.0);
        if slack_seen {
            let peak = self.recent.iter().fold(0.0f64, |m, &(a, _)| m.max(a));
            cap = cap.min(peak.ceil() as usize + 1);
        }

        if let Some(limit) = self.saturated_at {
            cap = cap.min(limit);
        }
        cap
    }
}

/// Whether a window carries enough work to draw a conclusion from.
///
/// The cost share is a ratio, and a window containing a handful of nanoseconds
/// of work on one side gives a ratio dominated by where the window boundaries
/// happened to fall relative to the counters' flush granularity. The CLR thread
/// pool's hill climber guards the same hazard from the other direction,
/// rejecting a sample window when in-flight work is large relative to completed
/// work.
///
/// The test is on the *total*, not on each side, so that a producer doing no
/// work at all is admitted rather than discarded. That case is real and needs
/// the right answer: a run whose inputs are already uncompressed has no decode
/// cost, and the honest reading of the window is a cost share of zero, which
/// solves to the producer floor and hands the whole budget to the consumer.
/// Requiring both sides to be busy would instead reject every window forever and
/// strand the initial guess.
///
/// Progress is required as well as busy time, because a window in which the
/// pipeline moved nothing says only that something is stalled, and which side
/// happened to be spinning during it is not evidence about their relative costs.
///
/// The floor is a twentieth of one thread-window, far below anything a working
/// stage produces and far above the residue of one that has stopped.
fn usable_window(dc: Work, dp: Work, window: Duration) -> bool {
    let floor = window.as_nanos() as u64 / 20;
    dc.items > 0 && dc.busy_nanos.saturating_add(dp.busy_nanos) > floor
}

fn rate_per_second(items: u64, window: Duration) -> f64 {
    let secs = window.as_secs_f64();
    if secs <= 0.0 {
        return 0.0;
    }
    items as f64 / secs
}

fn fmt_cap(cap: usize) -> String {
    if cap == usize::MAX {
        "none".into()
    } else {
        cap.to_string()
    }
}

struct PhaseName<'a>(&'a Phase);

impl std::fmt::Debug for PhaseName<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self.0 {
            Phase::Survey => "survey",
            Phase::Blackout { .. } => "blackout",
            Phase::Ratify { .. } => "ratify",
            Phase::Steady { .. } => "steady",
        })
    }
}
