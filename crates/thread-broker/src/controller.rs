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

use std::collections::HashMap;
use std::sync::Arc;
use std::thread::JoinHandle;
use std::time::{Duration, Instant};

use crate::{
    BrokerError, BrokerErrorKind, Consumer, OpeningPolicy, Producer, ProducerMeasurementMode,
    ProducerMeasurementStats, ProducerPressure, ResizeSide, SteadyStatePolicy, Stop, ThreadBroker,
    Work, idle_fraction, now,
};

/// What the broker did, for logging and for tests.
#[derive(Debug, Clone, Default, serde::Serialize)]
pub struct BrokerReport {
    /// Configured behavior after the first stable split.
    pub steady_state_policy: SteadyStatePolicy,
    /// True when the controller terminated itself at convergence instead of
    /// continuing to sample until the application called `finish`.
    pub monitoring_stopped_after_convergence: bool,
    /// Number of controller wakeups, including startup warm-up windows.
    pub controller_samples: u64,
    /// Wall time for which the controller remained alive.
    pub controller_elapsed: Duration,
    /// Controller thread CPU time across its complete lifetime. This is
    /// measured with only one thread-clock read at start and one at exit.
    pub controller_cpu_nanos: Option<u64>,
    /// Controller thread-clock reads that failed.
    pub controller_cpu_accounting_failures: usize,
    /// Effective cadence between responsive steady-state probes.
    pub steady_probe_interval: Duration,
    /// Cost and result of the optional startup-only opening bracket.
    pub opening_bracket: OpeningBracketReport,
    /// Consumer threads at each decision, in order.
    pub consumer_trajectory: Vec<usize>,
    /// Producer limit at each decision, in order.
    pub producer_trajectory: Vec<usize>,
    /// Actual live/active occupancy observed at every recorded sample.
    pub actual_consumer_trajectory: Vec<usize>,
    pub actual_producer_trajectory: Vec<usize>,
    /// Controller elapsed time for each retained trajectory point.
    pub trajectory_elapsed: Vec<Duration>,
    /// Times at which the inferred producer ceiling or its evidence changed.
    pub cap_trajectory: Vec<CapObservation>,
    /// Real pipeline-throughput windows retained while a decision was open.
    /// Steady-state probes are excluded, so this remains bounded startup or
    /// resurvey evidence rather than recurring telemetry overhead.
    pub throughput_trace: Vec<f64>,
    pub throughput_trace_dropped: usize,
    /// Oldest samples discarded after reaching `trace_capacity`.
    pub trace_samples_dropped: usize,
    /// Splits applied, including reverts.
    pub moves: usize,
    /// Moves undone because throughput got worse.
    pub reverts: usize,
    /// Guarded local candidates attempted after the opening rejected or could
    /// not distinguish the model's differing answer.
    pub nonlinear_probes: usize,
    /// Probe targets retained because their throughput confidence interval was
    /// strictly better than the previous split.
    pub nonlinear_probe_improvements: usize,
    /// Whether the final split is held by measured opening-bracket evidence
    /// rather than the allocation-independent cost model.
    pub nonlinear_override: bool,
    /// Ratifications whose uncertainty interval could not establish a
    /// regression or a clear non-regression.
    pub inconclusive_ratifications: usize,
    /// Times a settled split was re-opened by a detected regime change.
    pub resurveys: usize,
    /// Windows discarded as carrying too little work to estimate a rate from.
    pub rejected_windows: usize,
    /// Otherwise-usable windows rejected because producer buffers filled or
    /// drained materially, so the two stages did not process equivalent work.
    pub flow_transient_windows: usize,
    /// Samples where the run was source-bound rather than decode-bound.
    pub source_bound_samples: usize,
    /// Samples where the producer could not use more concurrency at all.
    pub inelastic_samples: usize,
    /// Model-requested producer shrinks vetoed because runnable producer work
    /// was still queued behind the current limit.
    pub pressure_vetoed_shrinks: usize,
    /// Time from the end of warm-up until the broker settled, if it did.
    pub time_to_converge: Option<Duration>,
    /// Final split.
    pub final_consumer_threads: usize,
    pub final_producer_limit: usize,
    /// Actual controlled occupancy observed when the controller stopped. Under
    /// freeze-after-convergence this is earlier than application shutdown.
    pub final_consumer_live: usize,
    pub final_producer_active: usize,
    /// Largest observed sum of live consumer workers and active producer slots.
    pub peak_controlled_slots: usize,
    /// Peak live threads explicitly outside the controlled slot budget: the
    /// controller itself plus producer-reported coordinators/optional samplers.
    pub peak_auxiliary_threads: usize,
    /// Producer measurement diagnostics. Native cumulative counters report
    /// zero samples and observation time; sampled adapters also report their
    /// cadence and sampler CPU.
    pub producer_measurement: Option<ProducerMeasurementStats>,
    /// Rejected targets and the workload epoch in which they failed.
    pub rejections: Vec<Rejection>,
    pub final_epoch: u64,
    pub final_phase: BrokerPhase,
    pub current_drift: Duration,
    pub terminal_error: Option<String>,
    /// Last solved model, for diagnosing a split that looks wrong.
    pub final_model: Option<Model>,
}

/// Auditable cost and result of [`OpeningPolicy::Bracket`].
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Serialize)]
pub struct OpeningBracketReport {
    /// Candidate allocations whose evidence horizon completed. The already
    /// measured opening is not counted.
    pub points_measured: usize,
    /// Clean throughput samples consumed by candidate evidence horizons.
    pub samples: usize,
    /// Wall time from first stable disagreement until the bracket completed.
    pub wall_nanos: u64,
    pub outcome: OpeningBracketOutcome,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum OpeningBracketOutcome {
    #[default]
    NotConfigured,
    /// The first stable model answer matched the opening, so no bracket work
    /// was performed.
    ModelAgreed,
    /// The model's differing answer was measured and retained.
    ModelSelected,
    /// The model answer and local alternatives failed to beat the opening.
    OpeningRetained,
    /// A measured local alternative beat both the model answer and opening.
    AlternativeSelected,
    /// The configured wall budget ended the experiment.
    BudgetExhausted,
    /// Mapping ended before the experiment reached a decision.
    Incomplete,
    /// The model-only freeze policy deliberately skipped calibration.
    SkippedBySteadyPolicy,
}

#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize)]
pub struct Rejection {
    pub epoch: u64,
    pub producer_target: usize,
    pub baseline_rate: f64,
    pub achieved_rate: f64,
    pub baseline_uncertainty: f64,
    pub achieved_uncertainty: f64,
    pub producer_cost_share: f64,
}

impl BrokerReport {
    /// Whether the broker settled and stayed settled.
    pub fn converged(&self) -> bool {
        self.time_to_converge.is_some()
            && self.final_phase == BrokerPhase::Steady
            && self.current_drift.is_zero()
            && self.terminal_error.is_none()
    }

    fn record_sample(
        &mut self,
        requested: Split,
        actual: Split,
        elapsed: Duration,
        capacity: usize,
    ) {
        let unchanged = self.consumer_trajectory.last() == Some(&requested.0)
            && self.producer_trajectory.last() == Some(&requested.1)
            && self.actual_consumer_trajectory.last() == Some(&actual.0)
            && self.actual_producer_trajectory.last() == Some(&actual.1);
        if unchanged {
            return;
        }
        if self.consumer_trajectory.len() == capacity {
            self.consumer_trajectory.remove(0);
            self.producer_trajectory.remove(0);
            self.actual_consumer_trajectory.remove(0);
            self.actual_producer_trajectory.remove(0);
            self.trajectory_elapsed.remove(0);
            self.trace_samples_dropped += 1;
        }
        self.consumer_trajectory.push(requested.0);
        self.producer_trajectory.push(requested.1);
        self.actual_consumer_trajectory.push(actual.0);
        self.actual_producer_trajectory.push(actual.1);
        self.trajectory_elapsed.push(elapsed);
    }

    fn record_cap(
        &mut self,
        elapsed: Duration,
        useful_cap: usize,
        reason: ProducerCapReason,
        capacity: usize,
    ) {
        let unchanged = self.cap_trajectory.last().is_some_and(|observation| {
            observation.useful_cap == useful_cap && observation.reason == reason
        });
        if unchanged {
            return;
        }
        if self.cap_trajectory.len() == capacity {
            self.cap_trajectory.remove(0);
        }
        self.cap_trajectory.push(CapObservation {
            elapsed,
            useful_cap,
            reason,
        });
    }

    fn record_throughput(&mut self, throughput: f64, capacity: usize) {
        if self.throughput_trace.len() == capacity {
            self.throughput_trace.remove(0);
            self.throughput_trace_dropped += 1;
        }
        self.throughput_trace.push(throughput);
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
pub struct CapObservation {
    pub elapsed: Duration,
    pub useful_cap: usize,
    pub reason: ProducerCapReason,
}

/// The solved model at one instant.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize)]
pub struct Model {
    /// Fraction of the pipeline's total per-item cost that falls on the
    /// producer: `s_p / (s_p + s_c)`. The split is this times the budget.
    pub producer_cost_share: f64,
    /// Approximate 95% half-width of the cost-share estimate across retained
    /// clean windows.
    pub producer_cost_share_uncertainty: f64,
    /// Producer slots the model asks for, before clamping.
    pub ideal_producer_slots: usize,
    /// Producer slots above which more concurrency has been observed not to
    /// help. `usize::MAX` when nothing has been observed to bound it.
    pub useful_cap: usize,
    /// Evidence that established `useful_cap`.
    pub useful_cap_reason: ProducerCapReason,
    /// Movement thresholds after scaling the configured absolute floors by the
    /// observed uncertainty at this budget.
    pub effective_deadband_threads: usize,
    pub effective_resurvey_distance: usize,
}

/// Evidence behind an empirical producer-concurrency ceiling.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum ProducerCapReason {
    #[default]
    None,
    /// The producer persistently left granted execution slots idle.
    Slack,
    /// The producer persistently reported source-bound or inelastic pressure.
    Source,
    /// Both signals independently support a ceiling.
    SlackAndSource,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum BrokerPhase {
    #[default]
    Starting,
    Survey,
    Draining,
    Blackout,
    Ratifying,
    Steady,
}

/// A broker sampling in the background. Dropping it stops the loop.
pub struct RunningBroker {
    stop: Arc<Stop>,
    join: Option<JoinHandle<Result<BrokerReport, BrokerError>>>,
}

impl RunningBroker {
    /// Stop sampling and collect the report.
    pub fn finish(mut self) -> Result<BrokerReport, BrokerError> {
        self.stop.set();
        let Some(join) = self.join.take() else {
            return Err(BrokerError::new(BrokerErrorKind::ThreadPanicked));
        };
        join.join()
            .map_err(|_| BrokerError::new(BrokerErrorKind::ThreadPanicked))?
    }
}

impl Drop for RunningBroker {
    fn drop(&mut self) {
        self.stop.set();
        if let Some(h) = self.join.take() {
            match h.join() {
                Ok(Ok(_)) => {}
                Ok(Err(error)) => tracing::error!(%error, "thread broker failed while dropping"),
                Err(_) => tracing::error!("thread-broker sampler panicked while dropping"),
            }
        }
    }
}

/// A split: consumer threads, producer slots.
type Split = (usize, usize);

/// Where the loop is.
#[derive(Debug, Clone)]
enum Phase {
    /// Gathering cost evidence at the current split.
    Survey,
    /// The shrinking side has accepted its new target; wait until its actual
    /// occupancy releases the slots before growing the other side.
    Draining {
        side: ResizeSide,
        to: Split,
        started: Instant,
        next: Box<Phase>,
    },
    /// A move was just applied; measurements are contaminated until it settles.
    Blackout { left: usize, next: Box<Phase> },
    /// The opening bracket exhausted its wall budget. Restore the last proven
    /// split before returning to ordinary steady-state control.
    OpeningBracketRestore { restore: Split },
    /// The move has settled. Is throughput actually better?
    Ratify {
        left: usize,
        from: Split,
        target: usize,
        baseline: RateEstimate,
        rates: Vec<f64>,
        kind: RatificationKind,
    },
    /// Settled. Sleeping until the model's answer moves a long way, and stays
    /// moved.
    Steady { drifted_for: Duration },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RatificationKind {
    Model,
    OpeningModel,
    NonlinearProbe,
}

impl Phase {
    fn public(&self) -> BrokerPhase {
        match self {
            Self::Survey => BrokerPhase::Survey,
            Self::Draining { .. } => BrokerPhase::Draining,
            Self::Blackout { .. } | Self::OpeningBracketRestore { .. } => BrokerPhase::Blackout,
            Self::Ratify { .. } => BrokerPhase::Ratifying,
            Self::Steady { .. } => BrokerPhase::Steady,
        }
    }

    fn drift_duration(&self) -> Duration {
        match self {
            Self::Steady { drifted_for } => *drifted_for,
            _ => Duration::ZERO,
        }
    }

    fn measurement_mode(&self) -> ProducerMeasurementMode {
        match self {
            Self::Survey | Self::Ratify { .. } | Self::OpeningBracketRestore { .. } => {
                ProducerMeasurementMode::Calibration
            }
            Self::Draining { .. } | Self::Blackout { .. } | Self::Steady { .. } => {
                ProducerMeasurementMode::Monitoring
            }
        }
    }
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
    pub fn start(self) -> Result<RunningBroker, BrokerError> {
        let floor = self.config.min_producer_slots.max(1);
        let ceiling = self.budget - self.config.min_consumer_threads;
        let coarse = (self.budget / 8).clamp(floor, ceiling);
        let initial_producer = self
            .initial_producer_slots
            .unwrap_or(coarse)
            .clamp(floor, ceiling);
        let initial_consumer = self.budget - initial_producer;

        self.consumer
            .set_threads(initial_consumer)
            .map_err(|source| {
                BrokerError::new(BrokerErrorKind::ResizeRefused {
                    side: ResizeSide::Consumer,
                    requested: initial_consumer,
                    source,
                })
            })?;
        self.producer
            .set_limit(initial_producer)
            .map_err(|source| {
                BrokerError::new(BrokerErrorKind::ResizeRefused {
                    side: ResizeSide::Producer,
                    requested: initial_producer,
                    source,
                })
            })?;
        let initial = (initial_consumer, initial_producer);

        let stop = Stop::new();
        let stop_thread = Arc::clone(&stop);
        let join = std::thread::Builder::new()
            .name("thread-broker".into())
            .spawn(move || self.run(stop_thread, initial))
            .map_err(|source| BrokerError::new(BrokerErrorKind::ThreadSpawn(source)))?;

        Ok(RunningBroker {
            stop,
            join: Some(join),
        })
    }

    fn run(self, stop: Arc<Stop>, initial: Split) -> Result<BrokerReport, BrokerError> {
        let controller_cpu = crate::ThreadCpuTimer::start();
        let cfg = self.config;
        let steady_probe_interval = cfg.steady_probe_interval.unwrap_or(cfg.sample_interval);
        let opening_bracket = match cfg.opening_policy {
            OpeningPolicy::Fixed => None,
            OpeningPolicy::Bracket(bracket) => Some(bracket),
        };
        let opening_bracket_enabled = opening_bracket.is_some()
            && matches!(
                self.steady_state_policy,
                SteadyStatePolicy::Responsive | SteadyStatePolicy::FreezeAfterFullCalibration
            );
        let mut report = BrokerReport {
            steady_state_policy: self.steady_state_policy,
            steady_probe_interval,
            opening_bracket: OpeningBracketReport {
                outcome: match (opening_bracket, opening_bracket_enabled) {
                    (None, _) => OpeningBracketOutcome::NotConfigured,
                    (Some(_), false) => OpeningBracketOutcome::SkippedBySteadyPolicy,
                    (Some(_), true) => OpeningBracketOutcome::Incomplete,
                },
                ..OpeningBracketReport::default()
            },
            ..BrokerReport::default()
        };

        let mut measurement_mode = ProducerMeasurementMode::Calibration;
        self.producer
            .set_monitoring_interval((steady_probe_interval / 4).max(Duration::from_nanos(1)));
        self.producer.set_measurement_mode(measurement_mode);
        self.consumer.set_measurement_mode(measurement_mode);

        let started = now();
        let mut prev_consumer = self.consumer.work();
        let mut prev_producer = self.producer.work();
        let mut prev_buffer = self.producer.buffered_items();
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
        let opening_bracket_samples = opening_bracket
            .map(|bracket| duration_samples(bracket.horizon, cfg.sample_interval))
            .unwrap_or(0);
        let mut recent = Recent::new(cfg.ratify_samples.max(opening_bracket_samples));
        let mut caps = Caps::new(cfg.cap_history, cfg.cap_persistence);
        // Targets that were tried and made throughput worse. FDP calls this the
        // tabu set; without it a model that keeps recomputing the same answer
        // re-applies a move that has already been measured as a regression.
        let mut rejected: HashMap<usize, f64> = HashMap::new();
        let mut epoch = 0u64;
        let mut phase = Phase::Survey;
        let mut warm_ended: Option<Instant> = None;
        let mut surveyed = 0usize;
        let mut settled_rate: Option<f64> = None;
        let mut pressure_vetoed_for = Duration::ZERO;
        // The full responsive policy pays a bounded startup probe when the
        // isolated cost model collapses to the producer floor. This catches
        // negative consumer scaling and allocation-dependent stage costs that
        // no one-point service-cost estimate can identify. The cheaper
        // convergence freeze is explicitly model-only; full-calibration freeze
        // runs the same bounded probe as responsive mode before stopping.
        let mut nonlinear_targets = std::collections::VecDeque::new();
        let mut nonlinear_probe_complete = !opening_bracket_enabled;
        let mut opening_bracket_started: Option<Instant> = None;
        let mut opening_bracket_points_started = 0usize;
        let mut nonlinear_override = false;
        // The achieved estimate at the last retained nonlinear point. Reusing
        // it makes the next probe compare adjacent points on the response
        // curve; `recent` still contains pre-probe windows and is only the
        // baseline for the first point.
        let mut nonlinear_rate: Option<RateEstimate> = None;
        // Best *proven* point in the current bounded exploration. An
        // inconclusive candidate is never made the rollback destination.
        let mut nonlinear_best_split: Option<Split> = None;
        let mut nonlinear_tried = std::collections::HashSet::new();
        // Model drift is meaningful only relative to the model that established
        // the settled split. Comparing only the newly rounded target makes a
        // stable share near a slot boundary alternate forever between adjacent
        // integers. Keep the settled evidence so reopening requires either a
        // material share change or a changed empirical cap.
        let mut settled_model: Option<Model> = None;

        while !stop.is_set() {
            // Establish convergence with one ordinary-cadence clean steady
            // window. Only an already-converged responsive broker adopts the
            // sparse cadence; otherwise a five-minute probe interval would also
            // delay convergence (and freeze) by five minutes.
            let mut interval =
                if matches!(phase, Phase::Steady { .. }) && report.time_to_converge.is_some() {
                    steady_probe_interval
                } else {
                    cfg.sample_interval
                };
            // A wall-bounded startup bracket must not inherit a long steady
            // sleep just as its deadline arrives. Wake on the deadline itself;
            // an in-flight drain remains governed by resize_timeout because it
            // is unsafe to grow the opposite pool before shrink acknowledgement.
            if !matches!(phase, Phase::Draining { .. })
                && !nonlinear_probe_complete
                && let (Some(bracket), Some(bracket_started)) =
                    (opening_bracket, opening_bracket_started)
            {
                let remaining = bracket
                    .total_budget
                    .saturating_sub(now().saturating_duration_since(bracket_started));
                interval = interval.min(remaining);
            }
            if stop.wait_timeout(interval) {
                break;
            }

            let sampled_at = now();
            report.controller_samples += 1;
            let consumer = self.consumer.work();
            let producer = self.producer.work();
            for (side, previous, observed) in [
                (ResizeSide::Consumer, prev_consumer, consumer),
                (ResizeSide::Producer, prev_producer, producer),
            ] {
                if observed.busy_nanos < previous.busy_nanos || observed.items < previous.items {
                    let actual = (self.consumer.live_threads(), self.producer.active_slots());
                    return Err(runtime_failure(
                        BrokerErrorKind::WorkCountersRegressed {
                            side,
                            previous,
                            observed,
                        },
                        report,
                        split,
                        actual,
                        phase.public(),
                        epoch,
                    ));
                }
            }
            let window = sampled_at.saturating_duration_since(prev_at);
            let dc = consumer.delta(prev_consumer);
            let dp = producer.delta(prev_producer);
            let buffer = self.producer.buffered_items();
            let flow_stable = match (prev_buffer, buffer) {
                (Some(previous), Some(current)) => {
                    let drift = previous.abs_diff(current);
                    let relative = (dp.items as f64 * cfg.max_buffer_drift_fraction).round() as u64;
                    drift <= cfg.max_buffer_drift_items.max(relative)
                }
                _ => true,
            };
            prev_consumer = consumer;
            prev_producer = producer;
            prev_buffer = buffer;
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

            let usable = usable_window(dc, dp, window);
            let resizing = matches!(phase, Phase::Draining { .. });
            let clean = usable && flow_stable && !resizing;
            // A zero-progress window is unusable for per-item service cost,
            // but it is essential throughput evidence. Dropping it conditions
            // the rate on batch completion and can make a bursty, slow split
            // look as fast as a smooth one. Keep rate and cost eligibility
            // deliberately separate.
            let rate_clean = flow_stable && !resizing && !window.is_zero();
            if rate_clean && !matches!(phase, Phase::Steady { .. }) {
                report.record_throughput(throughput, cfg.trace_capacity);
            }
            if clean {
                costs.push(dc.busy_nanos, dp.busy_nanos);
                caps.observe_slack(dp, window, split.1);
            } else if !resizing && !usable {
                report.rejected_windows += 1;
            } else if !resizing && !flow_stable {
                report.flow_transient_windows += 1;
            }
            // Source/elasticity pressure describes producer admission directly;
            // unlike busy-derived slack, it does not require equivalent work to
            // have crossed the downstream buffer. Requiring stable flow here
            // delayed a persistent source ceiling by 5--12 seconds on Gate G.
            if !resizing {
                caps.observe_pressure(window, split.1, pressure);
            }
            if rate_clean && matches!(phase, Phase::Survey | Phase::Steady { .. }) {
                recent.push(dc.items, window);
            }

            let live = self.consumer.live_threads();
            let active = self.producer.active_slots();
            report.peak_controlled_slots = report
                .peak_controlled_slots
                .max(live.saturating_add(active));
            report.peak_auxiliary_threads = report
                .peak_auxiliary_threads
                .max(self.producer.auxiliary_threads().saturating_add(1));
            report.producer_measurement = self.producer.measurement_stats();
            let idle = idle_fraction(dc.busy_nanos, live, window);

            if !matches!(phase, Phase::Survey) {
                pressure_vetoed_for = Duration::ZERO;
            }
            let opening_bracket_expired = !nonlinear_probe_complete
                && opening_bracket.is_some_and(|bracket| {
                    opening_bracket_started.is_some_and(|started| {
                        sampled_at.saturating_duration_since(started) >= bracket.total_budget
                    })
                });

            phase = match phase {
                Phase::OpeningBracketRestore { restore } => {
                    nonlinear_probe_complete = true;
                    nonlinear_targets.clear();
                    finish_opening_bracket(
                        &mut report,
                        &mut opening_bracket_started,
                        OpeningBracketOutcome::BudgetExhausted,
                        sampled_at,
                    );
                    settled_rate = nonlinear_rate.map(|rate| rate.mean);
                    recent.clear();
                    if split == restore {
                        Phase::Steady {
                            drifted_for: Duration::ZERO,
                        }
                    } else {
                        report.reverts += 1;
                        match self.begin_move(
                            split,
                            restore,
                            Phase::Blackout {
                                left: cfg.blackout_samples,
                                next: Box::new(Phase::Steady {
                                    drifted_for: Duration::ZERO,
                                }),
                            },
                        ) {
                            Ok((interim, phase)) => {
                                split = interim;
                                phase
                            }
                            Err(kind) => {
                                report.final_consumer_threads = split.0;
                                report.final_producer_limit = split.1;
                                report.final_consumer_live = live;
                                report.final_producer_active = active;
                                return Err(runtime_failure(
                                    kind,
                                    report,
                                    split,
                                    (live, active),
                                    BrokerPhase::Blackout,
                                    epoch,
                                ));
                            }
                        }
                    }
                }
                Phase::Draining {
                    side,
                    to,
                    started,
                    next,
                } => {
                    let (observed, target) = match side {
                        ResizeSide::Consumer => (live, to.0),
                        ResizeSide::Producer => (active, to.1),
                    };
                    if observed <= target {
                        let growth = match side {
                            ResizeSide::Consumer => {
                                self.producer.set_limit(to.1).map_err(|source| {
                                    BrokerErrorKind::ResizeRefused {
                                        side: ResizeSide::Producer,
                                        requested: to.1,
                                        source,
                                    }
                                })
                            }
                            ResizeSide::Producer => {
                                self.consumer.set_threads(to.0).map_err(|source| {
                                    BrokerErrorKind::ResizeRefused {
                                        side: ResizeSide::Consumer,
                                        requested: to.0,
                                        source,
                                    }
                                })
                            }
                        };
                        if let Err(kind) = growth {
                            report.final_consumer_threads = split.0;
                            report.final_producer_limit = split.1;
                            report.final_consumer_live = live;
                            report.final_producer_active = active;
                            return Err(runtime_failure(
                                kind,
                                report,
                                split,
                                (live, active),
                                BrokerPhase::Draining,
                                epoch,
                            ));
                        }
                        split = to;
                        report.moves += 1;
                        costs.clear();
                        *next
                    } else if sampled_at.saturating_duration_since(started) >= cfg.resize_timeout {
                        report.final_consumer_threads = split.0;
                        report.final_producer_limit = split.1;
                        report.final_consumer_live = live;
                        report.final_producer_active = active;
                        let kind = BrokerErrorKind::ResizeTimedOut {
                            side,
                            target,
                            observed,
                            timeout: cfg.resize_timeout,
                        };
                        return Err(runtime_failure(
                            kind,
                            report,
                            split,
                            (live, active),
                            BrokerPhase::Draining,
                            epoch,
                        ));
                    } else {
                        Phase::Draining {
                            side,
                            to,
                            started,
                            next,
                        }
                    }
                }

                Phase::Blackout { left, next } => {
                    if opening_bracket_expired {
                        Phase::OpeningBracketRestore {
                            restore: nonlinear_best_split.unwrap_or(initial),
                        }
                    } else if !rate_clean {
                        Phase::Blackout { left, next }
                    } else if left > 1 {
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
                    mut rates,
                    kind,
                } => {
                    if opening_bracket_expired {
                        Phase::OpeningBracketRestore {
                            restore: nonlinear_best_split.unwrap_or(initial),
                        }
                    } else if !rate_clean {
                        Phase::Ratify {
                            left,
                            from,
                            target,
                            baseline,
                            rates,
                            kind,
                        }
                    } else {
                        rates.push(rate_per_second(dc.items, window));
                        if matches!(
                            kind,
                            RatificationKind::OpeningModel | RatificationKind::NonlinearProbe
                        ) {
                            report.opening_bracket.samples += 1;
                        }
                        if left > 1 {
                            Phase::Ratify {
                                left: left - 1,
                                from,
                                target,
                                baseline,
                                rates,
                                kind,
                            }
                        } else {
                            let achieved = RateEstimate::from_rates(rates.iter().copied());
                            if matches!(
                                kind,
                                RatificationKind::OpeningModel | RatificationKind::NonlinearProbe
                            ) {
                                report.opening_bracket.points_measured += 1;
                            }
                            if kind == RatificationKind::NonlinearProbe {
                                let improved = nonlinear_probe_improved(
                                    baseline,
                                    achieved,
                                    cfg.regression_tolerance,
                                );
                                if improved {
                                    tracing::debug!(
                                        "thread-broker: nonlinear probe kept producer {} -> {} \
                                         ({:.0} [{:.0}, {:.0}] vs {:.0} [{:.0}, {:.0}] items/s)",
                                        from.1,
                                        target,
                                        achieved.mean,
                                        achieved.lower(),
                                        achieved.upper(),
                                        baseline.mean,
                                        baseline.lower(),
                                        baseline.upper(),
                                    );
                                    report.nonlinear_probe_improvements += 1;
                                    nonlinear_override = true;
                                    nonlinear_rate = Some(achieved);
                                    nonlinear_best_split = Some(split);
                                    settled_rate = Some(achieved.mean);
                                    costs.clear();
                                    // The model has already lost to the
                                    // opening. The first local point is on the
                                    // opposite side of that rejected answer;
                                    // once it proves an improvement, continuing
                                    // back toward the rejected model cannot
                                    // justify another startup move.
                                    nonlinear_probe_complete = true;
                                    nonlinear_targets.clear();
                                    finish_opening_bracket(
                                        &mut report,
                                        &mut opening_bracket_started,
                                        OpeningBracketOutcome::AlternativeSelected,
                                        sampled_at,
                                    );
                                    Phase::Steady {
                                        drifted_for: Duration::ZERO,
                                    }
                                } else {
                                    let comparison =
                                        compare_rates(baseline, achieved, cfg.regression_tolerance);
                                    tracing::debug!(
                                        "thread-broker: nonlinear probe did not retain producer {} \
                                         ({:.0} [{:.0}, {:.0}] vs {:.0} [{:.0}, {:.0}] items/s)",
                                        target,
                                        achieved.mean,
                                        achieved.lower(),
                                        achieved.upper(),
                                        baseline.mean,
                                        baseline.lower(),
                                        baseline.upper(),
                                    );
                                    if comparison == RatificationOutcome::Inconclusive {
                                        report.inconclusive_ratifications += 1;
                                    }
                                    // Every optional point has the higher
                                    // burden of proving improvement. A failure
                                    // does not discard a queued candidate on
                                    // the other side of the opening: that point
                                    // is what turns a one-sided guess into a
                                    // bracket.
                                    costs.clear();
                                    if nonlinear_targets.is_empty() {
                                        nonlinear_probe_complete = true;
                                        let restore = nonlinear_best_split.unwrap_or(from);
                                        let outcome = if nonlinear_override {
                                            OpeningBracketOutcome::AlternativeSelected
                                        } else {
                                            OpeningBracketOutcome::OpeningRetained
                                        };
                                        finish_opening_bracket(
                                            &mut report,
                                            &mut opening_bracket_started,
                                            outcome,
                                            sampled_at,
                                        );
                                        settled_rate = nonlinear_rate.map(|rate| rate.mean);
                                        nonlinear_best_split = None;
                                        recent.clear();
                                        if split == restore {
                                            Phase::Steady {
                                                drifted_for: Duration::ZERO,
                                            }
                                        } else {
                                            report.reverts += 1;
                                            match self.begin_move(
                                                split,
                                                restore,
                                                Phase::Blackout {
                                                    left: cfg.blackout_samples,
                                                    next: Box::new(Phase::Steady {
                                                        drifted_for: Duration::ZERO,
                                                    }),
                                                },
                                            ) {
                                                Ok((interim, phase)) => {
                                                    split = interim;
                                                    phase
                                                }
                                                Err(error_kind) => {
                                                    report.final_consumer_threads = split.0;
                                                    report.final_producer_limit = split.1;
                                                    report.final_consumer_live = live;
                                                    report.final_producer_active = active;
                                                    return Err(runtime_failure(
                                                        error_kind,
                                                        report,
                                                        split,
                                                        (live, active),
                                                        BrokerPhase::Ratifying,
                                                        epoch,
                                                    ));
                                                }
                                            }
                                        }
                                    } else {
                                        // Restore the last proven point before
                                        // measuring the queued candidate on the
                                        // other side of the opening. Otherwise
                                        // it inherits cold-start cost from a
                                        // failed endpoint while being compared
                                        // with a warm baseline.
                                        let restore = nonlinear_best_split.unwrap_or(from);
                                        if split != restore {
                                            match self.begin_move(
                                                split,
                                                restore,
                                                Phase::Blackout {
                                                    left: cfg.blackout_samples,
                                                    next: Box::new(Phase::Steady {
                                                        drifted_for: Duration::ZERO,
                                                    }),
                                                },
                                            ) {
                                                Ok((interim, phase)) => {
                                                    split = interim;
                                                    phase
                                                }
                                                Err(error_kind) => {
                                                    report.final_consumer_threads = split.0;
                                                    report.final_producer_limit = split.1;
                                                    report.final_consumer_live = live;
                                                    report.final_producer_active = active;
                                                    return Err(runtime_failure(
                                                        error_kind,
                                                        report,
                                                        split,
                                                        (live, active),
                                                        BrokerPhase::Ratifying,
                                                        epoch,
                                                    ));
                                                }
                                            }
                                        } else {
                                            Phase::Steady {
                                                drifted_for: Duration::ZERO,
                                            }
                                        }
                                    }
                                }
                            } else {
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
                                let comparison =
                                    compare_rates(baseline, achieved, cfg.regression_tolerance);
                                let opening_model_validation =
                                    kind == RatificationKind::OpeningModel;
                                let unresolved_opening_model = opening_model_validation
                                    && comparison == RatificationOutcome::Inconclusive;
                                if comparison == RatificationOutcome::Regressed
                                    || unresolved_opening_model
                                {
                                    tracing::debug!(
                                        "thread-broker: reverting producer {} -> {} after {} \
                                         ({:.0} vs {:.0} items/s)",
                                        target,
                                        from.1,
                                        if unresolved_opening_model {
                                            "inconclusive opening-model validation"
                                        } else {
                                            "measured regression"
                                        },
                                        achieved.mean,
                                        baseline.mean,
                                    );
                                    let rejected_share = report
                                        .final_model
                                        .map(|model| model.producer_cost_share)
                                        .unwrap_or_default();
                                    if comparison == RatificationOutcome::Regressed {
                                        rejected.insert(target, rejected_share);
                                        report.rejections.push(Rejection {
                                            epoch,
                                            producer_target: target,
                                            baseline_rate: baseline.mean,
                                            achieved_rate: achieved.mean,
                                            baseline_uncertainty: baseline.half_width,
                                            achieved_uncertainty: achieved.half_width,
                                            producer_cost_share: rejected_share,
                                        });
                                    } else {
                                        report.inconclusive_ratifications += 1;
                                    }
                                    // A failed model candidate makes the
                                    // opening the retained baseline. Spend the
                                    // remaining point budget on the adjacent
                                    // candidate away from the rejected model,
                                    // then the adjacent candidate toward it if
                                    // the first local test fails.
                                    // This brackets both t32 (model wins below
                                    // the opening) and t8 (the peak is just
                                    // above it) without a geometric sweep.
                                    if opening_model_validation {
                                        nonlinear_rate = Some(baseline);
                                        nonlinear_best_split = Some(from);
                                        nonlinear_tried.insert(from.1);
                                        nonlinear_tried.insert(target);
                                        let bracket = opening_bracket.expect(
                                            "opening model validation requires bracket config",
                                        );
                                        nonlinear_targets = opening_bracket_targets(
                                            self.budget,
                                            cfg.min_consumer_threads,
                                            cfg.min_producer_slots,
                                            from.1,
                                            target,
                                            bracket
                                                .max_points
                                                .saturating_sub(opening_bracket_points_started),
                                        );
                                    }
                                    settled_rate = Some(baseline.mean);
                                    report.reverts += 1;
                                    costs.clear();
                                    match self.begin_move(
                                        split,
                                        from,
                                        Phase::Blackout {
                                            left: cfg.blackout_samples,
                                            next: Box::new(Phase::Steady {
                                                drifted_for: Duration::ZERO,
                                            }),
                                        },
                                    ) {
                                        Ok((interim, phase)) => {
                                            split = interim;
                                            phase
                                        }
                                        Err(kind) => {
                                            report.final_consumer_threads = split.0;
                                            report.final_producer_limit = split.1;
                                            report.final_consumer_live = live;
                                            report.final_producer_active = active;
                                            return Err(runtime_failure(
                                                kind,
                                                report,
                                                split,
                                                (live, active),
                                                BrokerPhase::Ratifying,
                                                epoch,
                                            ));
                                        }
                                    }
                                } else {
                                    if comparison == RatificationOutcome::Inconclusive {
                                        report.inconclusive_ratifications += 1;
                                    }
                                    // The differing model answer was measured
                                    // and retained. The bracket has answered
                                    // its only question, so it must not add an
                                    // exploratory tax after agreement.
                                    if opening_model_validation {
                                        nonlinear_probe_complete = true;
                                        nonlinear_targets.clear();
                                        nonlinear_rate = None;
                                        nonlinear_best_split = None;
                                        finish_opening_bracket(
                                            &mut report,
                                            &mut opening_bracket_started,
                                            OpeningBracketOutcome::ModelSelected,
                                            sampled_at,
                                        );
                                    }
                                    // Kept. Re-survey once: the rates may have been
                                    // measured under a split that distorted them, and
                                    // DS2 reports convergence in at most a few such
                                    // rounds. The deadband and the tabu set are what
                                    // stop this from cycling.
                                    surveyed = 0;
                                    settled_rate = None;
                                    Phase::Survey
                                }
                            }
                        }
                    }
                }

                Phase::Steady { .. } if opening_bracket_expired => Phase::OpeningBracketRestore {
                    restore: nonlinear_best_split.unwrap_or(initial),
                },

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
                    let bracket = opening_bracket;
                    let bracket_budget_expired = bracket.is_some_and(|bracket| {
                        opening_bracket_started.is_some_and(|started| {
                            sampled_at.saturating_duration_since(started)
                                + bracket.horizon.saturating_mul(2)
                                > bracket.total_budget
                        })
                    });
                    let bracket_has_point_budget = bracket
                        .is_some_and(|bracket| opening_bracket_points_started < bracket.max_points);
                    let probe_target = if !nonlinear_probe_complete
                        && !bracket_budget_expired
                        && bracket_has_point_budget
                    {
                        let target = loop {
                            match nonlinear_targets.pop_front() {
                                Some(target)
                                    if target != split.1 && !nonlinear_tried.contains(&target) =>
                                {
                                    break Some(target);
                                }
                                Some(_) => {}
                                None => break None,
                            }
                        };
                        if target.is_none() {
                            nonlinear_probe_complete = true;
                            let outcome = if nonlinear_override {
                                OpeningBracketOutcome::AlternativeSelected
                            } else {
                                OpeningBracketOutcome::OpeningRetained
                            };
                            finish_opening_bracket(
                                &mut report,
                                &mut opening_bracket_started,
                                outcome,
                                sampled_at,
                            );
                        }
                        target
                    } else {
                        if !nonlinear_probe_complete {
                            nonlinear_probe_complete = true;
                            nonlinear_targets.clear();
                            let outcome = if bracket_budget_expired {
                                OpeningBracketOutcome::BudgetExhausted
                            } else if nonlinear_override {
                                OpeningBracketOutcome::AlternativeSelected
                            } else {
                                OpeningBracketOutcome::OpeningRetained
                            };
                            finish_opening_bracket(
                                &mut report,
                                &mut opening_bracket_started,
                                outcome,
                                sampled_at,
                            );
                        }
                        None
                    };

                    if let Some(target) = probe_target {
                        opening_bracket_points_started += 1;
                        report.nonlinear_probes += 1;
                        let from = split;
                        let to = (self.budget - target, target);
                        // A local bracket point is optional and must prove a real
                        // improvement, unlike a model-directed move which only
                        // has to avoid regression. Progress publishers must be
                        // fine grained enough that the ordinary ratification
                        // horizon is representative; stretching the horizon to
                        // compensate for bursty publication only delays the
                        // decision and hides the real measurement problem.
                        let probe_samples = opening_bracket_samples;
                        let probe_blackout_samples = cfg.blackout_samples;
                        let baseline = nonlinear_rate.unwrap_or_else(|| recent.estimate());
                        nonlinear_rate.get_or_insert(baseline);
                        nonlinear_best_split.get_or_insert(from);
                        nonlinear_tried.insert(from.1);
                        nonlinear_tried.insert(target);
                        tracing::debug!(
                            "thread-broker: probing opening bracket producer {} -> {}",
                            from.1,
                            target,
                        );
                        match self.begin_move(
                            from,
                            to,
                            Phase::Blackout {
                                left: probe_blackout_samples,
                                next: Box::new(Phase::Ratify {
                                    left: probe_samples,
                                    from,
                                    target,
                                    baseline,
                                    rates: Vec::with_capacity(probe_samples),
                                    kind: RatificationKind::NonlinearProbe,
                                }),
                            },
                        ) {
                            Ok((interim, phase)) => {
                                split = interim;
                                phase
                            }
                            Err(kind) => {
                                report.final_consumer_threads = split.0;
                                report.final_producer_limit = split.1;
                                report.final_consumer_live = live;
                                report.final_producer_active = active;
                                return Err(runtime_failure(
                                    kind,
                                    report,
                                    split,
                                    (live, active),
                                    BrokerPhase::Steady,
                                    epoch,
                                ));
                            }
                        }
                    } else {
                        let model_drifted = !nonlinear_override
                            && solved.as_ref().is_some_and(|s| {
                                let far_enough = s.target.abs_diff(split.1)
                                    >= s.snapshot.effective_resurvey_distance;
                                let changed_since_settling = settled_model.is_none_or(|settled| {
                                    let combined_uncertainty = settled
                                        .producer_cost_share_uncertainty
                                        + s.snapshot.producer_cost_share_uncertainty;
                                    let share_changed = (settled.producer_cost_share
                                        - s.snapshot.producer_cost_share)
                                        .abs()
                                        > combined_uncertainty.max(0.05);
                                    let cap_changed = settled.useful_cap != s.snapshot.useful_cap
                                        && (settled.useful_cap_reason != ProducerCapReason::None
                                            || s.snapshot.useful_cap_reason
                                                != ProducerCapReason::None);
                                    share_changed || cap_changed
                                });
                                let pressure_allows =
                                    !(s.target < split.1 && pressure == ProducerPressure::Starved);
                                let changed_since_rejection =
                                    rejected.get(&s.target).is_none_or(|rejected_share| {
                                        (s.snapshot.producer_cost_share - rejected_share).abs()
                                            > (2.0 * s.snapshot.producer_cost_share_uncertainty)
                                                .max(0.05)
                                    });
                                far_enough
                                    && changed_since_settling
                                    && pressure_allows
                                    && changed_since_rejection
                            });
                        let current_rate = recent.estimate();
                        // A common-mode rate change with an unchanged cost share
                        // does not imply a different split. Reopening on throughput
                        // alone made ordinary machine-load variation (and the tail
                        // of a finite input) look like a workload epoch. The rate
                        // detector remains useful only when it can expire rejection
                        // evidence: a previously bad target may become valid after
                        // the response surface changes even if the local cost share
                        // rounds to the same split.
                        let throughput_drifted = settled_rate.is_some_and(|baseline| {
                            if baseline <= 0.0 {
                                return false;
                            }
                            let regressed = current_rate.upper() < baseline * 0.9;
                            let improved = current_rate.lower() > baseline * 1.1;
                            // After measuring a nonlinear response surface, an
                            // upside-only rate change at the retained split is
                            // not evidence that another split is better. The
                            // first scATAC probe otherwise repeatedly reopened
                            // while its reduced consumer set warmed up. A real
                            // loss still reopens it. Model rejections retain the
                            // two-sided test so changed work can expire evidence
                            // against a previously rejected target.
                            (nonlinear_override && regressed)
                                || (!nonlinear_override
                                    && !rejected.is_empty()
                                    && (regressed || improved))
                        });
                        let next_drift =
                            consecutive_drift_duration(drifted_for, window, cfg.sample_interval);
                        match model_drifted || throughput_drifted {
                            // Distance alone is not enough: the noise in the solved
                            // target grows with the budget while a fixed distance
                            // does not, so at 64 threads ordinary jitter clears a
                            // band that 32 threads never approaches. Persistence is
                            // what tells the two apart -- noise does not survive a
                            // second of consecutive windows, a regime change does.
                            true if next_drift >= cfg.resurvey_persistence => {
                                tracing::debug!(
                                    "thread-broker: persistent workload change at producer split {} \
                                 (model target {}, throughput {:.0}/s); re-opening",
                                    split.1,
                                    solved.as_ref().map_or(split.1, |model| model.target),
                                    current_rate.mean,
                                );
                                report.resurveys += 1;
                                epoch += 1;
                                rejected.clear();
                                nonlinear_override = false;
                                nonlinear_rate = None;
                                nonlinear_best_split = None;
                                nonlinear_tried.clear();
                                // Opening calibration is deliberately
                                // startup-only. A later regime change reuses
                                // the ordinary model/ratification path and
                                // cannot re-arm the bracket or its fine
                                // measurement cadence.
                                nonlinear_targets.clear();
                                nonlinear_probe_complete = true;
                                report.time_to_converge = None;
                                settled_model = None;
                                // Persistent drift was itself measured from a
                                // full clean cost window. Reuse that evidence
                                // and require one fresh confirmation instead of
                                // paying another complete smoothing horizon;
                                // this keeps recovery under one second without
                                // weakening the 800 ms false-resurvey guard.
                                surveyed = cfg.smoothing_windows.saturating_sub(1);
                                settled_rate = None;
                                Phase::Survey
                            }
                            true => Phase::Steady {
                                drifted_for: next_drift,
                            },
                            false => {
                                if settled_model.is_none() {
                                    settled_model = solved.map(|model| model.snapshot);
                                }
                                if report.time_to_converge.is_none() {
                                    report.time_to_converge =
                                        Some(sampled_at.saturating_duration_since(warm_ended));
                                }
                                settled_rate.get_or_insert(current_rate.mean);
                                Phase::Steady {
                                    drifted_for: Duration::ZERO,
                                }
                            }
                        }
                    }
                }

                Phase::Survey => {
                    surveyed += clean as usize;
                    match costs.solve(self.budget, &cfg, &caps) {
                        // Not enough evidence yet, or one side has produced no
                        // work at all. Keep looking rather than guessing.
                        _ if surveyed < cfg.smoothing_windows => Phase::Survey,
                        None => Phase::Survey,
                        Some(model) => {
                            // A target that regressed is tabu only in this
                            // workload epoch. Keep the last ratified split
                            // until a proven regime change opens a new epoch,
                            // rather than walking through unsupported neighbors.
                            let mut target = if rejected.contains_key(&model.target) {
                                split.1
                            } else {
                                model.target
                            };
                            let opening_model_disagrees = opening_bracket_enabled
                                && !nonlinear_probe_complete
                                && split == initial
                                && model.target != initial.1;
                            // The cost model is intentionally primary, but it
                            // assumes allocation-independent service cost. At a
                            // low allocation, a decoder can under-report the
                            // thread-time a second slot would unlock. A live
                            // `Starved` signal may therefore veto *removing* a
                            // slot that still has runnable work behind it. It
                            // never chooses a larger target, preserving the
                            // isolated cost-share solution instead of reviving
                            // the old starvation hill climber.
                            let current_capacity_used =
                                producer_capacity_used(active, split.1, dp.busy_nanos, window);
                            let source_growth_wait = target > split.1
                                && matches!(
                                    pressure,
                                    ProducerPressure::SourceBound | ProducerPressure::Inelastic
                                )
                                && model.snapshot.useful_cap_reason == ProducerCapReason::None;
                            let pressure_veto = target < split.1
                                && pressure == ProducerPressure::Starved
                                && current_capacity_used
                                && !opening_model_disagrees;
                            if pressure_veto {
                                target = split.1;
                                report.pressure_vetoed_shrinks += 1;
                                pressure_vetoed_for += window;
                            } else {
                                pressure_vetoed_for = Duration::ZERO;
                            }
                            let distance = target.abs_diff(split.1);
                            let opening_disagrees = opening_model_disagrees;
                            if opening_bracket_enabled
                                && !nonlinear_probe_complete
                                && split == initial
                                && target == initial.1
                            {
                                finish_opening_bracket(
                                    &mut report,
                                    &mut opening_bracket_started,
                                    OpeningBracketOutcome::ModelAgreed,
                                    sampled_at,
                                );
                                nonlinear_probe_complete = true;
                            } else if opening_disagrees && opening_bracket_started.is_none() {
                                opening_bracket_started = Some(sampled_at);
                            }
                            // A single queued snapshot can be a short decode
                            // burst, especially for mapping-heavy modalities.
                            // Keep surveying until starvation has persisted as
                            // long as other cap evidence before accepting the
                            // current split as the bounded bracket correction.
                            if source_growth_wait {
                                // Direct source pressure can arrive before its
                                // duration-based cap is mature. Do not race a
                                // large growth move against that evidence; keep
                                // surveying until the signal either clears or
                                // establishes the current grant as the cap.
                                Phase::Survey
                            } else if pressure_veto && pressure_vetoed_for < cfg.cap_persistence {
                                Phase::Survey
                            } else if opening_disagrees && recent.len() < opening_bracket_samples {
                                // Compare like-sized blocks. Starting the model
                                // candidate with only the smoothing horizon in
                                // its baseline made short, bursty workloads
                                // needlessly inconclusive.
                                Phase::Survey
                            } else if opening_disagrees
                                && opening_bracket.is_some_and(|bracket| {
                                    opening_bracket_started.is_some_and(|started| {
                                        sampled_at.saturating_duration_since(started)
                                            + bracket.horizon.saturating_mul(2)
                                            > bracket.total_budget
                                    })
                                })
                            {
                                nonlinear_probe_complete = true;
                                finish_opening_bracket(
                                    &mut report,
                                    &mut opening_bracket_started,
                                    OpeningBracketOutcome::BudgetExhausted,
                                    sampled_at,
                                );
                                report.final_model = Some(model.snapshot);
                                settled_rate = Some(recent.estimate().mean);
                                Phase::Steady {
                                    drifted_for: Duration::ZERO,
                                }
                            } else if distance < model.snapshot.effective_deadband_threads
                                && !opening_disagrees
                            {
                                report.final_model = Some(model.snapshot);
                                settled_rate = Some(recent.estimate().mean);
                                Phase::Steady {
                                    drifted_for: Duration::ZERO,
                                }
                            } else {
                                let from = split;
                                let ratification_kind = if opening_disagrees {
                                    opening_bracket_points_started += 1;
                                    RatificationKind::OpeningModel
                                } else {
                                    RatificationKind::Model
                                };
                                let ratify_samples = if opening_disagrees {
                                    opening_bracket_samples
                                } else {
                                    cfg.ratify_samples
                                };
                                let blackout_samples = cfg.blackout_samples;
                                tracing::debug!(
                                    "thread-broker: solved producer {} -> {} \
                                     (cost share {:.2}, cap {})",
                                    from.1,
                                    target,
                                    model.snapshot.producer_cost_share,
                                    fmt_cap(model.snapshot.useful_cap),
                                );
                                report.final_model = Some(model.snapshot);
                                let to = (self.budget - target, target);
                                match self.begin_move(
                                    from,
                                    to,
                                    Phase::Blackout {
                                        left: blackout_samples,
                                        next: Box::new(Phase::Ratify {
                                            left: ratify_samples,
                                            from,
                                            target,
                                            baseline: recent.estimate(),
                                            rates: Vec::with_capacity(ratify_samples),
                                            kind: ratification_kind,
                                        }),
                                    },
                                ) {
                                    Ok((interim, phase)) => {
                                        split = interim;
                                        phase
                                    }
                                    Err(kind) => {
                                        report.final_consumer_threads = split.0;
                                        report.final_producer_limit = split.1;
                                        report.final_consumer_live = live;
                                        report.final_producer_active = active;
                                        return Err(runtime_failure(
                                            kind,
                                            report,
                                            split,
                                            (live, active),
                                            BrokerPhase::Survey,
                                            epoch,
                                        ));
                                    }
                                }
                            }
                        }
                    }
                }
            };

            // Keep high-resolution progress publication through the bounded
            // opening bracket, including its resize/blackout stages. A
            // slow split is precisely where coarse item batches quantize the
            // baseline. Once probing is complete, both stages can use their
            // settled cadence.
            let next_measurement_mode = if nonlinear_probe_complete {
                phase.measurement_mode()
            } else {
                ProducerMeasurementMode::Calibration
            };
            if next_measurement_mode != measurement_mode {
                self.producer.set_measurement_mode(next_measurement_mode);
                self.consumer.set_measurement_mode(next_measurement_mode);
                measurement_mode = next_measurement_mode;
            }

            report.record_sample(
                split,
                (live, active),
                sampled_at.saturating_duration_since(started),
                cfg.trace_capacity,
            );
            let (useful_cap, useful_cap_reason) = caps.useful();
            report.record_cap(
                sampled_at.saturating_duration_since(started),
                useful_cap,
                useful_cap_reason,
                cfg.trace_capacity,
            );

            tracing::trace!(
                "thread-broker: {:?} idle {:.1}% pressure {:?} {:.0} items/s \
                 consumer {} (live {live}) producer {} (active {active})",
                PhaseName(&phase),
                idle * 100.0,
                pressure,
                throughput,
                split.0,
                split.1,
            );

            // The low-overhead policy pays the full, guarded convergence cost
            // and one clean steady-state validation window, then ends both the
            // controller thread and (when it is owned only by this broker) the
            // producer's sampling adapter. There is no long timeout or dormant
            // polling loop: recurring cost after this point is exactly zero.
            if matches!(
                self.steady_state_policy,
                SteadyStatePolicy::FreezeAfterConvergence
                    | SteadyStatePolicy::FreezeAfterFullCalibration
            ) && report.time_to_converge.is_some()
                && matches!(
                    phase,
                    Phase::Steady {
                        drifted_for: Duration::ZERO
                    }
                )
            {
                report.monitoring_stopped_after_convergence = true;
                break;
            }
        }

        if report.opening_bracket.outcome == OpeningBracketOutcome::Incomplete {
            finish_opening_bracket(
                &mut report,
                &mut opening_bracket_started,
                OpeningBracketOutcome::Incomplete,
                now(),
            );
        }
        report.final_consumer_threads = split.0;
        report.final_producer_limit = split.1;
        report.final_consumer_live = self.consumer.live_threads();
        report.final_producer_active = self.producer.active_slots();
        report.peak_controlled_slots = report.peak_controlled_slots.max(
            report
                .final_consumer_live
                .saturating_add(report.final_producer_active),
        );
        report.final_epoch = epoch;
        report.final_phase = phase.public();
        report.nonlinear_override = nonlinear_override;
        report.current_drift = phase.drift_duration();
        report.controller_elapsed = now().saturating_duration_since(started);
        report.producer_measurement = self.producer.finish_measurement();
        report.controller_cpu_nanos = controller_cpu
            .elapsed()
            .and_then(|elapsed| u64::try_from(elapsed.as_nanos()).ok());
        report.controller_cpu_accounting_failures =
            usize::from(report.controller_cpu_nanos.is_none());
        Ok(report)
    }

    /// Begin a two-phase move by changing only the shrinking side.
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
    fn begin_move(
        &self,
        from: Split,
        to: Split,
        next: Phase,
    ) -> Result<(Split, Phase), BrokerErrorKind> {
        if to.1 > from.1 {
            self.consumer
                .set_threads(to.0)
                .map_err(|source| BrokerErrorKind::ResizeRefused {
                    side: ResizeSide::Consumer,
                    requested: to.0,
                    source,
                })?;
            Ok((
                (to.0, from.1),
                Phase::Draining {
                    side: ResizeSide::Consumer,
                    to,
                    started: now(),
                    next: Box::new(next),
                },
            ))
        } else {
            self.producer
                .set_limit(to.1)
                .map_err(|source| BrokerErrorKind::ResizeRefused {
                    side: ResizeSide::Producer,
                    requested: to.1,
                    source,
                })?;
            Ok((
                (from.0, to.1),
                Phase::Draining {
                    side: ResizeSide::Producer,
                    to,
                    started: now(),
                    next: Box::new(next),
                },
            ))
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

    fn estimate(&self) -> RateEstimate {
        RateEstimate::from_rates(
            self.ring
                .iter()
                .map(|&(items, elapsed)| rate_per_second(items, elapsed)),
        )
    }

    fn len(&self) -> usize {
        self.ring.len()
    }

    fn clear(&mut self) {
        self.ring.clear();
    }
}

fn duration_samples(horizon: Duration, interval: Duration) -> usize {
    let samples = horizon.as_nanos().div_ceil(interval.as_nanos().max(1));
    usize::try_from(samples).unwrap_or(usize::MAX).max(1)
}

fn finish_opening_bracket(
    report: &mut BrokerReport,
    started: &mut Option<Instant>,
    outcome: OpeningBracketOutcome,
    at: Instant,
) {
    if let Some(started) = started.take() {
        report.opening_bracket.wall_nanos =
            u64::try_from(at.saturating_duration_since(started).as_nanos()).unwrap_or(u64::MAX);
    }
    report.opening_bracket.outcome = outcome;
}

/// Candidate allocations adjacent to the opening, first away from the rejected
/// model and then toward it. The model point itself has already been measured.
/// Alternating around the opening distinguishes opposite response shapes with
/// the same model answer: at scATAC t32 the lower model point wins, whereas at
/// t8 the first point above the opening is the real peak.
fn opening_bracket_targets(
    budget: usize,
    min_consumer_threads: usize,
    min_producer_slots: usize,
    opening: usize,
    model: usize,
    maximum_points: usize,
) -> std::collections::VecDeque<usize> {
    let mut targets = std::collections::VecDeque::new();
    if maximum_points == 0 {
        return targets;
    }
    let floor = min_producer_slots;
    let ceiling = budget.saturating_sub(min_consumer_threads).max(floor);
    let toward_is_lower = model < opening;
    for distance in 1..budget {
        let lower = opening
            .checked_sub(distance)
            .filter(|target| *target >= floor);
        let upper = opening
            .checked_add(distance)
            .filter(|target| *target <= ceiling);
        let ordered = if toward_is_lower {
            [upper, lower]
        } else {
            [lower, upper]
        };
        for target in ordered.into_iter().flatten() {
            if target != model && !targets.contains(&target) {
                targets.push_back(target);
                if targets.len() == maximum_points {
                    return targets;
                }
            }
        }
    }
    targets
}

#[derive(Debug, Clone, Copy, Default)]
struct RateEstimate {
    mean: f64,
    half_width: f64,
    samples: usize,
    zero_samples: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RatificationOutcome {
    Kept,
    Inconclusive,
    Regressed,
}

fn compare_rates(
    baseline: RateEstimate,
    achieved: RateEstimate,
    tolerance: f64,
) -> RatificationOutcome {
    if baseline.mean <= 0.0 {
        return RatificationOutcome::Kept;
    }
    // One output flush or I/O pause is evidence that this block is noisy, not
    // evidence that the target is bad. Sustained zero progress remains visible
    // (and can still regress); only the isolated case is forced inconclusive.
    if achieved.samples > 1 && achieved.zero_samples == 1 {
        return RatificationOutcome::Inconclusive;
    }

    let boundary = baseline.mean * (1.0 - tolerance);
    // The blocks are independent, so uncertainty in their difference combines
    // in quadrature rather than by adding both 95% widths. Ratification asks a
    // directional question (is the new rate below the tolerated boundary?),
    // hence the one-sided 95% factor 1.645 rather than 1.96.
    let difference_half_width =
        (1.645 / 1.96) * (baseline.half_width * (1.0 - tolerance)).hypot(achieved.half_width);
    if achieved.mean + difference_half_width < boundary {
        RatificationOutcome::Regressed
    } else if achieved.mean - difference_half_width < boundary {
        RatificationOutcome::Inconclusive
    } else {
        RatificationOutcome::Kept
    }
}

/// A nonlinear probe has to be both statistically and materially better.
///
/// Applying the materiality tolerance to the baseline confidence endpoint
/// double-counts uncertainty: two already-disjoint 95% intervals can fail only
/// because the noisier endpoint is inflated a second time. Keep the questions
/// separate instead: intervals must not overlap, and means must differ by the
/// configured practical threshold.
fn nonlinear_probe_improved(
    baseline: RateEstimate,
    achieved: RateEstimate,
    tolerance: f64,
) -> bool {
    achieved.lower() > baseline.upper() && achieved.mean > baseline.mean * (1.0 + tolerance)
}

impl RateEstimate {
    fn from_rates(rates: impl Iterator<Item = f64>) -> Self {
        let rates: Vec<_> = rates.filter(|rate| rate.is_finite()).collect();
        if rates.is_empty() {
            return Self::default();
        }
        let mean = rates.iter().sum::<f64>() / rates.len() as f64;
        let half_width = if rates.len() > 1 {
            let variance = rates.iter().map(|rate| (rate - mean).powi(2)).sum::<f64>()
                / (rates.len() - 1) as f64;
            1.96 * variance.sqrt() / (rates.len() as f64).sqrt()
        } else {
            f64::INFINITY
        };
        Self {
            mean,
            half_width,
            samples: rates.len(),
            zero_samples: rates.iter().filter(|rate| **rate == 0.0).count(),
        }
    }

    fn lower(self) -> f64 {
        (self.mean - self.half_width).max(0.0)
    }

    fn upper(self) -> f64 {
        self.mean + self.half_width
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
        let window_shares: Vec<_> = self
            .ring
            .iter()
            .filter_map(|(consumer, producer)| {
                let total = consumer.saturating_add(*producer);
                (total > 0).then_some(*producer as f64 / total as f64)
            })
            .collect();
        let uncertainty = if window_shares.len() > 1 {
            let variance = window_shares
                .iter()
                .map(|sample| (sample - share).powi(2))
                .sum::<f64>()
                / (window_shares.len() - 1) as f64;
            (1.96 * variance.sqrt() / (window_shares.len() as f64).sqrt()).clamp(0.0, 0.5)
        } else {
            0.5
        };
        let ideal = (budget as f64 * share).round().max(1.0) as usize;

        let lo = cfg.min_producer_slots.max(1);
        let hi = budget.saturating_sub(cfg.min_consumer_threads).max(lo);
        let (useful, useful_cap_reason) = caps.useful();
        let target = ideal.clamp(lo, hi).min(useful.max(lo));
        let effective_deadband = cfg
            .deadband_threads
            .max((budget as f64 * uncertainty).ceil() as usize);
        let effective_resurvey = cfg
            .resurvey_distance
            .max((budget as f64 * uncertainty * 2.0).ceil() as usize)
            .max(effective_deadband);

        Some(Solved {
            target,
            snapshot: Model {
                producer_cost_share: share,
                producer_cost_share_uncertainty: uncertainty,
                ideal_producer_slots: ideal,
                useful_cap: useful,
                useful_cap_reason,
                effective_deadband_threads: effective_deadband,
                effective_resurvey_distance: effective_resurvey,
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
    history: Duration,
    persistence: Duration,
    slack_for: Duration,
    slack_age: Duration,
    slack_peak: f64,
    slack_limit: Option<usize>,
    slack_candidate: Option<(usize, usize)>,
    slack_cap: Option<usize>,
    saturated_for: Duration,
    saturated_age: Duration,
    saturated_at: Option<usize>,
}

impl Caps {
    fn new(history: Duration, persistence: Duration) -> Self {
        Self {
            history,
            persistence,
            slack_for: Duration::ZERO,
            slack_age: Duration::ZERO,
            slack_peak: 0.0,
            slack_limit: None,
            slack_candidate: None,
            slack_cap: None,
            saturated_for: Duration::ZERO,
            saturated_age: Duration::ZERO,
            saturated_at: None,
        }
    }

    /// Observe busy-derived slack only when producer and consumer progress is
    /// comparable. Buffer fill/drain can otherwise make ordinary pipeline lag
    /// look like unused producer concurrency.
    fn observe_slack(&mut self, dp: Work, window: Duration, producer_limit: usize) {
        // Concurrency the producer actually achieved this window, in
        // thread-equivalents. Busy time excludes blocking, so this is work done
        // rather than slots held.
        let achieved = dp.busy_nanos as f64 / window.as_nanos().max(1) as f64;
        // One quiet sample is a lull, not scalability evidence. Require the
        // same condition continuously for a duration independent of sampling
        // cadence before turning it into a cap.
        if self.slack_limit != Some(producer_limit) {
            self.slack_for = Duration::ZERO;
            self.slack_peak = 0.0;
            self.slack_limit = Some(producer_limit);
        }
        if producer_limit as f64 > achieved + 1.0 {
            self.slack_for += window;
            self.slack_age = Duration::ZERO;
            self.slack_peak = self.slack_peak.max(achieved);
            if self.slack_for >= self.persistence {
                let observed_cap = self.slack_peak.ceil() as usize + 1;
                if producer_limit >= observed_cap.saturating_mul(2) {
                    // A grant at least twice the sustained ceiling is already
                    // a decisive exploration. This lets an initially invisible
                    // source ceiling recover after one oversized model move.
                    self.slack_cap = Some(observed_cap);
                    self.slack_candidate = Some((producer_limit, observed_cap));
                } else {
                    match self.slack_candidate {
                        Some((prior_limit, prior_cap))
                            if prior_limit != producer_limit
                                && prior_limit > prior_cap
                                && producer_limit > observed_cap
                                && prior_cap.abs_diff(observed_cap) <= 1 =>
                        {
                            // One mildly under-filled grant can be a batching or
                            // barrier artifact: real scRNA used more concurrency
                            // and halved wall time when its grant grew. Comparable
                            // ceilings at two distinct grants, each with real
                            // headroom above its inferred ceiling, make it a
                            // hard cap. A grant equal to the inferred cap
                            // cannot confirm its own saturation.
                            self.slack_cap = Some(prior_cap.max(observed_cap));
                        }
                        Some((prior_limit, _)) if prior_limit != producer_limit => {
                            self.slack_cap = None;
                            self.slack_candidate = Some((producer_limit, observed_cap));
                        }
                        None => {
                            self.slack_candidate = Some((producer_limit, observed_cap));
                        }
                        _ => {}
                    }
                }
            }
        } else {
            self.slack_for = Duration::ZERO;
            self.slack_peak = 0.0;
            self.slack_age += window;
            if self.slack_age >= self.history {
                self.slack_candidate = None;
                self.slack_cap = None;
            }
        }
    }

    /// Observe the producer's direct source/elasticity classification.
    ///
    /// This evidence is independent of downstream buffer equivalence. It still
    /// requires duration-based persistence and expires after contradictory
    /// observations, so a lone source lull cannot become a lasting cap.
    fn observe_pressure(&mut self, window: Duration, producer_limit: usize, p: ProducerPressure) {
        // Cap history uses only the source/elasticity half of
        // `ProducerPressure`. It is a saturating signal -- it cannot tell "just
        // enough" from "far too much" -- so it may cap growth but must never
        // size the split. `Starved`'s separate bounded shrink veto lives at the
        // decision point, where it can retain but never increase the target.
        match p {
            ProducerPressure::SourceBound | ProducerPressure::Inelastic => {
                self.saturated_for += window;
                self.saturated_age = Duration::ZERO;
                if self.saturated_for >= self.persistence {
                    self.saturated_at = Some(match self.saturated_at {
                        Some(prev) => prev.min(producer_limit.max(1)),
                        None => producer_limit.max(1),
                    });
                }
            }
            _ => {
                self.saturated_for = Duration::ZERO;
                self.saturated_age += window;
                if self.saturated_age >= self.history {
                    self.saturated_at = None;
                }
            }
        }
    }

    #[cfg(test)]
    fn observe(&mut self, dp: Work, window: Duration, producer_limit: usize, p: ProducerPressure) {
        self.observe_slack(dp, window, producer_limit);
        self.observe_pressure(window, producer_limit, p);
    }

    /// The largest producer limit worth asking for, or `usize::MAX` if nothing
    /// has been observed to bound it.
    fn useful(&self) -> (usize, ProducerCapReason) {
        let mut cap = usize::MAX;

        // Slack: if the producer was handed more than a thread's worth of slack
        // and declined to use it, more slots are not what it is short of. The
        // ceiling is the best it managed while it had that slack, plus one so a
        // marginal case is not pinched.
        if let Some(slack) = self.slack_cap {
            cap = cap.min(slack);
        }

        if let Some(limit) = self.saturated_at {
            cap = cap.min(limit);
        }
        let reason = match (self.slack_cap.is_some(), self.saturated_at.is_some()) {
            (false, false) => ProducerCapReason::None,
            (true, false) => ProducerCapReason::Slack,
            (false, true) => ProducerCapReason::Source,
            (true, true) => ProducerCapReason::SlackAndSource,
        };
        (cap, reason)
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

/// Whether a `Starved` producer demonstrably used the allocation a shrink would
/// remove. Queue state alone can remain stale while already-produced buffers
/// drain, so it is not enough to invoke the nonlinear shrink veto.
fn producer_capacity_used(
    active_slots: usize,
    requested_slots: usize,
    busy_nanos: u64,
    window: Duration,
) -> bool {
    let capacity_nanos = (window.as_nanos() as u64).saturating_mul(requested_slots as u64);
    active_slots >= requested_slots
        || busy_nanos.saturating_mul(4) >= capacity_nanos.saturating_mul(3)
}

fn rate_per_second(items: u64, window: Duration) -> f64 {
    let secs = window.as_secs_f64();
    if secs <= 0.0 {
        return 0.0;
    }
    items as f64 / secs
}

/// Duration represented by consecutive drift observations.
///
/// The first sparse steady-state probe says only that drift exists *now*; it
/// cannot prove that drift was present throughout the entire sleep since the
/// previous clean probe. Credit it with one active-cadence window. A second
/// consecutive observation spans real elapsed time and may then satisfy the
/// wall-clock persistence guard.
fn consecutive_drift_duration(
    previous: Duration,
    observation_window: Duration,
    active_interval: Duration,
) -> Duration {
    if previous.is_zero() {
        observation_window.min(active_interval)
    } else {
        previous.saturating_add(observation_window)
    }
}

fn fmt_cap(cap: usize) -> String {
    if cap == usize::MAX {
        "none".into()
    } else {
        cap.to_string()
    }
}

fn runtime_failure(
    kind: BrokerErrorKind,
    mut report: BrokerReport,
    requested: Split,
    actual: Split,
    phase: BrokerPhase,
    epoch: u64,
) -> BrokerError {
    let error = BrokerError::new(kind);
    report.final_consumer_threads = requested.0;
    report.final_producer_limit = requested.1;
    report.final_consumer_live = actual.0;
    report.final_producer_active = actual.1;
    report.final_phase = phase;
    report.final_epoch = epoch;
    report.terminal_error = Some(error.to_string());
    error.with_report(report)
}

struct PhaseName<'a>(&'a Phase);

impl std::fmt::Debug for PhaseName<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self.0 {
            Phase::Survey => "survey",
            Phase::Draining { .. } => "draining",
            Phase::Blackout { .. } => "blackout",
            Phase::OpeningBracketRestore { .. } => "opening-bracket-restore",
            Phase::Ratify { .. } => "ratify",
            Phase::Steady { .. } => "steady",
        })
    }
}

#[cfg(test)]
mod tests {
    use super::{
        BrokerPhase, BrokerReport, Caps, Costs, Phase, ProducerCapReason, RateEstimate,
        RatificationKind, RatificationOutcome, Recent, compare_rates, consecutive_drift_duration,
        duration_samples, nonlinear_probe_improved, opening_bracket_targets,
        producer_capacity_used,
    };
    use crate::{BrokerConfig, ProducerMeasurementMode, ProducerPressure, ResizeSide, Work};
    use std::time::{Duration, Instant};

    struct Noise(u64);

    impl Noise {
        fn new(seed: u64) -> Self {
            Self(seed.max(1))
        }

        fn uniform(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1);
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }

        /// Irwin-Hall approximation to a standard normal. It is deterministic,
        /// portable, and sufficient for controller acceptance simulations.
        fn normal(&mut self) -> f64 {
            (0..12).map(|_| self.uniform()).sum::<f64>() - 6.0
        }
    }

    fn estimate(rate: f64, relative_sigma: f64, noise: &mut Noise) -> RateEstimate {
        RateEstimate::from_rates((0..4).map(|_| rate * (1.0 + relative_sigma * noise.normal())))
    }

    fn push_noisy_share(costs: &mut Costs, share: f64, noise: &mut Noise) {
        let observed = (share + 0.005 * noise.normal()).clamp(0.001, 0.999);
        let total = 1_000_000.0;
        costs.push(((1.0 - observed) * total) as u64, (observed * total) as u64);
    }

    fn work_at(concurrency: usize, window: Duration) -> Work {
        Work {
            busy_nanos: (window.as_nanos() as u64).saturating_mul(concurrency as u64),
            items: 1,
        }
    }

    #[test]
    fn one_slack_lull_does_not_create_a_cap() {
        let mut caps = Caps::new(Duration::from_secs(2), Duration::from_millis(300));
        let window = Duration::from_millis(100);
        caps.observe(work_at(2, window), window, 8, ProducerPressure::Satisfied);
        assert_eq!(caps.useful(), (usize::MAX, ProducerCapReason::None));
    }

    #[test]
    fn cap_confidence_uses_duration_and_expires_by_duration() {
        for window in [Duration::from_millis(25), Duration::from_millis(100)] {
            let mut caps = Caps::new(Duration::from_secs(1), Duration::from_millis(300));
            let samples = Duration::from_millis(300)
                .as_nanos()
                .div_ceil(window.as_nanos());
            for _ in 0..samples {
                caps.observe(work_at(6, window), window, 8, ProducerPressure::Satisfied);
            }
            assert_eq!(
                caps.useful(),
                (usize::MAX, ProducerCapReason::None),
                "one grant must remain provisional"
            );

            for _ in 0..samples {
                caps.observe(work_at(6, window), window, 12, ProducerPressure::Satisfied);
            }
            assert_eq!(caps.useful(), (7, ProducerCapReason::Slack));

            let expiry = Duration::from_secs(1)
                .as_nanos()
                .div_ceil(window.as_nanos());
            for _ in 0..expiry {
                caps.observe(work_at(12, window), window, 12, ProducerPressure::Starved);
            }
            assert_eq!(caps.useful(), (usize::MAX, ProducerCapReason::None));
        }
    }

    #[test]
    fn adjacent_grant_cannot_confirm_its_own_inferred_cap() {
        let mut caps = Caps::new(Duration::from_secs(1), Duration::from_millis(300));
        let window = Duration::from_millis(100);
        for _ in 0..3 {
            caps.observe(work_at(6, window), window, 8, ProducerPressure::Satisfied);
        }
        assert_eq!(caps.useful(), (usize::MAX, ProducerCapReason::None));
        for _ in 0..3 {
            caps.observe(work_at(6, window), window, 7, ProducerPressure::Satisfied);
        }
        assert_eq!(
            caps.useful(),
            (usize::MAX, ProducerCapReason::None),
            "a grant equal to the inferred cap has no headroom to confirm saturation"
        );
    }

    #[test]
    fn source_caps_arrive_and_lull_caps_recover_within_one_second() {
        for window in [Duration::from_millis(25), Duration::from_millis(100)] {
            let mut source = Caps::new(Duration::from_millis(800), Duration::from_millis(300));
            let evidence = Duration::from_millis(300)
                .as_nanos()
                .div_ceil(window.as_nanos());
            for _ in 0..evidence {
                source.observe(work_at(8, window), window, 8, ProducerPressure::SourceBound);
            }
            assert_eq!(source.useful(), (8, ProducerCapReason::Source));

            let mut lull = Caps::new(Duration::from_millis(800), Duration::from_millis(300));
            let lull_samples = Duration::from_secs(1)
                .as_nanos()
                .div_ceil(window.as_nanos());
            for _ in 0..lull_samples {
                lull.observe(work_at(0, window), window, 8, ProducerPressure::Satisfied);
            }
            assert_eq!(lull.useful(), (1, ProducerCapReason::Slack));

            let recovery_samples = Duration::from_millis(800)
                .as_nanos()
                .div_ceil(window.as_nanos());
            for _ in 0..recovery_samples {
                lull.observe(work_at(8, window), window, 8, ProducerPressure::Starved);
            }
            assert_eq!(
                lull.useful(),
                (usize::MAX, ProducerCapReason::None),
                "{window:?} sampling retained a lull cap after 800 ms of contrary evidence"
            );
        }
    }

    #[test]
    fn direct_pressure_cap_does_not_wait_for_flow_stability() {
        let mut caps = Caps::new(Duration::from_millis(800), Duration::from_millis(300));
        for _ in 0..3 {
            caps.observe_pressure(Duration::from_millis(100), 2, ProducerPressure::Inelastic);
        }
        assert_eq!(caps.useful(), (2, ProducerCapReason::Source));
    }

    #[test]
    fn interval_ratification_meets_seeded_error_gates() {
        let mut false_rejections = 0;
        let mut regressions_detected = 0;
        for seed in 1..=1_000 {
            let mut noise = Noise::new(seed);
            let baseline = estimate(100.0, 0.02, &mut noise);
            let equal = estimate(100.0, 0.02, &mut noise);
            false_rejections +=
                usize::from(compare_rates(baseline, equal, 0.05) == RatificationOutcome::Regressed);

            let baseline = estimate(100.0, 0.02, &mut noise);
            let regressed = estimate(90.0, 0.02, &mut noise);
            regressions_detected += usize::from(
                compare_rates(baseline, regressed, 0.05) == RatificationOutcome::Regressed,
            );
        }
        assert!(
            false_rejections < 10,
            "{false_rejections}/1000 equal-rate traces were confidently rejected"
        );
        assert!(
            regressions_detected >= 950,
            "only {regressions_detected}/1000 ten-percent regressions were detected"
        );
    }

    #[test]
    fn uncertainty_hysteresis_meets_seeded_budget_gates() {
        let cfg = BrokerConfig::default();
        let caps = Caps::new(cfg.cap_history, cfg.cap_persistence);

        for budget in [2usize, 4, 8, 16, 32, 64, 96] {
            let base_target = (budget / 3).clamp(1, budget - 1);
            let base_share = base_target as f64 / budget as f64;
            let mut stable_runs_without_resurvey = 0;
            let mut shifts_detected = 0;

            for seed in 1..=1_000 {
                let mut noise = Noise::new(seed ^ ((budget as u64) << 32));
                let mut costs = Costs::new(cfg.smoothing_windows);
                let mut drifted_for = Duration::ZERO;
                let mut false_resurvey = false;
                for _ in 0..300 {
                    push_noisy_share(&mut costs, base_share, &mut noise);
                    if let Some(solved) = costs.solve(budget, &cfg, &caps) {
                        if solved.target.abs_diff(base_target)
                            >= solved.snapshot.effective_resurvey_distance
                        {
                            drifted_for += Duration::from_millis(100);
                            false_resurvey |= drifted_for >= cfg.resurvey_persistence;
                        } else {
                            drifted_for = Duration::ZERO;
                        }
                    }
                }
                stable_runs_without_resurvey += usize::from(!false_resurvey);

                // A two-slot budget has only one legal split, so there is no
                // allocation change for any controller to detect.
                if budget > 2 {
                    let shift = ((budget as f64 * 0.10).ceil() as usize).max(1);
                    let changed_target = (base_target + shift).min(budget - 1);
                    let changed_share = changed_target as f64 / budget as f64;
                    drifted_for = Duration::ZERO;
                    let mut detected = false;
                    for _ in 0..20 {
                        push_noisy_share(&mut costs, changed_share, &mut noise);
                        if let Some(solved) = costs.solve(budget, &cfg, &caps) {
                            if solved.target.abs_diff(base_target)
                                >= solved.snapshot.effective_resurvey_distance
                            {
                                drifted_for += Duration::from_millis(100);
                                detected |= drifted_for >= cfg.resurvey_persistence;
                            } else {
                                drifted_for = Duration::ZERO;
                            }
                        }
                    }
                    shifts_detected += usize::from(detected);
                }
            }

            assert!(
                stable_runs_without_resurvey >= 990,
                "budget {budget}: only {stable_runs_without_resurvey}/1000 stable 30-second traces avoided resurvey"
            );
            if budget > 2 {
                assert!(
                    shifts_detected >= 990,
                    "budget {budget}: only {shifts_detected}/1000 ten-percent shifts reopened within two seconds"
                );
            }
        }
    }

    #[test]
    fn report_trace_is_bounded_and_convergence_is_terminal_state_aware() {
        let mut report = BrokerReport::default();
        for producer in 1..=10 {
            report.record_sample(
                (16 - producer, producer),
                (15 - producer, producer),
                Duration::from_millis(producer as u64),
                3,
            );
        }
        assert_eq!(report.producer_trajectory, vec![8, 9, 10]);
        assert_eq!(report.actual_producer_trajectory, vec![8, 9, 10]);
        assert_eq!(
            report.trajectory_elapsed,
            vec![
                Duration::from_millis(8),
                Duration::from_millis(9),
                Duration::from_millis(10),
            ]
        );
        assert_eq!(report.trace_samples_dropped, 7);

        report.time_to_converge = Some(Duration::from_secs(1));
        for phase in [
            BrokerPhase::Starting,
            BrokerPhase::Survey,
            BrokerPhase::Draining,
            BrokerPhase::Blackout,
            BrokerPhase::Ratifying,
        ] {
            report.final_phase = phase;
            assert!(!report.converged(), "{phase:?} was labeled converged");
        }
        report.final_phase = BrokerPhase::Steady;
        assert!(report.converged());
        report.current_drift = Duration::from_millis(100);
        assert!(!report.converged());
        report.current_drift = Duration::ZERO;
        report.terminal_error = Some("injected".into());
        assert!(!report.converged());
    }

    #[test]
    fn producer_measurement_is_dense_only_while_a_decision_is_open() {
        assert_eq!(
            Phase::Survey.measurement_mode(),
            ProducerMeasurementMode::Calibration
        );
        assert_eq!(
            Phase::Ratify {
                left: 1,
                from: (7, 1),
                target: 2,
                baseline: RateEstimate::from_rates([1.0].into_iter()),
                rates: Vec::new(),
                kind: RatificationKind::Model,
            }
            .measurement_mode(),
            ProducerMeasurementMode::Calibration
        );

        for phase in [
            Phase::Draining {
                side: ResizeSide::Producer,
                to: (6, 2),
                started: Instant::now(),
                next: Box::new(Phase::Survey),
            },
            Phase::Blackout {
                left: 1,
                next: Box::new(Phase::Survey),
            },
            Phase::Steady {
                drifted_for: Duration::ZERO,
            },
        ] {
            assert_eq!(
                phase.measurement_mode(),
                ProducerMeasurementMode::Monitoring
            );
        }
    }

    #[test]
    fn recent_throughput_keeps_zero_progress_windows() {
        let window = Duration::from_millis(100);
        let mut recent = Recent::new(3);
        recent.push(100, window);
        recent.push(0, window);
        recent.push(100, window);

        let estimate = recent.estimate();
        assert!((estimate.mean - (2_000.0 / 3.0)).abs() < 1e-9);
        assert!(estimate.lower() < estimate.mean);
        assert!(estimate.upper() > estimate.mean);
    }

    #[test]
    fn nonlinear_probe_separates_significance_from_materiality() {
        let baseline = RateEstimate {
            mean: 6_386.0,
            half_width: 3_066.0,
            ..RateEstimate::default()
        };
        let achieved = RateEstimate {
            mean: 12_453.0,
            half_width: 2_688.0,
            ..RateEstimate::default()
        };
        assert!(nonlinear_probe_improved(baseline, achieved, 0.05));

        // A large mean gain without separated intervals is not evidence.
        assert!(!nonlinear_probe_improved(
            baseline,
            RateEstimate {
                mean: 12_000.0,
                half_width: 4_000.0,
                ..RateEstimate::default()
            },
            0.05,
        ));
        // Separated intervals with a sub-threshold mean gain are not useful.
        assert!(!nonlinear_probe_improved(
            RateEstimate {
                mean: 100.0,
                half_width: 0.1,
                ..RateEstimate::default()
            },
            RateEstimate {
                mean: 104.0,
                half_width: 0.1,
                ..RateEstimate::default()
            },
            0.05,
        ));
    }

    #[test]
    fn opening_bracket_horizon_rounds_up_to_whole_samples() {
        assert_eq!(
            duration_samples(Duration::from_millis(300), Duration::from_millis(100)),
            3
        );
        assert_eq!(
            duration_samples(Duration::from_millis(301), Duration::from_millis(100)),
            4
        );
        assert_eq!(
            duration_samples(Duration::from_nanos(1), Duration::from_millis(100)),
            1
        );
    }

    #[test]
    fn opening_bracket_alternates_around_the_opening_and_is_bounded() {
        assert_eq!(
            opening_bracket_targets(32, 1, 1, 4, 1, 2)
                .into_iter()
                .collect::<Vec<_>>(),
            vec![5, 3]
        );
        assert_eq!(
            opening_bracket_targets(8, 1, 1, 4, 7, 3)
                .into_iter()
                .collect::<Vec<_>>(),
            vec![3, 5, 2]
        );
        assert!(opening_bracket_targets(8, 1, 1, 4, 1, 0).is_empty());
        assert_eq!(
            opening_bracket_targets(8, 1, 3, 4, 6, 4)
                .into_iter()
                .collect::<Vec<_>>(),
            vec![3, 5, 7],
            "candidate generation must preserve the configured safety floor"
        );
    }

    #[test]
    fn isolated_zero_progress_cannot_create_rejection_evidence() {
        let baseline = RateEstimate::from_rates([100.0; 10].into_iter());
        let mut achieved = [100.0; 10];
        achieved[4] = 0.0;
        let achieved = RateEstimate::from_rates(achieved.into_iter());
        assert_eq!(
            compare_rates(baseline, achieved, 0.05),
            RatificationOutcome::Inconclusive
        );
    }

    #[test]
    fn one_sparse_probe_cannot_satisfy_resurvey_persistence() {
        let active = Duration::from_millis(100);
        let sparse = Duration::from_secs(5);
        let first = consecutive_drift_duration(Duration::ZERO, sparse, active);
        assert_eq!(first, active);
        assert!(first < Duration::from_millis(800));

        // A second consecutive observation really does establish that the
        // drift persisted across the sparse interval.
        assert_eq!(
            consecutive_drift_duration(first, sparse, active),
            Duration::from_millis(5_100)
        );
    }

    #[test]
    fn shrink_veto_requires_evidence_that_current_capacity_was_used() {
        let window = Duration::from_millis(100);
        assert!(!producer_capacity_used(0, 2, 20_000_000, window));
        assert!(producer_capacity_used(2, 2, 0, window));
        assert!(producer_capacity_used(0, 2, 150_000_000, window));
    }
}
