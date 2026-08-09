//! Dividing one thread budget between decoding and the work that consumes it.
//!
//! # The problem
//!
//! A tool that reads compressed input has to split `-t` between threads that
//! decompress and threads that do the real work — mapping, aligning, counting.
//! The split matters enormously: measured on a 150 M read-pair FASTQ archive at
//! a fixed 64 threads, it was worth **5.75x**, while doubling the whole budget
//! bought 9%.
//!
//! It is also unpredictable. The right answer moves by an order of magnitude
//! with the workload — measured per-thread consumption ranged 0.064 GB/s on a
//! 96.3 M-k-mer index to 0.43 GB/s on a 1.5 M-k-mer one, a 6.7x spread — so no
//! constant spans it, and every attempt to *predict* it from a pre-run
//! measurement landed 7-60% off the optimum depending on the workload.
//!
//! # The approach
//!
//! Do not search for the split. Solve for it.
//!
//! Writing `s_p` and `s_c` for the two stages' per-item costs in thread-time, a
//! budget `N` finishes both sides together when `d / s_p = (N - d) / s_c`, so
//!
//! ```text
//!     d* = N * s_p / (s_p + s_c)
//! ```
//!
//! and neither cost has to be measured on its own. Over a window in which the
//! pipeline moved `X` items the producer spent `X * s_p` thread-nanoseconds and
//! the consumer `X * s_c`, so their **ratio of busy times is the ratio of costs
//! exactly** — the `X` cancels. The whole control law is
//!
//! ```text
//!     d* = N * busy_producer / (busy_producer + busy_consumer)
//! ```
//!
//! measured over the same window. No shared unit of work is needed, and a run
//! mixing compressed and uncompressed inputs is handled without a special case:
//! the producer simply does less work, and the split follows.
//!
//! The one thing that makes any of this valid is measuring **busy time, not wall
//! time**. Wall time contains blocking; blocking is a function of the split; so
//! a controller built on wall time is measuring its own actuation. Busy time is
//! a property of the work rather than of the current assignment, which is what
//! makes the equation solvable rather than circular. DS2 (Kalavri et al.,
//! OSDI'18) calls the distinction *true* versus *observed* rate. See [`Work`],
//! where the contract lives.
//!
//! An earlier version of this crate searched instead — grow whichever side looks
//! starved, one step at a time. It fails for a structural reason worth recording:
//! consumer starvation reaches zero at the balance point and *stays* zero above
//! it, so the target set is a half-line and the controller has no restoring force
//! on it. Measured, it walked to 47 producer slots of 64 and ran **44% slower
//! than the best split**, reporting near-zero error the whole way. See
//! [`ProducerPressure`] and `DESIGN-thread-broker.md`.
//!
//! A model can still be wrong — a producer reading from a saturated disk cannot
//! use the slots it is given — so the jump is guarded: capped against observed
//! saturation, ratified against measured throughput, and reverted if it made
//! things worse. Both requirements — resizing either side mid-run, and both
//! sides reporting their own busy time — are met by `paraseq` (consumer) and
//! `rapidgzip` (producer).
//!
//! Measured against swept optima on a 10x Flex dataset: **5.0% off at `-t 32`
//! (one move, settled in 1.4 s)** and **3.4% off at `-t 64` (three moves, 2.6 s)**.
//!
//! # What this crate is not
//!
//! It is not a thread pool. It owns no threads except the one that samples, and
//! it creates no workers. It reads two signals and calls two setters. The pools
//! themselves belong to the caller.
//!
//! # Wiring it up
//!
//! Implement [`Consumer`] for the side doing the real work and [`Producer`] for
//! the decompressor, then let the broker run for the duration of the job:
//!
//! ```no_run
//! # use thread_broker::*;
//! # use std::time::Duration;
//! # fn demo<C: Consumer + 'static, P: Producer + 'static>(mapper: C, decoders: P)
//! # -> Result<(), Box<dyn std::error::Error>> {
//! let broker = ThreadBroker::builder(mapper, decoders)
//!     .budget(32)
//!     .sample_interval(Duration::from_millis(100))
//!     .build()?;
//!
//! let running = broker.start()?;
//! // ... run the job ...
//! let report = running.finish()?;
//! tracing::info!(
//!     "decode settled at {} of 32 threads after {:?}",
//!     report.final_producer_limit,
//!     report.time_to_converge,
//! );
//! # Ok(())
//! # }
//! ```
//!
//! # Minimum viable adoption
//!
//! A first integration does **not** need to tune every [`BrokerConfig`] field.
//! The minimum contract is:
//!
//! 1. Define one real execution-slot budget shared by both pools. Construct the
//!    pools so every permitted split obeys that budget, choose any valid opening
//!    split, and set only `.budget(...)` plus real consumer/producer floors.
//! 2. Implement cumulative [`Work`] counters on both adapters. Busy time must
//!    exclude queue/source blocking; counters must be monotonic; progress must
//!    publish often enough that ordinary decision windows are not empty.
//! 3. Make resize acknowledgement truthful. [`Consumer::live_threads`] and
//!    [`Producer::active_slots`] must continue to include work admitted before a
//!    shrink until it actually retires. Batch geometry must let that happen
//!    inside [`BrokerConfig::resize_timeout`].
//! 4. Implement [`Producer::pressure`] as a directional classification, not a
//!    desired thread count. If the producer has native cumulative busy time,
//!    the measurement-mode and monitoring hooks can keep their defaults.
//! 5. Keep broker lifecycle failure separate from application correctness. A
//!    configuration error should fail before work starts; a runtime tuning
//!    failure normally leaves the last valid split in force and should not make
//!    successfully produced output disappear.
//!
//! Leave warm-up, sampling, smoothing, blackout, ratification, cap history,
//! deadband, resurvey, regression tolerance, and [`OpeningPolicy`] at their
//! defaults initially. Choose a [`SteadyStatePolicy`] only
//! from the application's workload contract: responsive for possible regime
//! changes, or a freeze policy for stable long runs where recurring work should
//! become zero.
//!
//! Tune cadence or batch/progress granularity only after measurement shows a
//! concrete failure: stale/empty decision windows, shrink acknowledgement near
//! the timeout, poor oracle-gap/convergence results, or material overhead.
//! Qualify a new adapter with canonical-output equality, a pinned-split sweep,
//! source-bound and inelastic cases, a mid-run regime change, and both absolute
//! administrative CPU and whole-process fractional overhead on representative
//! long runs. That is the adopter boundary; modality-specific knobs belong in
//! the application until repeated integrations demonstrate a generic policy.

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::{Duration, Instant};

mod controller;
pub use controller::{
    BrokerPhase, BrokerReport, Model, OpeningBracketOutcome, OpeningBracketReport,
    ProducerCapReason, Rejection, RunningBroker,
};

/// Whether engaging a *parallel* producer can pay for the threads it takes.
///
/// # The economics this encodes
///
/// A consumer that decodes its own input inline — each worker inflating under
/// its own stream's lock, then processing what it produced — is
/// **work-conserving**: it is never idle, because a thread with less decode work
/// simply does more consuming, rebalancing at batch granularity with no
/// coordination. It gives `streams` concurrent decoders for free, where `streams`
/// is the number of independent input streams.
///
/// A parallel producer gives `d` decoders but *takes* `d` threads away from
/// consuming, and those threads can only decode. So it pays exactly when
/// `d > streams` — when it can supply more concurrency than the consumer already
/// gets for nothing.
///
/// # Why a threshold rather than that rule
///
/// `d` is what the broker *converges to*, so it is unknown when the decision has
/// to be made. This is the budget-side proxy: with `budget` threads and
/// `streams` inputs, the producer can only plausibly exceed `streams` decoders
/// if the budget is several times the stream count.
///
/// # Provenance of the default
///
/// Measured on a single-cell Flex library (150 bp pairs, human index) and a bulk
/// SE library (gencode), sweeping `budget` over `{8, 32, 64}` and `streams` over
/// `{1..16}` — 16 cells. Inline decoding won every cell where `d <= streams` by
/// 8–28%, and the parallel producer won every cell where `d > streams` by up to
/// 4.2x. `budget >= 8 * streams` reproduces all 16 verdicts except one boundary
/// cell that costs ~8%.
///
/// The default is calibrated on the *most* decode-bound workload available, so it
/// errs toward inline decoding for anything less decode-heavy — the direction
/// that protects a consumer-bound baseline. The asymmetry is deliberate: missing
/// a parallel win costs up to 4.2x, taking one wrongly costs ~25%, but the
/// wrongly-taken case is the one a user notices as a regression against not
/// having the feature at all.
///
/// Override it when your workload's decode-to-consume ratio is known to differ.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct EngagementPolicy {
    /// Budget threads required per independent input stream before a parallel
    /// producer is engaged at all.
    ///
    /// `0` always engages; a very large value never does.
    pub min_threads_per_stream: usize,
}

impl Default for EngagementPolicy {
    fn default() -> Self {
        Self {
            min_threads_per_stream: 8,
        }
    }
}

impl EngagementPolicy {
    /// Engage a parallel producer for `streams` inputs under `budget` threads?
    ///
    /// `streams == 0` never engages: there is nothing for a parallel producer to
    /// do, so it could only cost threads.
    pub const fn engages(&self, budget: usize, streams: usize) -> bool {
        streams != 0 && budget >= self.min_threads_per_stream.saturating_mul(streams)
    }

    /// Always engage, leaving the decision entirely to the broker.
    ///
    /// There is deliberately no `never()` counterpart. This type expresses a
    /// *threshold*, and no threshold declines every budget — `usize::MAX` still
    /// engages at a `usize::MAX` budget, so such a constructor would be a lie in
    /// exactly the case someone would reach for it to be safe. A caller that
    /// wants the producer off should not construct it at all; piscem exposes
    /// that as `--decoder serial`.
    pub const fn always() -> Self {
        Self {
            min_threads_per_stream: 0,
        }
    }
}

/// Cumulative work done by one stage.
///
/// # Busy time, never wall time
///
/// `busy_nanos` must be time a worker spent *doing the work*, with blocking
/// excluded — a consumer must not count time waiting for input, a producer must
/// not count time blocked on a full output queue or on the source.
///
/// This is the single most important contract in the crate, and the whole
/// control law rests on it. Wall time contains blocking, and blocking is a
/// function of the split being controlled, so a controller built on wall time
/// feeds its own actuation back into its measurement. That endogeneity is what
/// makes utilisation-style signals degenerate — see [`ProducerPressure`].
///
/// Busy time has the opposite property: over a window in which the pipeline
/// moved `X` items, the two stages spent `X * s_producer` and `X * s_consumer`
/// thread-nanoseconds respectively, so their **ratio is the ratio of per-item
/// costs regardless of how many threads either side had**. That is the quantity
/// the split is solved from, and it is why no shared unit of work is needed:
/// the `X` cancels.
///
/// # Items are a progress measure, not a shared currency
///
/// `items` is used only to compare throughput at one split against throughput at
/// another — for ratifying a move and for noticing a regime change. It need only
/// be consistent over time within a stage. The two stages may count entirely
/// different things (records here, decompressed bytes there) without affecting
/// the split, because the split never divides one by the other.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Work {
    /// Nanoseconds of real work, summed over this stage's workers.
    pub busy_nanos: u64,
    /// Progress in whatever unit is natural for this stage.
    pub items: u64,
}

impl Work {
    fn delta(self, prev: Work) -> Work {
        Work {
            busy_nanos: self.busy_nanos.saturating_sub(prev.busy_nanos),
            items: self.items.saturating_sub(prev.items),
        }
    }
}

/// The side that consumes decoded bytes: mapping, aligning, counting.
///
/// Implementors resize a worker pool and report how much work it has done.
pub trait Consumer: Send + Sync {
    /// Ask for `n` worker threads. May take effect gradually.
    ///
    /// Returning success acknowledges the request, not its completion. The
    /// broker observes [`Self::live_threads`] before spending the released
    /// capacity on the producer.
    fn set_threads(&self, n: usize) -> Result<(), ResizeError>;

    /// Workers **actually running now**, which may lag what was asked for.
    ///
    /// Diagnostic only, and deliberately not "the target": the broker tracks
    /// what it asked for itself. The gap between the two is what you want in a
    /// log when a split looks wrong, because it separates "the broker decided
    /// badly" from "the broker decided well and the pool has not caught up".
    /// Reporting the target here would hide exactly that distinction.
    fn live_threads(&self) -> usize;

    /// Cumulative work. See [`Work`] for the two contracts that matter.
    ///
    /// Publish more often than once per batch. A counter updated only at batch
    /// boundaries cannot be sampled faster than a batch completes, and a window
    /// containing no completed batch reads as zero work — indistinguishable
    /// from a totally idle consumer, and reporting maximum starvation for one
    /// running flat out. [`BusyMeter`] handles this.
    fn work(&self) -> Work;

    /// Switch instrumentation between decision-quality and settled cadence.
    ///
    /// Most consumers can leave this as a no-op. Consumers for which publishing
    /// progress is itself measurable work may use a finer cadence only while a
    /// decision is open and a coarser cadence in steady state.
    fn set_measurement_mode(&self, _mode: ProducerMeasurementMode) {}
}

/// The side that produces decoded bytes: a decompressor or pool of them.
pub trait Producer: Send + Sync {
    /// Ask for an aggregate concurrency limit of `n`.
    ///
    /// Returning success acknowledges the request, not the completion of work
    /// already admitted under the previous limit.
    fn set_limit(&self, n: usize) -> Result<(), ResizeError>;

    /// The current limit.
    fn limit(&self) -> usize;

    /// Slots executing producer work now.
    ///
    /// This is the acknowledgement signal for a shrink. It must include work
    /// admitted before a lower limit was requested and exclude coordinator or
    /// I/O threads that do not consume controlled execution slots.
    fn active_slots(&self) -> usize;

    /// Current queued output in the same unit as [`Work::items`], when known.
    ///
    /// The controller rejects cost windows in which this buffer changes
    /// materially, because the two stages then processed different amounts of
    /// logical work and their busy-time ratio is not a per-item cost ratio.
    fn buffered_items(&self) -> Option<u64> {
        None
    }

    /// Live auxiliary threads outside the controlled execution-slot budget.
    fn auxiliary_threads(&self) -> usize {
        0
    }

    /// Select how aggressively an adapter should observe producer busy time.
    ///
    /// Calibration is requested only while the broker is gathering or
    /// ratifying cost evidence. Monitoring is requested while resizing,
    /// blacking out contaminated windows, or sitting at a settled split. The
    /// default is a no-op because producers with native cumulative counters do
    /// not need a sampling policy.
    fn set_measurement_mode(&self, _mode: ProducerMeasurementMode) {}

    /// Suggest an observation interval while the broker is settled.
    ///
    /// Producers with native cumulative counters may ignore this. An adapter
    /// that reconstructs busy time by polling should reduce its own observation
    /// rate when steady-state probes are far apart. It may clamp the request to
    /// the minimum cadence needed for a useful estimate.
    fn set_monitoring_interval(&self, _interval: Duration) {}

    /// Optional diagnostics for a sampled producer measurement.
    fn measurement_stats(&self) -> Option<ProducerMeasurementStats> {
        None
    }

    /// Stop any private measurement adapter and return its final diagnostics.
    ///
    /// The controller calls this exactly once when it exits. Native cumulative
    /// counters may use the default snapshot implementation. Sampled adapters
    /// should join their sampler here so lifetime CPU accounting is complete
    /// before the report is serialized.
    fn finish_measurement(&self) -> Option<ProducerMeasurementStats> {
        self.measurement_stats()
    }

    /// Why the producer is or is not keeping up.
    ///
    /// Used only as a directional guard; it never selects a target. Source-bound
    /// or inelastic evidence may cap growth. A live [`ProducerPressure::Starved`]
    /// signal may veto a model-requested shrink while the current limit is
    /// demonstrably occupied, which protects against allocation-dependent cost
    /// under-reporting. That exception can retain the current split but cannot
    /// grow it or choose a larger one. Pressure signals are deliberately not
    /// used to size the split because they saturate: they cannot distinguish
    /// "just enough" from "far too much".
    fn pressure(&self) -> ProducerPressure;

    /// Cumulative work. See [`Work`].
    fn work(&self) -> Work;
}

/// Requested fidelity for an adapter that reconstructs cumulative busy time
/// from instantaneous observations.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum ProducerMeasurementMode {
    /// Short, high-resolution observation while a cost decision is open.
    #[default]
    Calibration,
    /// Low-frequency observation sufficient to detect a possible regime change.
    Monitoring,
}

/// What the broker does after it has established a stable split.
///
/// All variants share warm-up, cost measurement, guarded model resizing,
/// blackout, ratification, and the first clean steady-state check. Responsive
/// may additionally run an opt-in opening bracket. The two freeze variants make
/// the calibration/overhead tradeoff explicit rather than overloading one
/// ambiguous "freeze" setting.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum SteadyStatePolicy {
    /// Continue low-frequency observation and re-open the controller after a
    /// persistent workload change. Also permits an opt-in opening bracket.
    #[default]
    Responsive,
    /// Use the model-only path, then stop the controller and release recurring
    /// measurement adapters after the first stable split is established.
    ///
    /// This makes the recurring broker cost zero for the rest of a stable run,
    /// but intentionally gives up adaptation to later workload regimes. Use it
    /// only when that tradeoff is known to fit the application.
    FreezeAfterConvergence,
    /// Complete the same opt-in opening calibration as
    /// [`Responsive`](Self::Responsive), then stop the controller and release
    /// recurring measurement adapters after the first stable split.
    ///
    /// This pays a bounded startup calibration cost in exchange for zero
    /// recurring broker work. It cannot react to a later workload regime.
    FreezeAfterFullCalibration,
}

/// How the broker treats the caller's initial split.
///
/// The default trusts the allocation-independent cost model. Applications
/// with measured allocation-dependent scaling can instead ask the broker to
/// confirm a disagreement between that model and the opening. The experiment
/// is startup-only, bounded in points and wall time, and does no work when the
/// first stable model answer agrees with the opening.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum OpeningPolicy {
    /// Treat the opening as an initial hint and follow ordinary model moves.
    #[default]
    Fixed,
    /// Confirm a model/opening disagreement with a bounded local bracket. The
    /// model answer is rejected only when a statistically separated, material
    /// throughput gain over the best measured opening/model rate is proven at
    /// an adjacent point. An inconclusive or regressed comparison may trigger
    /// that search unless independent empirical-cap evidence confirms the
    /// model target. The opening is never a fallback; deadline-limited evidence
    /// also retains the model.
    Bracket(OpeningBracketConfig),
}

/// Bounds for [`OpeningPolicy::Bracket`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OpeningBracketConfig {
    /// Maximum allocations measured by the bracket, including the model's
    /// answer but excluding the already-measured opening.
    pub max_points: usize,
    /// Throughput evidence gathered at each candidate allocation.
    pub horizon: Duration,
    /// Hard wall-time budget for the startup experiment.
    pub total_budget: Duration,
}

impl Default for OpeningBracketConfig {
    fn default() -> Self {
        Self {
            max_points: 3,
            horizon: Duration::from_millis(200),
            // Three 200 ms evidence horizons plus their ordinary resize
            // blackouts fit comfortably below four seconds on the measured
            // scATAC path, leaving margin for scheduling and acknowledgements.
            total_budget: Duration::from_secs(4),
        }
    }
}

/// Bounded diagnostics for producer-side busy-time measurement.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Serialize)]
pub struct ProducerMeasurementStats {
    /// Cumulative executing-worker time used by the controller. This is exact
    /// when `accounted_busy_nanos` is present and sampled otherwise.
    pub busy_nanos: u64,
    /// Exact event-time executing-worker integral, when the producer supplies
    /// one. Production adapters may use this as `busy_nanos`; retaining it here
    /// makes the signal source and sampled-fallback diagnostics explicit.
    pub accounted_busy_nanos: Option<u64>,
    /// CPU time from completed producer worker-thread lifetimes, when the
    /// adapter supports optional component accounting.
    pub completed_worker_cpu_nanos: Option<u64>,
    /// CPU time from completed producer auxiliary-thread lifetimes, when the
    /// adapter supports optional component accounting.
    pub completed_auxiliary_cpu_nanos: Option<u64>,
    /// Thread CPU-clock read failures, when component accounting is enabled.
    pub cpu_accounting_failures: Option<usize>,
    /// Lifetime CPU time of the measurement sampler itself. This requires only
    /// one thread-clock read at sampler start and one at sampler exit; no clock
    /// reads are added to its observation loop.
    pub sampler_cpu_nanos: Option<u64>,
    /// Sampler thread-clock reads that failed.
    pub sampler_cpu_accounting_failures: usize,
    pub calibration_samples: u64,
    pub monitoring_samples: u64,
    pub mode_changes: u64,
    pub final_mode: ProducerMeasurementMode,
    /// Wall time spent obtaining instantaneous observations, excluding sleeps.
    pub observation_nanos: u64,
    pub calibration_interval_micros: u64,
    pub monitoring_interval_micros: u64,
}

/// A two-read lifetime CPU timer for low-overhead administrative accounting.
///
/// It must be started and read on the same thread. Constructing it reads that
/// thread's CPU clock once, and [`Self::elapsed`] reads it once more. It does
/// not install periodic instrumentation or affect application hot paths.
#[doc(hidden)]
pub struct ThreadCpuTimer(Option<cpu_time::ThreadTime>);

impl ThreadCpuTimer {
    pub fn start() -> Self {
        Self(cpu_time::ThreadTime::try_now().ok())
    }

    pub fn elapsed(&self) -> Option<Duration> {
        self.0.as_ref()?.try_elapsed().ok()
    }
}

/// Whether more producer concurrency would buy anything.
///
/// # This guards direction; it does not size the split
///
/// [`Self::SourceBound`] and [`Self::Inelastic`] may cap how far the producer
/// grows. [`Self::Starved`] has one bounded exception: it may retain the current
/// limit when the cost model asks to shrink a producer that is demonstrably
/// using its present capacity. It cannot choose or increase the target, so the
/// cost-share model remains the only sizing mechanism. [`Self::Satisfied`] is
/// diagnostic. This one-sided restriction is deliberate and hard-won.
///
/// Every variant here is a *saturating* signal, and an earlier version of this
/// crate sized the split from one. It does not work, and the reason is
/// structural rather than a matter of tuning.
///
/// Let `starve(d)` be the fraction of consumer time blocked on an empty buffer.
/// It falls as `d` rises and reaches zero at the balance point — and then
/// **stays** zero for every larger `d`. So `{d : starve(d) = 0}` is a half-line,
/// not a point, and a controller driving `starve -> 0` has no restoring force
/// anywhere on it. Any bias walks it to the boundary while the error stays
/// genuinely near zero the whole way. Measured: it ran to 47 producer / 17
/// consumer where 32/32 was optimal, 44% slower, "converging" throughout.
///
/// DS2 (Kalavri et al., OSDI'18, section 2) classifies exactly that controller
/// -- threshold policies over utilisation, queue depth or observed rates -- and
/// gives the consequence as "incorrect provisioning, oscillations, and long
/// convergence times".
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProducerPressure {
    /// Work is queued behind the limit: more slots would be used.
    ///
    /// This variant carried a `wanted` count while the split was grown one
    /// decision at a time. Nothing reads it now — the split is solved, not
    /// stepped — and a demand estimate is a poor thing to keep alive with no
    /// consumer to keep it honest.
    ///
    /// May veto removal of an occupied current slot when allocation-dependent
    /// producer work makes the low-allocation cost estimate incomplete. The
    /// controller counts that event and requires persistent evidence; this
    /// variant never requests growth.
    Starved,
    /// Keeping up with the consumer. More slots would idle.
    Satisfied,
    /// Idle, but *not* because of the limit — the source is the constraint.
    ///
    /// Distinct from [`Self::Satisfied`] because the consumer may simultaneously
    /// be starving, and that combination is the one case where a starved
    /// consumer must **not** cause the producer to grow: the threads would sit
    /// waiting on storage. Without this variant a broker spends its whole budget
    /// on decode for an I/O-bound run and gains nothing.
    SourceBound,
    /// Cannot use more concurrency at any price.
    ///
    /// For a stream the decompressor can only process serially. Such an input
    /// reports total consumer starvation forever and can never be helped, so
    /// without this the broker would spend every thread chasing it.
    Inelastic,
}

/// Invalid [`ThreadBroker`] configuration.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BrokerConfigError(&'static str);

impl std::fmt::Display for BrokerConfigError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "invalid thread-broker configuration: {}", self.0)
    }
}

impl std::error::Error for BrokerConfigError {}

/// A pool refused a resize request.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResizeError(Box<str>);

impl ResizeError {
    pub fn new(message: impl Into<String>) -> Self {
        Self(message.into().into_boxed_str())
    }
}

impl std::fmt::Display for ResizeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for ResizeError {}

/// Which pipeline side was involved in a broker operation or failure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ResizeSide {
    Consumer,
    Producer,
}

/// Why a broker could not start or complete safely.
#[derive(Debug)]
pub enum BrokerErrorKind {
    ResizeRefused {
        side: ResizeSide,
        requested: usize,
        source: ResizeError,
    },
    ResizeTimedOut {
        side: ResizeSide,
        target: usize,
        observed: usize,
        timeout: Duration,
    },
    /// An adapter violated [`Work`]'s cumulative-counter contract.
    WorkCountersRegressed {
        side: ResizeSide,
        previous: Work,
        observed: Work,
    },
    ThreadSpawn(std::io::Error),
    ThreadPanicked,
}

/// A broker failure, optionally carrying the report accumulated before it.
#[derive(Debug)]
pub struct BrokerError {
    kind: BrokerErrorKind,
    report: Option<Box<BrokerReport>>,
}

impl BrokerError {
    pub(crate) fn new(kind: BrokerErrorKind) -> Self {
        Self { kind, report: None }
    }

    pub(crate) fn with_report(mut self, report: BrokerReport) -> Self {
        self.report = Some(Box::new(report));
        self
    }

    pub fn kind(&self) -> &BrokerErrorKind {
        &self.kind
    }

    pub fn report(&self) -> Option<&BrokerReport> {
        self.report.as_deref()
    }
}

impl std::fmt::Display for BrokerError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match &self.kind {
            BrokerErrorKind::ResizeRefused {
                side,
                requested,
                source,
            } => write!(
                f,
                "thread broker {side:?} resize to {requested} was refused: {source}"
            ),
            BrokerErrorKind::ResizeTimedOut {
                side,
                target,
                observed,
                timeout,
            } => write!(
                f,
                "thread broker timed out after {timeout:?} waiting for {side:?} to shrink to {target} (observed {observed})"
            ),
            BrokerErrorKind::WorkCountersRegressed {
                side,
                previous,
                observed,
            } => write!(
                f,
                "thread broker observed decreasing {side:?} work counters: {previous:?} -> {observed:?}"
            ),
            BrokerErrorKind::ThreadSpawn(source) => {
                write!(f, "could not spawn thread-broker sampler: {source}")
            }
            BrokerErrorKind::ThreadPanicked => f.write_str("thread-broker sampler panicked"),
        }
    }
}

impl std::error::Error for BrokerError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match &self.kind {
            BrokerErrorKind::ResizeRefused { source, .. } => Some(source),
            BrokerErrorKind::ThreadSpawn(source) => Some(source),
            _ => None,
        }
    }
}

/// Tuning, all of it overridable.
///
/// Every default below was arrived at by fixing an observed failure, so they
/// carry real information — but none is a law, and a caller with different
/// timings should move them.
#[derive(Debug, Clone, Copy)]
pub struct BrokerConfig {
    /// How often to sample and decide.
    ///
    /// 100 ms by default, with a 400 ms warm-up over it. An earlier pairing of
    /// 250 ms and 150 ms was *self-inconsistent*: it failed this type's own
    /// `warmup >= 2 * sample_interval` rule, so every broker built from the
    /// defaults failed to build. Nothing caught it, because every unit test
    /// passed an explicit config — the defaults were the one combination never
    /// exercised. There is now a test that builds from them.
    ///
    /// 250 ms was chosen to match `rapidgzip`'s worker-retirement hysteresis, on
    /// the grounds that instantaneous fields sampled at that period alias
    /// against their own duty cycle. That reasoning applies to a producer that
    /// reports instantaneous worker counts; the signals used here are an
    /// integrated busy counter and a demand estimate, neither of which aliases.
    pub sample_interval: Duration,
    /// How often a responsive broker probes after convergence.
    ///
    /// `None` (the default) uses [`Self::sample_interval`], preserving the fully
    /// responsive behavior. A longer interval changes only settled monitoring;
    /// warm-up, calibration, resizing, blackout, and ratification continue at
    /// `sample_interval`.
    pub steady_probe_interval: Option<Duration>,
    /// Time discarded before the first decision.
    ///
    /// **Load-bearing, not merely tuned.** Every workload measured reported ~99%
    /// consumer starvation in its first 75 ms, because the consumer's threads
    /// are still starting and nothing has been processed yet. With no warm-up,
    /// *every* input looks producer-bound and the broker's first move is always
    /// wrong.
    ///
    /// Must span at least two sampling intervals, or the first decision rests on
    /// a single window.
    pub warmup: Duration,
    /// Recent windows averaged into the rate estimates.
    ///
    /// Deciding on a single window makes the broker a noise amplifier: a real
    /// consumer throws the occasional idle sample from a batch boundary or a
    /// page fault. Averaging is what makes the *signal* steady, rather than
    /// trying to out-tune the noise downstream of it.
    ///
    /// One disables smoothing.
    pub smoothing_windows: usize,
    /// Minimum distance, in threads, between the model's answer and the current
    /// split before the broker will act.
    ///
    /// Moves are not free — each one rebuilds per-thread state on both sides,
    /// including a multi-MB record buffer per new consumer worker — so the
    /// broker should not spend one on a difference it cannot be sure about.
    ///
    /// # This value is not yet derived, and should be
    ///
    /// The measurements it was originally justified from are void. They were
    /// taken through a path that did not enforce the thread budget, which made a
    /// sharply peaked throughput surface look like a wide plateau, and the
    /// conclusion drawn — "precision near the answer is worth almost nothing" —
    /// does not survive re-measurement.
    ///
    /// What the corrected sweeps show, on a decode-bound single-cell workload:
    ///
    /// | budget | optimum | surface |
    /// |---:|---:|---|
    /// | 32 | 12 slots | 8 slots +29%, 16 +13%, 24 +115% |
    /// | 64 | 22 slots | 16 slots +17%, 28 +5%, 40 +34% |
    ///
    /// Two threads is defensible against both — the broker landed 5.0% and 3.4%
    /// off the best swept point — but that is an observation, not a derivation,
    /// and the two budgets differ enough in tolerance that a single constant is
    /// unlikely to be right for both. It matters most at small budgets, where
    /// two threads is a large fraction of the answer.
    pub deadband_threads: usize,
    /// Samples discarded after a move before the model is trusted again.
    ///
    /// Two distinct effects need to decay. The pools converge toward their new
    /// size over some interval, so a window spanning a move measures a split
    /// that no longer exists; and the move itself does work — a new consumer
    /// worker builds its thread state and allocates its record set — which
    /// inflates busy time without moving any items, biasing the very ratio the
    /// decision was made from.
    ///
    /// Must be at least `smoothing_windows`, or the smoothed estimate still
    /// contains pre-move windows when the broker next decides.
    pub blackout_samples: usize,
    /// Windows of throughput evidence gathered before a move is ratified.
    /// The default ten-window block spans one second at the default 100 ms
    /// cadence. Applications should validate its power against stationary
    /// windows from their real workload; a deliberately changing or pausing
    /// workload is not a statistically identifiable ratification population.
    pub ratify_samples: usize,
    /// Whether a disagreement between the opening and the first stable model
    /// answer should be confirmed by a bounded startup bracket.
    pub opening_policy: OpeningPolicy,
    /// Relative throughput loss that counts as a regression rather than noise.
    ///
    /// Only a *regression* reverts a move. Absence of improvement does not,
    /// because the surface is flat near the optimum: "no better" is the expected
    /// reading for a correct move, and treating it as failure would revert the
    /// broker to wherever it happened to start. Flink's equivalent check,
    /// `detectIneffectiveScaleUp`, compares against a model-*predicted* increase
    /// for this reason, and ships disabled by default.
    pub regression_tolerance: f64,
    /// How far the model's answer must drift from the settled split, in threads,
    /// before the broker re-opens a decision it had closed.
    ///
    /// Settling must not mean going deaf: a file ending or a different sample
    /// beginning genuinely changes the costs. But adapting continuously spends
    /// moves on noise forever, and each move rebuilds per-thread state on both
    /// sides. So the broker settles and then sleeps until surprised. This value
    /// is only an absolute floor: observed cost-share uncertainty can widen the
    /// effective threshold, and persistent evidence is required before the
    /// decision re-opens. Keeping the floor at one preserves sensitivity to a
    /// one-slot regime change on small budgets without making noisy large
    /// budgets unstable.
    ///
    /// # Why drift rather than a throughput CUSUM
    ///
    /// The obvious alternative is to watch throughput and re-arm when it falls —
    /// a one-sided CUSUM, and the form the autoscaling literature usually
    /// reaches for. It is genuinely one-sided, and that turns out to matter
    /// here: when a consumer suddenly gets *cheaper*, throughput does not fall
    /// at all — the pipeline simply stays pinned at the old producer-limited
    /// rate while a much better split goes unclaimed. A CUSUM sleeps through it.
    /// The solved target moves immediately, because it is computed from the cost
    /// ratio rather than from the outcome, so watching the model detects change
    /// in both directions with one mechanism instead of two.
    pub resurvey_distance: usize,
    /// Wall-clock time model drift must persist before a settled split is
    /// re-opened.
    ///
    /// Distance alone is not enough, and the reason is that the noise in the
    /// solved target scales with the budget while a fixed distance does not.
    /// Measured on a 116-second Flex run at `-t 64`: with distance as the only
    /// gate the broker re-opened 18 times and never reported convergence, while
    /// the same build at `-t 32` settled once in 1.4 s and never moved again —
    /// identical workload, identical noise in the *share*, twice the noise in
    /// threads.
    ///
    /// Requiring duration-based persistence separates the two cleanly, because
    /// it keys on something that actually distinguishes them: sampling noise
    /// does not survive a second of consecutive windows, and a real regime
    /// change does. It also costs nothing when there is no drift.
    pub resurvey_persistence: Duration,
    /// Maximum time to wait for the shrinking side to release capacity before
    /// aborting the broker. The other side is never grown before this
    /// acknowledgement arrives.
    pub resize_timeout: Duration,
    /// Absolute queue drift tolerated in one cost window.
    pub max_buffer_drift_items: u64,
    /// Queue drift tolerated as a fraction of producer progress in the window.
    pub max_buffer_drift_fraction: f64,
    /// Wall-clock history retained for producer scalability caps.
    pub cap_history: Duration,
    /// Continuous corroborating evidence required before slack or source
    /// saturation is allowed to constrain a move.
    pub cap_persistence: Duration,
    /// Maximum number of requested/actual samples retained in a report.
    pub trace_capacity: usize,
    /// The consumer never drops below this.
    pub min_consumer_threads: usize,
    /// Producer concurrency never drops below this, in aggregate.
    ///
    /// One is right for a producer whose limit is a pure runtime throttle, and
    /// that is worth checking rather than assuming. `rapidgzip` used to read
    /// this limit while *choosing a backend*, so anything below two committed a
    /// decoder to sequential decoding irreversibly — a caller on such a version,
    /// or on any producer that couples its limit to a one-time decision, needs
    /// this at two or more.
    ///
    /// Measured on the shared-pool decoder: an aggregate limit of 1 leaves the
    /// path at `DenseMembers` and still reports a demand of 12, so admission is
    /// genuinely decoupled and the floor can be a single slot. Note this is
    /// **aggregate, not per input** — scaling it by input count would reserve
    /// half a 32-thread budget on an 8-file run that needs no decode at all.
    pub min_producer_slots: usize,
}

impl Default for BrokerConfig {
    fn default() -> Self {
        Self {
            sample_interval: Duration::from_millis(100),
            steady_probe_interval: None,
            warmup: Duration::from_millis(400),
            smoothing_windows: 3,
            deadband_threads: 1,
            blackout_samples: 4,
            ratify_samples: 10,
            opening_policy: OpeningPolicy::Fixed,
            regression_tolerance: 0.05,
            resurvey_distance: 1,
            resurvey_persistence: Duration::from_millis(800),
            resize_timeout: Duration::from_secs(2),
            max_buffer_drift_items: 1 << 20,
            max_buffer_drift_fraction: 0.05,
            cap_history: Duration::from_millis(800),
            cap_persistence: Duration::from_millis(300),
            trace_capacity: 256,
            min_consumer_threads: 1,
            min_producer_slots: 1,
        }
    }
}

/// A configured, not-yet-running broker.
pub struct ThreadBroker<C: Consumer, P: Producer> {
    pub(crate) consumer: Arc<C>,
    pub(crate) producer: Arc<P>,
    pub(crate) budget: usize,
    pub(crate) config: BrokerConfig,
    pub(crate) initial_producer_slots: Option<usize>,
    pub(crate) steady_state_policy: SteadyStatePolicy,
}

/// Builder for [`ThreadBroker`].
pub struct ThreadBrokerBuilder<C: Consumer, P: Producer> {
    consumer: Arc<C>,
    producer: Arc<P>,
    budget: Option<usize>,
    config: BrokerConfig,
    initial_producer_slots: Option<usize>,
    steady_state_policy: SteadyStatePolicy,
}

impl<C: Consumer, P: Producer> ThreadBroker<C, P> {
    /// Start configuring a broker over these two sides, from the defaults.
    pub fn builder(consumer: C, producer: P) -> ThreadBrokerBuilder<C, P> {
        Self::builder_with(consumer, producer, BrokerConfig::default())
    }

    /// As [`Self::builder`], but starting from `config` rather than the
    /// defaults.
    ///
    /// A whole `BrokerConfig` is taken here, at construction, rather than
    /// offered as a `.config()` setter. A setter would silently discard every
    /// individual override called before it — which it did, wiping a floor of 4
    /// back to 1 and letting a test watch the broker breach a limit it had
    /// explicitly been given. Taking it up front leaves no ordering to get
    /// wrong.
    pub fn builder_with(
        consumer: C,
        producer: P,
        config: BrokerConfig,
    ) -> ThreadBrokerBuilder<C, P> {
        ThreadBrokerBuilder {
            consumer: Arc::new(consumer),
            producer: Arc::new(producer),
            budget: None,
            config,
            initial_producer_slots: None,
            steady_state_policy: SteadyStatePolicy::default(),
        }
    }
}

impl<C: Consumer, P: Producer> ThreadBrokerBuilder<C, P> {
    /// Total threads to divide. Required.
    pub fn budget(mut self, threads: usize) -> Self {
        self.budget = Some(threads);
        self
    }

    /// Start the producer at this many slots instead of the coarse default.
    ///
    /// The default start is deliberately low, because growth is cheap and many
    /// workloads never need decode threads at all. A caller that already knows
    /// better — from a previous run, or from an explicit user request — can skip
    /// the climb. Still clamped to the per-input floor and the budget.
    pub fn initial_producer_slots(mut self, slots: usize) -> Self {
        self.initial_producer_slots = Some(slots);
        self
    }

    /// Choose whether a settled broker keeps watching for regime changes.
    ///
    /// See [`SteadyStatePolicy`] for the responsiveness/overhead tradeoff.
    pub fn steady_state_policy(mut self, policy: SteadyStatePolicy) -> Self {
        self.steady_state_policy = policy;
        self
    }

    /// See [`BrokerConfig::sample_interval`].
    pub fn sample_interval(mut self, interval: Duration) -> Self {
        self.config.sample_interval = interval;
        self
    }

    /// Set the interval between responsive steady-state probes.
    ///
    /// This does not change the high-resolution convergence cadence. The
    /// producer is asked to use roughly four monitoring observations per probe,
    /// subject to any accuracy floor imposed by its adapter.
    pub fn steady_probe_interval(mut self, interval: Duration) -> Self {
        self.config.steady_probe_interval = Some(interval);
        self
    }

    /// See [`BrokerConfig::warmup`].
    pub fn warmup(mut self, warmup: Duration) -> Self {
        self.config.warmup = warmup;
        self
    }

    /// See [`BrokerConfig::smoothing_windows`].
    pub fn smoothing_windows(mut self, windows: usize) -> Self {
        self.config.smoothing_windows = windows;
        self
    }

    /// See [`BrokerConfig::deadband_threads`].
    pub fn deadband_threads(mut self, threads: usize) -> Self {
        self.config.deadband_threads = threads;
        self
    }

    /// See [`BrokerConfig::blackout_samples`].
    pub fn blackout_samples(mut self, samples: usize) -> Self {
        self.config.blackout_samples = samples;
        self
    }

    /// See [`BrokerConfig::ratify_samples`].
    pub fn ratify_samples(mut self, samples: usize) -> Self {
        self.config.ratify_samples = samples;
        self
    }

    /// See [`BrokerConfig::opening_policy`].
    pub fn opening_policy(mut self, policy: OpeningPolicy) -> Self {
        self.config.opening_policy = policy;
        self
    }

    /// See [`BrokerConfig::regression_tolerance`].
    pub fn regression_tolerance(mut self, fraction: f64) -> Self {
        self.config.regression_tolerance = fraction;
        self
    }

    /// See [`BrokerConfig::resurvey_distance`].
    pub fn resurvey_distance(mut self, threads: usize) -> Self {
        self.config.resurvey_distance = threads;
        self
    }

    /// See [`BrokerConfig::resurvey_persistence`].
    pub fn resurvey_persistence(mut self, duration: Duration) -> Self {
        self.config.resurvey_persistence = duration;
        self
    }

    /// See [`BrokerConfig::resize_timeout`].
    pub fn resize_timeout(mut self, timeout: Duration) -> Self {
        self.config.resize_timeout = timeout;
        self
    }

    /// See [`BrokerConfig::min_consumer_threads`].
    pub fn min_consumer_threads(mut self, threads: usize) -> Self {
        self.config.min_consumer_threads = threads.max(1);
        self
    }

    /// See [`BrokerConfig::min_producer_slots`].
    pub fn min_producer_slots(mut self, slots: usize) -> Self {
        self.config.min_producer_slots = slots.max(1);
        self
    }

    /// Validate and finish.
    ///
    /// # Errors
    ///
    /// Rejects a budget too small to split, a zero sampling interval, and a
    /// warm-up shorter than two sampling intervals — the last because a broker
    /// that decides on one sample is deciding on the startup transient, which
    /// is the failure the warm-up exists to prevent.
    pub fn build(self) -> Result<ThreadBroker<C, P>, BrokerConfigError> {
        let budget = self.budget.ok_or(BrokerConfigError("budget is required"))?;
        let c = &self.config;
        if budget < 2 {
            return Err(BrokerConfigError(
                "budget must be at least 2 to divide between two sides",
            ));
        }
        if c.sample_interval.is_zero() {
            return Err(BrokerConfigError("sample_interval must be non-zero"));
        }
        if c.steady_probe_interval
            .is_some_and(|interval| interval.is_zero())
        {
            return Err(BrokerConfigError(
                "steady_probe_interval must be non-zero when configured",
            ));
        }
        if c.warmup < c.sample_interval * 2 {
            return Err(BrokerConfigError(
                "warmup must span at least two sample intervals, or the first \
                 decision is made on the startup transient",
            ));
        }
        if !(0.0..1.0).contains(&c.regression_tolerance) {
            return Err(BrokerConfigError(
                "regression_tolerance must be in 0.0..1.0",
            ));
        }
        if c.smoothing_windows == 0 {
            return Err(BrokerConfigError("smoothing_windows must be non-zero"));
        }
        // A blackout shorter than the smoothing window does not do its job: the
        // smoothed estimate the next decision reads would still be carrying
        // windows from before the move, so the broker would partly re-decide on
        // the split it has already left.
        if c.blackout_samples < c.smoothing_windows {
            return Err(BrokerConfigError(
                "blackout_samples must be at least smoothing_windows, or the \
                 next decision still averages in pre-move windows",
            ));
        }
        if c.ratify_samples == 0 {
            return Err(BrokerConfigError("ratify_samples must be non-zero"));
        }
        if let OpeningPolicy::Bracket(bracket) = c.opening_policy {
            if bracket.max_points == 0 {
                return Err(BrokerConfigError(
                    "opening bracket max_points must be non-zero",
                ));
            }
            if bracket.horizon.is_zero() {
                return Err(BrokerConfigError(
                    "opening bracket horizon must be non-zero",
                ));
            }
            if bracket.total_budget.is_zero() || bracket.total_budget < bracket.horizon {
                return Err(BrokerConfigError(
                    "opening bracket total_budget must be at least its horizon",
                ));
            }
        }
        // Re-opening is made harder by persistence. Its distance may equal the
        // movement deadband so a one-slot regime change remains observable on
        // small budgets.
        if c.resurvey_persistence.is_zero() {
            return Err(BrokerConfigError("resurvey_persistence must be non-zero"));
        }
        if c.resize_timeout.is_zero() {
            return Err(BrokerConfigError("resize_timeout must be non-zero"));
        }
        if !(0.0..=1.0).contains(&c.max_buffer_drift_fraction) {
            return Err(BrokerConfigError(
                "max_buffer_drift_fraction must be in 0.0..=1.0",
            ));
        }
        if c.cap_history.is_zero()
            || c.cap_persistence.is_zero()
            || c.cap_persistence > c.cap_history
        {
            return Err(BrokerConfigError(
                "cap durations must be non-zero and cap_persistence must not exceed cap_history",
            ));
        }
        if c.trace_capacity == 0 {
            return Err(BrokerConfigError("trace_capacity must be non-zero"));
        }
        if c.resurvey_distance < c.deadband_threads {
            return Err(BrokerConfigError(
                "resurvey_distance must be at least deadband_threads",
            ));
        }
        if c.min_consumer_threads == 0 || c.min_producer_slots == 0 {
            return Err(BrokerConfigError(
                "min_consumer_threads and min_producer_slots must be non-zero",
            ));
        }
        if c.min_consumer_threads + c.min_producer_slots > budget {
            return Err(BrokerConfigError(
                "the two floors do not fit inside the budget",
            ));
        }
        Ok(ThreadBroker {
            consumer: self.consumer,
            producer: self.producer,
            budget,
            config: self.config,
            initial_producer_slots: self.initial_producer_slots,
            steady_state_policy: self.steady_state_policy,
        })
    }
}

/// Shared stop flag and wakeup, so even a very sparse steady probe remains
/// immediately interruptible at application shutdown.
pub(crate) struct Stop {
    value: AtomicBool,
    mutex: std::sync::Mutex<()>,
    wake: std::sync::Condvar,
}

impl Stop {
    pub(crate) fn new() -> Arc<Self> {
        Arc::new(Self {
            value: AtomicBool::new(false),
            mutex: std::sync::Mutex::new(()),
            wake: std::sync::Condvar::new(),
        })
    }
    pub(crate) fn is_set(&self) -> bool {
        self.value.load(Ordering::Acquire)
    }
    pub(crate) fn set(&self) {
        self.value.store(true, Ordering::Release);
        self.wake.notify_all();
    }
    pub(crate) fn wait_timeout(&self, timeout: Duration) -> bool {
        if self.is_set() {
            return true;
        }
        let Ok(guard) = self.mutex.lock() else {
            return self.is_set();
        };
        if self.is_set() {
            return true;
        }
        let _wait = self.wake.wait_timeout(guard, timeout);
        self.is_set()
    }
}

/// A window of consumer activity, normalised against the capacity available.
pub(crate) fn idle_fraction(busy_delta: u64, threads: usize, window: Duration) -> f64 {
    let capacity = (threads as u128).saturating_mul(window.as_nanos());
    if capacity == 0 {
        return 0.0;
    }
    let used = (busy_delta as f64 / capacity as f64).clamp(0.0, 1.0);
    1.0 - used
}

pub(crate) fn now() -> Instant {
    Instant::now()
}

// ---------------------------------------------------------------------------
// Helpers for implementing `Consumer::busy_nanos`
// ---------------------------------------------------------------------------

/// A shared counter of time spent doing real work.
///
/// Every [`Consumer`] needs one, and the obvious implementation of it is wrong
/// in a way that is hard to notice, so it is provided rather than left to each
/// caller. See [`BusyMeter::timer`].
#[derive(Debug, Default)]
pub struct BusyMeter {
    nanos: std::sync::atomic::AtomicU64,
    items: std::sync::atomic::AtomicU64,
}

impl BusyMeter {
    pub fn new() -> Self {
        Self::default()
    }

    /// Everything recorded so far, across every worker.
    pub fn work(&self) -> Work {
        Work {
            busy_nanos: self.nanos.load(Ordering::Relaxed),
            items: self.items.load(Ordering::Relaxed),
        }
    }

    /// Total nanoseconds recorded, across every worker.
    pub fn nanos(&self) -> u64 {
        self.nanos.load(Ordering::Relaxed)
    }

    /// Start timing a batch.
    ///
    /// The returned timer publishes periodically *within* the batch, not only
    /// at its end. That distinction is the whole reason this type exists.
    ///
    /// Publishing once per batch cannot be sampled faster than a batch
    /// completes. With a 16384-record batch at a realistic 5 us per record, a
    /// batch takes ~80 ms against a 250 ms sampling window — so windows
    /// routinely contain no completed batch, read a delta of zero, and report
    /// **maximum starvation for a consumer running flat out**. Measured
    /// directly before this was fixed: a deliberately expensive kernel reported
    /// 100.0% idle, the exact opposite of the truth, while cheaper kernels
    /// reported less.
    pub fn timer(&self) -> BatchTimer<'_> {
        BatchTimer {
            meter: self,
            single_writer_progress: None,
            single_writer_progress_total: 0,
            mark: Instant::now(),
            items_since_publish: 0,
            time_items_since_publish: 0,
            progress_every: DEFAULT_FLUSH_EVERY,
            time_every: DEFAULT_FLUSH_EVERY,
        }
    }
}

/// Records published every this many items by default.
///
/// Keeps the publish interval far below one sampling window even for slow
/// consumers, while costing one `Instant::now()` pair per 256 items — nothing
/// against a consumer doing real work, and still under a percent for a trivial
/// one.
pub const DEFAULT_FLUSH_EVERY: u64 = 256;

/// Accumulates work time at a granularity the broker can actually resolve.
///
/// Publishes on drop, so an early return or a `?` still records what was done,
/// and so the tail of a batch is never lost.
///
/// # Construct one per batch
///
/// The timer measures wall time between publishes, so anything that happens
/// while it is alive is charged as work. Constructing it *inside* the
/// per-batch callback is what keeps blocking out of the measurement: the
/// callback is entered only once a filled buffer is in hand, so the wait for
/// that buffer falls between two timers rather than inside one. A timer hoisted
/// to thread scope would charge every `fill` to the consumer's busy time — and
/// busy time that includes blocking is exactly the wall-time measurement this
/// crate exists to avoid.
pub struct BatchTimer<'a> {
    meter: &'a BusyMeter,
    single_writer_progress: Option<&'a std::sync::atomic::AtomicU64>,
    single_writer_progress_total: u64,
    mark: Instant,
    items_since_publish: u64,
    time_items_since_publish: u64,
    progress_every: u64,
    time_every: u64,
}

impl<'a> BatchTimer<'a> {
    /// Publish item progress to a single-writer caller-owned counter.
    ///
    /// A consumer with several hot workers can give each worker a cache-padded
    /// counter and sum those counters in [`Consumer::work`]. Busy time remains
    /// in this meter, while fine progress uses a relaxed cumulative store
    /// rather than a contended atomic read-modify-write. No other writer may
    /// update `counter`; it must remain cumulative across successive timers.
    pub fn single_writer_progress_counter(
        mut self,
        counter: &'a std::sync::atomic::AtomicU64,
    ) -> Self {
        self.single_writer_progress_total = counter.load(Ordering::Relaxed);
        self.single_writer_progress = Some(counter);
        self
    }

    /// Publish both progress and busy time every `n` items instead of the default.
    ///
    /// Prefer [`Self::progress_every`] when only the completed-item signal needs
    /// finer resolution. That avoids adding clock reads to a hot record loop.
    pub fn flush_every(mut self, n: u64) -> Self {
        let n = n.max(1);
        self.progress_every = n;
        self.time_every = n;
        self
    }

    /// Publish completed-item progress every `n` items.
    ///
    /// This cadence is independent of busy-time publication: making throughput
    /// visible to short controller windows therefore costs one relaxed counter
    /// update per interval, without forcing an `Instant::now()` at that cadence.
    pub fn progress_every(mut self, n: u64) -> Self {
        self.progress_every = n.max(1);
        self
    }

    /// Publish accumulated busy time every `n` items.
    ///
    /// Most consumers should retain [`DEFAULT_FLUSH_EVERY`]. This separate
    /// control primarily lets a converged freeze policy disable all in-batch
    /// measurement for new callbacks.
    pub fn time_every(mut self, n: u64) -> Self {
        self.time_every = n.max(1);
        self
    }

    /// Call once per item processed.
    ///
    /// The count is only ever compared against itself over time, so any
    /// consistent unit will do; it plays no part in sizing the split.
    #[inline]
    pub fn tick(&mut self) {
        self.items_since_publish += 1;
        self.time_items_since_publish += 1;
        if self.items_since_publish >= self.progress_every {
            self.publish_items();
        }
        if self.time_items_since_publish >= self.time_every {
            self.publish_time();
        }
    }

    /// Count one item when the returned guard leaves scope.
    ///
    /// Use this when a loop has early exits and throughput must mean completed,
    /// rather than started, items. Construct the guard at the top of an
    /// iteration; normal completion, `continue`, and `?` all publish the item
    /// only after that iteration's work has finished.
    #[inline]
    pub fn completed_item(&mut self) -> CompletedItem<'_, 'a> {
        CompletedItem { timer: self }
    }

    #[inline]
    fn publish_items(&mut self) {
        if let Some(counter) = self.single_writer_progress {
            self.single_writer_progress_total = self
                .single_writer_progress_total
                .saturating_add(self.items_since_publish);
            counter.store(self.single_writer_progress_total, Ordering::Relaxed);
        } else {
            self.meter
                .items
                .fetch_add(self.items_since_publish, Ordering::Relaxed);
        }
        self.items_since_publish = 0;
    }

    #[inline]
    fn publish_time(&mut self) {
        let now = Instant::now();
        self.meter
            .nanos
            .fetch_add((now - self.mark).as_nanos() as u64, Ordering::Relaxed);
        self.mark = now;
        self.time_items_since_publish = 0;
    }
}

/// Completion guard returned by [`BatchTimer::completed_item`].
pub struct CompletedItem<'timer, 'meter> {
    timer: &'timer mut BatchTimer<'meter>,
}

impl Drop for CompletedItem<'_, '_> {
    #[inline]
    fn drop(&mut self) {
        self.timer.tick();
    }
}

impl Drop for BatchTimer<'_> {
    fn drop(&mut self) {
        self.publish_items();
        self.publish_time();
    }
}
