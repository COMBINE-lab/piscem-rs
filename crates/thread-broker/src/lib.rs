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
//! # -> Result<(), BrokerConfigError> {
//! let broker = ThreadBroker::builder(mapper, decoders)
//!     .budget(32)
//!     .sample_interval(Duration::from_millis(100))
//!     .build()?;
//!
//! let running = broker.start();
//! // ... run the job ...
//! let report = running.finish();
//! tracing::info!(
//!     "decode settled at {} of 32 threads after {:?}",
//!     report.final_producer_limit,
//!     report.time_to_converge,
//! );
//! # Ok(())
//! # }
//! ```

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::{Duration, Instant};

mod controller;
pub use controller::{BrokerReport, Model, RunningBroker};

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
    fn set_threads(&self, n: usize);

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
}

/// The side that produces decoded bytes: a decompressor or pool of them.
pub trait Producer: Send + Sync {
    /// Ask for an aggregate concurrency limit of `n`.
    fn set_limit(&self, n: usize);

    /// The current limit.
    fn limit(&self) -> usize;

    /// Why the producer is or is not keeping up.
    ///
    /// Used only to *veto* growth — see [`ProducerPressure::SourceBound`] and
    /// [`ProducerPressure::Inelastic`]. It is deliberately not used to size the
    /// split, because pressure signals saturate: they cannot distinguish "just
    /// enough" from "far too much".
    fn pressure(&self) -> ProducerPressure;

    /// Cumulative work. See [`Work`].
    fn work(&self) -> Work;
}

/// Whether more producer concurrency would buy anything.
///
/// # This vetoes growth; it does not size the split
///
/// Only [`Self::SourceBound`] and [`Self::Inelastic`] change a decision, and only
/// by capping how far the producer may grow. The other two are diagnostic. That
/// restriction is deliberate and hard-won.
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
    pub ratify_samples: usize,
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
    /// sides. So the broker settles and then sleeps until surprised, with this
    /// as the surprise threshold. It should be comfortably wider than
    /// [`Self::deadband_threads`], which is the asymmetric hysteresis the
    /// load-balancing literature recommends: cheap to keep a decision, dear to
    /// re-open one.
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
    /// Consecutive samples the drift must persist before a settled split is
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
    /// Requiring persistence separates the two cleanly, because it keys on
    /// something that actually distinguishes them: sampling noise does not
    /// survive a second of consecutive windows, and a real regime change does.
    /// It also costs nothing when there is no drift.
    pub resurvey_samples: usize,
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
            warmup: Duration::from_millis(400),
            smoothing_windows: 3,
            deadband_threads: 2,
            blackout_samples: 4,
            ratify_samples: 4,
            regression_tolerance: 0.05,
            resurvey_distance: 6,
            resurvey_samples: 8,
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
}

/// Builder for [`ThreadBroker`].
pub struct ThreadBrokerBuilder<C: Consumer, P: Producer> {
    consumer: Arc<C>,
    producer: Arc<P>,
    budget: Option<usize>,
    config: BrokerConfig,
    initial_producer_slots: Option<usize>,
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

    /// See [`BrokerConfig::sample_interval`].
    pub fn sample_interval(mut self, interval: Duration) -> Self {
        self.config.sample_interval = interval;
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

    /// See [`BrokerConfig::resurvey_samples`].
    pub fn resurvey_samples(mut self, samples: usize) -> Self {
        self.config.resurvey_samples = samples;
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
        // Re-opening a settled decision must be strictly harder than making one,
        // or the broker can reach a split it immediately wants to leave and
        // oscillate between the two forever.
        if c.resurvey_samples == 0 {
            return Err(BrokerConfigError("resurvey_samples must be non-zero"));
        }
        if c.resurvey_distance <= c.deadband_threads {
            return Err(BrokerConfigError(
                "resurvey_distance must exceed deadband_threads, or a settled \
                 split can immediately re-open itself",
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
        })
    }
}

/// Shared stop flag, so the sampling thread can be told to wind up.
pub(crate) struct Stop(pub(crate) AtomicBool);

impl Stop {
    pub(crate) fn new() -> Arc<Self> {
        Arc::new(Self(AtomicBool::new(false)))
    }
    pub(crate) fn is_set(&self) -> bool {
        self.0.load(Ordering::Relaxed)
    }
    pub(crate) fn set(&self) {
        self.0.store(true, Ordering::Relaxed);
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
            mark: Instant::now(),
            since_flush: 0,
            flush_every: DEFAULT_FLUSH_EVERY,
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
    mark: Instant,
    since_flush: u64,
    flush_every: u64,
}

impl<'a> BatchTimer<'a> {
    /// Publish every `n` items instead of the default.
    pub fn flush_every(mut self, n: u64) -> Self {
        self.flush_every = n.max(1);
        self
    }

    /// Call once per item processed.
    ///
    /// The count is only ever compared against itself over time, so any
    /// consistent unit will do; it plays no part in sizing the split.
    #[inline]
    pub fn tick(&mut self) {
        self.since_flush += 1;
        if self.since_flush >= self.flush_every {
            self.publish();
        }
    }

    #[inline]
    fn publish(&mut self) {
        let now = Instant::now();
        self.meter
            .nanos
            .fetch_add((now - self.mark).as_nanos() as u64, Ordering::Relaxed);
        self.meter
            .items
            .fetch_add(self.since_flush, Ordering::Relaxed);
        self.mark = now;
        self.since_flush = 0;
    }
}

impl Drop for BatchTimer<'_> {
    fn drop(&mut self) {
        self.publish();
    }
}
