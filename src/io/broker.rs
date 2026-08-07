//! Adapters joining piscem's two resizable pools to the thread broker.
//!
//! The broker knows nothing about mapping or gzip. These translate:
//!
//! * `paraseq::parallel::ThreadPool` + [`MappingStats`] → [`Consumer`]
//! * `rapidgzip_core::DecoderPool` + its handles → [`Producer`]

use thread_broker::{Consumer, Work};
#[cfg(feature = "rapidgzip")]
use thread_broker::{Producer, ProducerPressure};

use crate::io::threads::MappingStats;

/// The mapping side.
///
/// Owns its meter rather than borrowing the stats it lives in: the broker runs
/// on a thread of its own for the length of the job, so anything it holds has
/// to outlive the stack frame that set it up.
pub struct MappingConsumer {
    pool: paraseq::parallel::ThreadPool,
    busy: std::sync::Arc<thread_broker::BusyMeter>,
}

impl MappingConsumer {
    pub fn new(pool: paraseq::parallel::ThreadPool, stats: &MappingStats) -> Self {
        Self {
            pool,
            busy: std::sync::Arc::clone(&stats.busy),
        }
    }
}

impl Consumer for MappingConsumer {
    fn set_threads(&self, n: usize) {
        self.pool.set_threads(n);
    }

    /// Workers really running, which lags `set_threads` by up to one batch.
    ///
    /// `ThreadPool::threads()` would return the target instead, and using that
    /// to normalise busy time while the pool is still converging divides real
    /// work by capacity that did not exist yet.
    fn live_threads(&self) -> usize {
        // The *aggregate*, not this share. `Collection` splits the pool across
        // readers, and the pool held here is the parent — whose own share never
        // runs anything, so `live()` reads zero for the whole run and every
        // consumer measures as 0% idle however starved it is.
        self.pool.total_live().max(1)
    }

    fn work(&self) -> Work {
        self.busy.work()
    }
}

/// The decode side.
///
/// Holds the handles as well as the pool because the pool's aggregate view
/// cannot answer the two questions that matter most: whether a decoder has
/// fallen back to sequential decoding, and how much more concurrency the
/// decoders would actually use.
///
/// # Size the pool's `workers` to the whole thread budget
///
/// `DecoderPool::workers` is an immutable maximum and `set_worker_limit` is
/// refused above it. Construct the pool with the *entire* budget and let the
/// broker move the mutable limit within it — do not size `workers` to the split
/// you expect.
///
/// Otherwise the broker can reach a state it cannot act on: the decoders are
/// genuinely starved, it grants more, the pool refuses, and the broker's own
/// view of the split silently diverges from the pool's. The broker tracks what
/// it asked for rather than reading it back, precisely so two decisions cannot
/// double-spend the same threads, and a refused grant breaks that.
/// [`DecodeProducer::set_limit`] logs any refusal at `WARN` for that reason: it
/// means the pool was built too small, which is a wiring bug rather than a
/// runtime condition.
#[cfg(feature = "rapidgzip")]
pub struct DecodeProducer {
    pool: rapidgzip_core::DecoderPool,
    handles: Vec<rapidgzip_core::DecoderHandle>,
    busy: std::sync::Mutex<BusyIntegrator>,
}

/// Turns the pool's instantaneous busy-slot count into cumulative busy time.
///
/// The broker needs decode time with blocking excluded, and the pool does not
/// publish a cumulative counter — but it does publish `busy_workers`, and the
/// pool's contract is precisely the one required: *"a task releases its slot
/// before any result or reader-channel operation that would block"*. So a busy
/// slot is a slot doing inflate work, never one waiting on the source or on a
/// slow consumer, and integrating that count over time gives exactly the
/// quantity the model wants.
///
/// # The aliasing this is exposed to
///
/// It is a left Riemann sum over an instantaneous sample, so it inherits that
/// sample's aliasing. The broker's own sampling interval is deliberately not a
/// multiple of the decoder's 250 ms worker-retirement hysteresis for this
/// reason. Two things then keep the error small: the sum is over the whole run
/// rather than a single window, so independent sampling error averages down; and
/// the split depends on the *ratio* of the two stages' busy times, so a
/// consistent proportional bias in this one cancels rather than accumulating.
#[cfg(feature = "rapidgzip")]
#[derive(Default)]
struct BusyIntegrator {
    nanos: u64,
    last: Option<std::time::Instant>,
}

#[cfg(feature = "rapidgzip")]
impl BusyIntegrator {
    fn accumulate(&mut self, busy_workers: usize) -> u64 {
        let now = std::time::Instant::now();
        if let Some(last) = self.last {
            let dt = now.saturating_duration_since(last).as_nanos() as u64;
            self.nanos += dt.saturating_mul(busy_workers as u64);
        }
        self.last = Some(now);
        self.nanos
    }
}

#[cfg(feature = "rapidgzip")]
impl DecodeProducer {
    pub fn new(
        pool: rapidgzip_core::DecoderPool,
        handles: Vec<rapidgzip_core::DecoderHandle>,
    ) -> Self {
        Self {
            pool,
            handles,
            busy: std::sync::Mutex::new(BusyIntegrator::default()),
        }
    }
}

#[cfg(feature = "rapidgzip")]
impl Producer for DecodeProducer {
    fn set_limit(&self, n: usize) {
        // A refusal means the pool was constructed smaller than the budget the
        // broker is dividing, so the broker's accounting has already diverged
        // from reality. That is a wiring bug and worth saying so out loud.
        if let Err(e) = self.pool.set_worker_limit(n.max(1)) {
            tracing::warn!(
                "decode pool refused a worker limit of {n}: {e}. The pool must be \
                 built with `workers` equal to the whole thread budget."
            );
        }
    }

    fn limit(&self) -> usize {
        self.pool.stats().worker_limit
    }

    /// Decode time with blocking excluded, plus bytes as a progress measure.
    ///
    /// Must be called at a roughly steady cadence -- the broker's sampling loop
    /// is the only caller -- since the busy time is integrated from the samples
    /// this method itself takes. See [`BusyIntegrator`].
    fn work(&self) -> Work {
        let pool = self.pool.stats();
        let busy_nanos = self
            .busy
            .lock()
            .map(|mut b| b.accumulate(pool.busy_workers))
            .unwrap_or(0);
        // Only ever compared against itself over time, and never against the
        // consumer's count, so decompressed bytes are a perfectly good unit even
        // though the consumer counts reads.
        let items = self
            .handles
            .iter()
            .map(|h| h.stats().decompressed_bytes)
            .sum();
        Work { busy_nanos, items }
    }

    /// Classify the decode side.
    ///
    /// # How ambiguity is resolved, and why
    ///
    /// Upstream deliberately does *not* report a source-bound state: its own
    /// docs say it cannot yet separate source waiting from inflate time
    /// reliably enough to promise as a contract. So [`ProducerPressure::SourceBound`]
    /// here is **inferred**, and the direction of that inference matters,
    /// because the two mistakes are not symmetric:
    ///
    /// * calling a genuinely starved decoder source-bound makes the broker
    ///   *reclaim* decode threads — the expensive error, worth 22-48% when
    ///   decode is really the constraint;
    /// * calling a genuinely source-bound decoder starved wastes some threads
    ///   on decode that will idle — measured at a few percent of CPU.
    ///
    /// So source-bound is claimed only on strong evidence: nothing queued,
    /// nothing executing, and our limit demonstrably not the constraint. A
    /// decoder with data to work on would be busy or have work queued; one held
    /// back by us would be waiting for a slot. Anything short of that resolves
    /// toward `Starved`.
    fn pressure(&self) -> ProducerPressure {
        use rapidgzip_core::{DecoderPath, DecoderPressure};

        let pool = self.pool.stats();
        let snaps: Vec<_> = self.handles.iter().map(|h| h.stats()).collect();
        if snaps.is_empty() {
            return ProducerPressure::Satisfied;
        }

        // Checked first: it dominates everything else. A decoder that committed
        // to the sequential backend reports `Idle` with no workers forever and
        // cannot be helped at any limit, so every other signal it emits is
        // misleading. Only when *every* decoder is in that state is the whole
        // producer inelastic -- one sequential file among several parallel ones
        // still leaves concurrency worth buying.
        // Finished decoders are excluded: a file that is done says nothing about
        // whether the ones still running want more concurrency. When they are
        // all finished there is nothing left to size.
        let considered: Vec<_> = snaps
            .iter()
            .filter(|s| !matches!(s.pressure, DecoderPressure::Finished))
            .collect();
        if considered.is_empty() {
            return ProducerPressure::Satisfied;
        }
        if considered.iter().all(|s| s.path == DecoderPath::Sequential) {
            return ProducerPressure::Inelastic;
        }

        let limited = considered.iter().any(|s| s.pool_limited)
            || pool.waiting_decoders > 0
            || considered
                .iter()
                .any(|s| matches!(s.pressure, DecoderPressure::DecoderBound { .. }));

        if limited {
            return ProducerPressure::Starved;
        }

        // Ahead of the consumer, or settled at a measured knee: more slots
        // would idle.
        if considered.iter().any(|s| {
            matches!(
                s.pressure,
                DecoderPressure::ConsumerBound { .. } | DecoderPressure::Converged { .. }
            )
        }) {
            return ProducerPressure::Satisfied;
        }

        // Nothing queued, nothing executing, and we are not the constraint:
        // the decoders are waiting on their input.
        let idle_everywhere = considered
            .iter()
            .all(|s| matches!(s.pressure, DecoderPressure::Idle))
            && pool.queued_tasks == 0
            && pool.busy_workers < pool.worker_limit;
        if idle_everywhere {
            return ProducerPressure::SourceBound;
        }

        let desired: usize = considered.iter().map(|s| s.desired_workers).sum();
        if desired > pool.worker_limit {
            ProducerPressure::Starved
        } else {
            ProducerPressure::Satisfied
        }
    }
}
