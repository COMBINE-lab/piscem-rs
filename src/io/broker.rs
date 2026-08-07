//! Adapters joining piscem's two resizable pools to the thread broker.
//!
//! The broker knows nothing about mapping or gzip. These translate:
//!
//! * `paraseq::parallel::ThreadPool` + [`MappingStats`] → [`Consumer`]
//! * `rapidgzip_core::DecoderPool` + its handles → [`thread_broker::Producer`]

use thread_broker::{Consumer, Work};
#[cfg(feature = "rapidgzip")]
use thread_broker::{Producer, ProducerPressure};

use crate::io::threads::MappingStats;

/// Validation/application hook for selecting the broker's settled behavior.
///
/// The public library API is
/// [`thread_broker::ThreadBrokerBuilder::steady_state_policy`]. Piscem exposes
/// the same choice through an environment variable so policy A/B runs can use
/// one identical release binary.
#[cfg(feature = "rapidgzip")]
pub const THREAD_BROKER_POLICY_ENV: &str = "PISCEM_THREAD_BROKER_POLICY";

/// Same-binary validation hook for responsive steady-state cadence.
#[cfg(feature = "rapidgzip")]
pub const THREAD_BROKER_PROBE_INTERVAL_ENV: &str = "PISCEM_THREAD_BROKER_PROBE_INTERVAL_MS";

#[cfg(feature = "rapidgzip")]
pub fn broker_policy_from_environment()
-> Result<thread_broker::SteadyStatePolicy, thread_broker::ResizeError> {
    parse_broker_policy(std::env::var(THREAD_BROKER_POLICY_ENV).ok().as_deref())
}

#[cfg(feature = "rapidgzip")]
fn parse_broker_policy(
    value: Option<&str>,
) -> Result<thread_broker::SteadyStatePolicy, thread_broker::ResizeError> {
    use thread_broker::SteadyStatePolicy::{
        FreezeAfterConvergence, FreezeAfterFullCalibration, Responsive,
    };

    match value {
        None | Some("responsive") => Ok(Responsive),
        Some("freeze-after-convergence") => Ok(FreezeAfterConvergence),
        Some("freeze-after-full-calibration") => Ok(FreezeAfterFullCalibration),
        Some(other) => Err(thread_broker::ResizeError::new(format!(
            "invalid {THREAD_BROKER_POLICY_ENV}={other:?}; expected responsive, \
             freeze-after-convergence, or freeze-after-full-calibration"
        ))),
    }
}

#[cfg(feature = "rapidgzip")]
pub fn broker_config_from_environment(
    mut config: thread_broker::BrokerConfig,
) -> Result<thread_broker::BrokerConfig, thread_broker::ResizeError> {
    let Some(value) = std::env::var(THREAD_BROKER_PROBE_INTERVAL_ENV).ok() else {
        return Ok(config);
    };
    let millis = value.parse::<u64>().map_err(|_| {
        thread_broker::ResizeError::new(format!(
            "invalid {THREAD_BROKER_PROBE_INTERVAL_ENV}={value:?}; expected a positive integer number of milliseconds"
        ))
    })?;
    if millis == 0 {
        return Err(thread_broker::ResizeError::new(format!(
            "invalid {THREAD_BROKER_PROBE_INTERVAL_ENV}=0; interval must be non-zero"
        )));
    }
    config.steady_probe_interval = Some(std::time::Duration::from_millis(millis));
    Ok(config)
}

/// The mapping side.
///
/// Owns its meter rather than borrowing the stats it lives in: the broker runs
/// on a thread of its own for the length of the job, so anything it holds has
/// to outlive the stack frame that set it up.
pub struct MappingConsumer {
    pool: paraseq::parallel::ThreadPool,
    busy: std::sync::Arc<thread_broker::BusyMeter>,
    progress_control: Option<MappingProgressControl>,
}

struct MappingProgressControl {
    flush_every: std::sync::Arc<std::sync::atomic::AtomicU64>,
    progress: std::sync::Arc<crate::io::threads::BrokerProgress>,
    calibration: u64,
    monitoring: u64,
}

impl MappingConsumer {
    pub fn new(pool: paraseq::parallel::ThreadPool, stats: &MappingStats) -> Self {
        Self {
            pool,
            busy: std::sync::Arc::clone(&stats.busy),
            progress_control: None,
        }
    }

    /// Use phase-aware in-batch progress publication.
    pub fn with_progress_cadence(
        mut self,
        stats: &MappingStats,
        calibration: u64,
        monitoring: u64,
    ) -> Self {
        self.progress_control = Some(MappingProgressControl {
            flush_every: stats.broker_progress_control(),
            progress: stats.broker_progress(),
            calibration: calibration.max(1),
            monitoring: monitoring.max(1),
        });
        self
    }
}

impl Drop for MappingConsumer {
    fn drop(&mut self) {
        if let Some(control) = &self.progress_control {
            // Freeze mode drops its controller-side consumer at convergence.
            // Existing callbacks retain their cadence until they finish; new
            // callbacks then publish only on drop, removing recurring clock and
            // atomic updates while preserving the per-record counter branch.
            control
                .flush_every
                .store(u64::MAX, std::sync::atomic::Ordering::Release);
        }
    }
}

impl Consumer for MappingConsumer {
    fn set_threads(&self, n: usize) -> Result<(), thread_broker::ResizeError> {
        self.pool.set_threads(n);
        Ok(())
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
        self.pool.total_live()
    }

    fn work(&self) -> Work {
        let mut work = self.busy.work();
        if let Some(control) = &self.progress_control {
            work.items = control.progress.items();
        }
        work
    }

    fn set_measurement_mode(&self, mode: thread_broker::ProducerMeasurementMode) {
        if let Some(control) = &self.progress_control {
            let flush_every = match mode {
                thread_broker::ProducerMeasurementMode::Calibration => control.calibration,
                thread_broker::ProducerMeasurementMode::Monitoring => control.monitoring,
            };
            control
                .flush_every
                .store(flush_every, std::sync::atomic::Ordering::Release);
        }
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
/// [`DecodeProducer::set_limit`] returns any refusal to the broker for that
/// reason: it means the pool was built too small, which is a wiring bug rather
/// than a runtime condition, and continuing would invalidate the budget.
#[cfg(feature = "rapidgzip")]
pub struct DecodeProducer {
    pool: rapidgzip_core::DecoderPool,
    handles: Vec<rapidgzip_core::DecoderHandle>,
    busy: BusySampler,
}

/// Nominal high-resolution polling period, active only while the broker is
/// gathering or ratifying cost evidence. Individual sleeps are deterministically
/// jittered around this value so periodic decoder work cannot phase-lock to it.
#[cfg(feature = "rapidgzip")]
const CALIBRATION_SAMPLE_INTERVAL: std::time::Duration = std::time::Duration::from_millis(3);

/// Settled monitoring needs only enough fidelity to notice a possible regime
/// change. The broker switches back to calibration before making a new move.
#[cfg(feature = "rapidgzip")]
const MONITORING_SAMPLE_INTERVAL: std::time::Duration = std::time::Duration::from_millis(25);

/// Turns instantaneous executing-worker counts into cumulative busy time on a
/// phase-aware sampler.
///
/// This deliberately sums the lock-free per-decoder `busy_workers` snapshots,
/// not pool-permit occupancy. A rapidgzip worker can retain a permit across a
/// nonblocking handoff while its executing flag is clear, so permit time can
/// conservatively exceed actual decode time. Per-decoder busy state changes at
/// the CPU-region boundary the broker intends to measure.
///
/// This remains a sampled estimate. A jittered cadence averaging three
/// milliseconds is used only while a decision is open; steady state polls every
/// 25 ms by default and may poll more sparsely under an explicit low-frequency
/// responsive policy. Both sample counts and time spent taking observations are
/// emitted in machine telemetry.
/// The sampler is one explicitly measured auxiliary thread outside the
/// controlled execution-slot budget.
#[cfg(feature = "rapidgzip")]
struct BusySampler {
    nanos: std::sync::Arc<std::sync::atomic::AtomicU64>,
    stop: std::sync::Arc<std::sync::atomic::AtomicBool>,
    mode: std::sync::Arc<std::sync::atomic::AtomicU8>,
    calibration_samples: std::sync::Arc<std::sync::atomic::AtomicU64>,
    monitoring_samples: std::sync::Arc<std::sync::atomic::AtomicU64>,
    mode_changes: std::sync::Arc<std::sync::atomic::AtomicU64>,
    observation_nanos: std::sync::Arc<std::sync::atomic::AtomicU64>,
    cpu_nanos: std::sync::Arc<std::sync::atomic::AtomicU64>,
    cpu_available: std::sync::Arc<std::sync::atomic::AtomicBool>,
    cpu_accounting_failures: std::sync::Arc<std::sync::atomic::AtomicUsize>,
    calibration_interval: std::time::Duration,
    monitoring_interval_nanos: std::sync::Arc<std::sync::atomic::AtomicU64>,
    join: std::sync::Mutex<Option<std::thread::JoinHandle<()>>>,
}

#[cfg(feature = "rapidgzip")]
impl BusySampler {
    const CALIBRATION: u8 = 0;
    const MONITORING: u8 = 1;

    fn start(
        handles: Vec<rapidgzip_core::DecoderHandle>,
    ) -> Result<Self, thread_broker::ResizeError> {
        Self::start_with(
            CALIBRATION_SAMPLE_INTERVAL,
            MONITORING_SAMPLE_INTERVAL,
            move || {
                handles
                    .iter()
                    .map(|handle| handle.stats().busy_workers)
                    .sum()
            },
        )
    }

    fn start_with(
        calibration_interval: std::time::Duration,
        monitoring_interval: std::time::Duration,
        busy_workers: impl Fn() -> usize + Send + 'static,
    ) -> Result<Self, thread_broker::ResizeError> {
        use std::sync::atomic::Ordering;

        let nanos = std::sync::Arc::new(std::sync::atomic::AtomicU64::new(0));
        let stop = std::sync::Arc::new(std::sync::atomic::AtomicBool::new(false));
        let mode = std::sync::Arc::new(std::sync::atomic::AtomicU8::new(Self::CALIBRATION));
        let calibration_samples = std::sync::Arc::new(std::sync::atomic::AtomicU64::new(0));
        let monitoring_samples = std::sync::Arc::new(std::sync::atomic::AtomicU64::new(0));
        let mode_changes = std::sync::Arc::new(std::sync::atomic::AtomicU64::new(0));
        let observation_nanos = std::sync::Arc::new(std::sync::atomic::AtomicU64::new(0));
        let cpu_nanos = std::sync::Arc::new(std::sync::atomic::AtomicU64::new(0));
        let cpu_available = std::sync::Arc::new(std::sync::atomic::AtomicBool::new(false));
        let cpu_accounting_failures = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));
        let monitoring_interval_nanos = std::sync::Arc::new(std::sync::atomic::AtomicU64::new(
            monitoring_interval.as_nanos() as u64,
        ));
        let thread_nanos = std::sync::Arc::clone(&nanos);
        let thread_stop = std::sync::Arc::clone(&stop);
        let thread_mode = std::sync::Arc::clone(&mode);
        let thread_calibration_samples = std::sync::Arc::clone(&calibration_samples);
        let thread_monitoring_samples = std::sync::Arc::clone(&monitoring_samples);
        let thread_observation_nanos = std::sync::Arc::clone(&observation_nanos);
        let thread_cpu_nanos = std::sync::Arc::clone(&cpu_nanos);
        let thread_cpu_available = std::sync::Arc::clone(&cpu_available);
        let thread_cpu_accounting_failures = std::sync::Arc::clone(&cpu_accounting_failures);
        let thread_monitoring_interval_nanos = std::sync::Arc::clone(&monitoring_interval_nanos);
        let join = std::thread::Builder::new()
            .name("decode-busy-sampler".into())
            .spawn(move || {
                let cpu_timer = thread_broker::ThreadCpuTimer::start();
                let mut sampled_at = std::time::Instant::now();
                let observation_started = std::time::Instant::now();
                let mut busy = busy_workers();
                let mut calibration_jitter_index = 0usize;
                thread_observation_nanos.fetch_add(
                    observation_started.elapsed().as_nanos() as u64,
                    Ordering::Relaxed,
                );
                loop {
                    let interval = if thread_mode.load(Ordering::Acquire) == Self::CALIBRATION {
                        // Mean multiplier is exactly one. The deliberately
                        // irregular seven-sample cycle breaks cadence aliasing
                        // without increasing the average polling rate.
                        const FIFTHS: [u32; 7] = [3, 7, 4, 6, 2, 8, 5];
                        let fifths = FIFTHS[calibration_jitter_index];
                        calibration_jitter_index = (calibration_jitter_index + 1) % FIFTHS.len();
                        calibration_interval.saturating_mul(fifths) / 5
                    } else {
                        std::time::Duration::from_nanos(
                            thread_monitoring_interval_nanos.load(Ordering::Acquire),
                        )
                    };
                    std::thread::park_timeout(interval);
                    let now = std::time::Instant::now();
                    let observation_started = std::time::Instant::now();
                    let observed_busy = busy_workers();
                    thread_observation_nanos.fetch_add(
                        observation_started.elapsed().as_nanos() as u64,
                        Ordering::Relaxed,
                    );
                    let elapsed = now.saturating_duration_since(sampled_at).as_nanos() as u64;
                    let endpoint_sum = busy.saturating_add(observed_busy) as u64;
                    thread_nanos
                        .fetch_add(elapsed.saturating_mul(endpoint_sum) / 2, Ordering::Relaxed);
                    sampled_at = now;
                    busy = observed_busy;
                    if thread_mode.load(Ordering::Acquire) == Self::CALIBRATION {
                        thread_calibration_samples.fetch_add(1, Ordering::Relaxed);
                    } else {
                        thread_monitoring_samples.fetch_add(1, Ordering::Relaxed);
                    }
                    if thread_stop.load(Ordering::Acquire) {
                        break;
                    }
                }
                if let Some(elapsed) = cpu_timer.elapsed()
                    && let Ok(nanos) = u64::try_from(elapsed.as_nanos())
                {
                    thread_cpu_nanos.store(nanos, Ordering::Relaxed);
                    thread_cpu_available.store(true, Ordering::Release);
                } else {
                    thread_cpu_accounting_failures.fetch_add(1, Ordering::Relaxed);
                }
            })
            .map_err(|error| {
                thread_broker::ResizeError::new(format!(
                    "could not spawn decode busy-time sampler: {error}"
                ))
            })?;
        Ok(Self {
            nanos,
            stop,
            mode,
            calibration_samples,
            monitoring_samples,
            mode_changes,
            observation_nanos,
            cpu_nanos,
            cpu_available,
            cpu_accounting_failures,
            calibration_interval,
            monitoring_interval_nanos,
            join: std::sync::Mutex::new(Some(join)),
        })
    }

    fn nanos(&self) -> u64 {
        use std::sync::atomic::Ordering;
        self.nanos.load(Ordering::Relaxed)
    }

    fn set_mode(&self, mode: thread_broker::ProducerMeasurementMode) {
        use std::sync::atomic::Ordering;

        let encoded = match mode {
            thread_broker::ProducerMeasurementMode::Calibration => Self::CALIBRATION,
            thread_broker::ProducerMeasurementMode::Monitoring => Self::MONITORING,
        };
        if self.mode.swap(encoded, Ordering::AcqRel) != encoded {
            self.mode_changes.fetch_add(1, Ordering::Relaxed);
            if let Ok(join) = self.join.lock()
                && let Some(join) = join.as_ref()
            {
                join.thread().unpark();
            }
        }
    }

    fn set_monitoring_interval(&self, interval: std::time::Duration) {
        use std::sync::atomic::Ordering;

        let nanos = u64::try_from(interval.as_nanos())
            .unwrap_or(u64::MAX)
            .max(1);
        self.monitoring_interval_nanos
            .store(nanos, Ordering::Release);
        if self.mode.load(Ordering::Acquire) == Self::MONITORING
            && let Ok(join) = self.join.lock()
            && let Some(join) = join.as_ref()
        {
            join.thread().unpark();
        }
    }

    fn stats(&self) -> thread_broker::ProducerMeasurementStats {
        use std::sync::atomic::Ordering;

        thread_broker::ProducerMeasurementStats {
            busy_nanos: self.nanos.load(Ordering::Relaxed),
            completed_worker_cpu_nanos: None,
            completed_auxiliary_cpu_nanos: None,
            cpu_accounting_failures: None,
            sampler_cpu_nanos: self
                .cpu_available
                .load(Ordering::Acquire)
                .then(|| self.cpu_nanos.load(Ordering::Relaxed)),
            sampler_cpu_accounting_failures: self.cpu_accounting_failures.load(Ordering::Relaxed),
            calibration_samples: self.calibration_samples.load(Ordering::Relaxed),
            monitoring_samples: self.monitoring_samples.load(Ordering::Relaxed),
            mode_changes: self.mode_changes.load(Ordering::Relaxed),
            final_mode: if self.mode.load(Ordering::Acquire) == Self::CALIBRATION {
                thread_broker::ProducerMeasurementMode::Calibration
            } else {
                thread_broker::ProducerMeasurementMode::Monitoring
            },
            observation_nanos: self.observation_nanos.load(Ordering::Relaxed),
            calibration_interval_micros: self.calibration_interval.as_micros() as u64,
            monitoring_interval_micros: self.monitoring_interval_nanos.load(Ordering::Acquire)
                / 1_000,
        }
    }

    fn stop_and_join(&self) {
        use std::sync::atomic::Ordering;

        self.stop.store(true, Ordering::Release);
        if let Ok(mut slot) = self.join.lock()
            && let Some(join) = slot.take()
        {
            join.thread().unpark();
            if join.join().is_err() {
                tracing::error!("decode busy-time sampler panicked");
                self.cpu_accounting_failures.fetch_add(1, Ordering::Relaxed);
            }
        }
    }
}

/// Same-binary control for measuring observation overhead at a fixed split.
///
/// This is intentionally an environment-only validation hook rather than a
/// user-facing scheduling control. With a pinned aggregate decode allocation,
/// unset/`off` runs no sampler, `calibration` runs the high-resolution cadence,
/// and `monitoring` runs the settled cadence. The mapping/decode allocation is
/// identical in all three cases.
#[cfg(feature = "rapidgzip")]
pub const FIXED_DECODE_MEASUREMENT_ENV: &str = "PISCEM_FIXED_DECODE_MEASUREMENT";

#[cfg(feature = "rapidgzip")]
pub struct FixedDecodeMeasurement {
    sampler: BusySampler,
    handles: Vec<rapidgzip_core::DecoderHandle>,
}

#[cfg(feature = "rapidgzip")]
impl FixedDecodeMeasurement {
    pub fn from_environment(
        handles: Vec<rapidgzip_core::DecoderHandle>,
        adaptive: bool,
    ) -> Result<Option<Self>, thread_broker::ResizeError> {
        let value = std::env::var(FIXED_DECODE_MEASUREMENT_ENV).ok();
        let Some(mode) = parse_fixed_measurement_mode(value.as_deref())? else {
            return Ok(None);
        };
        if adaptive {
            return Err(thread_broker::ResizeError::new(format!(
                "{FIXED_DECODE_MEASUREMENT_ENV} is a fixed-split overhead control and cannot be combined with adaptive decoding"
            )));
        }
        if handles.is_empty() {
            return Err(thread_broker::ResizeError::new(format!(
                "{FIXED_DECODE_MEASUREMENT_ENV} requires at least one parallel decoder handle"
            )));
        }
        let sampler = BusySampler::start(handles.clone())?;
        sampler.set_mode(mode);
        Ok(Some(Self { sampler, handles }))
    }

    pub fn finish(self) -> thread_broker::ProducerMeasurementStats {
        self.sampler.stop_and_join();
        measurement_stats_with_cpu(self.sampler.stats(), &self.handles)
    }
}

#[cfg(feature = "rapidgzip")]
fn measurement_stats_with_cpu(
    mut measurement: thread_broker::ProducerMeasurementStats,
    handles: &[rapidgzip_core::DecoderHandle],
) -> thread_broker::ProducerMeasurementStats {
    refresh_measurement_cpu(&mut measurement, handles);
    measurement
}

/// Refresh lifetime CPU counters after mapping has ended.
///
/// A freeze-after-convergence broker intentionally releases its sampler early,
/// so its report cannot yet contain CPU time from producer workers that finish
/// later. The rapidgzip handle counters remain available to the application and
/// are refreshed at pipeline shutdown.
#[cfg(feature = "rapidgzip")]
pub fn refresh_measurement_cpu(
    measurement: &mut thread_broker::ProducerMeasurementStats,
    handles: &[rapidgzip_core::DecoderHandle],
) {
    fn duration_sum(
        mut durations: impl Iterator<Item = Option<std::time::Duration>>,
    ) -> Option<u64> {
        durations.try_fold(0u64, |total, duration| {
            let nanos = u64::try_from(duration?.as_nanos()).ok()?;
            total.checked_add(nanos)
        })
    }

    measurement.completed_worker_cpu_nanos = duration_sum(
        handles
            .iter()
            .map(|handle| handle.stats().completed_worker_cpu_time),
    );
    measurement.completed_auxiliary_cpu_nanos = duration_sum(
        handles
            .iter()
            .map(|handle| handle.stats().completed_auxiliary_cpu_time),
    );
    measurement.cpu_accounting_failures = handles
        .iter()
        .map(|handle| handle.stats().cpu_accounting_failures)
        .try_fold(0usize, |total, failures| total.checked_add(failures?));
}

#[cfg(feature = "rapidgzip")]
fn parse_fixed_measurement_mode(
    value: Option<&str>,
) -> Result<Option<thread_broker::ProducerMeasurementMode>, thread_broker::ResizeError> {
    match value {
        None | Some("") | Some("off") => Ok(None),
        Some("calibration") => Ok(Some(thread_broker::ProducerMeasurementMode::Calibration)),
        Some("monitoring") => Ok(Some(thread_broker::ProducerMeasurementMode::Monitoring)),
        Some(value) => Err(thread_broker::ResizeError::new(format!(
            "{FIXED_DECODE_MEASUREMENT_ENV} must be off, calibration, or monitoring, not {value:?}"
        ))),
    }
}

#[cfg(all(test, feature = "rapidgzip"))]
mod busy_sampler_tests {
    use super::{
        BusySampler, CALIBRATION_SAMPLE_INTERVAL, MONITORING_SAMPLE_INTERVAL, MappingConsumer,
        parse_broker_policy, parse_fixed_measurement_mode,
    };
    use crate::io::threads::MappingStats;
    use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
    use std::sync::{Arc, Mutex};
    use std::time::{Duration, Instant};

    struct Truth {
        nanos: u64,
        at: Instant,
        busy: usize,
    }

    impl Truth {
        fn set(&mut self, now: Instant, busy: usize) {
            self.nanos = self.nanos.saturating_add(
                now.saturating_duration_since(self.at)
                    .as_nanos()
                    .saturating_mul(self.busy as u128) as u64,
            );
            self.at = now;
            self.busy = busy;
        }

        fn snapshot(&self, now: Instant) -> u64 {
            self.nanos.saturating_add(
                now.saturating_duration_since(self.at)
                    .as_nanos()
                    .saturating_mul(self.busy as u128) as u64,
            )
        }
    }

    #[test]
    fn calibration_sampling_tracks_bursty_ground_truth() {
        let busy = Arc::new(AtomicUsize::new(0));
        let truth = Arc::new(Mutex::new(Truth {
            nanos: 0,
            at: Instant::now(),
            busy: 0,
        }));
        let stop = Arc::new(AtomicBool::new(false));
        let sampler =
            BusySampler::start_with(CALIBRATION_SAMPLE_INTERVAL, MONITORING_SAMPLE_INTERVAL, {
                let busy = Arc::clone(&busy);
                move || busy.load(Ordering::Acquire)
            })
            .unwrap();
        let pulse = {
            let busy = Arc::clone(&busy);
            let truth = Arc::clone(&truth);
            let stop = Arc::clone(&stop);
            std::thread::spawn(move || {
                let mut active = false;
                while !stop.load(Ordering::Acquire) {
                    active = !active;
                    let value = if active { 4 } else { 0 };
                    let now = Instant::now();
                    truth.lock().unwrap().set(now, value);
                    busy.store(value, Ordering::Release);
                    std::thread::sleep(Duration::from_millis(5));
                }
            })
        };

        let started = Instant::now();
        let intervals = [25u64, 50, 73, 100, 137, 250];
        let mut windows = Vec::new();
        let mut previous_sampled = 0u64;
        let mut previous_truth = 0u64;
        let mut interval_index = 0usize;
        while started.elapsed() < Duration::from_secs(6) {
            std::thread::sleep(Duration::from_millis(intervals[interval_index]));
            let now = Instant::now();
            let sampled = sampler.nanos();
            let exact = truth.lock().unwrap().snapshot(now);
            windows.push((
                interval_index,
                sampled.saturating_sub(previous_sampled),
                exact.saturating_sub(previous_truth),
            ));
            previous_sampled = sampled;
            previous_truth = exact;
            interval_index = (interval_index + 1) % intervals.len();
        }
        stop.store(true, Ordering::Release);
        pulse.join().unwrap();

        let sampled_total: u64 = windows.iter().map(|window| window.1).sum();
        let exact_total: u64 = windows.iter().map(|window| window.2).sum();
        let whole_error = sampled_total.abs_diff(exact_total) as f64 / exact_total as f64;
        assert!(
            whole_error <= 0.03,
            "whole-run busy-time error was {:.1}%",
            whole_error * 100.0
        );

        let mut aggregate_errors: Vec<_> = windows
            .windows(3)
            .filter_map(|three| {
                let sampled: u64 = three.iter().map(|window| window.1).sum();
                let exact: u64 = three.iter().map(|window| window.2).sum();
                (exact > 0).then_some(sampled.abs_diff(exact) as f64 / exact as f64)
            })
            .collect();
        aggregate_errors.sort_by(f64::total_cmp);
        let p95 =
            aggregate_errors[(aggregate_errors.len() * 95 / 100).min(aggregate_errors.len() - 1)];
        assert!(
            p95 <= 0.10,
            "three-window p95 error was {:.1}%",
            p95 * 100.0
        );

        // Hold an exact consumer cost at three times producer cost. If the
        // producer counter changes the solved share by cadence, it violates the
        // controller's model even when the raw relative-error thresholds pass.
        for (interval_index, interval) in intervals.iter().enumerate() {
            let (sampled, exact) = windows
                .iter()
                .filter(|window| window.0 == interval_index)
                .fold((0u64, 0u64), |(sampled, exact), window| {
                    (sampled + window.1, exact + window.2)
                });
            let consumer = exact.saturating_mul(3);
            let solved_share = sampled as f64 / sampled.saturating_add(consumer) as f64;
            assert!(
                (solved_share - 0.25).abs() <= 0.02,
                "{} ms cadence shifted cost share to {:.1}%",
                interval,
                solved_share * 100.0,
            );
        }
    }

    #[test]
    fn settled_monitoring_reduces_the_poll_rate() {
        let sampler =
            BusySampler::start_with(Duration::from_millis(1), Duration::from_millis(50), || 0)
                .unwrap();
        std::thread::sleep(Duration::from_millis(200));
        let calibration_samples = sampler.stats().calibration_samples;
        sampler.set_mode(thread_broker::ProducerMeasurementMode::Monitoring);
        std::thread::sleep(Duration::from_millis(200));
        let stats = sampler.stats();

        assert_eq!(stats.mode_changes, 1);
        assert!(
            calibration_samples > stats.monitoring_samples.saturating_mul(5),
            "calibration={} monitoring={}",
            calibration_samples,
            stats.monitoring_samples,
        );
        assert!(stats.observation_nanos > 0);
    }

    #[test]
    fn mapping_progress_is_fine_only_while_a_decision_is_open() {
        use thread_broker::Consumer as _;
        use thread_broker::ProducerMeasurementMode::{Calibration, Monitoring};

        let stats = MappingStats::new();
        stats.set_broker_progress_flush_every(32);
        {
            let consumer =
                MappingConsumer::new(paraseq::parallel::ThreadPool::with_max(1, 1), &stats)
                    .with_progress_cadence(&stats, 32, 256);
            consumer.set_measurement_mode(Calibration);
            assert_eq!(stats.broker_progress_flush_every(), 32);
            consumer.set_measurement_mode(Monitoring);
            assert_eq!(stats.broker_progress_flush_every(), 256);
        }
        assert_eq!(stats.broker_progress_flush_every(), u64::MAX);
    }

    #[test]
    fn settled_monitoring_interval_is_runtime_configurable() {
        let sampler =
            BusySampler::start_with(Duration::from_millis(1), Duration::from_millis(25), || 0)
                .unwrap();
        sampler.set_mode(thread_broker::ProducerMeasurementMode::Monitoring);
        sampler.set_monitoring_interval(Duration::from_millis(200));
        std::thread::sleep(Duration::from_millis(450));
        sampler.stop_and_join();
        let stats = sampler.stats();

        assert_eq!(stats.monitoring_interval_micros, 200_000);
        assert!(
            (1..=6).contains(&stats.monitoring_samples),
            "unexpected monitoring sample count: {stats:?}"
        );
        #[cfg(target_os = "linux")]
        {
            assert!(stats.sampler_cpu_nanos.is_some());
            assert_eq!(stats.sampler_cpu_accounting_failures, 0);
        }
    }

    #[test]
    fn fixed_measurement_control_is_explicit() {
        use thread_broker::ProducerMeasurementMode::{Calibration, Monitoring};

        assert_eq!(parse_fixed_measurement_mode(None).unwrap(), None);
        assert_eq!(parse_fixed_measurement_mode(Some("off")).unwrap(), None);
        assert_eq!(
            parse_fixed_measurement_mode(Some("calibration")).unwrap(),
            Some(Calibration)
        );
        assert_eq!(
            parse_fixed_measurement_mode(Some("monitoring")).unwrap(),
            Some(Monitoring)
        );
        assert!(parse_fixed_measurement_mode(Some("yes")).is_err());
    }

    #[test]
    fn broker_policy_control_is_explicit() {
        use thread_broker::SteadyStatePolicy::{
            FreezeAfterConvergence, FreezeAfterFullCalibration, Responsive,
        };

        assert_eq!(parse_broker_policy(None).unwrap(), Responsive);
        assert_eq!(parse_broker_policy(Some("responsive")).unwrap(), Responsive);
        assert_eq!(
            parse_broker_policy(Some("freeze-after-convergence")).unwrap(),
            FreezeAfterConvergence
        );
        assert_eq!(
            parse_broker_policy(Some("freeze-after-full-calibration")).unwrap(),
            FreezeAfterFullCalibration
        );
        assert!(parse_broker_policy(Some("freeze")).is_err());
    }
}

#[cfg(feature = "rapidgzip")]
impl Drop for BusySampler {
    fn drop(&mut self) {
        self.stop_and_join();
    }
}

#[cfg(feature = "rapidgzip")]
impl DecodeProducer {
    pub fn new(
        pool: rapidgzip_core::DecoderPool,
        handles: Vec<rapidgzip_core::DecoderHandle>,
    ) -> Result<Self, thread_broker::ResizeError> {
        let busy = BusySampler::start(handles.clone())?;
        Ok(Self {
            busy,
            pool,
            handles,
        })
    }
}

#[cfg(feature = "rapidgzip")]
impl Producer for DecodeProducer {
    fn set_limit(&self, n: usize) -> Result<(), thread_broker::ResizeError> {
        self.pool.set_worker_limit(n.max(1)).map_err(|error| {
            thread_broker::ResizeError::new(format!(
                "decode pool refused worker limit {n}: {error}; its immutable worker ceiling must equal the whole execution budget"
            ))
        })
    }

    fn limit(&self) -> usize {
        self.pool.stats().worker_limit
    }

    fn active_slots(&self) -> usize {
        self.pool.stats().busy_workers
    }

    fn buffered_items(&self) -> Option<u64> {
        Some(
            self.handles
                .iter()
                .map(|handle| {
                    let stats = handle.stats();
                    stats
                        .decompressed_bytes
                        .saturating_sub(stats.consumed_bytes)
                })
                .sum(),
        )
    }

    fn auxiliary_threads(&self) -> usize {
        // rapidgzip coordinators/scanners plus this adapter's busy-time sampler.
        self.pool.stats().auxiliary_threads.saturating_add(1)
    }

    fn set_measurement_mode(&self, mode: thread_broker::ProducerMeasurementMode) {
        self.busy.set_mode(mode);
    }

    fn set_monitoring_interval(&self, interval: std::time::Duration) {
        // A shorter cadence did not improve the deterministic burst gate, while
        // 25 ms already passed the measured <=1% steady-state overhead gate.
        // Longer responsive probe intervals can safely lower the observation
        // rate while retaining roughly four samples per controller window.
        self.busy
            .set_monitoring_interval(interval.max(MONITORING_SAMPLE_INTERVAL));
    }

    fn measurement_stats(&self) -> Option<thread_broker::ProducerMeasurementStats> {
        Some(measurement_stats_with_cpu(self.busy.stats(), &self.handles))
    }

    fn finish_measurement(&self) -> Option<thread_broker::ProducerMeasurementStats> {
        self.busy.stop_and_join();
        Some(measurement_stats_with_cpu(self.busy.stats(), &self.handles))
    }

    /// Decode time with blocking excluded, plus bytes as a progress measure.
    ///
    /// Busy time is sampled independently of the broker's decision cadence; a
    /// slow or irregular caller therefore cannot alias the estimate.
    fn work(&self) -> Work {
        let busy_nanos = self.busy.nanos();
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
