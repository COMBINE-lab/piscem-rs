//! Process-wide decoder-thread budget across all input files.
//!
//! # Division of responsibility
//!
//! `rapidgzip-core` adapts *within* one decoder: its calibrator sizes the worker
//! population empirically, `configured_workers` is a lazy ceiling rather than an
//! eager allocation, and sustained backpressure at the final reader handoff
//! already drops task admission and retires workers. None of that needs help.
//!
//! What it deliberately does not provide (see upstream issue #4) is a shared
//! pool *across* decoders. A paired run opens two decoders and a multi-file run
//! opens 2N, each sizing itself in ignorance of the others, so the process can
//! hold far more decoder threads than the user's `-t` budget implies.
//!
//! This module owns exactly that gap: it holds one [`DecoderHandle`] per open
//! file and continuously divides a single decode budget between them according
//! to the pressure each reports. It does not try to second-guess the per-decoder
//! calibrator — it only ever sets a *ceiling*, and the calibrator remains free
//! to sit below it.
//!
//! # Why a ceiling and not a target
//!
//! `set_worker_limit` caps admission; it cannot force workers into existence.
//! Handing a starved decoder a higher ceiling lets its own controller grow if
//! the work is really there, and costs nothing if it is not. That keeps one
//! authority over "how many workers actually help" (upstream) and one over "how
//! many this process can afford" (here).

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::thread::JoinHandle;
use std::time::Duration;

use rapidgzip_core::{DecoderHandle, DecoderPressure};

/// How often the supervisor samples telemetry.
///
/// `stats()` is relaxed atomic loads, so sampling is cheap; the interval is set
/// by how fast a decoder's state can meaningfully change. Upstream's worker
/// retirement hysteresis is 250 ms, so sampling much faster than that mostly
/// re-observes the same state.
const SAMPLE_INTERVAL: Duration = Duration::from_millis(250);

/// Evidence needed before a decoder's share is changed, in net samples.
///
/// A *score*, not a consecutive-run counter: `Starved` adds one, `Slack`
/// subtracts one, and the share moves when the score crosses this. Consecutive
/// runs were tried first and do not work — a healthy pipeline alternates
/// between "workers queued" and "parser cannot drain" from sample to sample, so
/// a run-length rule either never fires or, if made short enough to fire,
/// oscillates. Accumulating evidence is robust to that alternation while still
/// ignoring single samples.
const EVIDENCE_THRESHOLD: i32 = 2;

/// Bound on the score, so a long steady phase cannot build up so much inertia
/// that the controller stops responding when conditions change.
const EVIDENCE_CLAMP: i32 = 4;

/// Workers each decoder starts with, before it has shown it needs more.
///
/// Deliberately small. Measured across three workloads, the cost asymmetry is
/// stark and one-sided: over-provisioning decode costs threads but never time
/// (every workload's curve is flat past its knee), while under-provisioning
/// costs 22% (bulk) to 48% (Flex). But the knee itself is not portable — it sat
/// at 16 workers/file for Flex, 4 for bulk, and 0 for PBMC-on-transcriptome,
/// because a large index makes mapping so much costlier per read that decode
/// demand nearly vanishes. No static starting point serves all three, so start
/// below all of them and let demonstrated starvation pay for growth.
const INITIAL_WORKERS_PER_FILE: usize = 2;

/// Smallest share change worth making, in workers. Avoids thrashing the
/// decoder for adjustments that cannot matter given the measured knee.
const MIN_ADJUSTMENT: usize = 2;

/// Summary of what the decoders actually did, for end-of-run reporting.
#[derive(Debug, Default, Clone)]
pub struct DecodeBudgetReport {
    /// Largest number of live decoder OS threads observed across all files.
    pub peak_worker_threads: usize,
    /// Largest number of live coordinator/scanner threads observed.
    pub peak_auxiliary_threads: usize,
    /// Largest number of workers observed actually executing decode tasks.
    pub peak_busy_workers: usize,
    /// Per-file worker count chosen by upstream calibration, where it settled.
    pub converged_workers: Vec<usize>,
    /// Decoder threads the budget permitted but that were never created.
    pub unused_budget: usize,
}

/// Divides one decoder-thread budget across every open input file.
pub struct DecodeBudget {
    stop: Arc<AtomicBool>,
    join: Option<JoinHandle<DecodeBudgetReport>>,
    /// Mirrors the report's peak so callers can read it without joining.
    peak_threads: Arc<AtomicUsize>,
}

impl DecodeBudget {
    /// Start supervising `handles` under a total ceiling of `budget` workers.
    ///
    /// Returns `None` when there is nothing to supervise (no gzip inputs), so
    /// callers pay nothing on the plain-FASTQ path.
    pub fn spawn(handles: Vec<DecoderHandle>, budget: usize) -> Option<Self> {
        if handles.is_empty() {
            return None;
        }
        let stop = Arc::new(AtomicBool::new(false));
        let peak_threads = Arc::new(AtomicUsize::new(0));

        // The fair-share limit is applied at open time (see `open_gz_rapidgzip`);
        // by the time this runs, paraseq has already read a first batch.

        let worker_stop = Arc::clone(&stop);
        let worker_peak = Arc::clone(&peak_threads);
        let join = std::thread::Builder::new()
            .name("decode-budget".into())
            .spawn(move || supervise(handles, budget, worker_stop, worker_peak))
            .ok()?;

        Some(Self {
            stop,
            join: Some(join),
            peak_threads,
        })
    }

    /// Peak live decoder OS threads observed so far.
    pub fn peak_threads(&self) -> usize {
        self.peak_threads.load(Ordering::Relaxed)
    }

    /// Stop supervising and collect the report.
    pub fn finish(mut self) -> DecodeBudgetReport {
        self.stop.store(true, Ordering::Relaxed);
        self.join
            .take()
            .and_then(|j| j.join().ok())
            .unwrap_or_default()
    }
}

impl Drop for DecodeBudget {
    fn drop(&mut self) {
        self.stop.store(true, Ordering::Relaxed);
        if let Some(j) = self.join.take() {
            let _ = j.join();
        }
    }
}

/// The standing verdict on one decoder, distilled from repeated observations.
///
/// Deliberately coarse. Upstream already runs a fast inner loop: its calibrator
/// sizes workers empirically and its handoff feedback reacts to a slow parser
/// within milliseconds. A second fast loop on top fights it. This one is slow,
/// mostly monotone, and only answers a question upstream cannot: how should a
/// *fixed shared budget* be divided between files.
#[derive(Clone, Copy, PartialEq, Debug)]
enum Verdict {
    /// Not enough evidence yet; leave the share alone.
    Undecided,
    /// Persistently starved for admission — worth a larger share if spare exists.
    Starved,
    /// Persistently unable to use what it has, or done.
    Slack,
    /// Calibration settled at a specific count; never exceed it.
    Settled(usize),
}

fn classify(pressure: &DecoderPressure) -> Verdict {
    match pressure {
        // Runnable work queued behind the admission limit is the one signal
        // that more workers could actually buy throughput.
        DecoderPressure::DecoderBound { .. } => Verdict::Starved,
        // Upstream measured the knee. Trust it — this is the reclamation that
        // matters, and it is monotone, so it cannot oscillate.
        DecoderPressure::Converged { at_workers } => Verdict::Settled(*at_workers),
        // The parser cannot drain what is produced, nothing is runnable, or the
        // stream is over. None of these are fixed by more workers.
        DecoderPressure::ConsumerBound { .. }
        | DecoderPressure::Idle
        | DecoderPressure::Finished => Verdict::Slack,
        // Still classifying the input.
        DecoderPressure::Starting => Verdict::Undecided,
        // `DecoderPressure` is `#[non_exhaustive]`: upstream deferred
        // `SourceBound` because it cannot yet separate source waiting from
        // inflate time reliably enough to promise as a contract. Treat any
        // future variant as slack — the conservative reading, since the
        // classifications we would most expect next are I/O-shaped, and those
        // argue for fewer workers rather than more.
        _ => Verdict::Slack,
    }
}

fn supervise(
    handles: Vec<DecoderHandle>,
    budget: usize,
    stop: Arc<AtomicBool>,
    peak_threads: Arc<AtomicUsize>,
) -> DecodeBudgetReport {
    let n = handles.len();
    let fair = (budget / n).max(1);
    let mut report = DecodeBudgetReport::default();
    let mut converged: Vec<Option<usize>> = vec![None; n];
    let mut last_verdict: Vec<Verdict> = vec![Verdict::Undecided; n];
    let mut score: Vec<i32> = vec![0; n];
    // Remember the last ceiling written so we only call into the decoder when
    // the decision actually changes.
    let mut applied: Vec<usize> = vec![INITIAL_WORKERS_PER_FILE.min(fair.max(1)); n];

    while !stop.load(Ordering::Relaxed) {
        let snapshots: Vec<_> = handles.iter().map(|h| h.stats()).collect();

        let live: usize = snapshots.iter().map(|s| s.spawned_workers).sum();
        let aux: usize = snapshots.iter().map(|s| s.auxiliary_threads).sum();
        let busy: usize = snapshots.iter().map(|s| s.busy_workers).sum();
        report.peak_worker_threads = report.peak_worker_threads.max(live);
        report.peak_auxiliary_threads = report.peak_auxiliary_threads.max(aux);
        report.peak_busy_workers = report.peak_busy_workers.max(busy);
        peak_threads.store(report.peak_worker_threads, Ordering::Relaxed);

        for (i, s) in snapshots.iter().enumerate() {
            if let DecoderPressure::Converged { at_workers } = s.pressure {
                converged[i] = Some(at_workers);
            }
        }

        if snapshots
            .iter()
            .all(|s| matches!(s.pressure, DecoderPressure::Finished))
        {
            break;
        }

        // Accumulate evidence rather than requiring consecutive agreement.
        for (i, s) in snapshots.iter().enumerate() {
            last_verdict[i] = classify(&s.pressure);
            match last_verdict[i] {
                Verdict::Starved => score[i] = (score[i] + 1).min(EVIDENCE_CLAMP),
                Verdict::Slack => score[i] = (score[i] - 1).max(-EVIDENCE_CLAMP),
                _ => {}
            }
        }

        // Grow only what has earned it. Everything a decoder is not currently
        // holding stays unspent, which is the whole point of starting low.
        let held: usize = applied.iter().sum();
        let spare = budget.saturating_sub(held);

        for (i, handle) in handles.iter().enumerate() {
            let configured = snapshots[i].configured_workers;
            if configured == 0 {
                continue;
            }

            let current = applied[i].max(1);
            let mut want = current;

            if score[i] >= EVIDENCE_THRESHOLD && spare > 0 {
                // Multiplicative growth: the knee can be 8x the starting point
                // (Flex wanted 16/file), and linear steps would spend most of a
                // short run ramping. Doubling reaches 16 from 2 in three
                // decisions, roughly a second at this sample interval.
                want = (current * 2).min(current + spare);
                score[i] = 0;
            }
            // Deliberately no `Slack`-driven shrink. It was tried and it
            // oscillates with a ~1 s period: at 4 workers the decoder is
            // `DecoderBound`, so the share grows to 8; at 8 it gets ahead of the
            // parser and reports `ConsumerBound`, so the share halves back to 4;
            // repeat forever. The bug is treating `ConsumerBound` as waste, when
            // it actually means the decoder is comfortably ahead — which is the
            // goal, not a fault.
            //
            // The measured asymmetry says to ratchet: past its knee every
            // workload's curve is flat, so surplus workers cost threads but no
            // time, whereas re-starving costs 22-48%. Growth is therefore
            // one-way within a run, bounded by `budget`, and given back only on
            // `Converged`, where upstream has actually measured the knee rather
            // than us inferring it from one sample.
            // Never exceed what upstream calibration measured as useful. This
            // is monotone reclamation and is the main source of savings.
            if let Some(at) = converged[i] {
                want = want.min(at);
            }
            let want = want.clamp(1, configured);

            // Ignore adjustments too small to matter; they only churn threads.
            if want.abs_diff(applied[i]) >= MIN_ADJUSTMENT || applied[i] == 0 {
                if handle.set_worker_limit(want).is_ok() {
                    tracing::debug!(
                        "decode-budget: file {} {} -> {} ({:?} score {}, spawned {})",
                        i,
                        applied[i],
                        want,
                        last_verdict[i],
                        score[i],
                        snapshots[i].spawned_workers
                    );
                    applied[i] = want;
                }
            }
        }

        std::thread::sleep(SAMPLE_INTERVAL);
    }

    report.converged_workers = converged.iter().filter_map(|c| *c).collect();
    report.unused_budget = budget.saturating_sub(report.peak_worker_threads);
    report
}
