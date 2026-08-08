//! Resolve whether gzip input can use the shared parallel decoder.
//!
//! This module no longer predicts a mapping/decode split before the run. Under
//! `auto`, seekable gzip input opens the parallel path and the live
//! `thread-broker` controller divides the execution-slot budget from cumulative
//! pipeline evidence. This module retains only decisions that can be made
//! without measurement: user preference, input seekability, compression type,
//! and the one-thread serial case.
//!
//! [`choose_decoder`] reports request conflicts, applies an explicit
//! [`DecoderPreference`], then applies [`forced_choice`]. If no rule forces a
//! result, it selects the adaptive parallel path. Non-regular inputs are never
//! opened merely to sniff compression: a FIFO or process substitution cannot
//! be rewound, and the parallel decoder cannot use positional reads on it.

use std::path::{Path, PathBuf};

/// One read group: the files a single logical record is assembled from.
///
/// `[r1]` for single-end, `[r1, r2]` for paired, `[r1, barcode, r2]` for scATAC.
pub type ReadGroup = Vec<PathBuf>;

/// What kind of input this is, which decides whether a restart is even possible.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputKind {
    /// A regular file: supports positional reads, and can be re-opened.
    Regular,
    /// A pipe, FIFO, socket or character device. Cannot be re-read, and
    /// `rapidgzip` degrades to sequential decoding on it.
    Stream,
}

/// Classify `path`. Anything that is not a regular file is a [`InputKind::Stream`].
///
/// Mirrors `rapidgzip`'s own `supports_positional_reads`, which gates its
/// parallel path on `file_type().is_file()`. Keeping the same test here means
/// this crate never hands the parallel decoder an input it would silently
/// decode sequentially.
///
/// Reads metadata only, so it is safe to call on a FIFO.
pub fn classify_input(path: &Path) -> InputKind {
    match std::fs::metadata(path) {
        Ok(md) if md.file_type().is_file() => InputKind::Regular,
        // Unreadable metadata is treated as a stream: the conservative answer,
        // since the serial path works on everything.
        _ => InputKind::Stream,
    }
}

/// Split an input set into the regular files and the ones that cannot be rewound.
///
/// The decoder choice is made **per file** — a FIFO among the inputs does not
/// cost the regular files their parallel decoder, because
/// [`crate::io::fastx::open_input`] re-checks each path and falls back on its
/// own. What the whole set decides is only whether the parallel decoder is worth
/// considering at all while preserving which paths must not be opened.
pub fn partition_inputs(paths: &[PathBuf]) -> (Vec<PathBuf>, Vec<PathBuf>) {
    paths
        .iter()
        .cloned()
        .partition(|p| classify_input(p) == InputKind::Regular)
}

/// Which decoder to use, and why.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Decision {
    pub parallel: bool,
    pub reason: Reason,
}

/// Why a decision was reached. Carried so callers can log it; users hitting an
/// unexpectedly slow run should be able to see which branch fired.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Reason {
    /// No input is a regular file. The parallel decoder cannot use positional
    /// reads and would decode sequentially anyway, while costing coordinator
    /// threads.
    NonSeekableInput,
    /// Only one mapping thread, so decode can never be the constraint.
    SingleThreaded,
    /// The user asked for this decoder explicitly.
    UserRequested,
    /// Nothing forces the answer, so the parallel path is opened and the thread
    /// broker sizes it from live evidence.
    Adaptive,
}

/// Decisions that follow from the input alone.
///
/// Returns `None` when the answer should be left to adaptive allocation during
/// the real run.
///
/// These are forcings rather than guesses:
///
/// * With no regular file among the inputs, nothing can be decoded in parallel.
/// * With one mapping thread there is nothing to starve.
///
/// A threads-per-file bound used to live here too, justified as a forcing on
/// the grounds that serial supply exceeded what so few threads could consume
/// "at the fastest per-thread rate ever measured". That reasoning embeds a
/// measured consumer rate, so it was empirical rather than logical. Consumer
/// cost is now observed by the live broker rather than predicted here.
///
/// A user's explicit `--decoder serial` is also honoured unconditionally, but
/// earlier: see [`preference_choice`], which runs before this.
pub fn forced_choice(kind: InputKind, map_threads: usize) -> Option<Decision> {
    if kind == InputKind::Stream {
        return Some(Decision {
            parallel: false,
            reason: Reason::NonSeekableInput,
        });
    }
    if map_threads <= 1 {
        return Some(Decision {
            parallel: false,
            reason: Reason::SingleThreaded,
        });
    }
    None
}

/// A user's `--decoder` request.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum DecoderPreference {
    /// Decide by forcing rules, then measurement.
    #[default]
    Auto,
    /// Always use the serial (libz-rs) path.
    Serial,
    /// Use the parallel decoder wherever it is possible at all, with an
    /// optional per-file worker ceiling.
    ///
    /// Still yields to [`Reason::NonSeekableInput`]: on a pipe the parallel
    /// decoder degrades to sequential decoding, so honouring the request
    /// literally would spend coordinator threads to obtain nothing. A
    /// preference cannot make an input seekable.
    Parallel { workers_per_file: Option<usize> },
}

impl DecoderPreference {
    /// Resolve `s` as given on a command line: `auto`, `serial`, `parallel`, or
    /// `parallel=N`.
    pub fn parse(s: &str) -> Result<Self, String> {
        let s = s.trim();
        let (head, tail) = match s.split_once('=') {
            Some((h, t)) => (h, Some(t)),
            None => (s, None),
        };
        match (head.to_ascii_lowercase().as_str(), tail) {
            ("auto", None) => Ok(Self::Auto),
            ("serial" | "libz" | "libdeflate", None) => Ok(Self::Serial),
            ("parallel" | "rapidgzip", None) => Ok(Self::Parallel {
                workers_per_file: None,
            }),
            ("parallel" | "rapidgzip", Some(n)) => n
                .parse::<usize>()
                .map(|w| Self::Parallel {
                    workers_per_file: Some(w.max(1)),
                })
                .map_err(|_| format!("`{n}` is not a worker count")),
            _ => Err(format!(
                "expected `auto`, `serial`, `parallel`, or `parallel=N`, got `{s}`"
            )),
        }
    }
}

/// How the input is compressed, as far as the decoder choice cares.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputCompression {
    /// gzip — the only format either decoder path competes over.
    Gzip,
    /// Plain text, or a codec (zstd, bz2) that only the serial path handles.
    NotGzip,
    /// Not determined, because determining it would consume the input.
    Unknown,
}

/// Sniff whether `path` is gzip.
///
/// Only ever called for [`InputKind::Regular`]. Sniffing a stream would consume
/// the bytes it read — exactly the double-open that made FIFO input hang — so a
/// non-regular input reports [`InputCompression::Unknown`] rather than being
/// opened speculatively. **Keep it that way.**
pub fn detect_compression(path: &Path, kind: InputKind) -> InputCompression {
    use std::io::Read;
    if kind != InputKind::Regular {
        return InputCompression::Unknown;
    }
    let Ok(mut f) = std::fs::File::open(path) else {
        return InputCompression::Unknown;
    };
    let mut magic = [0u8; 2];
    match f.read_exact(&mut magic) {
        Ok(()) if magic == [0x1f, 0x8b] => InputCompression::Gzip,
        Ok(()) => InputCompression::NotGzip,
        Err(_) => InputCompression::Unknown,
    }
}

/// A request that cannot be carried out as asked.
///
/// Surfaced so the user is told their flag was overridden, rather than quietly
/// getting something else. Someone who passes `--decoder parallel` and sees no
/// speedup deserves to know why — and, when only *part* of the input is
/// affected, which part.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PreferenceConflict {
    /// `parallel` asked for, but some inputs cannot support positional reads.
    /// Carries them, since naming the offender is the whole point on a mixed set.
    ParallelOnStream { streams: Vec<PathBuf>, total: usize },
    /// A decoder named for an input neither decoder path applies to.
    NotGzipInput,
}

impl std::fmt::Display for PreferenceConflict {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ParallelOnStream { streams, total } => {
                let names = streams
                    .iter()
                    .map(|p| p.display().to_string())
                    .collect::<Vec<_>>()
                    .join(", ");
                if streams.len() == *total {
                    write!(
                        f,
                        "--decoder parallel was requested, but no input is a regular file \
                         (these are pipes, FIFOs or process substitutions): {names}. The \
                         parallel decoder needs positional reads and would decode \
                         sequentially anyway, so the serial decoder is being used. To get \
                         parallel decoding, pass real files rather than streams."
                    )
                } else {
                    write!(
                        f,
                        "--decoder parallel was requested, but {} of {total} inputs are not \
                         regular files and will use the serial decoder: {names}. The \
                         remaining inputs are unaffected and still decode in parallel.",
                        streams.len()
                    )
                }
            }
            Self::NotGzipInput => write!(
                f,
                "--decoder was set, but the input is not gzip-compressed, so the setting \
                 has no effect: only gzip input has both a serial and a parallel decoder. \
                 Plain, zstd and bzip2 inputs are always read by the serial path."
            ),
        }
    }
}

/// Report any way in which `pref` cannot be honoured for this input set.
///
/// `non_regular` is the subset of inputs that cannot be rewound, and `total` the
/// number of inputs, so a partial downgrade can be described as one.
pub fn preference_conflict(
    pref: DecoderPreference,
    non_regular: &[PathBuf],
    total: usize,
    compression: InputCompression,
) -> Option<PreferenceConflict> {
    if pref == DecoderPreference::Auto {
        return None;
    }
    if compression == InputCompression::NotGzip {
        return Some(PreferenceConflict::NotGzipInput);
    }
    if matches!(pref, DecoderPreference::Parallel { .. }) && !non_regular.is_empty() {
        return Some(PreferenceConflict::ParallelOnStream {
            streams: non_regular.to_vec(),
            total,
        });
    }
    None
}

/// Apply an explicit preference, if it settles the question.
///
/// Returns `None` for [`DecoderPreference::Auto`], and for a parallel request on
/// an input set with nothing seekable in it, so the caller falls through to the
/// normal path in both cases.
///
/// `kind` is [`InputKind::Regular`] when *any* input is a regular file: the
/// decision is run-level ("is the parallel decoder worth considering"), and the
/// per-file downgrade happens later in `open_input`.
pub fn preference_choice(pref: DecoderPreference, kind: InputKind) -> Option<Decision> {
    match pref {
        DecoderPreference::Auto => None,
        DecoderPreference::Serial => Some(Decision {
            parallel: false,
            reason: Reason::UserRequested,
        }),
        DecoderPreference::Parallel { .. } if kind == InputKind::Stream => Some(Decision {
            parallel: false,
            reason: Reason::NonSeekableInput,
        }),
        DecoderPreference::Parallel { .. } => Some(Decision {
            parallel: true,
            reason: Reason::UserRequested,
        }),
    }
}

/// Resolve the decoder for a whole run.
///
/// `groups` is the run's input in the order `paraseq` receives it — one entry
/// per logical record, holding that record's files (`[r1]`, `[r1, r2]`, or
/// `[r1, barcode, r2]`).
///
/// # This no longer measures anything
///
/// It used to run a miniature pipeline before the real one, to predict how the
/// thread budget should be split. That is gone. Both sides of the split are now
/// resizable while the run is in flight, so the split is *converged to* rather
/// than predicted — see [`crate::io::broker`] and the `thread-broker` crate.
///
/// The prediction was worth removing on its own merits. Measured against the
/// optimum on our own workloads, the closed-form supply/demand model was **60%
/// off**, and a two-point measurement was 7-38% off depending on where the
/// second point was taken, with no placement good on both. A plain constant beat
/// both. Meanwhile the probe cost 215-500 ms before every run and could not, in
/// principle, notice a workload that changes partway through.
///
/// What is left here is only what follows from the input without measurement:
/// whether anything *can* be decoded in parallel, and what the user asked for.
pub fn choose_decoder(
    groups: &[ReadGroup],
    map_threads: usize,
    pref: DecoderPreference,
) -> Decision {
    let all: Vec<PathBuf> = groups.iter().flatten().cloned().collect();
    if all.is_empty() {
        return Decision {
            parallel: false,
            reason: Reason::NonSeekableInput,
        };
    }
    let (regular, streams) = partition_inputs(&all);

    // Sniffing is safe only on a regular file, so ask one of those.
    let compression = regular
        .first()
        .map(|p| detect_compression(p, InputKind::Regular))
        .unwrap_or(InputCompression::Unknown);

    // Volume matched to what was asked for. An explicit request that cannot be
    // honoured is a WARN, which the subscriber shows without `RUST_LOG` -- a
    // message saying a flag was overridden is useless if only visible to
    // someone who already suspected a problem. Under `auto` nothing was
    // overridden, so the same fact is merely informative.
    if let Some(conflict) = preference_conflict(pref, &streams, all.len(), compression) {
        tracing::warn!("{conflict}");
    } else if !streams.is_empty() {
        tracing::info!(
            "{} of {} inputs are not regular files and will use the serial decoder: {}",
            streams.len(),
            all.len(),
            streams
                .iter()
                .map(|p| p.display().to_string())
                .collect::<Vec<_>>()
                .join(", ")
        );
    }

    // Run-level: the parallel decoder is worth considering if anything is
    // seekable. Per-file downgrade happens in `open_input`.
    let kind = if regular.is_empty() {
        InputKind::Stream
    } else {
        InputKind::Regular
    };

    if let Some(chosen) = preference_choice(pref, kind) {
        tracing::info!(
            "decoder selected by request: {} ({:?})",
            if chosen.parallel {
                "parallel"
            } else {
                "serial"
            },
            chosen.reason
        );
        return chosen;
    }
    if let Some(forced) = forced_choice(kind, map_threads) {
        tracing::debug!("decoder choice forced: {:?}", forced.reason);
        return forced;
    }

    // Nothing forces the answer, so open the parallel path and let the broker
    // size it. Costing a run that does not need decode threads is bounded: the
    // broker reclaims them to the floor within a few sampling intervals, and an
    // input the decoder cannot parallelise reports itself inelastic and is left
    // alone entirely.
    Decision {
        parallel: true,
        reason: Reason::Adaptive,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn streams_never_take_the_parallel_path() {
        // The FIFO case: rapidgzip would decode it sequentially anyway, and a
        // compression sniff would consume bytes.
        let d = forced_choice(InputKind::Stream, 64).unwrap();
        assert!(!d.parallel);
        assert_eq!(d.reason, Reason::NonSeekableInput);
    }

    #[test]
    fn single_thread_is_forced_serial() {
        let d = forced_choice(InputKind::Regular, 1).unwrap();
        assert!(!d.parallel);
        assert_eq!(d.reason, Reason::SingleThreaded);
    }

    /// Only genuinely logical cases are forced. A threads-per-file bound used
    /// to live here; it encoded consumer cost as a constant, which the live
    /// broker now measures directly.
    #[test]
    fn thread_ratios_are_adapted_not_forced() {
        assert!(forced_choice(InputKind::Regular, 16).is_none());
        assert!(forced_choice(InputKind::Regular, 4).is_none());
        assert!(forced_choice(InputKind::Regular, 32).is_none());
    }

    #[test]
    fn preference_parses_its_command_line_forms() {
        use DecoderPreference as P;
        assert_eq!(P::parse("auto").unwrap(), P::Auto);
        assert_eq!(P::parse("serial").unwrap(), P::Serial);
        assert_eq!(
            P::parse("PARALLEL").unwrap(),
            P::Parallel {
                workers_per_file: None
            }
        );
        assert_eq!(
            P::parse("parallel=8").unwrap(),
            P::Parallel {
                workers_per_file: Some(8)
            }
        );
        // A zero worker request is a parallel request, not a serial one.
        assert_eq!(
            P::parse("parallel=0").unwrap(),
            P::Parallel {
                workers_per_file: Some(1)
            }
        );
        assert!(P::parse("sometimes").is_err());
        assert!(P::parse("parallel=lots").is_err());
    }

    #[test]
    fn a_preference_cannot_make_a_pipe_seekable() {
        // Honouring `parallel` with nothing seekable would spend coordinator
        // threads to get sequential decoding, so the forcing still wins.
        let d = preference_choice(
            DecoderPreference::Parallel {
                workers_per_file: Some(8),
            },
            InputKind::Stream,
        )
        .unwrap();
        assert!(!d.parallel);
        assert_eq!(d.reason, Reason::NonSeekableInput);

        // With a regular file present it is honoured.
        let d = preference_choice(
            DecoderPreference::Parallel {
                workers_per_file: None,
            },
            InputKind::Regular,
        )
        .unwrap();
        assert!(d.parallel);
        assert_eq!(d.reason, Reason::UserRequested);
    }

    #[test]
    fn auto_defers_to_the_normal_path() {
        assert!(preference_choice(DecoderPreference::Auto, InputKind::Regular).is_none());
        assert!(preference_choice(DecoderPreference::Auto, InputKind::Stream).is_none());
    }

    #[test]
    fn conflicts_are_reported_only_when_real() {
        use DecoderPreference as P;
        let fifo = vec![PathBuf::from("/tmp/a.fifo")];
        // Auto never conflicts, whatever the input.
        assert!(preference_conflict(P::Auto, &fifo, 2, InputCompression::NotGzip).is_none());
        // Parallel with a pipe in the set.
        assert!(matches!(
            preference_conflict(
                P::Parallel {
                    workers_per_file: None
                },
                &fifo,
                2,
                InputCompression::Unknown
            ),
            Some(PreferenceConflict::ParallelOnStream { .. })
        ));
        // Any explicit decoder on non-gzip input, including `serial`, since the
        // flag implies a choice that does not exist there.
        assert_eq!(
            preference_conflict(P::Serial, &[], 1, InputCompression::NotGzip),
            Some(PreferenceConflict::NotGzipInput)
        );
        // Honourable request: nothing to say.
        assert!(
            preference_conflict(
                P::Parallel {
                    workers_per_file: Some(4)
                },
                &[],
                2,
                InputCompression::Gzip
            )
            .is_none()
        );
    }

    /// A mixed set says which inputs were downgraded and that the rest were not,
    /// because "parallel was requested but you got serial" is actively
    /// misleading when only one file of eight was affected.
    #[test]
    fn a_partial_downgrade_says_so_and_names_names() {
        let c = PreferenceConflict::ParallelOnStream {
            streams: vec![PathBuf::from("/tmp/r2.fifo")],
            total: 4,
        };
        let msg = c.to_string();
        assert!(msg.contains("/tmp/r2.fifo"), "{msg}");
        assert!(msg.contains("1 of 4"), "{msg}");
        assert!(msg.contains("unaffected"), "{msg}");

        // When nothing is seekable it is a whole-run downgrade, worded as one.
        let c = PreferenceConflict::ParallelOnStream {
            streams: vec![PathBuf::from("/tmp/r1.fifo")],
            total: 1,
        };
        let msg = c.to_string();
        assert!(msg.contains("no input is a regular file"), "{msg}");
        assert!(!msg.contains("unaffected"), "{msg}");
    }

    /// The guard that matters most: decoder selection must never open a FIFO,
    /// because the bytes it reads are gone from the run.
    #[test]
    fn an_all_stream_set_is_forced_before_any_read() {
        let kind = InputKind::Stream;
        // Both the preference path and the forcing path refuse.
        assert!(
            !preference_choice(
                DecoderPreference::Parallel {
                    workers_per_file: Some(4)
                },
                kind
            )
            .unwrap()
            .parallel
        );
        assert!(!forced_choice(kind, 32).unwrap().parallel);
    }
}

/// Tests for the guard that matters most: decoder selection must never read a
/// non-rewindable input.
///
/// Compression sniffing reads a prefix and the run then re-opens from the start.
/// On a FIFO those bytes are gone, so sniffing one would silently drop reads.
///
/// These tests exploit a property of FIFOs to make the assertion airtight:
/// opening one for reading **blocks until a writer appears**. Every FIFO below
/// is created with no writer, so any code that opens one deadlocks and the test
/// times out. A test that completes is proof the path was never opened.
#[cfg(all(test, unix))]
mod non_rewindable_input_tests {
    use super::*;
    /// A FIFO with no writer. Reading it blocks forever, which is the point.
    fn fifo(dir: &std::path::Path, name: &str) -> PathBuf {
        let path = dir.join(name);
        let ok = std::process::Command::new("mkfifo")
            .arg(&path)
            .status()
            .expect("mkfifo should be available on unix");
        assert!(ok.success(), "mkfifo failed for {}", path.display());
        path
    }

    fn tmpdir(tag: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!(
            "piscem-fifo-{tag}-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    fn gzip_fixture(dir: &std::path::Path, name: &str) -> PathBuf {
        // A real gzip FASTQ, small but valid, so the regular-file path is
        // genuinely exercised rather than erroring out early.
        let path = dir.join(name);
        let mut body = Vec::new();
        for i in 0..64 {
            body.extend_from_slice(
                format!("@r{i}\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n").as_bytes(),
            );
        }
        let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        std::io::Write::write_all(&mut enc, &body).unwrap();
        std::fs::write(&path, enc.finish().unwrap()).unwrap();
        path
    }

    /// Metadata only -- must not block even on a writer-less FIFO.
    #[test]
    fn classification_never_opens_the_input() {
        let d = tmpdir("classify");
        let f = fifo(&d, "r1.fifo");
        assert_eq!(classify_input(&f), InputKind::Stream);
        // The compression sniff *does* open, so it must refuse on a stream.
        assert_eq!(
            detect_compression(&f, InputKind::Stream),
            InputCompression::Unknown
        );
        let _ = std::fs::remove_dir_all(&d);
    }

    /// Every input is a pipe: forced serial without opening it.
    #[test]
    fn an_all_fifo_run_is_forced_before_anything_is_read() {
        let d = tmpdir("allfifo");
        let groups = vec![vec![fifo(&d, "a.fifo"), fifo(&d, "b.fifo")]];
        let cal = choose_decoder(&groups, 32, DecoderPreference::Auto);
        assert!(!cal.parallel);
        assert_eq!(cal.reason, Reason::NonSeekableInput);
        let _ = std::fs::remove_dir_all(&d);
    }

    /// `--decoder parallel` cannot make a pipe seekable, and must not try.
    #[test]
    fn an_explicit_parallel_request_still_refuses_to_read_a_pipe() {
        let d = tmpdir("explicit");
        let groups = vec![vec![fifo(&d, "a.fifo")]];
        let cal = choose_decoder(
            &groups,
            32,
            DecoderPreference::Parallel {
                workers_per_file: Some(8),
            },
        );
        assert!(!cal.parallel);
        assert_eq!(cal.reason, Reason::NonSeekableInput);
        let _ = std::fs::remove_dir_all(&d);
    }

    /// A group with one regular file and one pipe retains parallel decoding for
    /// the regular file without ever inspecting the pipe's contents.
    #[test]
    fn a_split_group_is_skipped_rather_than_half_read() {
        let d = tmpdir("split");
        let groups = vec![vec![gzip_fixture(&d, "r1.fq.gz"), fifo(&d, "r2.fifo")]];
        let cal = choose_decoder(&groups, 32, DecoderPreference::Auto);
        // Adaptive -- seekable gzip exists and no data was read to decide.
        assert_eq!(cal.reason, Reason::Adaptive);
        // But not forced serial either: the regular file is still a candidate.
        assert!(cal.parallel);
        let _ = std::fs::remove_dir_all(&d);
    }

    /// Mixed groups: regular input remains adaptive and the FIFO stays untouched.
    #[test]
    fn a_regular_group_is_adaptive_while_a_fifo_group_is_left_alone() {
        let d = tmpdir("mixed");
        let groups = vec![
            vec![gzip_fixture(&d, "a1.fq.gz"), gzip_fixture(&d, "a2.fq.gz")],
            vec![fifo(&d, "b1.fifo"), fifo(&d, "b2.fifo")],
        ];
        // Reaching the next line at all is the real assertion: opening either
        // FIFO would block forever and this would never return.
        let cal = choose_decoder(&groups, 8, DecoderPreference::Auto);
        // The FIFO group must not have forced the run serial -- half the input
        // is perfectly seekable.
        assert_ne!(cal.reason, Reason::NonSeekableInput);
        let _ = std::fs::remove_dir_all(&d);
    }
}
