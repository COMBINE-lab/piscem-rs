//! User-supplied overrides for thread and decoder policy.
//!
//! Every value here has a measured default that suits the workloads piscem was
//! validated on. The file exists for the cases those measurements do not cover —
//! an unusual decode-to-map ratio, a filesystem with very different read
//! characteristics, or an experiment that wants to force a decision.
//!
//! # Format
//!
//! ```json
//! {
//!   "parallel_decode": {
//!     "min_threads_per_stream": 8
//!   }
//! }
//! ```
//!
//! Every field is optional and falls back to its default, so a file need only
//! name what it changes. Unknown fields are a hard error rather than a silent
//! no-op: a policy file that looks like it is doing something while doing
//! nothing is worse than one that refuses to load.
//!
//! # Path or inline JSON
//!
//! The same value may be given as a path to a JSON file *or* as an inline JSON
//! document. A value whose first non-whitespace character is `{` is parsed
//! directly; anything else is treated as a filesystem path. This lets a caller
//! that constructs a run programmatically (e.g. `simpleaf`) pass a policy without
//! materialising a temporary file, while a hand-written file path keeps working
//! unchanged.

use anyhow::{Context, Result};

/// Thread and decoder policy for a run.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct ThreadPolicy {
    /// When the parallel gzip decoder is engaged at all.
    ///
    /// This is the *engagement* decision, not the split. Once engaged, the
    /// thread broker sizes decode against mapping from live measurement; this
    /// only decides whether that machinery is worth starting.
    pub parallel_decode: thread_broker::EngagementPolicy,
}

impl ThreadPolicy {
    /// Resolve a policy from a `--thread-policy` argument, or use the defaults
    /// when none was given.
    ///
    /// `spec` is either an inline JSON document (first non-whitespace character
    /// `{`) or a path to a JSON file. Both routes parse to the same
    /// [`ThreadPolicy`] and reject unknown fields identically; the only
    /// difference is where the bytes come from and how a parse error names its
    /// source.
    pub fn load(spec: Option<&str>) -> Result<Self> {
        let Some(spec) = spec else {
            return Ok(Self::default());
        };
        let inline = spec.trim_start().starts_with('{');
        let text = if inline {
            spec.to_string()
        } else {
            std::fs::read_to_string(spec)
                .with_context(|| format!("could not read thread policy file {spec}"))?
        };
        let policy: Self = serde_json::from_str(&text).with_context(|| {
            if inline {
                format!("could not parse inline thread policy `{spec}`")
            } else {
                format!("could not parse thread policy file {spec}")
            }
        })?;
        tracing::info!(
            "thread policy from {}: parallel decoding requires {} thread(s) per gzip input",
            if inline {
                "inline JSON".to_string()
            } else {
                spec.to_string()
            },
            policy.parallel_decode.min_threads_per_stream,
        );
        Ok(policy)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write(dir: &std::path::Path, name: &str, body: &str) -> std::path::PathBuf {
        let p = dir.join(name);
        std::fs::write(&p, body).unwrap();
        p
    }

    /// A path is passed to `load` the same way the CLI passes it: as the string
    /// the user typed.
    fn load_path(p: &std::path::Path) -> Result<ThreadPolicy> {
        ThreadPolicy::load(Some(p.to_str().unwrap()))
    }

    #[test]
    fn absent_path_yields_the_measured_defaults() {
        let p = ThreadPolicy::load(None).unwrap();
        assert_eq!(p, ThreadPolicy::default());
        assert_eq!(p.parallel_decode.min_threads_per_stream, 8);
    }

    #[test]
    fn a_partial_file_overrides_only_what_it_names() {
        let dir = tempfile::tempdir().unwrap();
        let f = write(
            dir.path(),
            "p.json",
            r#"{"parallel_decode": {"min_threads_per_stream": 2}}"#,
        );
        let p = load_path(&f).unwrap();
        assert_eq!(p.parallel_decode.min_threads_per_stream, 2);
    }

    #[test]
    fn an_empty_object_is_the_default() {
        let dir = tempfile::tempdir().unwrap();
        let f = write(dir.path(), "p.json", "{}");
        assert_eq!(load_path(&f).unwrap(), ThreadPolicy::default());
    }

    /// A typo that silently does nothing is the failure mode a policy file is
    /// most likely to have, and the hardest to notice: the run looks configured
    /// and behaves as if it were not.
    #[test]
    fn a_misspelled_field_is_rejected_rather_than_ignored() {
        let dir = tempfile::tempdir().unwrap();
        let f = write(
            dir.path(),
            "p.json",
            r#"{"parallel_decode": {"min_threads_per_input": 2}}"#,
        );
        let err = load_path(&f).unwrap_err().to_string();
        assert!(err.contains("could not parse"), "{err}");
    }

    #[test]
    fn a_missing_file_names_itself() {
        let err = ThreadPolicy::load(Some("/nonexistent/p.json"))
            .unwrap_err()
            .to_string();
        assert!(err.contains("/nonexistent/p.json"), "{err}");
    }

    /// The `simpleaf` route: the policy arrives inline, never touching the disk.
    #[test]
    fn inline_json_is_parsed_directly() {
        let p = ThreadPolicy::load(Some(
            r#"{"parallel_decode": {"min_threads_per_stream": 4}}"#,
        ))
        .unwrap();
        assert_eq!(p.parallel_decode.min_threads_per_stream, 4);
    }

    /// Leading whitespace still selects the inline branch, so a value assembled
    /// with a stray space is not mistaken for a path.
    #[test]
    fn inline_json_with_leading_whitespace_is_still_inline() {
        let p = ThreadPolicy::load(Some(
            r#"  {"parallel_decode": {"min_threads_per_stream": 4}}"#,
        ))
        .unwrap();
        assert_eq!(p.parallel_decode.min_threads_per_stream, 4);
    }

    /// The equivalence the plan gates on: an inline document and the same JSON
    /// in a file resolve to an identical policy.
    #[test]
    fn inline_and_file_resolve_identically() {
        let body = r#"{"parallel_decode": {"min_threads_per_stream": 4}}"#;
        let dir = tempfile::tempdir().unwrap();
        let f = write(dir.path(), "p.json", body);
        assert_eq!(
            ThreadPolicy::load(Some(body)).unwrap(),
            load_path(&f).unwrap()
        );
    }

    /// Inline JSON gets the same `deny_unknown_fields` treatment as a file.
    #[test]
    fn a_misspelled_inline_field_is_rejected() {
        let err = ThreadPolicy::load(Some(
            r#"{"parallel_decode": {"min_threads_per_input": 2}}"#,
        ))
        .unwrap_err()
        .to_string();
        assert!(err.contains("could not parse inline"), "{err}");
    }

    /// Malformed inline JSON names itself as inline rather than as a path.
    #[test]
    fn malformed_inline_json_is_named_inline() {
        let err = ThreadPolicy::load(Some(r#"{"parallel_decode": "#))
            .unwrap_err()
            .to_string();
        assert!(err.contains("inline"), "{err}");
    }

    #[test]
    fn a_policy_round_trips() {
        let doc = ThreadPolicy {
            parallel_decode: thread_broker::EngagementPolicy {
                min_threads_per_stream: 3,
            },
        };
        let text = serde_json::to_string(&doc).unwrap();
        assert_eq!(serde_json::from_str::<ThreadPolicy>(&text).unwrap(), doc);
    }

    #[test]
    fn zero_always_engages_and_no_input_never_does() {
        assert!(thread_broker::EngagementPolicy::always().engages(1, 64));
        // Nothing for a parallel producer to do, so it could only cost threads.
        assert!(!ThreadPolicy::default().parallel_decode.engages(1024, 0));
    }
}

#[cfg(test)]
mod show {
    use super::*;
    #[test]
    #[ignore = "prints the canonical policy document"]
    fn canonical_form() {
        println!(
            "{}",
            serde_json::to_string_pretty(&ThreadPolicy::default()).unwrap()
        );
    }
}
