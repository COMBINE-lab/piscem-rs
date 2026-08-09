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

use std::path::Path;

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
    /// Read a policy file, or use the defaults when none was given.
    pub fn load(path: Option<&Path>) -> Result<Self> {
        let Some(path) = path else {
            return Ok(Self::default());
        };
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("could not read thread policy {}", path.display()))?;
        let policy: Self = serde_json::from_str(&text)
            .with_context(|| format!("could not parse thread policy {}", path.display()))?;
        tracing::info!(
            "thread policy from {}: parallel decoding requires {} thread(s) per gzip input",
            path.display(),
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
        let p = ThreadPolicy::load(Some(&f)).unwrap();
        assert_eq!(p.parallel_decode.min_threads_per_stream, 2);
    }

    #[test]
    fn an_empty_object_is_the_default() {
        let dir = tempfile::tempdir().unwrap();
        let f = write(dir.path(), "p.json", "{}");
        assert_eq!(
            ThreadPolicy::load(Some(&f)).unwrap(),
            ThreadPolicy::default()
        );
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
        let err = ThreadPolicy::load(Some(&f)).unwrap_err().to_string();
        assert!(err.contains("could not parse"), "{err}");
    }

    #[test]
    fn a_missing_file_names_itself() {
        let err = ThreadPolicy::load(Some(std::path::Path::new("/nonexistent/p.json")))
            .unwrap_err()
            .to_string();
        assert!(err.contains("/nonexistent/p.json"), "{err}");
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
