//! Memory-budget detection for sshash's minimizer sort.
//!
//! sshash decides between an in-memory minimizer sort and a disk-backed
//! external k-way merge sort by comparing an estimated footprint
//! (`total_kmers * 48` bytes) against [`BuildConfiguration::ram_limit_gib`].
//! Its default is a fixed 8 GiB, which puts the spill point at ~179M k-mers
//! regardless of how much memory the machine actually has.
//!
//! [`resolve_ram_limit_gib`] keeps an explicit `--ram-limit-gib` verbatim and
//! otherwise derives a budget from the memory this *process* may actually use.
//!
//! # Why not just physical RAM
//!
//! `sysinfo` reports total physical memory, not the ceiling enforced on this
//! process. On a scheduler-managed node the two differ by orders of magnitude:
//! a machine with 1.5 TiB installed may confine a job to 32 GiB via a cgroup.
//! Budgeting from physical RAM there would tell sshash to sort in memory and
//! get the job OOM-killed, turning a slow build into a failed one. Container
//! runtimes behave the same way.
//!
//! So the budget is derived from the *minimum* of every ceiling the platform
//! exposes:
//!
//! | source                          | Linux | macOS | Windows |
//! |---------------------------------|-------|-------|---------|
//! | physical RAM (`sysinfo`)        | yes   | yes   | yes     |
//! | `getrlimit(RLIMIT_AS)` (POSIX)  | yes   | yes   | no      |
//! | cgroup v2 / v1 memory limit     | yes   | no    | no      |
//!
//! No single mechanism is portable: `getrlimit` covers Linux and macOS, but
//! schedulers and container runtimes enforce through cgroups rather than
//! rlimits, so cgroups catch a *different* set of cases rather than subsuming
//! them. Each source is compiled only where it exists; a platform offering
//! none degrades to physical RAM, i.e. today's behaviour.

const GIB: u64 = 1 << 30;

/// sshash's own default, and the floor this module will not drop below when a
/// machine has at least this much usable memory.
pub const DEFAULT_RAM_LIMIT_GIB: usize = 8;

/// Fraction of usable memory to hand to the minimizer sort. The build needs
/// headroom for everything else it holds (contig table, MPHF construction), so
/// half is deliberately conservative.
const BUDGET_NUMERATOR: u64 = 1;
const BUDGET_DENOMINATOR: u64 = 2;

/// Total installed physical memory, in bytes.
fn physical_memory_bytes() -> Option<u64> {
    use sysinfo::{MemoryRefreshKind, RefreshKind, System};
    let sys = System::new_with_specifics(
        RefreshKind::nothing().with_memory(MemoryRefreshKind::nothing().with_ram()),
    );
    match sys.total_memory() {
        0 => None,
        n => Some(n),
    }
}

/// The `RLIMIT_AS` (maximum address space) soft limit, when it is finite.
///
/// `RLIMIT_AS` bounds virtual address space rather than resident memory, so it
/// is an upper bound on what the process can allocate, which is exactly the
/// question here. `RLIM_INFINITY` (the common case) yields `None`.
#[cfg(unix)]
fn rlimit_as_bytes() -> Option<u64> {
    // SAFETY: `getrlimit` only writes into the provided `rlimit` struct.
    let mut lim = unsafe { std::mem::zeroed::<libc::rlimit>() };
    let rc = unsafe { libc::getrlimit(libc::RLIMIT_AS, &mut lim) };
    if rc != 0 {
        return None;
    }
    // `rlim_t` is u64 on Linux/macOS but the width is not guaranteed, so widen
    // through u128 rather than assuming.
    let soft = u128::from(lim.rlim_cur);
    if soft == u128::from(libc::RLIM_INFINITY) || soft == 0 {
        None
    } else {
        Some(u64::try_from(soft).unwrap_or(u64::MAX))
    }
}

#[cfg(not(unix))]
fn rlimit_as_bytes() -> Option<u64> {
    None
}

/// The cgroup memory ceiling applied to this process, if any.
///
/// Tries cgroup v2 (`memory.max`) first, then v1 (`memory.limit_in_bytes`).
/// A scheduler nests the job in its own cgroup and the mount-root file usually
/// reads `max`, so the process's own path is resolved from `/proc/self/cgroup`
/// before falling back to the root.
#[cfg(target_os = "linux")]
fn cgroup_limit_bytes() -> Option<u64> {
    use std::path::{Path, PathBuf};

    /// Read a cgroup limit file. `max` / `-1` / absurdly large values mean
    /// "no limit" and yield `None`.
    fn read_limit(path: &Path) -> Option<u64> {
        let raw = std::fs::read_to_string(path).ok()?;
        let trimmed = raw.trim();
        if trimmed.is_empty() || trimmed == "max" || trimmed == "-1" {
            return None;
        }
        let bytes: u64 = trimmed.parse().ok()?;
        // cgroup v1 encodes "unlimited" as a saturated value (commonly
        // u64::MAX rounded down to the page size).
        if bytes == 0 || bytes >= u64::MAX / 2 {
            None
        } else {
            Some(bytes)
        }
    }

    /// The v2 (`0::<path>`) and v1 (`<id>:memory:<path>`) relative paths for
    /// this process, from `/proc/self/cgroup`.
    fn self_cgroup_paths() -> (Option<String>, Option<String>) {
        let Ok(content) = std::fs::read_to_string("/proc/self/cgroup") else {
            return (None, None);
        };
        let (mut v2, mut v1) = (None, None);
        for line in content.lines() {
            // hierarchy-ID : controller-list : cgroup-path
            let mut parts = line.splitn(3, ':');
            let (Some(id), Some(controllers), Some(path)) =
                (parts.next(), parts.next(), parts.next())
            else {
                continue;
            };
            if id == "0" && controllers.is_empty() {
                v2 = Some(path.to_string());
            } else if controllers.split(',').any(|c| c == "memory") {
                v1 = Some(path.to_string());
            }
        }
        (v2, v1)
    }

    /// Candidate files for one hierarchy: the process's own cgroup first, then
    /// each ancestor up to the mount root (an ancestor may carry the limit
    /// even when the leaf does not), then the mount root itself.
    fn candidates(mount: &str, rel: Option<&str>, file: &str) -> Vec<PathBuf> {
        let mut out = Vec::new();
        if let Some(rel) = rel {
            let rel = rel.trim_start_matches('/');
            let mut cur = PathBuf::from(mount);
            if !rel.is_empty() {
                for component in rel.split('/') {
                    cur.push(component);
                    out.push(cur.join(file));
                }
            }
            out.reverse(); // deepest (most specific) first
        }
        out.push(PathBuf::from(mount).join(file));
        out
    }

    let (v2_rel, v1_rel) = self_cgroup_paths();

    // cgroup v2 unified hierarchy.
    for path in candidates("/sys/fs/cgroup", v2_rel.as_deref(), "memory.max") {
        if let Some(bytes) = read_limit(&path) {
            return Some(bytes);
        }
    }

    // cgroup v1: the memory controller is mounted at its own subdirectory.
    // Under a rootless container the namespace may remap paths, in which case
    // the relative path does not resolve and the mount root is used.
    for mount in ["/sys/fs/cgroup/memory", "/sys/fs/cgroup"] {
        for path in candidates(mount, v1_rel.as_deref(), "memory.limit_in_bytes") {
            if let Some(bytes) = read_limit(&path) {
                return Some(bytes);
            }
        }
    }

    None
}

#[cfg(not(target_os = "linux"))]
fn cgroup_limit_bytes() -> Option<u64> {
    None
}

/// Memory this process may actually use: the minimum of every ceiling the
/// platform exposes. `None` when no source reports anything usable.
pub fn usable_memory_bytes() -> Option<u64> {
    [
        physical_memory_bytes(),
        rlimit_as_bytes(),
        cgroup_limit_bytes(),
    ]
    .into_iter()
    .flatten()
    .min()
}

/// Resolve the `ram_limit_gib` handed to sshash.
///
/// * `Some(n)` is used verbatim, including `Some(0)`, which is sshash's
///   "unlimited" (always sort in memory).
/// * `None` auto-detects: half of [`usable_memory_bytes`], floored to whole
///   GiB, raised to [`DEFAULT_RAM_LIMIT_GIB`] when the machine has at least
///   that much usable memory.
///
/// The floor is clamped to what exists. On a machine with less than 8 GiB
/// usable, the fixed 8 GiB default already budgets more than the machine has;
/// spilling to disk there is slower but survivable, being OOM-killed is not.
/// If nothing can be detected, the result is [`DEFAULT_RAM_LIMIT_GIB`], i.e.
/// sshash's own default.
pub fn resolve_ram_limit_gib(requested: Option<usize>) -> usize {
    if let Some(explicit) = requested {
        tracing::info!("sshash minimizer-sort RAM budget: {explicit} GiB (user-specified)");
        return explicit;
    }

    let Some(usable) = usable_memory_bytes() else {
        tracing::info!(
            "sshash minimizer-sort RAM budget: {DEFAULT_RAM_LIMIT_GIB} GiB \
             (could not detect usable memory)"
        );
        return DEFAULT_RAM_LIMIT_GIB;
    };

    let usable_gib = usable / GIB;
    let half_gib = usable * BUDGET_NUMERATOR / BUDGET_DENOMINATOR / GIB;
    let budget = if usable_gib < DEFAULT_RAM_LIMIT_GIB as u64 {
        // Never promise more than the machine has; keep at least 1 GiB so the
        // value stays distinct from sshash's "unlimited" sentinel of 0.
        usable_gib.max(1)
    } else {
        half_gib.max(DEFAULT_RAM_LIMIT_GIB as u64)
    };

    let budget = usize::try_from(budget).unwrap_or(DEFAULT_RAM_LIMIT_GIB);
    tracing::info!(
        "sshash minimizer-sort RAM budget: {budget} GiB (auto-detected from {usable_gib} GiB usable)"
    );
    budget
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn explicit_value_round_trips() {
        assert_eq!(resolve_ram_limit_gib(Some(64)), 64);
        assert_eq!(resolve_ram_limit_gib(Some(1)), 1);
    }

    #[test]
    fn explicit_zero_is_preserved_as_unlimited() {
        // 0 is sshash's "unlimited" sentinel and must not be rewritten.
        assert_eq!(resolve_ram_limit_gib(Some(0)), 0);
    }

    #[test]
    fn auto_never_returns_the_unlimited_sentinel() {
        // Machine-dependent, so assert the invariant rather than a value.
        assert!(resolve_ram_limit_gib(None) > 0);
    }

    #[test]
    fn auto_does_not_exceed_detected_usable_memory() {
        let budget = resolve_ram_limit_gib(None);
        if let Some(usable) = usable_memory_bytes() {
            let usable_gib = usize::try_from(usable / GIB).unwrap_or(usize::MAX);
            assert!(
                budget <= usable_gib.max(1),
                "budget {budget} GiB exceeds {usable_gib} GiB usable"
            );
        } else {
            assert_eq!(budget, DEFAULT_RAM_LIMIT_GIB);
        }
    }

    #[test]
    fn auto_matches_the_default_when_the_machine_is_large_enough() {
        if let Some(usable) = usable_memory_bytes()
            && usable / GIB >= DEFAULT_RAM_LIMIT_GIB as u64
        {
            assert!(resolve_ram_limit_gib(None) >= DEFAULT_RAM_LIMIT_GIB);
        }
    }

    #[test]
    fn usable_memory_does_not_exceed_physical_memory() {
        if let (Some(usable), Some(physical)) = (usable_memory_bytes(), physical_memory_bytes()) {
            assert!(usable <= physical);
        }
    }
}
