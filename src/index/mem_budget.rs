//! How much memory the sshash minimizer sort may use.
//!
//! # Why this is not just "half the RAM"
//!
//! sshash chooses its minimizer-sort strategy from a RAM budget: when the
//! estimated footprint (`num_kmers * 48` bytes) exceeds it, the sort spills to
//! a disk-backed external k-way merge; otherwise it sorts in memory. The
//! historical default is 8 GiB, which puts the spill point at ~179 M k-mers —
//! so a human genome index (chr1+chr2 alone estimates 17.7 GiB) spills even on
//! a machine with a terabyte free.
//!
//! Scaling the budget to the machine is therefore worth doing, but "detect
//! physical RAM and take half" is actively dangerous in the environments this
//! tool actually runs in. A scheduler-managed node commonly reports 1.5 TiB of
//! physical RAM while the job is confined to 32 GiB by a cgroup; budgeting
//! 755 GiB there tells sshash to sort in memory and the job is killed. Physical
//! RAM is an upper bound on what exists, not on what this process may use.
//!
//! So take the **minimum of every constraint the platform can show us**:
//!
//! | source | Linux | macOS | Windows |
//! |---|---|---|---|
//! | physical RAM (`sysinfo`) | yes | yes | yes |
//! | `getrlimit(RLIMIT_AS)` (POSIX) | yes | yes | no |
//! | cgroup v2 / v1 memory limit | yes | no | no |
//!
//! There is no single portable mechanism. `getrlimit` is POSIX and covers Linux
//! and macOS, but schedulers and container runtimes enforce through cgroups
//! rather than rlimits, so it catches a different set of cases rather than
//! subsuming them. Each source is compiled only where it exists; on a platform
//! that offers none, the answer degrades to physical RAM, which is the current
//! behaviour.

use tracing::{debug, info};

/// Budget floor, in GiB.
///
/// sshash's historical default. Auto-detection never resolves below this on a
/// machine that can afford it, so no existing build gets a smaller budget than
/// it has today.
const MIN_RAM_LIMIT_GIB: usize = 8;

/// Fraction of the visible limit to hand to the sort.
///
/// The sort is not the only live allocation during a build — the contig table,
/// the reference sequences and the MPHF partitions are all resident — so half
/// leaves room for the rest rather than inviting the OOM killer at the moment
/// the sort peaks.
const BUDGET_FRACTION: f64 = 0.5;

const BYTES_PER_GIB: u64 = 1024 * 1024 * 1024;

/// Total physical memory, in bytes.
fn physical_ram_bytes() -> Option<u64> {
    let mut sys = sysinfo::System::new();
    sys.refresh_memory();
    match sys.total_memory() {
        0 => None,
        n => Some(n),
    }
}

/// Address-space ceiling from `getrlimit(RLIMIT_AS)`, if one is set.
///
/// POSIX, so this is the one limit visible on both Linux and macOS. Usually
/// `RLIM_INFINITY`; some sites do set it, and when they do it binds this
/// process directly.
#[cfg(unix)]
fn rlimit_as_bytes() -> Option<u64> {
    // SAFETY: `getrlimit` writes a fully-initialized `rlimit` on success, and
    // we only read the struct when it returns 0.
    let mut lim = unsafe { std::mem::zeroed::<libc::rlimit>() };
    let rc = unsafe { libc::getrlimit(libc::RLIMIT_AS, &mut lim) };
    if rc != 0 || lim.rlim_cur == libc::RLIM_INFINITY {
        return None;
    }
    Some(lim.rlim_cur as u64)
}

#[cfg(not(unix))]
fn rlimit_as_bytes() -> Option<u64> {
    None
}

/// Memory ceiling imposed by the process's cgroup, if any.
///
/// Handles cgroup v2 (`memory.max`) and v1 (`memory.limit_in_bytes`). Reads the
/// process's *own* cgroup path from `/proc/self/cgroup` before falling back to
/// the mount root, because a scheduler places the job in a nested cgroup and
/// the root file usually reads `max`.
///
/// This is the case that matters: SLURM, LSF and container runtimes all limit
/// memory here rather than through rlimits.
#[cfg(target_os = "linux")]
fn cgroup_limit_bytes() -> Option<u64> {
    use std::fs;

    fn parse(s: &str) -> Option<u64> {
        let s = s.trim();
        if s == "max" {
            return None;
        }
        let v: u64 = s.parse().ok()?;
        // cgroup v1 spells "unlimited" as a huge sentinel rather than a word.
        if v >= i64::MAX as u64 / 2 { None } else { Some(v) }
    }

    // The process's own cgroup path, e.g. "0::/slurm/uid_1000/job_123".
    let own_path = fs::read_to_string("/proc/self/cgroup").ok().and_then(|s| {
        s.lines()
            .find_map(|l| l.strip_prefix("0::").map(str::to_string))
    });

    let mut candidates = Vec::new();
    if let Some(p) = own_path {
        let p = p.trim_start_matches('/');
        if !p.is_empty() {
            candidates.push(format!("/sys/fs/cgroup/{p}/memory.max"));
        }
    }
    candidates.push("/sys/fs/cgroup/memory.max".to_string());
    candidates.push("/sys/fs/cgroup/memory/memory.limit_in_bytes".to_string());

    candidates
        .iter()
        .filter_map(|p| fs::read_to_string(p).ok())
        .find_map(|s| parse(&s))
}

#[cfg(not(target_os = "linux"))]
fn cgroup_limit_bytes() -> Option<u64> {
    None
}

/// Memory this process may actually use, in bytes: the tightest limit visible.
fn usable_memory_bytes() -> Option<u64> {
    let physical = physical_ram_bytes();
    let rlimit = rlimit_as_bytes();
    let cgroup = cgroup_limit_bytes();

    debug!(
        "memory limits: physical={:?} rlimit_as={:?} cgroup={:?}",
        physical, rlimit, cgroup
    );

    [physical, rlimit, cgroup].into_iter().flatten().min()
}

/// Resolve the RAM budget, in GiB, handed to sshash's minimizer sort.
///
/// * `Some(0)` keeps sshash's "unlimited" semantics — always sort in memory.
/// * `Some(n)` is used verbatim; an explicit request is never second-guessed.
/// * `None` auto-detects.
///
/// Auto-detection is `max(MIN_RAM_LIMIT_GIB, usable/2)`, except on a machine
/// whose usable memory is *below* the floor, where it yields the usable amount
/// instead. That last case is a deliberate change from the fixed 8 GiB default,
/// which on a 6 GiB machine budgets more memory than exists; spilling to disk
/// there is slower but survivable, and being OOM-killed is not.
pub fn resolve_ram_limit_gib(requested: Option<usize>) -> usize {
    if let Some(n) = requested {
        if n == 0 {
            info!("sshash minimizer-sort RAM budget: unlimited (requested)");
        } else {
            info!("sshash minimizer-sort RAM budget: {n} GiB (requested)");
        }
        return n;
    }

    let Some(usable) = usable_memory_bytes() else {
        info!(
            "sshash minimizer-sort RAM budget: {MIN_RAM_LIMIT_GIB} GiB \
             (could not determine available memory)"
        );
        return MIN_RAM_LIMIT_GIB;
    };

    let usable_gib = (usable / BYTES_PER_GIB) as usize;
    let half_gib = ((usable as f64 * BUDGET_FRACTION) as u64 / BYTES_PER_GIB) as usize;

    let budget = if half_gib >= MIN_RAM_LIMIT_GIB {
        half_gib
    } else {
        // Never budget more than the machine has, even to reach the floor.
        MIN_RAM_LIMIT_GIB.min(usable_gib.max(1))
    };

    info!(
        "sshash minimizer-sort RAM budget: {budget} GiB \
         (auto-detected; {usable_gib} GiB usable)"
    );
    budget
}

#[cfg(test)]
mod tests {
    use super::*;

    /// An explicit request is authoritative, including the `0` that means
    /// "unlimited" — auto-detection must not quietly override either.
    #[test]
    fn explicit_requests_pass_through() {
        assert_eq!(resolve_ram_limit_gib(Some(64)), 64);
        assert_eq!(resolve_ram_limit_gib(Some(1)), 1);
        assert_eq!(resolve_ram_limit_gib(Some(0)), 0);
    }

    /// Auto-detection must never return 0, which sshash reads as "unlimited" —
    /// the single most dangerous value to arrive at by accident.
    #[test]
    fn auto_detection_never_yields_unlimited() {
        assert!(resolve_ram_limit_gib(None) > 0);
    }

    /// On any machine large enough to afford it, auto-detection must not
    /// regress below the historical default.
    #[test]
    fn auto_detection_holds_the_floor_on_a_capable_machine() {
        if usable_memory_bytes().is_some_and(|b| b >= 16 * BYTES_PER_GIB) {
            assert!(resolve_ram_limit_gib(None) >= MIN_RAM_LIMIT_GIB);
        }
    }

    /// The budget must never exceed what the process may actually use — the
    /// property the cgroup and rlimit layers exist to guarantee.
    #[test]
    fn budget_never_exceeds_usable_memory() {
        if let Some(usable) = usable_memory_bytes() {
            let budget = resolve_ram_limit_gib(None) as u64 * BYTES_PER_GIB;
            assert!(
                budget <= usable,
                "budget {budget} exceeds usable {usable}"
            );
        }
    }

    /// Whatever the platform, asking twice must give the same answer; a budget
    /// that drifts between calls would make builds irreproducible.
    #[test]
    fn detection_is_stable() {
        assert_eq!(resolve_ram_limit_gib(None), resolve_ram_limit_gib(None));
    }
}
