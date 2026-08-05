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
/// Handles cgroup v2 (`memory.max`) and v1 (`memory.limit_in_bytes`), resolving
/// the process's own path from `/proc/self/cgroup` and then taking the
/// **minimum over that cgroup and every ancestor**.
///
/// Walking the chain matters as much as finding the leaf. The effective ceiling
/// is the tightest limit anywhere above the process, and a scheduler routinely
/// caps an ancestor while the job's own cgroup reads `max` — so reading only
/// the leaf reports "unlimited" for a confined job, which is precisely the
/// failure this module exists to prevent. (Raised by @BenjaminDEMAILLE in
/// review of PR #4, against an earlier version that read only the leaf.)
///
/// This is the case that matters: SLURM, LSF and container runtimes all limit
/// memory here rather than through rlimits.
#[cfg(target_os = "linux")]
fn cgroup_limit_bytes() -> Option<u64> {
    use std::fs;
    use std::path::PathBuf;

    fn parse(s: &str) -> Option<u64> {
        let s = s.trim();
        if s == "max" {
            return None;
        }
        let v: u64 = s.parse().ok()?;
        // cgroup v1 spells "unlimited" as a huge sentinel rather than a word.
        if v >= i64::MAX as u64 / 2 { None } else { Some(v) }
    }

    /// Minimum limit over `dir` and every ancestor up to `root`.
    ///
    /// The effective ceiling is the tightest limit anywhere on the chain, not
    /// the one on the leaf: a scheduler routinely caps an ancestor while the
    /// job's own cgroup reads `max`. Reading only the leaf therefore reports
    /// "unlimited" for a confined job — the exact failure this module exists to
    /// prevent.
    fn min_along_chain(root: &str, mut dir: PathBuf, file: &str) -> Option<u64> {
        let mut best: Option<u64> = None;
        loop {
            if let Some(v) = fs::read_to_string(dir.join(file)).ok().and_then(|s| parse(&s)) {
                best = Some(best.map_or(v, |b: u64| b.min(v)));
            }
            if dir.as_os_str() == root || !dir.pop() {
                break;
            }
        }
        best
    }

    // The process's own cgroup path. v2 lines start "0::"; v1 lines name their
    // controllers, and we want the one carrying `memory`.
    let cgroup_file = fs::read_to_string("/proc/self/cgroup").unwrap_or_default();
    let v2_rel = cgroup_file
        .lines()
        .find_map(|l| l.strip_prefix("0::"))
        .map(|p| p.trim_start_matches('/').to_string());
    let v1_rel = cgroup_file.lines().find_map(|l| {
        let mut f = l.splitn(3, ':');
        let (_, ctrls, path) = (f.next()?, f.next()?, f.next()?);
        ctrls
            .split(',')
            .any(|c| c == "memory")
            .then(|| path.trim_start_matches('/').to_string())
    });

    let mut found: Option<u64> = None;
    let mut take = |v: Option<u64>| {
        if let Some(v) = v {
            found = Some(found.map_or(v, |b: u64| b.min(v)));
        }
    };

    // cgroup v2: unified hierarchy at the mount root.
    const V2_ROOT: &str = "/sys/fs/cgroup";
    take(min_along_chain(
        V2_ROOT,
        PathBuf::from(V2_ROOT).join(v2_rel.clone().unwrap_or_default()),
        "memory.max",
    ));

    // cgroup v1: the memory controller is mounted in its own subdirectory.
    const V1_ROOT: &str = "/sys/fs/cgroup/memory";
    take(min_along_chain(
        V1_ROOT,
        PathBuf::from(V1_ROOT).join(v1_rel.unwrap_or_default()),
        "memory.limit_in_bytes",
    ));

    if found.is_none() && v2_rel.is_some() {
        // Under a rootless container the path from /proc/self/cgroup may not
        // resolve inside the namespace, and we fall back to whatever the mount
        // root reports — possibly a larger limit than the real one. That is
        // today's behaviour rather than a new hazard, but it is silent, so say
        // so at debug level.
        debug!("no cgroup memory limit resolved; falling back to other sources");
    }
    found
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
