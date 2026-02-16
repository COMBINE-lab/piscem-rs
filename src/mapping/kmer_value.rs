//! Strong-typed canonical k-mer value.
//!
//! `CanonicalKmer` wraps a packed 2-bit k-mer representation, replacing bare
//! `u64` values in the poison table, canonical k-mer iteration, and unitig-end
//! cache. This provides type safety (a `u64` could be anything) and a single
//! point of change for upgrading to wider k-mer support.
//!
//! # Current state
//!
//! `CanonicalKmer` wraps `u64`, supporting k <= 31 (2 bits x 31 = 62 bits).
//! This matches the C++ piscem limitation (`static_assert(K <= 32)` in
//! Kmer.hpp).
//!
//! # Upgrade path for k <= 63
//!
//! sshash-rs `Kmer<K>` already supports k up to 63 via `u128` backing (the
//! `KmerBits` trait selects `u64` for K <= 31, `u128` for K > 31). To support
//! k > 31 in piscem-rs:
//!
//! 1. **Change backing**: Use `u128` always (simplest), an enum
//!    `Small(u64)/Large(u128)` (optimal size), or generic `CanonicalKmer<B>`
//!    (compile-time dispatch via `KmerBits`).
//!
//! 2. **PoisonTable**: `HashMap<CanonicalKmer, u64>` works with any backing
//!    that implements `Hash + Eq`.
//!
//! 3. **CanonicalKmerIter**: Rolling k-mer state (`fw_kmer`, `rc_kmer`,
//!    `mask`) must use the same backing type.
//!
//! 4. **UnitigEndCache**: `DashMap<CanonicalKmer, CachedLookup>` works
//!    directly.
//!
//! 5. **Serialization**: A new format version (e.g., `PPOIS02\0`) would write
//!    the wider type, with backward-compatible loading.
//!
//! 6. **Construct from `Kmer<K>`**: Add `CanonicalKmer::from_kmer()` using
//!    `KmerBits::to_u64()` (k <= 31) or `to_u128()` (k > 31).

/// A canonical k-mer packed into a `u64`.
///
/// The minimum of the forward and reverse-complement 2-bit encodings.
/// Supports k <= 31 (62 bits). See module docs for upgrade path.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CanonicalKmer(pub(crate) u64);

impl CanonicalKmer {
    /// Create from a raw packed value.
    #[inline]
    pub fn new(raw: u64) -> Self {
        Self(raw)
    }

    /// The raw packed `u64` value.
    #[inline]
    pub fn as_u64(self) -> u64 {
        self.0
    }
}

impl std::fmt::Debug for CanonicalKmer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "CanonicalKmer(0x{:016x})", self.0)
    }
}

impl std::fmt::Display for CanonicalKmer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "0x{:016x}", self.0)
    }
}
