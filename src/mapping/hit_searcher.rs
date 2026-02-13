//! HitSearcher — k-mer hit collection engine.
//!
//! Given a read sequence, scans k-mers, looks them up in the dictionary, and
//! produces a list of `(read_position, ProjectedHits)` pairs. The algorithm
//! matches the C++ `hit_searcher` implementation (same output).
//!
//! Two modes:
//! - **PERMISSIVE**: skips along contigs using SPSS verification, with binary
//!   search midpoint recovery.
//! - **STRICT**: queries every valid k-mer position, with contig extension via
//!   SPSS comparison.
//!
//! C++ reference: `piscem-cpp/include/hit_searcher.hpp`,
//! `piscem-cpp/src/hit_searcher.cpp`.

use sshash_lib::{Kmer, KmerBits};

use crate::index::reference_index::ReferenceIndex;
use crate::mapping::projected_hits::ProjectedHits;
use crate::mapping::streaming_query::PiscemStreamingQuery;

// ---------------------------------------------------------------------------
// SkippingStrategy
// ---------------------------------------------------------------------------

/// Strategy for how the hit searcher advances along a read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SkippingStrategy {
    /// Query every k-mer position; walk contigs via SPSS comparison.
    Strict,
    /// Skip along contigs using SPSS verification with binary search fallback.
    Permissive,
}

// ---------------------------------------------------------------------------
// KmerMatchType
// ---------------------------------------------------------------------------

/// Result of comparing a read k-mer against a reference k-mer from the SPSS.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(clippy::enum_variant_names)]
pub(crate) enum KmerMatchType {
    /// No match: k-mers are distinct.
    NoMatch,
    /// Identity match: ref k-mer equals the forward read k-mer.
    IdentityMatch,
    /// Twin match: ref k-mer equals the reverse complement of the read k-mer.
    TwinMatch,
}

/// Compare a read k-mer against a reference k-mer.
///
/// - `IdentityMatch`: `ref_kmer == read_kmer`
/// - `TwinMatch`: `ref_kmer == read_kmer.reverse_complement()`
/// - `NoMatch`: neither
#[inline]
fn kmer_match<const K: usize>(read_kmer: &Kmer<K>, ref_kmer: &Kmer<K>) -> KmerMatchType
where
    Kmer<K>: KmerBits,
{
    if ref_kmer == read_kmer {
        KmerMatchType::IdentityMatch
    } else if *ref_kmer == read_kmer.reverse_complement() {
        KmerMatchType::TwinMatch
    } else {
        KmerMatchType::NoMatch
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Check if a byte is a valid DNA base (A, C, G, T, case-insensitive).
#[inline]
fn is_valid_base(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
}

/// Check if a k-mer byte slice is a homopolymer (all bases identical).
#[inline]
pub(crate) fn is_homopolymer(kmer: &[u8]) -> bool {
    if kmer.is_empty() {
        return false;
    }
    let first = kmer[0];
    kmer.iter().all(|&b| b == first)
}

/// Parse a read k-mer from a byte slice. Returns `None` if the slice contains
/// invalid characters (should not happen if ReadKmerIter is used correctly).
#[inline]
fn parse_read_kmer<const K: usize>(bytes: &[u8]) -> Option<Kmer<K>>
where
    Kmer<K>: KmerBits,
{
    // Safety: ReadKmerIter guarantees only valid DNA bases in the window.
    let s = std::str::from_utf8(bytes).ok()?;
    Kmer::<K>::from_str(s).ok()
}

// ---------------------------------------------------------------------------
// ReadKmerIter
// ---------------------------------------------------------------------------

/// Iterator over valid k-mer positions in a read, handling N characters.
///
/// `Clone` for save/restore during PERMISSIVE skip-and-fallback (matching the
/// C++ `kit_tmp`/`kit_swap` pattern).
#[derive(Clone)]
pub(crate) struct ReadKmerIter<'a> {
    read: &'a [u8],
    k: usize,
    /// Current start position of the k-mer window. When exhausted,
    /// `pos + k > read.len()`.
    pos: i32,
    /// Index of the last non-ACGT byte seen so far (init: -1).
    last_invalid: i32,
}

impl<'a> ReadKmerIter<'a> {
    /// Create a new iterator. Positions at the first valid k-mer window.
    pub fn new(read: &'a [u8], k: usize) -> Self {
        let mut iter = Self {
            read,
            k,
            pos: 0,
            last_invalid: -1,
        };
        // Scan the first k-1 bases for invalid characters.
        let scan_end = k.saturating_sub(1).min(read.len());
        for (i, &base) in read.iter().enumerate().take(scan_end) {
            if !is_valid_base(base) {
                iter.last_invalid = i as i32;
            }
        }
        // Find the first valid window (scan rightmost base and beyond).
        iter.find_next_valid();
        iter
    }

    /// Whether the iterator is past the last valid k-mer position.
    #[inline]
    pub fn is_exhausted(&self) -> bool {
        self.pos as usize + self.k > self.read.len()
    }

    /// Current k-mer start position on the read.
    #[inline]
    pub fn pos(&self) -> i32 {
        self.pos
    }

    /// Byte slice of the current k-mer window. Only valid when not exhausted.
    #[inline]
    pub fn kmer_bytes(&self) -> &'a [u8] {
        let start = self.pos as usize;
        &self.read[start..start + self.k]
    }

    /// Advance by one base. Returns the actual distance moved (may be > 1
    /// if N characters are skipped).
    pub fn advance(&mut self) -> i32 {
        let old_pos = self.pos;
        self.pos += 1;
        // Check the new rightmost base entering the window.
        let new_right = self.pos as usize + self.k - 1;
        if new_right < self.read.len() && !is_valid_base(self.read[new_right]) {
            self.last_invalid = new_right as i32;
        }
        // If the window contains an invalid base, find the next valid window.
        if self.last_invalid >= self.pos {
            self.pos = self.last_invalid + 1;
            self.find_next_valid();
        }
        self.pos - old_pos
    }

    /// Advance by `n` bases. Returns the actual distance moved (may be > n
    /// if N characters are skipped).
    ///
    /// Matches C++ `operator+=(n)` semantics: performs `n` individual
    /// `advance()` steps so that each step independently handles N-skipping.
    pub fn advance_by(&mut self, n: i32) -> i32 {
        let old_pos = self.pos;
        for _ in 0..n {
            self.advance();
            if self.is_exhausted() {
                break;
            }
        }
        self.pos - old_pos
    }

    /// Scan forward from the current position to find the next valid k-mer
    /// window (one with no invalid bases).
    fn find_next_valid(&mut self) {
        let read_len = self.read.len() as i32;
        let k = self.k as i32;

        loop {
            if self.pos + k > read_len {
                break; // exhausted
            }
            // Jump past any known invalid base that's inside the current window.
            if self.last_invalid >= self.pos {
                self.pos = self.last_invalid + 1;
                if self.pos + k > read_len {
                    break;
                }
            }
            // Scan from our scan frontier to the end of the current window.
            let scan_from = ((self.last_invalid + 1).max(self.pos)) as usize;
            let scan_to = (self.pos + k - 1) as usize;

            let mut found_invalid = false;
            for i in scan_from..=scan_to {
                if !is_valid_base(self.read[i]) {
                    self.last_invalid = i as i32;
                    found_invalid = true;
                    break;
                }
            }
            if !found_invalid {
                return; // Window is valid.
            }
            // Found invalid: loop will jump past it on next iteration.
        }
    }
}

// ---------------------------------------------------------------------------
// HitSearcher
// ---------------------------------------------------------------------------

/// Core k-mer hit collection engine.
///
/// Given a read sequence, scans k-mers, looks them up in the dictionary, and
/// produces a list of `(read_position, ProjectedHits)` pairs.
pub struct HitSearcher<'idx> {
    index: &'idx ReferenceIndex,
    k: usize,
    alt_skip: u32,
    left_hits: Vec<(i32, ProjectedHits<'idx>)>,
    right_hits: Vec<(i32, ProjectedHits<'idx>)>,
}

impl<'idx> HitSearcher<'idx> {
    /// Create a new hit searcher for the given index.
    pub fn new(index: &'idx ReferenceIndex) -> Self {
        Self {
            k: index.k(),
            index,
            alt_skip: 3,
            left_hits: Vec::new(),
            right_hits: Vec::new(),
        }
    }

    /// Set the alternative skip amount (used on misses before the first hit).
    pub fn set_alt_skip(&mut self, n: u32) {
        self.alt_skip = n;
    }

    /// Clear all collected hits.
    pub fn clear(&mut self) {
        self.left_hits.clear();
        self.right_hits.clear();
    }

    /// Access the left-end hits.
    pub fn left_hits(&self) -> &[(i32, ProjectedHits<'idx>)] {
        &self.left_hits
    }

    /// Access the right-end hits.
    pub fn right_hits(&self) -> &[(i32, ProjectedHits<'idx>)] {
        &self.right_hits
    }

    /// K-mer size.
    pub fn k(&self) -> usize {
        self.k
    }

    /// Reference to the underlying index.
    pub fn index(&self) -> &'idx ReferenceIndex {
        self.index
    }

    /// Collect raw k-mer hits from `read` using the given strategy.
    ///
    /// Populates either `left_hits` or `right_hits` depending on `is_left`.
    /// Returns `true` if at least one hit was found.
    pub fn get_raw_hits_sketch<const K: usize>(
        &mut self,
        read: &[u8],
        query: &mut PiscemStreamingQuery<'_, K>,
        strategy: SkippingStrategy,
        is_left: bool,
    ) -> bool
    where
        Kmer<K>: KmerBits,
    {
        // Reset streaming query state for this new read.
        query.reset();

        let raw_hits = if is_left {
            &mut self.left_hits
        } else {
            &mut self.right_hits
        };

        let index = self.index;
        let k = self.k;

        match strategy {
            SkippingStrategy::Strict => {
                let read_end_pos = read.len() as i32 - k as i32;
                let mut iter = ReadKmerIter::new(read, k);
                Self::walk_safely_until::<K>(
                    index,
                    &mut iter,
                    query,
                    read,
                    read_end_pos,
                    raw_hits,
                );
            }
            SkippingStrategy::Permissive => {
                Self::collect_permissive::<K>(index, k, read, query, raw_hits);
            }
        }

        !raw_hits.is_empty()
    }

    // -----------------------------------------------------------------------
    // PERMISSIVE mode
    // -----------------------------------------------------------------------

    /// PERMISSIVE hit collection loop.
    ///
    /// Matches C++ `get_raw_hits_sketch` lines 914-1147.
    fn collect_permissive<const K: usize>(
        index: &'idx ReferenceIndex,
        k_usize: usize,
        read: &[u8],
        query: &mut PiscemStreamingQuery<'_, K>,
        raw_hits: &mut Vec<(i32, ProjectedHits<'idx>)>,
    ) where
        Kmer<K>: KmerBits,
    {
        let k = k_usize as i32;
        let read_len = read.len() as i32;
        let mut iter = ReadKmerIter::new(read, k_usize);

        while !iter.is_exhausted() {
            let kmer_bytes = iter.kmer_bytes();
            let read_pos = iter.pos();

            // Skip homopolymers.
            if is_homopolymer(kmer_bytes) {
                iter.advance();
                continue;
            }

            // Query the index.
            let kmer_str = unsafe { std::str::from_utf8_unchecked(kmer_bytes) };
            let result = query.lookup(kmer_str);

            if let Some(phit) = index.resolve_lookup(&result) {
                if phit.is_empty() {
                    iter.advance();
                    continue;
                }

                // --- Hit found ---
                let c_start_pos = phit.global_pos() as i64 - phit.contig_pos() as i64;
                let c_end_pos = c_start_pos + phit.contig_len() as i64;
                let c_curr_pos_initial = phit.global_pos() as i64;

                // Record hit if advancing past all previous.
                if raw_hits.is_empty() || read_pos > raw_hits.last().unwrap().0 {
                    let mut proj_hit = phit.clone();
                    proj_hit.set_resulted_from_open_search(true);
                    raw_hits.push((read_pos, proj_hit));
                }

                // Compute skip direction and distance.
                let direction: i32;
                let dist_to_contig_end: i64;
                if phit.hit_fw_on_contig() {
                    direction = 1;
                    dist_to_contig_end = c_end_pos - (c_curr_pos_initial + k as i64);
                } else {
                    direction = -1;
                    dist_to_contig_end = phit.contig_pos() as i64;
                }

                let dist_to_read_end = (read_len - k) as i64 - read_pos as i64;
                let skip_dist =
                    std::cmp::min(dist_to_read_end, dist_to_contig_end) as i32;

                if skip_dist > 1 {
                    // Save backup iterator and contig position.
                    let backup_iter = iter.clone();
                    let backup_cpos = c_curr_pos_initial;
                    let mut c_curr_pos = c_curr_pos_initial;

                    // --- Neighbor check (advance by 1) ---
                    let neighbor_dist = iter.advance();
                    if iter.is_exhausted() {
                        continue;
                    }

                    if neighbor_dist < skip_dist {
                        let found_match = Self::check_direct_match::<K>(
                            index,
                            &iter,
                            read,
                            direction,
                            neighbor_dist,
                            &mut c_curr_pos,
                            raw_hits,
                            false, // don't add hit
                        );
                        if !found_match {
                            // Go to top of loop; will do regular query at
                            // neighbor position.
                            continue;
                        } else {
                            // Restore backup for the actual skip.
                            iter = backup_iter.clone();
                            c_curr_pos = backup_cpos;
                        }
                    }

                    // --- Skip to target ---
                    let actual_dist = iter.advance_by(skip_dist);
                    if iter.is_exhausted() {
                        // Fall back to safe walk.
                        iter = backup_iter.clone();
                        iter.advance();
                        Self::walk_safely_until::<K>(
                            index,
                            &mut iter,
                            query,
                            read,
                            read_pos + skip_dist,
                            raw_hits,
                        );
                        continue;
                    }

                    let alt_iter = iter.clone();
                    let next_read_pos = iter.pos();

                    // --- Direct match at target ---
                    if actual_dist == skip_dist {
                        let found_match = Self::check_direct_match::<K>(
                            index,
                            &iter,
                            read,
                            direction,
                            skip_dist,
                            &mut c_curr_pos,
                            raw_hits,
                            true, // add hit if successful
                        );
                        if found_match {
                            iter.advance();
                            continue;
                        }
                    }

                    // --- Index query at target ---
                    let alt_kmer_bytes = iter.kmer_bytes();
                    let alt_str =
                        unsafe { std::str::from_utf8_unchecked(alt_kmer_bytes) };
                    let alt_result = query.lookup(alt_str);
                    let alt_phit = index.resolve_lookup(&alt_result);
                    let alt_found = alt_phit
                        .as_ref()
                        .is_some_and(|p| !p.is_empty());

                    if let Some(ref check_phit) = alt_phit {
                        if !check_phit.is_empty() {
                            let accept = check_phit.contig_id() == phit.contig_id()
                                && check_phit.hit_fw_on_contig()
                                    == phit.hit_fw_on_contig()
                                && if direction > 0 {
                                    check_phit.contig_pos() > phit.contig_pos()
                                } else {
                                    check_phit.contig_pos() < phit.contig_pos()
                                };

                            if accept {
                                let mut accepted = check_phit.clone();
                                accepted.set_resulted_from_open_search(false);
                                raw_hits.push((iter.pos(), accepted));
                                iter.advance();
                                continue;
                            }
                        }
                    }

                    // --- Binary search midpoint ---
                    let mut mid_acceptable = false;
                    if skip_dist > 4 {
                        let half_skip = skip_dist / 2;
                        let mut mid_iter = backup_iter.clone();
                        mid_iter.advance_by(half_skip);

                        if !mid_iter.is_exhausted() {
                            let mid_bytes = mid_iter.kmer_bytes();
                            let mid_str =
                                unsafe { std::str::from_utf8_unchecked(mid_bytes) };
                            let mid_result = query.lookup(mid_str);
                            if let Some(mid_phit) =
                                index.resolve_lookup(&mid_result)
                            {
                                if !mid_phit.is_empty() {
                                    if mid_phit.contig_id() == phit.contig_id() {
                                        // Matched first contig: add mid hit,
                                        // optionally add alt hit.
                                        let mut mp = mid_phit;
                                        mp.set_resulted_from_open_search(false);
                                        raw_hits.push((mid_iter.pos(), mp));
                                        if alt_found {
                                            if let Some(ref ap) = alt_phit {
                                                let mut ap_clone = ap.clone();
                                                ap_clone
                                                    .set_resulted_from_open_search(
                                                        true,
                                                    );
                                                raw_hits
                                                    .push((alt_iter.pos(), ap_clone));
                                            }
                                        }
                                        mid_acceptable = true;
                                    } else if alt_found {
                                        if let Some(ref ap) = alt_phit {
                                            if !ap.is_empty()
                                                && mid_phit.contig_id()
                                                    == ap.contig_id()
                                            {
                                                // Matched second contig.
                                                let mut ap_clone = ap.clone();
                                                ap_clone
                                                    .set_resulted_from_open_search(
                                                        true,
                                                    );
                                                raw_hits
                                                    .push((alt_iter.pos(), ap_clone));
                                                mid_acceptable = true;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                    if mid_acceptable {
                        // Advance past the alt position and continue.
                        let mut next_iter = alt_iter;
                        next_iter.advance();
                        iter = next_iter;
                        continue;
                    } else {
                        // Walk safely from backup+2 to the target.
                        iter = backup_iter;
                        iter.advance();
                        iter.advance();
                        Self::walk_safely_until::<K>(
                            index,
                            &mut iter,
                            query,
                            read,
                            next_read_pos,
                            raw_hits,
                        );
                        continue;
                    }
                }

                // skip_dist <= 1: just advance by 1.
                iter.advance();
            } else {
                // Miss: advance by 1.
                iter.advance();
            }
        }
    }

    // -----------------------------------------------------------------------
    // walk_safely_until
    // -----------------------------------------------------------------------

    /// Walk safely one k-mer at a time, querying the index and extending along
    /// contigs via SPSS comparison.
    ///
    /// Matches C++ `walk_safely_until` (lines 548-692).
    fn walk_safely_until<const K: usize>(
        index: &'idx ReferenceIndex,
        iter: &mut ReadKmerIter<'_>,
        query: &mut PiscemStreamingQuery<'_, K>,
        _read: &[u8],
        end_read_pos: i32,
        raw_hits: &mut Vec<(i32, ProjectedHits<'idx>)>,
    ) where
        Kmer<K>: KmerBits,
    {
        let k = index.k() as i32;

        while !iter.is_exhausted() && iter.pos() <= end_read_pos {
            let kmer_bytes = iter.kmer_bytes();

            // Skip homopolymers.
            if is_homopolymer(kmer_bytes) {
                iter.advance();
                continue;
            }

            // Query the index.
            let kmer_str = unsafe { std::str::from_utf8_unchecked(kmer_bytes) };
            let result = query.lookup(kmer_str);

            if let Some(phit_info) = index.resolve_lookup(&result) {
                if phit_info.is_empty() {
                    iter.advance();
                    continue;
                }

                let read_pos = iter.pos();
                let initial_search_pos = read_pos;

                let c_start_pos =
                    phit_info.global_pos() as i64 - phit_info.contig_pos() as i64;
                let c_end_pos = c_start_pos + phit_info.contig_len() as i64;
                let mut c_curr_pos = phit_info.global_pos() as i64;

                // Add hit if advancing past previous.
                if raw_hits.is_empty() || read_pos > raw_hits.last().unwrap().0 {
                    let mut proj_hit = phit_info.clone();
                    proj_hit.set_resulted_from_open_search(true);
                    raw_hits.push((read_pos, proj_hit));
                }

                // Compute direction and distance to contig end.
                let direction: i32;
                let mut dist_to_contig_end: i64;
                if phit_info.hit_fw_on_contig() {
                    direction = 1;
                    dist_to_contig_end = c_end_pos - (c_curr_pos + k as i64);
                } else {
                    direction = -1;
                    dist_to_contig_end = phit_info.contig_pos() as i64;
                }

                // Walk along the contig, comparing SPSS k-mers to read k-mers.
                let mut matches = true;
                let mut ended_on_match = false;
                let mut last_valid_hit = raw_hits.last().unwrap().clone();

                while !iter.is_exhausted() && matches && dist_to_contig_end > 0 {
                    let inc_amt = iter.advance();
                    dist_to_contig_end -= inc_amt as i64;

                    if dist_to_contig_end >= 0 {
                        let inc_offset = direction * inc_amt;
                        c_curr_pos += inc_offset as i64;
                        debug_assert!(c_curr_pos >= 0, "cCurrPos must be >= 0");

                        // Read ref k-mer from SPSS and compare to read k-mer.
                        if iter.is_exhausted() {
                            break;
                        }
                        let ref_kmer: Kmer<K> =
                            index.dict().kmer_at_pos(c_curr_pos as usize);
                        let read_kmer_bytes = iter.kmer_bytes();
                        if let Some(read_kmer) = parse_read_kmer::<K>(read_kmer_bytes)
                        {
                            let match_type = kmer_match(&read_kmer, &ref_kmer);
                            matches = match_type != KmerMatchType::NoMatch;

                            if matches {
                                let hit_fw =
                                    match_type == KmerMatchType::IdentityMatch;
                                let phit = &mut last_valid_hit.1;
                                phit.set_resulted_from_open_search(false);
                                phit.set_contig_orientation(hit_fw);
                                phit.set_global_pos(
                                    (phit.global_pos() as i64 + inc_offset as i64)
                                        as u64,
                                );
                                phit.set_contig_pos(
                                    (phit.contig_pos() as i32 + inc_offset) as u32,
                                );
                                last_valid_hit.0 = iter.pos();
                                ended_on_match = dist_to_contig_end == 0;
                            } else {
                                break;
                            }
                        } else {
                            break;
                        }
                    } else {
                        // Overshot the end of the contig.
                        break;
                    }
                }

                // If we walked past the initial hit, record the last valid hit.
                if last_valid_hit.0 > raw_hits.last().unwrap().0 {
                    raw_hits.push(last_valid_hit);
                }

                // If we ended on a contig boundary or didn't advance at all,
                // increment to move past.
                if ended_on_match || iter.pos() == initial_search_pos {
                    iter.advance();
                }
            } else {
                // Miss: advance by 1.
                iter.advance();
            }
        }
    }

    // -----------------------------------------------------------------------
    // check_direct_match
    // -----------------------------------------------------------------------

    /// Move forward on the current unitig by the proposed amount and check if
    /// the read k-mer matches the reference k-mer from the SPSS.
    ///
    /// If `add_hit` is true and there is a match with consistent orientation,
    /// the hit is added to `raw_hits`.
    ///
    /// Matches C++ `check_direct_match<add_hit_if_successful>` (lines 699-744).
    ///
    /// **Side effect**: `c_curr_pos` is always modified (+=direction*dist),
    /// regardless of whether the match succeeds.
    #[allow(clippy::too_many_arguments)]
    fn check_direct_match<const K: usize>(
        index: &'idx ReferenceIndex,
        iter: &ReadKmerIter<'_>,
        _read: &[u8],
        direction: i32,
        dist: i32,
        c_curr_pos: &mut i64,
        raw_hits: &mut Vec<(i32, ProjectedHits<'idx>)>,
        add_hit: bool,
    ) -> bool
    where
        Kmer<K>: KmerBits,
    {
        let inc_offset = direction as i64 * dist as i64;
        *c_curr_pos += inc_offset;

        // Read reference k-mer from SPSS at the new position.
        let ref_kmer: Kmer<K> = index.dict().kmer_at_pos(*c_curr_pos as usize);

        // Get the previous hit's orientation for consistency check.
        let prev_phit = &raw_hits.last().unwrap().1;
        let prev_hit_fw = prev_phit.hit_fw_on_contig();

        // Prepare the projected hit (copy of previous, with offset applied).
        // We'll only use this if add_hit is true and the match succeeds.
        let mut direct_phit = prev_phit.clone();
        if add_hit {
            direct_phit.set_resulted_from_open_search(false);
            direct_phit.set_global_pos(
                (direct_phit.global_pos() as i64 + inc_offset) as u64,
            );
            direct_phit.set_contig_pos(
                (direct_phit.contig_pos() as i32 + inc_offset as i32) as u32,
            );
        }

        // Compare read k-mer to ref k-mer.
        let read_kmer_bytes = iter.kmer_bytes();
        let read_kmer = match parse_read_kmer::<K>(read_kmer_bytes) {
            Some(k) => k,
            None => return false,
        };
        let match_type = kmer_match(&read_kmer, &ref_kmer);
        let matches = match_type != KmerMatchType::NoMatch;

        let hit_fw = match_type == KmerMatchType::IdentityMatch;
        if matches && (hit_fw == prev_hit_fw) {
            if add_hit {
                direct_phit.set_contig_orientation(hit_fw);
                raw_hits.push((iter.pos(), direct_phit));
            }
            true
        } else {
            false
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- ReadKmerIter tests ----

    #[test]
    fn test_kmer_iter_all_valid() {
        let read = b"ACGTACGTACGT";
        let iter = ReadKmerIter::new(read, 5);
        assert!(!iter.is_exhausted());
        assert_eq!(iter.pos(), 0);

        // Count all valid positions: 0..=(12-5) = 0..=7
        let mut it = iter;
        let mut positions = Vec::new();
        while !it.is_exhausted() {
            positions.push(it.pos());
            it.advance();
        }
        assert_eq!(positions, vec![0, 1, 2, 3, 4, 5, 6, 7]);
    }

    #[test]
    fn test_kmer_iter_with_n() {
        // "ACGTNACGTACGT" — N at position 4
        let read = b"ACGTNACGTACGT";
        let iter = ReadKmerIter::new(read, 5);
        // First valid window starts at position 5 (past the N).
        assert_eq!(iter.pos(), 5);
        assert!(!iter.is_exhausted());
        assert_eq!(iter.kmer_bytes(), b"ACGTA");
    }

    #[test]
    fn test_kmer_iter_leading_n() {
        let read = b"NNACGTACGT";
        let iter = ReadKmerIter::new(read, 5);
        // N at 0,1 → first valid at 2.
        assert_eq!(iter.pos(), 2);
        assert_eq!(iter.kmer_bytes(), b"ACGTA");
    }

    #[test]
    fn test_kmer_iter_trailing_n() {
        let read = b"ACGTACGTN";
        let iter = ReadKmerIter::new(read, 5);
        assert_eq!(iter.pos(), 0);

        // Walk until exhausted.
        let mut it = iter;
        let mut positions = Vec::new();
        while !it.is_exhausted() {
            positions.push(it.pos());
            it.advance();
        }
        // Valid positions: 0,1,2,3 (window 4 would include the N at index 8).
        assert_eq!(positions, vec![0, 1, 2, 3]);
    }

    #[test]
    fn test_kmer_iter_all_n() {
        let read = b"NNNNN";
        let iter = ReadKmerIter::new(read, 5);
        assert!(iter.is_exhausted());
    }

    #[test]
    fn test_kmer_iter_short_read() {
        let read = b"ACG";
        let iter = ReadKmerIter::new(read, 5);
        assert!(iter.is_exhausted());
    }

    #[test]
    fn test_kmer_iter_exact_k() {
        let read = b"ACGTG";
        let iter = ReadKmerIter::new(read, 5);
        assert!(!iter.is_exhausted());
        assert_eq!(iter.pos(), 0);
        assert_eq!(iter.kmer_bytes(), b"ACGTG");

        let mut it = iter;
        it.advance();
        assert!(it.is_exhausted());
    }

    #[test]
    fn test_kmer_iter_advance_by() {
        let read = b"ACGTACGTACGTACGT"; // 16 bases
        let mut iter = ReadKmerIter::new(read, 5);
        assert_eq!(iter.pos(), 0);

        let dist = iter.advance_by(5);
        assert_eq!(dist, 5);
        assert_eq!(iter.pos(), 5);
        assert_eq!(iter.kmer_bytes(), b"CGTAC");
    }

    #[test]
    fn test_kmer_iter_advance_by_with_n() {
        // 18 bases, N's at positions 8 and 9.
        let read = b"ACGTACGTNNACGTACGT";
        let mut iter = ReadKmerIter::new(read, 5);
        assert_eq!(iter.pos(), 0); // Window [0,4] = "ACGTA" valid

        // advance_by(5) does 5 individual advance() calls (matching C++
        // operator+=(n) semantics):
        //   pos 0→1→2→3→10(skip N)→11
        let dist = iter.advance_by(5);
        assert_eq!(iter.pos(), 11);
        assert_eq!(dist, 11); // jumped from 0 to 11
        assert_eq!(iter.kmer_bytes(), b"CGTAC");
    }

    #[test]
    fn test_kmer_iter_clone_restore() {
        let read = b"ACGTACGTACGT";
        let mut iter = ReadKmerIter::new(read, 5);
        assert_eq!(iter.pos(), 0);

        let backup = iter.clone();
        iter.advance();
        iter.advance();
        assert_eq!(iter.pos(), 2);

        // Restore from backup.
        iter = backup;
        assert_eq!(iter.pos(), 0);
    }

    // ---- Logic tests ----

    #[test]
    fn test_is_homopolymer() {
        assert!(is_homopolymer(b"AAAA"));
        assert!(is_homopolymer(b"CCCC"));
        assert!(is_homopolymer(b"TTTTTTTT"));
        assert!(!is_homopolymer(b"ACGT"));
        assert!(!is_homopolymer(b"AACG"));
        assert!(!is_homopolymer(b""));
    }

    #[test]
    fn test_kmer_match_identity() {
        let kmer = Kmer::<5>::from_str("ACGTG").unwrap();
        let same = Kmer::<5>::from_str("ACGTG").unwrap();
        assert_eq!(kmer_match(&kmer, &same), KmerMatchType::IdentityMatch);
    }

    #[test]
    fn test_kmer_match_twin() {
        let kmer = Kmer::<5>::from_str("ACGTG").unwrap();
        let rc = kmer.reverse_complement();
        assert_eq!(kmer_match(&kmer, &rc), KmerMatchType::TwinMatch);
    }

    #[test]
    fn test_kmer_match_none() {
        let kmer1 = Kmer::<5>::from_str("ACGTG").unwrap();
        let kmer2 = Kmer::<5>::from_str("TTTTT").unwrap();
        assert_eq!(kmer_match(&kmer1, &kmer2), KmerMatchType::NoMatch);
    }

    #[test]
    fn test_skip_distance_fw() {
        // FW on contig: dist_to_contig_end = contig_len - (contig_pos + k)
        let contig_len: i64 = 100;
        let contig_pos: i64 = 10;
        let k: i64 = 31;
        let dist = contig_len - (contig_pos + k);
        assert_eq!(dist, 59);
    }

    #[test]
    fn test_skip_distance_rc() {
        // RC on contig: dist_to_contig_end = contig_pos
        let contig_pos: i64 = 10;
        let dist = contig_pos;
        assert_eq!(dist, 10);
    }

    // ---- HitSearcher struct tests ----

    // These tests use a minimal mock. Full integration tests require a real
    // dictionary (deferred to integration test suite).

    #[test]
    fn test_is_valid_base() {
        assert!(is_valid_base(b'A'));
        assert!(is_valid_base(b'C'));
        assert!(is_valid_base(b'G'));
        assert!(is_valid_base(b'T'));
        assert!(is_valid_base(b'a'));
        assert!(is_valid_base(b'c'));
        assert!(is_valid_base(b'g'));
        assert!(is_valid_base(b't'));
        assert!(!is_valid_base(b'N'));
        assert!(!is_valid_base(b'n'));
        assert!(!is_valid_base(b'X'));
        assert!(!is_valid_base(b' '));
    }

    #[test]
    fn test_parse_read_kmer() {
        let kmer: Option<Kmer<5>> = parse_read_kmer(b"ACGTG");
        assert!(kmer.is_some());

        let bad: Option<Kmer<5>> = parse_read_kmer(b"ACNGT");
        assert!(bad.is_none());
    }
}
