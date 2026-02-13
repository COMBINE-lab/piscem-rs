//! Chained hit tracking — optional structural constraint mode.
//!
//! Port of C++ `chain_state` + `sketch_hit_info` (the structural constraint
//! variant). Off by default; enabled with `--struct-constraints`.
//!
//! Maintains chains of consecutive hits that satisfy positional constraints,
//! ensuring that mapping positions are consistent across the read.

use smallvec::SmallVec;

use crate::mapping::hits::{HitDirection, SimpleHit, SketchHitInfo, MAX_DISTORTION};

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Maximum number of chains before we give up structural constraints.
const MAX_NUM_CHAINS: usize = 8;

/// Maximum positional stretch allowed for chain extension.
const MAX_CHAIN_STRETCH: i32 = 31;

// ---------------------------------------------------------------------------
// ChainState
// ---------------------------------------------------------------------------

/// State of a single hit chain for structural constraints.
///
/// Corresponds to C++ `chain_state`.
#[derive(Debug, Clone)]
pub struct ChainState {
    pub read_start_pos: i32,
    pub prev_pos: i32,
    pub curr_pos: i32,
    pub num_hits: u8,
    pub min_distortion: u8,
}

impl Default for ChainState {
    fn default() -> Self {
        Self {
            read_start_pos: -1,
            prev_pos: -1,
            curr_pos: -1,
            num_hits: 0,
            min_distortion: MAX_DISTORTION,
        }
    }
}

// ---------------------------------------------------------------------------
// SketchHitInfoChained
// ---------------------------------------------------------------------------

/// Hit accumulator with structural constraints (optional mode).
///
/// Maintains chains of positionally consistent hits. When too many chains
/// accumulate (overflow), falls back to simple counting (same as
/// `SketchHitInfoSimple`).
///
/// Corresponds to C++ `sketch_hit_info`.
#[derive(Debug)]
pub struct SketchHitInfoChained {
    last_read_pos_fw: i32,
    last_read_pos_rc: i32,
    last_ref_pos_fw: i32,
    last_ref_pos_rc: i32,
    approx_pos_fw: i32,
    approx_pos_rc: i32,
    approx_end_pos_rc: i32,
    first_read_pos_rc: i32,
    rightmost_bound_rc: i32,
    fw_hits: u32,
    rc_hits: u32,
    fw_score: f32,
    rc_score: f32,
    ignore_struct_constraints_fw: bool,
    ignore_struct_constraints_rc: bool,
    fw_rank: i32,
    rc_rank: i32,
    fw_chains: SmallVec<[ChainState; MAX_NUM_CHAINS]>,
    rc_chains: SmallVec<[ChainState; MAX_NUM_CHAINS]>,
}

impl Default for SketchHitInfoChained {
    fn default() -> Self {
        Self {
            last_read_pos_fw: -1,
            last_read_pos_rc: -1,
            last_ref_pos_fw: -1,
            last_ref_pos_rc: i32::MAX,
            approx_pos_fw: -1,
            approx_pos_rc: -1,
            approx_end_pos_rc: -1,
            first_read_pos_rc: -1,
            rightmost_bound_rc: i32::MAX,
            fw_hits: 0,
            rc_hits: 0,
            fw_score: 0.0,
            rc_score: 0.0,
            ignore_struct_constraints_fw: false,
            ignore_struct_constraints_rc: false,
            fw_rank: -1,
            rc_rank: -1,
            fw_chains: SmallVec::new(),
            rc_chains: SmallVec::new(),
        }
    }
}

impl SketchHitInfoChained {
    /// Remove chains that don't meet the hit count threshold, then prepare
    /// surviving chains for the next rank by shifting curr_pos → prev_pos.
    fn compact_chains(chains: &mut SmallVec<[ChainState; MAX_NUM_CHAINS]>, required_hits: u32) {
        chains.retain_mut(|s| {
            if (s.num_hits as u32) < required_hits {
                return false;
            }
            s.prev_pos = s.curr_pos;
            s.curr_pos = -1;
            true
        });
    }

    /// Process the first k-mer hit (rank 0) for a chain orientation.
    fn process_rank0_hit(
        approx_map_pos: i32,
        hit_pos: i32,
        chains: &mut SmallVec<[ChainState; MAX_NUM_CHAINS]>,
        approx_pos_out: &mut i32,
        ignore_struct_constraints: &mut bool,
        num_hits: &mut u32,
    ) -> bool {
        if chains.len() == MAX_NUM_CHAINS {
            // Too many chains — disable structural constraints.
            *ignore_struct_constraints = true;
            *approx_pos_out = chains[0].read_start_pos;
            *num_hits = chains[0].num_hits as u32;
            return true;
        }
        chains.push(ChainState {
            read_start_pos: approx_map_pos,
            prev_pos: -1,
            curr_pos: hit_pos,
            num_hits: 1,
            min_distortion: MAX_DISTORTION,
        });
        *num_hits = 1;
        true
    }

    /// Process a hit at rank > 0 by extending an existing chain.
    fn process_hit(
        is_fw_hit: bool,
        read_start_pos: i32,
        next_hit_pos: i32,
        chains: &mut SmallVec<[ChainState; MAX_NUM_CHAINS]>,
        num_hits: &mut u32,
    ) -> bool {
        if chains.is_empty() {
            return false;
        }

        // Binary search for the chain matching this read_start_pos.
        let probe_idx = chains
            .binary_search_by(|c| c.read_start_pos.cmp(&read_start_pos))
            .unwrap_or_else(|i| i);

        // For FW hits, start at the match or just before it and scan backward.
        // For RC hits, start at the match or insertion point and scan forward.
        let start_idx = if is_fw_hit {
            // Start at probe_idx itself (if valid) or the last element
            probe_idx.min(chains.len() - 1)
        } else {
            probe_idx
        };

        let mut added = false;
        let mut tries = 0;
        let mut idx = start_idx;

        while idx < chains.len() && !added && tries < 2 {
            let stretch = (chains[idx].read_start_pos - read_start_pos)
                .abs()
                .min(MAX_DISTORTION as i32);

            if chains[idx].curr_pos == -1 {
                if stretch < MAX_CHAIN_STRETCH {
                    chains[idx].curr_pos = next_hit_pos;
                    chains[idx].min_distortion = stretch as u8;
                    chains[idx].num_hits += 1;
                    *num_hits = (*num_hits).max(chains[idx].num_hits as u32);
                    added = true;
                }
            } else if stretch < chains[idx].min_distortion as i32 {
                chains[idx].min_distortion = stretch as u8;
                added = true;
            }

            if is_fw_hit {
                if idx == 0 {
                    break;
                }
                idx -= 1;
            } else {
                idx += 1;
            }
            tries += 1;
        }

        added
    }
}

impl SketchHitInfo for SketchHitInfoChained {
    fn add_fw(
        &mut self,
        ref_pos: i32,
        read_pos: i32,
        _rl: i32,
        _k: i32,
        _max_stretch: i32,
        score_inc: f32,
    ) -> bool {
        let approx_map_pos = ref_pos - read_pos;

        // If structural constraints are disabled, fall through to simple counting.
        if self.ignore_struct_constraints_fw {
            if read_pos > self.last_read_pos_fw {
                if self.last_read_pos_fw == -1 {
                    self.approx_pos_fw = approx_map_pos;
                }
                self.last_ref_pos_fw = ref_pos;
                self.last_read_pos_fw = read_pos;
                self.fw_score += score_inc;
                self.fw_hits += 1;
                return true;
            }
            return false;
        }

        // New rank of k-mer
        if read_pos > self.last_read_pos_fw {
            if self.last_read_pos_fw == -1 {
                self.approx_pos_fw = approx_map_pos;
            } else {
                Self::compact_chains(&mut self.fw_chains, self.fw_hits);
            }
            self.last_read_pos_fw = read_pos;
            self.fw_rank += 1;
        }

        if self.fw_rank == 0 {
            Self::process_rank0_hit(
                approx_map_pos,
                ref_pos,
                &mut self.fw_chains,
                &mut self.approx_pos_fw,
                &mut self.ignore_struct_constraints_fw,
                &mut self.fw_hits,
            )
        } else {
            Self::process_hit(
                true,
                approx_map_pos,
                ref_pos,
                &mut self.fw_chains,
                &mut self.fw_hits,
            )
        }
    }

    fn add_rc(
        &mut self,
        ref_pos: i32,
        read_pos: i32,
        rl: i32,
        k: i32,
        _max_stretch: i32,
        score_inc: f32,
    ) -> bool {
        let approx_map_pos = ref_pos - (rl - (read_pos + k));

        // If structural constraints are disabled, fall through to simple counting.
        if self.ignore_struct_constraints_rc {
            if read_pos > self.last_read_pos_rc {
                self.approx_pos_rc = approx_map_pos;
                if self.last_read_pos_rc == -1 {
                    self.approx_end_pos_rc = ref_pos + read_pos;
                    self.first_read_pos_rc = read_pos;
                }
                self.rc_score += score_inc;
                self.rc_hits += 1;
                self.rightmost_bound_rc = self.last_ref_pos_rc;
                self.last_ref_pos_rc = ref_pos;
                self.last_read_pos_rc = read_pos;
                return true;
            }
            return false;
        }

        // New rank of k-mer
        if read_pos > self.last_read_pos_rc {
            self.approx_pos_rc = approx_map_pos;
            if self.last_read_pos_rc == -1 {
                self.approx_end_pos_rc = ref_pos + read_pos;
                self.first_read_pos_rc = read_pos;
            } else {
                Self::compact_chains(&mut self.rc_chains, self.rc_hits);
            }
            self.rc_rank += 1;
            self.rightmost_bound_rc = self.last_ref_pos_rc;
            self.last_ref_pos_rc = ref_pos;
            self.last_read_pos_rc = read_pos;
        }

        let mut added = if self.rc_rank == 0 {
            Self::process_rank0_hit(
                approx_map_pos,
                ref_pos,
                &mut self.rc_chains,
                &mut self.approx_pos_rc,
                &mut self.ignore_struct_constraints_rc,
                &mut self.rc_hits,
            )
        } else {
            Self::process_hit(
                false,
                approx_map_pos,
                ref_pos,
                &mut self.rc_chains,
                &mut self.rc_hits,
            )
        };

        if added {
            self.approx_pos_rc = approx_map_pos;
        }
        // For RC rank 0, the C++ always sets added = true via process_rank0_hit
        // and the "added" variable above handles that.
        let _ = &mut added; // suppress unused-mut if needed
        added
    }

    #[inline]
    fn inc_fw_hits(&mut self) {
        self.fw_hits += 1;
    }

    #[inline]
    fn inc_rc_hits(&mut self) {
        self.rc_hits += 1;
    }

    #[inline]
    fn max_hits_for_target(&self) -> u32 {
        self.fw_hits.max(self.rc_hits)
    }

    #[inline]
    fn best_hit_direction(&self) -> HitDirection {
        let diff = self.fw_hits as i32 - self.rc_hits as i32;
        if diff > 0 {
            HitDirection::Fw
        } else if diff < 0 {
            HitDirection::Rc
        } else {
            HitDirection::Both
        }
    }

    fn get_fw_hit(&self) -> SimpleHit {
        SimpleHit {
            is_fw: true,
            mate_is_fw: false,
            pos: self.approx_pos_fw,
            score: self.fw_score,
            num_hits: self.fw_hits,
            tid: u32::MAX,
            ..SimpleHit::default()
        }
    }

    fn get_rc_hit(&self) -> SimpleHit {
        SimpleHit {
            is_fw: false,
            mate_is_fw: false,
            pos: self.approx_pos_rc,
            score: self.rc_score,
            num_hits: self.rc_hits,
            tid: u32::MAX,
            ..SimpleHit::default()
        }
    }

    fn get_best_hit(&self) -> SimpleHit {
        let dir = self.best_hit_direction();
        if dir != HitDirection::Rc {
            self.get_fw_hit()
        } else {
            self.get_rc_hit()
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compact_chains() {
        let mut chains: SmallVec<[ChainState; MAX_NUM_CHAINS]> = SmallVec::new();
        chains.push(ChainState {
            read_start_pos: 100,
            prev_pos: -1,
            curr_pos: 50,
            num_hits: 3,
            min_distortion: 5,
        });
        chains.push(ChainState {
            read_start_pos: 200,
            prev_pos: -1,
            curr_pos: 60,
            num_hits: 1,
            min_distortion: 10,
        });
        chains.push(ChainState {
            read_start_pos: 300,
            prev_pos: -1,
            curr_pos: 70,
            num_hits: 3,
            min_distortion: 2,
        });

        SketchHitInfoChained::compact_chains(&mut chains, 3);
        assert_eq!(chains.len(), 2); // chain with num_hits=1 removed
        // curr_pos should have shifted to prev_pos
        assert_eq!(chains[0].prev_pos, 50);
        assert_eq!(chains[0].curr_pos, -1);
        assert_eq!(chains[1].prev_pos, 70);
        assert_eq!(chains[1].curr_pos, -1);
    }

    #[test]
    fn test_rank0_overflow_disables() {
        let mut info = SketchHitInfoChained::default();
        // Fill up to MAX_NUM_CHAINS hits at rank 0 for FW
        for i in 0..MAX_NUM_CHAINS {
            info.add_fw(100 + i as i32, 0, 100, 31, 100, 1.0);
        }
        assert!(!info.ignore_struct_constraints_fw);

        // One more should trigger the overflow
        info.add_fw(100 + MAX_NUM_CHAINS as i32, 0, 100, 31, 100, 1.0);
        assert!(info.ignore_struct_constraints_fw);
    }

    #[test]
    fn test_chain_extension() {
        let mut info = SketchHitInfoChained::default();
        // rank 0: first hit
        assert!(info.add_fw(100, 0, 100, 31, 100, 1.0));
        assert_eq!(info.fw_hits, 1);

        // rank 1: extend the chain
        assert!(info.add_fw(105, 5, 100, 31, 100, 1.0));
        assert_eq!(info.fw_hits, 2);
    }

    #[test]
    fn test_chained_matches_simple_when_disabled() {
        // When structural constraints are disabled (via overflow),
        // behavior should match SketchHitInfoSimple.
        let mut info = SketchHitInfoChained::default();
        info.ignore_struct_constraints_fw = true;
        info.ignore_struct_constraints_rc = true;

        assert!(info.add_fw(100, 0, 100, 31, 100, 1.0));
        assert!(info.add_fw(105, 5, 100, 31, 100, 1.0));
        assert!(!info.add_fw(110, 5, 100, 31, 100, 1.0)); // same read_pos, rejected
        assert_eq!(info.fw_hits, 2);

        assert!(info.add_rc(200, 0, 100, 31, 100, 1.0));
        assert!(info.add_rc(210, 5, 100, 31, 100, 1.0));
        assert_eq!(info.rc_hits, 2);
    }

    #[test]
    fn test_chained_best_direction() {
        let mut info = SketchHitInfoChained::default();
        info.ignore_struct_constraints_fw = true;
        info.ignore_struct_constraints_rc = true;

        info.add_fw(100, 0, 100, 31, 100, 1.0);
        info.add_fw(105, 5, 100, 31, 100, 1.0);
        info.add_rc(200, 0, 100, 31, 100, 1.0);
        assert_eq!(info.best_hit_direction(), HitDirection::Fw);
        assert_eq!(info.max_hits_for_target(), 2);
    }

    #[test]
    fn test_chained_get_hits() {
        let mut info = SketchHitInfoChained::default();
        info.ignore_struct_constraints_fw = true;
        info.ignore_struct_constraints_rc = true;

        info.add_fw(100, 0, 100, 31, 100, 1.0);
        info.add_rc(200, 0, 100, 31, 100, 2.0);

        let fw = info.get_fw_hit();
        assert!(fw.is_fw);
        assert_eq!(fw.num_hits, 1);

        let rc = info.get_rc_hit();
        assert!(!rc.is_fw);
        assert_eq!(rc.num_hits, 1);
        assert_eq!(rc.score, 2.0);

        let best = info.get_best_hit();
        // Tied → FW wins
        assert!(best.is_fw);
    }
}
