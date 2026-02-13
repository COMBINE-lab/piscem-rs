//! No-constraint hit tracking — the default mode.
//!
//! Port of C++ `sketch_hit_info_no_struct_constraint`. Simply counts unique
//! read positions that produce hits in each orientation, without structural
//! chaining constraints.

use crate::mapping::hits::{HitDirection, SimpleHit, SketchHitInfo};

// ---------------------------------------------------------------------------
// SketchHitInfoSimple
// ---------------------------------------------------------------------------

/// Hit accumulator without structural constraints (the default mode).
///
/// Tracks the number of unique read positions producing hits in each
/// orientation (FW / RC), along with approximate mapping positions.
///
/// Corresponds to C++ `sketch_hit_info_no_struct_constraint`.
#[derive(Debug)]
pub struct SketchHitInfoSimple {
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
}

impl Default for SketchHitInfoSimple {
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
        }
    }
}

impl SketchHitInfo for SketchHitInfoSimple {
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
        false
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
        false
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
    fn test_add_fw_basic() {
        let mut info = SketchHitInfoSimple::default();
        assert!(info.add_fw(100, 0, 100, 31, 100, 1.0));
        assert_eq!(info.fw_hits, 1);
        assert_eq!(info.approx_pos_fw, 100); // ref_pos - read_pos = 100 - 0
    }

    #[test]
    fn test_add_fw_dedup_same_pos() {
        let mut info = SketchHitInfoSimple::default();
        assert!(info.add_fw(100, 5, 100, 31, 100, 1.0));
        // Same read_pos again — should not count
        assert!(!info.add_fw(200, 5, 100, 31, 100, 1.0));
        assert_eq!(info.fw_hits, 1);
    }

    #[test]
    fn test_add_fw_multiple() {
        let mut info = SketchHitInfoSimple::default();
        assert!(info.add_fw(100, 0, 100, 31, 100, 1.0));
        assert!(info.add_fw(105, 5, 100, 31, 100, 1.0));
        assert!(info.add_fw(110, 10, 100, 31, 100, 1.0));
        assert_eq!(info.fw_hits, 3);
        assert_eq!(info.fw_score, 3.0);
    }

    #[test]
    fn test_add_rc_basic() {
        let mut info = SketchHitInfoSimple::default();
        assert!(info.add_rc(200, 0, 100, 31, 100, 1.0));
        assert_eq!(info.rc_hits, 1);
        assert_eq!(info.first_read_pos_rc, 0);
    }

    #[test]
    fn test_add_rc_updates_tracking() {
        let mut info = SketchHitInfoSimple::default();
        assert!(info.add_rc(200, 0, 100, 31, 100, 1.0));
        assert!(info.add_rc(210, 5, 100, 31, 100, 1.0));
        assert_eq!(info.rc_hits, 2);
        // first_read_pos_rc should stay as 0 (only set on first hit)
        assert_eq!(info.first_read_pos_rc, 0);
        // rightmost_bound_rc should be updated to previous last_ref_pos_rc
        assert_eq!(info.rightmost_bound_rc, 200);
    }

    #[test]
    fn test_best_direction_fw() {
        let mut info = SketchHitInfoSimple::default();
        info.add_fw(100, 0, 100, 31, 100, 1.0);
        info.add_fw(105, 5, 100, 31, 100, 1.0);
        info.add_rc(200, 0, 100, 31, 100, 1.0);
        assert_eq!(info.best_hit_direction(), HitDirection::Fw);
    }

    #[test]
    fn test_best_direction_rc() {
        let mut info = SketchHitInfoSimple::default();
        info.add_fw(100, 0, 100, 31, 100, 1.0);
        info.add_rc(200, 0, 100, 31, 100, 1.0);
        info.add_rc(210, 5, 100, 31, 100, 1.0);
        assert_eq!(info.best_hit_direction(), HitDirection::Rc);
    }

    #[test]
    fn test_best_direction_tie() {
        let mut info = SketchHitInfoSimple::default();
        info.add_fw(100, 0, 100, 31, 100, 1.0);
        info.add_rc(200, 0, 100, 31, 100, 1.0);
        assert_eq!(info.best_hit_direction(), HitDirection::Both);
    }

    #[test]
    fn test_get_fw_hit_fields() {
        let mut info = SketchHitInfoSimple::default();
        info.add_fw(100, 0, 100, 31, 100, 1.0);
        info.add_fw(110, 10, 100, 31, 100, 2.0);
        let hit = info.get_fw_hit();
        assert!(hit.is_fw);
        assert_eq!(hit.pos, 100); // approx_pos_fw set on first hit
        assert_eq!(hit.num_hits, 2);
        assert_eq!(hit.score, 3.0);
    }

    #[test]
    fn test_get_rc_hit_fields() {
        let mut info = SketchHitInfoSimple::default();
        info.add_rc(200, 0, 100, 31, 100, 1.5);
        info.add_rc(210, 5, 100, 31, 100, 1.5);
        let hit = info.get_rc_hit();
        assert!(!hit.is_fw);
        assert_eq!(hit.num_hits, 2);
        assert_eq!(hit.score, 3.0);
    }

    #[test]
    fn test_max_hits_for_target() {
        let mut info = SketchHitInfoSimple::default();
        info.add_fw(100, 0, 100, 31, 100, 1.0);
        info.add_fw(105, 5, 100, 31, 100, 1.0);
        info.add_rc(200, 0, 100, 31, 100, 1.0);
        assert_eq!(info.max_hits_for_target(), 2);
    }

    #[test]
    fn test_inc_hits() {
        let mut info = SketchHitInfoSimple::default();
        info.inc_fw_hits();
        info.inc_fw_hits();
        info.inc_rc_hits();
        assert_eq!(info.fw_hits, 2);
        assert_eq!(info.rc_hits, 1);
        assert_eq!(info.max_hits_for_target(), 2);
    }
}
