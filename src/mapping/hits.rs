//! Core mapping types — enums, structs, and traits used throughout the mapping
//! pipeline.
//!
//! Corresponds to the shared types in the C++ `mapping/utils.hpp`.

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Sentinel value for an invalid fragment length.
pub const INVALID_FRAG_LEN: i32 = i32::MIN;

/// Sentinel value for an invalid mate position.
pub const INVALID_MATE_POS: i32 = i32::MIN;

/// Maximum distortion allowed in chain-based structural constraints.
pub const MAX_DISTORTION: u8 = u8::MAX;

// ---------------------------------------------------------------------------
// MappingType
// ---------------------------------------------------------------------------

/// How a read (or read pair) mapped.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum MappingType {
    #[default]
    Unmapped = 0,
    SingleMapped = 1,
    MappedFirstOrphan = 2,
    MappedSecondOrphan = 3,
    MappedPair = 4,
}

// ---------------------------------------------------------------------------
// HitDirection
// ---------------------------------------------------------------------------

/// Best-scoring orientation for a target.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum HitDirection {
    Fw = 0,
    Rc = 1,
    Both = 2,
}

// ---------------------------------------------------------------------------
// FragmentEnd
// ---------------------------------------------------------------------------

/// Which end of a read pair we are processing.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum FragmentEnd {
    #[default]
    Left = 0,
    Right = 1,
}

// ---------------------------------------------------------------------------
// SimpleHit
// ---------------------------------------------------------------------------

/// A resolved mapping hit on a reference target.
///
/// Corresponds to the C++ `simple_hit` struct.
#[derive(Debug, Clone, Copy)]
pub struct SimpleHit {
    pub is_fw: bool,
    pub mate_is_fw: bool,
    pub pos: i32,
    pub score: f32,
    pub num_hits: u32,
    pub tid: u32,
    pub mate_pos: i32,
    pub fragment_length: i32,
    /// Bin ID for scATAC binned mapping. `u64::MAX` = invalid/unused.
    pub bin_id: u64,
}

impl Default for SimpleHit {
    fn default() -> Self {
        Self {
            is_fw: false,
            mate_is_fw: false,
            pos: -1,
            score: 0.0,
            num_hits: 0,
            tid: u32::MAX,
            mate_pos: INVALID_MATE_POS,
            fragment_length: INVALID_FRAG_LEN,
            bin_id: u64::MAX,
        }
    }
}

impl SimpleHit {
    /// Whether this hit has a valid mate position.
    #[inline]
    pub fn has_mate(&self) -> bool {
        self.mate_pos != INVALID_MATE_POS
    }

    /// Fragment length, or 0 if invalid.
    #[inline]
    pub fn frag_len(&self) -> i32 {
        if self.fragment_length != INVALID_FRAG_LEN {
            self.fragment_length
        } else {
            0
        }
    }
}

// ---------------------------------------------------------------------------
// SketchHitInfo trait
// ---------------------------------------------------------------------------

/// Generic interface for accumulating k-mer hits on a target.
///
/// Implemented by both `SketchHitInfoSimple` (no structural constraints,
/// the default) and `SketchHitInfoChained` (structural constraint variant).
///
/// The trait is used as a generic parameter in `map_read<K, S>()` for static
/// dispatch via monomorphization.
pub trait SketchHitInfo: Default {
    /// Record a forward-strand hit.
    fn add_fw(
        &mut self,
        ref_pos: i32,
        read_pos: i32,
        rl: i32,
        k: i32,
        max_stretch: i32,
        score_inc: f32,
    ) -> bool;

    /// Record a reverse-complement hit.
    fn add_rc(
        &mut self,
        ref_pos: i32,
        read_pos: i32,
        rl: i32,
        k: i32,
        max_stretch: i32,
        score_inc: f32,
    ) -> bool;

    /// Increment forward hit count (used by EC filtering).
    fn inc_fw_hits(&mut self);

    /// Increment RC hit count (used by EC filtering).
    fn inc_rc_hits(&mut self);

    /// Maximum hit count across both orientations.
    fn max_hits_for_target(&self) -> u32;

    /// Best scoring orientation.
    fn best_hit_direction(&self) -> HitDirection;

    /// Retrieve the forward-strand hit summary.
    fn get_fw_hit(&self) -> SimpleHit;

    /// Retrieve the RC-strand hit summary.
    fn get_rc_hit(&self) -> SimpleHit;

    /// Retrieve the best-scoring hit summary.
    fn get_best_hit(&self) -> SimpleHit;
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_hit_defaults() {
        let hit = SimpleHit::default();
        assert!(!hit.is_fw);
        assert!(!hit.mate_is_fw);
        assert_eq!(hit.pos, -1);
        assert_eq!(hit.score, 0.0);
        assert_eq!(hit.num_hits, 0);
        assert_eq!(hit.tid, u32::MAX);
        assert_eq!(hit.mate_pos, INVALID_MATE_POS);
        assert_eq!(hit.fragment_length, INVALID_FRAG_LEN);
        assert!(!hit.has_mate());
        assert_eq!(hit.frag_len(), 0);
    }

    #[test]
    fn test_mapping_type_repr() {
        assert_eq!(MappingType::Unmapped as u8, 0);
        assert_eq!(MappingType::SingleMapped as u8, 1);
        assert_eq!(MappingType::MappedFirstOrphan as u8, 2);
        assert_eq!(MappingType::MappedSecondOrphan as u8, 3);
        assert_eq!(MappingType::MappedPair as u8, 4);
        assert_eq!(MappingType::default(), MappingType::Unmapped);
    }

    #[test]
    fn test_hit_direction_from_counts() {
        // FW > RC → FW
        let fw: i32 = 5;
        let rc: i32 = 3;
        let diff = fw - rc;
        let dir = if diff > 0 {
            HitDirection::Fw
        } else if diff < 0 {
            HitDirection::Rc
        } else {
            HitDirection::Both
        };
        assert_eq!(dir, HitDirection::Fw);

        // RC > FW → RC
        let fw: i32 = 2;
        let rc: i32 = 7;
        let diff = fw - rc;
        let dir = if diff > 0 {
            HitDirection::Fw
        } else if diff < 0 {
            HitDirection::Rc
        } else {
            HitDirection::Both
        };
        assert_eq!(dir, HitDirection::Rc);

        // FW == RC → Both
        let fw: i32 = 4;
        let rc: i32 = 4;
        let diff = fw - rc;
        let dir = if diff > 0 {
            HitDirection::Fw
        } else if diff < 0 {
            HitDirection::Rc
        } else {
            HitDirection::Both
        };
        assert_eq!(dir, HitDirection::Both);
    }

    #[test]
    fn test_simple_hit_with_mate() {
        let hit = SimpleHit {
            is_fw: true,
            mate_is_fw: false,
            pos: 100,
            score: 5.0,
            num_hits: 5,
            tid: 42,
            mate_pos: 200,
            fragment_length: 150,
            ..SimpleHit::default()
        };
        assert!(hit.has_mate());
        assert_eq!(hit.frag_len(), 150);
    }
}
