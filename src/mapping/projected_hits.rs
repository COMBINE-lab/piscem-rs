//! Projected hits — k-mer hits projected onto reference coordinates.
//!
//! When a k-mer is found in the SSHash dictionary, it maps to a position on a
//! unitig (contig). The contig table then tells us which reference sequences
//! that unitig occurs in, and at what positions. `ProjectedHits` holds the
//! intermediate state needed to decode each reference occurrence into a concrete
//! `RefPos` (position + orientation on a reference).
//!
//! Corresponds to the C++ `projected_hits` struct in `projected_hits.hpp`.

use crate::index::contig_table::{ContigSpan, EntryEncoding};

// ---------------------------------------------------------------------------
// RefPos
// ---------------------------------------------------------------------------

/// A decoded reference position and orientation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RefPos {
    /// Position on the reference sequence.
    pub pos: u32,
    /// Whether the k-mer is in the forward orientation on the reference.
    pub is_fw: bool,
}

// ---------------------------------------------------------------------------
// ProjectedHits
// ---------------------------------------------------------------------------

/// K-mer hit projected onto reference coordinates.
///
/// Created by `ReferenceIndex::resolve_lookup()` from an sshash-rs
/// `LookupResult`. Holds the contig-level coordinates and a `ContigSpan`
/// over the contig table entries for that unitig. Call `decode_hit()` on
/// each entry to get the final `RefPos` on the reference.
#[derive(Clone)]
pub struct ProjectedHits<'a> {
    contig_id: u32,
    /// Position of the k-mer within the contig (0-based, in bases).
    contig_pos: u32,
    /// Whether the k-mer is in forward orientation on the contig.
    contig_orientation: bool,
    /// Length of the contig in bases.
    contig_len: u32,
    /// Global k-mer position in the dictionary.
    global_pos: u64,
    /// K-mer size.
    k: u32,
    /// Reference occurrences for this contig.
    ref_range: ContigSpan<'a>,
    /// Whether this hit resulted from a fresh dictionary search ("open search")
    /// rather than an extension along a known contig.
    ///
    /// Used by poison scanning: when the right bound of a hit interval came
    /// from an open search, the poison check uses the broader `key_exists()`
    /// instead of the constrained `key_occurs_in_unitig_between()`.
    resulted_from_open_search: bool,
}

impl<'a> ProjectedHits<'a> {
    /// Create a new `ProjectedHits`.
    ///
    /// `resulted_from_open_search` defaults to `false`; set it via
    /// `set_resulted_from_open_search()` if the hit came from an open search.
    pub fn new(
        contig_id: u32,
        contig_pos: u32,
        contig_orientation: bool,
        contig_len: u32,
        global_pos: u64,
        k: u32,
        ref_range: ContigSpan<'a>,
    ) -> Self {
        Self {
            contig_id,
            contig_pos,
            contig_orientation,
            contig_len,
            global_pos,
            k,
            ref_range,
            resulted_from_open_search: false,
        }
    }

    /// Whether this contig has no reference occurrences.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.ref_range.is_empty()
    }

    /// The contig (unitig) ID.
    #[inline]
    pub fn contig_id(&self) -> u32 {
        self.contig_id
    }

    /// Position of the k-mer within the contig.
    #[inline]
    pub fn contig_pos(&self) -> u32 {
        self.contig_pos
    }

    /// Whether the k-mer is in the forward orientation on the contig.
    #[inline]
    pub fn hit_fw_on_contig(&self) -> bool {
        self.contig_orientation
    }

    /// Length of the contig in bases.
    #[inline]
    pub fn contig_len(&self) -> u32 {
        self.contig_len
    }

    /// Global k-mer position in the dictionary.
    #[inline]
    pub fn global_pos(&self) -> u64 {
        self.global_pos
    }

    /// K-mer size.
    #[inline]
    pub fn k(&self) -> u32 {
        self.k
    }

    /// Reference occurrences for this contig.
    #[inline]
    pub fn ref_range(&self) -> &ContigSpan<'a> {
        &self.ref_range
    }

    /// Number of reference occurrences.
    #[inline]
    pub fn num_hits(&self) -> usize {
        self.ref_range.len()
    }

    /// Whether this hit resulted from an open dictionary search.
    #[inline]
    pub fn resulted_from_open_search(&self) -> bool {
        self.resulted_from_open_search
    }

    /// Mark whether this hit resulted from an open dictionary search.
    #[inline]
    pub fn set_resulted_from_open_search(&mut self, v: bool) {
        self.resulted_from_open_search = v;
    }

    /// Set the global k-mer position.
    #[inline]
    pub fn set_global_pos(&mut self, pos: u64) {
        self.global_pos = pos;
    }

    /// Set the contig position.
    #[inline]
    pub fn set_contig_pos(&mut self, pos: u32) {
        self.contig_pos = pos;
    }

    /// Set the contig orientation.
    #[inline]
    pub fn set_contig_orientation(&mut self, fw: bool) {
        self.contig_orientation = fw;
    }

    /// Decode a packed contig table entry into a reference position and
    /// orientation.
    ///
    /// The 4-case logic depends on two booleans:
    /// - `contigFW`: whether the contig is stored forward on the reference
    ///   (from the packed entry's orientation bit)
    /// - `contigOri`: whether the k-mer is forward on the contig
    ///   (from the dictionary lookup)
    ///
    /// | contigFW | contigOri | rpos                                    | rfw   |
    /// |----------|-----------|-----------------------------------------|-------|
    /// | true     | true      | pos(v) + contig_pos                     | true  |
    /// | true     | false     | pos(v) + contig_pos                     | false |
    /// | false    | true      | pos(v) + contig_len - (contig_pos + k)  | false |
    /// | false    | false     | pos(v) + contig_len - (contig_pos + k)  | true  |
    ///
    /// Pattern: `rfw = contigFW == contigOri`. Position depends only on
    /// `contigFW`.
    #[inline]
    pub fn decode_hit(&self, entry: u64, encoding: &EntryEncoding) -> RefPos {
        let contig_fw = encoding.orientation(entry);
        let base_pos = encoding.pos(entry);

        let rpos = if contig_fw {
            base_pos + self.contig_pos
        } else {
            base_pos + self.contig_len - (self.contig_pos + self.k)
        };

        RefPos {
            pos: rpos,
            is_fw: contig_fw == self.contig_orientation,
        }
    }
}

impl std::fmt::Debug for ProjectedHits<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ProjectedHits")
            .field("contig_id", &self.contig_id)
            .field("contig_pos", &self.contig_pos)
            .field("contig_orientation", &self.contig_orientation)
            .field("contig_len", &self.contig_len)
            .field("global_pos", &self.global_pos)
            .field("k", &self.k)
            .field("num_hits", &self.ref_range.len())
            .field("resulted_from_open_search", &self.resulted_from_open_search)
            .finish()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::contig_table::ContigTableBuilder;

    /// Build a small test contig table and return it along with its encoding.
    ///
    /// Layout:
    /// - Contig 0: ref 0, pos 100, FW; ref 1, pos 200, RC
    /// - Contig 1: (empty — no reference occurrences)
    /// - Contig 2: ref 0, pos 500, FW
    fn make_test_table() -> (crate::index::contig_table::ContigTable, EntryEncoding) {
        let mut builder = ContigTableBuilder::new(3, 10_000, 5);
        builder.add_occurrence(0, 0, 100, true);  // contig 0: ref0 pos100 FW
        builder.add_occurrence(0, 1, 200, false); // contig 0: ref1 pos200 RC
        builder.add_occurrence(2, 0, 500, true);  // contig 2: ref0 pos500 FW
        let table = builder.build();
        let enc = table.encoding();
        (table, enc)
    }

    #[test]
    fn test_decode_hit_case1_contig_fw_kmer_fw() {
        // contigFW=true, contigOri=true → rfw=true, rpos=pos+contigPos
        let (table, enc) = make_test_table();
        let span = table.contig_entries(0);
        let entry = span.get(0); // ref 0, pos 100, FW

        let hits = ProjectedHits::new(
            0,    // contig_id
            10,   // contig_pos (k-mer at position 10 within contig)
            true, // contig_orientation (k-mer is FW on contig)
            50,   // contig_len
            999,  // global_pos
            31,   // k
            span,
        );

        let rp = hits.decode_hit(entry, &enc);
        assert_eq!(rp.pos, 100 + 10); // pos(entry) + contig_pos
        assert!(rp.is_fw);            // contigFW == contigOri → true
    }

    #[test]
    fn test_decode_hit_case2_contig_fw_kmer_rc() {
        // contigFW=true, contigOri=false → rfw=false, rpos=pos+contigPos
        let (table, enc) = make_test_table();
        let span = table.contig_entries(0);
        let entry = span.get(0); // ref 0, pos 100, FW

        let hits = ProjectedHits::new(
            0,     // contig_id
            10,    // contig_pos
            false, // contig_orientation (k-mer is RC on contig)
            50,    // contig_len
            999,   // global_pos
            31,    // k
            span,
        );

        let rp = hits.decode_hit(entry, &enc);
        assert_eq!(rp.pos, 100 + 10); // pos(entry) + contig_pos
        assert!(!rp.is_fw);           // contigFW != contigOri → false
    }

    #[test]
    fn test_decode_hit_case3_contig_rc_kmer_fw() {
        // contigFW=false, contigOri=true → rfw=false, rpos=pos+contigLen-(contigPos+k)
        let (table, enc) = make_test_table();
        let span = table.contig_entries(0);
        let entry = span.get(1); // ref 1, pos 200, RC

        let hits = ProjectedHits::new(
            0,    // contig_id
            10,   // contig_pos
            true, // contig_orientation (k-mer is FW on contig)
            50,   // contig_len
            999,  // global_pos
            31,   // k
            span,
        );

        let rp = hits.decode_hit(entry, &enc);
        // pos(entry) + contig_len - (contig_pos + k) = 200 + 50 - (10 + 31) = 209
        assert_eq!(rp.pos, 200 + 50 - (10 + 31));
        assert!(!rp.is_fw); // contigFW != contigOri → false
    }

    #[test]
    fn test_decode_hit_case4_contig_rc_kmer_rc() {
        // contigFW=false, contigOri=false → rfw=true, rpos=pos+contigLen-(contigPos+k)
        let (table, enc) = make_test_table();
        let span = table.contig_entries(0);
        let entry = span.get(1); // ref 1, pos 200, RC

        let hits = ProjectedHits::new(
            0,     // contig_id
            10,    // contig_pos
            false, // contig_orientation (k-mer is RC on contig)
            50,    // contig_len
            999,   // global_pos
            31,    // k
            span,
        );

        let rp = hits.decode_hit(entry, &enc);
        assert_eq!(rp.pos, 200 + 50 - (10 + 31));
        assert!(rp.is_fw); // contigFW == contigOri → true
    }

    #[test]
    fn test_projected_hits_accessors() {
        let (table, _enc) = make_test_table();
        let span = table.contig_entries(0);

        let mut hits = ProjectedHits::new(42, 7, true, 100, 12345, 31, span);
        assert_eq!(hits.contig_id(), 42);
        assert_eq!(hits.contig_pos(), 7);
        assert!(hits.hit_fw_on_contig());
        assert_eq!(hits.contig_len(), 100);
        assert_eq!(hits.global_pos(), 12345);
        assert_eq!(hits.k(), 31);
        assert_eq!(hits.num_hits(), 2); // contig 0 has 2 entries
        assert!(!hits.is_empty());
        assert!(!hits.resulted_from_open_search());

        hits.set_resulted_from_open_search(true);
        assert!(hits.resulted_from_open_search());
    }

    #[test]
    fn test_projected_hits_empty_contig() {
        let (table, _enc) = make_test_table();
        let span = table.contig_entries(1); // contig 1 has no entries

        let hits = ProjectedHits::new(1, 0, true, 50, 0, 31, span);
        assert!(hits.is_empty());
        assert_eq!(hits.num_hits(), 0);
    }

    #[test]
    fn test_decode_hit_boundary_pos_zero() {
        // k-mer at position 0 within the contig
        let (table, enc) = make_test_table();
        let span = table.contig_entries(0);
        let entry_fw = span.get(0); // ref 0, pos 100, FW
        let entry_rc = span.get(1); // ref 1, pos 200, RC

        let hits = ProjectedHits::new(0, 0, true, 50, 0, 31, span);

        // FW entry: pos = 100 + 0 = 100
        let rp_fw = hits.decode_hit(entry_fw, &enc);
        assert_eq!(rp_fw.pos, 100);
        assert!(rp_fw.is_fw);

        // RC entry: pos = 200 + 50 - (0 + 31) = 219
        let rp_rc = hits.decode_hit(entry_rc, &enc);
        assert_eq!(rp_rc.pos, 200 + 50 - 31);
        assert!(!rp_rc.is_fw);
    }

    #[test]
    fn test_decode_hit_boundary_pos_max() {
        // k-mer at last valid position: contig_len - k
        let k: u32 = 31;
        let contig_len: u32 = 50;
        let contig_pos = contig_len - k; // 19

        let (table, enc) = make_test_table();
        let span = table.contig_entries(0);
        let entry_fw = span.get(0); // ref 0, pos 100, FW
        let entry_rc = span.get(1); // ref 1, pos 200, RC

        let hits = ProjectedHits::new(0, contig_pos, true, contig_len, 0, k, span);

        // FW entry: pos = 100 + 19 = 119
        let rp_fw = hits.decode_hit(entry_fw, &enc);
        assert_eq!(rp_fw.pos, 100 + contig_pos);
        assert!(rp_fw.is_fw);

        // RC entry: pos = 200 + 50 - (19 + 31) = 200 + 0 = 200
        let rp_rc = hits.decode_hit(entry_rc, &enc);
        assert_eq!(rp_rc.pos, 200 + contig_len - (contig_pos + k));
        assert_eq!(rp_rc.pos, 200); // exactly at the start
        assert!(!rp_rc.is_fw);
    }
}
