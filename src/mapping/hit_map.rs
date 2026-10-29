//! Per-read target accumulator map with stable slots.
//!
//! [`HitMap`] maps a target id to its hit accumulator, like the
//! `FixedMap<u32, S>` it replaces, but stores the accumulators in a `Vec` and
//! keeps only `tid → slot` in the hash table. A slot stays valid until the
//! next [`HitMap::clear`], so hit collection can resolve a contig's
//! occurrence list to slots once and reuse them for every further raw hit on
//! the same contig, instead of hashing each `(hit, occurrence)` pair again.
//!
//! Iteration visits targets in the hash table's order, which depends only on
//! the hasher, the capacity and the sequence of inserted keys — never on the
//! value type. Callers that insert the same keys in the same order therefore
//! see exactly the order the old `FixedMap<u32, S>` produced, which matters
//! because that order reaches the RAD output (see [`crate::hash`]).

use std::collections::hash_map::Entry;

use crate::hash::{FixedMap, fixed_map_with_capacity};

/// Target id → hit accumulator, with slot indices that are stable until
/// [`clear`](Self::clear).
#[derive(Debug)]
pub struct HitMap<S> {
    slots: FixedMap<u32, u32>,
    targets: Vec<S>,
}

impl<S: Default> HitMap<S> {
    /// An empty map with room for `capacity` targets.
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            slots: fixed_map_with_capacity(capacity),
            targets: Vec::with_capacity(capacity),
        }
    }

    /// The slot holding `tid`'s accumulator, inserting a default one first if
    /// `tid` is absent.
    #[inline]
    pub fn slot_or_insert(&mut self, tid: u32) -> u32 {
        match self.slots.entry(tid) {
            Entry::Occupied(e) => *e.get(),
            Entry::Vacant(e) => {
                let slot = self.targets.len() as u32;
                self.targets.push(S::default());
                e.insert(slot);
                slot
            }
        }
    }

    /// `tid`'s accumulator, inserting a default one first if absent.
    #[inline]
    pub fn entry_or_default(&mut self, tid: u32) -> &mut S {
        let slot = self.slot_or_insert(tid);
        &mut self.targets[slot as usize]
    }

    /// Insert or replace `tid`'s accumulator.
    pub fn insert(&mut self, tid: u32, target: S) {
        let slot = self.slot_or_insert(tid);
        self.targets[slot as usize] = target;
    }
}

impl<S> HitMap<S> {
    /// The accumulator in `slot`, as returned by
    /// [`slot_or_insert`](HitMap::slot_or_insert) since the last clear.
    #[inline]
    pub fn slot_mut(&mut self, slot: u32) -> &mut S {
        &mut self.targets[slot as usize]
    }

    /// `tid`'s accumulator, if present.
    #[inline]
    pub fn get(&self, tid: &u32) -> Option<&S> {
        self.slots.get(tid).map(|&s| &self.targets[s as usize])
    }

    /// `tid`'s accumulator, if present.
    #[inline]
    pub fn get_mut(&mut self, tid: &u32) -> Option<&mut S> {
        self.slots.get(tid).map(|&s| &mut self.targets[s as usize])
    }

    /// Remove every target; previously returned slots become invalid.
    #[inline]
    pub fn clear(&mut self) {
        self.slots.clear();
        self.targets.clear();
    }

    /// Number of targets.
    #[inline]
    pub fn len(&self) -> usize {
        self.targets.len()
    }

    /// Whether there are no targets.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.targets.is_empty()
    }

    /// `(tid, accumulator)` pairs in hash-table order.
    pub fn iter(&self) -> Iter<'_, S> {
        Iter {
            slots: self.slots.iter(),
            targets: &self.targets,
        }
    }

    /// Accumulators in hash-table order.
    pub fn values(&self) -> impl Iterator<Item = &S> {
        self.iter().map(|(_, t)| t)
    }
}

/// Iterator over a [`HitMap`]'s `(tid, accumulator)` pairs.
pub struct Iter<'a, S> {
    slots: std::collections::hash_map::Iter<'a, u32, u32>,
    targets: &'a [S],
}

impl<'a, S> Iterator for Iter<'a, S> {
    type Item = (&'a u32, &'a S);

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        self.slots
            .next()
            .map(|(tid, &s)| (tid, &self.targets[s as usize]))
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.slots.size_hint()
    }
}

impl<'a, S> IntoIterator for &'a HitMap<S> {
    type Item = (&'a u32, &'a S);
    type IntoIter = Iter<'a, S>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slots_are_stable_and_shared_per_tid() {
        let mut m: HitMap<u32> = HitMap::with_capacity(4);
        let a = m.slot_or_insert(7);
        let b = m.slot_or_insert(9);
        assert_ne!(a, b);
        assert_eq!(m.slot_or_insert(7), a);
        *m.slot_mut(a) += 3;
        assert_eq!(m.get(&7), Some(&3));
        assert_eq!(m.len(), 2);
        m.clear();
        assert!(m.is_empty());
        assert_eq!(m.get(&7), None);
    }

    #[test]
    fn iteration_order_matches_a_plain_fixed_map() {
        // Same keys, same insertion order, same hasher and capacity: the
        // order must be the one `FixedMap<u32, S>` gives, whatever `S` is.
        let keys: Vec<u32> = (0..500u32)
            .map(|i| i.wrapping_mul(2_654_435_761) % 100_003)
            .collect();
        let mut plain: FixedMap<u32, [u64; 7]> = fixed_map_with_capacity(256);
        let mut slotted: HitMap<[u64; 7]> = HitMap::with_capacity(256);
        for round in 0..3 {
            plain.clear();
            slotted.clear();
            for &k in keys.iter().skip(round * 37) {
                plain.entry(k).or_default()[0] += 1;
                slotted.entry_or_default(k)[0] += 1;
            }
            let a: Vec<(u32, u64)> = plain.iter().map(|(&k, v)| (k, v[0])).collect();
            let b: Vec<(u32, u64)> = slotted.iter().map(|(&k, v)| (k, v[0])).collect();
            assert_eq!(a, b);
        }
    }
}
