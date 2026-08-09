//! hashbrown vs PHast for TinyDictionary's lookup workload.
//!
//! Both sides answer the same question — "is this k-mer present, and if so
//! what is its packed value?" — because that is what `TinyDictionary::lookup`
//! needs. That framing matters: an MPHF is *not* a dictionary. `phast::
//! Function::get` returns a `usize` for any input whatsoever, so for a k-mer
//! that is not in the reference it hands back some other key's slot. To answer
//! "present?" the PHast side must store the key at that slot and compare it.
//! Anything cheaper is answering an easier question than tiny-dict is asking.
//!
//! Structures compared:
//!   hashbrown: HashMap<u64, u64, rapidhash::fast::SeedableState>  (tiny-dict's)
//!   phast:     Function<Bits8, SeedOnly, DefaultCompressedArray, BuildRapidHash>
//!              + Vec<(key, value)> indexed by the MPHF slot   (sshash's MPHF)
//!
//! The PHast slot array is interleaved (key next to value) so that verification
//! and payload share one cache line — the layout you would actually build.

use hashbrown::HashMap;
use ph::seedable_hash::BuildRapidHash;
use ph::seeds::Bits8;
use ph::{GetSize, phast};
use rapidhash::fast::SeedableState;
use std::hint::black_box;
use std::time::Instant;

type TinyBuildHasher = SeedableState<'static>;
type Mphf = phast::Function<Bits8, phast::SeedOnly, phast::DefaultCompressedArray, BuildRapidHash>;

// ---------------------------------------------------------------------------
// deterministic RNG (splitmix64) — no rand dependency, reproducible runs
// ---------------------------------------------------------------------------

struct Rng(u64);
impl Rng {
    #[inline]
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
}

// ---------------------------------------------------------------------------
// The PHast-backed dictionary
// ---------------------------------------------------------------------------

struct PhastDict {
    f: Mphf,
    /// (key, value) at the MPHF's slot for that key. Key stored for membership.
    slots: Vec<(u64, u64)>,
}

impl PhastDict {
    fn build(keys: &[u64], vals: &[u64]) -> Self {
        let params = phast::Params::new(Bits8, phast::bits_per_seed_to_100_bucket_size(8));
        let f = Mphf::with_slice_p_hash_sc(keys, &params, BuildRapidHash, phast::SeedOnly);
        let mut slots = vec![(0u64, 0u64); keys.len()];
        for (k, v) in keys.iter().zip(vals.iter()) {
            let i = f.get(k);
            assert!(i < slots.len(), "MPHF slot {i} out of range");
            slots[i] = (*k, *v);
        }
        Self { f, slots }
    }

    #[inline]
    fn get(&self, k: u64) -> Option<u64> {
        // Bounds-checked: `get` is only a bijection over the key set, so a
        // non-key may land anywhere. This branch is perfectly predicted.
        match self.slots.get(self.f.get(&k)) {
            Some(&(sk, sv)) if sk == k => Some(sv),
            _ => None,
        }
    }

    fn size_bytes(&self) -> usize {
        self.f.size_bytes() + self.slots.len() * std::mem::size_of::<(u64, u64)>()
    }
}

// ---------------------------------------------------------------------------
// harness
// ---------------------------------------------------------------------------

const REPS: usize = 7;
const NQ: usize = 2_000_000;

/// Minimum over reps, in ns/lookup. Min rather than mean: the noise here is
/// one-sided (scheduling, migration), so the floor is the honest estimate.
fn time_lookups(name: &str, queries: &[u64], mut f: impl FnMut(u64) -> Option<u64>) -> f64 {
    let mut best = f64::MAX;
    let mut hits_seen = 0usize;
    for _ in 0..REPS {
        let t0 = Instant::now();
        let mut acc = 0u64;
        let mut hits = 0usize;
        for &q in queries {
            if let Some(v) = f(black_box(q)) {
                acc ^= v;
                hits += 1;
            }
        }
        let ns = t0.elapsed().as_nanos() as f64 / queries.len() as f64;
        black_box(acc);
        hits_seen = hits;
        best = best.min(ns);
    }
    let _ = (name, hits_seen);
    best
}

fn make_queries(keys: &[u64], hit_pct: usize, rng: &mut Rng) -> Vec<u64> {
    (0..NQ)
        .map(|_| {
            if (rng.next() % 100) < hit_pct as u64 {
                keys[(rng.next() % keys.len() as u64) as usize]
            } else {
                // Random 2-bit-packed k-mer space; collision with the key set is
                // negligible and would only make the miss path look slower.
                rng.next()
            }
        })
        .collect()
}

fn main() {
    let sizes: Vec<usize> = std::env::args()
        .skip(1)
        .filter_map(|a| a.parse().ok())
        .collect();
    let sizes = if sizes.is_empty() {
        vec![100_000, 1_000_000, 4_000_000, 8_000_000]
    } else {
        sizes
    };

    println!(
        "{:>10} {:>7} {:>12} {:>12} {:>8}   {:>10} {:>10} {:>7}",
        "keys", "hit%", "hashbrown", "phast", "ratio", "hb MB", "phast MB", "bits/k"
    );
    println!("{}", "-".repeat(88));

    for &n in &sizes {
        let mut rng = Rng(0x5EED_1234_ABCD_0001);

        // Distinct keys.
        let mut keys = Vec::with_capacity(n);
        let mut seen: HashMap<u64, (), TinyBuildHasher> =
            HashMap::with_capacity_and_hasher(n, SeedableState::fixed());
        while keys.len() < n {
            let k = rng.next();
            if seen.insert(k, ()).is_none() {
                keys.push(k);
            }
        }
        drop(seen);
        let vals: Vec<u64> = (0..n).map(|_| rng.next()).collect();

        // --- build both ---
        let t0 = Instant::now();
        let mut hb: HashMap<u64, u64, TinyBuildHasher> =
            HashMap::with_capacity_and_hasher(n, SeedableState::fixed());
        for (k, v) in keys.iter().zip(vals.iter()) {
            hb.insert(*k, *v);
        }
        let hb_build = t0.elapsed();

        let t0 = Instant::now();
        let ph = PhastDict::build(&keys, &vals);
        let ph_build = t0.elapsed();

        // --- correctness: they must agree on every query ---
        let check = make_queries(&keys, 50, &mut Rng(0xC0FFEE));
        for &q in check.iter().take(200_000) {
            assert_eq!(hb.get(&q).copied(), ph.get(q), "disagreement on {q:#x}");
        }

        // --- memory ---
        let raw_buckets = ((n * 8) / 7).next_power_of_two();
        let hb_bytes = raw_buckets * (std::mem::size_of::<(u64, u64)>() + 1);
        let ph_bytes = ph.size_bytes();
        let mphf_bits_per_key = ph.f.size_bytes() as f64 * 8.0 / n as f64;

        for hit_pct in [100usize, 50, 0] {
            let queries = make_queries(&keys, hit_pct, &mut rng);
            let hb_ns = time_lookups("hb", &queries, |q| hb.get(&q).copied());
            let ph_ns = time_lookups("ph", &queries, |q| ph.get(q));
            println!(
                "{:>10} {:>7} {:>10.2}ns {:>10.2}ns {:>8.2}   {:>10.1} {:>10.1} {:>7.2}",
                n,
                hit_pct,
                hb_ns,
                ph_ns,
                ph_ns / hb_ns,
                hb_bytes as f64 / (1 << 20) as f64,
                ph_bytes as f64 / (1 << 20) as f64,
                mphf_bits_per_key,
            );
        }
        println!(
            "{:>10} build: hashbrown {:?}, phast {:?} ({:.1}x slower)",
            n,
            hb_build,
            ph_build,
            ph_build.as_secs_f64() / hb_build.as_secs_f64()
        );
    }
    println!("\nratio > 1 means PHast is slower.");
}
