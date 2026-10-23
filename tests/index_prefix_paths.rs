//! Index file naming, end to end through the `piscem-rs` binary.
//!
//! Every index file is the full prefix plus `.<ext>`. These tests build small
//! synthetic indexes (random references, each its own unitig, so no
//! cuttlefish is needed) and check that dotted prefixes neither collide nor
//! lose their decoys, that a dotted index in the pre-0.11.1 layout (core files
//! at the replaced-suffix paths) still loads and maps identically, with the
//! decoys 0.11.0 applied to it, and that a prefix naming a directory is
//! rejected before anything is written.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use piscem_rs::index::paths::{
    BUILD_ARTIFACT_EXTENSIONS, DICT_EXT, POISON_ARTIFACT_EXTENSIONS, POISON_EXT,
    index_artifact_path,
};
use tempfile::TempDir;

const READ_LEN: usize = 100;

/// Deterministic xorshift64* stream, so every run builds the same data.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_f491_4f6c_dd1d)
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
    fn seq(&mut self, len: usize) -> Vec<u8> {
        (0..len).map(|_| b"ACGT"[self.below(4)]).collect()
    }
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        })
        .collect()
}

struct Refs {
    names: Vec<String>,
    seqs: Vec<Vec<u8>>,
}

fn random_refs(seed: u64, tag: &str, n: usize, len: usize) -> Refs {
    let mut rng = Rng(seed);
    Refs {
        names: (0..n).map(|i| format!("{tag}{i}")).collect(),
        seqs: (0..n).map(|_| rng.seq(len)).collect(),
    }
}

/// Cuttlefish-style output in which every reference is one unitig.
fn write_cfish(prefix: &Path, refs: &Refs) {
    let mut seg = String::new();
    let mut seq = String::new();
    for (i, (name, s)) in refs.names.iter().zip(&refs.seqs).enumerate() {
        seg += &format!("{i}\t{}\n", String::from_utf8_lossy(s));
        seq += &format!("Reference:{i}_Sequence:{name}\t{i}+\n");
    }
    std::fs::write(index_artifact_path(prefix, "cf_seg"), seg).unwrap();
    std::fs::write(index_artifact_path(prefix, "cf_seq"), seq).unwrap();
    std::fs::write(index_artifact_path(prefix, "json"), "{}").unwrap();
}

/// The signature file the piscem wrapper writes next to an index.
fn write_sigs(prefix: &Path, num_records: usize) {
    let json = format!(
        r#"{{"sha256_names":"aa","sha256_seqs":"bb","sha512_names":"cc","sha512_seqs":"dd","num_records":{num_records},"num_trimmed_polya":0}}"#
    );
    std::fs::write(index_artifact_path(prefix, "sigs.json"), json).unwrap();
}

/// `n` reads drawn from `seqs`, half of them reverse-complemented.
fn write_reads(path: &Path, seqs: &[Vec<u8>], n: usize, seed: u64) {
    let mut rng = Rng(seed);
    let mut fq = String::new();
    for i in 0..n {
        let s = &seqs[rng.below(seqs.len())];
        let start = rng.below(s.len() - READ_LEN + 1);
        let mut read = s[start..start + READ_LEN].to_vec();
        if rng.below(2) == 1 {
            read = revcomp(&read);
        }
        fq += &format!(
            "@r{i}\n{}\n+\n{}\n",
            String::from_utf8_lossy(&read),
            "I".repeat(READ_LEN)
        );
    }
    std::fs::write(path, fq).unwrap();
}

/// Copies of `seqs` with a substitution every 50 bases: reads from them keep
/// enough reference k-mers to map, and the k-mers spanning each substitution
/// are decoy-only, so `build-poison` finds them.
fn decoys_of(seqs: &[Vec<u8>]) -> Vec<Vec<u8>> {
    seqs.iter()
        .map(|s| {
            let mut d = s.clone();
            for i in (25..d.len()).step_by(50) {
                d[i] = match d[i] {
                    b'A' => b'C',
                    b'C' => b'G',
                    b'G' => b'T',
                    _ => b'A',
                };
            }
            d
        })
        .collect()
}

fn write_fasta(path: &Path, seqs: &[Vec<u8>]) {
    let mut fa = String::new();
    for (i, s) in seqs.iter().enumerate() {
        fa += &format!(">decoy{i}\n{}\n", String::from_utf8_lossy(s));
    }
    std::fs::write(path, fa).unwrap();
}

fn piscem(args: &[&std::ffi::OsStr]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_piscem-rs"))
        .args(args)
        .output()
        .expect("failed to run piscem-rs")
}

fn log_of(out: &Output) -> String {
    format!(
        "{}{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    )
}

fn ok(out: Output, what: &str) -> String {
    let log = log_of(&out);
    assert!(out.status.success(), "{what} failed:\n{log}");
    log
}

struct Work {
    dir: TempDir,
}

impl Work {
    fn new() -> Self {
        Self {
            dir: tempfile::tempdir().unwrap(),
        }
    }

    fn path(&self, rel: &str) -> PathBuf {
        self.dir.path().join(rel)
    }

    /// Build an index at `out` from `refs` (`-t 1`, EC table and, by the
    /// `auto` rule for small indexes, the Tiny artifacts).
    fn build(&self, out: &Path, refs: &Refs) {
        let cfish = self.path(&format!(
            "cfish_{}",
            out.file_name().unwrap().to_string_lossy()
        ));
        write_cfish(&cfish, refs);
        let tmp = self.path("sshash_tmp");
        let log = ok(
            piscem(&[
                "build".as_ref(),
                "-i".as_ref(),
                cfish.as_os_str(),
                "-o".as_ref(),
                out.as_os_str(),
                "-k".as_ref(),
                "31".as_ref(),
                "-m".as_ref(),
                "19".as_ref(),
                "-t".as_ref(),
                "1".as_ref(),
                "--build-ec-table".as_ref(),
                "--tmp-dir".as_ref(),
                tmp.as_os_str(),
            ]),
            "build",
        );
        assert!(!log.contains("WARN"), "unexpected warning:\n{log}");
    }

    fn build_poison(&self, index: &Path, decoys: &Path) {
        ok(
            piscem(&[
                "build-poison".as_ref(),
                "-i".as_ref(),
                index.as_os_str(),
                "-d".as_ref(),
                decoys.as_os_str(),
                "-t".as_ref(),
                "1".as_ref(),
            ]),
            "build-poison",
        );
    }

    /// Map single-end `reads` with `-t 1`; returns the RAD bytes, the
    /// map_info JSON and the log.
    fn map(
        &self,
        index: &Path,
        reads: &Path,
        out: &str,
        extra: &[&str],
    ) -> (Vec<u8>, serde_json::Value, String) {
        let stem = self.path(out);
        let mut args: Vec<&std::ffi::OsStr> = vec![
            "map-bulk".as_ref(),
            "-i".as_ref(),
            index.as_os_str(),
            "-r".as_ref(),
            reads.as_os_str(),
            "-o".as_ref(),
            stem.as_os_str(),
            "-t".as_ref(),
            "1".as_ref(),
        ];
        args.extend(extra.iter().map(std::ffi::OsStr::new));
        let log = ok(piscem(&args), "map-bulk");
        let rad = std::fs::read(index_artifact_path(&stem, "rad")).unwrap();
        let info: serde_json::Value = serde_json::from_slice(
            &std::fs::read(index_artifact_path(&stem, "map_info.json")).unwrap(),
        )
        .unwrap();
        (rad, info, log)
    }
}

fn num(info: &serde_json::Value, key: &str) -> u64 {
    info[key]
        .as_u64()
        .unwrap_or_else(|| panic!("map_info lacks {key}: {info}"))
}

fn dir_listing(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = std::fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    names
}

/// The files `build` wrote for `prefix` (EC table and Tiny artifacts on).
fn built_files(prefix: &str) -> Vec<String> {
    let mut v: Vec<String> = BUILD_ARTIFACT_EXTENSIONS
        .iter()
        .map(|e| format!("{prefix}.{e}"))
        .collect();
    v.sort();
    v
}

#[test]
fn dotted_prefixes_in_one_directory_do_not_collide() {
    let w = Work::new();
    let refs_a = random_refs(1, "a", 20, 1500);
    let refs_b = random_refs(2, "b", 30, 1500);
    std::fs::create_dir(w.path("idx")).unwrap();
    let v49 = w.path("idx/gencode.v49");
    let v44 = w.path("idx/gencode.v44");
    w.build(&v49, &refs_a);
    w.build(&v44, &refs_b);

    let mut expected = built_files("gencode.v49");
    expected.extend(built_files("gencode.v44"));
    expected.sort();
    assert_eq!(dir_listing(&w.path("idx")), expected);

    let reads_a = w.path("a.fq");
    let reads_b = w.path("b.fq");
    write_reads(&reads_a, &refs_a.seqs, 400, 11);
    write_reads(&reads_b, &refs_b.seqs, 400, 12);

    let (rad_a, info, _) = w.map(&v49, &reads_a, "a_v49", &[]);
    assert_eq!(num(&info, "num_refs"), 20);
    assert_eq!(num(&info, "num_mapped"), 400);
    let (_, info, _) = w.map(&v44, &reads_b, "b_v44", &[]);
    assert_eq!(num(&info, "num_refs"), 30);
    assert_eq!(num(&info, "num_mapped"), 400);
    let (_, info, _) = w.map(&v44, &reads_a, "a_v44", &[]);
    assert_eq!(
        num(&info, "num_mapped"),
        0,
        "v44 must not hold v49's sequences"
    );

    // The same index at an undotted prefix: identical files, identical
    // mapping output.
    let plain = w.path("idx_plain/gencode_v49");
    w.build(&plain, &refs_a);
    for ext in BUILD_ARTIFACT_EXTENSIONS {
        assert_eq!(
            std::fs::read(index_artifact_path(&v49, ext)).unwrap(),
            std::fs::read(index_artifact_path(&plain, ext)).unwrap(),
            ".{ext} differs between the dotted and undotted builds"
        );
    }
    let (rad_plain, _, _) = w.map(&plain, &reads_a, "a_plain", &[]);
    assert_eq!(rad_a, rad_plain);
    for dict in ["sshash", "tiny"] {
        let (rad, _, _) = w.map(&v49, &reads_a, &format!("a_v49_{dict}"), &["--dict", dict]);
        let (rad_p, _, _) = w.map(&plain, &reads_a, &format!("a_p_{dict}"), &["--dict", dict]);
        assert_eq!(rad, rad_p, "--dict {dict}");
    }
}

#[test]
fn dotted_prefix_applies_its_poison_table() {
    let w = Work::new();
    let refs = random_refs(3, "t", 20, 1500);
    let decoys = decoys_of(&refs.seqs[..8]);
    let decoy_fa = w.path("decoys.fa");
    write_fasta(&decoy_fa, &decoys);
    let reads = w.path("decoy_reads.fq");
    write_reads(&reads, &decoys, 300, 13);

    let dotted = w.path("idx/ref.v1");
    let plain = w.path("idx/refv1");
    for prefix in [&dotted, &plain] {
        w.build(prefix, &refs);
        w.build_poison(prefix, &decoy_fa);
        for ext in POISON_ARTIFACT_EXTENSIONS {
            assert!(index_artifact_path(prefix, ext).exists(), "{ext}");
        }
    }

    let (rad_plain, info_plain, _) = w.map(&plain, &reads, "plain", &[]);
    let (rad_dotted, info_dotted, _) = w.map(&dotted, &reads, "dotted", &[]);
    let (_, info_off, _) = w.map(&dotted, &reads, "dotted_off", &["--no-poison"]);
    assert!(num(&info_plain, "num_poisoned") > 0, "{info_plain}");
    assert_eq!(
        num(&info_dotted, "num_poisoned"),
        num(&info_plain, "num_poisoned")
    );
    assert_eq!(
        num(&info_dotted, "num_mapped"),
        num(&info_plain, "num_mapped")
    );
    assert_eq!(rad_dotted, rad_plain);
    assert_eq!(num(&info_off, "num_poisoned"), 0);
    assert!(num(&info_off, "num_mapped") > num(&info_dotted, "num_mapped"));
}

/// Move a freshly built dotted index's core files to the paths piscem-rs
/// 0.11.0 gave them (`set_extension`), leaving the poison table and
/// signatures at the appended paths, as 0.11.0's `build-poison` and the
/// wrapper wrote them.
fn to_legacy_layout(prefix: &Path) {
    for ext in BUILD_ARTIFACT_EXTENSIONS {
        let from = index_artifact_path(prefix, ext);
        if from.exists() {
            std::fs::rename(&from, prefix.with_extension(ext)).unwrap();
        }
    }
}

#[test]
fn legacy_layout_index_loads_with_a_warning_and_maps_identically() {
    let w = Work::new();
    let refs = random_refs(4, "t", 20, 1500);
    let decoys = decoys_of(&refs.seqs[..8]);
    let decoy_fa = w.path("decoys.fa");
    write_fasta(&decoy_fa, &decoys);
    let reads = w.path("reads.fq");
    let mut seqs = refs.seqs.clone();
    seqs.extend(decoys);
    write_reads(&reads, &seqs, 600, 14);

    let prefix = w.path("idx/gencode.v49");
    w.build(&prefix, &refs);
    w.build_poison(&prefix, &decoy_fa);
    write_sigs(&prefix, refs.names.len());

    let modes: [&[&str]; 3] = [&[], &["--dict", "sshash"], &["--dict", "tiny"]];
    let before: Vec<_> = modes
        .iter()
        .enumerate()
        .map(|(i, m)| w.map(&prefix, &reads, &format!("new{i}"), m))
        .collect();
    assert!(num(&before[0].1, "num_poisoned") > 0);

    to_legacy_layout(&prefix);
    let mut expected = built_files("gencode");
    expected.extend(
        POISON_ARTIFACT_EXTENSIONS
            .iter()
            .chain(&["sigs.json"])
            .map(|e| format!("gencode.v49.{e}")),
    );
    expected.sort();
    assert_eq!(dir_listing(&w.path("idx")), expected);

    for (i, m) in modes.iter().enumerate() {
        let (rad, info, log) = w.map(&prefix, &reads, &format!("legacy{i}"), m);
        assert!(log.contains("pre-0.11.1 file layout"), "{m:?}:\n{log}");
        assert!(!log.contains("probably"), "{m:?}:\n{log}");
        assert_eq!(rad, before[i].0, "{m:?}: RAD differs after the move");
        for key in ["num_mapped", "num_poisoned", "num_refs"] {
            assert_eq!(num(&info, key), num(&before[i].1, key), "{m:?} {key}");
        }
        assert_eq!(info["signatures"], before[i].1["signatures"], "{m:?}");
    }

    // `stats` reads the same files.
    let stats = ok(
        piscem(&["stats".as_ref(), "-i".as_ref(), prefix.as_os_str()]),
        "stats",
    );
    let json: serde_json::Value =
        serde_json::from_str(&stats[stats.find("{\n").expect("stats JSON")..]).unwrap();
    for label in ["dictionary", "contig_table", "poison_table", "signatures"] {
        assert!(
            json["size_on_disk_bytes"][label].as_u64().is_some(),
            "{label}: {json}"
        );
    }

    // Core files that another index overwrote: the signature file disagrees.
    write_sigs(&prefix, 21);
    let (_, _, log) = w.map(&prefix, &reads, "mismatch", &[]);
    assert!(log.contains("probably overwritten"), "{log}");
}

/// A 0.11.0 dotted index whose decoys were built at the stem
/// (`build-poison -i idx/gencode`, the natural call once `ls` shows
/// `idx/gencode.ssi`, ...). 0.11.0 read that table, and so does this release.
/// A later `build-poison` at the dotted prefix writes the appended path,
/// which then wins, and leaves the stem's table alone.
#[test]
fn legacy_layout_keeps_decoys_built_at_the_stem() {
    let w = Work::new();
    let refs = random_refs(9, "t", 20, 1500);
    let decoys = decoys_of(&refs.seqs[..8]);
    let decoy_fa = w.path("decoys.fa");
    write_fasta(&decoy_fa, &decoys);
    let fewer_fa = w.path("fewer_decoys.fa");
    write_fasta(&fewer_fa, &decoys[..4]);
    let reads = w.path("reads.fq");
    let mut seqs = refs.seqs.clone();
    seqs.extend(decoys);
    write_reads(&reads, &seqs, 600, 17);

    let prefix = w.path("idx/gencode.v49");
    let stem = w.path("idx/gencode");
    w.build(&prefix, &refs);
    w.build_poison(&prefix, &decoy_fa);
    let modes: [&[&str]; 3] = [&[], &["--dict", "sshash"], &["--dict", "tiny"]];
    let before: Vec<_> = modes
        .iter()
        .enumerate()
        .map(|(i, m)| w.map(&prefix, &reads, &format!("new{i}"), m))
        .collect();
    let all_poisoned = num(&before[0].1, "num_poisoned");
    assert!(all_poisoned > 0);

    to_legacy_layout(&prefix);
    for ext in POISON_ARTIFACT_EXTENSIONS {
        std::fs::rename(
            index_artifact_path(&prefix, ext),
            index_artifact_path(&stem, ext),
        )
        .unwrap();
    }
    let mut expected = built_files("gencode");
    expected.extend(
        POISON_ARTIFACT_EXTENSIONS
            .iter()
            .map(|e| format!("gencode.{e}")),
    );
    expected.sort();
    assert_eq!(dir_listing(&w.path("idx")), expected);

    for (i, m) in modes.iter().enumerate() {
        let (rad, info, log) = w.map(&prefix, &reads, &format!("legacy{i}"), m);
        for key in ["num_poisoned", "num_mapped", "num_refs"] {
            assert_eq!(num(&info, key), num(&before[i].1, key), "{m:?} {key}");
        }
        assert_eq!(
            rad, before[i].0,
            "{m:?}: RAD differs with the stem's poison"
        );
        assert!(log.contains("pre-0.11.1 file layout"), "{m:?}:\n{log}");
        assert!(log.contains("with the decoy table"), "{m:?}:\n{log}");
    }
    let (_, info_off, _) = w.map(&prefix, &reads, "legacy_off", &["--no-poison"]);
    assert_eq!(num(&info_off, "num_poisoned"), 0);

    // `stats` reports the table the loader used.
    let stem_poison = std::fs::read(index_artifact_path(&stem, POISON_EXT)).unwrap();
    let stats = ok(
        piscem(&["stats".as_ref(), "-i".as_ref(), prefix.as_os_str()]),
        "stats",
    );
    let json: serde_json::Value =
        serde_json::from_str(&stats[stats.find("{\n").expect("stats JSON")..]).unwrap();
    assert_eq!(
        json["size_on_disk_bytes"]["poison_table"].as_u64(),
        Some(stem_poison.len() as u64),
        "{json}"
    );

    // New decoys for the dotted prefix go to its appended path and win; the
    // stem's table is untouched. They give what the same decoys give a
    // current-layout index.
    w.build_poison(&prefix, &fewer_fa);
    for ext in POISON_ARTIFACT_EXTENSIONS {
        assert!(index_artifact_path(&prefix, ext).is_file(), "{ext}");
    }
    assert_eq!(
        std::fs::read(index_artifact_path(&stem, POISON_EXT)).unwrap(),
        stem_poison
    );
    let (rad, info, log) = w.map(&prefix, &reads, "legacy_fewer", &[]);
    assert!(log.contains("pre-0.11.1 file layout"), "{log}");
    assert!(!log.contains("with the decoy table"), "{log}");
    let fewer_poisoned = num(&info, "num_poisoned");
    assert!(
        fewer_poisoned > 0 && fewer_poisoned < all_poisoned,
        "{fewer_poisoned} of {all_poisoned}"
    );
    let fresh = w.path("fresh/gencode.v49");
    w.build(&fresh, &refs);
    w.build_poison(&fresh, &fewer_fa);
    let (rad_fresh, info_fresh, _) = w.map(&fresh, &reads, "fresh_fewer", &[]);
    assert_eq!(num(&info_fresh, "num_poisoned"), fewer_poisoned);
    assert_eq!(rad, rad_fresh);
}

/// `-i dir/` looks for `dir.ssi` beside the directory, as it always has;
/// when nothing is there, the error says so and names the index inside.
#[test]
fn load_from_a_directory_prefix_names_the_index_inside() {
    let w = Work::new();
    let refs = random_refs(10, "t", 5, 500);
    std::fs::create_dir(w.path("d")).unwrap();
    w.build(&w.path("d/piscem_idx"), &refs);
    let reads = w.path("reads.fq");
    write_reads(&reads, &refs.seqs, 10, 18);
    let mut dir = w.path("d").into_os_string();
    dir.push("/");
    let dir = PathBuf::from(dir);
    for dict in ["sshash", "tiny"] {
        let res = piscem(&[
            "map-bulk".as_ref(),
            "-i".as_ref(),
            dir.as_os_str(),
            "-r".as_ref(),
            reads.as_os_str(),
            "-o".as_ref(),
            w.path("out").as_os_str(),
            "-t".as_ref(),
            "1".as_ref(),
            "--dict".as_ref(),
            dict.as_ref(),
        ]);
        let log = log_of(&res);
        assert!(!res.status.success(), "{dict}: mapping succeeded:\n{log}");
        assert!(log.contains("ends in a path separator"), "{dict}:\n{log}");
        let inside = format!("{:?}", dir.join("piscem_idx"));
        assert!(log.contains(&inside), "{dict}: no {inside}:\n{log}");
    }
}

/// A build at the new layout next to a complete legacy set at the same stem
/// is never mixed with it.
#[test]
fn current_layout_wins_over_legacy_files_at_the_stem() {
    let w = Work::new();
    let refs_old = random_refs(5, "old", 25, 1500);
    let refs_new = random_refs(6, "new", 20, 1500);
    let prefix = w.path("idx/gencode.v49");
    w.build(&prefix, &refs_old);
    to_legacy_layout(&prefix);
    w.build(&prefix, &refs_new);

    let reads = w.path("reads.fq");
    write_reads(&reads, &refs_new.seqs, 200, 15);
    let (_, info, log) = w.map(&prefix, &reads, "out", &[]);
    assert!(!log.contains("pre-0.11.1"), "{log}");
    assert_eq!(num(&info, "num_refs"), 20);
    assert_eq!(num(&info, "num_mapped"), 200);
}

/// `-i idx/plain.ssi` loaded `idx/plain` (and its poison table) in 0.11.0.
#[test]
fn prefix_naming_an_index_file_still_loads_that_index() {
    let w = Work::new();
    let refs = random_refs(7, "t", 20, 1500);
    let decoys = decoys_of(&refs.seqs[..8]);
    let decoy_fa = w.path("decoys.fa");
    write_fasta(&decoy_fa, &decoys);
    let reads = w.path("reads.fq");
    write_reads(&reads, &decoys, 200, 16);

    let prefix = w.path("idx/plain");
    w.build(&prefix, &refs);
    w.build_poison(&prefix, &decoy_fa);
    let (rad, info, _) = w.map(&prefix, &reads, "prefix", &[]);
    let (rad_file, info_file, log) =
        w.map(&index_artifact_path(&prefix, DICT_EXT), &reads, "file", &[]);
    assert!(!log.contains("WARN"), "{log}");
    assert!(num(&info, "num_poisoned") > 0);
    assert_eq!(num(&info_file, "num_poisoned"), num(&info, "num_poisoned"));
    assert_eq!(rad_file, rad);
}

#[test]
fn build_rejects_a_prefix_naming_a_directory() {
    let w = Work::new();
    let refs = random_refs(8, "t", 5, 500);
    let cfish = w.path("cfish");
    write_cfish(&cfish, &refs);
    std::fs::create_dir(w.path("myidx")).unwrap();
    for bad in ["myidx/", "myidx/.", "myidx/.."] {
        let out = w.dir.path().join(bad);
        let res = piscem(&[
            "build".as_ref(),
            "-i".as_ref(),
            cfish.as_os_str(),
            "-o".as_ref(),
            out.as_os_str(),
            "-k".as_ref(),
            "31".as_ref(),
            "-m".as_ref(),
            "19".as_ref(),
        ]);
        let log = log_of(&res);
        assert!(!res.status.success(), "{bad}: build succeeded:\n{log}");
        assert!(log.contains("must end in a file name"), "{bad}:\n{log}");
        let suggestion = format!("{:?}", out.join("index"));
        assert!(log.contains(&suggestion), "{bad}: no {suggestion}:\n{log}");
    }
    // Nothing written, inside the directory or beside it.
    assert!(dir_listing(&w.path("myidx")).is_empty());
    assert_eq!(
        dir_listing(w.dir.path()),
        ["cfish.cf_seg", "cfish.cf_seq", "cfish.json", "myidx"]
    );
}
