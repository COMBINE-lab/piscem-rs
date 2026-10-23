//! Where an index's files live on disk.
//!
//! An index is a set of files sharing a path prefix: `piscem build -o
//! idx/gencode.v49` writes `idx/gencode.v49.ssi`, `idx/gencode.v49.ctab`, and
//! so on. Every file name is the **full prefix** with `.<ext>` appended, so a
//! dot inside the prefix's file name is part of the name, never an extension
//! to replace. [`index_artifact_path`] is the one place that rule lives; the
//! piscem wrapper and other callers should use it too rather than
//! concatenating strings, so that all the files of an index agree.
//!
//! # The pre-0.11.1 layout
//!
//! piscem-rs 0.11.0 and earlier built the core files (everything `build`
//! writes, [`BUILD_ARTIFACT_EXTENSIONS`]) with [`PathBuf::set_extension`], which
//! *replaces* the prefix's last dotted component: `idx/gencode.v49` got
//! `idx/gencode.ssi`, `idx/gencode.ctab`, ... . The files written by others
//! were appended to the prefix they were given, as now: `build-poison`'s
//! `.poison` / `.poison.json` and the wrapper's `.sigs.json`. For a prefix
//! without a dot in its file name the two rules give the same paths, so only
//! dotted prefixes are affected.
//!
//! [`IndexPaths::resolve`] still loads such an index. It picks one layout per
//! prefix, never mixing core files from both. The signatures are read from
//! the appended path, where the wrapper wrote them. The poison table is read
//! from the appended path if one is there (`build-poison -i idx/gencode.v49`
//! wrote it there, and 0.11.0 never read it), and otherwise from the stem
//! (`idx/gencode.poison`, written by `build-poison -i idx/gencode` and read
//! by 0.11.0), so an old index keeps the decoys 0.11.0 applied.

use std::path::{Path, PathBuf};

use anyhow::{Result, bail};

/// sshash dictionary. sshash writes its MPHF beside it as [`DICT_MPHF_EXT`].
pub const DICT_EXT: &str = "ssi";
/// sshash MPHF container, derived by sshash from the [`DICT_EXT`] path.
pub const DICT_MPHF_EXT: &str = "ssi.mphf";
/// Contig table.
pub const CTAB_EXT: &str = "ctab";
/// Reference names and lengths.
pub const REFINFO_EXT: &str = "refinfo";
/// Equivalence class table (`build --build-ec-table`).
pub const ECTAB_EXT: &str = "ectab";
/// Tiny dictionary (`build --dict tiny`, or `auto` below the size threshold).
pub const TDCT_EXT: &str = "tdct";
/// Tiny contig table, written alongside [`TDCT_EXT`].
pub const TCT_EXT: &str = "tct";
/// Poison (decoy) k-mer table, written by `build-poison`.
pub const POISON_EXT: &str = "poison";
/// Poison table statistics, written by `build-poison`.
pub const POISON_JSON_EXT: &str = "poison.json";
/// Reference signatures. piscem-rs only reads this file; the piscem wrapper
/// writes it when it builds from FASTA.
pub const SIGS_EXT: &str = "sigs.json";

/// Extensions of the files `build` may write for a prefix. The optional ones
/// ([`ECTAB_EXT`], [`TDCT_EXT`], [`TCT_EXT`]) are present only when enabled;
/// a caller cleaning up before a rebuild should remove whichever exist.
pub const BUILD_ARTIFACT_EXTENSIONS: &[&str] = &[
    DICT_EXT,
    DICT_MPHF_EXT,
    CTAB_EXT,
    REFINFO_EXT,
    ECTAB_EXT,
    TDCT_EXT,
    TCT_EXT,
];

/// Extensions of the files `build-poison` writes for a prefix. A poison table
/// left over from an earlier build would be applied to a rebuilt index, so a
/// rebuild without decoys should remove these too.
pub const POISON_ARTIFACT_EXTENSIONS: &[&str] = &[POISON_EXT, POISON_JSON_EXT];

/// The path of an index file: `prefix` with `.ext` appended.
///
/// `idx/gencode.v49` with `ctab` is `idx/gencode.v49.ctab`, and
/// `idx/piscem_idx` with `ssi.mphf` is `idx/piscem_idx.ssi.mphf`. A trailing
/// separator is dropped (`myidx/` with `ssi` is `myidx.ssi`), as
/// [`Path::with_added_extension`] does.
///
/// A prefix that names no file (empty, `/`, `..`) is returned unchanged;
/// check prefixes with [`validate_index_prefix`] before writing to them.
pub fn index_artifact_path(prefix: &Path, ext: &str) -> PathBuf {
    prefix.with_added_extension(ext)
}

/// Check that `prefix` can name an index, before anything is written.
///
/// Index files are named `<prefix>.<ext>`, so the prefix must end in a file
/// name. Rejected, with a message suggesting `<prefix>/index`:
/// - an empty prefix;
/// - a prefix ending in a path separator (`myidx/`), which would put the
///   files beside the directory (`myidx.ssi`), not inside it;
/// - a prefix whose last component is `.` or `..`.
///
/// Every other prefix is accepted, dots and all.
pub fn validate_index_prefix(prefix: &Path) -> Result<()> {
    let raw = prefix.as_os_str().as_encoded_bytes();
    let last = last_component_bytes(prefix);
    let names_a_dir = raw.is_empty()
        || last.is_empty() // ends with a separator (or is one)
        || last == b"."
        || last == b".."
        || prefix.file_name().is_none();
    if names_a_dir {
        bail!(
            "invalid index prefix {:?}: it must end in a file name, because the index is \
             written as <prefix>.ssi, <prefix>.ctab, ...; to put the index inside that \
             directory, pass {:?}",
            prefix,
            prefix.join("index")
        );
    }
    Ok(())
}

/// The bytes after the last path separator of `prefix`, as typed (empty if
/// it ends in one). Unlike [`Path::file_name`], this sees a trailing `/`.
fn last_component_bytes(prefix: &Path) -> &[u8] {
    prefix
        .as_os_str()
        .as_encoded_bytes()
        .rsplit(|&b| std::path::is_separator(b as char))
        .next()
        .unwrap_or_default()
}

/// For a load that found no index at a prefix ending in a path separator
/// (`-i dir/`, which looks for `dir.ssi` beside the directory, as it always
/// has): a note saying so, naming the indexes inside the directory if there
/// are any. `None` for any other prefix.
pub(crate) fn directory_prefix_hint(prefix: &Path) -> Option<String> {
    let raw = prefix.as_os_str().as_encoded_bytes();
    if raw.is_empty() || !last_component_bytes(prefix).is_empty() {
        return None;
    }
    let suffix = format!(".{REFINFO_EXT}");
    let mut inside: Vec<PathBuf> = std::fs::read_dir(prefix)
        .into_iter()
        .flatten()
        .filter_map(|e| {
            let name = e.ok()?.file_name();
            let stem = name.to_str()?.strip_suffix(&suffix)?;
            (!stem.is_empty()).then(|| prefix.join(stem))
        })
        .collect();
    inside.sort();
    let what = match inside.as_slice() {
        [] => format!(
            "pass the prefix the index was built with, e.g. {:?}",
            prefix.join("index")
        ),
        found => {
            let names: Vec<String> = found.iter().take(3).map(|p| format!("{p:?}")).collect();
            let noun = if found.len() == 1 { "index" } else { "indexes" };
            format!("the directory holds the {noun} {}", names.join(", "))
        }
    };
    Some(format!(
        "the index prefix {prefix:?} ends in a path separator, so its files are looked for \
         beside that directory, not inside it; {what}"
    ))
}

/// Which naming rule an index's core files follow.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum IndexLayout {
    /// `<prefix>.<ext>` for every file (piscem-rs 0.11.1 and later, and any
    /// prefix without a dot in its file name).
    Current,
    /// Core files at `<prefix>` with its last extension replaced, as
    /// piscem-rs 0.11.0 and earlier wrote them for a dotted prefix.
    Legacy,
}

/// The on-disk files of one index, with the layout decided once.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IndexPaths {
    /// The index's prefix; the signature file is always appended to it.
    prefix: PathBuf,
    /// Prefix the core files are appended to.
    core_prefix: PathBuf,
    /// Prefix the poison table and its stats are appended to.
    poison_prefix: PathBuf,
    /// The prefix as given, when it differs from `prefix`.
    given: Option<PathBuf>,
}

impl IndexPaths {
    /// The paths of an index at `prefix` in the current layout, without
    /// looking at the disk (use this to write an index).
    pub fn new(prefix: &Path) -> Self {
        Self {
            prefix: prefix.to_path_buf(),
            core_prefix: prefix.to_path_buf(),
            poison_prefix: prefix.to_path_buf(),
            given: None,
        }
    }

    /// Find the index at `prefix` on disk.
    ///
    /// In order:
    /// 1. If the core files are at `<prefix>.<ext>`, that is the index.
    /// 2. If `prefix` names one of an index's files (`idx/gencode.ssi`, which
    ///    0.11.0 accepted in place of the prefix) and the index at the shorter
    ///    prefix is complete, that index is used, poison table included.
    /// 3. If `prefix` has an extension, no core file exists at its appended
    ///    paths, and a complete set exists at the replaced paths, the index is
    ///    in the [`IndexLayout::Legacy`] layout: its core files are read
    ///    from there, and its signatures from `<prefix>.sigs.json`, where the
    ///    wrapper wrote them. The poison table and its stats are read as a
    ///    pair from `<prefix>.poison` if that exists (`build-poison` given the
    ///    dotted prefix wrote it there), and otherwise from the replaced path
    ///    `<stem>.poison` (`build-poison` given the stem wrote it there, and
    ///    0.11.0 read it from there). Either was built against these same
    ///    core files.
    /// 4. Otherwise the current layout, so a missing index is reported at the
    ///    paths `build` would have written.
    ///
    /// Only regular files (or links to them) count. This only inspects the
    /// disk; [`IndexPaths::log_resolution`] reports cases 2 and 3.
    pub fn resolve(prefix: &Path) -> Self {
        let current = Self::new(prefix);
        if current.has_complete_core() {
            return current;
        }
        if prefix.extension().is_none() {
            return current;
        }
        // `with_extension("")` drops the last extension: `idx/gencode.v49`
        // -> `idx/gencode`, whose appended paths are exactly the paths
        // `set_extension` gave (`idx/gencode.ctab`).
        let stem = prefix.with_extension("");
        if stem.file_name().is_none() {
            // `odd/..v50` -> `odd/..`: no file can be named after it.
            return current;
        }
        let named_file = prefix
            .extension()
            .and_then(|e| e.to_str())
            .is_some_and(|e| FILE_NAME_EXTENSIONS.contains(&e));
        if named_file {
            let at_stem = Self::new(&stem);
            if at_stem.has_complete_core() {
                return Self {
                    given: Some(prefix.to_path_buf()),
                    ..at_stem
                };
            }
        }
        let mut legacy = Self {
            prefix: prefix.to_path_buf(),
            core_prefix: stem.clone(),
            poison_prefix: prefix.to_path_buf(),
            given: None,
        };
        if current.has_any_core() || !legacy.has_complete_core() {
            return current;
        }
        if !index_artifact_path(prefix, POISON_EXT).is_file()
            && index_artifact_path(&stem, POISON_EXT).is_file()
        {
            legacy.poison_prefix = stem;
        }
        legacy
    }

    /// The prefix of the index (after resolving a prefix that named a file).
    pub fn prefix(&self) -> &Path {
        &self.prefix
    }

    /// Which layout the core files follow.
    pub fn layout(&self) -> IndexLayout {
        if self.core_prefix == self.prefix {
            IndexLayout::Current
        } else {
            IndexLayout::Legacy
        }
    }

    /// Path the index file with extension `ext` (e.g. [`CTAB_EXT`]) is read
    /// from.
    ///
    /// A file written for an existing index, such as a new poison table,
    /// belongs at `index_artifact_path(paths.prefix(), ext)`. For a
    /// [`IndexLayout::Legacy`] index that can differ from this path: a new
    /// poison table goes to `<prefix>.poison`, which a later load prefers,
    /// and the stem's table is left alone.
    pub fn file(&self, ext: &str) -> PathBuf {
        match ext {
            POISON_EXT | POISON_JSON_EXT => index_artifact_path(&self.poison_prefix, ext),
            SIGS_EXT => index_artifact_path(&self.prefix, ext),
            _ => index_artifact_path(&self.core_prefix, ext),
        }
    }

    /// Whether both Tiny artifacts (`.tdct` + `.tct`) exist.
    pub fn tiny_artifacts_exist(&self) -> bool {
        self.file(TDCT_EXT).is_file() && self.file(TCT_EXT).is_file()
    }

    /// Log how the prefix was resolved, if it was not used as given: a warning
    /// for a [`IndexLayout::Legacy`] index, a note for a prefix that named
    /// an index file. Call once per load.
    pub fn log_resolution(&self) {
        if let Some(given) = &self.given {
            tracing::info!(
                "{} names an index file; loading the index at prefix {}",
                given.display(),
                self.prefix.display()
            );
        }
        if self.layout() == IndexLayout::Legacy {
            let core = self.core_prefix.display();
            // Name the stem's poison table when it is the one used: it is not
            // where a current-layout index keeps it.
            let decoys = if self.poison_prefix == self.core_prefix {
                format!(
                    ", with the decoy table {}, as piscem-rs <= 0.11.0 did (there is no {})",
                    self.file(POISON_EXT).display(),
                    index_artifact_path(&self.prefix, POISON_EXT).display(),
                )
            } else {
                String::new()
            };
            tracing::warn!(
                "index {prefix} uses the pre-0.11.1 file layout (piscem-rs <= 0.11.0): its \
                 core files are {core}.ssi, {core}.ctab, ... rather than {prefix}.ssi, ...; \
                 files at these paths are shared with the index at {core} and with any \
                 {core}.* index built by piscem-rs <= 0.11.0. Loading them{decoys}; rebuild \
                 the index to move it to the current layout",
                prefix = self.prefix.display(),
            );
        }
    }

    /// Reference info plus one complete dictionary backend.
    fn has_complete_core(&self) -> bool {
        let present = |ext| self.file(ext).is_file();
        present(REFINFO_EXT)
            && ((present(DICT_EXT) && present(DICT_MPHF_EXT) && present(CTAB_EXT))
                || (present(TDCT_EXT) && present(TCT_EXT)))
    }

    fn has_any_core(&self) -> bool {
        BUILD_ARTIFACT_EXTENSIONS
            .iter()
            .any(|ext| self.file(ext).is_file())
    }
}

/// Extensions that, ending a prefix, name one of the index's own files
/// (`-i idx/gencode.ssi`).
const FILE_NAME_EXTENSIONS: &[&str] = &[
    DICT_EXT,
    CTAB_EXT,
    REFINFO_EXT,
    ECTAB_EXT,
    TDCT_EXT,
    TCT_EXT,
    POISON_EXT,
];

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    fn p(s: &str) -> &Path {
        Path::new(s)
    }

    #[test]
    fn artifact_path_appends_to_the_full_prefix() {
        assert_eq!(
            index_artifact_path(p("idx/gencode.v49"), "ctab"),
            p("idx/gencode.v49.ctab")
        );
        assert_eq!(
            index_artifact_path(p("idx/gencode.v49"), DICT_MPHF_EXT),
            p("idx/gencode.v49.ssi.mphf")
        );
        assert_eq!(
            index_artifact_path(p("a.b.c"), POISON_JSON_EXT),
            p("a.b.c.poison.json")
        );
        assert_eq!(index_artifact_path(p("foo."), "ssi"), p("foo..ssi"));
        assert_eq!(index_artifact_path(p(".hidden"), "ssi"), p(".hidden.ssi"));
        assert_eq!(index_artifact_path(p("myidx/"), "ssi"), p("myidx.ssi"));
    }

    /// For a prefix without a dot in its file name the appended path is the
    /// one 0.11.0 wrote with `set_extension`, so those indexes keep their
    /// file names. Dots in directory names never mattered.
    #[test]
    fn undotted_prefix_paths_match_set_extension() {
        for prefix in ["piscem_idx", "out/piscem_idx", "/abs/v1.2/index", "./x"] {
            for ext in BUILD_ARTIFACT_EXTENSIONS
                .iter()
                .chain(POISON_ARTIFACT_EXTENSIONS)
                .chain(&[SIGS_EXT])
            {
                assert_eq!(
                    index_artifact_path(p(prefix), ext),
                    p(prefix).with_extension(ext),
                    "{prefix} + {ext}"
                );
            }
        }
    }

    #[cfg(unix)]
    #[test]
    fn artifact_path_keeps_non_utf8_bytes() {
        use std::os::unix::ffi::OsStrExt;
        let prefix = Path::new(std::ffi::OsStr::from_bytes(b"idx/g\xffn.v1"));
        let got = index_artifact_path(prefix, "ctab");
        assert_eq!(got.as_os_str().as_bytes(), b"idx/g\xffn.v1.ctab");
    }

    #[test]
    fn validate_accepts_prefixes_that_name_a_file() {
        for ok in [
            "index",
            "piscem_idx",
            "idx/gencode.v49",
            "./idx",
            "../idx/x.y",
            "/abs/index",
            ".hidden",
            "foo.",
        ] {
            validate_index_prefix(p(ok)).unwrap_or_else(|e| panic!("{ok}: {e:#}"));
        }
    }

    #[test]
    fn validate_rejects_prefixes_that_name_a_directory() {
        for (bad, suggestion) in [
            ("", "index"),
            ("myidx/", "myidx/index"),
            ("myidx//", "myidx//index"),
            ("/", "/index"),
            (".", "./index"),
            ("..", "../index"),
            ("out/.", "out/./index"),
            ("out/..", "out/../index"),
        ] {
            let err = validate_index_prefix(p(bad)).expect_err(bad);
            let msg = format!("{err:#}");
            assert!(msg.contains("must end in a file name"), "{bad}: {msg}");
            assert!(msg.contains(&format!("{suggestion:?}")), "{bad}: {msg}");
        }
    }

    fn touch(dir: &Path, names: &[&str]) {
        for n in names {
            std::fs::write(dir.join(n), b"x").unwrap();
        }
    }

    const CORE: &[&str] = &["ssi", "ssi.mphf", "ctab", "refinfo"];

    fn touch_core(dir: &Path, stem: &str) {
        let names: Vec<String> = CORE.iter().map(|e| format!("{stem}.{e}")).collect();
        touch(dir, &names.iter().map(String::as_str).collect::<Vec<_>>());
    }

    #[test]
    fn resolve_current_layout() {
        let d = tempdir().unwrap();
        touch_core(d.path(), "g.v49");
        touch(d.path(), &["g.v49.poison", "g.v49.tdct", "g.v49.tct"]);
        let paths = IndexPaths::resolve(&d.path().join("g.v49"));
        assert_eq!(paths.layout(), IndexLayout::Current);
        assert_eq!(paths.file(CTAB_EXT), d.path().join("g.v49.ctab"));
        assert_eq!(paths.file(POISON_EXT), d.path().join("g.v49.poison"));
        assert!(paths.tiny_artifacts_exist());
    }

    #[test]
    fn resolve_prefers_current_over_legacy() {
        let d = tempdir().unwrap();
        touch_core(d.path(), "g.v49");
        touch_core(d.path(), "g");
        let paths = IndexPaths::resolve(&d.path().join("g.v49"));
        assert_eq!(paths.layout(), IndexLayout::Current);
        assert_eq!(paths.file(DICT_EXT), d.path().join("g.v49.ssi"));
    }

    /// Poison at both paths: the appended one (`build-poison -i g.v49`) wins.
    #[test]
    fn resolve_legacy_layout_prefers_the_appended_poison_table() {
        let d = tempdir().unwrap();
        touch_core(d.path(), "g");
        touch(d.path(), &["g.tdct", "g.tct", "g.v49.poison", "g.poison"]);
        let paths = IndexPaths::resolve(&d.path().join("g.v49"));
        assert_eq!(paths.layout(), IndexLayout::Legacy);
        assert_eq!(paths.prefix(), d.path().join("g.v49"));
        for ext in BUILD_ARTIFACT_EXTENSIONS {
            assert_eq!(
                paths.file(ext),
                d.path().join("g.v49").with_extension(ext),
                "{ext} at the 0.11.0 path"
            );
        }
        assert!(paths.tiny_artifacts_exist());
        assert_eq!(paths.file(POISON_EXT), d.path().join("g.v49.poison"));
        assert_eq!(
            paths.file(POISON_JSON_EXT),
            d.path().join("g.v49.poison.json")
        );
        assert_eq!(paths.file(SIGS_EXT), d.path().join("g.v49.sigs.json"));
    }

    /// Poison only at the stem (`build-poison -i g`, read by 0.11.0): kept,
    /// as a pair; the signatures stay at the appended path.
    #[test]
    fn resolve_legacy_layout_falls_back_to_the_stem_poison_table() {
        let d = tempdir().unwrap();
        touch_core(d.path(), "g");
        touch(d.path(), &["g.poison", "g.poison.json", "g.sigs.json"]);
        let paths = IndexPaths::resolve(&d.path().join("g.v49"));
        assert_eq!(paths.layout(), IndexLayout::Legacy);
        assert_eq!(paths.file(POISON_EXT), d.path().join("g.poison"));
        assert_eq!(paths.file(POISON_JSON_EXT), d.path().join("g.poison.json"));
        assert_eq!(paths.file(SIGS_EXT), d.path().join("g.v49.sigs.json"));
    }

    /// No poison table anywhere: reported at the appended path.
    #[test]
    fn resolve_legacy_layout_without_poison_uses_the_appended_path() {
        let d = tempdir().unwrap();
        touch_core(d.path(), "g");
        let paths = IndexPaths::resolve(&d.path().join("g.v49"));
        assert_eq!(paths.layout(), IndexLayout::Legacy);
        assert_eq!(paths.file(POISON_EXT), d.path().join("g.v49.poison"));
        assert_eq!(
            paths.file(POISON_JSON_EXT),
            d.path().join("g.v49.poison.json")
        );
    }

    /// The poison fallback only applies to a legacy index: a current-layout
    /// index never reads the stem's table.
    #[test]
    fn resolve_current_layout_never_reads_the_stem_poison_table() {
        let d = tempdir().unwrap();
        touch_core(d.path(), "g.v49");
        touch(d.path(), &["g.poison"]);
        let paths = IndexPaths::resolve(&d.path().join("g.v49"));
        assert_eq!(paths.layout(), IndexLayout::Current);
        assert_eq!(paths.file(POISON_EXT), d.path().join("g.v49.poison"));
    }

    /// Directories are not index files, and a stem that names no file
    /// (`odd/..v50` -> `odd/..`) is never a legacy index.
    #[test]
    fn resolve_ignores_directories() {
        let d = tempdir().unwrap();
        let odd = d.path().join("odd");
        std::fs::create_dir(&odd).unwrap();
        let paths = IndexPaths::resolve(&odd.join("..v50"));
        assert_eq!(paths, IndexPaths::new(&odd.join("..v50")));
        assert!(!paths.tiny_artifacts_exist());

        for e in CORE {
            std::fs::create_dir(d.path().join(format!("g.{e}"))).unwrap();
        }
        let paths = IndexPaths::resolve(&d.path().join("g.v49"));
        assert_eq!(paths.layout(), IndexLayout::Current);
    }

    #[test]
    fn resolve_legacy_tiny_only_layout() {
        let d = tempdir().unwrap();
        touch(d.path(), &["g.refinfo", "g.tdct", "g.tct"]);
        let paths = IndexPaths::resolve(&d.path().join("g.v49"));
        assert_eq!(paths.layout(), IndexLayout::Legacy);
    }

    /// A partial current-layout index is never completed from legacy files.
    #[test]
    fn resolve_never_mixes_layouts() {
        let d = tempdir().unwrap();
        touch_core(d.path(), "g");
        touch(d.path(), &["g.v49.tdct"]);
        let paths = IndexPaths::resolve(&d.path().join("g.v49"));
        assert_eq!(paths.layout(), IndexLayout::Current);
        assert_eq!(paths.file(DICT_EXT), d.path().join("g.v49.ssi"));
    }

    #[test]
    fn resolve_incomplete_legacy_reports_current_paths() {
        let d = tempdir().unwrap();
        touch(d.path(), &["g.ssi", "g.ctab"]); // no .refinfo / .ssi.mphf
        let paths = IndexPaths::resolve(&d.path().join("g.v49"));
        assert_eq!(paths.layout(), IndexLayout::Current);
        assert_eq!(paths.file(REFINFO_EXT), d.path().join("g.v49.refinfo"));
    }

    #[test]
    fn resolve_undotted_prefix_never_falls_back() {
        let d = tempdir().unwrap();
        let paths = IndexPaths::resolve(&d.path().join("piscem_idx"));
        assert_eq!(paths, IndexPaths::new(&d.path().join("piscem_idx")));
    }

    /// `-i idx/gencode.ssi` loaded `idx/gencode`, poison table included, in
    /// 0.11.0; it still does.
    #[test]
    fn resolve_prefix_naming_an_index_file() {
        let d = tempdir().unwrap();
        touch_core(d.path(), "g");
        touch(d.path(), &["g.poison"]);
        for given in ["g.ssi", "g.ctab", "g.refinfo"] {
            let paths = IndexPaths::resolve(&d.path().join(given));
            assert_eq!(paths.layout(), IndexLayout::Current, "{given}");
            assert_eq!(paths.prefix(), d.path().join("g"), "{given}");
            assert_eq!(paths.file(POISON_EXT), d.path().join("g.poison"));
        }
    }

    #[test]
    fn directory_prefix_hint_only_for_a_trailing_separator() {
        let d = tempdir().unwrap();
        let dir = d.path().join("idxdir");
        std::fs::create_dir(&dir).unwrap();
        assert_eq!(directory_prefix_hint(&dir), None);
        assert_eq!(directory_prefix_hint(&dir.join("x.v1")), None);
        assert_eq!(directory_prefix_hint(p("")), None);

        let mut slash = dir.clone().into_os_string();
        slash.push("/");
        let slash = PathBuf::from(slash);
        let hint = directory_prefix_hint(&slash).unwrap();
        assert!(hint.contains("ends in a path separator"), "{hint}");
        assert!(
            hint.contains(&format!("{:?}", slash.join("index"))),
            "{hint}"
        );

        touch(&dir, &["piscem_idx.refinfo", "piscem_idx.ssi", "notes.txt"]);
        let hint = directory_prefix_hint(&slash).unwrap();
        assert!(
            hint.contains(&format!("{:?}", slash.join("piscem_idx"))),
            "{hint}"
        );
    }

    #[test]
    fn new_is_the_current_layout() {
        let paths = IndexPaths::new(p("idx/gencode.v49"));
        assert_eq!(paths.layout(), IndexLayout::Current);
        for ext in BUILD_ARTIFACT_EXTENSIONS
            .iter()
            .chain(POISON_ARTIFACT_EXTENSIONS)
        {
            assert_eq!(
                paths.file(ext),
                index_artifact_path(p("idx/gencode.v49"), ext)
            );
        }
    }
}
