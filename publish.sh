#!/usr/bin/env bash
set -euo pipefail

die() {
    echo "error: $*" >&2
    exit 1
}

usage() {
    cat <<'EOF'
Usage:
  ./publish.sh <crate> <version> [--publish] [--dry-run]

Arguments:
  <crate>    Crate to release: seq_geom_parser or piscem-rs
  <version>  New version (X.Y.Z format)

Options:
  --publish  Publish to crates.io after bumping and committing
  --dry-run  Show what would be done without modifying anything
  -h, --help Show this help message

Examples:
  ./publish.sh seq_geom_parser 1.0.0 --publish
  ./publish.sh piscem-rs 0.3.0 --publish
  ./publish.sh piscem-rs 0.3.0 --dry-run
EOF
}

print_cmd() {
    printf '+'
    printf ' %q' "$@"
    printf '\n'
}

run() {
    print_cmd "$@"
    if [[ "$DRY_RUN" == true ]]; then
        return 0
    fi
    "$@"
}

CRATE=""
VERSION=""
PUBLISH=false
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --publish)
            PUBLISH=true
            ;;
        --dry-run)
            DRY_RUN=true
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        -*)
            die "unknown option: $1"
            ;;
        *)
            if [[ -z "$CRATE" ]]; then
                CRATE="$1"
            elif [[ -z "$VERSION" ]]; then
                VERSION="$1"
            else
                die "too many positional arguments"
            fi
            ;;
    esac
    shift
done

[[ -n "$CRATE" ]] || { usage; exit 1; }
[[ -n "$VERSION" ]] || { usage; exit 1; }

if ! [[ "$VERSION" =~ ^[0-9]+\.[0-9]+\.[0-9]+([+-][0-9A-Za-z.-]+)*$ ]]; then
    die "version must look like X.Y.Z, optionally with prerelease/build suffixes"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Determine the Cargo.toml path for the crate
case "$CRATE" in
    piscem-rs)
        CARGO_TOML="Cargo.toml"
        TAG="${CRATE}-v${VERSION}"
        PUBLISH_ARGS=""
        ;;
    seq_geom_parser)
        CARGO_TOML="crates/seq_geom_parser/Cargo.toml"
        TAG="${CRATE}-v${VERSION}"
        PUBLISH_ARGS="-p seq_geom_parser"
        ;;
    *)
        die "unknown crate: $CRATE (expected: piscem-rs or seq_geom_parser)"
        ;;
esac

LOCKFILE="Cargo.lock"

[[ -f "$CARGO_TOML" ]] || die "not found: $CARGO_TOML"
[[ -f "$LOCKFILE" ]] || die "not found: $LOCKFILE"

CURRENT_VERSION="$(sed -n 's/^version = "\(.*\)"/\1/p' "$CARGO_TOML" | head -1)"
[[ -n "$CURRENT_VERSION" ]] || die "could not determine current crate version from $CARGO_TOML"

if [[ "$CURRENT_VERSION" == "$VERSION" ]]; then
    die "crate version is already $VERSION"
fi

if git rev-parse "$TAG" >/dev/null 2>&1; then
    die "tag $TAG already exists"
fi

if [[ -n "$(git status --porcelain)" ]]; then
    die "working tree is not clean; commit or stash existing changes first"
fi

echo "Crate                 : $CRATE"
echo "Cargo.toml            : $CARGO_TOML"
echo "Current version       : $CURRENT_VERSION"
echo "New version           : $VERSION"
echo "Tag                   : $TAG"
if [[ "$PUBLISH" == true ]]; then
    echo "Publish crate         : yes"
else
    echo "Publish crate         : no"
fi
if [[ "$DRY_RUN" == true ]]; then
    echo "Dry-run               : yes"
else
    echo "Dry-run               : no"
fi
echo

echo "Updating $CARGO_TOML"
echo "  version: $CURRENT_VERSION -> $VERSION"

if [[ "$DRY_RUN" == false ]]; then
    sed -i.bak "1,/^version = /s/^version = \".*\"/version = \"${VERSION}\"/" "$CARGO_TOML"
    rm -f "${CARGO_TOML}.bak"
fi

UPDATED_VERSION="$(sed -n 's/^version = "\(.*\)"/\1/p' "$CARGO_TOML" | head -1)"

if [[ "$DRY_RUN" == false ]]; then
    [[ "$UPDATED_VERSION" == "$VERSION" ]] || die "crate version update failed"
else
    echo "Dry-run: would rewrite $CARGO_TOML and refresh $LOCKFILE"
fi

run cargo check -q
run git add "$CARGO_TOML" "$LOCKFILE"
run git commit -m "chore(release): bump ${CRATE} to v${VERSION}"

if [[ "$PUBLISH" == true ]]; then
    run cargo publish $PUBLISH_ARGS --allow-dirty
fi

run git tag -a "$TAG" -m "${CRATE} v${VERSION}"
run git push origin HEAD
run git push origin "$TAG"

if [[ "$DRY_RUN" == true ]]; then
    echo
    echo "Dry-run complete"
else
    echo
    echo "Release bump complete for ${CRATE} v${VERSION}"
fi
