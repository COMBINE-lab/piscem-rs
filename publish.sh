#!/usr/bin/env bash
set -euo pipefail

die() {
    echo "error: $*" >&2
    exit 1
}

usage() {
    cat <<'EOF'
Usage:
  ./publish.sh <version> [--publish] [--dry-run]
  ./publish.sh [--publish] [--dry-run] <version>

Options:
  --publish  Publish piscem-rs after bumping and committing
  --dry-run  Show what would be done without modifying files, creating commits, tags, or publishing
  -h, --help Show this help message
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
            if [[ -n "$VERSION" ]]; then
                die "version specified more than once"
            fi
            VERSION="$1"
            ;;
    esac
    shift
done

[[ -n "$VERSION" ]] || {
    usage
    exit 1
}

if ! [[ "$VERSION" =~ ^[0-9]+\.[0-9]+\.[0-9]+([+-][0-9A-Za-z.-]+)*$ ]]; then
    die "version must look like X.Y.Z, optionally with prerelease/build suffixes"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

ROOT_CARGO="Cargo.toml"
LOCKFILE="Cargo.lock"
TAG="piscem-rs-v${VERSION}"

[[ -f "$ROOT_CARGO" ]] || die "not found: $ROOT_CARGO"
[[ -f "$LOCKFILE" ]] || die "not found: $LOCKFILE"

CURRENT_VERSION="$(sed -n 's/^version = "\(.*\)"/\1/p' "$ROOT_CARGO" | head -1)"
[[ -n "$CURRENT_VERSION" ]] || die "could not determine current crate version from $ROOT_CARGO"

if [[ "$CURRENT_VERSION" == "$VERSION" ]]; then
    die "crate version is already $VERSION"
fi

if git rev-parse "$TAG" >/dev/null 2>&1; then
    die "tag $TAG already exists"
fi

if [[ -n "$(git status --porcelain)" ]]; then
    die "working tree is not clean; commit or stash existing changes first"
fi

echo "Current crate version : $CURRENT_VERSION"
echo "New crate version     : $VERSION"
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

echo "Updating $ROOT_CARGO"
echo "  version: $CURRENT_VERSION -> $VERSION"

if [[ "$DRY_RUN" == false ]]; then
    sed -i.bak "1,/^version = /s/^version = \".*\"/version = \"${VERSION}\"/" "$ROOT_CARGO"
    rm -f "${ROOT_CARGO}.bak"
fi

UPDATED_VERSION="$(sed -n 's/^version = "\(.*\)"/\1/p' "$ROOT_CARGO" | head -1)"

if [[ "$DRY_RUN" == false ]]; then
    [[ "$UPDATED_VERSION" == "$VERSION" ]] || die "crate version update failed"
else
    echo "Dry-run: would rewrite $ROOT_CARGO and refresh $LOCKFILE"
fi

run cargo check -q
run git add "$ROOT_CARGO" "$LOCKFILE"
run git commit -m "chore(release): bump piscem-rs to v${VERSION}"

if [[ "$PUBLISH" == true ]]; then
    run cargo publish --allow-dirty
fi

run git tag -a "$TAG" -m "piscem-rs v${VERSION}"
run git push origin HEAD
run git push origin "$TAG"

if [[ "$DRY_RUN" == true ]]; then
    echo
    echo "Dry-run complete"
else
    echo
    echo "Release bump complete for v${VERSION}"
fi
