#!/usr/bin/env bash
# publish.sh — publish piscem-rs to crates.io and create a git tag.
#
# Usage:
#   ./publish.sh              # dry-run (no publish, no tag)
#   ./publish.sh --publish    # publish to crates.io, commit tag, push
#
# The version is read from Cargo.toml — bump it before running this script.

set -euo pipefail

# ── helpers ──────────────────────────────────────────────────────────────────

die() { echo "error: $*" >&2; exit 1; }

# ── args ─────────────────────────────────────────────────────────────────────

PUBLISH=false
[[ "${1-}" == "--publish" ]] && PUBLISH=true

# ── setup ────────────────────────────────────────────────────────────────────

cd "$(dirname "$0")"   # run from repo root regardless of cwd

[[ -f "Cargo.toml" ]] || die "not found: Cargo.toml"

VERSION=$(grep '^version' Cargo.toml | head -1 | sed 's/.*"\(.*\)"/\1/')
TAG="piscem-rs-v${VERSION}"

echo "piscem-rs version: $VERSION"
echo "tag:               $TAG"
echo ""

# ── pre-flight checks ───────────────────────────────────────────────────────

if git rev-parse "$TAG" &>/dev/null; then
    die "tag $TAG already exists — bump the version in Cargo.toml first"
fi

if ! git diff --quiet; then
    die "working tree has unstaged changes — commit or stash them first"
fi

if ! git diff --cached --quiet; then
    die "index has staged changes — commit or stash them first"
fi

# ── dry run ──────────────────────────────────────────────────────────────────

echo "Running cargo publish --dry-run ..."
echo ""
cargo publish --dry-run
echo ""
echo "Dry run passed."

if ! $PUBLISH; then
    echo ""
    echo "This was a dry run. To publish for real:"
    echo "  ./publish.sh --publish"
    exit 0
fi

# ── publish ──────────────────────────────────────────────────────────────────

echo ""
echo "Publishing piscem-rs $VERSION to crates.io ..."
cargo publish

# ── tag and push ─────────────────────────────────────────────────────────────

git tag -a "$TAG" -m "piscem-rs v${VERSION}"
echo "Created tag $TAG"

git push origin "$TAG"
echo "Pushed tag to origin."
