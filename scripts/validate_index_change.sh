#!/usr/bin/env bash
# End-to-end validation of the hashbrown/eq_classes change against the commit
# immediately before it (7ea5b84). Every stage compares NEW against OLD.
#
# Stage 0 is the control: if two builds from the SAME binary differ, then
# byte-comparison cannot isolate anything and later stages are meaningless.
set -uo pipefail

V=${VALIDATE_DIR:-/tmp/validate-index-change}
mkdir -p "$V"
NEW=$V/piscem-new
OLD=$V/piscem-old
REF=/scratch1/rob/flex_bench/human_ref/probes.fa
G='1{b[16]u[12]x[0-3]f[TTGCTAGGACCG]s[10]x:}2{r:}'
R1=/scratch2/rob/read_data/10x_flexv2/23C235LT4_8_0470676772_flexv2_gex_S2_L008_R1_001.fastq.gz
R2=/scratch2/rob/read_data/10x_flexv2/23C235LT4_8_0470676772_flexv2_gex_S2_L008_R2_001.fastq.gz
T1=/scratch1/rob/flex_bench/tiny/R1.fq.gz
T2=/scratch1/rob/flex_bench/tiny/R2.fq.gz

# Data artifacts only: build_time.txt and index_cfish.json embed paths/timings.
ARTIFACTS="index.ctab index.ectab index.refinfo index.ssi index.ssi.mphf index.tct index.tdct"

FAIL=0
note() { echo; echo "======== $* ========"; }
ok()   { echo "  PASS  $*"; }
bad()  { echo "  FAIL  $*"; FAIL=1; }

cmp_idx() { # dirA dirB label
  local a=$1 b=$2 label=$3 diffs=0
  for f in $ARTIFACTS; do
    if [[ ! -f $a/$f ]]; then bad "$label: $f missing from $a"; diffs=1; continue; fi
    if [[ ! -f $b/$f ]]; then bad "$label: $f missing from $b"; diffs=1; continue; fi
    if cmp -s "$a/$f" "$b/$f"; then
      printf "    %-16s identical (%s bytes)\n" "$f" "$(stat -c%s "$a/$f")"
    else
      printf "    %-16s DIFFERS\n" "$f"; diffs=1
    fi
  done
  [[ $diffs -eq 0 ]] && ok "$label: all artifacts byte-identical" || bad "$label: artifacts differ"
}

# Flags match exactly what the piscem wrapper passes to `piscem_rs::cli::build`
# (main.rs:291-303): canonical, ec-table on, seed 1, partitioned MPHF, dict auto.
# --build-ec-table is what invokes the rewritten eq_classes code; without it this
# whole validation would exercise nothing.
CF=$V/cfish/index_cfish
build_idx() { # binary outdir
  rm -rf "$2"; mkdir -p "$2"
  "$1" build -i "$CF" -o "$2/index" -k 23 -m 19 -t 32 \
       --build-ec-table --canonical --seed 1 --dict auto >"$2/build.log" 2>&1 \
    || { bad "index build failed: $1 -> $2 (see $2/build.log)"; return 1; }
  [[ -f $2/index.tdct ]] || bad "no Tiny artifacts emitted in $2 -- TinyDict not exercised"
  [[ -f $2/index.ectab ]] || bad "no ectab in $2 -- eq_classes not exercised"
}

map_run() { # binary index outdir threads r1 r2 label
  local bin=$1 idx=$2 out=$3 t=$4 r1=$5 r2=$6 label=$7
  rm -rf "$out"; mkdir -p "$out"
  /usr/bin/time -f "%e %U %M" -o "$out/.time" \
    "$bin" map-scrna -i "$idx/index" -g "$G" -o "$out" -t "$t" --dict auto \
      --skipping-strategy permissive --max-ec-card 4096 -1 "$r1" -2 "$r2" \
      >"$out/out.log" 2>"$out/err.log"
  local rc=$?
  if [[ $rc -ne 0 ]]; then bad "$label: map exited $rc (see $out/err.log)"; return 1; fi
  [[ -s $out/map_info.json ]] || { bad "$label: no map_info.json"; return 1; }
  ok "$label completed: $(tr '\n' ' ' < "$out/.time" | awk '{printf "wall=%ss user=%ss rssMB=%d", $1,$2,$3/1024}')"
}

info() { python3 -c "
import json,sys
d=json.load(open('$1/map_info.json'))
print(' '.join(f'{k}={d[k]}' for k in ('num_reads','num_mapped','percent_mapped','num_poisoned','num_refs') if k in d))
"; }

echo "NEW = $($NEW --version 2>&1 | head -1)   [$(git -C /scratch1/rob/alevin-fry-ecosystem/piscem-rs log --oneline -1 | cut -c1-7)]"
echo "OLD = $($OLD --version 2>&1 | head -1)   [7ea5b84, parent of the eq_classes change]"

# ---------------------------------------------------------------- Stage 0
note "Stage 0 CONTROL: same binary twice -- is index build even deterministic?"
build_idx "$OLD" "$V/idx_old_a" && build_idx "$OLD" "$V/idx_old_b"
cmp_idx "$V/idx_old_a" "$V/idx_old_b" "control (OLD vs OLD)"

# ---------------------------------------------------------------- Stage 1
note "Stage 1: index built by NEW vs OLD (isolates the eq_classes rewrite)"
build_idx "$NEW" "$V/idx_new"
cmp_idx "$V/idx_old_a" "$V/idx_new" "index (OLD vs NEW)"
echo "  ectab is the direct product of the changed code:"
ls -l "$V/idx_old_a/index.ectab" "$V/idx_new/index.ectab" 2>/dev/null | sed 's/^/    /'

# ---------------------------------------------------------------- Stage 2
note "Stage 2: 10k-read fixture at -t 1 -- RAD output must be byte-identical"
map_run "$OLD" "$V/idx_old_a" "$V/tiny_old" 1 "$T1" "$T2" "tiny OLD -t1"
map_run "$NEW" "$V/idx_new"   "$V/tiny_new" 1 "$T1" "$T2" "tiny NEW -t1"
echo "  OLD: $(info "$V/tiny_old")"
echo "  NEW: $(info "$V/tiny_new")"
if cmp -s "$V/tiny_old/map.rad" "$V/tiny_new/map.rad"; then
  ok "map.rad byte-identical ($(stat -c%s "$V/tiny_new/map.rad") bytes)"
else
  bad "map.rad DIFFERS at -t 1"
fi
cmp -s "$V/tiny_old/unmapped_bc_count.bin" "$V/tiny_new/unmapped_bc_count.bin" \
  && ok "unmapped_bc_count.bin byte-identical" || bad "unmapped_bc_count.bin differs"

# ---------------------------------------------------------------- Stage 3
if [[ "${1:-all}" == "quick" ]]; then
  note "Stage 3 SKIPPED (quick mode)"
  [[ $FAIL -eq 0 ]] && echo "  stages 0-2 PASSED" || echo "  *** stages 0-2 FAILED ***"
  exit $FAIL
fi
note "Stage 3: FULL flexv2 (281 GB gz, R1+R2) at -t 64, both binaries"
map_run "$OLD" "$V/idx_old_a" "$V/full_old" 64 "$R1" "$R2" "full OLD -t64"
echo "  OLD: $(info "$V/full_old")"
rm -f "$V/full_old/map.rad"
map_run "$NEW" "$V/idx_new" "$V/full_new" 64 "$R1" "$R2" "full NEW -t64"
echo "  NEW: $(info "$V/full_new")"
rm -f "$V/full_new/map.rad"

if [[ -s $V/full_old/map_info.json && -s $V/full_new/map_info.json ]]; then
  if [[ "$(info "$V/full_old")" == "$(info "$V/full_new")" ]]; then
    ok "full-run mapping counts identical"
  else
    bad "full-run mapping counts DIFFER"
  fi
fi
echo "  broker telemetry (NEW):"
grep -iE "decoder|broker|split|converge" "$V/full_new/err.log" 2>/dev/null | head -6 | sed 's/^/    /'

note "RESULT"
[[ $FAIL -eq 0 ]] && echo "  ALL STAGES PASSED" || echo "  *** ONE OR MORE STAGES FAILED ***"
exit $FAIL
