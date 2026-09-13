#!/bin/sh
# garlic test suite.  POSIX sh, no test framework, no network: it needs only a
# built garlic, a C++ compiler for the unit tests, and the bundled example data.
#
#   sh test/run_tests.sh                 # from the repository root
#   GARLIC=/path/to/garlic sh test/run_tests.sh
#   sh test/run_tests.sh --bless         # rewrite the golden checksums
#
# Exits 0 if everything passes, 1 otherwise.

set -u

ROOT=$(cd "$(dirname "$0")/.." && pwd)
GARLIC=${GARLIC:-$ROOT/src/garlic}
EX=$ROOT/example
GOLDEN=$ROOT/test/golden
WORK=${TMPDIR:-/tmp}/garlic-tests.$$
BLESS=no
[ "${1:-}" = "--bless" ] && BLESS=yes

pass=0; fail=0

# md5 is named differently on macOS and Linux.
if command -v md5sum >/dev/null 2>&1; then
    sum() { md5sum "$1" | cut -d' ' -f1; }
elif command -v md5 >/dev/null 2>&1; then
    sum() { md5 -q "$1"; }
else
    echo "ERROR: neither md5sum nor md5 found"; exit 1
fi
# Compare gzipped output by content, not by bytes: gzip embeds an mtime.
sumgz() { gzip -cd "$1" | { if command -v md5sum >/dev/null 2>&1; then md5sum; else md5 -q; fi; } | cut -d' ' -f1; }

ok()   { pass=$((pass+1)); }
bad()  { fail=$((fail+1)); echo "  FAIL  $1"; }

# ---------------------------------------------------------------------------
# 1. Unit tests
# ---------------------------------------------------------------------------
unit_tests() {
    echo "== unit tests =="
    CXX=${CXX:-c++}
    LIBDIR=""
    for d in macos-arm macos linux win32; do
        [ -f "$ROOT/lib/$d/libgsl.a" ] && { LIBDIR=$ROOT/lib/$d; break; }
    done
    OBJ=""
    for o in "$ROOT"/src/*.o; do
        case "$o" in *garlic-main.o|*countFeatures*) continue;; esac
        [ -f "$o" ] && OBJ="$OBJ $o"
    done
    if [ -z "$OBJ" ] || [ -z "$LIBDIR" ]; then
        echo "  SKIP  no object files or GSL archive found (run make in src/ first)"
        return
    fi
    # shellcheck disable=SC2086
    if $CXX -O1 -std=c++11 -I"$ROOT/include" -I"$ROOT/src" \
            "$ROOT/test/unit_tests.cpp" $OBJ -lz "$LIBDIR/libgsl.a" "$LIBDIR/libgslcblas.a" \
            -o "$WORK/unit_tests" 2>"$WORK/unit_build.log"; then
        if "$WORK/unit_tests"; then ok; else bad "unit tests reported failures"; fi
    else
        bad "unit tests did not compile (see $WORK/unit_build.log)"
    fi
}

# ---------------------------------------------------------------------------
# 2. Golden-output regression
#
# Every case below is fully specified -- explicit --lod-cutoff and
# --size-bounds where they exist, so no value is estimated and the output is a
# deterministic function of the input.  The auto-cutoff cases are included on
# purpose: --kde-subsample now defaults to 0, which removed the last consumer
# of random numbers on that path, so those are deterministic too and a
# regression there would reintroduce seed dependence.
# ---------------------------------------------------------------------------
# Every input named here is TRACKED in git, so a fresh clone can run this.
# The chr21.* files and example.GQ.tgls.gz exist in the author's working copy
# but are not in the repository (.gitignore has *.gz and only the example.*
# files were force-added), so cases needing them are listed as optional below
# and skipped rather than failed when absent.
#
# Every case is fully specified where it can be -- explicit --lod-cutoff and
# --size-bounds -- so no value is estimated and the output is a deterministic
# function of the input.  The auto-cutoff cases are included on purpose:
# --kde-subsample defaults to 0, which removed the last consumer of random
# numbers on that path, so a regression there would reintroduce seed
# dependence.
CASES="
unweighted|--tped $EX/example.tped.gz --tfam $EX/example.tfam --build hg18 --winsize 60 --error 0.001 --lod-cutoff 2.5 --size-bounds 500000 1000000|roh.bed freq.gz
autocutoff|--tped $EX/example.tped.gz --tfam $EX/example.tfam --build hg18 --winsize 60 --error 0.001|roh.bed
autowinsize|--tped $EX/example.tped.gz --tfam $EX/example.tfam --build hg18 --auto-winsize --winsize 30 --error 0.001|roh.bed
weighted|--tped $EX/example.tped.gz --tfam $EX/example.tfam --map $EX/example.map.gz --weighted --build hg18 --winsize 60 --error 0.001 --lod-cutoff 2.5 --size-bounds 500000 1000000|roh.bed
cm|--tped $EX/example.tped.gz --tfam $EX/example.tfam --map $EX/example.map.gz --cm --build hg18 --winsize 60 --error 0.001 --lod-cutoff 2.5 --size-bounds 0.5 1.0|roh.bed
phased|--tped $EX/example.tped.gz --tfam $EX/example.tfam --phased --build hg18 --winsize 60 --error 0.001 --lod-cutoff 2.5 --size-bounds 500000 1000000|roh.bed
tgls_gl|--tped $EX/example.tped.gz --tfam $EX/example.tfam --tgls $EX/example.tgls.gz --gl-type GL --build hg18 --winsize 60 --lod-cutoff 2.5 --size-bounds 500000 1000000|roh.bed
rawlod|--tped $EX/example.tped.gz --tfam $EX/example.tfam --build hg18 --winsize 60 --error 0.001 --lod-cutoff 2.5 --size-bounds 500000 1000000 --raw-lod|chr22.raw.lod.windows.gz
freqonly|--tped $EX/example.tped.gz --tfam $EX/example.tfam --build hg18 --error 0.001 --freq-only|freq.gz
froh|--tped $EX/example.tped.gz --tfam $EX/example.tfam --build hg18 --winsize 60 --error 0.001 --lod-cutoff 2.5 --size-bounds 500000 1000000 --froh|roh.bed froh.tsv
chr_subset|--tped $EX/example.tped.gz --tfam $EX/example.tfam --build hg18 --winsize 60 --error 0.001 --lod-cutoff 2.5 --size-bounds 500000 1000000 --chr chr21 chr22|roh.bed
multiclass|--tped $EX/example.tped.gz --tfam $EX/example.tfam --build hg18 --winsize 60 --error 0.001 --lod-cutoff 2.5 --size-bounds 1e5 2e5 3e5 5e5 1e6 2e6 3e6 4e6 5e6 6e6 7e6|roh.bed
tgls_gq|--tped $EX/example.tped.gz --tfam $EX/example.tfam --tgls $EX/example.GQ.tgls.gz --gl-type GQ --build hg18 --winsize 60 --lod-cutoff 2.5 --size-bounds 500000 1000000|roh.bed
"

golden() {
    echo "== golden output =="
    echo "$CASES" | while IFS='|' read -r name args outs; do
        [ -z "$name" ] && continue
        out=$WORK/$name
        # Skip rather than fail when an input is not in the repository.
        missing=
        for tok in $args; do
            case "$tok" in
                */*.gz|*/*.tfam|*/*.tped) [ -f "$tok" ] || missing="$missing $tok";;
            esac
        done
        if [ -n "$missing" ]; then
            echo "  SKIP  $name: input not present:$missing"
            continue
        fi
        # shellcheck disable=SC2086
        $GARLIC $args --out "$out" --quiet --force >"$WORK/$name.stdout" 2>"$WORK/$name.stderr"
        rc=$?
        if [ "$rc" -ne 0 ]; then
            echo "  FAIL  $name: garlic exited $rc (see $WORK/$name.stderr)"
            echo x >>"$WORK/failures"; continue
        fi
        for suffix in $outs; do
            f=$out.$suffix
            if [ ! -f "$f" ]; then
                echo "  FAIL  $name: expected output $name.$suffix was not created"
                echo x >>"$WORK/failures"; continue
            fi
            case "$suffix" in *.gz) got=$(sumgz "$f");; *) got=$(sum "$f");; esac
            g=$GOLDEN/$name.$suffix.md5
            if [ "$BLESS" = yes ]; then
                echo "$got" >"$g"
                continue
            fi
            if [ ! -f "$g" ]; then
                echo "  FAIL  $name.$suffix: no golden checksum (run with --bless)"
                echo x >>"$WORK/failures"; continue
            fi
            want=$(cat "$g")
            if [ "$got" != "$want" ]; then
                echo "  FAIL  $name.$suffix: $got != $want"
                echo "        output kept at $f"
                echo x >>"$WORK/failures"
            fi
        done
    done
    if [ "$BLESS" = yes ]; then echo "  blessed $(ls "$GOLDEN" | wc -l | tr -d ' ') checksums"; ok; return; fi
    if [ -s "$WORK/failures" ]; then fail=$((fail + $(wc -l <"$WORK/failures"))); : >"$WORK/failures"; else ok; fi
}

# ---------------------------------------------------------------------------
# 3. Determinism: the auto-cutoff path must not depend on the drawn seed
# ---------------------------------------------------------------------------
determinism() {
    echo "== determinism =="
    for i in 1 2 3; do
        $GARLIC --tped "$EX/example.tped.gz" --tfam "$EX/example.tfam" --build hg18 \
                --winsize 60 --error 0.001 --out "$WORK/det$i" --quiet --force >/dev/null 2>&1
    done
    a=$(sum "$WORK/det1.roh.bed"); b=$(sum "$WORK/det2.roh.bed"); c=$(sum "$WORK/det3.roh.bed")
    if [ "$a" = "$b" ] && [ "$b" = "$c" ]; then
        s1=$(grep -c . "$WORK/det1.params.json" 2>/dev/null || echo 0)
        [ "$s1" -gt 0 ] || bad "no params.json written"
        ok
    else
        bad "three auto-cutoff runs disagree ($a $b $c) -- the KDE path is seed dependent again"
    fi
}

# ---------------------------------------------------------------------------
# 4. Exit-code contract: 0 success, 1 usage error, 2 runtime error
# ---------------------------------------------------------------------------
expect_exit() {
    want=$1; what=$2; shift 2
    "$@" >/dev/null 2>&1
    got=$?
    if [ "$got" -eq "$want" ]; then ok; else bad "$what: exit $got, expected $want"; fi
}

exit_codes() {
    echo "== exit codes =="
    expect_exit 0 "--version" "$GARLIC" --version
    expect_exit 0 "--help" "$GARLIC" --help
    expect_exit 0 "-h" "$GARLIC" -h
    expect_exit 1 "no arguments" "$GARLIC"
    expect_exit 1 "unrecognised flag" "$GARLIC" --not-a-real-flag
    expect_exit 1 "missing tped" "$GARLIC" --tped /nonexistent.tped.gz --tfam "$EX/example.tfam" --build hg18 --winsize 60 --out "$WORK/e1"
    expect_exit 1 "no build and no centromere file" "$GARLIC" --tped "$EX/example.tped.gz" --tfam "$EX/example.tfam" --winsize 60 --out "$WORK/e2"
    expect_exit 1 "winsize below the minimum" "$GARLIC" --tped "$EX/example.tped.gz" --tfam "$EX/example.tfam" --build hg18 --winsize 1 --out "$WORK/e3"
    expect_exit 1 "--chr naming a chromosome not in the data" "$GARLIC" --tped "$EX/example.tped.gz" --tfam "$EX/example.tfam" --build hg18 --winsize 60 --lod-cutoff 2.5 --chr chr99 --out "$WORK/e4"
    expect_exit 1 "--load-params on a missing file" "$GARLIC" --load-params /nonexistent.json --out "$WORK/e5"
    expect_exit 1 "--weighted without a map" "$GARLIC" --tped "$EX/example.tped.gz" --tfam "$EX/example.tfam" --weighted --build hg18 --winsize 60 --out "$WORK/e6"
    expect_exit 2 "cutoff selection fails (--mode-smooth-span too wide)" "$GARLIC" --tped "$EX/example.tped.gz" --tfam "$EX/example.tfam" --build hg18 --winsize 60 --error 0.001 --mode-smooth-span 400 --out "$WORK/e7"

    # A rejected command line must not leave output files behind.
    rm -f "$WORK/noclobber".*
    $GARLIC --tped "$EX/example.tped.gz" --not-a-real-flag --out "$WORK/noclobber" >/dev/null 2>&1
    if ls "$WORK/noclobber".* >/dev/null 2>&1; then
        bad "a rejected command line created output files: $(ls "$WORK/noclobber".*)"
    else ok; fi

    # --force is required to overwrite existing calls.
    $GARLIC --tped "$EX/example.tped.gz" --tfam "$EX/example.tfam" --build hg18 --winsize 60 \
            --error 0.001 --lod-cutoff 2.5 --size-bounds 500000 1000000 --out "$WORK/clob" --quiet >/dev/null 2>&1
    expect_exit 1 "second run without --force refuses to clobber" "$GARLIC" --tped "$EX/example.tped.gz" --tfam "$EX/example.tfam" --build hg18 --winsize 60 --error 0.001 --lod-cutoff 2.5 --size-bounds 500000 1000000 --out "$WORK/clob" --quiet

    # A successful run must not create an empty .error file.
    if [ -f "$WORK/clob.error" ]; then bad "a successful run created a .error file"; else ok; fi
}

# ---------------------------------------------------------------------------
# 5. Round trip through --load-params
# ---------------------------------------------------------------------------
params_roundtrip() {
    echo "== params round trip =="
    $GARLIC --tped "$EX/example.tped.gz" --tfam "$EX/example.tfam" --map "$EX/example.map.gz" \
            --weighted --build hg18 --winsize 60 --error 0.001 --froh \
            --out "$WORK/rt1" --quiet --force >/dev/null 2>&1
    $GARLIC --load-params "$WORK/rt1.params.json" --out "$WORK/rt2" --quiet --force >/dev/null 2>&1
    if [ "$(sum "$WORK/rt1.roh.bed")" = "$(sum "$WORK/rt2.roh.bed")" ] &&
       [ "$(sum "$WORK/rt1.froh.tsv")" = "$(sum "$WORK/rt2.froh.tsv")" ]; then ok
    else bad "--load-params did not reproduce the run it recorded"; fi

    # An explicit flag must win over the recorded value.
    $GARLIC --load-params "$WORK/rt1.params.json" --lod-cutoff 3.0 --out "$WORK/rt3" --quiet --force >/dev/null 2>&1
    if grep -q '"lod_cutoff": 3' "$WORK/rt3.params.json" 2>/dev/null; then ok
    else bad "a command line flag did not override the loaded params file"; fi
}

# ---------------------------------------------------------------------------
# 6. BED conformance
# ---------------------------------------------------------------------------
bed_format() {
    echo "== bed format =="
    f=$WORK/unweighted.roh.bed
    [ -f "$f" ] || { bad "no bed file to check (golden stage did not run)"; return; }
    bad_rows=$(grep -v '^track' "$f" | awk '$3-$2 != $5 || $2 < 0 {n++} END {print n+0}')
    if [ "$bad_rows" -eq 0 ]; then ok
    else bad "$bad_rows rows where chromEnd-chromStart != length, or chromStart < 0"; fi
}

# ---------------------------------------------------------------------------

mkdir -p "$WORK" "$GOLDEN"; : >"$WORK/failures"
echo "garlic test suite"
echo "  binary: $GARLIC"
[ -x "$GARLIC" ] || { echo "ERROR: $GARLIC is not executable. Run make in src/ first."; exit 1; }
echo "  version: $("$GARLIC" --version 2>&1 | head -1)"
echo

for f in "$EX/example.tped.gz" "$EX/example.tfam" "$EX/example.map.gz"; do
    [ -f "$f" ] || { echo "ERROR: required example data missing: $f"; exit 1; }
done

unit_tests
golden
determinism
exit_codes
params_roundtrip
bed_format

echo
echo "$pass checks passed, $fail failure(s)"
if [ "$fail" -eq 0 ]; then rm -rf "$WORK"; echo "OK"; exit 0; fi
echo "outputs kept in $WORK for inspection"
exit 1
