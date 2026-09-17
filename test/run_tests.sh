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

# gzip -cd is spelled the same on macOS and Linux; gzcat/zcat are not.
gz() { gzip -cd "$1"; }

# Builds a .tgls fixture from a tracked .tped: same four leading columns, the
# same number of genotype columns, every value set to $2.  Derived rather than
# committed so there is no new binary fixture to keep in step with the data.
mk_tgls() {
    gz "$EX/chr21.tped.gz" | awk -v v="$2" '{
        printf "%s\t%s\t%s\t%s", $1, $2, $3, $4
        n = (NF - 4) / 2
        for (i = 0; i < n; i++) printf "\t%s", v
        printf "\n"
    }' | gzip > "$1"
}

# As mk_tgls, but puts $4 at locus $2 (1-based) for individual $3 and $5
# everywhere else -- for the guards that fire on a single bad value.
mk_tgls_one() {
    gz "$EX/chr21.tped.gz" | awk -v L="$2" -v S="$3" -v bad="$4" -v good="$5" '{
        printf "%s\t%s\t%s\t%s", $1, $2, $3, $4
        n = (NF - 4) / 2
        for (i = 0; i < n; i++) printf "\t%s", (NR == L && i == S) ? bad : good
        printf "\n"
    }' | gzip > "$1"
}
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
    OBJ=""
    for o in "$ROOT"/src/*.o; do
        case "$o" in *garlic-main.o|*countFeatures*) continue;; esac
        [ -f "$o" ] && OBJ="$OBJ $o"
    done
    if [ -z "$OBJ" ]; then
        echo "  SKIP  no object files found (run make in src/ first)"
        return
    fi
    # zlib is the only external library left; garlic-math.o replaced GSL.
    # shellcheck disable=SC2086
    if $CXX -O1 -std=c++11 -I"$ROOT/include" -I"$ROOT/src" \
            "$ROOT/test/unit_tests.cpp" $OBJ -lz \
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
# The four chr21.* files were force-added for exactly that reason (.gitignore
# has *.gz).  example.GQ.tgls.gz is still untracked, so the tgls_gq case is
# skipped rather than failed when it is absent.
#
# Most cases run on chr21 (8,599 loci) rather than genome-wide (577,489), so
# the suite is fast.  Four deliberately stay genome-wide because nothing
# smaller covers them: unweighted (22 chromosomes, per-chromosome centromere
# handling, genome-wide allele frequencies), chr_subset (needs more than one
# chromosome to select from), multiclass (needs a spread of ROH lengths to
# reach twelve size classes), and tgls_gl.
#
# TWO VACUITY TRAPS, both of which produced silently meaningless cases here
# before being caught:
#   - chr21.tgls.gz is a GL-convention file (-0.0004 everywhere, plus 9 exact
#     zeros).  Read as GQ, BOTH of those values convert to an error rate of 1,
#     so every LOD score is exactly 0: the raw matrix held 371,700 zeros and
#     2,655 NA, two distinct values, and the run called zero ROH.  The earlier
#     c21_tgls_rawlod case pinned that matrix -- non-empty, and still a
#     checksum of nothing.  It is replaced by c21_tgls_gq30, on a generated
#     GQ-30 fixture, whose raw matrix holds 274,055 distinct values; the
#     likelihood_guards stage asserts it equals the --error 0.001 run exactly,
#     which is the equivalence that makes the case meaningful rather than
#     merely non-empty.  chr21.tgls.gz is now rejected under every --gl-type:
#     negative values are invalid as GQ and PL, and its zeros are invalid as GL.
#   - example.GQ.tgls.gz is GQ 30 for every genotype, and GQ 30 converts to an
#     error rate of exactly 0.001, so tgls_gq reproduces the unweighted run
#     byte for byte.  That is expected, not a redundant case.
# The golden stage asserts every .roh.bed case called at least one ROH, so a
# case cannot silently become vacuous again.
#
# Every case is fully specified where it can be -- explicit --lod-cutoff and
# --size-bounds -- so no value is estimated and the output is a deterministic
# function of the input.  The auto-cutoff cases are included on purpose:
# --kde-subsample defaults to 0, which removed the last consumer of random
# numbers on that path, so a regression there would reintroduce seed
# dependence.
CASES="
unweighted|--tped $EX/example.tped.gz --tfam $EX/example.tfam --build hg18 --winsize 60 --error 0.001 --lod-cutoff 2.5 --size-bounds 500000 1000000|roh.bed freq.gz
chr_subset|--tped $EX/example.tped.gz --tfam $EX/example.tfam --build hg18 --winsize 60 --error 0.001 --lod-cutoff 2.5 --size-bounds 500000 1000000 --chr chr21 chr22|roh.bed
multiclass|--tped $EX/example.tped.gz --tfam $EX/example.tfam --build hg18 --winsize 60 --error 0.001 --lod-cutoff 2.5 --size-bounds 1e5 2e5 3e5 5e5 1e6 2e6 3e6 4e6 5e6 6e6 7e6|roh.bed
tgls_gl|--tped $EX/example.tped.gz --tfam $EX/example.tfam --tgls $EX/example.tgls.gz --gl-type GL --build hg18 --winsize 60 --lod-cutoff 2.5 --size-bounds 500000 1000000|roh.bed
c21_unweighted|--tped $EX/chr21.tped.gz --tfam $EX/chr21.tfam.gz --build hg18 --winsize 60 --error 0.001 --lod-cutoff 2.5 --size-bounds 500000 1000000|roh.bed freq.gz
c21_autocutoff|--tped $EX/chr21.tped.gz --tfam $EX/chr21.tfam.gz --build hg18 --winsize 60 --error 0.001|roh.bed
c21_autowinsize|--tped $EX/chr21.tped.gz --tfam $EX/chr21.tfam.gz --build hg18 --auto-winsize --winsize 30 --error 0.001|roh.bed
c21_weighted|--tped $EX/chr21.tped.gz --tfam $EX/chr21.tfam.gz --map $EX/chr21.map.gz --weighted --build hg18 --winsize 60 --error 0.001 --lod-cutoff 2.5 --size-bounds 500000 1000000|roh.bed
c21_cm|--tped $EX/chr21.tped.gz --tfam $EX/chr21.tfam.gz --map $EX/chr21.map.gz --cm --build hg18 --winsize 60 --error 0.001 --lod-cutoff 2.5 --size-bounds 0.5 1.0|roh.bed
c21_phased|--tped $EX/chr21.tped.gz --tfam $EX/chr21.tfam.gz --phased --build hg18 --winsize 60 --error 0.001 --lod-cutoff 2.5 --size-bounds 500000 1000000|roh.bed
c21_froh|--tped $EX/chr21.tped.gz --tfam $EX/chr21.tfam.gz --build hg18 --winsize 60 --error 0.001 --lod-cutoff 2.5 --size-bounds 500000 1000000 --froh|roh.bed froh.tsv
c21_rawlod|--tped $EX/chr21.tped.gz --tfam $EX/chr21.tfam.gz --build hg18 --winsize 60 --error 0.001 --lod-cutoff 2.5 --size-bounds 500000 1000000 --raw-lod|chr21.raw.lod.windows.gz
c21_freqonly|--tped $EX/chr21.tped.gz --tfam $EX/chr21.tfam.gz --build hg18 --error 0.001 --freq-only|freq.gz
c21_tgls_gq30|--tped $EX/chr21.tped.gz --tfam $EX/chr21.tfam.gz --tgls $WORK/gq30.tgls.gz --gl-type GQ --build hg18 --winsize 60 --lod-cutoff 2.5 --size-bounds 500000 1000000 --raw-lod|chr21.raw.lod.windows.gz
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
            # A case that calls no ROH compares two empty files and proves
            # nothing.  This is how two vacuous cases got into this suite.
            case "$suffix" in
                roh.bed)
                    # NOT 'grep -vc ... || echo 0': grep -vc prints 0 AND
                    # exits 1 when nothing matches, so the || fires too and
                    # tracts becomes "0\n0", which makes [ -eq ] a syntax
                    # error and this guard silently never fires.  wc always
                    # exits 0.
                    tracts=$(grep -v '^track' "$f" 2>/dev/null | wc -l | tr -d ' ')
                    if [ "$tracts" -eq 0 ]; then
                        echo "  FAIL  $name: called 0 ROH -- this case is vacuous"
                        echo x >>"$WORK/failures"; continue
                    fi
                    ;;
            esac
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
        $GARLIC --tped "$EX/chr21.tped.gz" --tfam "$EX/chr21.tfam.gz" --build hg18 \
                --winsize 60 --error 0.001 --out "$WORK/det$i" --quiet --force >/dev/null 2>&1
    done
    a=$(sum "$WORK/det1.roh.bed"); b=$(sum "$WORK/det2.roh.bed"); c=$(sum "$WORK/det3.roh.bed")
    if [ "$a" = "$b" ] && [ "$b" = "$c" ]; then
        s1=$(wc -l <"$WORK/det1.params.json" 2>/dev/null | tr -d ' ')
        [ -n "$s1" ] || s1=0
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
# 4b. Genotype-likelihood guards
# ---------------------------------------------------------------------------
# A tgls file states P(genotype CORRECT) for PL and GL, which is NOT what a
# VCF's PL/GL fields carry: VCF normalises them so the CALLED genotype is
# exactly 0.  Feeding those in used to run to completion with an error rate of
# 0 clamped to 1e-16, which makes every heterozygote contribute -16 to the
# window LOD -- 113 ROH instead of 171 on this data, with no other symptom.
# These are the checks that reject such input instead of computing with it.
#
# The fixtures are generated from chr21.tped.gz so the row and column counts
# match; the guards fire on the first offending value, but the file still has
# to be the right shape to reach them.
likelihood_guards() {
    echo "== genotype likelihood guards =="

    # An exact 0 is P(genotype correct) = 1, which no measurement supports and
    # which is exactly what VCF normalisation writes.  Field-specific: a GQ of
    # 0 legitimately means "no confidence".
    mk_tgls_one "$WORK/zero_gl.tgls.gz" 500 7 0 -0.0004
    mk_tgls_one "$WORK/zero_pl.tgls.gz" 500 7 0 0.00435
    expect_exit 2 "GL of exactly 0 is rejected" "$GARLIC" --tped "$EX/chr21.tped.gz" \
        --tfam "$EX/chr21.tfam.gz" --tgls "$WORK/zero_gl.tgls.gz" --gl-type GL \
        --build hg18 --winsize 60 --lod-cutoff 2.5 --size-bounds 500000 1000000 --out "$WORK/g1" --force
    expect_exit 2 "PL of exactly 0 is rejected" "$GARLIC" --tped "$EX/chr21.tped.gz" \
        --tfam "$EX/chr21.tfam.gz" --tgls "$WORK/zero_pl.tgls.gz" --gl-type PL \
        --build hg18 --winsize 60 --lod-cutoff 2.5 --size-bounds 500000 1000000 --out "$WORK/g2" --force

    # No INTEGER PL expresses a realistic error rate under this convention:
    # they map to 0, 0.206, 0.369, 0.499, ... while 0.001 needs PL = 0.00435.
    # So an all-integer PL file came from a VCF however it was produced.
    mk_tgls "$WORK/int_pl.tgls.gz" 3
    expect_exit 2 "all-integer PL file is rejected" "$GARLIC" --tped "$EX/chr21.tped.gz" \
        --tfam "$EX/chr21.tfam.gz" --tgls "$WORK/int_pl.tgls.gz" --gl-type PL \
        --build hg18 --winsize 60 --lod-cutoff 2.5 --size-bounds 500000 1000000 --out "$WORK/g3" --force
    # ... and the same file as GQ is fine: integer GQ is the native scale.
    expect_exit 0 "the same file as GQ is accepted" "$GARLIC" --tped "$EX/chr21.tped.gz" \
        --tfam "$EX/chr21.tfam.gz" --tgls "$WORK/int_pl.tgls.gz" --gl-type GQ \
        --build hg18 --winsize 60 --lod-cutoff 2.5 --size-bounds 500000 1000000 --out "$WORK/g4" --force

    # Values impossible on their own scale.
    mk_tgls_one "$WORK/neg_gq.tgls.gz" 3 0 -1 30
    mk_tgls_one "$WORK/pos_gl.tgls.gz" 3 0 0.5 -0.0004
    mk_tgls_one "$WORK/junk.tgls.gz"   9 3 NA 30
    expect_exit 2 "negative GQ is rejected" "$GARLIC" --tped "$EX/chr21.tped.gz" \
        --tfam "$EX/chr21.tfam.gz" --tgls "$WORK/neg_gq.tgls.gz" --gl-type GQ \
        --build hg18 --winsize 60 --lod-cutoff 2.5 --size-bounds 500000 1000000 --out "$WORK/g5" --force
    expect_exit 2 "positive GL is rejected" "$GARLIC" --tped "$EX/chr21.tped.gz" \
        --tfam "$EX/chr21.tfam.gz" --tgls "$WORK/pos_gl.tgls.gz" --gl-type GL \
        --build hg18 --winsize 60 --lod-cutoff 2.5 --size-bounds 500000 1000000 --out "$WORK/g6" --force
    # An unparseable token used to leave the value at 0 and fail the stream, so
    # the rest of the line silently became zeros.
    expect_exit 2 "an unparseable value is rejected" "$GARLIC" --tped "$EX/chr21.tped.gz" \
        --tfam "$EX/chr21.tfam.gz" --tgls "$WORK/junk.tgls.gz" --gl-type GQ \
        --build hg18 --winsize 60 --lod-cutoff 2.5 --size-bounds 500000 1000000 --out "$WORK/g7" --force

    # A fractional PL file is what the convention actually asks for, and runs.
    mk_tgls "$WORK/frac_pl.tgls.gz" 0.00435
    expect_exit 0 "a fractional PL file is accepted" "$GARLIC" --tped "$EX/chr21.tped.gz" \
        --tfam "$EX/chr21.tfam.gz" --tgls "$WORK/frac_pl.tgls.gz" --gl-type PL \
        --build hg18 --winsize 60 --lod-cutoff 2.5 --size-bounds 500000 1000000 --out "$WORK/g8" --force

    # The shipped chr21.tgls.gz is now invalid under every --gl-type: it holds
    # negative values (invalid as GQ and PL) and exact zeros (invalid as GL).
    # This is the file the replaced c21_tgls_rawlod case used to read as GQ.
    if [ -f "$EX/chr21.tgls.gz" ]; then
        expect_exit 2 "shipped chr21.tgls.gz rejected as GQ" "$GARLIC" --tped "$EX/chr21.tped.gz" \
            --tfam "$EX/chr21.tfam.gz" --tgls "$EX/chr21.tgls.gz" --gl-type GQ \
            --build hg18 --winsize 60 --lod-cutoff 2.5 --size-bounds 500000 1000000 --out "$WORK/g9" --force
        expect_exit 2 "shipped chr21.tgls.gz rejected as GL" "$GARLIC" --tped "$EX/chr21.tped.gz" \
            --tfam "$EX/chr21.tfam.gz" --tgls "$EX/chr21.tgls.gz" --gl-type GL \
            --build hg18 --winsize 60 --lod-cutoff 2.5 --size-bounds 500000 1000000 --out "$WORK/g10" --force
    else
        echo "  SKIP shipped chr21.tgls.gz: input not present"
    fi

    # The equivalence that makes c21_tgls_gq30 a real case rather than a
    # non-empty one: GQ 30 IS an error rate of 0.001, so supplying it per
    # genotype must reproduce --error 0.001 exactly.  If the GQ conversion ever
    # drifts, these two stop matching.
    $GARLIC --tped "$EX/chr21.tped.gz" --tfam "$EX/chr21.tfam.gz" --tgls "$WORK/gq30.tgls.gz" \
        --gl-type GQ --build hg18 --winsize 60 --lod-cutoff 2.5 --size-bounds 500000 1000000 \
        --raw-lod --out "$WORK/eqgq" --quiet --force >/dev/null 2>&1
    $GARLIC --tped "$EX/chr21.tped.gz" --tfam "$EX/chr21.tfam.gz" --error 0.001 \
        --build hg18 --winsize 60 --lod-cutoff 2.5 --size-bounds 500000 1000000 \
        --raw-lod --out "$WORK/eqfl" --quiet --force >/dev/null 2>&1
    if [ "$(sum "$WORK/eqgq.chr21.raw.lod.windows.gz")" = "$(sum "$WORK/eqfl.chr21.raw.lod.windows.gz")" ]; then
        ok
    else
        bad "GQ 30 per genotype did not reproduce --error 0.001"
    fi
    # ... and that case must not be vacuous the way the one it replaced was.
    n=$(gz "$WORK/eqgq.chr21.raw.lod.windows.gz" | tr -s ' \t' '\n' | sort -u | grep -c .)
    if [ "$n" -gt 1000 ]; then ok; else bad "c21_tgls_gq30 raw LOD has only $n distinct values"; fi
}

# ---------------------------------------------------------------------------
# 4c. Individual metadata
# ---------------------------------------------------------------------------
# The duplicate-ID error and the pooled-population warning were moved out of
# scanIndData3, which only --tfam calls, into checkIndData over the assembled
# IndData, so they apply to every input path.  These cases pin the TFAM path's
# behaviour across that move: same messages, same exit codes.  The pooled-
# population warning in particular is what catches a TFAM with population 0
# for every sample, which is what the deprecated vcf2tped.pl writes.
ind_metadata() {
    echo "== individual metadata =="

    # A duplicate individual ID.
    gz "$EX/chr21.tfam.gz" | awk 'NR==5{$2="HGDP00607"} {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > "$WORK/dup.tfam"
    expect_exit 2 "a duplicate individual ID is rejected" "$GARLIC" --tped "$EX/chr21.tped.gz" \
        --tfam "$WORK/dup.tfam" --build hg18 --winsize 60 --error 0.001 --lod-cutoff 2.5 \
        --size-bounds 500000 1000000 --out "$WORK/i1" --force
    if "$GARLIC" --tped "$EX/chr21.tped.gz" --tfam "$WORK/dup.tfam" --build hg18 --winsize 60 \
            --error 0.001 --lod-cutoff 2.5 --size-bounds 500000 1000000 --out "$WORK/i1" --force \
            2>&1 | grep -q "Found duplicate individual ID"; then ok
    else bad "the duplicate-ID error message is missing"; fi

    # Two populations in one file: allowed, but it must say so, because
    # frequencies are pooled across all of them.
    gz "$EX/chr21.tfam.gz" | awk 'NR>20{$1="OTHERPOP"} {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > "$WORK/twopop.tfam"
    expect_exit 0 "two populations in one TFAM still runs" "$GARLIC" --tped "$EX/chr21.tped.gz" \
        --tfam "$WORK/twopop.tfam" --build hg18 --winsize 60 --error 0.001 --lod-cutoff 2.5 \
        --size-bounds 500000 1000000 --out "$WORK/i2" --force
    if "$GARLIC" --tped "$EX/chr21.tped.gz" --tfam "$WORK/twopop.tfam" --build hg18 --winsize 60 \
            --error 0.001 --lod-cutoff 2.5 --size-bounds 500000 1000000 --out "$WORK/i2" --force \
            2>&1 | grep -q "Found multiple population IDs"; then ok
    else bad "the pooled-population warning is missing"; fi

    # A single-population TFAM must NOT warn.
    if "$GARLIC" --tped "$EX/chr21.tped.gz" --tfam "$EX/chr21.tfam.gz" --build hg18 --winsize 60 \
            --error 0.001 --lod-cutoff 2.5 --size-bounds 500000 1000000 --out "$WORK/i3" --force \
            2>&1 | grep -q "Found multiple population IDs"; then
        bad "a single-population TFAM produced the pooled-population warning"
    else ok; fi
    # The sex-chromosome warning must never state a number of affected
    # individuals the metadata cannot support.  Sex is optional in a TFAM, so
    # all three coverage cases are reachable.  The partial case used to print
    # only the male count: 3 of 45 coded male with 42 unknown printed
    # "individuals coded male in the TFAM: 3", and 3 reads as the answer.
    gz "$EX/chr21.tped.gz" | sed 's/^21/chrX/' | gzip > "$WORK/chrX.tped.gz"
    gz "$EX/chr21.tfam.gz" | awk '{print $1"\t"$2}' > "$WORK/nosex.tfam"
    gz "$EX/chr21.tfam.gz" | awk 'NR<=3{print $1"\t"$2"\t0\t0\t1\t0"} NR>3{print $1"\t"$2}' > "$WORK/partial.tfam"
    sexwarn() {   # $1 = tfam, $2 = pattern that must appear
        if "$GARLIC" --tped "$WORK/chrX.tped.gz" --tfam "$1" --build hg18 --winsize 60 \
                --error 0.001 --lod-cutoff 2.5 --size-bounds 500000 1000000 \
                --out "$WORK/sx" --force 2>&1 | grep -q "$2"; then ok
        else bad "sex warning on $(basename "$1"): expected /$2/"; fi
    }
    sexwarn "$WORK/nosex.tfam"  "sex is not recorded for any of the 45"
    sexwarn "$WORK/partial.tfam" "3 male, 0 female, 42 unknown"
    sexwarn "$WORK/partial.tfam" "may be as high as 45"
    sexwarn "$EX/chr21.tfam.gz" "26 male, 19 female"
    sexwarn "$EX/chr21.tfam.gz" "affected individuals: 26"
    # ... and a partial-coverage run must NOT print a bare affected count, which
    # is the misleading form it used to print.
    if "$GARLIC" --tped "$WORK/chrX.tped.gz" --tfam "$WORK/partial.tfam" --build hg18 --winsize 60 \
            --error 0.001 --lod-cutoff 2.5 --size-bounds 500000 1000000 --out "$WORK/sx" --force 2>&1 \
            | grep -q "affected individuals:"; then
        bad "a partial-coverage run stated a bare affected count"
    else ok; fi
}

# ---------------------------------------------------------------------------
# 5. Round trip through --load-params
# ---------------------------------------------------------------------------
params_roundtrip() {
    echo "== params round trip =="
    $GARLIC --tped "$EX/chr21.tped.gz" --tfam "$EX/chr21.tfam.gz" --map "$EX/chr21.map.gz" \
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
    f=$WORK/c21_unweighted.roh.bed
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

mk_tgls "$WORK/gq30.tgls.gz" 30

unit_tests
golden
determinism
likelihood_guards
ind_metadata
exit_codes
params_roundtrip
bed_format

echo
echo "$pass checks passed, $fail failure(s)"
if [ "$fail" -eq 0 ]; then rm -rf "$WORK"; echo "OK"; exit 0; fi
echo "outputs kept in $WORK for inspection"
exit 1
