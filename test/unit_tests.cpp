// Unit tests for garlic's computational helpers.
//
// Deliberately free of any test framework: the whole point is that this runs
// in CI on a bare toolchain with nothing to install.  Link against the same
// object files the program uses (see test/Makefile).
//
// Every case here is either a hand-computed value or an invariant that must
// hold by construction.  Cases marked REGRESSION reproduce a bug that was
// actually shipped.  Two of them were checked against the pre-fix objects and
// do fail there -- get_arg_max returns -1 instead of 3, and getMapInfo returns
// 8.5 instead of 1.5 -- which is the only evidence that a passing test here
// means anything.  The ones marked as previously inline cannot be linked
// against the old code, so they pin behaviour rather than prove discrimination.

#include "garlic-roh.h"
#include "garlic-data.h"
#include "garlic-kde.h"
#include "garlic-centromeres.h"
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include <set>

using namespace std;

static int failures = 0;
static int checks   = 0;

static void ck(bool ok, const char *what)
{
    checks++;
    if (!ok) { failures++; printf("  FAIL  %s\n", what); }
}

static void ckd(double got, double want, double tol, const char *what)
{
    checks++;
    if (!(fabs(got - want) <= tol))
    {
        failures++;
        printf("  FAIL  %s\n          got %.17g  want %.17g  (tol %g)\n", what, got, want, tol);
    }
}

static void cks(const string &got, const string &want, const char *what)
{
    checks++;
    if (got != want)
    {
        failures++;
        printf("  FAIL  %s\n          got '%s'  want '%s'\n", what, got.c_str(), want.c_str());
    }
}

// ---------------------------------------------------------------- lod() ----
// lod(g, p, e) = log10( P(g | autozygous) / P(g | not autozygous) ) with a
// per-genotype error rate e.  The three branches are hand-computable.
static void test_lod()
{
    const double e = 0.001;

    // A heterozygote under autozygosity can only arise from error, and the
    // 2p(1-p) term cancels, so the score is exactly log10(e) at every
    // frequency.  This is the sharpest invariant the function has.
    for (double p = 0.05; p < 1.0; p += 0.05)
        ckd(lod(1, p, e), log10(e), 1e-12, "lod(het) == log10(error), independent of freq");

    // A monomorphic site carries no information either way.
    ckd(lod(0, 0.0, e), 0.0, 0.0, "lod(hom-ref, freq=0) == 0 exactly");
    ckd(lod(2, 1.0, e), 0.0, 0.0, "lod(hom-alt, freq=1) == 0 exactly");
    ckd(lod(1, 0.0, e), 0.0, 0.0, "lod(het, freq=0) == 0 exactly");

    // Missing genotypes must be neutral, not scored.
    ckd(lod(MISSING, 0.3, e), 0.0, 0.0, "lod(MISSING) == 0 exactly");
    ckd(lod(-1, 0.3, e), 0.0, 0.0, "lod(unknown code) == 0 exactly");

    // Hand-computed: p = 0.5, e = 0.001, g = 0.
    //   nonAutozygous = (1-p)^2                   = 0.25
    //   autozygous    = (1-e)(1-p) + e(1-p)^2     = 0.999*0.5 + 0.001*0.25 = 0.49975
    //   lod           = log10(0.49975 / 0.25)     = log10(1.999)
    ckd(lod(0, 0.5, e), log10(1.999), 1e-15, "lod(hom-ref, p=0.5, e=1e-3) == log10(1.999)");
    ckd(lod(2, 0.5, e), log10(1.999), 1e-15, "lod(hom-alt, p=0.5, e=1e-3) == log10(1.999)");

    // The two homozygote branches are mirror images under p -> 1-p.
    for (double p = 0.05; p < 1.0; p += 0.05)
        ckd(lod(0, p, e), lod(2, 1.0 - p, e), 1e-12, "lod(hom-ref, p) == lod(hom-alt, 1-p)");

    // A rarer homozygote is stronger evidence of autozygosity, so the score
    // must increase monotonically as that genotype becomes less common.
    double prev = -1e300;
    bool increasing = true;
    for (double p = 0.95; p > 0.02; p -= 0.05)
    {
        double v = lod(2, p, e);
        if (v < prev) increasing = false;
        prev = v;
    }
    ck(increasing, "lod(hom-alt) increases as that homozygote gets rarer");

    // With no error, a heterozygote is impossible under autozygosity.
    ck(lod(1, 0.5, 0.0) == -INFINITY || lod(1, 0.5, 0.0) < -300,
       "lod(het, error=0) is -inf (an impossible observation)");
}

// -------------------------------------------------------- interpolate() ----
static void test_interpolate()
{
    // y = 2x + 1 through (1,3) and (5,11)
    ckd(interpolate(1, 3, 5, 11, 1), 3.0,  1e-12, "interpolate returns y0 at x0");
    ckd(interpolate(1, 3, 5, 11, 5), 11.0, 1e-12, "interpolate returns y1 at x1");
    ckd(interpolate(1, 3, 5, 11, 3), 7.0,  1e-12, "interpolate at the midpoint");
    // Extrapolation is linear too (the caller is responsible for bracketing).
    ckd(interpolate(1, 3, 5, 11, 7), 15.0, 1e-12, "interpolate extrapolates linearly");
    // A flat interval must return the constant.
    ckd(interpolate(0, 4, 10, 4, 6), 4.0, 1e-12, "interpolate on a flat interval");
}

// ------------------------------------------------------------- inGap() -----
static void test_inGap()
{
    // gap spans [1000, 2000]
    ck(!inGap(10,   100,  1000, 2000), "query entirely before the gap");
    ck(!inGap(3000, 4000, 1000, 2000), "query entirely after the gap");
    ck( inGap(1500, 3000, 1000, 2000), "query starts inside the gap");
    ck( inGap(10,   1500, 1000, 2000), "query ends inside the gap");
    ck( inGap(500,  3000, 1000, 2000), "gap entirely inside the query");
    ck( inGap(1200, 1800, 1000, 2000), "query entirely inside the gap");
    ck( inGap(2000, 3000, 1000, 2000), "query starts exactly at the gap end");
    ck( inGap(10,   1000, 1000, 2000), "query ends exactly at the gap start");
    ck(!inGap(10,   999,  1000, 2000), "query ends one base before the gap");
    ck(!inGap(2001, 3000, 1000, 2000), "query starts one base after the gap");
}

// -------------------------------------- automatic window / overlap fits ----
static void test_auto_fits()
{
    // size = round(8.3235*ln(d) + 138.0521), floored at 10.
    int want = int(8.3235 * log(100.0) + 138.0521 + 0.5);
    ck(selectWinsizeWeighted(100.0) == want, "selectWinsizeWeighted matches its formula");
    // The fit goes negative at very low density; the floor must hold.
    ck(selectWinsizeWeighted(1e-12) == 10, "selectWinsizeWeighted floors at 10");
    ck(selectWinsizeWeighted(1e-30) == 10, "selectWinsizeWeighted floors at 10 (extreme)");

    // frac = (6.375*ln(d) + 63.888)/100, clamped to (0,1]; <=0 becomes 1/winsize.
    double f = selectOverlapFrac(100.0, 60);
    ckd(f, (6.375 * log(100.0) + 63.888) / 100.0, 1e-12, "selectOverlapFrac matches its formula");
    ck(selectOverlapFrac(1e9, 60) == 1.0, "selectOverlapFrac clamps at 1");
    ckd(selectOverlapFrac(1e-12, 60), 1.0 / 60.0, 1e-12,
        "selectOverlapFrac falls back to 1/winsize when the fit goes non-positive");
}

// -------------------------------------------- size class labels/colours ----
static void test_class_labels()
{
    cks(sizeClassLabel(0),  "A",  "sizeClassLabel(0)");
    cks(sizeClassLabel(1),  "B",  "sizeClassLabel(1)");
    cks(sizeClassLabel(25), "Z",  "sizeClassLabel(25)");
    // REGRESSION (C5.19): the old code did char('A'+k) inline at the call
    // site, which runs past 'Z' into punctuation for more than 26 classes.
    // Being inline, it cannot be linked against, so unlike the B5 and B7
    // cases below this test cannot be run against the old code to prove it
    // discriminates -- it pins the fixed behaviour instead.
    cks(sizeClassLabel(26), "AA", "sizeClassLabel(26) is AA, not '['");
    cks(sizeClassLabel(27), "AB", "sizeClassLabel(27)");
    cks(sizeClassLabel(51), "AZ", "sizeClassLabel(51)");
    cks(sizeClassLabel(52), "BA", "sizeClassLabel(52)");

    // REGRESSION (C5.19): only nine colours existed and the index was clamped
    // at 8, so every class past the ninth was drawn the same grey.  Also
    // previously inline; see the note above.
    vector<string> c = makeClassColors(30);
    ck(c.size() == 30, "makeClassColors returns one colour per class");
    set<string> uniq(c.begin(), c.end());
    ck(uniq.size() == 30, "makeClassColors returns 30 DISTINCT colours");
    cks(c[0], "228,26,28", "makeClassColors keeps the original first colour");
    cks(c[8], "153,153,153", "makeClassColors keeps the original ninth colour");

    // Every colour must be three bytes a genome browser will accept.
    bool wellFormed = true;
    for (size_t i = 0; i < c.size(); i++)
    {
        int r = -1, g = -1, b = -1;
        if (sscanf(c[i].c_str(), "%d,%d,%d", &r, &g, &b) != 3) { wellFormed = false; break; }
        if (r < 0 || r > 255 || g < 0 || g > 255 || b < 0 || b > 255) { wellFormed = false; break; }
    }
    ck(wellFormed, "every generated colour is three values in 0..255");
}

// ------------------------------------------------------ checkChrName() -----
static void test_chr_names()
{
    cks(checkChrName("1"),    "chr1",  "checkChrName prefixes a bare number");
    cks(checkChrName("chr1"), "chr1",  "checkChrName leaves chr1 alone");
    cks(checkChrName("X"),    "chrX",  "checkChrName prefixes X");
    cks(checkChrName("22"),   "chr22", "checkChrName prefixes 22");
    // Documenting actual behaviour, not endorsing it: the test is a 'starts
    // with c' check, so names that already begin with c are passed through
    // and an upper-case Chr1 is NOT recognised.  Both sides of a comparison
    // go through this function, so garlic stays self-consistent; a map using
    // 'Chr1' against a TPED using 'chr1' would not match.
    cks(checkChrName("ctg7"), "ctg7", "checkChrName passes through names starting with c");
    cks(checkChrName("Chr1"), "chrChr1", "checkChrName is case sensitive (documented wart)");
}

// ------------------------------------------------------- glToError() ------
// The GQ/GL/PL to per-genotype-error conversion.  Tested here rather than
// end to end because both TGLS files bundled with garlic are constant -- GQ
// is 30 for every genotype and GL is -0.0004 for every genotype -- so no
// golden-output case can exercise per-genotype variation in this conversion.
static void test_glToError()
{
    // GQ is Phred scaled: error = 10^(-GQ/10).
    ckd(glToError(30.0, "GQ"), 1e-3, 1e-15, "glToError(GQ 30) == 0.001");
    ckd(glToError(20.0, "GQ"), 1e-2, 1e-15, "glToError(GQ 20) == 0.01");
    ckd(glToError(10.0, "GQ"), 1e-1, 1e-15, "glToError(GQ 10) == 0.1");
    ckd(glToError(0.0,  "GQ"), 1.0,  1e-15, "glToError(GQ 0) == 1 (no confidence)");
    // The exponent is floored at -10, so quality beyond 100 stops helping.
    ckd(glToError(100.0, "GQ"), 1e-10, 1e-20, "glToError(GQ 100) == 1e-10");
    ckd(glToError(500.0, "GQ"), 1e-10, 1e-20, "glToError(GQ 500) is floored at 1e-10");
    // This is why the bundled GQ file reproduces --error 0.001 exactly.
    ckd(glToError(30.0, "GQ"), 0.001, 1e-15,
        "GQ 30 is exactly equivalent to --error 0.001");

    // GL is a log10 likelihood: error = 1 - 10^GL.
    ckd(glToError(-0.0004, "GL"), 1.0 - pow(10, -0.0004), 1e-15,
        "glToError(GL -0.0004) matches 1 - 10^GL");
    ckd(glToError(-1.0, "GL"), 0.9, 1e-15, "glToError(GL -1) == 0.9");
    ckd(glToError(-3.0, "GL"), 0.999, 1e-15, "glToError(GL -3) == 0.999");
    // A likelihood of 1 means no error, which the clamp turns into a tiny
    // positive value rather than zero (a zero error rate makes a heterozygote
    // impossible and the LOD score -inf).
    ck(glToError(0.0, "GL") > 0.0, "glToError(GL 0) is clamped strictly above zero");
    ckd(glToError(0.0, "GL"), 1e-16, 1e-17, "glToError(GL 0) == 1e-16");
    ckd(glToError(-50.0, "GL"), 1.0 - pow(10, -10), 1e-15, "glToError(GL) floors its exponent at -10");

    // PL is Phred-scaled likelihood: error = 1 - 10^(-PL/10).
    ckd(glToError(10.0, "PL"), 0.9,   1e-15, "glToError(PL 10) == 0.9");
    ckd(glToError(30.0, "PL"), 0.999, 1e-15, "glToError(PL 30) == 0.999");
    ck(glToError(0.0, "PL") > 0.0, "glToError(PL 0) is clamped strictly above zero");

    // The result is a probability in every branch.
    const char *types[3] = {"GQ", "GL", "PL"};
    bool inRange = true;
    for (int t = 0; t < 3; t++)
        for (double v = -60.0; v <= 200.0; v += 0.5)
        {
            double e = glToError(v, types[t]);
            if (!(e > 0.0 && e <= 1.0)) inRange = false;
        }
    ck(inRange, "glToError always returns a value in (0, 1]");

    // Monotonic: higher GQ means lower error.
    bool mono = true;
    double prev = 2.0;
    for (double q = 0.0; q <= 90.0; q += 5.0)
    {
        double e = glToError(q, "GQ");
        if (e > prev) mono = false;
        prev = e;
    }
    ck(mono, "glToError(GQ) decreases as quality increases");
}

// -------------------------------------------------- KDE mode finding ------
static void test_kde_helpers()
{
    // REGRESSION (B7): get_arg_max seeded its running maximum with
    // numeric_limits<double>::min(), which is the smallest POSITIVE double,
    // so an all-negative array returned index 0 regardless of the values.
    double allNeg[5] = {-5.0, -2.0, -9.0, -1.0, -7.0};
    ck(get_arg_max(allNeg, 5) == 3, "get_arg_max works on an all-negative array");
    double allPos[5] = {1.0, 4.0, 2.0, 9.0, 3.0};
    ck(get_arg_max(allPos, 5) == 3, "get_arg_max on positives");
    ck(get_arg_min(allNeg, 5) == 2, "get_arg_min on negatives");
    ck(get_arg_min(allPos, 5) == 0, "get_arg_min on positives");
    double one[1] = {42.0};
    ck(get_arg_max(one, 1) == 0 && get_arg_min(one, 1) == 0, "arg_max/arg_min on a single element");

    // A synthetic bimodal density with its trough at a known place.  Two
    // Gaussian bumps at x = -2 and x = +2; the minimum between the modes is
    // at x = 0 by symmetry.
    const int N = 401;
    double x[N], y[N];
    for (int i = 0; i < N; i++)
    {
        x[i] = -5.0 + 10.0 * double(i) / double(N - 1);
        y[i] = exp(-(x[i] + 2.0) * (x[i] + 2.0) / 0.5) + exp(-(x[i] - 2.0) * (x[i] - 2.0) / 0.5);
    }
    double cut = get_min_btw_modes(x, y, N, 20);
    ckd(cut, 0.0, 0.25, "get_min_btw_modes finds the trough of a symmetric bimodal density");
}

// ----------------------------------------------------- getMapInfo() -------
static void test_getMapInfo()
{
    // REGRESSION (B5): getMapInfo used scaffold->currentIndex as a
    // forward-only cursor.  A query behind the cursor fell through the scan
    // loop without ever assigning startIndex/endIndex, and interpolate() was
    // called on uninitialised stack values.  The map below is deliberately
    // NON-linear (a hotspot between 2000 and 3000) so that interpolating on
    // the wrong interval gives a visibly wrong answer rather than the right
    // one by accident.
    GenMapScaffold *sc = initGenMapScaffold(5);
    double gp[5] = {1.0, 2.0, 10.0, 11.0, 12.0};
    for (int i = 0; i < 5; i++)
    {
        sc->physicalPos[i] = 1000 * (i + 1);
        sc->geneticPos[i]  = gp[i];
        sc->ppos2index[1000 * (i + 1)] = i;
    }
    sc->chr = "chr1";
    sc->currentIndex = 0;
    int count = 0;

    // Ascending queries: the easy direction.
    ckd(getMapInfo(1500, sc, count), 1.5,  1e-12, "getMapInfo interpolates 1500 -> 1.5 cM");
    ckd(getMapInfo(2500, sc, count), 6.0,  1e-12, "getMapInfo interpolates 2500 across the hotspot");
    ckd(getMapInfo(4500, sc, count), 11.5, 1e-12, "getMapInfo interpolates 4500 -> 11.5 cM");

    // Now go BACKWARDS, behind the cursor.  Verified to discriminate: built
    // against the pre-fix objects this exact sequence returns 8.5 instead of
    // 1.5 (interpolation on a stale/uninitialised interval).
    ckd(getMapInfo(1500, sc, count), 1.5, 1e-12,
        "getMapInfo is correct for a query BEHIND the cursor (B5)");
    ckd(getMapInfo(2500, sc, count), 6.0, 1e-12,
        "getMapInfo is correct for a second backward query (B5)");

    // Exact scaffold positions must come back exactly.
    ckd(getMapInfo(3000, sc, count), 10.0, 1e-12, "getMapInfo returns an exact scaffold position");
    releaseGenMapScaffold(sc);
}

int main()
{
    printf("garlic unit tests\n");
    test_lod();
    test_interpolate();
    test_inGap();
    test_auto_fits();
    test_class_labels();
    test_chr_names();
    test_glToError();
    test_kde_helpers();
    test_getMapInfo();
    printf("%d checks, %d failures\n", checks, failures);
    return failures == 0 ? 0 : 1;
}
