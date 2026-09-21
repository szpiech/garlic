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
#include "garlic-features.h"
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <cstring>
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

// Several cases below exercise paths whose whole job is to report and refuse.
// Those messages are asserted in test/run_tests.sh, where the wording is what
// the user sees; here they are noise that buries the one line that matters.
// Swapping cerr's buffer is the portable way to mute them -- freopen("/dev/null")
// is spelled differently on Windows.
static streambuf *savedCerrBuf = NULL;
static ostringstream cerrSink;
static void quietErrors(bool on)
{
    if (on && savedCerrBuf == NULL) savedCerrBuf = cerr.rdbuf(cerrSink.rdbuf());
    else if (!on && savedCerrBuf != NULL)
    {
        cerr.rdbuf(savedCerrBuf);
        savedCerrBuf = NULL;
        cerrSink.str("");
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

// ------------------------------------------- enumeratePopulations() -----
// The ORDER is the contract, not just the contents: it decides which
// population is analysed first and, once per-population seeding exists, which
// seed each one gets.  A map-keyed implementation would return sorted order
// for free, which is why the labels below are deliberately chosen so that
// sorted and first-appearance differ.
static IndData *makeInd(const char **pops, int n)
{
    IndData *d = initIndData(n);
    for (int i = 0; i < n; i++)
    {
        char buf[32];
        snprintf(buf, sizeof(buf), "ind%d", i);
        d->indID[i] = buf;
        d->pop[i]   = pops[i];
    }
    return d;
}

static void test_enumeratePopulations()
{
    // "zeta" first, "alpha" second: sorted order would invert these.
    const char *mixed[] = {"zeta", "zeta", "alpha", "zeta", "alpha", "mid"};
    IndData *d = makeInd(mixed, 6);
    vector< pair<string, int> > got = enumeratePopulations(d);
    ck(got.size() == 3, "enumeratePopulations: three distinct labels");
    if (got.size() == 3)
    {
        cks(got[0].first, "zeta",  "enumeratePopulations: first appearance wins, not sort order");
        cks(got[1].first, "alpha", "enumeratePopulations: second by first appearance");
        cks(got[2].first, "mid",   "enumeratePopulations: third by first appearance");
        ck(got[0].second == 3, "enumeratePopulations: count for the first label");
        ck(got[1].second == 2, "enumeratePopulations: count for the second label");
        ck(got[2].second == 1, "enumeratePopulations: count for the third label");
    }
    int total = 0;
    for (unsigned int i = 0; i < got.size(); i++) total += got[i].second;
    ck(total == d->nind, "enumeratePopulations: counts sum to nind");
    releaseIndData(d);

    const char *one[] = {"36", "36", "36"};
    IndData *s1 = makeInd(one, 3);
    vector< pair<string, int> > g1 = enumeratePopulations(s1);
    ck(g1.size() == 1 && g1[0].second == 3, "enumeratePopulations: single population");
    releaseIndData(s1);

    // Labels are compared as strings, so these are three populations, not one.
    const char *casey[] = {"pop", "POP", "Pop"};
    IndData *s2 = makeInd(casey, 3);
    ck(enumeratePopulations(s2).size() == 3, "enumeratePopulations: labels are case sensitive");
    releaseIndData(s2);

    ck(enumeratePopulations(NULL).empty(), "enumeratePopulations: NULL gives an empty list");
}

// ---------------------------------------------- subsetDataByIndex() -----
// The gather behind both --kde-subsample and, shortly, per-population
// selection.  Two properties matter beyond "it copies something":
//
//   ORDER follows keepInd, not the natural order -- the caller controls the
//   individual ordering of the result, which is what lets a population be
//   gathered in file order.
//
//   It DEEP copies.  subsetData used to alias IndData::pop, which leaked one
//   array and double-freed another (the --auto-winsize abort, B1).  That exact
//   bug is no longer expressible -- IndData holds std::vector and Matrix<T>
//   since the ownership refactor, and both deep copy on assignment -- so the
//   assertions below are a guard against a future reintroduction of raw
//   pointer sharing, not a live reproduction.  Injecting the original
//   `newIndData->pop = indData->pop` now fails the ORDER check instead, which
//   is worth knowing: it is the gather checks that discriminate today.
static void test_subsetDataByIndex()
{
    const int nind = 5, nloci = 4;
    IndData *ind = initIndData(nind);
    const char *ids[]  = {"i0", "i1", "i2", "i3", "i4"};
    const char *pops[] = {"A",  "B",  "A",  "C",  "B"};
    for (int i = 0; i < nind; i++)
    {
        ind->indID[i] = ids[i];
        ind->pop[i]   = pops[i];
        ind->sex[i]   = (i % 2) + 1;
    }

    //Heap allocated because releaseHapData(vector<HapData*>*) deletes the
    //vector itself, as it does for the real per-chromosome vectors.
    vector< HapData * > *hv = new vector< HapData * >;
    HapData *hap = initHapData(nind, nloci, true);
    for (int l = 0; l < nloci; l++)
        for (int i = 0; i < nind; i++)
        {
            hap->data[l][i]      = geno_t(10 * l + i);
            hap->firstCopy[l][i] = (unsigned char)(i % 2);
        }
    hv->push_back(hap);

    // the "B" individuals, in file order: 1 then 4
    vector<int> keep;
    keep.push_back(1); keep.push_back(4);

    vector< HapData * > *subHap = NULL;
    IndData *subInd = NULL;
    subsetDataByIndex(hv, NULL, ind, keep, &subHap, NULL, &subInd, false, true);

    ck(subInd != NULL && subInd->nind == 2, "subsetDataByIndex: subset has the requested size");
    cks(subInd->indID[0], "i1", "subsetDataByIndex: first individual follows keepInd");
    cks(subInd->indID[1], "i4", "subsetDataByIndex: second individual follows keepInd");
    cks(subInd->pop[0],   "B",  "subsetDataByIndex: population label follows");
    ck(subInd->sex[1] == 1,     "subsetDataByIndex: sex follows");
    ck(subHap->at(0)->nind == 2 && subHap->at(0)->nloci == nloci,
       "subsetDataByIndex: genotype matrix is nloci x |keepInd|");

    bool geno = true, fc = true;
    for (int l = 0; l < nloci; l++)
    {
        if (subHap->at(0)->data[l][0] != geno_t(10 * l + 1)) geno = false;
        if (subHap->at(0)->data[l][1] != geno_t(10 * l + 4)) geno = false;
        if (subHap->at(0)->firstCopy[l][0] != 1) fc = false;
        if (subHap->at(0)->firstCopy[l][1] != 0) fc = false;
    }
    ck(geno, "subsetDataByIndex: genotypes are the chosen individuals' own");
    ck(fc,   "subsetDataByIndex: firstCopy follows the same individuals");

    // reversed order must reverse the result, not re-sort it
    vector<int> rev;
    rev.push_back(4); rev.push_back(1);
    vector< HapData * > *revHap = NULL;
    IndData *revInd = NULL;
    subsetDataByIndex(hv, NULL, ind, rev, &revHap, NULL, &revInd, false, true);
    cks(revInd->indID[0], "i4", "subsetDataByIndex: keepInd order is honoured, not sorted");
    ck(revHap->at(0)->data[0][0] == geno_t(4),
       "subsetDataByIndex: genotypes follow the reversed order too");
    releaseHapData(revHap); releaseIndData(revInd);

    // Mutate the subset, then check the WHOLE original.  Checking one cell is
    // not enough: an alias shifts which element a given subset index lands on,
    // so a single-cell check can read an untouched one and pass.
    bool origOK = true;
    subInd->pop[0]   = "MUTATED";
    subInd->indID[0] = "MUTATED";
    subHap->at(0)->data[0][0] = geno_t(99);
    for (int i = 0; i < nind; i++)
        if (ind->pop[i] != pops[i] || ind->indID[i] != ids[i]) origOK = false;
    for (int l = 0; l < nloci; l++)
        for (int i = 0; i < nind; i++)
            if (hap->data[l][i] != geno_t(10 * l + i)) origOK = false;
    ck(origOK, "subsetDataByIndex: mutating the subset leaves the whole original intact");

    releaseHapData(subHap); releaseIndData(subInd);
    origOK = true;
    for (int i = 0; i < nind; i++)
        if (ind->pop[i] != pops[i] || ind->indID[i] != ids[i]) origOK = false;
    for (int l = 0; l < nloci; l++)
        for (int i = 0; i < nind; i++)
            if (hap->data[l][i] != geno_t(10 * l + i)) origOK = false;
    ck(origOK, "subsetDataByIndex: releasing the subset leaves the whole original intact");

    // an out-of-range index is an error, not an out-of-bounds read
    vector<int> bad;
    bad.push_back(0); bad.push_back(nind);
    vector< HapData * > *badHap = NULL;
    IndData *badInd = NULL;
    bool threw = false;
    try { subsetDataByIndex(hv, NULL, ind, bad, &badHap, NULL, &badInd, false, true); }
    catch (...) { threw = true; }
    ck(threw, "subsetDataByIndex: an out-of-range index is rejected");

    releaseHapData(hv);
    releaseIndData(ind);
}

// ------------------------------------------------ allele frequency -----
// alleleFrequency() is the rule four call sites had each written out:
// loadTPEDData, loadVCFData, freqOnly, freqOnlyVCF.
static void test_alleleFrequency()
{
    ckd(alleleFrequency(3, 6, 0, NULL), 0.5,  1e-12, "alleleFrequency: 3 of 6");
    ckd(alleleFrequency(5, 5, 0, NULL), 1.0,  1e-12, "alleleFrequency: fixed");
    ckd(alleleFrequency(0, 8, 0, NULL), 0.0,  1e-12, "alleleFrequency: absent");
    // total == 0 is the no-data case and must not divide; every caller relies
    // on this branch for a locus where nobody was called.
    ckd(alleleFrequency(0, 0, 0, NULL), 0.0,  1e-12, "alleleFrequency: no observations gives 0, not NaN");
    // r is documented as ignorable when nresample <= 0; NULL above proves it.
    ckd(alleleFrequency(1, 4, -1, NULL), 0.25, 1e-12, "alleleFrequency: negative nresample is no resampling");
}

// ------------------------------------------ calcFreqDataForIndices() ----
// Per-population frequencies, from the loaded genotype matrix.  The property
// that matters is not just "it divides correctly" but that the SAME locus
// gets DIFFERENT frequencies in different subsets -- that is the whole point
// of analysing populations separately.
static void test_calcFreqDataForIndices()
{
    const int nind = 4, nloci = 3;
    vector< HapData * > *hv = new vector< HapData * >;
    HapData *hap = initHapData(nind, nloci, false);
    //                 ind:   0            1            2            3
    geno_t g[3][4] = {{ 0,           1,           2,           GENO_MISSING },
                      { 2,           2,           2,           2            },
                      { GENO_MISSING, GENO_MISSING, 0,          1            }};
    for (int l = 0; l < nloci; l++)
        for (int i = 0; i < nind; i++) hap->data[l][i] = g[l][i];
    hv->push_back(hap);

    vector<int> all;
    for (int i = 0; i < nind; i++) all.push_back(i);
    vector< FreqData * > *fAll = calcFreqDataForIndices(hv, all, 0);
    ck(fAll->size() == 1 && fAll->at(0)->nloci == nloci,
       "calcFreqDataForIndices: one FreqData per chromosome, nloci long");
    ckd(fAll->at(0)->freq[0], 3.0 / 6.0, 1e-12, "calcFreqDataForIndices: missing individual excluded from both numerator and denominator");
    ckd(fAll->at(0)->freq[1], 1.0,       1e-12, "calcFreqDataForIndices: fixed locus");
    ckd(fAll->at(0)->freq[2], 1.0 / 4.0, 1e-12, "calcFreqDataForIndices: two missing individuals excluded");

    // the same loci, two different halves of the cohort
    vector<int> firstTwo, lastTwo;
    firstTwo.push_back(0); firstTwo.push_back(1);
    lastTwo.push_back(2);  lastTwo.push_back(3);
    vector< FreqData * > *fA = calcFreqDataForIndices(hv, firstTwo, 0);
    vector< FreqData * > *fB = calcFreqDataForIndices(hv, lastTwo,  0);

    ckd(fA->at(0)->freq[0], 1.0 / 4.0, 1e-12, "calcFreqDataForIndices: locus 0 in the first half");
    ckd(fB->at(0)->freq[0], 1.0,       1e-12, "calcFreqDataForIndices: locus 0 in the second half");
    ck(fA->at(0)->freq[0] != fB->at(0)->freq[0],
       "calcFreqDataForIndices: a locus gets different frequencies in different subsets");
    // every individual of the first half is missing at locus 2
    ckd(fA->at(0)->freq[2], 0.0,       1e-12, "calcFreqDataForIndices: all-missing subset gives 0, not NaN");
    ckd(fB->at(0)->freq[2], 1.0 / 4.0, 1e-12, "calcFreqDataForIndices: locus 2 in the second half");
    // a fixed locus is fixed in any subset
    ckd(fA->at(0)->freq[1], 1.0, 1e-12, "calcFreqDataForIndices: fixed locus is fixed in a subset too");

    releaseFreqData(fAll); releaseFreqData(fA); releaseFreqData(fB);
    releaseHapData(hv);
}

// ------------------------------------------------ addAlleleCounts() -----
// What each genotype code contributes to an allele frequency.  This is the
// rule that makes the TPED and VCF paths agree about half calls, so it is
// pinned rather than left implicit in four accumulation loops.
static void test_addAlleleCounts()
{
    double na, t;
    #define ACC(g) (na = 0, t = 0, addAlleleCounts((g), na, t))
    ACC(0); ck(na == 0 && t == 2, "addAlleleCounts: 0 contributes 0 of 2");
    ACC(1); ck(na == 1 && t == 2, "addAlleleCounts: 1 contributes 1 of 2");
    ACC(2); ck(na == 2 && t == 2, "addAlleleCounts: 2 contributes 2 of 2");
    ACC(GENO_HALF_COUNTED); ck(na == 1 && t == 1, "addAlleleCounts: a half call of the counted allele contributes 1 of 1");
    ACC(GENO_HALF_OTHER);   ck(na == 0 && t == 1, "addAlleleCounts: a half call of the other allele contributes 0 of 1");
    ACC(GENO_MISSING);      ck(na == 0 && t == 0, "addAlleleCounts: a fully missing call contributes nothing");

    // the two half codes together are exactly one heterozygote's worth, which
    // is why a half call cannot be treated as missing without biasing the
    // frequency towards whatever the called individuals carry
    na = 0; t = 0;
    addAlleleCounts(GENO_HALF_COUNTED, na, t);
    addAlleleCounts(GENO_HALF_OTHER,   na, t);
    ck(na == 1 && t == 2, "addAlleleCounts: two opposite half calls sum to one heterozygote");
    #undef ACC
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

// ----------------------------------------------------- canonChrKey() ------
// The MATCH key, as distinct from the display name above.  Everything that
// compares one chromosome name against another goes through this: --chr, the
// --build centromere table, the frequency file, the sex-chromosome detector.
// Display names are untouched, which is what keeps existing output
// byte-identical.
static void test_canonChrKey()
{
    // Case and the chr prefix are both ignored, in every combination.
    cks(canonChrKey("chrX"), "x", "canonChrKey chrX");
    cks(canonChrKey("chrx"), "x", "canonChrKey chrx");
    cks(canonChrKey("CHRX"), "x", "canonChrKey CHRX");
    cks(canonChrKey("ChrX"), "x", "canonChrKey ChrX");
    cks(canonChrKey("X"),    "x", "canonChrKey X");
    cks(canonChrKey("x"),    "x", "canonChrKey x");

    // Leading zeros go only when what is left is all digits: 01 and 1 are the
    // same chromosome in any file, scaffold_007 is not scaffold_7.
    cks(canonChrKey("chr1"),  "1", "canonChrKey chr1");
    cks(canonChrKey("01"),    "1", "canonChrKey 01");
    cks(canonChrKey("chr01"), "1", "canonChrKey chr01");
    cks(canonChrKey("0"),     "0", "canonChrKey 0 does not become empty");
    cks(canonChrKey("00"),    "0", "canonChrKey 00");
    cks(canonChrKey("scaffold_007"), "scaffold_007", "canonChrKey keeps zeros in a non-numeric name");

    // Non-standard but real names: chicken 4A, chimp 2A, linkage groups.
    cks(canonChrKey("chr4a"), "4a", "canonChrKey chr4a");
    cks(canonChrKey("4a"),    "4a", "canonChrKey 4a");
    cks(canonChrKey("Chr4A"), "4a", "canonChrKey Chr4A");
    cks(canonChrKey("LG12"),  "lg12", "canonChrKey LG12");
    cks(canonChrKey("NC_000023.11"), "nc_000023.11", "canonChrKey leaves an accession alone");

    // Only ONE leading prefix is stripped, and only when something follows:
    // a chromosome named "chr" keeps its name rather than becoming empty.
    cks(canonChrKey("chr"),       "chr",  "canonChrKey chr stays chr");
    cks(canonChrKey("chrchr1"),   "chr1", "canonChrKey strips one prefix only");
    cks(canonChrKey("contig7"),   "contig7", "canonChrKey does not eat a c that is not chr");

    // The lowercasing is an explicit 'A'..'Z' test rather than tolower(),
    // which is locale dependent, or std::tolower(char), which is undefined
    // for negative char values.  Every byte of a UTF-8 name above 0x7F is
    // negative on a signed-char platform; it must pass through untouched.
    string utf8 = "chr\xc3\xa9";           // "chré"
    string wantUtf8 = "\xc3\xa9";
    cks(canonChrKey(utf8), wantUtf8, "canonChrKey passes high bytes through unchanged");
}

// ------------------------------------------------------- SexModel --------
// What a chromosome IS.  The role is independent of the sex determination
// system -- a chromosome named X or Z is the shared one either way -- and the
// system contributes only which sex CODE is heterogametic, which matters to
// the calling step rather than to this table.
static void test_sex_model()
{
    ck(parseSexSystem("xy") == SEX_SYSTEM_XY,     "parseSexSystem xy");
    ck(parseSexSystem("zw") == SEX_SYSTEM_ZW,     "parseSexSystem zw");
    ck(parseSexSystem("none") == SEX_SYSTEM_NONE, "parseSexSystem none");
    ck(parseSexSystem("XY") < 0,                  "parseSexSystem rejects XY (lower case only)");
    ck(parseSexSystem("xo") < 0,                  "parseSexSystem rejects an unknown system");

    ck(isDetectedSexChrKey("x") && isDetectedSexChrKey("y"), "detector fires on x and y");
    ck(isDetectedSexChrKey("z") && isDetectedSexChrKey("w"), "detector fires on z and w");
    ck(isDetectedSexChrKey("m") && isDetectedSexChrKey("mt"), "detector fires on m and mt");
    ck(!isDetectedSexChrKey("23"), "detector does not fire on a bare number");
    ck(!isDetectedSexChrKey("4a") && !isDetectedSexChrKey("lg12"),
       "detector does not fire on an ordinary chromosome name");

    // The numbers are AMBIGUOUS, not autosomal: sex-linked under PLINK's
    // human coding, ordinary autosomes in a species with that many
    // chromosomes.  This is the distinction that keeps chicken chr24 safe.
    ck(isAmbiguousNumericChrKey("23") && isAmbiguousNumericChrKey("26"),
       "23 and 26 are ambiguous");
    ck(!isAmbiguousNumericChrKey("22") && !isAmbiguousNumericChrKey("27"),
       "22 and 27 are not ambiguous");
    ck(isHumanBuild("hg19") && isHumanBuild("t2t-chm13"), "human builds recognised");
    ck(!isHumanBuild("none"), "no build is not a human build");

    // A 3-chromosome fixture: one autosome, one X, one Y.
    const char *names[3] = {"chr1", "chrX", "chrY"};
    vector< MapData * > map3;
    for (int i = 0; i < 3; i++)
    {
        MapData *md = initMapData(1);
        md->chr = names[i];
        md->physicalPos[0] = 1000;
        map3.push_back(md);
    }
    IndData *ind = initIndData(2);
    ind->indID[0] = "a"; ind->pop[0] = "P"; ind->sex[0] = 1;
    ind->indID[1] = "b"; ind->pop[1] = "P"; ind->sex[1] = 2;

    vector<string> none;
    SexModel m;
    quietErrors(true);   //every refusal below logs; the wording is asserted in run_tests.sh

    // Unset system with a sex chromosome present: refuses rather than
    // guessing.  This is the whole point of the flag being unset by default.
    ck(buildSexModel(m, &map3, ind, SEX_SYSTEM_UNSET, none, none, none, false, false) < 0,
       "unset system refuses when chrX is present");

    // Declared XY: conventional names resolve, and nothing else moves.
    ck(buildSexModel(m, &map3, ind, SEX_SYSTEM_XY, none, none, none, false, false) == 0,
       "xy resolves chrX and chrY");
    ck(m.role[0] == CHR_AUTOSOME,       "chr1 is an autosome under xy");
    ck(m.role[1] == CHR_SEX_SHARED,     "chrX is the shared sex chromosome under xy");
    ck(m.role[2] == CHR_SEX_DEGENERATE, "chrY is the degenerate sex chromosome under xy");

    // Declared ZW against X/Y data: the conventional ZW names are z and w, so
    // chrX and chrY are left unaccounted for and the run stops.  A silent
    // reinterpretation here would invert who is hemizygous.
    ck(buildSexModel(m, &map3, ind, SEX_SYSTEM_ZW, none, none, none, false, false) < 0,
       "zw refuses X/Y data rather than reinterpreting it");

    // none asserts everything is diploid, and says so about every chromosome.
    ck(buildSexModel(m, &map3, ind, SEX_SYSTEM_NONE, none, none, none, false, false) == 0,
       "none accepts chrX");
    ck(m.role[1] == CHR_AUTOSOME, "none makes chrX an autosome");

    // --autosomes-only drops what is not an autosome, so it can use the
    // conventional names of both systems without knowing which sex is
    // heterogametic -- it is not calling anything on them.
    ck(buildSexModel(m, &map3, ind, SEX_SYSTEM_UNSET, none, none, none, false, true) == 0,
       "--autosomes-only resolves conventional names with no system");
    ck(m.role[1] == CHR_SEX_SHARED && m.role[2] == CHR_SEX_DEGENERATE,
       "--autosomes-only assigns both sex-chromosome roles");

    // A declaration cannot be made without the system, because the system is
    // the part the data cannot supply.
    vector<string> justX; justX.push_back("chrX");
    ck(buildSexModel(m, &map3, ind, SEX_SYSTEM_UNSET, justX, none, none, false, false) < 0,
       "--sex-chr without --sex-system is refused");
    // ...and naming something absent is an error, as with --chr.
    vector<string> absent; absent.push_back("chrQ");
    ck(buildSexModel(m, &map3, ind, SEX_SYSTEM_XY, absent, none, none, false, false) < 0,
       "--sex-chr naming a chromosome not in the data is refused");

    for (int i = 0; i < 3; i++) releaseMapData(map3[i]);
    releaseIndData(ind);

    // Numbers: ambiguous without a human build, PLINK's human codes with one.
    const char *nums[4] = {"chr23", "chr24", "chr25", "chr26"};
    vector< MapData * > map4;
    for (int i = 0; i < 4; i++)
    {
        MapData *md = initMapData(1);
        md->chr = nums[i];
        md->physicalPos[0] = 1000;
        map4.push_back(md);
    }
    IndData *ind2 = initIndData(1);
    ind2->indID[0] = "a"; ind2->pop[0] = "P"; ind2->sex[0] = 1;

    ck(buildSexModel(m, &map4, ind2, SEX_SYSTEM_XY, none, none, none, false, false) < 0,
       "bare 23-26 are refused without a human build");
    ck(buildSexModel(m, &map4, ind2, SEX_SYSTEM_XY, none, none, none, true, false) == 0,
       "a human build licenses PLINK's numbering");
    ck(m.role[0] == CHR_SEX_SHARED,     "23 is X under a human build");
    ck(m.role[1] == CHR_SEX_DEGENERATE, "24 is Y under a human build");
    ck(m.role[2] == CHR_PAR,            "25 is the PAR under a human build");
    ck(m.role[3] == CHR_HAPLOID,        "26 is the mitochondrion under a human build");

    // Same numbers in a species that simply has 26 chromosomes: the user says
    // so, and they stay autosomes.
    ck(buildSexModel(m, &map4, ind2, SEX_SYSTEM_NONE, none, none, none, false, false) == 0,
       "--sex-system none accepts 23-26 as autosomes");
    ck(m.role[0] == CHR_AUTOSOME && m.role[3] == CHR_AUTOSOME,
       "23 and 26 stay autosomes under none");

    // A non-conventional name, declared explicitly: the ZW species whose Z is
    // called LG12.  Nothing about the role depends on the system.
    map4[0]->chr = "LG12";
    vector<string> lg; lg.push_back("lg12");     // matched by key, any spelling
    ck(buildSexModel(m, &map4, ind2, SEX_SYSTEM_ZW, lg, none, none, false, false) < 0,
       "declaring LG12 does not excuse the other ambiguous numbers");
    ck(m.role[0] == CHR_SEX_SHARED, "LG12 is declared the shared sex chromosome");

    quietErrors(false);
    for (int i = 0; i < 4; i++) releaseMapData(map4[i]);
    releaseIndData(ind2);
}

// ------------------------------------------- winsizeForChr ---------------
// The window size stops being one number for the run once --sexchr-winsize
// exists.  Every stage that takes one -- the LOD windows, the LD matrix, the
// overlap threshold -- asks this, so a wrong answer here is a silent change
// of window size on the wrong chromosome.
static void test_winsize_for_chr()
{
    vector<ChrRole> role;
    role.push_back(CHR_AUTOSOME);
    role.push_back(CHR_SEX_SHARED);
    role.push_back(CHR_SEX_DEGENERATE);

    ck(winsizeForChr(60, 30, &role, 0) == 60, "an autosome keeps the run's window size");
    ck(winsizeForChr(60, 30, &role, 1) == 30, "the shared sex chromosome takes its own");
    ck(winsizeForChr(60, 30, &role, 2) == 60,
       "a degenerate sex chromosome is not the one --sexchr-winsize names");

    // Unset is the default and has to mean "the run's size", not zero.
    ck(winsizeForChr(60, 0, &role, 1) == 60, "no override leaves the sex chromosome alone");
    ck(winsizeForChr(60, -1, &role, 1) == 60, "a nonsensical override is ignored, not applied");

    // The role vector is absent on paths that never built one, and an index
    // past its end is reachable while chromosomes are being filtered.
    ck(winsizeForChr(60, 30, NULL, 1) == 60, "no role table means the run's size");
    ck(winsizeForChr(60, 30, &role, 7) == 60, "an index past the table means the run's size");
}

// ------------------------------------------- LD with half calls ---------
// Both LD statistics skip a genotype that is not called, and a half call is
// not called.  On the shared sex chromosome that removes every heterogametic
// individual, so a subsample can contain nobody with genotypes at both loci
// -- and 0/0 used to leave as a NaN, which the weighted LOD turns into
// 1/NaN and a silently uncalled stretch of chromosome.
static void test_ld_no_shared_genotypes()
{
    HapData *h = initHapData(6, 2, false);
    // 0-3 diploid and a mix of homozygous and heterozygous, so the genotype
    // frequency is strictly between 0 and 1 and the statistic is evaluated
    // rather than short-circuited; 4 and 5 hemizygous at both loci.
    geno_t L0[6] = {0, 1, 2, 1, GENO_HALF_COUNTED, GENO_HALF_OTHER};
    geno_t L1[6] = {2, 1, 0, 1, GENO_HALF_OTHER,   GENO_HALF_COUNTED};
    for (int i = 0; i < 6; i++) { h->data[0][i] = L0[i]; h->data[1][i] = L1[i]; }

    GenoFreqData *g = calculateGenoFreq(h);
    ck(g->homFreq[0] > 0 && g->homFreq[0] < 1,
       "half calls are left out of the genotype frequency");

    FreqData *f = initFreqData(2);
    f->freq[0] = 0.5; f->freq[1] = 0.5;

    int all[6] = {0, 1, 2, 3, 4, 5};
    int hemi[2] = {4, 5};
    ck(hr2(h, g, 0, 1, all, 6) == hr2(h, g, 0, 1, all, 6), "hr2 over the cohort is a number");
    double bh = hr2(h, g, 0, 1, hemi, 2);
    double br = r2(h, f, 0, 1, hemi, 2);
    ck(bh == bh && bh == 0, "hr2 with no shared genotypes is 0, not NaN");
    ck(br == br && br == 0, "r2 with no shared genotypes is 0, not NaN");

    releaseFreqData(f);
    releaseGenoFreq(g);
    releaseHapData(h);
}

// ------------------------------------------- ChrLengths ------------------
// The lengths are matched to whatever the data calls a chromosome, which is
// the whole reason a VCF's decoys and alt contigs need no filtering: nothing
// looks them up.  The lookup therefore has to be the same canonicalisation
// the rest of the program uses, not a string compare.
static void test_chr_lengths()
{
    ChrLengths L;
    L.add("chr1", 249250621);
    L.add("23", 155270560);

    ck(L.get("chr1") == 249250621, "a length is found under the name it was added with");
    ck(L.get("1") == 249250621,    "and under the same chromosome spelled without chr");
    ck(L.get("CHR1") == 249250621, "and case-insensitively");
    ck(L.get("chr23") == 155270560, "and the other way round, chr23 finding 23");
    //0 is the "not here" answer, and a real length is never 0.
    ck(L.get("chr2") == 0,          "a chromosome not in the table has no length");
    ck(L.get("chrUn_KI270302v1") == 0, "nor does an unplaced contig");
    ck(!L.empty() && L.chromosomes().size() == 2, "two chromosomes, in insertion order");

    ChrLengths E;
    ck(E.empty() && E.get("chr1") == 0, "an empty table answers nothing");

    //The built-in table is the one --build supplies, and the hg19 autosomes
    //have to sum to the published denominator or the point of the mode is
    //lost.
    ChrLengths B;
    ck(builtinChrLengths("hg19", B), "hg19 has a built-in length table");
    long long total = 0;
    for (int c = 1; c <= 22; c++)
    {
        ostringstream ss; ss << "chr" << c;
        total += (long long)B.get(ss.str());
    }
    ck(total == 2881033286LL, "the hg19 autosomes sum to the published 2,881,033,286 bp");
    ck(B.get("chrX") == 155270560, "and chrX is the published 155,270,560 bp");
    //Y is a degenerate sex chromosome; it is in no denominator, so carrying a
    //length for it would be a number nothing reads.
    ck(B.get("chrY") == 0, "no length is carried for the degenerate sex chromosome");
    ck(B.source().compare("build hg19") == 0, "the table says where it came from");

    //An unknown build must leave the caller's table untouched and say so,
    //rather than reporting whatever was already in it as its own.
    ChrLengths P;
    P.add("chr1", 12345);
    ck(!builtinChrLengths("none", P), "an unknown build has no table");
    ck(P.get("chr1") == 12345 && P.source().empty(),
       "and does not claim the table it was handed");
}

// ------------------------------------------- ExcludedRegions -------------
// A chromosome has a SET of excluded regions, not one.  Humans already have
// two pseudoautosomal regions; the container the centromere table uses is one
// interval per chromosome and would silently keep only the last.
static void test_excluded_regions()
{
    ExcludedRegions e;
    e.add("chrX", 100, 200);
    e.add("X",    500, 600);          // same chromosome, any spelling
    ck(e.finalise() == 0, "two disjoint regions on one chromosome are accepted");
    const vector<Interval> *v = e.get("chrx");
    ck(v != NULL && v->size() == 2, "both regions survive; neither overwrites the other");
    ck(v != NULL && v->size() == 2 && v->at(0).start == 100 && v->at(1).start == 500,
       "regions come back sorted by start");

    ck(e.contains("chrX", 100) && e.contains("chrX", 200),
       "contains() is inclusive at both ends of a region");
    ck(!e.contains("chrX", 99) && !e.contains("chrX", 201),
       "contains() excludes the bases either side");
    ck(!e.contains("chrX", 350), "contains() is false in the gap between two regions");
    ck(!e.contains("chr1", 150), "regions do not leak to another chromosome");

    // Overlap is what the FROH denominator subtracts, so it must be clipped
    // and must not double-count.
    ck(e.overlap("chrX", 1, 1000) == 202, "overlap sums both regions");
    ck(e.overlap("chrX", 150, 1000) == 152, "overlap clips to the query start");
    ck(e.overlap("chrX", 1, 150) == 51, "overlap clips to the query end");
    ck(e.overlap("chrX", 250, 400) == 0, "overlap is zero between the regions");
    ck(e.overlap("chr1", 1, 1000) == 0, "overlap is zero on another chromosome");

    // Overlapping and abutting regions merge, so nothing is subtracted twice.
    ExcludedRegions m;
    m.add("chrX", 100, 200);
    m.add("chrX", 150, 300);
    m.add("chrX", 301, 400);
    ck(m.finalise() == 0, "overlapping regions are accepted");
    const vector<Interval> *mv = m.get("chrX");
    ck(mv != NULL && mv->size() == 1, "overlapping and abutting regions merge into one");
    ck(mv != NULL && mv->size() == 1 && mv->at(0).start == 100 && mv->at(0).end == 400,
       "the merged region spans all three");
    ck(m.overlap("chrX", 1, 1000) == 301, "a merged region is counted once");

    ExcludedRegions bad;
    bad.add("chrX", 500, 100);
    quietErrors(true);
    ck(bad.finalise() < 0, "a region that ends before it starts is refused");
    quietErrors(false);

    // Both input spellings reach the same place.
    ExcludedRegions p;
    vector<string> specs;
    specs.push_back("chrX:60001-2699520,chrX:154931044-155260560");   // one token
    ck(parsePARSpecs(specs, p) == 0 && p.finalise() == 0, "a comma-separated list parses");
    const vector<Interval> *pv = p.get("chrX");
    ck(pv != NULL && pv->size() == 2, "a comma-separated list gives two regions");
    ck(pv != NULL && pv->size() == 2 && pv->at(0).end == 2699520 && pv->at(1).start == 154931044,
       "the parsed coordinates are the ones given");

    ExcludedRegions q;
    vector<string> two;
    two.push_back("chrX:1-10"); two.push_back("chrX:20-30");          // two tokens
    ck(parsePARSpecs(two, q) == 0 && q.finalise() == 0, "a space-separated list parses");
    ck(q.get("chrX") != NULL && q.get("chrX")->size() == 2, "both tokens are kept");

    ExcludedRegions junk;
    vector<string> bads;
    bads.push_back("chrX-100-200");
    quietErrors(true);
    ck(parsePARSpecs(bads, junk) < 0, "a spec without chr:start-end is refused");
    quietErrors(false);
}

// ------------------------------------------- checkChrKeyCollisions() ------
// Two display names that differ only by case or a chr prefix are two
// independent MapData entries, and every matcher now treats them as one name.
// The same display name twice means one chromosome whose rows are not
// contiguous, which the readers split silently.
static void test_chr_key_collisions()
{
    vector< MapData * > m;
    const char *ok[2] = {"chr1", "chr2"};
    for (int i = 0; i < 2; i++)
    {
        MapData *md = initMapData(1);
        md->chr = ok[i];
        m.push_back(md);
    }
    ck(checkChrKeyCollisions(&m) == 0, "distinct chromosomes pass the collision check");

    quietErrors(true);
    m[1]->chr = "CHR1";
    ck(checkChrKeyCollisions(&m) < 0, "chr1 and CHR1 collide");

    m[1]->chr = "chr1";
    ck(checkChrKeyCollisions(&m) < 0, "the same name twice is a non-contiguous chromosome");

    m[1]->chr = "chr01";
    ck(checkChrKeyCollisions(&m) < 0, "chr1 and chr01 collide");
    quietErrors(false);

    for (unsigned int i = 0; i < m.size(); i++) releaseMapData(m[i]);
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

// ------------------------------------------------------- keepSites() -------
// The site-retention predicate.  This is the case the review's D4 filtering
// item was about: the predicate used to be re-derived inline in ten separate
// functions and could not be tested at all -- reaching it meant running the
// whole pipeline and inferring the outcome from the locus count of the output.
// Now it is one function taking a hand-built scaffold.
static void test_keepSites()
{
    // Six sites at 1000..6000.  Frequencies make site 0 monomorphic (0.0) and
    // site 5 fixed (1.0); both must be dropped in every mode.
    const int N = 6;
    double freqs[N] = {0.0, 0.25, 0.5, 0.5, 0.75, 1.0};

    MapData *md = initMapData(N);
    FreqData *fd = initFreqData(N);
    for (int i = 0; i < N; i++)
    {
        md->physicalPos[i] = 1000 * (i + 1);
        md->geneticPos[i]  = i;
        fd->freq[i]        = freqs[i];
    }
    md->chr = "chr1";

    // --- NULL scaffold: monomorphic filter only -------------------------
    vector<int> keep = keepSites(md, fd, NULL);
    ck(keep.size() == 4, "keepSites drops monomorphic and fixed sites");
    ck(keep.size() == 4 && keep[0] == 1 && keep[1] == 2 && keep[2] == 3 && keep[3] == 4,
       "keepSites returns ORIGINAL indices of the retained sites, in order");

    // --- scaffold spanning 2500..5500, no centromere ---------------------
    // Site 1 (1000) is below the map; sites 2,3,4 are inside.
    GenMapScaffold *sc = initGenMapScaffold(2);
    sc->physicalPos[0] = 2500;
    sc->physicalPos[1] = 5500;
    sc->geneticPos[0]  = 0.0;
    sc->geneticPos[1]  = 1.0;
    sc->chr = "chr1";
    sc->centroStart = 0;
    sc->centroEnd   = 0;

    keep = keepSites(md, fd, sc);
    ck(keep.size() == 3 && keep[0] == 2 && keep[1] == 3 && keep[2] == 4,
       "keepSites drops sites below the scaffold's first position");

    // Shrink the map to 2500..3500: sites 3 (4000) and 4 (5000) are now both
    // above it, leaving only site 2 (3000).
    sc->physicalPos[1] = 3500;
    keep = keepSites(md, fd, sc);
    ck(keep.size() == 1 && keep[0] == 2,
       "keepSites drops sites above the scaffold's last position");

    // --- centromere gap ---------------------------------------------------
    sc->physicalPos[1] = 5500;
    sc->centroStart = 2500;
    sc->centroEnd   = 3500;
    keep = keepSites(md, fd, sc);
    ck(keep.size() == 2 && keep[0] == 3 && keep[1] == 4,
       "keepSites drops sites strictly inside the centromere gap");

    // The gap test is STRICT on both sides (pos > start && pos < end), so a
    // site sitting exactly on an edge is retained.  Widen the map to 0..99999
    // so only the gap decides, and put the gap edges exactly on sites 1 (2000)
    // and 3 (4000): both survive, and site 2 (3000) between them does not.
    sc->physicalPos[0] = 0;
    sc->physicalPos[1] = 99999;
    sc->centroStart = 2000;
    sc->centroEnd   = 4000;
    keep = keepSites(md, fd, sc);
    ck(keep.size() == 3 && keep[0] == 1 && keep[1] == 3 && keep[2] == 4,
       "keepSites RETAINS sites exactly at both gap boundaries (predicate is strict)");

    // --- every site dropped ----------------------------------------------
    for (int i = 0; i < N; i++) fd->freq[i] = 0.0;
    keep = keepSites(md, fd, NULL);
    ck(keep.empty(), "keepSites returns an empty list when nothing is retained");

    // --- every site retained ---------------------------------------------
    for (int i = 0; i < N; i++) fd->freq[i] = 0.5;
    sc->centroStart = 0; sc->centroEnd = 0;
    sc->physicalPos[0] = 0; sc->physicalPos[1] = 99999;
    keep = keepSites(md, fd, sc);
    ck(keep.size() == (size_t)N, "keepSites retains every site when nothing excludes any");

    // The invariant the old code could violate: the count that sizes the
    // destination and the indices that fill it come from the same list, so
    // the largest index is always addressable within keep.size() rows.
    bool monotone = true;
    for (size_t k = 1; k < keep.size(); k++) if (keep[k] <= keep[k - 1]) monotone = false;
    ck(monotone, "keepSites indices are strictly increasing");

    // The scaffold is bounded by freqData->nloci, not mapData->nloci, exactly
    // as all ten predecessors were.  A shorter FreqData must shorten the scan.
    fd->nloci = 3;
    keep = keepSites(md, fd, NULL);
    ck(keep.size() == 3, "keepSites scans freqData->nloci, not mapData->nloci");
    fd->nloci = N;

    releaseGenMapScaffold(sc);
    releaseMapData(md);
    releaseFreqData(fd);
}

// ---------------------------------------------------------- parseGT() -------
// The whole correctness surface of VCF genotype reading.  Made a named
// function so it can be tested directly rather than inferred from ROH counts
// -- the same move as keepSites() and glToError().
//
// Two conventions are load-bearing and are asserted here rather than left to
// the reader of the code:
//   - a HALF call (0/. or ./1) is MISSING, not a partial dosage.  That is what
//     loadTPEDData does: it accumulates -9 per missing allele and then clamps
//     anything negative to -9.
//   - dosage counts altIndex, and firstCopy says whether the FIRST haplotype
//     carries that allele, matching loadTPEDData's
//     firstCopy[i] = (alleleStr1 == oneAllele).
static void test_parseGT()
{
    int d, pl; bool fc, ph;
    const char *s;

    #define GT(str) (s = (str), parseGT(s, s + strlen(s), 0, 1, d, fc, pl, ph))
    #define GTI(str, gi, ai) (s = (str), parseGT(s, s + strlen(s), (gi), (ai), d, fc, pl, ph))

    // --- unphased biallelic, the ordinary cases ---
    ck(GT("0/0") && d == 0 && pl == 2 && !ph, "parseGT 0/0 -> dosage 0, unphased");
    ck(GT("0/1") && d == 1 && pl == 2 && !ph, "parseGT 0/1 -> dosage 1");
    ck(GT("1/0") && d == 1 && pl == 2,        "parseGT 1/0 -> dosage 1");
    ck(GT("1/1") && d == 2 && pl == 2,        "parseGT 1/1 -> dosage 2");

    // --- phased forms: same dosage, phased flag set, firstCopy meaningful ---
    ck(GT("0|0") && d == 0 && ph, "parseGT 0|0 -> dosage 0, phased");
    ck(GT("0|1") && d == 1 && ph, "parseGT 0|1 -> dosage 1, phased");
    ck(GT("1|0") && d == 1 && ph, "parseGT 1|0 -> dosage 1, phased");
    ck(GT("1|1") && d == 2 && ph, "parseGT 1|1 -> dosage 2, phased");

    // firstCopy distinguishes the two heterozygotes; dosage does not.
    ck(GT("1|0") && fc,  "parseGT 1|0 -> firstCopy true (hap 1 carries ALT)");
    ck(GT("0|1") && !fc, "parseGT 0|1 -> firstCopy false");
    ck(GT("1|1") && fc,  "parseGT 1|1 -> firstCopy true");
    ck(GT("0|0") && !fc, "parseGT 0|0 -> firstCopy false");

    // --- missing and half calls ---
    ck(GT("./.") && d == GENO_MISSING && pl == 2, "parseGT ./. -> missing");
    ck(GT(".|.") && d == GENO_MISSING,            "parseGT .|. -> missing");
    // A HALF call keeps the allele that WAS observed.  These three used to
    // return GENO_MISSING, throwing that allele away; loadTPEDData has always
    // counted it towards the allele frequency while storing the genotype as
    // missing, and the VCF path now does the same.  The genotype is still
    // unusable -- every negative code is missing to lod() -- but the observed
    // allele reaches the frequency.
    ck(GT("0/.") && d == GENO_HALF_OTHER,   "parseGT 0/. -> half call, observed allele is not the counted one");
    ck(GT("./1") && d == GENO_HALF_COUNTED, "parseGT ./1 -> half call, observed allele IS the counted one");
    ck(GT("1/.") && d == GENO_HALF_COUNTED, "parseGT 1/. -> half call, order does not matter");
    ck(GT("1|.") && d == GENO_HALF_COUNTED, "parseGT 1|. -> half call on a phased genotype");
    ck(!genoIsCalled(GENO_HALF_COUNTED) && !genoIsCalled(GENO_HALF_OTHER) && !genoIsCalled(GENO_MISSING),
       "parseGT: no half code counts as a called genotype");
    ck(genoIsCalled(0) && genoIsCalled(1) && genoIsCalled(2),
       "parseGT: 0, 1 and 2 are called genotypes");
    // lod() must not distinguish them: all negatives take its default branch
    ck(lod(GENO_HALF_COUNTED, 0.3, 0.001) == lod(GENO_MISSING, 0.3, 0.001) &&
       lod(GENO_HALF_OTHER,   0.3, 0.001) == lod(GENO_MISSING, 0.3, 0.001),
       "parseGT: a half call scores exactly as missing in lod()");
    // at a multiallelic site the counted allele is whichever altIndex names
    ck(GTI("2/.", 0, 2) && d == GENO_HALF_COUNTED, "parseGT 2/. with altIndex 2 -> counted");
    ck(GTI("1/.", 0, 2) && d == GENO_HALF_OTHER,   "parseGT 1/. with altIndex 2 -> not counted");
    ck(GT(".")   && d == GENO_MISSING && pl == 1, "parseGT . -> missing, ploidy 1");

    // --- ploidy is reported, not rejected: the caller names site and sample ---
    ck(GT("0")     && d == 0 && pl == 1, "parseGT haploid 0 -> dosage 0, ploidy 1");
    ck(GT("1")     && d == 1 && pl == 1, "parseGT haploid 1 -> dosage 1, ploidy 1");
    ck(GT("0/0/1") && pl == 3,           "parseGT polyploid 0/0/1 -> ploidy 3 for the caller to reject");
    ck(GT("0|1|1") && pl == 3 && d == 2, "parseGT polyploid 0|1|1 -> ploidy 3, dosage 2");

    // --- multiallelic: dosage counts altIndex and nothing else ---
    ck(GTI("1/2", 0, 1) && d == 1, "parseGT 1/2 with altIndex 1 -> dosage 1");
    ck(GTI("2/2", 0, 1) && d == 0, "parseGT 2/2 with altIndex 1 -> dosage 0");
    ck(GTI("2/2", 0, 2) && d == 2, "parseGT 2/2 with altIndex 2 -> dosage 2");
    ck(GTI("1/2", 0, 2) && d == 1, "parseGT 1/2 with altIndex 2 -> dosage 1");
    // a two-digit allele index must not be read as two alleles
    ck(GTI("10/10", 0, 10) && d == 2 && pl == 2, "parseGT 10/10 -> one two-digit index, ploidy 2");

    // --- GT inside a colon-separated sample field ---
    ck(GT("0/1:35:99")   && d == 1 && pl == 2, "parseGT GT:DP:GQ suffix ignored");
    ck(GT("1|1:0,30,60") && d == 2 && ph,      "parseGT GT followed by a PL array");
    ck(GTI("35:0/1:99", 1, 1) && d == 1,       "parseGT GT as the SECOND FORMAT field");
    ck(GTI("35:99:1|1", 2, 1) && d == 2 && ph, "parseGT GT as the THIRD FORMAT field");

    // --- input that cannot be interpreted ---
    ck(!GTI("35:99", 2, 1), "parseGT fails when FORMAT promises more fields than exist");
    ck(!GT(""),             "parseGT fails on an empty sample field");
    ck(!GTI("0/1:", 1, 1),  "parseGT fails on an empty GT sub-field");
    ck(!GT("A/T"),          "parseGT fails on a non-digit allele (a TPED-style genotype)");
    ck(!GT("0-1"),          "parseGT fails on an unexpected separator");

    #undef GT
    #undef GTI
}

// -------------------------------------------------------- plToError() ------
// The correct conversion from a VCF PL/GL array, and the reason Commit 4
// exists.  A VCF normalises PL so the CALLED genotype is exactly 0, so the
// single value at the called genotype is 0 for every call, confident or not --
// the information is in how much worse the alternatives are.
//
// The reference values are the GQ the same confidence would be written as:
// error = 10^(-GQ/10).  plToError reproducing them is what establishes that
// the array form is right and the single-value --tgls PL path was not.
static void test_plToError()
{
    vector<double> pl;
    #define PL3(a,b,c) pl.clear(); pl.push_back(a); pl.push_back(b); pl.push_back(c)

    PL3(0,30,60);  ckd(plToError(pl,0), 0.001,     1e-6, "plToError (0,30,60) -> 0.001  == GQ 30");
    PL3(0,10,20);  ckd(plToError(pl,0), 0.0990991, 1e-6, "plToError (0,10,20) -> 0.0991 == GQ 10");
    PL3(0,3,6);    ckd(plToError(pl,0), 0.429346,  1e-5, "plToError (0,3,6)   -> 0.429  == GQ 3.7");
    PL3(0,99,99);  ck(plToError(pl,0) > 2.5e-10 && plToError(pl,0) < 2.6e-10, "plToError (0,99,99) -> 2.5e-10");

    // Which genotype was called changes the answer: the called index selects
    // the numerator.  A single value could not express this.
    PL3(30,0,30);  ckd(plToError(pl,1), 0.00199601,1e-7, "plToError (30,0,30) called=het -> 0.002");
    PL3(60,30,0);  ckd(plToError(pl,2), 0.001,     1e-6, "plToError (60,30,0) called=hom-alt -> 0.001");
    PL3(0,30,60);  ck(plToError(pl,2) > 0.999,                    "plToError (0,30,60) called=hom-alt -> nearly 1");

    // A flat array is no information: 2 of 3 genotypes are wrong.
    PL3(0,0,0);    ckd(plToError(pl,0), 2.0/3.0,   1e-9, "plToError (0,0,0) -> 2/3, no information");

    // A large PL must not underflow the SUM to zero.  Subtracting the minimum
    // before exponentiating is what prevents it; the result hits the 1e-16
    // floor glToError also uses.
    PL3(0,2000,2000); ck(plToError(pl,0) > 0 && plToError(pl,0) <= 1e-16,   "plToError (0,2000,2000) does not underflow to 0");

    // Degenerate inputs return maximum uncertainty rather than garbage.
    pl.clear();      ckd(plToError(pl,0), 1.0, 1e-12, "plToError on an empty array -> 1");
    PL3(0,30,60);    ckd(plToError(pl,-1), 1.0, 1e-12, "plToError with calledIndex -1 -> 1");
    PL3(0,30,60);    ckd(plToError(pl,7), 1.0, 1e-12, "plToError with calledIndex past the end -> 1");

    // Every result is a usable error rate.
    PL3(0,30,60);
    ck(plToError(pl,0) > 0 && plToError(pl,0) <= 1, "plToError stays in (0,1]");

    #undef PL3
}

// ------------------------------------------------ feature counting --------
// Where a classified site falls, for one individual.  The distinction these
// cases protect is NONE against UNASSESSED: the Perl this replaces folded
// them together, so every site on a chromosome the user had excluded from the
// analysis was reported as "not in a run of homozygosity", which inflates the
// outside-ROH count by however much of the genome was left out.

static void test_bucketAt()
{
    ROHIndex::ChrInfo ci;
    ci.lo = 1000; ci.hi = 9000;
    ci.gapLo = 4000; ci.gapHi = 4500;
    ci.role = CHR_AUTOSOME;
    ci.excluded.push_back(Interval(8000, 8100));

    ROHIndex::ChrTracts t;
    t.start.push_back(1500); t.stop.push_back(2000); t.cls.push_back(0);
    t.start.push_back(5000); t.stop.push_back(6000); t.cls.push_back(2);

    // A chromosome that was never analysed.  Not NONE: nothing here is
    // evidence either way.
    ck(ROHIndex::bucketAt(&t, NULL, ZYG_UNKNOWN, 1600) == BUCKET_UNASSESSED,
       "bucketAt: no chromosome info -> UNASSESSED");

    // Outside the span the markers cover, in either direction.
    ck(ROHIndex::bucketAt(&t, &ci, ZYG_UNKNOWN, 999)  == BUCKET_UNASSESSED,
       "bucketAt: before the first analysed marker -> UNASSESSED");
    ck(ROHIndex::bucketAt(&t, &ci, ZYG_UNKNOWN, 9001) == BUCKET_UNASSESSED,
       "bucketAt: past the last analysed marker -> UNASSESSED");
    ck(ROHIndex::bucketAt(&t, &ci, ZYG_UNKNOWN, 1000) != BUCKET_UNASSESSED,
       "bucketAt: the first analysed marker itself is in the span");
    ck(ROHIndex::bucketAt(&t, &ci, ZYG_UNKNOWN, 9000) != BUCKET_UNASSESSED,
       "bucketAt: the last analysed marker itself is in the span");

    // The assembly gap and an excluded region, both inclusive at their ends
    // exactly as inGap and ExcludedRegions are.
    ck(ROHIndex::bucketAt(&t, &ci, ZYG_UNKNOWN, 4000) == BUCKET_UNASSESSED,
       "bucketAt: assembly gap start -> UNASSESSED");
    ck(ROHIndex::bucketAt(&t, &ci, ZYG_UNKNOWN, 4500) == BUCKET_UNASSESSED,
       "bucketAt: assembly gap end -> UNASSESSED");
    ck(ROHIndex::bucketAt(&t, &ci, ZYG_UNKNOWN, 4250) == BUCKET_UNASSESSED,
       "bucketAt: inside the assembly gap -> UNASSESSED");
    ck(ROHIndex::bucketAt(&t, &ci, ZYG_UNKNOWN, 8050) == BUCKET_UNASSESSED,
       "bucketAt: inside an excluded region -> UNASSESSED");
    ck(ROHIndex::bucketAt(&t, &ci, ZYG_UNKNOWN, 8000) == BUCKET_UNASSESSED,
       "bucketAt: excluded region start -> UNASSESSED");
    ck(ROHIndex::bucketAt(&t, &ci, ZYG_UNKNOWN, 8100) == BUCKET_UNASSESSED,
       "bucketAt: excluded region end -> UNASSESSED");

    // Inside a run, at both endpoints and in the middle, returning the class
    // the run was assigned rather than a fixed letter.
    ck(ROHIndex::bucketAt(&t, &ci, ZYG_UNKNOWN, 1500) == 0,
       "bucketAt: first base of a run -> its class");
    ck(ROHIndex::bucketAt(&t, &ci, ZYG_UNKNOWN, 2000) == 0,
       "bucketAt: last base of a run -> its class");
    ck(ROHIndex::bucketAt(&t, &ci, ZYG_UNKNOWN, 1750) == 0,
       "bucketAt: inside a run -> its class");
    ck(ROHIndex::bucketAt(&t, &ci, ZYG_UNKNOWN, 5500) == 2,
       "bucketAt: inside the second run -> its own class");

    // In the callable span but in no run.
    ck(ROHIndex::bucketAt(&t, &ci, ZYG_UNKNOWN, 1400) == BUCKET_NONE,
       "bucketAt: before the first run -> NONE");
    ck(ROHIndex::bucketAt(&t, &ci, ZYG_UNKNOWN, 2001) == BUCKET_NONE,
       "bucketAt: one past a run -> NONE");
    ck(ROHIndex::bucketAt(&t, &ci, ZYG_UNKNOWN, 3000) == BUCKET_NONE,
       "bucketAt: between two runs -> NONE");
    ck(ROHIndex::bucketAt(&t, &ci, ZYG_UNKNOWN, 7000) == BUCKET_NONE,
       "bucketAt: after the last run -> NONE");
    ck(ROHIndex::bucketAt(NULL, &ci, ZYG_UNKNOWN, 3000) == BUCKET_NONE,
       "bucketAt: individual with no runs on this chromosome -> NONE");

    // An individual with one copy of the chromosome.  Distinct from
    // UNASSESSED because the caller must also refuse to read the genotype: a
    // hemizygous call written as a doubled allele looks like a homozygote.
    ROHIndex::ChrInfo sx = ci;
    sx.role = CHR_SEX_SHARED;
    ck(ROHIndex::bucketAt(&t, &sx, ZYG_HETEROGAMETIC, 1750) == BUCKET_INELIGIBLE,
       "bucketAt: heterogametic on the shared sex chromosome -> INELIGIBLE");
    ck(ROHIndex::bucketAt(&t, &sx, ZYG_HOMOGAMETIC, 1750) == 0,
       "bucketAt: homogametic on the shared sex chromosome is called normally");

    // A chromosome nobody is diploid on is ineligible whatever the zygosity.
    ROHIndex::ChrInfo dg = ci;
    dg.role = CHR_SEX_DEGENERATE;
    ck(ROHIndex::bucketAt(&t, &dg, ZYG_HOMOGAMETIC, 1750) == BUCKET_INELIGIBLE,
       "bucketAt: degenerate sex chromosome -> INELIGIBLE");
}

// The index built from one population's calls: individual lookup by ID, the
// per-chromosome resolution the streaming readers use, and the refusal that
// the Perl never made.
static void test_ROHIndex()
{
    vector< MapData * > *mapDataByChr = new vector< MapData * >;
    MapData *m1 = initMapData(3);
    m1->chr = "chr21";
    m1->physicalPos[0] = 1000; m1->physicalPos[1] = 5000; m1->physicalPos[2] = 9000;
    mapDataByChr->push_back(m1);
    // A chromosome carrying one locus is not analysed -- writeFROH leaves it
    // out of the denominator for the same reason -- so its sites must not
    // come back NONE.
    MapData *m2 = initMapData(1);
    m2->chr = "chr22";
    m2->physicalPos[0] = 2000;
    mapDataByChr->push_back(m2);

    IndData *indData = initIndData(2);
    indData->indID[0] = "IND1"; indData->pop[0] = "POPA";
    indData->indID[1] = "IND2"; indData->pop[1] = "POPA";

    vector< ROHData * > *rohDataByInd = new vector< ROHData * >;
    ROHData *r0 = new ROHData;
    r0->indID = "IND1";
    // Deliberately out of order: bucketAt binary-searches, so addPopulation
    // has to sort rather than trust the order it was handed.
    r0->chr.push_back(0); r0->start.push_back(6000); r0->stop.push_back(7000); r0->length.push_back(1001);
    r0->chr.push_back(0); r0->start.push_back(1500); r0->stop.push_back(2000); r0->length.push_back(501);
    rohDataByInd->push_back(r0);
    ROHData *r1 = new ROHData;
    r1->indID = "IND2";
    rohDataByInd->push_back(r1);

    vector<double> bounds;
    bounds.push_back(1000.0);     // < 1000 is class A, else class B

    ROHIndex idx;
    ck(idx.addPopulation("POPA", rohDataByInd, mapDataByChr, bounds,
                         indData, NULL, NULL, NULL) == 0,
       "ROHIndex: a population is added");
    ck(idx.nind() == 2, "ROHIndex: both individuals are in the index");
    ck(idx.npop() == 1, "ROHIndex: one population");
    ck(idx.nclass(0) == 2, "ROHIndex: one boundary gives two size classes");
    ck(idx.indexOf("IND1") == 0 && idx.indexOf("IND2") == 1,
       "ROHIndex: individuals are found by ID");
    // The check the Perl never made: an ID that is not in the call set.
    ck(idx.indexOf("NOT_A_SAMPLE") == -1,
       "ROHIndex: an unknown ID is reported, not silently treated as callable");

    vector<const ROHIndex::ChrTracts *> tracts;
    vector<const ROHIndex::ChrInfo *> info;
    // Resolved by canonical key, so the counting file may spell the
    // chromosome differently from the calling file.
    idx.resolveChr("21", tracts, info);
    ck(info[0] != NULL && info[1] != NULL, "ROHIndex: chr21 resolves from the key 21");
    ck(tracts[0] != NULL && tracts[0]->start.size() == 2,
       "ROHIndex: the individual's two runs are indexed");
    ck(tracts[0]->start[0] == 1500 && tracts[0]->start[1] == 6000,
       "ROHIndex: runs are sorted by start");
    ck(tracts[1] == NULL, "ROHIndex: an individual with no runs has no tracts");
    ck(info[0]->lo == 1000 && info[0]->hi == 9000,
       "ROHIndex: the analysed span is the first and last marker");

    ck(ROHIndex::bucketAt(tracts[0], info[0], idx.zygoOf(0), 1750) == 0,
       "ROHIndex: a 501 bp run is class A");
    ck(ROHIndex::bucketAt(tracts[0], info[0], idx.zygoOf(0), 6500) == 1,
       "ROHIndex: a 1001 bp run is class B");
    ck(ROHIndex::bucketAt(tracts[1], info[1], idx.zygoOf(1), 6500) == BUCKET_NONE,
       "ROHIndex: an individual with no runs is NONE inside the span");

    idx.resolveChr("22", tracts, info);
    ck(info[0] == NULL, "ROHIndex: a chromosome with one locus is not analysed");
    idx.resolveChr("7", tracts, info);
    ck(info[0] == NULL, "ROHIndex: a chromosome absent from the run is not analysed");

    // The same ID in two populations would make one population's counts
    // overwrite the other's.  Nothing else in garlic looks for it: duplicate
    // IDs are refused within a file, not across the population split.
    quietErrors(true);
    ck(idx.addPopulation("POPB", rohDataByInd, mapDataByChr, bounds,
                         indData, NULL, NULL, NULL) == -1,
       "ROHIndex: an individual in two populations is refused");
    quietErrors(false);

    delete r0; delete r1; delete rohDataByInd;
    releaseIndData(indData);
    releaseMapData(mapDataByChr);
}

// Inverse of sizeClassLabel, which is how a class letter read back out of a
// .roh.bed becomes a bucket again.  Round-tripping is the property that
// matters: --nclust is unbounded, so the letters go past Z.
static void test_sizeClassIndexFromLabel()
{
    for (int k = 0; k < 60; k++)
        ck(ROHIndex::sizeClassIndexFromLabel(sizeClassLabel(k)) == k,
           "sizeClassIndexFromLabel round-trips sizeClassLabel");
    ck(ROHIndex::sizeClassIndexFromLabel("A")  == 0,  "label A is class 0");
    ck(ROHIndex::sizeClassIndexFromLabel("Z")  == 25, "label Z is class 25");
    ck(ROHIndex::sizeClassIndexFromLabel("AA") == 26, "label AA is class 26");
    ck(ROHIndex::sizeClassIndexFromLabel("")   == -1, "an empty label is not a class");
    ck(ROHIndex::sizeClassIndexFromLabel("a")  == -1, "a lower-case label is not one garlic writes");
    ck(ROHIndex::sizeClassIndexFromLabel("A1") == -1, "a digit is not a class label");
    // Any run of capitals is arithmetically a label, so the guard is on size
    // rather than on spelling: garlic never writes one this long, and a file
    // that carries one is not a .roh.bed.
    ck(ROHIndex::sizeClassIndexFromLabel("ABCDEFGH") == -1,
       "an absurdly long label is refused rather than becoming a huge class index");
}

// The feature file.  Return codes only: the wording of each refusal is
// asserted in test/run_tests.sh, where it is what the user actually sees.
static string tmpPath(const char *name)
{
    const char *d = getenv("TMPDIR");
    if (d == NULL || d[0] == '\0') d = getenv("TMP");
    if (d == NULL || d[0] == '\0') d = "/tmp";
    string p(d);
    if (!p.empty() && p[p.size() - 1] != '/' && p[p.size() - 1] != '\\') p += "/";
    return p + name;
}

static bool writeLines(const string &path, const char *text)
{
    ofstream f(path.c_str());
    if (f.fail()) return false;
    f << text;
    f.close();
    return true;
}

static void test_FeatureTable()
{
    const string good = tmpPath("garlic_ut_features_good.txt");
    // Comments, a blank line, mixed chromosome spellings, a lower-case
    // allele, and one site carrying two classes.
    if (!writeLines(good,
        "# chr pos allele class\n"
        "chr21\t100\tT\tnonsynonymous\n"
        "\n"
        "chr21 100 T probably_damaging\n"
        "21 200 a benign\n"
        "  # indented comment\n"
        "Chr21 300 G benign\n"))
    {
        ck(false, "FeatureTable: could not write the test file");
        return;
    }

    FeatureTable ft;
    ck(ft.read(good) == 0, "FeatureTable: a well-formed file reads");
    ck(ft.nrows() == 4, "FeatureTable: four classified rows");
    ck(ft.nsites() == 3, "FeatureTable: three distinct sites");
    ck(ft.nclass() == 3, "FeatureTable: three distinct classes");
    // Sorted, so the same annotation scheme gives the same columns whatever
    // order the rows were written in.
    ck(ft.classes()[0] == "benign" && ft.classes()[1] == "nonsynonymous" &&
       ft.classes()[2] == "probably_damaging",
       "FeatureTable: classes come back in sorted order");

    // All three spellings landed on one chromosome.
    const map<pos_t, vector<FeatureAllele> > *c = ft.chromosome("21");
    ck(c != NULL, "FeatureTable: chr21, 21 and Chr21 are one chromosome");
    if (c != NULL)
    {
        ck(c->size() == 3, "FeatureTable: three sites on that chromosome");
        map<pos_t, vector<FeatureAllele> >::const_iterator it = c->find(100);
        ck(it != c->end() && it->second.size() == 2,
           "FeatureTable: a site may carry two classes");
        it = c->find(200);
        ck(it != c->end() && it->second.size() == 1 && it->second[0].allele == "A",
           "FeatureTable: alleles are upper-cased for comparison");
    }
    ck(ft.chromosome("22") == NULL, "FeatureTable: an absent chromosome is NULL");
    remove(good.c_str());

    quietErrors(true);

    const string dup = tmpPath("garlic_ut_features_dup.txt");
    writeLines(dup, "chr21 100 T benign\nchr21 100 t benign\n");
    FeatureTable f2;
    ck(f2.read(dup) == -1, "FeatureTable: the same class twice at one site is refused");
    remove(dup.c_str());

    // The format the Perl scripts took.  Recognised so the message can say so.
    const string old = tmpPath("garlic_ut_features_old.txt");
    writeLines(old, "chr21:100 C T probably_damaging\n");
    FeatureTable f3;
    ck(f3.read(old) == -1, "FeatureTable: the old chr:pos format is refused");
    remove(old.c_str());

    const string badpos = tmpPath("garlic_ut_features_badpos.txt");
    writeLines(badpos, "chr21 1e5 T benign\n");
    FeatureTable f4;
    ck(f4.read(badpos) == -1, "FeatureTable: an unparsable position is refused");
    remove(badpos.c_str());

    const string zeropos = tmpPath("garlic_ut_features_zeropos.txt");
    writeLines(zeropos, "chr21 0 T benign\n");
    FeatureTable f5;
    ck(f5.read(zeropos) == -1, "FeatureTable: position 0 is refused (positions are 1-based)");
    remove(zeropos.c_str());

    const string ncol = tmpPath("garlic_ut_features_ncol.txt");
    writeLines(ncol, "chr21 100 T benign extra\n");
    FeatureTable f6;
    ck(f6.read(ncol) == -1, "FeatureTable: five columns are refused");
    remove(ncol.c_str());

    const string blank = tmpPath("garlic_ut_features_blank.txt");
    writeLines(blank, "# nothing but a comment\n\n");
    FeatureTable f7;
    ck(f7.read(blank) == -1, "FeatureTable: a file with no classified sites is refused");
    remove(blank.c_str());

    FeatureTable f8;
    ck(f8.read(tmpPath("garlic_ut_features_absent.txt")) == -1,
       "FeatureTable: a missing file is refused");

    quietErrors(false);
}

// --------------------------------------------- logger overload coverage ----
// A COMPILE-TIME check, never called.  int64_t is long long under LLP64
// (macOS arm64, Windows) and long under LP64 (Linux, most BSDs), so an errlog
// that overloads on only one of the two 64-bit types gives an exact match on
// one ABI and none on the other -- where the argument then converts equally
// well to int, long long, double, bool and char, and the call is ambiguous.
//
// That is exactly how the first CI run failed: garlic-data.cpp:1353 logs a
// pos_t, which compiled on macOS and was ambiguous on linux-x86_64.  Nothing
// in a macOS build can notice, so the check has to be structural rather than
// a test that runs.  If an overload goes missing this file stops compiling,
// and the unit-test stage runs on every platform in CI.
static void compile_time_logger_overload_coverage(errlog &L)
{
    int       i  = 0;
    long      l  = 0;
    long long ll = 0;
    int64_t   i64 = 0;
    pos_t     p  = 0;
    double    d  = 0;
    char      c  = 'x';
    bool      b  = false;
    std::string s = "s";

    L.err("x", i, false);   L.log("x", i, false);
    L.err("x", l, false);   L.log("x", l, false);
    L.err("x", ll, false);  L.log("x", ll, false);
    L.err("x", i64, false); L.log("x", i64, false);
    L.err("x", p, false);   L.log("x", p, false);
    L.err("x", d, false);   L.log("x", d, false);
    L.err("x", c, false);   L.log("x", c, false);
    L.err("x", b, false);   L.log("x", b, false);
    L.err("x", s, false);   L.log("x", s, false);
}

int main()
{
    //Referenced, never invoked: taking the address requires the function to
    //compile (which is the whole point) while keeping it out of the run.
    (void)&compile_time_logger_overload_coverage;
    printf("garlic unit tests\n");
    test_enumeratePopulations();
    test_subsetDataByIndex();
    test_alleleFrequency();
    test_addAlleleCounts();
    test_calcFreqDataForIndices();
    test_lod();
    test_interpolate();
    test_inGap();
    test_auto_fits();
    test_class_labels();
    test_chr_names();
    test_canonChrKey();
    test_sex_model();
    test_winsize_for_chr();
    test_ld_no_shared_genotypes();
    test_chr_lengths();
    test_chr_key_collisions();
    test_excluded_regions();
    test_glToError();
    test_plToError();
    test_kde_helpers();
    test_getMapInfo();
    test_keepSites();
    test_parseGT();
    test_bucketAt();
    test_ROHIndex();
    test_FeatureTable();
    test_sizeClassIndexFromLabel();
    printf("%d checks, %d failures\n", checks, failures);
    return failures == 0 ? 0 : 1;
}
