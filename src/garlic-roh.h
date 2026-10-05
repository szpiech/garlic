#ifndef __GARLIC_ROH_H__
#define __GARLIC_ROH_H__
#include "garlic-cli.h"
#include "garlic-data.h"
#include "garlic-centromeres.h"
#include "param_t.h"
#include <thread>
#include <cmath>
#include <iostream>
#include "garlic-errlog.h"
#include "garlic-kde.h"
#include "gmm.h"
#include "garlic-math.h"
#include "BoundFinder.h"
#include <limits>
#include "garlic-pbar.h"

using namespace std;

struct WLOD_work_order_t
{
  MapData *mapData;
  HapData *hapData;
  FreqData *freqData;
  GenoLikeData *GLData;
  LDData *LD;
  WinData *winData;
  int cStart;
  int cEnd;
  int winsize;
  double error;
  int MAX_GAP;
  bool USE_GL;
  double mu;
  int M;
  int start;
  int stop;
  Bar *bar;
  int numThreads;
};

struct ROHData
{
  string indID;
  vector<int> chr;
  //Physical positions in both bp and --cm mode; only length becomes a genetic
  //distance under --cm, so it stays a double.
  vector<pos_t> start;
  vector<pos_t> stop;
  vector<double> length;
};

struct ROHLength
{
  //string pop;
  vector<double> length;
  double size;
};

//extern double **LD;
//extern double **LD1;

int selectWinsizeWeighted(double density);

void parallelwLOD(WLOD_work_order_t *p);

void setLODThreads(int n);
//Both were hardcoded: the smoothness threshold at which the --auto-winsize
//search stops (0.50, duplicated in two functions), and the coefficients of the
//empirical winsize-vs-density fit used with --weighted.
void setAutoWinsizeThreshold(double t);
void setAutoWinsizeCoef(double slope, double intercept);
void setGMMParams(int maxIter, double tol);

//Per-individual autozygous total and fraction, by size class.  Long format so
//the number of classes (which --nclust controls) does not change the columns.
//popLabel is empty for a single-population run and pooled marks
//--pool-populations, so a table can be told apart from one produced by a file
//that really had a single population.
//Per-individual information report, <out>.ind.tsv.
//
//What a window's LOD means depends on how many of its sites carried
//information for THIS individual: a missing genotype (or a flat likelihood
//triplet) contributes exactly 0, so two individuals with the same window size
//can be summing very different amounts of evidence.  This records, per
//individual: the fraction of sites with information, the median and 5th
//percentile of informative sites per window, the observed/expected
//heterozygote ratio over the sites in windows judged non-autozygous (window
//LOD < 0, a cutoff-free rule so this can run before the genotype matrices are
//released), and, for likelihood input only, the allelic-dropout estimate
//d-hat: the maximiser of the profile marginal likelihood over those sites
//under a confusion model where a true heterozygote is observed as a
//homozygote with probability d.  d-hat is 0 for an unbiased individual and
//rises with reference bias; it is a diagnostic, and nothing corrects for it.
//The het ratio from called genotypes is NOT a low-depth estimator -- under
//low depth it reads a large deficit with no bias at all -- which the manual
//says; it is the ordinary sample-quality statistic for called data.
struct IndReportRow
{
    long long nSites;        //analysed sites (autosomes + shared sex chromosome where diploid)
    long long nInformative;  //called genotype / non-flat triplet
    double    snpsPerWindowMedian;
    double    snpsPerWindowP05;
    double    hetExpected;   //sum of 2p(1-p) over the non-autozygous sites used: the sample size behind the two statistics below
    double    hetObsExp;     //NaN when no expected hets
    double    dropoutHat;    //NaN for non-likelihood input
};

//Stratified LOD cutoff (--stratify-cutoff).
//
//The automatic cutoff is the valley between the two modes of the window-LOD
//density POOLED over individuals.  That is self-calibrating when every
//individual carries the same amount of evidence per window, and wrong when
//they do not: an individual at low depth has both modes compressed toward
//zero, so a cutoff placed for the cohort calls most of that individual's
//genome.  Two things were tried and rejected before this: standardising each
//individual's windows by (mode, spread) and pooling -- the non-autozygous
//spread is set by allele-frequency heterogeneity across windows while the
//autozygous mode is set by information per site, so the two do not scale
//together and the standardised density has no common autozygous mode -- and
//a valley per individual, which is too noisy on one individual's windows
//(10-26% more genome called on the bundled uniform-depth example, and worse
//with a shared bandwidth).
//
//What works is to notice that the location of an individual's
//non-autozygous mode is a precise statistic of how much evidence its windows
//carry (sd 1.3 LOD across the bundled cohort, against 18 and 29 LOD between
//2x, 4x and 30x), and to apply the POOLED rule within strata of individuals
//whose modes agree: sorted by mode, two neighbours join a stratum when their
//modes differ by less than a quarter of the smaller spread.  A uniform cohort
//is one stratum and reproduces the pooled cutoff exactly; a mixed cohort
//splits by depth and each stratum gets the valley of its own pooled density,
//with tens of thousands of windows behind it.  A stratum whose density has
//no valley borrows the others' median standardised valley through its own
//median mode and spread.
struct IndStandardization
{
    double mode;        //location of the individual's non-autozygous mode, LOD units
    double spread;      //1.4826 * median |w - mode| over windows <= mode
    double cutoff;      //the cutoff applied to this individual
    int    stratum;     //0-based, in order of increasing mode; -1 = too few windows
    int    stratumSize; //individuals in it
    bool   own;         //the stratum's cutoff is its own valley (else borrowed)
    bool   small;       //stratum too small for a valley; took the largest stratum's cutoff
    bool   hasScale;    //mode and spread were estimated
    bool   unresolved;  //enough windows, but no density peak at LOD < 0: the non-autozygous mode could not be located
    IndStandardization() : mode(0), spread(1), cutoff(0), stratum(-1), stratumSize(0), own(false), small(false), hasScale(false), unresolved(false) {}
};

//How many evidence strata the cohort falls into, with each stratum's median
//non-autozygous mode and size -- for the pooled selector to warn when it is
//about to place one cutoff for individuals that do not carry the same
//evidence per window.
int countEvidenceStrata(vector< WinData * > *winDataByChr, IndData *indData, int step,
                        const vector<ChrRole> *role, ChrRole keep, vector<double> &stratumModes,
                        vector<int> &stratumSizes);

//Returns the number of strata; perInd receives one entry per individual with
//the cutoff actually applied.  The KDE written to kdeoutfile is of the
//largest stratum's pooled windows.
//`noAutozygousMode` is set when no stratum's density has an autozygous mode:
//every individual then receives an infinite cutoff (no calls) and ok stays
//true, as for the pooled selector.
int selectLODCutoffStratified(vector< WinData * > *winDataByChr, IndData *indData,
                              string kdeoutfile, int step, int wsize, bool &ok,
                              const vector<ChrRole> *role, ChrRole keep,
                              vector<IndStandardization> &perInd, bool *noAutozygousMode = NULL);

vector<IndReportRow> computeIndReport(vector< HapData * > *hapDataByChr,
                                      vector< FreqData * > *freqDataByChr,
                                      vector< GenoLikeData * > *GLDataByChr,
                                      vector< WinData * > *winDataByChr,
                                      IndData *indData,
                                      int winsize, double error, bool USE_GL,
                                      const vector<ChrRole> *role);

//perInd may be NULL; when given, more columns record the stratum and the
//cutoff each individual was called with.
void writeIndReport(string outfile, const vector<IndReportRow> &rows, IndData *indData,
                    int winsize, double lodCutoff, bool likelihoodInput,
                    const vector<IndStandardization> *perInd = NULL);

void writeFROH(string outfile,
               vector< ROHData * > *rohDataByInd,
               vector< MapData * > *mapDataByChr,
               vector< double > bounds,
               IndData *indData,
               centromere *centro,
               bool CM,
               const string &popLabel,
               bool pooled,
               const vector<ChrRole> *role = NULL,
               //Excluded regions are already gone from the map, so a TERMINAL
               //one needs no handling: dropping its loci moved the first or
               //last marker inward and the span shrank by itself.  An INTERIOR
               //one does not, and stays in the denominator unless it is
               //subtracted here.
               const ExcludedRegions *excluded = NULL,
               //One of the FROHDenominator values, and the lengths it needs.
               //chrLengths may be NULL under FROH_ANALYZED, which reads none.
               int denomMode = FROH_ANALYZED,
               const ChrLengths *chrLengths = NULL);
void calcLOD(MapData *mapData,
             HapData *hapData, FreqData *freqData,
             GenoLikeData *GLData,
             WinData *winData, centromere *centro,
             int winsize, double error, int MAX_GAP, bool USE_GL);

void calcwLOD(MapData *mapData,
              HapData *hapData, FreqData *freqData,
              GenoLikeData *GLData,
              LDData *LD,
              WinData *winData, centromere *centro,
              int winsize, double error, int MAX_GAP, bool USE_GL,
              double mu, int M, int numThreads);

double nomut(double M, double mu, double interval);
double norec(double M, double interval);

double lod(const int genotype, const double &freq, const double &error);

//The same emission marginalised over a genotype-likelihood triplet L[0..2]
//(copies of the counted allele 0, 1, 2; any scale):
//    log10( sum_g L[g] A(g) / sum_g L[g] N(g) )
//with A and N exactly the two mixtures lod() forms.  A degenerate triplet
//reproduces lod(g) bit for bit; the flat triplet gives exactly 0, like a
//missing genotype; and both sums lie in (0, 1], so no value needs clamping.
//`error` is the run-wide --error: the autozygous state's allowance for a
//heterozygote, not a per-genotype quality.
double lodGL(const float *L, const double &freq, const double &error);

/*
KDEResult *automaticallyChooseWindowSize(vector< HapData * > *hapDataByChr, vector< FreqData * > *freqDataByChr,
    vector< MapData * > *mapDataByChr, IndData *indData,
    centromere *centro, int &winsize, double error, int MAX_GAP,
    int KDE_SUBSAMPLE, int numThreads, bool WINSIZE_EXPLORE, double AUTO_WINSIZE_THRESHOLD, string outfile);

KDEWinsizeReport *calculateLODOverWinsizeRange(vector< HapData * > *hapDataByChr, vector< FreqData * > *freqDataByChr,
    vector< MapData * > *mapDataByChr, IndData *indData,
    centromere *centro, vector<int> *multiWinsizes, double error, int MAX_GAP,
    int KDE_SUBSAMPLE, int numThreads, bool WINSIZE_EXPLORE, string outfile);
*/
vector< WinData * > *calcLODWindows(vector< HapData * > *hapDataByChr,
                                    vector< FreqData * > *freqDataByChr,
                                    vector< MapData * > *mapDataByChr,
                                    vector< GenoLikeData * > *GLDataByChr,
                                    centromere *centro,
                                    int winsize, double error,
                                    int MAX_GAP, bool USE_GL,
                                    int sexWinsize = 0,
                                    const vector<ChrRole> *role = NULL);

vector< WinData * > *calcwLODWindows(vector< HapData * > *hapDataByChr,
                                     vector< FreqData * > *freqDataByChr,
                                     vector< MapData * > *mapDataByChr,
                                     vector< GenoLikeData * > *GLDataByChr,
                                     vector< LDData * > *ldDataByChr,
                                     centromere *centro,
                                     int winsize, double error,
                                     int MAX_GAP, bool USE_GL, 
                                     int M, double mu, int numThreads,
                                     int sexWinsize = 0,
                                     const vector<ChrRole> *role = NULL);

vector< ROHData * > *assembleROHWindows(vector< WinData * > *winDataByChr,
                                        vector< MapData * > *mapDataByChr,
                                        IndData *indData,
                                        centromere *centro,
                                        double lodScoreCutoff,
                                        ROHLength **rohLength,
                                        int winSize,
                                        int MAX_GAP,
                                        double OVERLAP_FRAC,
                                        bool CM,
                                        //Tracts from these chromosomes, and only
                                        //these, feed the size-class GMM.  Every
                                        //tract is still written; what a 5 Mb run
                                        //means differs between an autosome and a
                                        //sex chromosome, and one chromosome's
                                        //worth of tracts cannot support a
                                        //three-component fit anyway.
                                        const vector<ChrRole> *role = NULL,
                                        //The shared sex chromosome's own
                                        //cutoff and window size, when it was
                                        //given one.  A cutoff belongs to the
                                        //window size it was estimated at, so
                                        //these travel together.
                                        double sexCutoff = 0, bool haveSexCutoff = false,
                                        int sexWinsize = 0,
                                        //--stratify-cutoff: one cutoff per
                                        //individual, replacing lodScoreCutoff
                                        //wherever the autosomal cutoff applies.
                                        const vector<double> *indCutoff = NULL);

ROHLength *initROHLength(int size);
void releaseROHLength(ROHLength *rohLength);

vector< ROHData * > *initROHData(IndData *indData);
void writeROHData(string outfile,
                  vector< ROHData * > *rohDataByInd,
                  vector< MapData * > *mapDataByChr,
                  vector< double > bounds,
                  const vector<string> &pop,
                  string version,
                  bool CM);
void releaseROHData(vector< ROHData * > *rohDataByInd);

string makeROHFilename(string outfile);
string sizeClassLabel(int k);
vector<string> makeClassColors(int nclass);
//Which size class a run of this length falls in, given the boundaries in
//force.  Shared with the feature counter so a classified genotype is bucketed
//by exactly the class letter the .roh.bed carries for the same run.
int rohSizeClassIndex(double size, const vector<double> &bounds);

//Automatic cutoff from the pooled LOD density.  `modes` receives what the
//density had (see identifyLODModes): on a bimodal density the return is the
//valley; on a unimodal-null density the return is +infinity (no window
//passes) and ok stays true -- the caller decides how to report it; on a
//density with no non-autozygous mode ok is false.
double selectLODCutoff(KDEResult *kdeResult, int wisize, bool &ok, LODModes *modes = NULL);
double selectLODCutoff(vector< WinData * > *winDataByChr, IndData *indData, int KDE_SUBSAMPLE, string kdeoutfile, int step, int wisize, bool &ok,
                       const vector<ChrRole> *role = NULL, ChrRole keep = CHR_AUTOSOME, LODModes *modes = NULL);

//The cutoff the SEX CHROMOSOME's own windows would have given, logged and not
//used.  Cotter et al. (2024) reused the autosomal cutoff on the X for want of
//enough X tracts to locate a second minimum, and noted that doing so may
//inflate X ROH; this reports how far apart the two are for the data in hand.
//Returns false, quietly, when there is no minimum between modes to find --
//which is the normal outcome on one chromosome and is not an error.
bool reportSexChrLODCutoff(vector< WinData * > *winDataByChr, IndData *indData,
                           const vector<ChrRole> *role, int step, int wsize,
                           double appliedCutoff, bool cutoffIsItsOwn = false);

//Sets MISSING on every window of an individual that cannot carry a run of
//homozygosity on that chromosome.  convertWinData2DoubleData already skips
//MISSING and parallelAssembleROH's `>= cutoff` test already fails it, so this
//removes those windows from both the density and the calls without either
//having to know why.  Necessary rather than incidental: a hemizygous window
//scores exactly 0, not MISSING, because lod() returns log10(1/1) for a
//genotype outside {0,1,2} -- so left alone it would pile a spike at zero into
//the density and would be CALLED outright under any negative cutoff.
long long maskIneligibleWindows(vector< WinData * > *winDataByChr, IndData *indData,
                                const vector<ChrRole> *role);

void exploreWinsizes(vector< HapData * > *hapDataByChr,
                     vector< FreqData * > *freqDataByChr,
                     vector< MapData * > *mapDataByChr,
                     IndData *indData,
                     centromere *centro,
                     vector<int> &multiWinsizes,
                     double error,
                     vector< GenoLikeData * > *GLDataByChr,
                     vector< GenoFreqData * > *genoFreqDataByChr, bool USE_GL,
                     int MAX_GAP, int KDE_SUBSAMPLE, string outfile,
                     bool WEIGHTED, int M, double mu, int numThreads, bool PHASED, int thinStep, int LD_SUBSAMPLE,
                     const vector<ChrRole> *role = NULL);

KDEResult *selectWinsizeFromList(vector< HapData * > *hapDataByChr,
                                 vector< FreqData * > *freqDataByChr,
                                 vector< MapData * > *mapDataByChr,
                                 IndData *indData, centromere *centro,
                                 vector<int> *multiWinsizes, int &winsize, double error,
                                 vector< GenoLikeData * > *GLDataByChr, bool USE_GL,
                                 int MAX_GAP, int KDE_SUBSAMPLE, string outfile,
                                 bool WEIGHTED, vector< GenoFreqData * > *genoFreqDataByChr, bool PHASED, int thinStep,
                                 const vector<ChrRole> *role = NULL);

KDEResult *selectWinsize(vector< HapData * > *hapDataByChr,
                         vector< FreqData * > *freqDataByChr,
                         vector< MapData * > *mapDataByChr,
                         IndData *indData, centromere *centro,
                         int &winsize, int step, double error,
                         vector< GenoLikeData * > *GLDataByChr, bool USE_GL,
                         int MAX_GAP, int KDE_SUBSAMPLE, string outfile,
                         bool WEIGHTED, vector< GenoFreqData * > *genoFreqDataByChr, bool PHASED, int thinStep,
                         int MAX_WINSIZE,
                         const vector<ChrRole> *role = NULL);

//int selectWinsize(KDEWinsizeReport *winsizeReport, double AUTO_WINSIZE_THRESHOLD);

vector<int> *getWinsizeList(int lastWinsize, int stepSize, int numThreads);

bool inGap(pos_t qStart, pos_t qEnd, pos_t targetStart, pos_t targetEnd);

//int_pair_t selectSizeClasses(ROHLength *rohLength);
vector<double> selectSizeClasses(ROHLength *rohLength, int NCLUST);





/*
vector< vector< WinData * >* > *calcLODWindowsSinglePop(vector< vector< HapData * >* > *hapDataByPopByChr,
    vector< vector< FreqData * >* > *freqDataByPopByChr,
    vector< MapData * > *mapDataByChr,
    vector< IndData * > *indDataByPop, centromere *centro,
    int* winsize, double error, int MAX_GAP, int numThreads, int pop);
*/
//void scan(void *work_order);
//void scanSinglePop(void *order);

#endif
