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
void writeFROH(string outfile,
               vector< ROHData * > *rohDataByInd,
               vector< MapData * > *mapDataByChr,
               vector< double > bounds,
               const vector<string> &pop,
               centromere *centro,
               bool CM,
               const string &popLabel,
               bool pooled);
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
                                    int MAX_GAP, bool USE_GL);

vector< WinData * > *calcwLODWindows(vector< HapData * > *hapDataByChr,
                                     vector< FreqData * > *freqDataByChr,
                                     vector< MapData * > *mapDataByChr,
                                     vector< GenoLikeData * > *GLDataByChr,
                                     vector< LDData * > *ldDataByChr,
                                     centromere *centro,
                                     int winsize, double error,
                                     int MAX_GAP, bool USE_GL, 
                                     int M, double mu, int numThreads);

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
                                        const vector<ChrRole> *role = NULL);

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

double selectLODCutoff(KDEResult *kdeResult, int wisize, bool &ok);
double selectLODCutoff(vector< WinData * > *winDataByChr, IndData *indData, int KDE_SUBSAMPLE, string kdeoutfile, int step, int wisize, bool &ok,
                       const vector<ChrRole> *role = NULL, ChrRole keep = CHR_AUTOSOME);

//The cutoff the SEX CHROMOSOME's own windows would have given, logged and not
//used.  Cotter et al. (2024) reused the autosomal cutoff on the X for want of
//enough X tracts to locate a second minimum, and noted that doing so may
//inflate X ROH; this reports how far apart the two are for the data in hand.
//Returns false, quietly, when there is no minimum between modes to find --
//which is the normal outcome on one chromosome and is not an error.
bool reportSexChrLODCutoff(vector< WinData * > *winDataByChr, IndData *indData,
                           const vector<ChrRole> *role, int step, int wsize,
                           double autosomalCutoff);

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
