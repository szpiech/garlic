#ifndef __GARLIC_ROH_H__
#define __GARLIC_ROH_H__
#include "garlic-cli.h"
#include "garlic-data.h"
#include "garlic-centromeres.h"
#include "param_t.h"
//#include <pthread.h>
#include <cmath>
#include <iostream>
#include "garlic-errlog.h"
#include "garlic-kde.h"
#include "gmm.h"
#include "gsl/gsl_statistics.h"
#include "gsl/gsl_sort.h"
#include "BoundFinder.h"
#include <limits>

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
};

struct ROHData
{
  string indID;
  vector<int> chr;
  vector<int> start;
  vector<int> stop;
};

struct ROHLength
{
  string pop;
  double *length;
  int size;
};

//extern double **LD;
//extern double **LD1;

void parallelwLOD(void *order);

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

double lod(const short &genotype, const double &freq, const double &error);

KDEResult *automaticallyChooseWindowSize(vector< HapData * > *hapDataByChr, vector< FreqData * > *freqDataByChr,
    vector< MapData * > *mapDataByChr, IndData *indData,
    centromere *centro, int &winsize, double error, int MAX_GAP,
    int KDE_SUBSAMPLE, int numThreads, bool WINSIZE_EXPLORE, double AUTO_WINSIZE_THRESHOLD, string outfile);

KDEWinsizeReport *calculateLODOverWinsizeRange(vector< HapData * > *hapDataByChr, vector< FreqData * > *freqDataByChr,
    vector< MapData * > *mapDataByChr, IndData *indData,
    centromere *centro, vector<int> *multiWinsizes, double error, int MAX_GAP,
    int KDE_SUBSAMPLE, int numThreads, bool WINSIZE_EXPLORE, string outfile);

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
                                        double OVERLAP_FRAC);

ROHLength *initROHLength(int size, string pop);
void releaseROHLength(ROHLength *rohLength);

vector< ROHData * > *initROHData(IndData *indData);
void writeROHData(string outfile,
                  vector< ROHData * > *rohDataByInd,
                  vector< MapData * > *mapDataByChr,
                  int_pair_t bounds,
                  string popName,
                  string version);
void releaseROHData(vector< ROHData * > *rohDataByInd);

string makeROHFilename(string outfile);

double selectLODCutoff(KDEResult *kdeResult);
double selectLODCutoff(vector< WinData * > *winDataByChr, IndData *indData, int KDE_SUBSAMPLE, string kdeoutfile);

void exploreWinsizes(vector< HapData * > *hapDataByChr,
                     vector< FreqData * > *freqDataByChr,
                     vector< MapData * > *mapDataByChr,
                     IndData *indData,
                     centromere *centro,
                     vector<int> &multiWinsizes,
                     double error,
                     vector< GenoLikeData * > *GLDataByChr, bool USE_GL,
                     int MAX_GAP, int KDE_SUBSAMPLE, string outfile,
                     bool WEIGHTED, vector< GenoFreqData * > *genoFreqDataByChr, bool PHASED);

KDEResult *selectWinsizeFromList(vector< HapData * > *hapDataByChr,
                                 vector< FreqData * > *freqDataByChr,
                                 vector< MapData * > *mapDataByChr,
                                 IndData *indData, centromere *centro,
                                 vector<int> *multiWinsizes, int &winsize, double error,
                                 vector< GenoLikeData * > *GLDataByChr, bool USE_GL,
                                 int MAX_GAP, int KDE_SUBSAMPLE, string outfile,
                                 bool WEIGHTED, vector< GenoFreqData * > *genoFreqDataByChr, bool PHASED);

KDEResult *selectWinsize(vector< HapData * > *hapDataByChr,
                         vector< FreqData * > *freqDataByChr,
                         vector< MapData * > *mapDataByChr,
                         IndData *indData, centromere *centro,
                         int &winsize, int step, double error,
                         vector< GenoLikeData * > *GLDataByChr, bool USE_GL,
                         int MAX_GAP, int KDE_SUBSAMPLE, string outfile,
                         bool WEIGHTED, vector< GenoFreqData * > *genoFreqDataByChr, bool PHASED);

//int selectWinsize(KDEWinsizeReport *winsizeReport, double AUTO_WINSIZE_THRESHOLD);

vector<int> *getWinsizeList(int lastWinsize, int stepSize, int numThreads);

bool inGap(int qStart, int qEnd, int targetStart, int targetEnd);

int_pair_t selectSizeClasses(ROHLength *rohLength);

void compute(void *order);

//extern pthread_mutex_t io_mutex;

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
