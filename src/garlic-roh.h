#ifndef __GARLIC_ROH_H__
#define __GARLIC_ROH_H__
#include "garlic-cli.h"
#include "garlic-data.h"
#include "garlic-centromeres.h"
#include "param_t.h"
#include <pthread.h>
#include <cmath>
#include <iostream>
#include "garlic-errlog.h"
#include "garlic-kde.h"
#include "gmm.h"
#include "gsl/gsl_statistics.h"
#include "gsl/gsl_sort.h"
#include "BoundFinder.h"

using namespace std;

struct work_order_t
{
  int id;
  int first_index;
  int last_index;

  int* winsize;
  double error;
  int MAX_GAP;

  vector<int_pair_t> *popChrPairs;

  vector< IndData * > *indDataByPop;
  vector< MapData * > *mapDataByChr;

  vector< vector< HapData * >* > *hapDataByPopByChr;
  vector< vector< FreqData * >* > *freqDataByPopByChr;

  vector< vector< WinData * >* > *winDataByPopByChr;

  centromere *centro;

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


void calcLOD(IndData *indData, MapData *mapData,
             HapData *hapData, FreqData *freqData,
             WinData *winData, centromere *centro,
             int winsize, double error, int MAX_GAP);

double lod(const short &genotype, const double &freq, const double &error);

vector< WinData * > *calcLODWindows(vector< HapData * > *hapDataByChr,
                                    vector< FreqData * > *freqDataByChr,
                                    vector< MapData * > *mapDataByChr,
                                    IndData *indData,
                                    centromere *centro,
                                    int winsize, double error, int MAX_GAP);

vector< ROHData * > *assembleROHWindows(vector< WinData * > *winDataByChr,
                                        vector< MapData * > *mapDataByChr,
                                        IndData *indData,
                                        centromere *centro,
                                        double lodScoreCutoff,
                                        ROHLength **rohLength,
                                        int winSize,
                                        int MAX_GAP);

ROHLength *initROHLength(int size, string pop);
void releaseROHLength(ROHLength *rohLength);
void releaseROHLength(vector< ROHLength * > *rohLengthByPop);

vector< ROHData * > *initROHData(IndData *indData);
void writeROHData(string outfile,
                  vector< ROHData * > *rohDataByInd,
                  vector< MapData * > *mapDataByChr,
                  int_pair_t bounds,
                  string popName,
                  string version);

string makeROHFilename(string outfile);

double selectLODCutoff(vector< WinData * > *winDataByChr, int KDE_SUBSAMPLE, string kdeoutfile);
void exploreWinsizes(vector< HapData * > *hapDataByChr,
                     vector< FreqData * > *freqDataByChr,
                     vector< MapData * > *mapDataByChr,
                     IndData *indData,
                     centromere *centro,
                     vector<int> &multiWinsizes,
                     double error,
                     int MAX_GAP, int KDE_SUBSAMPLE, string outfile);

int selectWinsize(vector< HapData * > *hapDataByChr,
                  vector< FreqData * > *freqDataByChr,
                  vector< MapData * > *mapDataByChr,
                  IndData *indData, centromere *centro,
                  int winsize, double error,
                  int MAX_GAP, int KDE_SUBSAMPLE);

bool inGap(int qStart, int qEnd, int targetStart, int targetEnd);

int_pair_t selectSizeClasses(ROHLength *rohLength);

extern pthread_mutex_t cerr_mutex;

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
