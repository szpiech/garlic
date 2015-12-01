#ifndef __GARLIC_ROH_H__
#define __GARLIC_ROH_H__
#include "garlic-data.h"
#include "garlic-centromeres.h"
#include "param_t.h"
#include <pthread.h>
#include <cmath>
#include <iostream>

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

void scan(void *work_order);
void scanSinglePop(void *order);
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
                                    int* winsize, double error, int MAX_GAP);
/*
vector< vector< WinData * >* > *calcLODWindowsSinglePop(vector< vector< HapData * >* > *hapDataByPopByChr,
    vector< vector< FreqData * >* > *freqDataByPopByChr,
    vector< MapData * > *mapDataByChr,
    vector< IndData * > *indDataByPop, centromere *centro,
    int* winsize, double error, int MAX_GAP, int numThreads, int pop);
*/
vector< ROHData * > *initROHData(IndData *indData);

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

void writeROHData(string outfile,
                  vector< ROHData * > *rohDataByInd,
                  vector< MapData * > *mapDataByChr,
                  double shortMedBound,
                  double medLongBound,
                  string popName,
                  string version);

bool inGap(int qStart, int qEnd, int targetStart, int targetEnd);

extern pthread_mutex_t cerr_mutex;


#endif
