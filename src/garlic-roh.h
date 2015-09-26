#ifndef __GARLIC_ROH_H__
#define __GARLIC_ROH_H__
#include "garlic-data.h"
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
             WinData *winData, int winsize, double error, int MAX_GAP);
double lod(const short &genotype, const double &freq, const double &error);

vector< vector< WinData * >* > *calcLODWindows(vector< vector< HapData * >* > *hapDataByPopByChr,
        vector< vector< FreqData * >* > *freqDataByPopByChr,
        vector< MapData * > *mapDataByChr,
        vector< IndData * > *indDataByPop,
        int* winsize, double error, int MAX_GAP, int numThreads);
vector< vector< WinData * >* > *calcLODWindowsSinglePop(vector< vector< HapData * >* > *hapDataByPopByChr,
        vector< vector< FreqData * >* > *freqDataByPopByChr,
        vector< MapData * > *mapDataByChr,
        vector< IndData * > *indDataByPop,
        int* winsize, double error, int MAX_GAP, int numThreads, int pop);

vector< vector< ROHData * >* > *initROHData(vector< IndData * > *indDataByPop);
vector< vector< ROHData * >* > *assembleROHWindows(vector< vector< WinData * >* > *winDataByPopByChr,
        vector< MapData * > *mapDataByChr,
        vector< IndData * > *indDataByPop,
        double *lodScoreCutoffByPop,
        vector< ROHLength * > **rohLengthByPop,
        int winSize,
        int MAX_GAP);
ROHLength *initROHLength(int size, string pop);
void releaseROHLength(ROHLength *rohLength);
void releaseROHLength(vector< ROHLength * > *rohLengthByPop);

void writeROHData(string outfile,
                  vector< vector< ROHData * >* > *rohDataByPopByInd,
                  vector< MapData * > *mapDataByChr,
                  double *shortMedBound,
                  double *medLongBound,
                  map<string, string> &ind2pop,
                  string version);

extern pthread_mutex_t cerr_mutex;


#endif
