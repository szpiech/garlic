#ifndef __ROHSCAN_ROH_H__
#define __ROHSCAN_ROH_H__
#include "rohscan-data.h"
#include "param_t.h"
#include <cmath>
#include <iostream>

using namespace std;

struct work_order_t
{
  int id;
  int first_index;
  int last_index;

  int winsize;
  double error;
  int MAX_GAP;

  vector<int_pair_t> *popChrPairs;

  vector< IndData* >* indDataByPop;
  vector< MapData* >* mapDataByChr;

  vector< vector< HapData* >* >* hapDataByPopByChr;
  vector< vector< FreqData* >* >* freqDataByPopByChr;

  vector< vector< WinData* >* >* winDataByPopByChr;

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
  int* length;
  int size;
};

void scan(void* work_order);
void calcLOD(IndData* indData, MapData* mapData,
	     HapData* hapData, FreqData* freqData,
	     WinData* winData, int winsize, double error, int MAX_GAP);
double lod(const short &genotype, const double &freq, const double &error);

vector< vector< ROHData* >* >* initROHData(vector< IndData* >* indDataByPop);
vector< vector< ROHData* >* >* assembleROHWindows(vector< vector< WinData* >* >* winDataByPopByChr,
						  vector< MapData* >* mapDataByChr,
						  vector< IndData* >* indDataByPop, 
						  double *lodScoreCutoffByPop,
						  vector< ROHLength* >* rohLengthByPop);
ROHLength* initROHLength(int size, string pop);
void releaseROHLength(ROHLength *rohLength);
void releaseROHLength(vector< ROHLength* >* rohLengthByPop);

extern pthread_mutex_t cerr_mutex;


#endif
