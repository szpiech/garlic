#ifndef __GARLIC_DATA_H__
#define __GARLIC_DATA_H__
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <sstream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <cctype>
#include <map>
#include <pthread.h>
#include "gzstream.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "garlic-errlog.h"
#include "garlic-centromeres.h"
#include "garlic-pbar.h"

using namespace std;

const int MISSING = -9999;

struct int_pair_t
{
  int first;
  int second;
};

struct HapData
{
  bool **firstCopy;
  short **data;
  int nind;
  int nloci;
};

struct GenMapScaffold {
  int *physicalPos;
  double *geneticPos;
  map<int, int> ppos2index;
  int nloci;
  string chr;
  int centroStart;
  int centroEnd;
  int currentIndex;
};

struct MapData
{
  double *physicalPos;
  double *geneticPos;
  string *locusName;
  char *allele;
  //char *allele0;
  int nloci;
  string chr;
};

struct IndData
{
  string pop;
  string *indID;
  int nind;
};

struct FreqData
{
  double *freq;
  int nloci;
};

struct GenoFreqData
{
  double *homFreq;
  int nloci;
};

struct WinData
{
  double **data;
  int nind;
  int nloci;
  //int nmiss;
};

struct GenoLikeData
{
  double **data;
  int nind;
  int nloci;
  //int nmiss;
};

struct DoubleData
{
  double *data;
  int size;
};

struct LDData
{
  double **LD;
  int nloci;
  int winsize;
};

struct HR2_work_order_t
{
  //int id;
  int start;
  int stop;
  int winsize;
  HapData *hapData;
  GenoFreqData *genoFreqData;
  LDData *LD;
  Bar *bar;
};

struct R2_work_order_t
{
  //int id;
  int start;
  int stop;
  int winsize;
  HapData *hapData;
  FreqData *freqData;
  LDData *LD;
  Bar *bar;
};

double selectOverlapFrac(double variantDensity);

FreqData *initFreqData(const vector<double> &freq, int nloci);

HapData *initHapData(const vector< short * > &hap,
                     const vector< bool * > &fc,
                     int nloci, int nind, bool PHASED);

MapData *initMapData(const vector<double> &geneticPos,
                     const vector<double> &physicalPos,
                     const vector<string> &locusNames,
                     const vector<char> &allele,
                     int nloci, string chr);

void loadTPEDData(string tpedfile, int &numLoci, int &numInd,
                  vector< HapData * > **hapDataByChr,
                  vector< MapData * > **mapDataByChr,
                  vector< FreqData * > **freqDataByChr,
                  char TPED_MISSING, int nresample, bool PHASED, bool AUTO_FREQ);

void freqOnly(string tpedfile, string outfile, int nresample, char TPED_MISSING);

double calcDensity(int numLoci, vector< MapData * > *mapDataByChr, centromere *centro);

void parallelHR2(void *order);
void parallelR2(void *order);

unsigned int *make_thread_partition(int &num_threads, int nloci);

void ldHR2(LDData *LD, HapData *hapData, GenoFreqData *genoFreqData, int site, int start, int end);
void ldR2(LDData *LD, HapData *hapData, FreqData *freqData, int site, int start, int end);

LDData *calcHR2LD(HapData *hapData, GenoFreqData *genoFreqData, int winsize, int numThreads);
LDData *calcR2LD(HapData *hapData, FreqData *freqData, int winsize, int numThreads);

//double ld(HapData *hapData, GenoFreqData *genoFreqData, int site, int start, int end, int ind);
double hr2(HapData *hapData, GenoFreqData *genoFreqData, int i, int j);
double r2(HapData *hapData, FreqData *freqData, int i, int j);

vector< LDData * > *calcLDData(vector< HapData * > *hapDataByChr, 
                               vector< FreqData * > *freqDataByChr,
                               vector< MapData * > *mapDataByChr,
                               vector< GenoFreqData * > *genoFreqDataByChr,
                               centromere *centro,
                               int winsize,
                               int MAX_GAP,
                               bool PHASED,
                               int numThreads);

LDData *initLDData(int nloci, int winsize);
void releaseLDData(LDData *data);
void releaseLDData(vector< LDData * > *ldDataByChr);

GenoFreqData *initGenoFreq(int nloci);
void releaseGenoFreq(GenoFreqData *genoFreqData);
void releaseGenoFreq(vector< GenoFreqData * > *genoFreqDataByChr);

GenoFreqData *calculateGenoFreq(HapData *hapData);
vector< GenoFreqData * > *calculateGenoFreq(vector <HapData *> *hapDataByChr);

double getMapInfo(int queryPos, GenMapScaffold *scaffold, int &count);
double interpolate(double x0, double y0, double x1, double y1, double query);
int interpolateGeneticmap(vector< MapData * > **mapDataByChr, vector< GenMapScaffold * > *scaffoldMapByChr);
int interpolateGeneticmap(MapData *mapData, GenMapScaffold *scaffold);
vector< GenMapScaffold *> *loadMapScaffold(string mapfile, centromere *centro);

GenMapScaffold *initGenMapScaffold(int nloci);
void releaseGenMapScaffold(GenMapScaffold *scaffoldMap);
void releaseGenMapScaffold(vector< GenMapScaffold * > *scaffoldMapByChr);

int filterMonomorphicSites(vector< MapData * > **mapDataByChr,
                           vector< HapData * > **hapDataByChr,
                           vector< FreqData * > **freqDataByChr,
                           vector< GenoLikeData * > **GLDataByChr,
                           bool USE_GL, bool PHASED);

int filterMonomorphicAndOOBSites(vector< MapData * > **mapDataByChr,
                                 vector< HapData * > **hapDataByChr,
                                 vector< FreqData * > **freqDataByChr,
                                 vector< GenoLikeData * > **GLDataByChr,
                                 vector< GenMapScaffold * > *scaffoldMapByChr,
                                 bool USE_GL, bool PHASED);

MapData *filterMonomorphicSites(MapData *mapData, FreqData *freqData, int &newLoci);
HapData *filterMonomorphicSites(HapData *hapData, FreqData *freqData, int &newLoci, bool PHASED);
GenoLikeData *filterMonomorphicSites(GenoLikeData *GLData, FreqData *freqData, int &newLoci);
FreqData *filterMonomorphicSites(FreqData *freqData, int &newLoci);

MapData *filterMonomorphicAndOOBSites(MapData *mapData, FreqData *freqData, GenMapScaffold *scaffold, int &newLoci);
HapData *filterMonomorphicAndOOBSites(HapData *hapData, MapData *mapData, FreqData *freqData, GenMapScaffold *scaffold, int &newLoci, bool PHASED);
GenoLikeData *filterMonomorphicAndOOBSites(GenoLikeData *GLData, MapData *mapData, FreqData *freqData, GenMapScaffold *scaffold, int &newLoci);
FreqData *filterMonomorphicAndOOBSites(FreqData *freqData, MapData *mapData, GenMapScaffold *scaffold, int &newLoci);

string getPost(int num);
bool goodDouble(string str);

FreqData *initFreqData(int nloci);
void releaseFreqData(FreqData *data);
void releaseFreqData(vector< FreqData * > *freqDataByChr);

FreqData *calcFreqData(HapData *data, int nresample, const gsl_rng *r);
vector< FreqData * > *calcFreqData2(vector< HapData * > *hapDataByChr, int nresample);
void writeFreqData(string freqOutfile, string popName,
                   vector< FreqData * > *freqDataByChr,
                   vector< MapData * > *mapDataByChr,
                   IndData *indData);

vector< FreqData * > *readFreqData(string freqfile, string popName,
                                   vector< MapData * > *mapDataByChr);

vector< GenoLikeData * > *readTGLSData(string filename,
                                       int expectedLoci,
                                       int expectedInd,
                                       vector< MapData * > *mapDataByChr,
                                       string GL_TYPE);

MapData *initMapData(int nloci);
void releaseMapData(MapData *data);
void releaseMapData(vector< MapData * > *mapDataByChr);

void scanIndData3(string filename, int &numInd, string &popName);
IndData *readIndData3(string filename, int numInd);
IndData *initIndData(int nind);
void releaseIndData(IndData *data);

HapData *initHapData(unsigned int nind, unsigned int nloci, bool PHASED);
void releaseHapData(HapData *data);
void releaseHapData(vector< HapData * > *hapDataByChr);

GenoLikeData *initGLData(unsigned int nind, unsigned int nloci);
void releaseGLData(GenoLikeData *data);
void releaseGLData(vector< GenoLikeData * > *GLDataByChr);

void subsetData(vector< HapData * > *hapDataByChr,
                vector< GenoLikeData *> *GLDataByChr,
                IndData *indData,
                vector< HapData * > **subsetHapDataByChr,
                vector< GenoLikeData *> **subsetGLDataByChr,
                IndData **subsetIndData,
                int subsample, bool USE_GL, bool PHASED);


WinData *initWinData(unsigned int nind, unsigned int nloci);
vector< WinData * > *initWinData(vector< MapData * > *mapDataByChr, int nind);
void releaseWinData(WinData *data);
void releaseWinData(vector< WinData * > *winDataByChr);
void writeWinData(vector< WinData * > *winDataByChr,
                  IndData *indData,
                  vector< MapData * > *mapDataByChr,
                  string outfile);

DoubleData *initDoubleData(int n);
DoubleData *convertWinData2DoubleData(vector< WinData * > *winDataByChr, int step);
DoubleData *convertSubsetWinData2DoubleData(vector< WinData * > *winDataByChr, IndData *indData, int subsample, int step);
void releaseDoubleData(DoubleData *data);
void writeDoubleData(vector < DoubleData * > *rawWinDataByPop, vector< MapData * > *mapDataByChr, vector< IndData * > *indDataByPop);

//counts the number of "fields" in a string
//where a field is defined as a contiguous set of non whitespace
//characters and fields are delimited by whitespace
int countFields(const string &str);
string lc(string str);
string checkChrName(string chr);
#endif
