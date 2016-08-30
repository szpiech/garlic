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
#include "gzstream.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "garlic-errlog.h"
#include "garlic-centromeres.h"

using namespace std;

const int MISSING = -9999;

struct int_pair_t
{
  int first;
  int second;
};

struct HapData
{
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
  int *physicalPos;
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
                           vector< GenoLikeData * > **GLDataByChr);

int filterMonomorphicAndOOBSites(vector< MapData * > **mapDataByChr,
                                 vector< HapData * > **hapDataByChr,
                                 vector< FreqData * > **freqDataByChr,
                                 vector< GenoLikeData * > **GLDataByChr,
                                 vector< GenMapScaffold * > *scaffoldMapByChr);

MapData *filterMonomorphicSites(MapData *mapData, FreqData *freqData, int &newLoci);
HapData *filterMonomorphicSites(HapData *hapData, FreqData *freqData, int &newLoci);
GenoLikeData *filterMonomorphicSites(GenoLikeData *GLData, FreqData *freqData, int &newLoci);
FreqData *filterMonomorphicSites(FreqData *freqData, int &newLoci);

MapData *filterMonomorphicAndOOBSites(MapData *mapData, FreqData *freqData, GenMapScaffold *scaffold, int &newLoci);
HapData *filterMonomorphicAndOOBSites(HapData *hapData, MapData *mapData, FreqData *freqData, GenMapScaffold *scaffold, int &newLoci);
GenoLikeData *filterMonomorphicAndOOBSites(GenoLikeData *GLData, MapData *mapData, FreqData *freqData, GenMapScaffold *scaffold, int &newLoci);
FreqData *filterMonomorphicAndOOBSites(FreqData *freqData, MapData *mapData, GenMapScaffold *scaffold, int &newLoci);

string getPost(int num);
bool goodDouble(string str);

//map<string, double> readLODCutoff(string lodCutoffFile, map<string, int> &pop2size);
//void readBoundSizes(string boundSizeFile, map<string, double> &pop2SMbound, map<string, double> &pop2MLbound, map<string, int> &pop2size);

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
                                   vector< int_pair_t > *chrCoordList,
                                   vector< MapData * > *mapDataByChr);

vector< int_pair_t > *scanTFAMData(string filename, int &numInd);
vector< IndData * > *readTFAMData(string filename, vector< int_pair_t > *indCoordList);

vector< int_pair_t > *scanTPEDMapData(string filename, int &numLoci, int &numCols);
vector< MapData * > *readTPEDMapData(string filename, int numCols, vector< int_pair_t > *chrCoordList, char TPED_MISSING);
/*
void writeTPEDDataByPop(string outfile,
                        vector< vector< HapData * >* > *hapDataByPopByChr,
                        vector< MapData * > *mapDataByChr,
                        map<string, int> &pop2index);
void writeTFAMDataByPop(string outfile, vector< IndData * > *indDataByPop, map<string, int> &pop2index);
*/

vector< HapData * > *readTPEDHapData3(string filename,
                                      int expectedLoci,
                                      int expectedInd,
                                      char TPED_MISSING,
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

HapData *initHapData(unsigned int nind, unsigned int nloci);
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
                int subsample, bool USE_GL);


WinData *initWinData(unsigned int nind, unsigned int nloci);
vector< WinData * > *initWinData(vector< MapData * > *mapDataByChr, IndData *indData);
void releaseWinData(WinData *data);
void releaseWinData(vector< WinData * > *winDataByChr);
void writeWinData(vector< WinData * > *winDataByChr,
                  IndData *indData,
                  vector< MapData * > *mapDataByChr,
                  string outfile);

DoubleData *initDoubleData(int n);
DoubleData *convertWinData2DoubleData(vector< WinData * > *winDataByChr);
DoubleData *convertSubsetWinData2DoubleData(vector< WinData * > *winDataByChr, IndData *indData, int subsample);
void releaseDoubleData(DoubleData *data);
void writeDoubleData(vector < DoubleData * > *rawWinDataByPop, vector< MapData * > *mapDataByChr, vector< IndData * > *indDataByPop);

//counts the number of "fields" in a string
//where a field is defined as a contiguous set of non whitespace
//characters and fields are delimited by whitespace
int countFields(const string &str);
string lc(string str);
string checkChrName(string chr);
#endif
