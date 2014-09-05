#ifndef __ROHSCAN_DATA_H__
#define __ROHSCAN_DATA_H__
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <ctime>
#include "gzstream.h"
#include "gsl/gsl_rng.h"

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
    int nhaps;
    int nloci;
};

struct MapData
{
    int *physicalPos;
    double *geneticPos;
    string *locusName;
    int nloci;
    int chr;
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

struct DoubleData
{
    double *data;
    int size;
};

string getPost(int num);

FreqData *initFreqData(int nloci);
void releaseFreqData(FreqData *data);
void releaseFreqData(vector< vector< FreqData * >* > *freqDataByPopByChr);
FreqData *calcFreqData(HapData *data, int nresample, const gsl_rng *r);
vector< vector< FreqData * >* > *calcFreqData(vector< vector< HapData * >* > *hapDataByPopByChr, int nresample);

vector< int_pair_t > *scanTFAMData(string filename, int &numInd);
vector< IndData * > *readTFAMData(string filename, vector< int_pair_t > *indCoordList);

vector< int_pair_t > *scanTPEDMapData(string filename, int &numLoci);
vector< MapData * > *readTPEDMapData(string filename, vector< int_pair_t > *chrCoordList);

vector< vector< HapData * >* > *readTPEDHapData(string filename,
        int expectedLoci,
        int expectedInd,
        vector< int_pair_t > *chrCoordList,
        vector< int_pair_t > *indCoordList);

MapData *initMapData(int nloci);
void releaseMapData(MapData *data);
void releaseMapData(vector< MapData * > *mapDataByChr);
vector< MapData * > *readMapData(string filename, vector< int_pair_t > *chrCoordList);
vector< int_pair_t > *scanMapData(string filename, int &numLoci);

vector< int_pair_t > *scanIndData(string filename, int &numInd);
vector< IndData * > *readIndData(string filename, vector< int_pair_t > *indCoordList);
IndData *initIndData(int nind);
void releaseIndData(vector< IndData * > *indDataByPop);
void releaseIndData(IndData *data);

HapData *initHapData(unsigned int nhaps, unsigned int nloci);
void releaseHapData(HapData *data);
void releaseHapData(vector< vector< HapData * >* > *hapDataByPopByChr);
vector< vector< HapData * >* > *readHapData(string filename,
        int expectedLoci,
        int expectedInd,
        vector< int_pair_t > *chrCoordList,
        vector< int_pair_t > *indCoordList);

vector< vector< WinData * >* > *initWinData(vector< MapData * > *mapDataByChr, vector< IndData * > *indDataByPop);
WinData *initWinData(unsigned int nind, unsigned int nloci);
void releaseWinData(WinData *data);
void releaseWinData(vector< vector< WinData * >* > *winDataByPopByChr);
void writeWinData(vector< vector< WinData * >* > *winDataByPopByChr,
                  vector< IndData * > *indDataByPop,
                  vector< MapData * > *mapDataByChr,
                  string outfile);

DoubleData *initDoubleData(int n);
vector < DoubleData * > *convertWinData2DoubleData(vector< vector< WinData * >* > *winDataByPopByChr);
void releaseDoubleData(DoubleData *data);
void releaseDoubleData(vector < DoubleData * > *rawWinDataByPop);
void writeDoubleData(vector < DoubleData * > *rawWinDataByPop, vector< MapData * > *mapDataByChr, vector< IndData * > *indDataByPop);

//counts the number of "fields" in a string
//where a field is defined as a contiguous set of non whitespace
//characters and fields are delimited by whitespace
int countFields(const string &str);

#endif
