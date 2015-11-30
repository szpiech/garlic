#ifndef __GARLIC_DATA_H__
#define __GARLIC_DATA_H__
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <ctime>
#include <cstdlib>
#include <cctype>
#include <map>
#include "gzstream.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

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

struct DoubleData
{
    double *data;
    int size;
};

string getPost(int num);
bool goodDouble(string str);

map<string, double> readLODCutoff(string lodCutoffFile, map<string, int> &pop2size);
void readBoundSizes(string boundSizeFile, map<string, double> &pop2SMbound, map<string, double> &pop2MLbound, map<string, int> &pop2size);

FreqData *initFreqData(int nloci);
void releaseFreqData(FreqData *data);
void releaseFreqData(vector< vector< FreqData * >* > *freqDataByPopByChr);
FreqData *calcFreqData(HapData *data, int nresample, const gsl_rng *r);
vector< vector< FreqData * >* > *calcFreqData(vector< vector< HapData * >* > *hapDataByPopByChr, int nresample);
/*
void writeFreqData(string freqOutfile,
                   vector< vector< FreqData * >* > *freqDataByPopByChr,
                   vector< MapData * > *mapDataByChr,
                   vector< IndData * > *indDataByPop);
*/
void writeFreqData(string freqOutfile, string popName,
                   vector< FreqData * > *freqDataByChr,
                   vector< MapData * > *mapDataByChr,
                   IndData *indData);
/*
vector< vector< FreqData * >* > *readFreqData(string freqfile,
        vector< int_pair_t > *chrCoordList,
        vector< MapData * > *mapDataByChr,
        map<string, int> &pop2index);
*/
vector< FreqData * > *readFreqData(string freqfile, string popName,
        vector< int_pair_t > *chrCoordList,
        vector< MapData * > *mapDataByChr);

vector< int_pair_t > *scanTFAMData(string filename, int &numInd);
vector< IndData * > *readTFAMData(string filename, vector< int_pair_t > *indCoordList);

vector< int_pair_t > *scanTPEDMapData(string filename, int &numLoci, int &numCols);
vector< MapData * > *readTPEDMapData(string filename, int numCols, vector< int_pair_t > *chrCoordList, char TPED_MISSING);
void writeTPEDDataByPop(string outfile,
                        vector< vector< HapData * >* > *hapDataByPopByChr,
                        vector< MapData * > *mapDataByChr,
                        map<string, int> &pop2index);
void writeTFAMDataByPop(string outfile, vector< IndData * > *indDataByPop, map<string, int> &pop2index);

/*
vector< vector< HapData * >* > *readTPEDHapData(string filename,
        int expectedLoci,
        int expectedInd,
        vector< int_pair_t > *chrCoordList,
        vector< int_pair_t > *indCoordList);
*/
/*
vector< vector< HapData * >* > *readTPEDHapData2(string filename,
        int expectedLoci,
        int expectedInd,
        vector< int_pair_t > *chrCoordList,
        string *indList,
        map<string, string> &ind2pop,
        map<string, int> &pop2size,
        map<string, int> &pop2index,
        char TPED_MISSING,
        vector< MapData * > *mapDataByChr);
*/
vector< HapData * > *readTPEDHapData3(string filename,
        int expectedLoci,
        int expectedInd,
        char TPED_MISSING,
        vector< MapData * > *mapDataByChr);

MapData *initMapData(int nloci);
void releaseMapData(MapData *data);
void releaseMapData(vector< MapData * > *mapDataByChr);
//vector< MapData * > *readMapData(string filename, vector< int_pair_t > *chrCoordList);
//vector< int_pair_t > *scanMapData(string filename, int &numLoci);

//vector< int_pair_t > *scanIndData(string filename, int &numInd);
//void scanIndData2(string filename, int &numInd, map<string, string> &ind2pop, map<string, int> &pop2size);
void scanIndData3(string filename, int &numInd, string &popName);
//vector< IndData * > *readIndData(string filename, vector< int_pair_t > *indCoordList);
/*
vector< IndData * > *readIndData2(string filename, int numInd,
                                  map<string, string> &ind2pop,
                                  map<string, int> &pop2size,
                                  string *indList,
                                  map<string, int> &pop2index);
*/
IndData *readIndData3(string filename, int numInd, string *indList);

IndData *initIndData(int nind);
void releaseIndData(vector< IndData * > *indDataByPop);
void releaseIndData(IndData *data);

HapData *initHapData(unsigned int nind, unsigned int nloci);
void releaseHapData(HapData *data);
void releaseHapData(vector< vector< HapData * >* > *hapDataByPopByChr);

void subsetData(vector< vector< HapData * >* > *hapDataByPopByChr,
                vector< IndData * > *indDataByPop,
                vector< vector< HapData * >* > **subsetHapDataByPopByChr,
                vector< IndData * > **subsetIndDataByPop, int subsample);
/*
vector< vector< HapData * >* > *readHapData(string filename,
        int expectedLoci,
        int expectedInd,
        vector< int_pair_t > *chrCoordList,
        vector< int_pair_t > *indCoordList);
*/
/*
vector< vector< HapData * >* > *readHapData2(string filename,
        int expectedLoci,
        int expectedInd,
        vector< int_pair_t > *chrCoordList,
        string *indList,
        map<string, string> &ind2pop,
        map<string, int> &pop2size,
        map<string, int> &pop2index,
        vector< MapData * > *mapDataByChr);
*/
WinData *initWinData(unsigned int nind, unsigned int nloci);
vector< vector< WinData * >* > *initWinData(vector< MapData * > *mapDataByChr, vector< IndData * > *indDataByPop);
vector< vector< WinData * >* > *initWinData(vector< MapData * > *mapDataByChr, vector< IndData * > *indDataByPop, int pop);
void releaseWinData(WinData *data);
void releaseWinData(vector< vector< WinData * >* > *winDataByPopByChr);
void writeWinData(vector< vector< WinData * >* > *winDataByPopByChr,
                  vector< IndData * > *indDataByPop,
                  vector< MapData * > *mapDataByChr,
                  string outfile);

DoubleData *initDoubleData(int n);
vector < DoubleData * > *convertWinData2DoubleData(vector< vector< WinData * >* > *winDataByPopByChr);
vector < DoubleData * > *convertSubsetWinData2DoubleData(vector< vector< WinData * >* > *winDataByPopByChr, int subsample);
void releaseDoubleData(DoubleData *data);
void releaseDoubleData(vector < DoubleData * > *rawWinDataByPop);
void writeDoubleData(vector < DoubleData * > *rawWinDataByPop, vector< MapData * > *mapDataByChr, vector< IndData * > *indDataByPop);

//counts the number of "fields" in a string
//where a field is defined as a contiguous set of non whitespace
//characters and fields are delimited by whitespace
int countFields(const string &str);
string lc(string str);
string checkChrName(string chr);
#endif
