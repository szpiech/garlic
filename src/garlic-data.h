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
#include "garlic-pos.h"
#include "garlic-matrix.h"
#include "garlic-math.h"
#include <pthread.h>
#include "gzstream.h"
#include "garlic-errlog.h"
#include "garlic-centromeres.h"
#include "garlic-pbar.h"

using namespace std;

const int MISSING = -9999;

//Single process-wide RNG, seeded once from main so that runs are reproducible.
//Previously every consumer allocated its own generator and seeded it from
//time(NULL), which made results non-reproducible and -- because time(NULL) has
//one-second resolution -- frequently gave several "independent" generators the
//same seed.  All consumers are in serial code; do not call these from a worker
//thread without adding synchronisation.
void initRNG(unsigned long int seed);
GarlicRNG *getRNG();
void freeRNG();
//Draws a nondeterministic seed, for when the user did not supply one.
unsigned long int drawRandomSeed();

struct int_pair_t
{
  int first;
  int second;
};


//Genotypes take only the values {0, 1, 2} and MISSING(-9), so one signed byte
//is enough; they used to be stored as short, twice the size for no benefit.
//(2.0 GB -> 1.0 GB for a 1000-individual, 1M-locus callset.)
typedef signed char geno_t;

//The value HapData is filled with before the parser overwrites every element.
//MISSING is -9999, which does not fit in a geno_t and silently wrapped to -15.
//That was inert -- nothing compares a genotype against MISSING (the parser
//encodes a missing allele by accumulating -9 per allele), and the LOD lookup
//maps any genotype outside {0,1,2} to the same default branch lod() used -- but
//a sentinel that is not the value it claims to be is worth not having.
const geno_t GENO_UNSET = -9;



//All four of the bulk structures below are allocated as ONE contiguous block
//with the row pointers indexing into it, rather than one allocation per locus
//(577,489 separate new[] calls on the bundled example).  The release functions
//therefore free row 0 and the pointer array, not every row.
struct HapData
{
  //Matrix<unsigned char> rather than Matrix<bool>: vector<bool> is the proxy
  //specialisation and has no bool* to hand out (see garlic-matrix.h).  Only
  //populated under --phased; empty() otherwise.
  Matrix<unsigned char> firstCopy;
  Matrix<geno_t> data;
  int nind;
  int nloci;
};

struct GenMapScaffold {
  //Raw arrays over new[]/delete[] became vectors: the releaseX function no
  //longer has to enumerate them, which is how B11-style leaks happened.

  vector<pos_t> physicalPos;
  vector<double> geneticPos;
  map<pos_t, int> ppos2index;
  int nloci;
  string chr;
  pos_t centroStart;
  pos_t centroEnd;
  int currentIndex;
};

struct MapData
{
  //Raw arrays over new[]/delete[] became vectors: the releaseX function no
  //longer has to enumerate them, which is how B11-style leaks happened.

  vector<pos_t> physicalPos;
  vector<double> geneticPos;
  vector<string> locusName;
  vector<char> allele;
  //char *allele0;
  int nloci;
  string chr;
};

struct IndData
{
  //These were string*/int* over new[]/delete[].  B1 was exactly that: subsetData
  //did `newIndData->pop = indData->pop`, which aliased the caller's array and
  //then freed it twice, aborting --auto-winsize.  As vectors that assignment
  //copies and the bug cannot be written.
  vector<string> pop;
  vector<string> indID;
  //TFAM column 5, PLINK coding: 1 male, 2 female, 0 or -9 unknown.  Only used
  //to make the sex-chromosome warning specific about how many individuals
  //would be affected; garlic does not condition any calculation on it.
  vector<int> sex;
  int nind;
};

struct FreqData
{
  vector<double> freq;
  int nloci;
};

struct GenoFreqData
{
  vector<double> homFreq;
  int nloci;
};

struct WinData
{
  //Row per INDIVIDUAL here, unlike HapData/GenoLikeData which are row per locus.
  Matrix<double> data;
  int nind;
  int nloci;
  //int nmiss;
};

struct GenoLikeData
{
  Matrix<double> data;
  int nind;
  int nloci;
  //int nmiss;
};

struct DoubleData
{
  vector<double> data;
  int size;
};

struct LDData
{
  Matrix<double> LD;
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
  int *indIndex;
  int ldSubsample;
  double *band;
};

//Banded pairwise LD: band[i*winsize + d] holds the LD between locus i and
//locus i+d (d = 0 is the self pair and is 1 by definition).  Every pair that
//any window needs lies in this band, so each pair is evaluated exactly once
//instead of once per window that contains it.
struct BAND_work_order_t
{
  int start;
  int stop;
  int winsize;
  int nloci;
  double *band;
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
  int *indIndex;
  int ldSubsample;
  double *band;
};

double selectOverlapFrac(double variantDensity, int winsize);
//Coefficients of the empirical overlap-vs-density fit used by
//--auto-overlap-frac; hardcoded as 6.375 / 63.888.
void setAutoOverlapCoef(double slope, double intercept);

//Restrict the analysis to the named chromosomes.  Pruning happens after the
//files are read rather than during: the TGLS file has no chromosome column and
//is consumed positionally against the TPED, so dropping TPED rows at read time
//would desynchronise it.  Peak memory is therefore unchanged; runtime after
//the read is not.
//True for X, Y, chrX, chrY, 23, 24, chr23, chr24 (PLINK numbers X as 23 and
//Y as 24).  Names arrive already normalised by checkChrName.
bool isSexChromosome(const string &chr);

//Hemizygous male genotypes on chrX are written as homozygous calls in a TPED
//and are indistinguishable from true autozygosity, so a male X chromosome
//looks like one chromosome-length run.  Warns, and reports how many
//individuals are coded male if the TFAM said.  Returns true if any sex
//chromosome is present.
bool warnSexChromosomes(vector< MapData * > *mapDataByChr, IndData *indData);

int filterChromosomes(vector<string> &keep,
                      vector< MapData * > **mapDataByChr,
                      vector< HapData * > **hapDataByChr,
                      vector< FreqData * > **freqDataByChr,
                      vector< GenoLikeData * > **GLDataByChr,
                      bool USE_GL);

FreqData *initFreqData(const vector<double> &freq, int nloci);

HapData *initHapData(const vector< geno_t * > &hap,
                     const vector< bool * > &fc,
                     int nloci, int nind, bool PHASED);

MapData *initMapData(const vector<double> &geneticPos,
                     const vector<pos_t> &physicalPos,
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
void parallelLDFromBand(void *order);
void ldRowsFromBand(double *band, LDData *LD, int nloci, int winsize, int start, int stop, Bar *bar);

//Returns by value: the previous version handed back a new[] array that every
//caller had to remember to delete, and an early return skipped it.
vector<unsigned int> make_thread_partition(int &num_threads, int nloci);

void ldHR2(LDData *LD, HapData *hapData, GenoFreqData *genoFreqData, int site, int start, int end, int *indIndex, int ldSubsample);
void ldR2(LDData *LD, HapData *hapData, FreqData *freqData, int site, int start, int end, int *indIndex, int ldSubsample);

LDData *calcHR2LD(HapData *hapData, GenoFreqData *genoFreqData, int winsize, int numThreads, int *indIndex, int ldSubsample);
LDData *calcR2LD(HapData *hapData, FreqData *freqData, int winsize, int numThreads, int *indIndex, int ldSubsample);

//double ld(HapData *hapData, GenoFreqData *genoFreqData, int site, int start, int end, int ind);
double hr2(HapData *hapData, GenoFreqData *genoFreqData, int i, int j, int *indIndex, int ldSubsample);
double r2(HapData *hapData, FreqData *freqData, int i, int j, int *indIndex, int ldSubsample);

vector< LDData * > *calcLDData(vector< HapData * > *hapDataByChr, 
                               vector< FreqData * > *freqDataByChr,
                               vector< MapData * > *mapDataByChr,
                               vector< GenoFreqData * > *genoFreqDataByChr,
                               centromere *centro,
                               int winsize,
                               int MAX_GAP,
                               bool PHASED,
                               int numThreads,
                               int ldSubsample);

LDData *initLDData(int nloci, int winsize);
void releaseLDData(LDData *data);
void releaseLDData(vector< LDData * > *ldDataByChr);

GenoFreqData *initGenoFreq(int nloci);
void releaseGenoFreq(GenoFreqData *genoFreqData);
void releaseGenoFreq(vector< GenoFreqData * > *genoFreqDataByChr);

GenoFreqData *calculateGenoFreq(HapData *hapData);
vector< GenoFreqData * > *calculateGenoFreq(vector <HapData *> *hapDataByChr);

double getMapInfo(pos_t queryPos, GenMapScaffold *scaffold, int &count);
double interpolate(double x0, double y0, double x1, double y1, double query);
int interpolateGeneticmap(vector< MapData * > **mapDataByChr, vector< GenMapScaffold * > *scaffoldMapByChr);
//Reorders scaffoldMapByChr so entry i is the scaffold for mapDataByChr->at(i).
//Downstream code zips the two vectors positionally; without this, a map file
//whose chromosomes appear in a different order than the TPED silently applies
//the wrong chromosome's genetic map.  Returns false on any mismatch.
bool alignMapScaffold(vector< GenMapScaffold * > *scaffoldMapByChr, vector< MapData * > *mapDataByChr);
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

//The site-retention predicate, exposed so it can be unit-tested.  NULL
//scaffold = monomorphic only; non-NULL also drops out-of-map and in-centromere
//sites.  Returns the indices of the retained sites in the original arrays.
vector<int> keepSites(const MapData *mapData, const FreqData *freqData,
                      const GenMapScaffold *scaffold);



string getPost(int num);
bool goodDouble(string str);

FreqData *initFreqData(int nloci);
void releaseFreqData(FreqData *data);
void releaseFreqData(vector< FreqData * > *freqDataByChr);

void writeFreqData(string freqOutfile,
                   vector< FreqData * > *freqDataByChr,
                   vector< MapData * > *mapDataByChr,
                   IndData *indData);

vector< FreqData * > *readFreqData(string freqfile,
                                   vector< MapData * > *mapDataByChr);

vector< GenoLikeData * > *readTGLSData(string filename,
                                       int expectedLoci,
                                       int expectedInd,
                                       vector< MapData * > *mapDataByChr,
                                       string GL_TYPE);

MapData *initMapData(int nloci);
void releaseMapData(MapData *data);
void releaseMapData(vector< MapData * > *mapDataByChr);

void scanIndData3(string filename, int &numInd);
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

//Convert one genotype-quality/likelihood value into the per-genotype error
//rate the LOD calculation uses.  Extracted from readTGLSData so it can be
//tested directly: both TGLS files bundled with garlic are constant (GQ 30
//everywhere, GL -0.0004 everywhere), so no end-to-end test can exercise
//per-genotype variation in this conversion.
//  GQ: error = 10^(-GQ/10)        a Phred-scaled quality
//  GL: error = 1 - 10^(GL)        GL is a log10 likelihood
//  PL: error = 1 - 10^(-PL/10)    PL is Phred-scaled likelihood
//The exponent is floored at -10 in every branch, and the result clamped to
//(0, 1].
double glToError(double value, string glType);
#endif
