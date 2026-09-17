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

//The same value, named for the other thing it means: a genotype that WAS
//parsed and found missing, as opposed to one not yet parsed.  The LOD lookup
//table maps anything outside {0,1,2} to its default slot, so -9 needs no
//special case downstream.
const geno_t GENO_MISSING = -9;



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

//Reads a VCF (plain or gzipped) into the same structures loadTPEDData
//produces, so nothing downstream knows which reader ran.  sampleIDs receives
//the sample names from the #CHROM line, in file order; main builds IndData
//from them, which --pop can then relabel.
void loadVCFData(string vcffile, int &numLoci, int &numInd,
                 vector< HapData * > **hapDataByChr,
                 vector< MapData * > **mapDataByChr,
                 vector< FreqData * > **freqDataByChr,
                 vector< GenoLikeData * > **GLDataByChr,
                 int nresample, bool PHASED, bool AUTO_FREQ, bool PASS_ONLY,
                 string GL_TYPE, vector<string> &sampleIDs);

void freqOnly(string tpedfile, string outfile, int nresample, char TPED_MISSING);

//The --freq-only pass for VCF input.  Streams like freqOnly rather than going
//through loadVCFData, which is the point of --freq-only: the genotypes are
//never held.  Writes the same five columns, with the ALT allele in the ALLELE
//column, so readFreqData's orientation check makes the file interchangeable
//with one written from a TPED.
void freqOnlyVCF(string vcffile, string outfile, int nresample, bool PASS_ONLY);

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

//indData is used only to name the individual in a diagnostic; the values are
//read positionally as before.
vector< GenoLikeData * > *readTGLSData(string filename,
                                       int expectedLoci,
                                       int expectedInd,
                                       vector< MapData * > *mapDataByChr,
                                       string GL_TYPE,
                                       IndData *indData);

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

//Row-consuming form, mirroring initHapData: a reader that allocates a row per
//locus, because the locus count is not known in advance, hands them here to be
//gathered into one contiguous block and freed.  It was already implemented in
//garlic-data.cpp but never declared, so nothing outside that file could reach
//it and nothing called it; loadVCFData is its first caller.  readTGLSData
//knows nloci up front and keeps using the dimensioned form above.
GenoLikeData *initGLData(const vector< double * > &GL, int nloci, int nind);
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

//Applies a --pop file to already-assembled IndData, matching BY SAMPLE ID.
//
//Format: whitespace-separated, '#' comments and blank lines skipped, read
//through igzstream so a .gz works like every other input.
//
//    <sample_id>  <population>  [sex]
//
//Two required columns in that order -- note this is the OPPOSITE order to a
//TFAM, whose first column is the population.  The optional third column is
//sex in PLINK coding (1 male, 2 female, 0 or -9 unknown); unlike a TFAM, which
//has to tolerate legacy variants, anything else is an error, because in a
//purpose-built file an unrecognised value is a typo.
//
//Every sample in the data must have a row.  Extra rows are allowed and
//counted: a cohort-wide population file used against a subset of the samples
//is a normal workflow.
void applyPopFile(const string &filename, IndData *indData);

//The duplicate-individual-ID error and the pooled-population warning.  Both
//used to live inside scanIndData3, which only the --tfam path calls, so under
//any other input path they silently stopped applying -- and the pooled-
//population warning is the one that catches a TFAM with population 0 for every
//sample, which is what vcf2tped.pl writes.  They are properties of the
//assembled metadata rather than of a TFAM's text, so they belong here, where
//every input path reaches them.  `source` names the file in the messages.
void checkIndData(IndData *indData, const string &source);

//True when a VCF REF or ALT field is a single unambiguous nucleotide, which is
//garlic's whole domain: a biallelic SNV coded 0/1/2.  Case-insensitive ACGT
//only -- 'N' is excluded deliberately, since a site whose REF or ALT is
//unknown cannot be coded as a dosage, as are the symbolic ALTs a VCF permits
//('*', '.', '<NON_REF>').
bool isSNV(const string &s);

//Position of a key within a VCF FORMAT string, or -1 if absent.  Shared by
//loadVCFData and freqOnlyVCF so the two readers cannot disagree about where a
//field is -- the failure mode the filterMonomorphic* overloads had.
int formatIndexOf(const string &fmt, const string &key);
int gtIndexOf(const string &fmt);   //formatIndexOf(fmt, "GT")

//Parses one sample's VCF genotype column.
//
//  sample..sampleEnd  the whole colon-separated sample field, e.g. "0|1:35:99"
//  gtIndex            position of GT within FORMAT (0 when GT is first, which
//                     the spec requires, but the reader does not assume)
//  altIndex           which ALT allele to count; 1 for a biallelic site
//
//Outputs:
//  dosage     number of copies of altIndex, or GENO_MISSING if ANY allele of
//             the call is '.'.  A half call such as 0/. is therefore missing,
//             which is what the TPED path does: loadTPEDData accumulates -9 per
//             missing allele and then clamps anything negative to -9.
//  firstCopy  whether the FIRST haplotype carries the counted allele, matching
//             loadTPEDData's firstCopy[i] = (alleleStr1 == oneAllele).
//  ploidy     number of alleles in the call.  parseGT does NOT reject ploidy
//             other than 2; the caller does, so the message can name the site
//             and the sample.
//  phased     false if any separator was '/'.  Vacuously true for a haploid
//             call, which has no separator.
//
//Returns false only for input that cannot be interpreted: no such sub-field,
//an empty GT, a non-digit allele, or an unexpected separator.
bool parseGT(const char *sample, const char *sampleEnd, int gtIndex, int altIndex,
             int &dosage, bool &firstCopy, int &ploidy, bool &phased);

//The per-genotype error rate from a whole PL array, which is the only correct
//way to use a VCF's PL or GL field.
//
//A VCF normalises PL so the CALLED genotype is exactly 0, so the single value
//at the called genotype carries no information at all -- it is 0 for every
//call, confident or not.  The information is in the OTHER entries: how much
//worse the alternatives are.  So
//
//    P(g) proportional to 10^(-PL_g/10),   error = 1 - P(called)/sum(P)
//
//which is exactly the definition of the VCF GQ field, and reproduces GQ to
//the digit where both are present.  For GL, PL = -10*GL.
//
//calledIndex is the genotype's index in the VCF ordering; for a biallelic
//diploid site that is the ALT dosage (0/0 -> 0, 0/1 -> 1, 1/1 -> 2).
//
//Computed by subtracting the minimum PL before exponentiating, so a large PL
//cannot underflow the sum to zero.
double plToError(const vector<double> &pl, int calledIndex);

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
