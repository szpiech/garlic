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
#include <thread>
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
//Genotype encoding.  0, 1 and 2 are the number of copies of the counted
//allele.  Every negative value means "not a usable genotype" -- lod() returns
//0 for any of them and the LOD lookup table maps them all to its missing slot
//-- but they differ in what they contribute to an ALLELE FREQUENCY:
//
//   GENO_MISSING       both alleles unobserved        contributes nothing
//   GENO_HALF_COUNTED  one allele observed, and it IS the counted allele
//                                                     contributes 1 of 1
//   GENO_HALF_OTHER    one allele observed, and it is not the counted allele
//                                                     contributes 0 of 1
//
//A half call ("A 0" in a TPED, "0/." in a VCF) is unusable as a genotype, but
//the allele that WAS observed is real data.  loadTPEDData has always counted
//it towards the frequency while storing the genotype as missing; the VCF path
//discarded it entirely.  These codes let both paths agree, and let a frequency
//be recomputed for a subset of individuals from the matrix alone -- which the
//old encoding could not, because it threw the half call away.
const geno_t GENO_MISSING      = -9;
const geno_t GENO_HALF_COUNTED = -1;
const geno_t GENO_HALF_OTHER   = -2;

//Is this a usable genotype?  Any test for "not missing" must use this rather
//than `!= GENO_MISSING`, which the half codes pass.
inline bool genoIsCalled(geno_t g) { return g >= 0; }

//One genotype's contribution to an allele frequency.  Centralised for the same
//reason as alleleFrequency: four sites accumulate this, and they must agree
//about the half codes or the input paths diverge again.
inline void addAlleleCounts(geno_t g, double &nalleles, double &total)
{
    if (g >= 0)                      { nalleles += double(g); total += 2; }
    else if (g == GENO_HALF_COUNTED) { nalleles += 1;         total += 1; }
    else if (g == GENO_HALF_OTHER)   {                        total += 1; }
    //GENO_MISSING contributes nothing.
}



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
  //TFAM column 5, PLINK coding: 1 male, 2 female, 0 or -9 unknown.  What the
  //metadata SAYS; nothing conditions a calculation on it directly.
  vector<int> sex;
  //What the individual IS, on a shared sex chromosome: one of the Zygo
  //values.  Derived from sex and --sex-system, checked against observed
  //heterozygosity, and inferred from it when sex was not recorded -- so it is
  //a separate field rather than a reading of the one above.  ZYG_UNKNOWN
  //everywhere until the sex check runs, and on a run with no shared sex
  //chromosome it stays that way, because nothing needs it.
  vector<int> zygo;
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
//---- chromosome identity ----------------------------------------------------
//
//Case- and prefix-insensitive key for MATCHING one chromosome name against
//another, against a --build table, against a --chr list or against the
//sex-chromosome detector.  It is NOT a display name: output keeps whatever
//checkChrName produced, so existing files are byte-identical.
//
//  chrX  chrx  CHRX  ChrX  X  x      -> x
//  chr1  1  01  chr01              -> 1
//  chr4a  4a  Chr4A               -> 4a
//  contig7  Contig7               -> contig7
//
//Lowercasing is done by hand over 'A'..'Z' rather than with tolower(), which
//is locale-dependent, or std::tolower(char), which is undefined for negative
//char values -- and every byte of a UTF-8 scaffold name above 0x7F is one.
string canonChrKey(const string &name);

//Refuses when two DISTINCT display names in the data share a key ("X" and
//"chrX", or "1" and "01").  Those are two independent MapData entries today,
//and matching them permissively without saying so would silently merge or
//silently pick one.  This also catches a chromosome whose rows are not
//contiguous in the input, because the reader closes a chromosome whenever
//column 1 changes and so produces two entries with the same name.
//Returns 0, or -1 after logging.
int checkChrKeyCollisions(vector< MapData * > *mapDataByChr);

//---- sex chromosomes --------------------------------------------------------
//
//garlic does not decide what a chromosome is.  The user declares it, and one
//table drives every path: there is no XY code path and no ZW code path.  The
//role says what a chromosome is; the zygosity says what an individual is; the
//expected ploidy of an individual at a locus is the product of the two.
//
//  role \ zygosity   homogametic   heterogametic
//  AUTOSOME               2             2
//  SEX_SHARED             2             1        (X under XY, Z under ZW)
//  SEX_DEGENERATE         0             1        (Y under XY, W under ZW)
//  HAPLOID                1             1        (mitochondrion, chloroplast)
//  PAR                    2             2        (diploid in both, not called)
//
//Only AUTOSOME and SEX_SHARED can carry a run of homozygosity in anyone, so
//the other three are dropped before the pipeline starts.
enum ChrRole
{
    CHR_AUTOSOME = 0,
    CHR_SEX_SHARED = 1,
    CHR_SEX_DEGENERATE = 2,
    CHR_HAPLOID = 3,
    CHR_PAR = 4
};

//What an individual is on a shared sex chromosome.  Two copies, one copy, or
//not established -- which is a real third state, because sex is optional in a
//TFAM and in a --pop file and heterozygosity does not always resolve it.
enum Zygo
{
    ZYG_UNKNOWN = 0,
    ZYG_HOMOGAMETIC = 1,
    ZYG_HETEROGAMETIC = 2
};

//PLINK sex coding is 1 male / 2 female whatever the species' system is, so in
//a ZW species the HOMOGAMETIC sex is the one coded 1.  Which code is
//heterogametic is the one bit of this that no amount of data can supply, and
//--sex-system is how the user supplies it.
const int SEX_SYSTEM_UNSET = 0;   //nothing said: detection refuses rather than guessing
const int SEX_SYSTEM_XY    = 1;   //code 1 (male) is heterogametic
const int SEX_SYSTEM_ZW    = 2;   //code 2 (female) is heterogametic
const int SEX_SYSTEM_NONE  = 3;   //asserted: every chromosome diploid in everyone

//-1 if the string is not one of xy / zw / none.
int parseSexSystem(const string &s);
string sexSystemName(int system);

//Does this key name something that is probably not an autosome?
//
//  always              x y z w m mt
//  human --build only  23 24 25 26   (PLINK's human codes for X, Y, PAR, MT)
//
//The gate on the numbers is the whole answer to "chr24 is a real autosome in
//most species": a bare number means a sex chromosome only under PLINK's human
//convention, and a human --build is the user asserting human coordinates.
//Without it the numbers are AMBIGUOUS rather than autosomal -- see
//isAmbiguousNumericChrKey -- and garlic refuses instead of picking a reading.
bool isDetectedSexChrKey(const string &key);
bool isAmbiguousNumericChrKey(const string &key);
bool isHumanBuild(const string &build);

struct SexModel
{
    int system;
    vector<ChrRole> role;          //parallel to mapDataByChr
    map<string, ChrRole> byKey;    //survives chromosome filtering; see rebuild()
    SexModel() : system(SEX_SYSTEM_UNSET) {}

    bool anyOfRole(ChrRole r) const
    {
        for (unsigned int i = 0; i < role.size(); i++) if (role[i] == r) return true;
        return false;
    }

    //role is positional, so any filtering of mapDataByChr invalidates it.
    //byKey does not move, so the vector is rebuilt from it rather than
    //recomputed -- which also means a role can never be re-derived
    //differently after a filter than it was before one.
    void rebuild(vector< MapData * > *mapDataByChr);
};

//Which zygosity a recorded sex implies under a declared system.  This is the
//ONLY place the xy/zw distinction has any consequence: PLINK codes 1 male and
//2 female whatever the species does, so the system decides which of those two
//codes is the heterogametic one and nothing else about the analysis differs.
int zygosityForSex(int system, int sex);

//May this individual carry a run of homozygosity on this chromosome?  An
//autosome, yes; the shared sex chromosome, only if homogametic; anything else
//has expected ploidy of at most one in everybody and never reaches here.
bool eligibleForCalling(ChrRole role, int zygo);

//Per-individual heterozygosity on the shared sex chromosome(s), used to check
//the recorded sex and to infer it where it was not recorded, then to recode
//the genotypes that the resulting zygosity makes impossible:
//
//  heterogametic  0 -> GENO_HALF_OTHER, 2 -> GENO_HALF_COUNTED, 1 -> MISSING
//  unknown        everything -> MISSING
//
//The half codes are what make the allele frequency come out right with no
//change to the frequency code: a hemizygous call contributes its one observed
//allele and no genotype, which is exactly what it is.  A heterozygous call
//where the individual has one copy of the chromosome is impossible by
//construction, so it is counted and discarded rather than believed.
//
//Writes <outfile>.sexcheck.tsv.  Returns 0, or -1 after logging.
int runSexCheck(vector< HapData * > *hapDataByChr,
                vector< MapData * > *mapDataByChr,
                IndData *indData,
                const SexModel &model,
                double hetLo, double hetHi,
                const string &outfile,
                bool quietCheck);

//Recomputes the allele frequency of the shared sex chromosome(s) after the
//recode above, and of nothing else.  The readers compute frequencies while
//reading, before anything knows what a chromosome is, so the sex chromosome's
//are counted as though every individual were diploid.  Recomputing only the
//chromosomes whose genotypes actually changed keeps every autosome's
//frequency bit-for-bit what it was -- which matters because --resample draws
//random numbers, so a recomputation is not a no-op there.
void recomputeFreqForRole(vector< HapData * > *hapDataByChr,
                          vector< FreqData * > *freqDataByChr,
                          const SexModel &model,
                          ChrRole role,
                          int nresample);

//Fills model.role from the declarations, applying the conventional names for
//the declared system, and refuses (returns -1, having logged) when a detected
//chromosome is left with no role -- which is the "detected but undeclared is
//an error" rule.  humanBuild licenses the numeric codes.  Explicit names that
//are not in the data are an error, as with --chr.
int buildSexModel(SexModel &model,
                  vector< MapData * > *mapDataByChr,
                  IndData *indData,
                  int system,
                  const vector<string> &sexChr,
                  const vector<string> &degenerateChr,
                  const vector<string> &haploidChr,
                  bool humanBuild,
                  bool autosomesOnly);

//How many individuals each sex code covers, for messages that would otherwise
//state a number the metadata cannot support.
void countSexCodes(IndData *indData, int &males, int &females, int &unknown);

//---- pseudoautosomal regions ------------------------------------------------
//
//A PAR is diploid in both sexes, so it is neither hemizygous nor worth calling
//at typical window sizes, and garlic drops the loci inside one.
//
//There is always a SET of them.  Humans have two (PAR1 and PAR2), other
//species have more, and a system with several shared sex chromosomes can have
//them on each -- so nothing here may assume one interval, or that an interval
//is at the end of a chromosome.  This is deliberately NOT the container the
//centromere table uses: that is one map<string,pos_t> pair keyed by chromosome
//and a second row for the same chromosome silently overwrites the first, which
//would lose PAR1 without an error.
//
//A PAR is identified by COORDINATES on a declared sex chromosome, never by a
//chromosome code: PLINK's 25 collides with a real autosome in any species with
//25 or more chromosomes, which is most of them, while coordinates are read in
//the same assembly as the input and cannot collide.
struct Interval
{
    pos_t start, end;   //inclusive at both ends, as PAR coordinates are written
    Interval() : start(0), end(0) {}
    Interval(pos_t s, pos_t e) : start(s), end(e) {}
};

class ExcludedRegions
{
public:
    void add(const string &chr, pos_t s, pos_t e);
    //Sorts and merges overlapping or abutting intervals per chromosome, and
    //rejects an inverted or empty one.  Returns 0, or -1 after logging.
    int finalise();
    bool empty() const { return byChr.empty(); }
    //Intervals of a chromosome, in the order they were merged into; NULL when
    //it has none.  Keyed by canonChrKey, so any spelling finds them.
    const vector<Interval> *get(const string &chr) const;
    bool contains(const string &chr, pos_t p) const;
    //Length of the part of [lo, hi] that this chromosome's intervals cover.
    pos_t overlap(const string &chr, pos_t lo, pos_t hi) const;
    //Every chromosome named, in insertion order, for validation messages.
    const vector<string> &chromosomes() const { return order; }

    //The spelling the user gave, for messages: the key is lower case with any
    //"chr" stripped, which is not what they typed.
    string spelling(const string &key) const;

private:
    map<string, vector<Interval> > byChr;
    map<string, string> asGiven;
    vector<string> order;
};

//Parses "chr:start-end" items, comma or whitespace separated, and the three
//column <chr> <start> <end> file form -- in which, unlike --centromere,
//repeated rows for a chromosome ADD a region rather than replacing one.
//Return 0, or -1 after logging.
int parsePARSpecs(const vector<string> &specs, ExcludedRegions &par);
int readPARFile(const string &filename, ExcludedRegions &par);

//Drops every locus inside an excluded region and reports how many each
//INTERVAL removed -- a total is diagnostic enough for one PAR, but with
//several it is the per-interval counts that reveal a coordinate given in the
//wrong assembly.  Refuses an interval on a chromosome that is not a declared
//shared sex chromosome: a PAR is defined relative to a sex-chromosome pair,
//and an interval on chr7 is a mistake rather than a region to drop.
//Returns the number of loci left, or -1 after logging.
int dropExcludedSites(vector< MapData * > **mapDataByChr,
                      vector< HapData * > **hapDataByChr,
                      vector< FreqData * > **freqDataByChr,
                      vector< GenoLikeData * > **GLDataByChr,
                      const ExcludedRegions &par,
                      const SexModel &model,
                      bool USE_GL, bool PHASED);

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

//popOfInd carries one population label per individual, in file order.  Empty,
//or all one label, gives the single FREQ column every earlier version wrote.
//With several labels the output is the wide format -- one column per
//population -- which is what --freq-file reads back.
void freqOnly(string tpedfile, string outfile, int nresample, char TPED_MISSING,
              const vector<string> &popOfInd);

//The --freq-only pass for VCF input.  Streams like freqOnly rather than going
//through loadVCFData, which is the point of --freq-only: the genotypes are
//never held.  Writes the same five columns, with the ALT allele in the ALLELE
//column, so readFreqData's orientation check makes the file interchangeable
//with one written from a TPED.
//A VCF carries no population labels, so they can only come from --pop.  The
//file is applied to the sample IDs on the #CHROM line; pass DEFAULT_POP (or
//an empty string) for the pooled single-column output.
void freqOnlyVCF(string vcffile, string outfile, int nresample, bool PASS_ONLY,
                 const string &popfile);

double calcDensity(int numLoci, vector< MapData * > *mapDataByChr, centromere *centro);

void parallelHR2(HR2_work_order_t *p);
void parallelR2(R2_work_order_t *p);
void parallelLDFromBand(BAND_work_order_t *p);
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

//The multi-population counterpart of writeFreqData: one frequency column per
//population, headed by its label.  This is the format readFreqData accepts, so
//a run's own frequency file can be fed back with --freq-file to reproduce it.
void writeFreqDataWide(string freqOutfile,
                       const vector< vector< FreqData * >* > &freqByPop,
                       const vector<string> &popNames,
                       vector< MapData * > *mapDataByChr);

void writeFreqData(string freqOutfile,
                   vector< FreqData * > *freqDataByChr,
                   vector< MapData * > *mapDataByChr,
                   IndData *indData);

//Read a frequency file, returning one set of per-chromosome frequencies for
//each requested population, in the order of `populations`.
//
//The file is WIDE: one row per locus, one frequency column per population,
//with the population names in the header after ALLELE.
//
//    CHR   SNP     POS    ALLELE   YRI    CEU
//    chr1  rs1234  1005   A        0.31   0.62
//
//Wide rather than one row per (locus, population), for two reasons.  ALLELE
//orients the frequency -- readFreqData flips to 1-f when it disagrees with
//mapData->allele -- and the counted allele is a property of the SITE, so one
//ALLELE column makes it impossible for two populations to disagree about the
//orientation of the same locus.  And the reader walks loci in lockstep with
//mapData, one row each, which a long layout would break.
//
//Backward compatible: a file with exactly ONE frequency column applies that
//column to every requested population, whatever the column is called.  That
//is every frequency file garlic has ever written.
//
//`populations` empty means "one pooled set", which is what a pooled analysis
//wants; a file with several columns is then an error, because there is no
//correct way to choose between them.
vector< vector< FreqData * >* > *readFreqData(string freqfile,
                                              vector< MapData * > *mapDataByChr,
                                              const vector<string> &populations);


//indData is used only to name the individual in a diagnostic; the values are
//read positionally as before.
vector< GenoLikeData * > *readTGLSData(string filename,
                                       int expectedLoci,
                                       int expectedInd,
                                       vector< MapData * > *mapDataByChr,
                                       string GL_TYPE,
                                       IndData *indData);

MapData *initMapData(int nloci);

//A deep copy of the per-chromosome map.
//
//Site filtering prunes the map alongside the genotypes, and each population
//filters on its OWN allele frequencies, so each needs a map it can prune
//without disturbing the next population's.  MapData is all vectors, so the
//copy is the vectors copying themselves.
vector< MapData * > *cloneMapData(const vector< MapData * > *mapDataByChr);
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

//Gather an explicit set of individuals into freshly allocated structures.
//
//Everything is deep copied: subsetData used to alias IndData::pop instead,
//which leaked one array and double-freed another (the --auto-winsize abort).
//keepInd may be any subset in any order; the output is in keepInd's order, so
//the caller controls the individual ordering of the result.
//
//This is the gather subsetData has always done, separated from the choice of
//WHICH individuals, so that a population can be selected as easily as a random
//subsample.
void subsetDataByIndex(vector< HapData * > *hapDataByChr,
                       vector< GenoLikeData *> *GLDataByChr,
                       IndData *indData,
                       const vector<int> &keepInd,
                       vector< HapData * > **subsetHapDataByChr,
                       vector< GenoLikeData *> **subsetGLDataByChr,
                       IndData **subsetIndData,
                       bool USE_GL, bool PHASED);

//Random subsample of `subsample` individuals (all of them if subsample >= n).
//A thin wrapper over subsetDataByIndex that draws the indices.
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
//The LOD score cutoff and the window size are estimated from the AUTOSOMES
//and applied to the sex chromosome, so the density these build has to be able
//to exclude a chromosome while the calling step still uses its windows.  That
//cannot be done by writing MISSING into them, which is how ineligible
//INDIVIDUALS are excluded -- hence a filter here rather than another mask.
//
//role may be NULL, which includes every chromosome and is what every caller
//did before there was such a thing as a role.
DoubleData *convertWinData2DoubleData(vector< WinData * > *winDataByChr, int step,
                                      const vector<ChrRole> *role = NULL,
                                      ChrRole keep = CHR_AUTOSOME);
DoubleData *convertSubsetWinData2DoubleData(vector< WinData * > *winDataByChr, IndData *indData,
                                            int subsample, int step,
                                            const vector<ChrRole> *role = NULL,
                                            ChrRole keep = CHR_AUTOSOME);
void releaseDoubleData(DoubleData *data);
void writeDoubleData(vector < DoubleData * > *rawWinDataByPop, vector< MapData * > *mapDataByChr, vector< IndData * > *indDataByPop);

//counts the number of "fields" in a string
//where a field is defined as a contiguous set of non whitespace
//characters and fields are delimited by whitespace
int countFields(const string &str);
string lc(string str);
string checkChrName(string chr);

//One allele frequency from a count of counted alleles and a total of observed
//alleles, with the optional --resample step.
//
//Four call sites computed this identically and independently: loadTPEDData,
//loadVCFData, freqOnly and freqOnlyVCF.  That is the shape the filterMonomorphic*
//family had before it was collapsed -- one rule, re-derived in several places,
//free to drift.  Per-population frequencies add a fifth caller, so it is
//centralised here first.
//
//r may be NULL when nresample <= 0, since no draw is made then.
double alleleFrequency(double nalleles, double total, int nresample, GarlicRNG *r);

//Allele frequencies for a chosen set of individuals, computed from an ALREADY
//LOADED genotype matrix rather than from the file.
//
//This is how a population gets its own frequencies: the population labels are
//not known while the file is being read (readIndData3 runs after
//loadTPEDData, and --pop later still), so the frequencies are recomputed once
//the labels are in hand.
//
//NOTE a half call ("A 0" in a TPED, "0/." in a VCF) counts as MISSING here,
//contributing to neither numerator nor denominator, because that is what the
//genotype matrix holds.  loadVCFData already does exactly this; loadTPEDData
//does NOT -- it counts the present allele of a half call towards the
//frequency while storing the genotype as missing.  So on TPED data containing
//half calls these frequencies differ from loadTPEDData's, and agree with
//loadVCFData's.  The two input paths already disagreed; this follows the VCF
//one.  Neither bundled example contains a missing call of either kind.
vector< FreqData * > *calcFreqDataForIndices(vector< HapData * > *hapDataByChr,
                                             const vector<int> &keepInd,
                                             int nresample);

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

//The distinct population labels, in ORDER OF FIRST APPEARANCE, each paired
//with the number of individuals carrying it.
//
//First appearance rather than sorted: it is the order the input presents, it
//is stable under relabelling, and it is the order in which populations will be
//analysed once they are analysed separately -- which also fixes the order the
//per-population RNG seeds are derived in, so a run stays reproducible.
vector< pair<string, int> > enumeratePopulations(IndData *indData);

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
