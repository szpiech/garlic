#ifndef __GARLIC_FEATURES_H__
#define __GARLIC_FEATURES_H__

//Counting classified genotypes inside and outside runs of homozygosity.
//
//This replaces count_features_in_roh.pl and count_damaging_genotypes_in_roh.pl.
//Two things made the Perl unsafe rather than merely slow, and both are
//structural rather than bugs to be patched:
//
//  - an individual absent from the .roh.bed autovivified an empty interval
//    list, so EVERY homozygote was silently reported as outside ROH.  A sample
//    ID mismatch therefore produced a complete, plausible, wrong table.  Here
//    the individuals come from the run that made the calls, and a counting
//    file whose IDs do not match is reported (see ROHIndex::indexOf and its
//    callers).
//  - the size classes were hardcoded A/B/C/NONE, while --nclust is unbounded
//    and sizeClassLabel() runs A..Z,AA..  Calls in a fourth class were counted
//    into no column at all.  Here the buckets come from the boundaries the
//    population was actually analysed with.
//
//Nothing here reads garlic's genotype matrix.  By the time ROH exist that
//matrix has been filtered for monomorphic sites (per population), pruned by
//--chr and --autosomes-only, had excluded regions dropped and had hemizygous
//genotypes recoded -- and it has been released.  Counting from it would drop
//exactly the rare functional sites this analysis is about.  The counting pass
//streams a genotype file and reads raw calls.

#include <string>
#include <vector>
#include <map>
#include "garlic-pos.h"
#include "garlic-data.h"
#include "garlic-roh.h"

using namespace std;

//---- the classification table ----------------------------------------------
//
//A feature file is four whitespace-separated columns:
//
//    #chr   pos        allele  class
//    chr21  15481365   T       probably_damaging
//    chr21  15481365   T       nonsynonymous
//    21     15502312   A       benign
//
//'#' comments and blank lines are skipped.  The chromosome is matched through
//canonChrKey(), so 21, chr21 and Chr21 all resolve to the same chromosome.
//The allele is compared as a string, case-insensitively, so a VCF indel allele
//works and a single-character TPED allele compares the same way.
//
//A site may carry SEVERAL labels: two annotation schemes can classify the same
//variant, and the categories are not required to be mutually exclusive.  An
//exactly repeated (chr, pos, allele, class) is an error rather than the Perl's
//silent last-one-wins.

struct FeatureAllele
{
    string allele;   //upper-cased at read time; compare against upper-cased data
    int    cls;      //index into FeatureTable::classes()
};

class FeatureTable
{
public:
    FeatureTable() : rows(0) {}

    //Returns 0, or -1 after logging.
    int read(const string &filename);

    //The sites of one chromosome, or NULL.  The streaming readers hold this
    //across a chromosome rather than looking the name up once per record.
    const map<pos_t, vector<FeatureAllele> > *chromosome(const string &chrKey) const;

    const vector<string> &classes() const { return classNames; }
    int       nclass() const { return int(classNames.size()); }
    long long nrows()  const { return rows; }
    long long nsites() const;
    bool      empty()  const { return byChr.empty(); }

private:
    map<string, map<pos_t, vector<FeatureAllele> > > byChr;
    vector<string> classNames;      //in first-appearance order; sorted by read()
    map<string, int> classIndex;
    long long rows;
};

//---- where a site falls, for one individual --------------------------------
//
//Non-negative values are size class indices, as rohSizeClassIndex() produces
//them.  The three negative codes are distinct on purpose: folding UNASSESSED
//into NONE is what inflates an outside-ROH count by however much of the genome
//the user excluded from the analysis.

const int BUCKET_NONE       = -1;  //inside the callable span, not in a run
const int BUCKET_UNASSESSED = -2;  //no run could have been called here at all
const int BUCKET_INELIGIBLE = -3;  //this individual is not diploid here

//---- the calls, indexed for lookup -----------------------------------------
//
//Built from the ROHData of each population as that population finishes, so the
//counting file is streamed ONCE for the whole run however many populations it
//holds.  The tracts are small (a few hundred per individual); the counting
//file may be a whole-genome VCF.

class ROHIndex
{
public:
    //One individual's runs on one chromosome.  Parallel vectors, sorted by
    //start and disjoint, which is how assembleROHWindows emits them.
    struct ChrTracts
    {
        vector<pos_t> start, stop;
        vector<int>   cls;
    };

    //What was analysed on one chromosome, for one population.  Everything
    //needed to tell NONE from UNASSESSED without consulting the genotype data
    //again.
    struct ChrInfo
    {
        pos_t lo, hi;          //first and last analysed marker
        pos_t gapLo, gapHi;    //assembly gap (centromere); 0,0 when none
        ChrRole role;
        vector<Interval> excluded;   //--par and friends, inclusive at both ends
        ChrInfo() : lo(0), hi(0), gapLo(0), gapHi(0), role(CHR_AUTOSOME) {}
    };

    //Called once per population, after assembleROHWindows and before
    //releaseROHData.  mapDataByChr, bounds and indData are that population's
    //own: filtering, and therefore the analysed span, is per population.
    //Returns 0, or -1 after logging (a duplicate individual ID across
    //populations, which nothing else would catch).
    int addPopulation(const string &popLabel,
                      const vector< ROHData * > *rohDataByInd,
                      const vector< MapData * > *mapDataByChr,
                      const vector<double> &bounds,
                      const IndData *indData,
                      centromere *centro,
                      const ExcludedRegions *excluded,
                      const vector<ChrRole> *role);

    int  nind() const { return int(inds.size()); }
    int  npop() const { return int(pops.size()); }
    //-1 when the ID is not in the index.  This is the check the Perl never
    //made.
    int  indexOf(const string &id) const;
    const string &indID(int i)    const { return inds[i].id; }
    int  popOf(int i)             const { return inds[i].pop; }
    int  zygoOf(int i)            const { return inds[i].zygo; }
    const string &popName(int p)  const { return pops[p].label; }
    const vector<double> &bounds(int p) const { return pops[p].bounds; }
    //Stored rather than derived from the boundaries: when the calls come
    //from a .roh.bed the boundaries are not in the file and the class count
    //is whatever letters it used.
    int  nclass(int p)            const { return pops[p].nclassN; }
    //True when this index was built from a .roh.bed rather than from calls
    //made in this run, which limits what the table can say.
    bool fromBed()                const { return builtFromBed; }
    //Individuals of one population, in the order they were added.
    const vector<int> &indsOf(int p) const { return pops[p].members; }

    //Resolved once per chromosome by the streaming readers rather than once
    //per genotype: a map lookup per (individual, site) is 10^8 lookups on a
    //whole-genome file.  Both vectors are sized nind(); an entry is NULL where
    //that individual has nothing on this chromosome.
    void resolveChr(const string &chrKey,
                    vector<const ChrTracts *> &tracts,
                    vector<const ChrInfo *> &info) const;

    //Builds the index from a .roh.bed instead of from a run's own calls, for
    //counting against calls that already exist.  Strictly less is knowable
    //this way and the table says so: a .roh.bed records where the runs are,
    //not which chromosomes were analysed, how far the markers reached, or
    //what the size class letters mean -- so nothing can be UNASSESSED and
    //every site outside a run is NONE.
    //
    //coordNote comes back describing which coordinate convention was used;
    //it goes in the output header.  Returns 0, or -1 after logging.
    int addFromBed(const string &bedfile, string &coordNote);

    //One of the BUCKET_ codes, or a size class index.  Static and taking only
    //what it reads, so it is testable without building an index.
    static int bucketAt(const ChrTracts *t, const ChrInfo *ci, int zygo, pos_t pos);

    //Inverse of sizeClassLabel: A->0, Z->25, AA->26.  -1 when the label is
    //not one garlic would have written.
    static int sizeClassIndexFromLabel(const string &label);

    //The individual's own population label, which is what goes in the table's
    //pop column.  Not the same string as popName(): that one is empty for a
    //single-population run, because it also names the output file.
    const string &popDisplay(int i) const { return inds[i].popDisplay; }

private:
    struct IndEntry
    {
        string id;
        string popDisplay;
        int    pop;
        int    zygo;
        map<string, ChrTracts> byChr;
        IndEntry() : pop(0), zygo(ZYG_UNKNOWN) {}
    };

    struct PopEntry
    {
        string label;
        vector<double> bounds;
        int nclassN;
        map<string, ChrInfo> chrInfo;
        vector<int> members;
        PopEntry() : nclassN(1) {}
    };

    vector<IndEntry> inds;
    vector<PopEntry> pops;
    map<string, int> indOf;
    bool builtFromBed;

public:
    ROHIndex() : builtFromBed(false) {}
};

//---- the counts ------------------------------------------------------------
//
//Three numbers per (individual, class, bucket).  hom is what the Perl
//counted; n is the one that was missing, and without it "more deleterious
//homozygotes inside runs" cannot be told apart from "more sites inside runs".
//
//Bucket columns are 0..nclass-1 for the size classes, then NONE, then
//UNASSESSED.  ALL is their sum and is produced by the writer rather than
//stored.

struct FeatureCounts
{
    //[individual][class][bucket]
    vector< vector< vector<long long> > > hom, het, n;
    //Per individual, counted per classified ROW so that
    //  sum(n) + missing + hemi == rowsSeen
    //holds for every individual that the counting file carried.
    vector<long long> missing, hemi;
    //Whether the counting file had a column for this individual at all.
    vector<bool> seen;

    //Run-level accounting.
    long long rowsSeen;        //classified rows found in the genotype file
    long long rowsUnusable;    //rows garlic could not evaluate at their site
    long long sitesNotFound;   //feature sites the genotype file never reached

    int nclass;
    int nbucket;
    int colNone() const       { return nclass; }
    int colUnassessed() const { return nclass + 1; }

    void init(int nind, int nclasses, int nsizeclass);
};

//Streams a TPED, counting classified genotypes against the calls in index.
//Reads raw calls: no site filter, no frequency estimation, no recoding.
//Returns 0, or -1 after logging.
int countFeaturesTPED(const string &tpedfile,
                      const string &tfamfile,
                      char TPED_MISSING,
                      const FeatureTable &features,
                      const ROHIndex &index,
                      FeatureCounts &counts);

//Streams a VCF, counting classified genotypes against the calls in index.
//Sample names come from the #CHROM line, so there is no TFAM.
//
//Unlike the reader that CALLS runs of homozygosity, this one does not skip
//indels or multiallelic sites: a classified variant is frequently one or the
//other, and nothing here depends on the site being biallelic.  The classified
//allele is matched against REF and every ALT, and a row naming an allele that
//is at neither is reported rather than counted as absent.
//Returns 0, or -1 after logging.
int countFeaturesVCF(const string &vcffile,
                     bool PASS_ONLY,
                     const FeatureTable &features,
                     const ROHIndex &index,
                     FeatureCounts &counts);

//What the run did, for the table's header.  A counts file is not
//self-describing without it: which classification, which genotypes, and --
//when the calls came from a file rather than from this run -- which file and
//how its coordinates were read.
struct FeatureRunInfo
{
    string featureFile;
    string genotypeSource;
    string rohSource;    //empty when the calls were made by this run
    string coordNote;    //how a .roh.bed's coordinate convention was settled
    bool   pooled;
    FeatureRunInfo() : pooled(false) {}
};

//One table per population, named like that population's other outputs.
//outfile is the full path.  Returns 0, or -1 after logging.
int writeFeatureCounts(const string &outfile,
                       const FeatureTable &features,
                       const ROHIndex &index,
                       int pop,
                       const FeatureCounts &counts,
                       const FeatureRunInfo &runInfo);

#endif
