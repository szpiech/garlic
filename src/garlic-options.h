#ifndef __GARLIC_OPTIONS_H__
#define __GARLIC_OPTIONS_H__

#include <string>
#include <vector>
#include "param_t.h"

using namespace std;

//Everything the pipeline needs from the command line, in one place.
//
//main() used to read and validate all 54 flags inline -- 269 lines of
//interleaved reading, validation and side effects -- with the argerr
//accumulator tested after some checks and not others (31 assignments, 27
//tests), so a few failures fell through to the next statement.  Collecting the
//values here makes the validation, and in particular the mutual-exclusion
//rules between --auto-winsize, --winsize-multi and --winsize, reachable from a
//test without running the pipeline.
struct GarlicOptions
{
    string          outfile;
    string          tpedfile;
    string          tfamfile;
    string          tglsfile;
    string          popfile;
    string          vcffile;
    bool            VCF_PASS_ONLY;
    char            TPED_MISSING;
    string          GL_TYPE;
    bool            WEIGHTED;
    string          mapfile;
    bool            CM;
    //--pool-populations.  Read in main, carried here so the per-population
    //analysis can mark its outputs without another parameter.
    bool            POOLED;
    string          BUILD;
    string          centromereFile;
    bool            NO_CENTROMERE;
    int             nresample;
    string          freqfile;
    bool            AUTO_FREQ;
    vector<int>     multiWinsizes;
    bool            WINSIZE_EXPLORE;
    bool            AUTO_WINSIZE;
    int             AUTO_WINSIZE_STEP;
    int             winsize;
    double          LOD_CUTOFF;
    bool            AUTO_CUTOFF;
    vector<double>  boundSizes;
    bool            AUTO_BOUNDS;
    int             numThreads;
    double          error;
    int             MAX_GAP;
    double          OVERLAP_FRAC;
    bool            AUTO_OVERLAP_FRAC;
    double          mu;
    int             M;
    int             NCLUST;
    int             KDE_SUBSAMPLE;
    int             LD_SUBSAMPLE;
    bool            RAW_LOD;
    bool            PHASED;
    int             KDE_THIN_STEP;
    int             MAX_WINSIZE;
    //The shared sex chromosome's own cutoff and window size, when the user
    //gave them.  Unset means "whatever the autosomes use", which is the
    //default because a cutoff estimated on one chromosome's windows from part
    //of a cohort is usually worse than the autosomal one, not better.
    //--froh-denominator, as one of the FROHDenominator values, and the file
    //that can supply the lengths it needs.
    int             FROH_DENOM;
    string          CHR_LENGTHS_FILE;
    double          SEXCHR_LOD_CUTOFF;
    bool            SEXCHR_CUTOFF_SET;
    int             SEXCHR_WINSIZE;
    //The seed actually used: --seed if given, otherwise the one drawn.  Recorded
    //in <out>.params.json so a replay reuses it.
    unsigned long int SEED;
    //--features and the genotypes to count from.  Empty counting paths mean
    //"the file this run was called from", which is the usual case: the
    //classified variants are normally in the same data.
    string          featurefile;
    string          countTpedfile;
    string          countTfamfile;
    string          countVcffile;
};

//Return values of configureFromCommandLine.
const int OPTIONS_OK = 0;           //carry on
const int OPTIONS_USAGE_ERROR = 1;  //a flag was missing, malformed or contradictory
const int OPTIONS_DONE = 2;         //the run finished inside here (--freq-only)
const int OPTIONS_RUNTIME_ERROR = 3;//the command line was fine; something failed

//Reads every flag into opt, validates it, and applies the settings that are
//global rather than per-run (thread counts, RNG seed, the KDE/GMM/auto-winsize
//constants).  Also commits the log, which is held in memory until the command
//line is known to be valid so that a rejected invocation leaves no files.
//
//--freq-only completes its work here and returns OPTIONS_DONE; the caller owns
//`params` and must delete it on every path.
//argc/argv are taken only so the invocation can be written into the log.
int configureFromCommandLine(param_t *params, GarlicOptions &opt, int argc, char *argv[]);

#endif
