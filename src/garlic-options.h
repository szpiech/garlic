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
    char            TPED_MISSING;
    string          GL_TYPE;
    bool            WEIGHTED;
    string          mapfile;
    bool            CM;
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
    //The seed actually used: --seed if given, otherwise the one drawn.  Recorded
    //in <out>.params.json so a replay reuses it.
    unsigned long int SEED;
};

//Return values of configureFromCommandLine.
const int OPTIONS_OK = 0;           //carry on
const int OPTIONS_USAGE_ERROR = 1;  //a flag was missing, malformed or contradictory
const int OPTIONS_DONE = 2;         //the run finished inside here (--freq-only)

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
