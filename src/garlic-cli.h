#ifndef __GARLIC_CLI_H__
#define __GARLIC_CLI_H__

#include <iostream>
#include <string>
#include "param_t.h"
#include "garlic-errlog.h"

using namespace std;

extern const string VERSION;
extern const string PREAMBLE;

extern const string ARG_OVERLAP_FRAC;
extern const double DEFAULT_OVERLAP_FRAC;
extern const string HELP_OVERLAP_FRAC;

extern const string ARG_AUTO_OVERLAP_FRAC;
extern const bool DEFAULT_AUTO_OVERLAP_FRAC;
extern const string HELP_AUTO_OVERLAP_FRAC;

extern const string ARG_OUTFILE;
extern const string DEFAULT_OUTFILE;
extern const string HELP_OUTFILE;

extern const string ARG_THREADS;
extern const int DEFAULT_THREADS;
extern const string HELP_THREADS;

extern const string ARG_ERROR;
extern const double DEFAULT_ERROR;
extern const string HELP_ERROR;

extern const string ARG_WINSIZE;
extern const int DEFAULT_WINSIZE;
extern const string HELP_WINSIZE;

extern const string ARG_WINSIZE_MULTI;
extern const int DEFAULT_WINSIZE_MULTI;
extern const string HELP_WINSIZE_MULTI;

extern const string ARG_AUTO_WINSIZE;
extern const bool DEFAULT_AUTO_WINSIZE;
extern const string HELP_AUTO_WINSIZE;

extern const string ARG_AUTO_WINSIZE_STEP;
extern const int DEFAULT_AUTO_WINSIZE_STEP;
extern const string HELP_AUTO_WINSIZE_STEP;

extern const string ARG_MAX_GAP;
extern const int DEFAULT_MAX_GAP;
extern const string HELP_MAX_GAP;

extern const string ARG_RESAMPLE;
extern const int DEFAULT_RESAMPLE;
extern const string HELP_RESAMPLE;

extern const string ARG_TPED;
extern const string DEFAULT_TPED;
extern const string HELP_TPED;

extern const string ARG_TFAM;
extern const string DEFAULT_TFAM;
extern const string HELP_TFAM;

extern const string ARG_VCF;
extern const string DEFAULT_VCF;
extern const string HELP_VCF;

extern const string ARG_VCF_PASS_ONLY;
extern const bool DEFAULT_VCF_PASS_ONLY;
extern const string HELP_VCF_PASS_ONLY;

extern const string ARG_POP;
extern const string DEFAULT_POP;
extern const string HELP_POP;

extern const string ARG_TGLS;
extern const string DEFAULT_TGLS;
extern const string HELP_TGLS;

extern const string ARG_GL_TYPE;
extern const string DEFAULT_GL_TYPE;
extern const string HELP_GL_TYPE;

extern const string ARG_MAP;
extern const string DEFAULT_MAP;
extern const string HELP_MAP;

extern const string ARG_WEIGHTED;
extern const bool DEFAULT_WEIGHTED;
extern const string HELP_WEIGHTED;

extern const string ARG_RAW_LOD;
extern const bool DEFAULT_RAW_LOD;
extern const string HELP_RAW_LOD;

extern const string ARG_LOD_CUTOFF;
extern const double DEFAULT_LOD_CUTOFF;
extern const string HELP_LOD_CUTOFF;

extern const string ARG_BOUND_SIZE;
extern const double DEFAULT_BOUND_SIZE;
extern const string HELP_BOUND_SIZE;

extern const string ARG_TPED_MISSING;
extern const char DEFAULT_TPED_MISSING;
extern const string HELP_TPED_MISSING;

extern const string ARG_FREQ_FILE;
extern const string DEFAULT_FREQ_FILE;
extern const string HELP_FREQ_FILE;

extern const string ARG_FREQ_ONLY;
extern const bool DEFAULT_FREQ_ONLY;
extern const string HELP_FREQ_ONLY;

extern const string ARG_KDE_SUBSAMPLE;
extern const int DEFAULT_KDE_SUBSAMPLE;
extern const string HELP_KDE_SUBSAMPLE;

extern const string ARG_LD_SUBSAMPLE;
extern const int DEFAULT_LD_SUBSAMPLE;
extern const string HELP_LD_SUBSAMPLE;

extern const string ARG_BUILD;
extern const string DEFAULT_BUILD;
extern const string HELP_BUILD;

extern const string ARG_CENTROMERE_FILE;
extern const string DEFAULT_CENTROMERE_FILE;
extern const string HELP_CENTROMERE_FILE;

extern const string ARG_FEATURES;
extern const string DEFAULT_FEATURES;
extern const string HELP_FEATURES;

extern const string ARG_FEATURE_TPED;
extern const string DEFAULT_FEATURE_TPED;
extern const string HELP_FEATURE_TPED;

extern const string ARG_FEATURE_TFAM;
extern const string DEFAULT_FEATURE_TFAM;
extern const string HELP_FEATURE_TFAM;

extern const string ARG_M;
extern const int DEFAULT_M;
extern const string HELP_M;

extern const string ARG_MU;
extern const double DEFAULT_MU;
extern const string HELP_MU;

extern const string ARG_PHASED;
extern const bool DEFAULT_PHASED;
extern const string HELP_PHASED;

extern const string ARG_NCLUST;
extern const int DEFAULT_NCLUST;
extern const string HELP_NCLUST;

extern const string ARG_CM;
extern const bool DEFAULT_CM;
extern const string HELP_CM;

extern const string ARG_KDE_THINNING;
extern const bool DEFAULT_KDE_THINNING;
extern const string HELP_KDE_THINNING;

extern const string ARG_MAX_WINSIZE;
extern const int DEFAULT_MAX_WINSIZE;
extern const string HELP_MAX_WINSIZE;

extern const string ARG_DUMP_DOCS;
extern const string DEFAULT_DUMP_DOCS;
extern const string HELP_DUMP_DOCS;

extern const string ARG_LOAD_PARAMS;
extern const string DEFAULT_LOAD_PARAMS;
extern const string HELP_LOAD_PARAMS;

extern const string ARG_QUIET;
extern const bool DEFAULT_QUIET;
extern const string HELP_QUIET;

extern const string ARG_VERBOSE;
extern const bool DEFAULT_VERBOSE;
extern const string HELP_VERBOSE;

extern const string ARG_NO_CENTROMERE;
extern const bool DEFAULT_NO_CENTROMERE;
extern const string HELP_NO_CENTROMERE;

extern const string ARG_AUTOSOMES_ONLY;
extern const bool DEFAULT_AUTOSOMES_ONLY;
extern const string HELP_AUTOSOMES_ONLY;

extern const string ARG_SEX_SYSTEM;
extern const string DEFAULT_SEX_SYSTEM;
extern const string HELP_SEX_SYSTEM;

extern const string ARG_SEX_CHR;
extern const string HELP_SEX_CHR;

extern const string ARG_SEX_CHR_DEGENERATE;
extern const string HELP_SEX_CHR_DEGENERATE;

extern const string ARG_HAPLOID_CHR;
extern const string HELP_HAPLOID_CHR;

extern const string ARG_HET_RATE_BOUNDS;
extern const string HELP_HET_RATE_BOUNDS;

extern const string ARG_FROH_DENOM;
extern const string DEFAULT_FROH_DENOM;
extern const string HELP_FROH_DENOM;

extern const string ARG_CHR_LENGTHS;
extern const string DEFAULT_CHR_LENGTHS;
extern const string HELP_CHR_LENGTHS;

extern const string ARG_SEXCHR_LOD_CUTOFF;
extern const double DEFAULT_SEXCHR_LOD_CUTOFF;
extern const string HELP_SEXCHR_LOD_CUTOFF;

extern const string ARG_SEXCHR_WINSIZE;
extern const int DEFAULT_SEXCHR_WINSIZE;
extern const string HELP_SEXCHR_WINSIZE;

extern const string ARG_PAR;
extern const string HELP_PAR;

extern const string ARG_PAR_FILE;
extern const string DEFAULT_PAR_FILE;
extern const string HELP_PAR_FILE;

extern const string ARG_CHR;
extern const string HELP_CHR;

extern const string ARG_OUTDIR;
extern const string DEFAULT_OUTDIR;
extern const string HELP_OUTDIR;

extern const string ARG_POOL;
extern const bool DEFAULT_POOL;
extern const string HELP_POOL;

extern const string ARG_FROH;
extern const bool DEFAULT_FROH;
extern const string HELP_FROH;

extern const string ARG_AUTO_WINSIZE_THRESHOLD;
extern const double DEFAULT_AUTO_WINSIZE_THRESHOLD;
extern const string HELP_AUTO_WINSIZE_THRESHOLD;

extern const string ARG_KDE_POINTS;
extern const int DEFAULT_KDE_POINTS;
extern const string HELP_KDE_POINTS;

extern const string ARG_KDE_CUT;
extern const double DEFAULT_KDE_CUT;
extern const string HELP_KDE_CUT;

extern const string ARG_MODE_SPAN;
extern const int DEFAULT_MODE_SPAN;
extern const string HELP_MODE_SPAN;

extern const string ARG_AUTO_WINSIZE_COEF;
extern const string HELP_AUTO_WINSIZE_COEF;

extern const string ARG_AUTO_OVERLAP_COEF;
extern const string HELP_AUTO_OVERLAP_COEF;

extern const string ARG_GMM_MAX_ITER;
extern const int DEFAULT_GMM_MAX_ITER;
extern const string HELP_GMM_MAX_ITER;

extern const string ARG_GMM_TOL;
extern const double DEFAULT_GMM_TOL;
extern const string HELP_GMM_TOL;

extern const string ARG_KDE_THIN_STEP;
extern const int DEFAULT_KDE_THIN_STEP;
extern const string HELP_KDE_THIN_STEP;

extern const string ARG_VERSION;
extern const bool DEFAULT_VERSION;
extern const string HELP_VERSION;

extern const string ARG_FORCE;
extern const bool DEFAULT_FORCE;
extern const string HELP_FORCE;

extern const string ARG_SEED;
extern const int DEFAULT_SEED;
extern const string HELP_SEED;


//status is set to one of the PARAM_* codes in param_t.h; on PARAM_OK the
//returned pointer is valid, otherwise it is NULL.
param_t *getCLI(int argc, char *argv[], int &status);
bool checkOutfileClobber(string outfile, bool force);
bool checkSeed(int seed);
bool checkMaxWinsize(int maxWinsize, int winsize);
bool checkPopFile(string popfile, string tpedfile, string vcffile);
bool checkBuild(string BUILD);
bool checkBuildAndCentromereFile(string BUILD, string centromereFile, bool NO_CENTROMERE);
bool checkMultiWinsizes(vector<int> &multiWinsizes, bool &WINSIZE_EXPLORE, bool wasSet);
bool checkAutoFreq(string freqfile, bool FREQ_ONLY, bool &AUTO_FREQ);
bool checkAutoWinsizeStep(int auto_winsize_step);
bool checkAutoWinsize(bool WINSIZE_EXPLORE, bool AUTO_WINSIZE);
bool checkAutoCutoff(double LOD_CUTOFF, bool &AUTO_CUTOFF, bool wasSet);

//--sexchr-lod-cutoff and --sexchr-winsize.  A cutoff is a statement about
//windows of a given size, so a sex chromosome given its own size must be given
//its own cutoff too; the autosomal one no longer means anything there.
//--froh-denominator, and the file that can feed it.  Sets denom to one of the
//FROHDenominator values.  Returns true after logging when the combination
//cannot be honoured.
bool checkFROHDenominator(const string &mode, int &denom, bool CM, bool FROH,
                          bool chrLengthsSet);

bool checkSexChrEstimates(double cutoff, bool cutoffSet, int winsize, bool winsizeSet,
                          bool WINSIZE_EXPLORE, bool FREQ_ONLY);
bool checkBoundSizes(vector<double> &boundSizes, bool &AUTO_BOUNDS, bool wasSet);
bool checkRequiredFiles(string tpedfile, string tfamfile, string vcffile, string tglsfile);
bool checkMapFile(string mapfile, bool WEIGHTED);
bool checkThreads(int numThreads);
bool checkError(double error, string tglsfile, bool wasSet, bool haveVCFLikelihoods);
bool checkKDEThinStep(int step);
bool checkAutoWinsizeThreshold(double t);
bool checkKDEPoints(int m);
bool checkKDECut(double c);
bool checkModeSpan(int s);
bool checkCoefPair(vector<double> &coef, string flag);
bool checkGMMParams(int maxIter, double tol);
bool makeOutdir(string dir);
//Writes <out>.params.json: the effective value of every flag plus the values
//resolved during the run (auto-selected window size, LOD cutoff, size class
//boundaries, seed).  'resolved' holds pre-formatted JSON values.
void writeParamsJSON(string file, param_t *params, vector< pair<string,string> > &resolved);
void warnBoundsOverridesNclust(bool boundsSet, bool nclustSet);
bool checkGLType(string TYPE, string tglsfile);
bool checkWinsize(int winsize, bool WINSIZE_EXPLORE, bool AUTO_WINSIZE, bool WEIGHTED, bool FREQ_ONLY);
bool checkMaxGap(int MAX_GAP);
bool checkOverlapFrac(double OVERLAP_FRAC);
bool checkM(int M);
bool checkMU(double mu);
bool checkNCLUST(int nclust);
bool checkCM(string mapfile, bool CM);
#endif