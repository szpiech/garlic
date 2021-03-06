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

/*
extern const string ARG_FEATURE_TPED;
extern const string DEFAULT_FEATURE_TPED;
extern const string HELP_FEATURE_TPED;

extern const string ARG_FEATURE_TFAM;
extern const string DEFAULT_FEATURE_TFAM;
extern const string HELP_FEATURE_TFAM;

extern const string ARG_FEATURES;
extern const string DEFAULT_FEATURES;
extern const string HELP_FEATURES;
*/

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


param_t *getCLI(int argc, char *argv[]);
bool checkBuild(string BUILD);
bool checkBuildAndCentromereFile(string BUILD, string centromereFile);
bool checkMultiWinsizes(vector<int> &multiWinsizes, bool &WINSIZE_EXPLORE);
bool checkAutoFreq(string freqfile, bool FREQ_ONLY, bool &AUTO_FREQ);
bool checkAutoWinsizeStep(int auto_winsize_step);
bool checkAutoWinsize(bool WINSIZE_EXPLORE, bool AUTO_WINSIZE);
bool checkAutoCutoff(double LOD_CUTOFF, bool &AUTO_CUTOFF);
bool checkBoundSizes(vector<double> &boundSizes, bool &AUTO_BOUNDS);
bool checkRequiredFiles(string tpedfile, string tfamfile);
bool checkMapFile(string mapfile, bool WEIGHTED);
bool checkThreads(int numThreads);
bool checkError(double error, string tglsfile);
bool checkGLType(string TYPE, string tglsfile);
bool checkWinsize(int winsize, bool WINSIZE_EXPLORE, bool AUTO_WINSIZE, bool WEIGHTED);
bool checkMaxGap(int MAX_GAP);
bool checkOverlapFrac(double OVERLAP_FRAC);
bool checkM(int M);
bool checkMU(double mu);
bool checkNCLUST(int nclust);
bool checkCM(string mapfile, bool CM);
#endif