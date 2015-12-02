#ifndef __GARLIC_CLI_H__
#define __GARLIC_CLI_H__

#include <string>
#include "param_t.h"

using namespace std;

const string VERSION = "1.0.0";

const string ARG_POP_SPLIT = "--split-pops";
const bool DEFAULT_POP_SPLIT = false;
const string HELP_POP_SPLIT = "Read in data, calculate allele frequencies, and output\n\
data in multiple files (one for each population).";

const string ARG_OUTFILE = "--out";
const string DEFAULT_OUTFILE = "outfile";
const string HELP_OUTFILE = "The base name for all output files.";

const string ARG_THREADS = "--threads";
const int DEFAULT_THREADS = 1;
const string HELP_THREADS = "The number of threads to spawn during calculations.";

const string ARG_ERROR = "--error";
const double DEFAULT_ERROR = 0.001;
const string HELP_ERROR = "The assumed genotyping error rate.";

const string ARG_WINSIZE = "--winsize";
const int DEFAULT_WINSIZE = 60;
const string HELP_WINSIZE = "The window size in # of SNPs in which to calculate LOD scores.";

const string ARG_WINSIZE_MULTI = "--winsize-multi";
const int DEFAULT_WINSIZE_MULTI = -1;
const string HELP_WINSIZE_MULTI = "Provide several window sizes (in # of SNPs) to calculate LOD scores.\n\
\tLOD score KDEs for each window size will be output for inspection.";

const string ARG_AUTO_WINSIZE = "--auto-winsize";
const bool DEFAULT_AUTO_WINSIZE = 0;
const string HELP_AUTO_WINSIZE = "Initiates an ad hoc method for automatically selecting the # of SNPs in which to\n\
\tcalculate LOD scores. Starts at the value specified by --winsize and increases\n\
\tby 10 SNPs until finished.";

/*
const string ARG_POINTS = "--kde-points";
const int DEFAULT_POINTS = 512;
const string HELP_POINTS = "The number of equally spaced points at which to do KDE.";
*/
/*
const string ARG_BW = "--kde-bw";
const double DEFAULT_BW = -1;
const string HELP_BW = "Manually set the bandwidth for the KDE of lod scores.\n\
\tBy default, the nrd0 rule of thumb is used.";
*/
const string ARG_MAX_GAP = "--max-gap";
const int DEFAULT_MAX_GAP = 200000;
const string HELP_MAX_GAP = "A LOD score window is not calculated if the gap (in bps)\n\
\tbetween two loci is greater than this value.";

const string ARG_RESAMPLE = "--resample";
const int DEFAULT_RESAMPLE = 0;
const string HELP_RESAMPLE = "Number of resamples for estimating allele frequencies.\n\
\tWhen set to 0 (default), garlic will use allele\n\
\tfrequencies as calculated from the data.";

const string ARG_TPED = "--tped";
const string DEFAULT_TPED = "__tpedfile";
const string HELP_TPED = "A tped formatted file containing map and genotype information.";

const string ARG_TFAM = "--tfam";
const string DEFAULT_TFAM = "__tfamfile";
const string HELP_TFAM = "A tfam formatted file containing population and individual IDs.";

const string ARG_RAW_LOD = "--raw-lod";
const bool DEFAULT_RAW_LOD = false;
const string HELP_RAW_LOD = "If set, LOD scores will be output to gzip compressed files.";

const string ARG_LOD_CUTOFF = "--lod-cutoff";
const double DEFAULT_LOD_CUTOFF = -999999;
const string HELP_LOD_CUTOFF = "For LOD based ROH calling, specify a single LOD score cutoff\n\
\tabove which ROH are called in all populations.  By default, this is chosen\n\
\tautomatically per population with KDE.";

/*
const string ARG_LOD_CUTOFF_FILE = "--lod-cutoff-file";
const string DEFAULT_LOD_CUTOFF_FILE = "_none";
const string HELP_LOD_CUTOFF_FILE = "For LOD based ROH calling, specify a file with LOD score cutoffs\n\
\tabove which ROH are called for each population.\n\
\tFile format is <pop ID> <cutoff>.\n\
\tBy default, these cutoffs are chosen automatically per population with KDE.";
*/
const string ARG_BOUND_SIZE = "--size-bounds";
const double DEFAULT_BOUND_SIZE = -1;
const string HELP_BOUND_SIZE = "Specify the short/medium and medium/long\n\
\tROH boundaries.  By default, this is chosen automatically\n\
\twith a 3-component GMM.  Must provide 2 numbers.";
/*
const string ARG_BOUND_SIZE_FILE = "--size-bounds-file";
const string DEFAULT_BOUND_SIZE_FILE = "_none";
const string HELP_BOUND_SIZE_FILE = "A file specifying the short/medium and medium/long\n\
\tROH boundaries per population.\n\
\tFile format <pop ID> <short/medium boundary> <medium/long boundary>\n\
\tBy default, this is chosen automatically\n\
\twith a 3-component GMM.  Must provide 2 numbers.";
*/
const string ARG_TPED_MISSING = "--tped-missing";
const char DEFAULT_TPED_MISSING = '0';
const string HELP_TPED_MISSING = "Missing data code for TPED files.";

const string ARG_FREQ_FILE = "--freq-file";
const string DEFAULT_FREQ_FILE = "_none";
const string HELP_FREQ_FILE = "A file specifying allele frequencies for\n\
\teach population for all variants. File format:\n\
\tSNP\tALLELE\t<pop1 ID> <pop2 ID> ...\n\
\t<locus ID> <allele> <pop1 freq> <pop2 freq> ...\n\
\tBy default, this is calculated automatically\n\
\tfrom the provided data.";

const string ARG_FREQ_ONLY = "--freq-only";
const bool DEFAULT_FREQ_ONLY = false;
const string HELP_FREQ_ONLY = "If set, calculates a freq file from provided data and then exits.";

const string ARG_KDE_SUBSAMPLE = "--kde-subsample";
const int DEFAULT_KDE_SUBSAMPLE = 10;
const string HELP_KDE_SUBSAMPLE = "The number of individuals to randomly sample for LOD score KDE. If there\n\
\tare fewer individuals in the population all are used.\n\
Set <= 0 to use all individuals (may use large amounts of RAM).";

const string ARG_BUILD = "--build";
const string DEFAULT_BUILD = "none";
const string HELP_BUILD = "Choose which genome build to use for centromere locations (hg18, hg19, or hg38).\n";

const string ARG_CENTROMERE_FILE = "--centromere";
const string DEFAULT_CENTROMERE_FILE = "__none";
const string HELP_CENTROMERE_FILE = "Provide custom centromere boundaries. Format <chr> <start> <end>.\n";

param_t *getCLI(int argc, char *argv[]);

#endif