#include "garlic-cli.h"
#include <iostream>
#include <fstream>
#include <utility>
#include <sys/stat.h>
#include <cerrno>

const string VERSION = "1.1.6a";

const string PREAMBLE = "\ngarlic v" + VERSION + " -- a program to call runs of homozygosity in genetic data.\n\
Source code and binaries can be found at <https://www.github.com/szpiech/garlic>.\n\
\n\
Citations:\n\
\n\
A Blant, et al. (2017) bioRxiv, doi: 10.1101/177352\n\
ZA Szpiech, et al. (2017) Bioinformatics, doi: 10.1093/bioinformatics/btx102\n\
TJ Pemberton, et al. (2012) AJHG, 91: 275–292\n";

const string ARG_OVERLAP_FRAC = "--overlap-frac";
const double DEFAULT_OVERLAP_FRAC = 0.25;
const string HELP_OVERLAP_FRAC = "The minimum fraction of overlapping windows above the LOD cutoff required\n\
\tto begin constructing a run. ROH will have a lower bound size threshold of WINSIZE*OVERLAP_FRAC.\n\
\tIf set to 0, GARLIC sets the value to the lowest sensible value: 1/winsize.";

const string ARG_AUTO_OVERLAP_FRAC = "--auto-overlap-frac";
const bool DEFAULT_AUTO_OVERLAP_FRAC = false;
const string HELP_AUTO_OVERLAP_FRAC = "If set, GARLIC attempts to guess based on marker density.";

const string ARG_OUTFILE = "--out";
const string DEFAULT_OUTFILE = "outfile";
const string HELP_OUTFILE = "The base name for all output files.";

const string ARG_THREADS = "--threads";
const int DEFAULT_THREADS = 1;
const string HELP_THREADS = "The number of threads to use. Applies to LOD score calculation, the LD\n\
\tcalculations under --weighted, the LOD score KDE, and ROH assembly.";

const string ARG_ERROR = "--error";
const double DEFAULT_ERROR = -1;
const string HELP_ERROR = "The assumed genotyping error rate.";

const string ARG_WINSIZE = "--winsize";
const int DEFAULT_WINSIZE = 0;
const string HELP_WINSIZE = "The window size in # of SNPs in which to calculate LOD scores.";

const string ARG_WINSIZE_MULTI = "--winsize-multi";
const int DEFAULT_WINSIZE_MULTI = -1;
const string HELP_WINSIZE_MULTI = "Provide several window sizes (in # of SNPs) to calculate LOD scores.\n\
\tLOD score KDEs for each window size will be output for inspection.";

const string ARG_AUTO_WINSIZE = "--auto-winsize";
const bool DEFAULT_AUTO_WINSIZE = false;
const string HELP_AUTO_WINSIZE = "If --weighted is set, guesses the best window size based on SNP density, otherwise\n\
\tinitiates an ad hoc method for automatically selecting the # of SNPs in which to\n\
\tcalculate LOD scores. Starts at the value specified by --winsize and increases\n\
\tby <step size> SNPs until finished.";

const string ARG_AUTO_WINSIZE_STEP = "--auto-winsize-step";
const int DEFAULT_AUTO_WINSIZE_STEP = 10;
const string HELP_AUTO_WINSIZE_STEP = "Step size for automatic window selection algorithm.";

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
const string DEFAULT_TPED = "none";
const string HELP_TPED = "A tped formatted file containing map and genotype information.";

const string ARG_TFAM = "--tfam";
const string DEFAULT_TFAM = "none";
const string HELP_TFAM = "A tfam formatted file containing population and individual IDs.";

const string ARG_TGLS = "--tgls";
const string DEFAULT_TGLS = "none";
const string HELP_TGLS = "A tgls file containing per-genotype likelihoods.";

const string ARG_GL_TYPE = "--gl-type";
const string DEFAULT_GL_TYPE = "none";
const string HELP_GL_TYPE = "Specify the form of the genotype likelihood data: GQ, GL, or PL.\n\
\tGQ is a phred-scaled likelihood of the genotype being incorrect.\n\
\tPL is a phred-scaled likelihood of the genotype being correct.\n\
\tGL is a log10-scaled likelihood of the genotype being correct.";

const string ARG_MAP = "--map";
const string DEFAULT_MAP = "none";
const string HELP_MAP = "Provide a scaffold genetic map, sites that aren't present within this file are interpolated.\n\
\tSites outside the bounds are filtered. This is required for wLOD calculations\n\
\tand any runs for which you wish to report ROH in units of cM.";

const string ARG_WEIGHTED = "--weighted";
const bool DEFAULT_WEIGHTED = false;
const string HELP_WEIGHTED = "Compute LOD scores weighted by LD and probability of mutation.";

const string ARG_RAW_LOD = "--raw-lod";
const bool DEFAULT_RAW_LOD = false;
const string HELP_RAW_LOD = "If set, LOD scores will be output to gzip compressed files.";

const string ARG_LOD_CUTOFF = "--lod-cutoff";
const double DEFAULT_LOD_CUTOFF = -999999;
const string HELP_LOD_CUTOFF = "For LOD based ROH calling, specify a single LOD score cutoff\n\
\tabove which ROH are called in all populations.  By default, this is chosen\n\
\tautomatically with KDE.";

const string ARG_BOUND_SIZE = "--size-bounds";
const double DEFAULT_BOUND_SIZE = -1;
const string HELP_BOUND_SIZE = "Specify the size class boundaries\n\
\tROH boundaries.  By default, this is chosen automatically\n\
\twith a 3-component GMM.  Must provide numbers in increasing order.";

const string ARG_TPED_MISSING = "--tped-missing";
const char DEFAULT_TPED_MISSING = '0';
const string HELP_TPED_MISSING = "Single character missing data code for TPED files.";

const string ARG_FREQ_FILE = "--freq-file";
const string DEFAULT_FREQ_FILE = "none";
const string HELP_FREQ_FILE = "A file specifying allele frequencies for\n\
\teach population for all variants. File format:\n\
\tCHR SNP POS ALLELE FREQ\n\
\t<chr> <locus ID> <allele> <freq>\n\
\tBy default, this is calculated automatically\n\
\tfrom the provided data.";

const string ARG_FREQ_ONLY = "--freq-only";
const bool DEFAULT_FREQ_ONLY = false;
const string HELP_FREQ_ONLY = "If set, calculates a freq file from provided data and then exits. Uses minimal RAM.";

const string ARG_KDE_SUBSAMPLE = "--kde-subsample";
const int DEFAULT_KDE_SUBSAMPLE = 0;
const string HELP_KDE_SUBSAMPLE = "The number of individuals to randomly sample for LOD score KDE. If there\n\
\tare fewer individuals in the population all are used.\n\
\tThe default, 0, uses every individual, which makes the selected LOD cutoff a\n\
\tdeterministic function of the data alone. A positive value subsamples, which\n\
\tsaves memory but makes the cutoff depend on --seed.";

extern const string ARG_LD_SUBSAMPLE = "--ld-subsample";
extern const int DEFAULT_LD_SUBSAMPLE = 0;
extern const string HELP_LD_SUBSAMPLE = "The number of individuals to randomly sample for LD calculation during wLOD. If there\n\
\tare fewer individuals in the population all are used.\n\
Set <= 0 to use all individuals (will increase runtime).";

const string ARG_BUILD = "--build";
const string DEFAULT_BUILD = "none";
const string HELP_BUILD = "Choose which genome build to use for centromere locations (hg18, hg19, or hg38).";

const string ARG_CENTROMERE_FILE = "--centromere";
const string DEFAULT_CENTROMERE_FILE = "none";
const string HELP_CENTROMERE_FILE = "Provide custom centromere boundaries. Format <chr> <start> <end>.";

const string ARG_M = "--M";
const int DEFAULT_M = 7;
const string HELP_M = "The expected number of meioses since a recent common ancestor for --weighted calculation.";

const string ARG_MU = "--mu";
const double DEFAULT_MU = 1e-9;
const string HELP_MU = "Mutation rate per bp per generation for --weighted calculation.";

const string ARG_PHASED = "--phased";
const bool DEFAULT_PHASED = false;
const string HELP_PHASED = "Set if data are phased and you want to calculate r2 instead of hr2 while --weighted is set.\n\
\tUses extra RAM. Has no effect on computations without --weighted.";

const string ARG_NCLUST = "--nclust";
const int DEFAULT_NCLUST = 3;
const string HELP_NCLUST = "Set number of clusters for GMM classification of ROH lengths.";

const string ARG_CM = "--cm";
const bool DEFAULT_CM = false;
const string HELP_CM = "Measure ROH lengths in genetic distance units. This requires a mapfile.";

const string ARG_KDE_THINNING = "--no-kde-thinning";
const bool DEFAULT_KDE_THINNING = false;
const string HELP_KDE_THINNING = "Set this flag to send all LOD score data to KDE function. This may dramatically\n\
\tincrease runtime.";

const string ARG_MAX_WINSIZE = "--max-winsize";
const int DEFAULT_MAX_WINSIZE = 1000;
const string HELP_MAX_WINSIZE = "Upper bound on the window size that --auto-winsize will try. The search\n\
\tpreviously had no bound and would grow past the number of loci if the\n\
\tsmoothness criterion was never met.";

const string ARG_DUMP_DOCS = "--dump-docs";
const string DEFAULT_DUMP_DOCS = "__none";
const string HELP_DUMP_DOCS = "Write the command line reference to stdout in a documentation format and\n\
\texit: txt for the README block, tex for a LaTeX description list. Used by\n\
\t'make docs' so the documentation cannot drift from the program.";

const string ARG_LOAD_PARAMS = "--load-params";
const string DEFAULT_LOAD_PARAMS = "__none";
const string HELP_LOAD_PARAMS = "Read flag values from a <out>.params.json written by a previous run. Flags\n\
\tgiven on the command line take precedence over the file, so a run can be\n\
\trepeated with one parameter changed.";

const string ARG_QUIET = "--quiet";
const bool DEFAULT_QUIET = false;
const string HELP_QUIET = "Suppress the progress bar and informational messages. Errors and warnings\n\
\tstill go to stderr and to <out>.error.";

const string ARG_VERBOSE = "--verbose";
const bool DEFAULT_VERBOSE = false;
const string HELP_VERBOSE = "Show the progress bar even when stderr is not a terminal. By default it is\n\
\tdrawn only on a terminal, because it works by emitting backspaces.";

const string ARG_CHR = "--chr";
const string HELP_CHR = "Analyse only these chromosomes, e.g. --chr chr1 chr2 chrX. Names are matched\n\
\tafter normalisation, so '1' and 'chr1' are the same. It is an error to name a\n\
\tchromosome that is not in the data.\n\tDefault: all chromosomes in the input";

const string ARG_OUTDIR = "--outdir";
const string DEFAULT_OUTDIR = "";
const string HELP_OUTDIR = "Write all output files into this directory, which is created if it does not\n\
\texist. Equivalent to prefixing --out with the path.";

const string ARG_FROH = "--froh";
const bool DEFAULT_FROH = false;
const string HELP_FROH = "Also write <out>.froh.tsv: per-individual autozygous total and fraction,\n\
\tbroken down by ROH size class. The denominator is stated in the file header.";

const string ARG_AUTO_WINSIZE_THRESHOLD = "--auto-winsize-threshold";
const double DEFAULT_AUTO_WINSIZE_THRESHOLD = 0.50;
const string HELP_AUTO_WINSIZE_THRESHOLD = "Smoothness criterion at which the --auto-winsize search stops. Lower values\n\
\tdemand a smoother LOD score density and therefore a larger window.";

const string ARG_KDE_POINTS = "--kde-points";
const int DEFAULT_KDE_POINTS = 512;
const string HELP_KDE_POINTS = "Number of grid points at which the LOD score density is evaluated. Also sets\n\
\tthe resolution at which the between-mode minimum (the LOD cutoff) is located.";

const string ARG_KDE_CUT = "--kde-cut";
const double DEFAULT_KDE_CUT = 3;
const string HELP_KDE_CUT = "Extend the density grid this many bandwidths beyond the range of the data.";

const string ARG_MODE_SPAN = "--mode-smooth-span";
const int DEFAULT_MODE_SPAN = 20;
const string HELP_MODE_SPAN = "Number of adjacent grid points used to smooth the density when locating its\n\
\tmodes. Must be smaller than --kde-points.";

const string ARG_AUTO_WINSIZE_COEF = "--auto-winsize-coef";
const string HELP_AUTO_WINSIZE_COEF = "<slope> <intercept>: coefficients of winsize = slope*log(density) + intercept,\n\
\tused by --auto-winsize with --weighted. The defaults (8.3235 138.0521) are an\n\
\tempirical fit to human SNP-array data; override them for other ascertainment\n\
\tschemes or species.\n\tDefault: 8.3235 138.0521";

const string ARG_AUTO_OVERLAP_COEF = "--auto-overlap-coef";
const string HELP_AUTO_OVERLAP_COEF = "<slope> <intercept>: coefficients of overlap% = slope*log(density) + intercept,\n\
\tused by --auto-overlap-frac. The defaults (6.375 63.888) are an empirical fit\n\
\tto human SNP-array data.\n\tDefault: 6.375 63.888";

const string ARG_GMM_MAX_ITER = "--gmm-max-iter";
const int DEFAULT_GMM_MAX_ITER = 1000;
const string HELP_GMM_MAX_ITER = "Maximum EM iterations when fitting ROH size classes.";

const string ARG_GMM_TOL = "--gmm-tol";
const double DEFAULT_GMM_TOL = 1e-5;
const string HELP_GMM_TOL = "Convergence tolerance for the ROH size-class EM fit.";

const string ARG_KDE_THIN_STEP = "--kde-thin-step";
const int DEFAULT_KDE_THIN_STEP = 0;
const string HELP_KDE_THIN_STEP = "Take every Nth LOD score window when building the KDE. 0 (the default) uses\n\
\tthe window size, which keeps the sampled windows non-overlapping; 1 uses every\n\
\twindow. Replaces --no-kde-thinning, which is equivalent to 1 and is kept as a\n\
\tdeprecated alias.";

const string ARG_VERSION = "--version";
const bool DEFAULT_VERSION = false;
const string HELP_VERSION = "Print the version and exit.";

const string ARG_FORCE = "--force";
const bool DEFAULT_FORCE = false;
const string HELP_FORCE = "Overwrite existing output files. Without this, garlic refuses to clobber an\n\
\texisting <out>.roh.bed.";

const string ARG_SEED = "--seed";
const int DEFAULT_SEED = 0;
const string HELP_SEED = "Random seed, for reproducible runs. Affects --kde-subsample,\n\
\t--ld-subsample and --resample. If 0 (default), a nondeterministic seed is drawn\n\
\tand written to the log file so the run can be repeated exactly.";

/*
const string ARG_FEATURE_TPED = "--tped-counting";
const string DEFAULT_FEATURE_TPED = "_none";
const string HELP_FEATURE_TPED = "A TPED formatted file containing genotypes that are classified in the feature file.\n\
Sites not in the feature file are ignored.";

const string ARG_FEATURE_TFAM = "--tfam-counting";
const string DEFAULT_FEATURE_TFAM = "_none";
const string HELP_FEATURE_TFAM = "A TFAM formatted file containing individuals listed in the corresponding TPED file.\n\
Individuals without ROH calls are ignored.";

const string ARG_FEATURES = "--features";
const string DEFAULT_FEATURES = "_none";
const string HELP_FEATURES = "A feature file giving classifications";
*/


//Written by the Makefile so the stamp cannot go stale.  Absent when building
//outside the Makefile or from a tarball with no .git, hence the guard.
#if defined(__has_include)
#  if __has_include("garlic-version.h")
#    include "garlic-version.h"
#  endif
#endif
#ifndef GARLIC_GIT_SHA
#define GARLIC_GIT_SHA "unknown"
#endif

param_t *getCLI(int argc, char *argv[], int &status)
{
	param_t *params = new param_t;
	params->setPreamble(PREAMBLE);
	params->addFlag(ARG_OVERLAP_FRAC, DEFAULT_OVERLAP_FRAC, "", HELP_OVERLAP_FRAC);
	params->addFlag(ARG_AUTO_OVERLAP_FRAC, DEFAULT_AUTO_OVERLAP_FRAC, "", HELP_AUTO_OVERLAP_FRAC);
	params->addFlag(ARG_OUTFILE, DEFAULT_OUTFILE, "", HELP_OUTFILE);
	params->addFlag(ARG_THREADS, DEFAULT_THREADS, "", HELP_THREADS);
	params->addFlag(ARG_ERROR, DEFAULT_ERROR, "", HELP_ERROR);
	params->addFlag(ARG_WINSIZE, DEFAULT_WINSIZE, "", HELP_WINSIZE);
	params->addFlag(ARG_MAX_GAP, DEFAULT_MAX_GAP, "", HELP_MAX_GAP);
	params->addFlag(ARG_RESAMPLE, DEFAULT_RESAMPLE, "", HELP_RESAMPLE);
	params->addFlag(ARG_TPED, DEFAULT_TPED, "", HELP_TPED);
	params->addFlag(ARG_TFAM, DEFAULT_TFAM, "", HELP_TFAM);
	params->addFlag(ARG_TGLS, DEFAULT_TGLS, "", HELP_TGLS);
	params->addFlag(ARG_GL_TYPE, DEFAULT_GL_TYPE, "", HELP_GL_TYPE);
	params->addFlag(ARG_MAP, DEFAULT_MAP, "", HELP_MAP);
	params->addFlag(ARG_WEIGHTED, DEFAULT_WEIGHTED, "", HELP_WEIGHTED);
	params->addFlag(ARG_RAW_LOD, DEFAULT_RAW_LOD, "", HELP_RAW_LOD);
	params->addListFlag(ARG_BOUND_SIZE, DEFAULT_BOUND_SIZE, "", HELP_BOUND_SIZE);
	params->addFlag(ARG_LOD_CUTOFF, DEFAULT_LOD_CUTOFF, "", HELP_LOD_CUTOFF);
	params->addFlag(ARG_TPED_MISSING, DEFAULT_TPED_MISSING, "", HELP_TPED_MISSING);
	params->addFlag(ARG_FREQ_FILE, DEFAULT_FREQ_FILE, "", HELP_FREQ_FILE);
	params->addFlag(ARG_FREQ_ONLY, DEFAULT_FREQ_ONLY, "", HELP_FREQ_ONLY);
	params->addListFlag(ARG_WINSIZE_MULTI, DEFAULT_WINSIZE_MULTI, "", HELP_WINSIZE_MULTI);
	params->addFlag(ARG_KDE_SUBSAMPLE, DEFAULT_KDE_SUBSAMPLE , "", HELP_KDE_SUBSAMPLE);
	params->addFlag(ARG_LD_SUBSAMPLE, DEFAULT_LD_SUBSAMPLE , "", HELP_LD_SUBSAMPLE);
	params->addFlag(ARG_AUTO_WINSIZE, DEFAULT_AUTO_WINSIZE, "", HELP_AUTO_WINSIZE);
	params->addFlag(ARG_AUTO_WINSIZE_STEP, DEFAULT_AUTO_WINSIZE_STEP, "", HELP_AUTO_WINSIZE_STEP);
	params->addFlag(ARG_BUILD, DEFAULT_BUILD, "", HELP_BUILD);
	params->addFlag(ARG_CENTROMERE_FILE, DEFAULT_CENTROMERE_FILE, "", HELP_CENTROMERE_FILE);
	params->addFlag(ARG_M, DEFAULT_M, "", HELP_M);
	params->addFlag(ARG_MU, DEFAULT_MU, "", HELP_MU);
	params->addFlag(ARG_PHASED, DEFAULT_PHASED, "", HELP_PHASED);
	params->addFlag(ARG_NCLUST, DEFAULT_NCLUST, "", HELP_NCLUST);
	params->addFlag(ARG_CM, DEFAULT_CM, "", HELP_CM);
	params->addFlag(ARG_KDE_THINNING, DEFAULT_KDE_THINNING, "", HELP_KDE_THINNING);
	params->addFlag(ARG_SEED, DEFAULT_SEED, "", HELP_SEED);
	params->addFlag(ARG_MAX_WINSIZE, DEFAULT_MAX_WINSIZE, "", HELP_MAX_WINSIZE);
	params->addFlag(ARG_VERSION, DEFAULT_VERSION, "", HELP_VERSION);
	params->addFlag(ARG_FORCE, DEFAULT_FORCE, "", HELP_FORCE);
	params->addFlag(ARG_KDE_THIN_STEP, DEFAULT_KDE_THIN_STEP, "", HELP_KDE_THIN_STEP);
	params->addFlag(ARG_DUMP_DOCS, DEFAULT_DUMP_DOCS, "", HELP_DUMP_DOCS);
	params->addFlag(ARG_LOAD_PARAMS, DEFAULT_LOAD_PARAMS, "", HELP_LOAD_PARAMS);
	params->addFlag(ARG_QUIET, DEFAULT_QUIET, "", HELP_QUIET);
	params->addFlag(ARG_VERBOSE, DEFAULT_VERBOSE, "", HELP_VERBOSE);
	params->addListFlag(ARG_CHR, "_ALL", "", HELP_CHR);
	params->addFlag(ARG_OUTDIR, DEFAULT_OUTDIR, "", HELP_OUTDIR);
	params->addFlag(ARG_FROH, DEFAULT_FROH, "", HELP_FROH);
	params->addFlag(ARG_AUTO_WINSIZE_THRESHOLD, DEFAULT_AUTO_WINSIZE_THRESHOLD, "", HELP_AUTO_WINSIZE_THRESHOLD);
	params->addFlag(ARG_KDE_POINTS, DEFAULT_KDE_POINTS, "", HELP_KDE_POINTS);
	params->addFlag(ARG_KDE_CUT, DEFAULT_KDE_CUT, "", HELP_KDE_CUT);
	params->addFlag(ARG_MODE_SPAN, DEFAULT_MODE_SPAN, "", HELP_MODE_SPAN);
	params->addListFlag(ARG_AUTO_WINSIZE_COEF, 0.0, "", HELP_AUTO_WINSIZE_COEF);
	params->addListFlag(ARG_AUTO_OVERLAP_COEF, 0.0, "", HELP_AUTO_OVERLAP_COEF);
	params->addFlag(ARG_GMM_MAX_ITER, DEFAULT_GMM_MAX_ITER, "", HELP_GMM_MAX_ITER);
	params->addFlag(ARG_GMM_TOL, DEFAULT_GMM_TOL, "", HELP_GMM_TOL);

	//A bare invocation is a usage error, not a successful no-op run.
	if (argc < 2)
	{
		params->printHelp();
		delete params;
		status = PARAM_ERROR;
		return NULL;
	}

	status = params->parseCommandLine(argc, argv);

	if (status == PARAM_OK && params->isFlagSet(ARG_DUMP_DOCS))
	{
		bool ok = params->writeHelpDoc(cout, params->getStringFlag(ARG_DUMP_DOCS));
		delete params;
		status = ok ? PARAM_HELP : PARAM_ERROR;
		return NULL;
	}

	//Applied after parsing so the command line wins, and before main reads
	//any value.
	if (status == PARAM_OK && params->isFlagSet(ARG_LOAD_PARAMS))
	{
		if (!params->loadFlagsJSON(params->getStringFlag(ARG_LOAD_PARAMS)))
		{
			delete params;
			status = PARAM_ERROR;
			return NULL;
		}
	}

	if (status == PARAM_OK && params->getBoolFlag(ARG_VERSION))
	{
		cout << "garlic v" << VERSION << " (" << GARLIC_GIT_SHA << ")\n";
		delete params;
		status = PARAM_HELP;
		return NULL;
	}

	if (status != PARAM_OK)
	{
		delete params;
		return NULL;
	}
	return params;
}

void writeParamsJSON(string file, param_t *params, vector< pair<string,string> > &resolved)
{
	ofstream out(file.c_str());
	if (out.fail())
	{
		LOG.err("ERROR: Failed to open", file);
		return;
	}
	out << "{\n";
	out << "  \"garlic_version\": \"" << VERSION << "\",\n";
	out << "  \"git_sha\": \"" << GARLIC_GIT_SHA << "\",\n";
	out << "  \"resolved\": {\n";
	for (unsigned int i = 0; i < resolved.size(); i++)
	{
		out << "    \"" << resolved[i].first << "\": " << resolved[i].second;
		if (i + 1 < resolved.size()) out << ",";
		out << "\n";
	}
	out << "  },\n";
	vector<string> given = params->setFlags();
	out << "  \"set\": [";
	for (unsigned int i = 0; i < given.size(); i++)
	{
		if (i) out << ", ";
		out << "\"" << given[i] << "\"";
	}
	out << "],\n";
	out << "  \"flags\": {\n";
	params->writeFlagsJSON(out, "    ");
	out << "  }\n";
	out << "}\n";
	out.close();
	LOG.log("Effective parameters:", file);
	return;
}

//Refuse to clobber a previous run's calls unless asked to.
bool checkOutfileClobber(string outfile, bool force)
{
	if (force) return false;
	string roh = outfile + ".roh.bed";
	ifstream probe(roh.c_str());
	if (probe.good())
	{
		probe.close();
		LOG.err("ERROR: Output file already exists:", roh);
		LOG.err("\tPass --force to overwrite, or choose another --out.");
		return true;
	}
	return false;
}

bool checkSeed(int seed){
	if(seed < 0){
		LOG.err("ERROR: Random seed must be >= 0 (0 means draw one automatically).");
		return true;
	}
	return false;
}

bool checkMaxWinsize(int maxWinsize, int winsize){
	if(maxWinsize < winsize){
		LOG.err("ERROR: --max-winsize must be >= the starting --winsize.");
		return true;
	}
	return false;
}

bool checkCM(string mapfile, bool CM){
	if(CM && mapfile.compare(DEFAULT_MAP) == 0){
		LOG.err("ERROR: Must provide mapfile if you wish to construct ROH in genetic map units.");
		return true;
	}
	else return false;
}
bool checkNCLUST(int nclust){
	if(nclust <= 0){
		LOG.err("ERROR: Must choose positive number for number of GMM clusters.");
		return true;
	}
	else return false;
}

bool checkM(int M){
	if(M <= 0){
		LOG.err("ERROR: M must be an integer > 0.");
		return true;
	}
	else return false;
}

bool checkMU(double mu){
	if(mu <= 0 || mu >= 1){
		LOG.err("ERROR: mu must be between 0 and 1.");
		return true;
	}
	else return false;
}


bool checkBuild(string BUILD)
{
	if (BUILD.compare("hg18") != 0 &&
	        BUILD.compare("hg19") != 0 &&
	        BUILD.compare("hg38") != 0 &&
	        BUILD.compare(DEFAULT_BUILD) != 0) {
		//cerr << "ERROR: Must choose hg18/hg19/hg38/none for build version.\n";
		LOG.err("ERROR: Must choose hg18/hg19/hg38 for build version or provide a custom centromere file.");
		return true;
	}
	return false;
}

bool checkBuildAndCentromereFile(string BUILD, string centromereFile) {
	if (BUILD.compare(DEFAULT_BUILD) == 0 && centromereFile.compare(DEFAULT_CENTROMERE_FILE) == 0) {
		LOG.err("ERROR: Must choose hg18/hg19/hg38 for build version or provide a custom centromere file.");
		return true;
	}
	return false;
}

bool checkMultiWinsizes(vector<int> &multiWinsizes, bool &WINSIZE_EXPLORE, bool wasSet)
{
	if (wasSet)
	{
		for (unsigned int i = 0; i < multiWinsizes.size(); i++)
		{
			if (multiWinsizes[i] <= 0)
			{
				//cerr << "ERROR: SNP window sizes must be > 1.\n";
				LOG.err("ERROR: SNP window sizes must be > 1.");
				return true;
			}
		}
		WINSIZE_EXPLORE = true;
	}
	return false;
}

bool checkAutoFreq(string freqfile, bool FREQ_ONLY, bool &AUTO_FREQ)
{
	if (freqfile.compare(DEFAULT_FREQ_FILE) != 0)
	{
		AUTO_FREQ = false;
		if (FREQ_ONLY)
		{
			//cerr << "ERROR: Specifying both " << ARG_FREQ_ONLY << " and " << ARG_FREQ_FILE << " accomplishes nothing useful.\n";
			LOG.err("ERROR: Specifying both", ARG_FREQ_ONLY, false);
			LOG.err(" and", ARG_FREQ_FILE, false);
			LOG.err(" accomplishes nothing useful.");
			return true;
		}
	}
	return false;
}

bool checkAutoWinsizeStep(int auto_winsize_step) {
	if (auto_winsize_step <= 0) {
		LOG.err("ERROR: Step size for automatic window selection must be positive.");
		return true;
	}
	return false;
}

bool checkAutoWinsize(bool WINSIZE_EXPLORE, bool AUTO_WINSIZE)
{
	//Check if both AUTO_WINSIZE and WINSIZE_EXPLORE are set
	//If so, exit with error.
	if (WINSIZE_EXPLORE && AUTO_WINSIZE)
	{
		//cerr << "ERROR: Must set only one of " << ARG_WINSIZE_MULTI << " and " << ARG_AUTO_WINSIZE << ".\n";
		LOG.err("ERROR: Must set only one of", ARG_WINSIZE_MULTI, false);
		LOG.err(" and", ARG_AUTO_WINSIZE);
		return true;
	}
	return false;
}

bool checkAutoCutoff(double LOD_CUTOFF, bool &AUTO_CUTOFF, bool wasSet)
{
	//Was keyed on LOD_CUTOFF != -999999, so -999999 was an unusable value.
	(void)LOD_CUTOFF;
	if (wasSet) {
		AUTO_CUTOFF = false;
	}
	return false;
}

bool checkBoundSizes(vector<double> &boundSizes, bool &AUTO_BOUNDS, bool wasSet){

	if(!wasSet){
		return false;
	}
	else {
		AUTO_BOUNDS = false;

		for(unsigned int i = 0; i < boundSizes.size(); i++){		
			if (boundSizes[i] <= 0){
				LOG.err("ERROR: User provided size boundaries must be positive.");
				return true;
			}
			if(i > 0){
				if (boundSizes[i] <= boundSizes[i-1]){
					LOG.err("ERROR: User provided size boundaries must be in strictly increasing order.");
					return true;
				}
			}
		}
	}
	return false;
}

bool checkRequiredFiles(string tpedfile, string tfamfile)
{
	if (tpedfile.compare(DEFAULT_TPED) == 0 || tfamfile.compare(DEFAULT_TFAM) == 0)
	{
		//cerr << "ERROR: Must provide both a tped and a tfam file.\n";
		LOG.err("ERROR: Must provide both a tped and a tfam file.");
		return true;
	}
	return false;
}

bool checkMapFile(string mapfile, bool WEIGHTED){
	if(mapfile.compare(DEFAULT_MAP) == 0 && WEIGHTED){
		LOG.err("ERROR: Weighted LOD score method requires a map file.");
		return true;
	}
	return false;
}

bool checkThreads(int numThreads)
{
	if (numThreads <= 0)
	{
		//cerr << "ERROR: Number of threads must be > 0.\n";
		LOG.err("ERROR: Number of threads must be > 0.");
		return true;
	}
	return false;
}

bool checkError(double error, string tglsfile, bool wasSet)
{
	if (!wasSet)
	{
		if (tglsfile.compare(DEFAULT_TGLS) == 0) {
			LOG.err("ERROR: --error must be given, or a TGLS file must be provided.");
			return true;
		}
		return false;
	}
	if (error <= 0 || error >= 1)
	{
		LOG.err("ERROR: Genotype error rate must be > 0 and < 1.");
		return true;
	}
	return false;
}

//mkdir -p, so --outdir can name a nested path.
bool makeOutdir(string dir){
	if(dir.empty()) return false;
	string partial;
	for(unsigned int i = 0; i < dir.size(); i++){
		partial += dir[i];
		if(dir[i] == '/' || i + 1 == dir.size()){
			if(partial == "/" || partial == "./" || partial == "../") continue;
			string p = partial;
			if(p.size() > 1 && p[p.size()-1] == '/') p.erase(p.size()-1);
			if(mkdir(p.c_str(), 0777) != 0 && errno != EEXIST){
				LOG.err("ERROR: Could not create output directory:", p);
				return true;
			}
		}
	}
	return false;
}

bool checkAutoWinsizeThreshold(double t){
	if(t <= 0){ LOG.err("ERROR: --auto-winsize-threshold must be > 0."); return true; }
	return false;
}

bool checkKDEPoints(int m){
	if(m < 16){ LOG.err("ERROR: --kde-points must be >= 16."); return true; }
	return false;
}

bool checkKDECut(double c){
	if(c <= 0){ LOG.err("ERROR: --kde-cut must be > 0."); return true; }
	return false;
}

bool checkModeSpan(int s){
	if(s < 2){ LOG.err("ERROR: --mode-smooth-span must be >= 2."); return true; }
	return false;
}

bool checkCoefPair(vector<double> &coef, string flag){
	if(coef.size() != 2){
		LOG.err("ERROR: " + flag + " takes exactly two values: <slope> <intercept>.");
		return true;
	}
	return false;
}

bool checkGMMParams(int maxIter, double tol){
	if(maxIter < 1){ LOG.err("ERROR: --gmm-max-iter must be >= 1."); return true; }
	if(tol <= 0){ LOG.err("ERROR: --gmm-tol must be > 0."); return true; }
	return false;
}

bool checkKDEThinStep(int step){
	if(step < 0){
		LOG.err("ERROR: --kde-thin-step must be >= 0 (0 means use the window size).");
		return true;
	}
	return false;
}

//--size-bounds and --nclust are alternatives: supplying bounds skips the GMM
//that --nclust configures.  This used to be silent.
void warnBoundsOverridesNclust(bool boundsSet, bool nclustSet){
	if(boundsSet && nclustSet){
		LOG.err("WARNING: --size-bounds was given, so --nclust is ignored.");
		LOG.err("\tThe GMM that --nclust configures only runs when size boundaries are");
		LOG.err("\tchosen automatically.");
	}
}

bool checkGLType(string TYPE, string tglsfile)
{
	if ( TYPE.compare("GQ") != 0 && TYPE.compare("GL") != 0 && TYPE.compare("PL") != 0 && tglsfile.compare(DEFAULT_TGLS) != 0 ) {
		LOG.err("ERROR: Must choose GQ/GL/PL for genotype likelihood format or provide a single error rate with --error.");
		return true;
	}
	return false;
}

bool checkWinsize(int winsize, bool WINSIZE_EXPLORE, bool AUTO_WINSIZE, bool WEIGHTED, bool FREQ_ONLY)
{
	//--freq-only never forms a window, so do not make the user invent a
	//--winsize just to get past this check.
	if (FREQ_ONLY) return false;
	if (winsize <= 1){
		if(!WINSIZE_EXPLORE && !(AUTO_WINSIZE && WEIGHTED)){
			LOG.err("ERROR: SNP window size must be > 1. If using --auto-winsize, this is the starting value.");
			return true;
		}
	}
	return false;
}

bool checkMaxGap(int MAX_GAP)
{
	if (MAX_GAP < 0){
		LOG.err("ERROR: Max gap must be > 0.");
		return true;
	}
	else if (MAX_GAP < 1000){
		LOG.err("WARNING: max gap set very low:", MAX_GAP);
	}
	return false;
}

bool checkOverlapFrac(double OVERLAP_FRAC){
	if(OVERLAP_FRAC < 0 || OVERLAP_FRAC > 1){
		LOG.err("ERROR: Overlap fraction must be >= 0 and <= 1.");
		return true;
	}
	return false;
}