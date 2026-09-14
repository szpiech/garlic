#include "garlic-options.h"
#include "garlic-cli.h"
#include "garlic-errlog.h"
#include "garlic-data.h"
#include "garlic-roh.h"
#include "garlic-kde.h"
#include "garlic-centromeres.h"
#include <iostream>
#include <streambuf>
#include <unistd.h>

string getCommandLineString(int argc, char *argv[])
{
    string str = argv[0];
    for (int i = 1; i < argc; i++) {
        str += " " + string(argv[i]);
    }
    return str;
}

//Discards whatever is written to it; used to silence stdout under --quiet
//without touching every cout site.
class NullBuf : public std::streambuf
{
protected:
    int overflow(int c) { return c; }
};
static NullBuf GARLIC_NULLBUF;

int configureFromCommandLine(param_t *params, GarlicOptions &opt, int argc, char *argv[])
{

    opt.outfile = params->getStringFlag(ARG_OUTFILE);
    bool QUIET = params->getBoolFlag(ARG_QUIET);
    bool VERBOSE = params->getBoolFlag(ARG_VERBOSE);
    if (QUIET && VERBOSE) {
        LOG.err("ERROR: --quiet and --verbose are mutually exclusive.");
        return OPTIONS_USAGE_ERROR;
    }
    LOG.setVerbosity(QUIET, VERBOSE);
    //The bar writes backspaces, so it is only useful on a terminal.
    setProgressEnabled(!QUIET && (VERBOSE || isatty(STDERR_FILENO)));
    if (QUIET) cout.rdbuf(&GARLIC_NULLBUF);

    string outdir = params->getStringFlag(ARG_OUTDIR);
    if (!outdir.empty())
    {
        if (makeOutdir(outdir)) return OPTIONS_USAGE_ERROR;
        if (outdir[outdir.size() - 1] != '/') outdir += "/";
        opt.outfile = outdir + opt.outfile;
    }

    LOG.init(opt.outfile);
    LOG.log(getCommandLineString(argc, argv));
    LOG.log("Output file basename:", opt.outfile);

    //Accumulator for validation failures.  It is monotone -- every check does
    //`argerr = argerr || check(...)` -- so a check whose result is not tested
    //until a few flags later defers the report rather than losing it, and
    //nothing reaches the pipeline with a failed check outstanding.  Three
    //checks (--gl-type, --map with --weighted/--cm, --mode-smooth-span) are in
    //that position; each was verified to still exit 1 on its own.  The tests
    //are deliberately left where they are: moving them earlier would change
    //which error a user sees first when more than one flag is wrong.
    bool argerr = false;

    opt.tpedfile = params->getStringFlag(ARG_TPED);
    opt.tfamfile = params->getStringFlag(ARG_TFAM);
    opt.tglsfile = params->getStringFlag(ARG_TGLS);
    argerr = argerr || checkRequiredFiles(opt.tpedfile, opt.tfamfile);
    if (argerr) return OPTIONS_USAGE_ERROR;
    LOG.log("TPED file:", opt.tpedfile);

    opt.TPED_MISSING = params->getCharFlag(ARG_TPED_MISSING);
    LOG.log("TPED missing data code:", opt.TPED_MISSING);
    LOG.log("TFAM file:", opt.tfamfile);
    LOG.log("TGLS file:", opt.tglsfile);

    opt.GL_TYPE = params->getStringFlag(ARG_GL_TYPE);
    argerr = argerr || checkGLType(opt.GL_TYPE, opt.tglsfile);
    LOG.log("Genotype likelihood format:", opt.GL_TYPE);

    opt.WEIGHTED = params->getBoolFlag(ARG_WEIGHTED);
    opt.mapfile = params->getStringFlag(ARG_MAP);
    opt.CM = params->getBoolFlag(ARG_CM);
    argerr = argerr || checkCM(opt.mapfile, opt.CM);
    if (argerr) return OPTIONS_USAGE_ERROR;
    LOG.log("Measure ROH in genetic distance units:", opt.CM);
    argerr = argerr || checkMapFile(opt.mapfile, opt.WEIGHTED || opt.CM);
    LOG.log("Weighted LOD:", opt.WEIGHTED);
    if (opt.WEIGHTED) {
        LOG.log("Map file:", opt.mapfile);
    }

    opt.BUILD = params->getStringFlag(ARG_BUILD);
    argerr = argerr || checkBuild(opt.BUILD);
    if (argerr) return OPTIONS_USAGE_ERROR;
    LOG.log("Genome build:", opt.BUILD);

    opt.centromereFile = params->getStringFlag(ARG_CENTROMERE_FILE);
    opt.NO_CENTROMERE = params->getBoolFlag(ARG_NO_CENTROMERE);
    argerr = argerr || checkBuildAndCentromereFile(opt.BUILD, opt.centromereFile, opt.NO_CENTROMERE);
    if (argerr) return OPTIONS_USAGE_ERROR;
    LOG.log("User defined centromere file:", opt.centromereFile);

    opt.nresample = params->getIntFlag(ARG_RESAMPLE);
    opt.freqfile = params->getStringFlag(ARG_FREQ_FILE);
    bool FREQ_ONLY = params->getBoolFlag(ARG_FREQ_ONLY);
    opt.AUTO_FREQ = true;
    argerr = argerr || checkAutoFreq(opt.freqfile, FREQ_ONLY, opt.AUTO_FREQ);
    if (argerr) return OPTIONS_USAGE_ERROR;
    LOG.log("Calculate allele frequencies only:", FREQ_ONLY);
    LOG.log("Calculate allele frequencies from data:", opt.AUTO_FREQ);
    if (!opt.AUTO_FREQ) LOG.log("Allele frequencies file:", opt.freqfile);
    else
    {
        if (opt.nresample <= 0) LOG.log("Allele frequencies resampled: FALSE");
        else LOG.log("Allele frequencies resampled:", opt.nresample);
    }

    opt.multiWinsizes = params->getIntListFlag(ARG_WINSIZE_MULTI);
    opt.WINSIZE_EXPLORE = false;
    argerr = argerr || checkMultiWinsizes(opt.multiWinsizes, opt.WINSIZE_EXPLORE, params->isFlagSet(ARG_WINSIZE_MULTI));
    if (argerr) return OPTIONS_USAGE_ERROR;
    LOG.log("Explore window sizes:", opt.WINSIZE_EXPLORE);
    if (opt.WINSIZE_EXPLORE) LOG.logv("User defined window sizes:", opt.multiWinsizes);

    opt.AUTO_WINSIZE = params->getBoolFlag(ARG_AUTO_WINSIZE);
    argerr = argerr || checkAutoWinsize(opt.WINSIZE_EXPLORE, opt.AUTO_WINSIZE);
    if (argerr) return OPTIONS_USAGE_ERROR;
    LOG.log("Automatic window size:", opt.AUTO_WINSIZE);

    opt.AUTO_WINSIZE_STEP = params->getIntFlag(ARG_AUTO_WINSIZE_STEP);
    argerr = argerr || checkAutoWinsizeStep(opt.AUTO_WINSIZE_STEP);
    if (argerr) return OPTIONS_USAGE_ERROR;
    LOG.log("Automatic window step size:", opt.AUTO_WINSIZE_STEP);

    opt.winsize = params->getIntFlag(ARG_WINSIZE);
    argerr = argerr || checkWinsize(opt.winsize, opt.WINSIZE_EXPLORE, opt.AUTO_WINSIZE, opt.WEIGHTED, FREQ_ONLY);
    if (argerr) return OPTIONS_USAGE_ERROR;
    if (!opt.WINSIZE_EXPLORE && !opt.AUTO_WINSIZE) LOG.log("User defined window size:", opt.winsize);

    opt.LOD_CUTOFF = params->getDoubleFlag(ARG_LOD_CUTOFF);
    opt.AUTO_CUTOFF = true;
    argerr = argerr || checkAutoCutoff(opt.LOD_CUTOFF, opt.AUTO_CUTOFF, params->isFlagSet(ARG_LOD_CUTOFF));
    if (argerr) return OPTIONS_USAGE_ERROR;
    LOG.log("Choose LOD score cutoff automatically:", opt.AUTO_CUTOFF);
    if (!opt.AUTO_CUTOFF) LOG.log("User defined LOD score cutoff:", opt.LOD_CUTOFF);

    opt.boundSizes = params->getDoubleListFlag(ARG_BOUND_SIZE);
    opt.AUTO_BOUNDS = true;
    argerr = argerr || checkBoundSizes(opt.boundSizes, opt.AUTO_BOUNDS, params->isFlagSet(ARG_BOUND_SIZE));
    if (argerr) return OPTIONS_USAGE_ERROR;
    LOG.log("Choose ROH class thresholds automatically:", opt.AUTO_BOUNDS);
    if (!opt.AUTO_BOUNDS) LOG.logv("User defined ROH class thresholds:", opt.boundSizes);

    opt.numThreads = params->getIntFlag(ARG_THREADS);
    argerr = argerr || checkThreads(opt.numThreads);
    if (argerr) return OPTIONS_USAGE_ERROR;
    LOG.log("Threads:", opt.numThreads);
    //--threads used to affect only the weighted LD stage; the KDE targets are
    //independent, so give it the same budget.
    setKDEThreads(opt.numThreads);
    setLODThreads(opt.numThreads);

    opt.error = params->getDoubleFlag(ARG_ERROR);
    argerr = argerr || checkError(opt.error, opt.tglsfile, params->isFlagSet(ARG_ERROR));
    if (argerr) return OPTIONS_USAGE_ERROR;
    LOG.log("Genotyping error:", opt.error);

    opt.MAX_GAP = params->getIntFlag(ARG_MAX_GAP);
    argerr = argerr || checkMaxGap(opt.MAX_GAP);
    if (argerr) return OPTIONS_USAGE_ERROR;
    LOG.log("Max gap:", opt.MAX_GAP);

    opt.OVERLAP_FRAC = params->getDoubleFlag(ARG_OVERLAP_FRAC);
    argerr = argerr || checkOverlapFrac(opt.OVERLAP_FRAC);
    if (argerr) return OPTIONS_USAGE_ERROR;
    opt.AUTO_OVERLAP_FRAC = params->getBoolFlag(ARG_AUTO_OVERLAP_FRAC);
    if(opt.AUTO_OVERLAP_FRAC) LOG.log("Overlap fraction: automatic");
    else if(opt.OVERLAP_FRAC != 0) LOG.log("Overlap fraction:", opt.OVERLAP_FRAC);
    else LOG.log("Overlap fraction: 1/winsize");

    opt.mu = params->getDoubleFlag(ARG_MU);
    argerr = argerr || checkMU(opt.mu);
    if (argerr) return OPTIONS_USAGE_ERROR;
    LOG.log("mu:", opt.mu);

    opt.M = params->getIntFlag(ARG_M);
    argerr = argerr || checkM(opt.M);
    if (argerr) return OPTIONS_USAGE_ERROR;
    LOG.log("M:", opt.M);

    opt.NCLUST = params->getIntFlag(ARG_NCLUST);
    argerr = argerr || checkNCLUST(opt.NCLUST);
    if (argerr) return OPTIONS_USAGE_ERROR;
    LOG.log("# GMM clusters:", opt.NCLUST);
    warnBoundsOverridesNclust(params->isFlagSet(ARG_BOUND_SIZE), params->isFlagSet(ARG_NCLUST));

    opt.KDE_SUBSAMPLE = params->getIntFlag(ARG_KDE_SUBSAMPLE);
    if (opt.KDE_SUBSAMPLE <= 0) LOG.log("# of rand individuals for KDE: ALL");
    else LOG.log("# of rand individuals for KDE:", opt.KDE_SUBSAMPLE);

    opt.LD_SUBSAMPLE = params->getIntFlag(ARG_LD_SUBSAMPLE);
    if (opt.LD_SUBSAMPLE <= 0) LOG.log("# of rand individuals for LD: ALL");
    else LOG.log("# of rand individuals for LD:", opt.LD_SUBSAMPLE);

    opt.RAW_LOD = params->getBoolFlag(ARG_RAW_LOD);
    LOG.log("Output raw LOD scores:", opt.RAW_LOD);

    opt.PHASED = params->getBoolFlag(ARG_PHASED);
    LOG.log("Use r2 for weighting phased data:", opt.PHASED);

    //Values that determine a scientific result and used to be unreachable
    //constants in the source.
    double AUTO_WINSIZE_THRESHOLD = params->getDoubleFlag(ARG_AUTO_WINSIZE_THRESHOLD);
    argerr = argerr || checkAutoWinsizeThreshold(AUTO_WINSIZE_THRESHOLD);
    if (argerr) return OPTIONS_USAGE_ERROR;
    LOG.log("Auto window size smoothness threshold:", AUTO_WINSIZE_THRESHOLD);
    setAutoWinsizeThreshold(AUTO_WINSIZE_THRESHOLD);

    int KDE_POINTS = params->getIntFlag(ARG_KDE_POINTS);
    double KDE_CUT = params->getDoubleFlag(ARG_KDE_CUT);
    argerr = argerr || checkKDEPoints(KDE_POINTS) || checkKDECut(KDE_CUT);
    if (argerr) return OPTIONS_USAGE_ERROR;
    LOG.log("KDE grid points:", KDE_POINTS);
    LOG.log("KDE range extension (bandwidths):", KDE_CUT);
    setKDEGrid(KDE_POINTS, KDE_CUT);

    int MODE_SPAN = params->getIntFlag(ARG_MODE_SPAN);
    argerr = argerr || checkModeSpan(MODE_SPAN);
    if (MODE_SPAN >= KDE_POINTS) {
        LOG.err("ERROR: --mode-smooth-span must be smaller than --kde-points.");
        argerr = true;
    }
    if (argerr) return OPTIONS_USAGE_ERROR;
    LOG.log("Mode smoothing span:", MODE_SPAN);
    setModeSpan(MODE_SPAN);

    vector<double> awCoef = params->getDoubleListFlag(ARG_AUTO_WINSIZE_COEF);
    if (params->isFlagSet(ARG_AUTO_WINSIZE_COEF)) {
        argerr = argerr || checkCoefPair(awCoef, ARG_AUTO_WINSIZE_COEF);
        if (argerr) return OPTIONS_USAGE_ERROR;
        LOG.logv("Auto window size coefficients (slope intercept):", awCoef);
        setAutoWinsizeCoef(awCoef[0], awCoef[1]);
    }

    vector<double> aoCoef = params->getDoubleListFlag(ARG_AUTO_OVERLAP_COEF);
    if (params->isFlagSet(ARG_AUTO_OVERLAP_COEF)) {
        argerr = argerr || checkCoefPair(aoCoef, ARG_AUTO_OVERLAP_COEF);
        if (argerr) return OPTIONS_USAGE_ERROR;
        LOG.logv("Auto overlap fraction coefficients (slope intercept):", aoCoef);
        setAutoOverlapCoef(aoCoef[0], aoCoef[1]);
    }

    int GMM_MAX_ITER = params->getIntFlag(ARG_GMM_MAX_ITER);
    double GMM_TOL = params->getDoubleFlag(ARG_GMM_TOL);
    argerr = argerr || checkGMMParams(GMM_MAX_ITER, GMM_TOL);
    if (argerr) return OPTIONS_USAGE_ERROR;
    LOG.log("GMM max iterations:", GMM_MAX_ITER);
    LOG.log("GMM tolerance:", GMM_TOL);
    setGMMParams(GMM_MAX_ITER, GMM_TOL);

    //0 means "follow the window size"; --no-kde-thinning is the old spelling of 1.
    opt.KDE_THIN_STEP = params->getIntFlag(ARG_KDE_THIN_STEP);
    argerr = argerr || checkKDEThinStep(opt.KDE_THIN_STEP);
    if (argerr) return OPTIONS_USAGE_ERROR;
    if (!params->isFlagSet(ARG_KDE_THIN_STEP) && params->getBoolFlag(ARG_KDE_THINNING)) opt.KDE_THIN_STEP = 1;
    LOG.log("KDE thinning step (0 = window size):", opt.KDE_THIN_STEP);
    //double AUTO_WINSIZE_THRESHOLD = 0.5;

    opt.MAX_WINSIZE = params->getIntFlag(ARG_MAX_WINSIZE);
    argerr = argerr || checkMaxWinsize(opt.MAX_WINSIZE, opt.winsize);
    if (argerr) return OPTIONS_USAGE_ERROR;

    int seedFlag = params->getIntFlag(ARG_SEED);
    argerr = argerr || checkSeed(seedFlag);
    if (argerr) return OPTIONS_USAGE_ERROR;
    opt.SEED = (seedFlag == 0) ? drawRandomSeed() : (unsigned long int)(seedFlag);
    {   //logged as a string: errlog has no unsigned long overload
        stringstream seedss;
        seedss << opt.SEED;
        LOG.log("Random seed:", seedss.str());
        if (seedFlag == 0) LOG.log("\t(drawn automatically; pass --seed with this value to reproduce this run)");
    }
    //So the flags block of <out>.params.json is a command line that reproduces
    //this run, rather than one that draws a fresh seed.
    params->setIntFlag(ARG_SEED, int(opt.SEED));
    initRNG(opt.SEED);

    //All argument validation has passed.  Refuse to clobber a previous run's
    //calls, then materialise <out>.log -- everything logged above has been held
    //in memory so that a rejected command line leaves no files behind.
    bool FORCE = params->getBoolFlag(ARG_FORCE);
    if (checkOutfileClobber(opt.outfile, FORCE)) return OPTIONS_USAGE_ERROR;
    LOG.commit();

    if (FREQ_ONLY){//calculated on the fly, as the file is read, to save RAM
        freqOnly(opt.tpedfile,opt.outfile,opt.nresample,opt.TPED_MISSING);
        freeRNG();
        return OPTIONS_DONE;
    }

    return OPTIONS_OK;
}
