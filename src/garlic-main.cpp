#include "garlic-errlog.h"
#include "garlic-cli.h"
#include <iostream>
#include <cstdio>
#include <unistd.h>
#include <fstream>
#include <sstream>
#include <cmath>
#include <pthread.h>
#include "garlic-data.h"
#include "garlic-roh.h"
#include "garlic-kde.h"
#include "param_t.h"
#include "garlic-centromeres.h"
#include <gsl/gsl_errno.h>

using namespace std;

string getCommandLineString(int argc, char *argv[])
{
    string str = argv[0];
    for (int i = 1; i < argc; i++) {
        str += " " + string(argv[i]);
    }
    return str;
}

//GSL's default handler calls abort(), so a domain error deep in the GMM killed
//the process with SIGABRT after all the work was done.  Turn it into the throw
//the surrounding code already handles.
static void garlicGSLError(const char *reason, const char *file, int line, int gsl_errno)
{
    LOG.err("ERROR: numerical failure in GSL:", string(reason));
    LOG.err("\tat", string(file), false);
    LOG.err(":", line);
    LOG.err("\tThis usually means degenerate input to the size-class GMM (for example");
    LOG.err("\ta LOD cutoff so low that every window is called, giving near-identical");
    LOG.err("\tROH lengths). Pass --size-bounds to set the boundaries explicitly.");
    (void)gsl_errno;
    throw 0;
}

//Discards whatever is written to it; used to silence stdout under --quiet
//without touching every cout site.
class NullBuf : public std::streambuf
{
protected:
    int overflow(int c) { return c; }
};
static NullBuf GARLIC_NULLBUF;

int main(int argc, char *argv[])
{
    gsl_set_error_handler(&garlicGSLError);

    #ifdef PTW32_STATIC_LIB
        pthread_win32_process_attach_np();
    #endif
//++++++++++CLI handling++++++++++
    int cliStatus = PARAM_OK;
    param_t *params = getCLI(argc, argv, cliStatus);
    //0 success, 1 usage error, 2 runtime error.
    if (params == NULL) return (cliStatus == PARAM_HELP) ? 0 : 1;

    string outfile = params->getStringFlag(ARG_OUTFILE);
    bool QUIET = params->getBoolFlag(ARG_QUIET);
    bool VERBOSE = params->getBoolFlag(ARG_VERBOSE);
    if (QUIET && VERBOSE) {
        LOG.err("ERROR: --quiet and --verbose are mutually exclusive.");
        return 1;
    }
    LOG.setVerbosity(QUIET, VERBOSE);
    //The bar writes backspaces, so it is only useful on a terminal.
    setProgressEnabled(!QUIET && (VERBOSE || isatty(STDERR_FILENO)));
    if (QUIET) cout.rdbuf(&GARLIC_NULLBUF);

    string outdir = params->getStringFlag(ARG_OUTDIR);
    if (!outdir.empty())
    {
        if (makeOutdir(outdir)) return 1;
        if (outdir[outdir.size() - 1] != '/') outdir += "/";
        outfile = outdir + outfile;
    }

    LOG.init(outfile);
    LOG.log(getCommandLineString(argc, argv));
    LOG.log("Output file basename:", outfile);

    bool argerr = false;

    string tpedfile = params->getStringFlag(ARG_TPED);
    string tfamfile = params->getStringFlag(ARG_TFAM);
    string tglsfile = params->getStringFlag(ARG_TGLS);
    argerr = argerr || checkRequiredFiles(tpedfile, tfamfile);
    if (argerr) return 1;
    LOG.log("TPED file:", tpedfile);

    char TPED_MISSING = params->getCharFlag(ARG_TPED_MISSING);
    LOG.log("TPED missing data code:", TPED_MISSING);
    LOG.log("TFAM file:", tfamfile);
    LOG.log("TGLS file:", tglsfile);

    string GL_TYPE = params->getStringFlag(ARG_GL_TYPE);
    argerr = argerr || checkGLType(GL_TYPE, tglsfile);
    LOG.log("Genotype likelihood format:", GL_TYPE);

    bool WEIGHTED = params->getBoolFlag(ARG_WEIGHTED);
    string mapfile = params->getStringFlag(ARG_MAP);
    bool CM = params->getBoolFlag(ARG_CM);
    argerr = argerr || checkCM(mapfile, CM);
    if (argerr) return 1;
    LOG.log("Measure ROH in genetic distance units:", CM);
    argerr = argerr || checkMapFile(mapfile, WEIGHTED || CM);
    LOG.log("Weighted LOD:", WEIGHTED);
    if (WEIGHTED) {
        LOG.log("Map file:", mapfile);
    }

    string BUILD = params->getStringFlag(ARG_BUILD);
    argerr = argerr || checkBuild(BUILD);
    if (argerr) return 1;
    LOG.log("Genome build:", BUILD);

    string centromereFile = params->getStringFlag(ARG_CENTROMERE_FILE);
    argerr = argerr || checkBuildAndCentromereFile(BUILD, centromereFile);
    if (argerr) return 1;
    LOG.log("User defined centromere file:", centromereFile);

    int nresample = params->getIntFlag(ARG_RESAMPLE);
    string freqfile = params->getStringFlag(ARG_FREQ_FILE);
    bool FREQ_ONLY = params->getBoolFlag(ARG_FREQ_ONLY);
    bool AUTO_FREQ = true;
    argerr = argerr || checkAutoFreq(freqfile, FREQ_ONLY, AUTO_FREQ);
    if (argerr) return 1;
    LOG.log("Calculate allele frequencies only:", FREQ_ONLY);
    LOG.log("Calculate allele frequencies from data:", AUTO_FREQ);
    if (!AUTO_FREQ) LOG.log("Allele frequencies file:", freqfile);
    else
    {
        if (nresample <= 0) LOG.log("Allele frequencies resampled: FALSE");
        else LOG.log("Allele frequencies resampled:", nresample);
    }

    vector<int> multiWinsizes = params->getIntListFlag(ARG_WINSIZE_MULTI);
    bool WINSIZE_EXPLORE = false;
    argerr = argerr || checkMultiWinsizes(multiWinsizes, WINSIZE_EXPLORE, params->isFlagSet(ARG_WINSIZE_MULTI));
    if (argerr) return 1;
    LOG.log("Explore window sizes:", WINSIZE_EXPLORE);
    if (WINSIZE_EXPLORE) LOG.logv("User defined window sizes:", multiWinsizes);

    bool AUTO_WINSIZE = params->getBoolFlag(ARG_AUTO_WINSIZE);
    argerr = argerr || checkAutoWinsize(WINSIZE_EXPLORE, AUTO_WINSIZE);
    if (argerr) return 1;
    LOG.log("Automatic window size:", AUTO_WINSIZE);

    int AUTO_WINSIZE_STEP = params->getIntFlag(ARG_AUTO_WINSIZE_STEP);
    argerr = argerr || checkAutoWinsizeStep(AUTO_WINSIZE_STEP);
    if (argerr) return 1;
    LOG.log("Automatic window step size:", AUTO_WINSIZE_STEP);

    int winsize = params->getIntFlag(ARG_WINSIZE);
    argerr = argerr || checkWinsize(winsize, WINSIZE_EXPLORE, AUTO_WINSIZE, WEIGHTED, FREQ_ONLY);
    if (argerr) return 1;
    if (!WINSIZE_EXPLORE && !AUTO_WINSIZE) LOG.log("User defined window size:", winsize);

    double LOD_CUTOFF = params->getDoubleFlag(ARG_LOD_CUTOFF);
    bool AUTO_CUTOFF = true;
    argerr = argerr || checkAutoCutoff(LOD_CUTOFF, AUTO_CUTOFF, params->isFlagSet(ARG_LOD_CUTOFF));
    if (argerr) return 1;
    LOG.log("Choose LOD score cutoff automatically:", AUTO_CUTOFF);
    if (!AUTO_CUTOFF) LOG.log("User defined LOD score cutoff:", LOD_CUTOFF);

    vector<double> boundSizes = params->getDoubleListFlag(ARG_BOUND_SIZE);
    bool AUTO_BOUNDS = true;
    argerr = argerr || checkBoundSizes(boundSizes, AUTO_BOUNDS, params->isFlagSet(ARG_BOUND_SIZE));
    if (argerr) return 1;
    LOG.log("Choose ROH class thresholds automatically:", AUTO_BOUNDS);
    if (!AUTO_BOUNDS) LOG.logv("User defined ROH class thresholds:", boundSizes);

    int numThreads = params->getIntFlag(ARG_THREADS);
    argerr = argerr || checkThreads(numThreads);
    if (argerr) return 1;
    LOG.log("Threads:", numThreads);
    //--threads used to affect only the weighted LD stage; the KDE targets are
    //independent, so give it the same budget.
    setKDEThreads(numThreads);
    setLODThreads(numThreads);

    double error = params->getDoubleFlag(ARG_ERROR);
    argerr = argerr || checkError(error, tglsfile, params->isFlagSet(ARG_ERROR));
    if (argerr) return 1;
    LOG.log("Genotyping error:", error);

    int MAX_GAP = params->getIntFlag(ARG_MAX_GAP);
    argerr = argerr || checkMaxGap(MAX_GAP);
    if (argerr) return 1;
    LOG.log("Max gap:", MAX_GAP);

    double OVERLAP_FRAC = params->getDoubleFlag(ARG_OVERLAP_FRAC);
    argerr = argerr || checkOverlapFrac(OVERLAP_FRAC);
    if (argerr) return 1;
    bool AUTO_OVERLAP_FRAC = params->getBoolFlag(ARG_AUTO_OVERLAP_FRAC);
    if(AUTO_OVERLAP_FRAC) LOG.log("Overlap fraction: automatic");
    else if(OVERLAP_FRAC != 0) LOG.log("Overlap fraction:", OVERLAP_FRAC);
    else LOG.log("Overlap fraction: 1/winsize");

    double mu = params->getDoubleFlag(ARG_MU);
    argerr = argerr || checkMU(mu);
    if (argerr) return 1;
    LOG.log("mu:", mu);

    int M = params->getIntFlag(ARG_M);
    argerr = argerr || checkM(M);
    if (argerr) return 1;
    LOG.log("M:", M);

    int NCLUST = params->getIntFlag(ARG_NCLUST);
    argerr = argerr || checkNCLUST(NCLUST);
    if (argerr) return 1;
    LOG.log("# GMM clusters:", NCLUST);
    warnBoundsOverridesNclust(params->isFlagSet(ARG_BOUND_SIZE), params->isFlagSet(ARG_NCLUST));

    int KDE_SUBSAMPLE = params->getIntFlag(ARG_KDE_SUBSAMPLE);
    if (KDE_SUBSAMPLE <= 0) LOG.log("# of rand individuals for KDE: ALL");
    else LOG.log("# of rand individuals for KDE:", KDE_SUBSAMPLE);

    int LD_SUBSAMPLE = params->getIntFlag(ARG_LD_SUBSAMPLE);
    if (LD_SUBSAMPLE <= 0) LOG.log("# of rand individuals for LD: ALL");
    else LOG.log("# of rand individuals for LD:", LD_SUBSAMPLE);

    bool RAW_LOD = params->getBoolFlag(ARG_RAW_LOD);
    LOG.log("Output raw LOD scores:", RAW_LOD);

    bool PHASED = params->getBoolFlag(ARG_PHASED);
    LOG.log("Use r2 for weighting phased data:", PHASED);

    //Values that determine a scientific result and used to be unreachable
    //constants in the source.
    double AUTO_WINSIZE_THRESHOLD = params->getDoubleFlag(ARG_AUTO_WINSIZE_THRESHOLD);
    argerr = argerr || checkAutoWinsizeThreshold(AUTO_WINSIZE_THRESHOLD);
    if (argerr) return 1;
    LOG.log("Auto window size smoothness threshold:", AUTO_WINSIZE_THRESHOLD);
    setAutoWinsizeThreshold(AUTO_WINSIZE_THRESHOLD);

    int KDE_POINTS = params->getIntFlag(ARG_KDE_POINTS);
    double KDE_CUT = params->getDoubleFlag(ARG_KDE_CUT);
    argerr = argerr || checkKDEPoints(KDE_POINTS) || checkKDECut(KDE_CUT);
    if (argerr) return 1;
    LOG.log("KDE grid points:", KDE_POINTS);
    LOG.log("KDE range extension (bandwidths):", KDE_CUT);
    setKDEGrid(KDE_POINTS, KDE_CUT);

    int MODE_SPAN = params->getIntFlag(ARG_MODE_SPAN);
    argerr = argerr || checkModeSpan(MODE_SPAN);
    if (MODE_SPAN >= KDE_POINTS) {
        LOG.err("ERROR: --mode-smooth-span must be smaller than --kde-points.");
        argerr = true;
    }
    if (argerr) return 1;
    LOG.log("Mode smoothing span:", MODE_SPAN);
    setModeSpan(MODE_SPAN);

    vector<double> awCoef = params->getDoubleListFlag(ARG_AUTO_WINSIZE_COEF);
    if (params->isFlagSet(ARG_AUTO_WINSIZE_COEF)) {
        argerr = argerr || checkCoefPair(awCoef, ARG_AUTO_WINSIZE_COEF);
        if (argerr) return 1;
        LOG.logv("Auto window size coefficients (slope intercept):", awCoef);
        setAutoWinsizeCoef(awCoef[0], awCoef[1]);
    }

    vector<double> aoCoef = params->getDoubleListFlag(ARG_AUTO_OVERLAP_COEF);
    if (params->isFlagSet(ARG_AUTO_OVERLAP_COEF)) {
        argerr = argerr || checkCoefPair(aoCoef, ARG_AUTO_OVERLAP_COEF);
        if (argerr) return 1;
        LOG.logv("Auto overlap fraction coefficients (slope intercept):", aoCoef);
        setAutoOverlapCoef(aoCoef[0], aoCoef[1]);
    }

    int GMM_MAX_ITER = params->getIntFlag(ARG_GMM_MAX_ITER);
    double GMM_TOL = params->getDoubleFlag(ARG_GMM_TOL);
    argerr = argerr || checkGMMParams(GMM_MAX_ITER, GMM_TOL);
    if (argerr) return 1;
    LOG.log("GMM max iterations:", GMM_MAX_ITER);
    LOG.log("GMM tolerance:", GMM_TOL);
    setGMMParams(GMM_MAX_ITER, GMM_TOL);

    //0 means "follow the window size"; --no-kde-thinning is the old spelling of 1.
    int KDE_THIN_STEP = params->getIntFlag(ARG_KDE_THIN_STEP);
    argerr = argerr || checkKDEThinStep(KDE_THIN_STEP);
    if (argerr) return 1;
    if (!params->isFlagSet(ARG_KDE_THIN_STEP) && params->getBoolFlag(ARG_KDE_THINNING)) KDE_THIN_STEP = 1;
    LOG.log("KDE thinning step (0 = window size):", KDE_THIN_STEP);
    //double AUTO_WINSIZE_THRESHOLD = 0.5;

    int MAX_WINSIZE = params->getIntFlag(ARG_MAX_WINSIZE);
    argerr = argerr || checkMaxWinsize(MAX_WINSIZE, winsize);
    if (argerr) return 1;

    int seedFlag = params->getIntFlag(ARG_SEED);
    argerr = argerr || checkSeed(seedFlag);
    if (argerr) return 1;
    unsigned long int SEED = (seedFlag == 0) ? drawRandomSeed() : (unsigned long int)(seedFlag);
    {   //logged as a string: errlog has no unsigned long overload
        stringstream seedss;
        seedss << SEED;
        LOG.log("Random seed:", seedss.str());
        if (seedFlag == 0) LOG.log("\t(drawn automatically; pass --seed with this value to reproduce this run)");
    }
    //So the flags block of <out>.params.json is a command line that reproduces
    //this run, rather than one that draws a fresh seed.
    params->setIntFlag(ARG_SEED, int(SEED));
    initRNG(SEED);

    //All argument validation has passed.  Refuse to clobber a previous run's
    //calls, then materialise <out>.log -- everything logged above has been held
    //in memory so that a rejected command line leaves no files behind.
    bool FORCE = params->getBoolFlag(ARG_FORCE);
    if (checkOutfileClobber(outfile, FORCE)) return 1;
    LOG.commit();

    if (FREQ_ONLY){//calculated on the fly, as the file is read, to save RAM
        freqOnly(tpedfile,outfile,nresample,TPED_MISSING);
        freeRNG();
        return 0;
    }

//++++++++++Datafile reading++++++++++
    centromere *centro;
    centro = new centromere(BUILD, centromereFile, DEFAULT_CENTROMERE_FILE);

    int numLoci, numInd;
    double variantDensity = -1;;
    //vector< int_pair_t > *chrCoordList = NULL;
    vector< MapData * > *mapDataByChr = NULL;
    //string popName;
    IndData *indData = NULL;
    vector< HapData * > *hapDataByChr = NULL;
    vector< FreqData * > *freqDataByChr = NULL;
    vector< GenoFreqData * > *genoFreqDataByChr = NULL;
    vector< WinData * > *winDataByChr = NULL;
    vector< GenoLikeData * > *GLDataByChr = NULL;
    vector< GenMapScaffold *> *scaffoldMapByChr = NULL;
    vector< LDData * > *ldDataByChr = NULL;
    KDEResult *kdeResult = NULL;
    bool USE_GL = false;
    try
    {
        hapDataByChr = new vector< HapData * >;
        mapDataByChr = new vector< MapData * >;
        if(AUTO_FREQ) freqDataByChr = new vector< FreqData * >;

        loadTPEDData(tpedfile, numLoci, numInd,
                     &hapDataByChr, &mapDataByChr, &freqDataByChr,
                     TPED_MISSING, nresample, PHASED, AUTO_FREQ);

        LOG.log("Total loci:", numLoci);

        scanIndData3(tfamfile, numInd);
        indData = readIndData3(tfamfile, numInd);

        //LOG.log("Population:", popName);
        LOG.log("Total diploid individuals:", numInd);

        if (tglsfile.compare(DEFAULT_TGLS) != 0) {
            GLDataByChr = readTGLSData(tglsfile, numLoci, numInd, mapDataByChr, GL_TYPE);
            USE_GL = true;
        }

        if (WEIGHTED || CM) {
            scaffoldMapByChr = loadMapScaffold(mapfile, centro);
            if (scaffoldMapByChr->size() != mapDataByChr->size()) {
                LOG.err("ERROR: Scaffold genetic map does not have the same number of chromosomes as data.");
                return 2;
            }
            //Match scaffolds to data by chromosome NAME.  Everything downstream
            //zips the two vectors positionally, so a map file sorted
            //chr1, chr10, chr11, ... against a TPED sorted chr1, chr2, ...
            //used to apply the wrong chromosome's map to every site, silently.
            if (!alignMapScaffold(scaffoldMapByChr, mapDataByChr)) return 2;
        }

    }
    catch (...) { return 2; }

//++++++++++Allele frequencies++++++++++
    if (AUTO_FREQ)
    {
        //cout << "Calculating allele frequencies\n";
        //freqDataByChr = calcFreqData2(hapDataByChr, nresample);

        string freqOutfile = outfile;
        freqOutfile += ".freq";
        writeFreqData(freqOutfile, freqDataByChr, mapDataByChr, indData);
    }
    else //(!AUTO_FREQ)
    {
        cout << "Loading user provided allele frequencies from " << freqfile << "\n";
        try { freqDataByChr = readFreqData(freqfile, mapDataByChr); }
        catch (...) { return 2; }
    }

//Filter data based on frequency data.
//Remove all monomorphic sites.
//If a frequency file is provided that reports
//a frequency in (0,1) the site will be retained
//even if it appears monomorphic in the sample.

    int newLoci;

    if (WEIGHTED || CM) {
        newLoci = filterMonomorphicAndOOBSites(&mapDataByChr, &hapDataByChr, &freqDataByChr, &GLDataByChr, scaffoldMapByChr, USE_GL, PHASED);
        LOG.log("Monomorphic or out of bounds loci filtered:", numLoci - newLoci);
        int numInterpolated = interpolateGeneticmap(&mapDataByChr, scaffoldMapByChr);

        LOG.log("Number of genetic map locations interpolated:", numInterpolated);
        releaseGenMapScaffold(scaffoldMapByChr);
        if(!PHASED && WEIGHTED) genoFreqDataByChr = calculateGenoFreq(hapDataByChr);
    }
    else {
        newLoci = filterMonomorphicSites(&mapDataByChr, &hapDataByChr, &freqDataByChr, &GLDataByChr, USE_GL, PHASED);
        LOG.log("Monomorphic loci filtered:", numLoci - newLoci);
    }

    LOG.log("Total loci used for analysis:", newLoci);

    numLoci = newLoci;

    if((AUTO_WINSIZE && WEIGHTED) || AUTO_OVERLAP_FRAC){
        variantDensity = calcDensity(numLoci, mapDataByChr, centro);
    }

    //chrCoordList->clear();
    //delete chrCoordList;

    //A hemizygous male genotype is indistinguishable from a homozygous call,
    //so warn whenever a sex chromosome is present and offer to drop it.
    bool haveSexChr = warnSexChromosomes(mapDataByChr, indData);

    if (params->getBoolFlag(ARG_AUTOSOMES_ONLY))
    {
        if (!haveSexChr)
        {
            LOG.log("--autosomes-only: no sex chromosomes in the data, nothing to drop.");
        }
        else
        {
            vector<string> autosomes;
            for (unsigned int chr = 0; chr < mapDataByChr->size(); chr++)
                if (!isSexChromosome(mapDataByChr->at(chr)->chr))
                    autosomes.push_back(mapDataByChr->at(chr)->chr);
            if (autosomes.empty())
            {
                LOG.err("ERROR: --autosomes-only leaves no data: every chromosome is a sex chromosome.");
                return 1;
            }
            int nkept = filterChromosomes(autosomes, &mapDataByChr, &hapDataByChr,
                                          &freqDataByChr, &GLDataByChr, USE_GL);
            LOG.log("--autosomes-only: kept", nkept);
            LOG.log("chromosomes after dropping sex chromosomes:", int(mapDataByChr->size()));
        }
    }

    vector<string> keepChr = params->getStringListFlag(ARG_CHR);
    if (params->isFlagSet(ARG_CHR))
    {
        int nkept = filterChromosomes(keepChr, &mapDataByChr, &hapDataByChr,
                                      &freqDataByChr, &GLDataByChr, USE_GL);
        if (nkept < 0) return 1;
        LOG.loga("Restricted to chromosomes:", &(keepChr[0]), int(keepChr.size()));
        LOG.log("Chromosomes analysed:", nkept);
        if (!PHASED && WEIGHTED)
        {
            releaseGenoFreq(genoFreqDataByChr);
            genoFreqDataByChr = calculateGenoFreq(hapDataByChr);
        }
    }

//++++++++++Pipeline begins++++++++++
    if (WINSIZE_EXPLORE && AUTO_WINSIZE && !WEIGHTED)
    {
        kdeResult = selectWinsizeFromList(hapDataByChr, freqDataByChr, mapDataByChr,
                                          indData, centro, &multiWinsizes, winsize, error,
                                          GLDataByChr, USE_GL,
                                          MAX_GAP, KDE_SUBSAMPLE, outfile, WEIGHTED, genoFreqDataByChr, PHASED, KDE_THIN_STEP);
    }
    else if (WINSIZE_EXPLORE)
    {
        /*
        KDEWinsizeReport *winsizeReport =  calculateLODOverWinsizeRange(hapDataByChr, freqDataByChr,
                                           mapDataByChr, indData, centro, &multiWinsizes, error, MAX_GAP,
                                           KDE_SUBSAMPLE, numThreads, WINSIZE_EXPLORE, outfile);
        releaseKDEWinsizeReport(winsizeReport);
        */

        exploreWinsizes(hapDataByChr, freqDataByChr, mapDataByChr,
                        indData, centro, multiWinsizes, error,
                        GLDataByChr, genoFreqDataByChr, USE_GL,
                        MAX_GAP, KDE_SUBSAMPLE, outfile, WEIGHTED, M, mu, numThreads, PHASED, KDE_THIN_STEP, LD_SUBSAMPLE);

        freeRNG();
        return 0;
    }
    else if (AUTO_WINSIZE)
    {
        if(!WEIGHTED){
            try{
                kdeResult = selectWinsize(hapDataByChr, freqDataByChr, mapDataByChr,
                                          indData, centro, winsize, AUTO_WINSIZE_STEP, error,
                                          GLDataByChr, USE_GL,
                                          MAX_GAP, KDE_SUBSAMPLE, outfile, WEIGHTED, genoFreqDataByChr, PHASED, KDE_THIN_STEP,
                                          MAX_WINSIZE);
            }
            catch (...){
                return 2;
            }
        }
        else{
            winsize = selectWinsizeWeighted(variantDensity);
        }
        
        LOG.log("Selected window size:", winsize);
    }

    cout << "Window size: " << winsize << endl;

    if(AUTO_OVERLAP_FRAC){
        OVERLAP_FRAC = selectOverlapFrac(variantDensity, winsize);
        LOG.log("Selected overlap fraction:", OVERLAP_FRAC);
    }


    if(WEIGHTED){
        cerr << "Calculating LD matrix.\n";
        ldDataByChr = calcLDData(hapDataByChr, freqDataByChr, mapDataByChr, genoFreqDataByChr, centro, winsize, MAX_GAP, PHASED, numThreads, LD_SUBSAMPLE);
        if(!PHASED) releaseGenoFreq(genoFreqDataByChr);
        winDataByChr = calcwLODWindows(hapDataByChr, freqDataByChr, mapDataByChr,
                                       GLDataByChr, ldDataByChr,
                                       centro, winsize, error,
                                       MAX_GAP, USE_GL, M, mu, numThreads);
        releaseLDData(ldDataByChr);
    }
    else{
        winDataByChr = calcLODWindows(hapDataByChr, freqDataByChr, mapDataByChr,
                                      GLDataByChr,
                                      centro, winsize, error,
                                      MAX_GAP, USE_GL);
    }
    releaseHapData(hapDataByChr);
    releaseFreqData(freqDataByChr);
    if (USE_GL) releaseGLData(GLDataByChr);

    if (RAW_LOD){
        //Output raw windows
        try { writeWinData(winDataByChr, indData, mapDataByChr, outfile); }
        catch (...) { return 2; }
    }

    if (AUTO_CUTOFF){
        //if ((!AUTO_WINSIZE && !WINSIZE_EXPLORE) || (AUTO_WINSIZE && WEIGHTED) )
        bool cutoffOK = true;
        if(kdeResult == NULL)
        {
            LOD_CUTOFF = selectLODCutoff(winDataByChr, indData, KDE_SUBSAMPLE, makeKDEFilename(outfile, winsize), (KDE_THIN_STEP > 0 ? KDE_THIN_STEP : winsize), winsize, cutoffOK);
        }
        else LOD_CUTOFF = selectLODCutoff(kdeResult, winsize, cutoffOK);

        //A failed KDE used to return -1, which main then used as the cutoff and
        //happily wrote a complete, plausible-looking .roh.bed from it.
        if (!cutoffOK)
        {
            LOG.err("ERROR: Could not select a LOD score cutoff automatically. Stopping.");
            return 2;
        }

        LOG.log("Selected LOD score cutoff:", LOD_CUTOFF);
    }
    else cout << "User defined LOD score cutoff: " << LOD_CUTOFF << "\n";

    cout << "Assembling ROH windows\n";
    //Assemble ROH for each individual in each pop
    ROHLength *rohLength;
    vector< ROHData * > *rohDataByInd = assembleROHWindows(winDataByChr, mapDataByChr, indData,
                                        centro, LOD_CUTOFF, &rohLength, winsize, MAX_GAP, OVERLAP_FRAC, CM);

    releaseWinData(winDataByChr);
    
    if (AUTO_BOUNDS){
        cout << "Fitting " << NCLUST << "-component GMM for size classification\n";
        //A GMM with more components than observations has no solution; the GMM
        //used to throw an int that nothing caught, so the process died with
        //SIGABRT after having already done all the work.  Report it instead.
        if (rohLength == NULL || rohLength->size < NCLUST)
        {
            LOG.err("ERROR: Cannot fit a", NCLUST, false);
            LOG.err("-component GMM to", int(rohLength == NULL ? 0 : rohLength->size), false);
            LOG.err(" ROH.");
            LOG.err("\tNo (or too few) ROH were called, so size classes cannot be chosen automatically.");
            LOG.err("\tCheck --lod-cutoff, or pass --size-bounds to set the class boundaries yourself.");
            return 2;
        }
        try {
            boundSizes = selectSizeClasses(rohLength, NCLUST);
        }
        catch (...) {
            LOG.err("ERROR: GMM size-class fitting failed.");
            LOG.err("\tPass --size-bounds to set the ROH size class boundaries explicitly.");
            return 2;
        }
        LOG.logv("Selected ROH size boundaries = (", boundSizes, false);
        LOG.log(" )");
    }
    else{
        LOG.logv("User provided ROH size boundaries = (", boundSizes, false);
        LOG.log(" )");
    }
    //Output ROH calls to file, one for each individual
    //includes A/B/C/etc size classifications
    cout << "Writing ROH tracts.\n";
    writeROHData(makeROHFilename(outfile), rohDataByInd, mapDataByChr, boundSizes, indData->pop, VERSION, CM);

    if (params->getBoolFlag(ARG_FROH))
    {
        writeFROH(outfile + ".froh.tsv", rohDataByInd, mapDataByChr, boundSizes,
                  indData->pop, centro, CM);
    }

    //Machine-readable record of what this run actually did.  The auto-selected
    //values were previously only prose lines in the .log, which is what the
    //documented workflow for reusing a cutoff asks users to parse by hand.
    {
        vector< pair<string,string> > resolved;
        ostringstream v;
        v << SEED;                        resolved.push_back(make_pair("seed", v.str()));
        v.str(""); v << winsize;          resolved.push_back(make_pair("winsize", v.str()));
        v.str(""); v << OVERLAP_FRAC;     resolved.push_back(make_pair("overlap_frac", v.str()));
        v.str(""); v << LOD_CUTOFF;       resolved.push_back(make_pair("lod_cutoff", v.str()));
        v.str(""); v << KDE_THIN_STEP;    resolved.push_back(make_pair("kde_thin_step", v.str()));
        v.str(""); v << mapDataByChr->size(); resolved.push_back(make_pair("chromosomes_analysed", v.str()));
        v.str(""); v << "[";
        for (unsigned int i = 0; i < boundSizes.size(); i++) { if (i) v << ", "; v << boundSizes[i]; }
        v << "]";                         resolved.push_back(make_pair("size_bounds", v.str()));
        v.str(""); v << (AUTO_CUTOFF ? "true" : "false");  resolved.push_back(make_pair("cutoff_was_automatic", v.str()));
        v.str(""); v << (AUTO_BOUNDS ? "true" : "false");  resolved.push_back(make_pair("bounds_were_automatic", v.str()));
        writeParamsJSON(outfile + ".params.json", params, resolved);
    }

    //centro is read by writeFROH; it used to be deleted before the writers ran.
    delete centro;

    releaseIndData(indData);
    releaseROHLength(rohLength);
    releaseROHData(rohDataByInd);
    if(kdeResult != NULL) releaseKDEResult(kdeResult);
    releaseMapData(mapDataByChr);
    delete params;
    freeRNG();
    cout << "Finished.\n";
    
    #ifdef PTW32_STATIC_LIB
        pthread_win32_process_detach_np();
    #endif
    
    return 0;
}
