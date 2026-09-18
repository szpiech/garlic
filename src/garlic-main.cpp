#include "garlic-errlog.h"
#include "garlic-cli.h"
#include "garlic-options.h"
#include <iostream>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <cmath>

#include "garlic-data.h"
#include "garlic-roh.h"
#include "garlic-kde.h"
#include "param_t.h"
#include "garlic-centromeres.h"

using namespace std;




int main(int argc, char *argv[])
{
//++++++++++CLI handling++++++++++
    int cliStatus = PARAM_OK;
    param_t *params = getCLI(argc, argv, cliStatus);
    //0 success, 1 usage error, 2 runtime error.
    if (params == NULL) return (cliStatus == PARAM_HELP) ? 0 : 1;

    GarlicOptions opt;
    int optStatus = configureFromCommandLine(params, opt, argc, argv);
    if (optStatus != OPTIONS_OK)
    {
        //main owns params on every path now.  The --freq-only exit used to
        //return without deleting it, leaking 249 allocations / 31 KB.
        delete params;
        if (optStatus == OPTIONS_USAGE_ERROR) return 1;
        if (optStatus == OPTIONS_RUNTIME_ERROR) return 2;
        return 0;
    }

    //References rather than copies, so the pipeline below reads and writes the
    //same objects it always did (winsize, LOD_CUTOFF and boundSizes are all
    //reassigned further down when their automatic modes are in use).
    string &outfile = opt.outfile;
    string &tpedfile = opt.tpedfile;
    string &tfamfile = opt.tfamfile;
    string &tglsfile = opt.tglsfile;
    char &TPED_MISSING = opt.TPED_MISSING;
    string &GL_TYPE = opt.GL_TYPE;
    bool &WEIGHTED = opt.WEIGHTED;
    string &mapfile = opt.mapfile;
    bool &CM = opt.CM;
    string &BUILD = opt.BUILD;
    string &centromereFile = opt.centromereFile;
    bool &NO_CENTROMERE = opt.NO_CENTROMERE;
    int &nresample = opt.nresample;
    string &freqfile = opt.freqfile;
    bool &AUTO_FREQ = opt.AUTO_FREQ;
    vector<int> &multiWinsizes = opt.multiWinsizes;
    bool &WINSIZE_EXPLORE = opt.WINSIZE_EXPLORE;
    bool &AUTO_WINSIZE = opt.AUTO_WINSIZE;
    int &AUTO_WINSIZE_STEP = opt.AUTO_WINSIZE_STEP;
    int &winsize = opt.winsize;
    double &LOD_CUTOFF = opt.LOD_CUTOFF;
    bool &AUTO_CUTOFF = opt.AUTO_CUTOFF;
    vector<double> &boundSizes = opt.boundSizes;
    bool &AUTO_BOUNDS = opt.AUTO_BOUNDS;
    int &numThreads = opt.numThreads;
    double &error = opt.error;
    int &MAX_GAP = opt.MAX_GAP;
    double &OVERLAP_FRAC = opt.OVERLAP_FRAC;
    bool &AUTO_OVERLAP_FRAC = opt.AUTO_OVERLAP_FRAC;
    double &mu = opt.mu;
    int &M = opt.M;
    int &NCLUST = opt.NCLUST;
    int &KDE_SUBSAMPLE = opt.KDE_SUBSAMPLE;
    int &LD_SUBSAMPLE = opt.LD_SUBSAMPLE;
    bool &RAW_LOD = opt.RAW_LOD;
    bool &PHASED = opt.PHASED;
    int &KDE_THIN_STEP = opt.KDE_THIN_STEP;
    int &MAX_WINSIZE = opt.MAX_WINSIZE;
    unsigned long int &SEED = opt.SEED;

//++++++++++Datafile reading++++++++++
    centromere *centro;
    if (NO_CENTROMERE)
    {
        //Empty gap table: every centromereStart/End returns 0, and the
        //missing-chromosome warning is suppressed because the absence is
        //deliberate rather than a chromosome-name mismatch.
        centro = new centromere();
        centro->suppressMissingWarnings();
        LOG.log("--no-centromere: treating every chromosome as having no assembled gap.");
    }
    else centro = new centromere(BUILD, centromereFile, DEFAULT_CENTROMERE_FILE);

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

        //Two readers, one set of structures.  Everything after this point is
        //identical for both paths, which is the property worth preserving: a
        //VCF-specific branch anywhere downstream would be a second pipeline.
        string metaSource;
        if (opt.vcffile.compare(DEFAULT_VCF) != 0)
        {
            vector<string> sampleIDs;
            //--gl-type with --vcf reads the FORMAT field directly, so the VCF
            //path fills GLDataByChr itself and readTGLSData is not involved.
            //This is the only way PL and GL can be correct: a VCF normalises
            //them so the CALLED genotype is exactly 0, which the one-value
            //--tgls format cannot represent (see plToError).
            USE_GL = (GL_TYPE.compare(DEFAULT_GL_TYPE) != 0);
            if (USE_GL) GLDataByChr = new vector< GenoLikeData * >;
            loadVCFData(opt.vcffile, numLoci, numInd,
                        &hapDataByChr, &mapDataByChr, &freqDataByChr, &GLDataByChr,
                        nresample, PHASED, AUTO_FREQ, opt.VCF_PASS_ONLY,
                        USE_GL ? GL_TYPE : string("none"), sampleIDs);

            LOG.log("Total loci:", numLoci);

            //A VCF carries sample names but no population and no sex, so the
            //labels start as a single placeholder population and --pop
            //replaces them.  Without --pop every sample shares one label,
            //which is honest -- frequencies really are pooled over all of
            //them -- but it cannot be what a multi-population cohort wants,
            //so it warns.
            indData = initIndData(numInd);
            for (int i = 0; i < numInd; i++)
            {
                indData->indID[i] = sampleIDs[i];
                indData->pop[i]   = "unknown";
                indData->sex[i]   = 0;
            }
            metaSource = opt.vcffile;
            if (opt.popfile.compare(DEFAULT_POP) == 0)
                LOG.err("WARNING: no --pop given, so every sample is labelled 'unknown'. Allele");
            if (opt.popfile.compare(DEFAULT_POP) == 0)
                LOG.err("WARNING: frequencies are pooled over all of them; see --pop.");
        }
        else
        {
            loadTPEDData(tpedfile, numLoci, numInd,
                         &hapDataByChr, &mapDataByChr, &freqDataByChr,
                         TPED_MISSING, nresample, PHASED, AUTO_FREQ);

            LOG.log("Total loci:", numLoci);

            scanIndData3(tfamfile, numInd);
            indData = readIndData3(tfamfile, numInd);
            metaSource = tfamfile;
        }

        //--pop replaces the labels before they are checked, so the pooled-
        //population warning is computed on the labels actually used.
        if (opt.popfile.compare(DEFAULT_POP) != 0)
        {
            applyPopFile(opt.popfile, indData);
            metaSource = opt.popfile;
        }
        checkIndData(indData, metaSource);

        //LOG.log("Population:", popName);
        LOG.log("Total diploid individuals:", numInd);

        if (opt.vcffile.compare(DEFAULT_VCF) == 0 && tglsfile.compare(DEFAULT_TGLS) != 0) {
            GLDataByChr = readTGLSData(tglsfile, numLoci, numInd, mapDataByChr, GL_TYPE, indData);
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
    catch (...) { logCurrentException("reading the input files"); return 2; }

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
        catch (...) { logCurrentException("reading the allele frequency file"); return 2; }
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

        //This return used to skip main's entire cleanup block, so
        //--winsize-multi leaked every structure it had loaded: measured at 420
        //blocks / 1,074,816 bytes on chr21 with leaks(1), and LeakSanitizer in
        //CI reported 441 allocations / 968,474 bytes on the same path.  It is
        //the only early return that reaches here with a full dataset loaded.
        //
        //rohLength, rohDataByInd and kdeResult do not exist yet on this path,
        //so the list is main's normal exit minus those three.  genoFreqDataByChr
        //is NULL unless (!PHASED && WEIGHTED), and releaseGenoFreq dereferences
        //its argument, so it needs the guard.
        releaseHapData(hapDataByChr);
        releaseFreqData(freqDataByChr);
        if (USE_GL) releaseGLData(GLDataByChr);
        if (genoFreqDataByChr != NULL) releaseGenoFreq(genoFreqDataByChr);
        releaseMapData(mapDataByChr);
        releaseIndData(indData);
        delete centro;
        delete params;
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
                logCurrentException("automatic window size selection");
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
        catch (...) { logCurrentException("writing the raw LOD windows"); return 2; }
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
            logCurrentException("GMM size-class fitting");
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
    //The writers throw 0 when they cannot open their output, and nothing
    //caught it: the calls run after every other stage, so an output path that
    //became unwritable turned a complete analysis into SIGABRT with the
    //results discarded.  Reported and carried to the exit status instead, so
    //the cleanup below still runs.
    int writeStatus = 0;
    try
    {
        writeROHData(makeROHFilename(outfile), rohDataByInd, mapDataByChr, boundSizes, indData->pop, VERSION, CM);

        if (params->getBoolFlag(ARG_FROH))
        {
            writeFROH(outfile + ".froh.tsv", rohDataByInd, mapDataByChr, boundSizes,
                      indData->pop, centro, CM);
        }
    }
    catch (...) { logCurrentException("writing the ROH calls"); writeStatus = 2; }

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
        try { writeParamsJSON(outfile + ".params.json", params, resolved); }
        catch (...) { logCurrentException("writing the parameter record"); writeStatus = 2; }
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

    //The PTW32_STATIC_LIB process attach/detach calls that used to bracket
    //main are gone: they are pthreads-win32 specific, and std::thread needs no
    //per-process initialisation on any platform.

    return writeStatus;
}
