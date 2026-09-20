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




//----------------------------------------------------------------------------
// One population's analysis.
//
// Everything from the window size through the written calls, for one set of
// individuals.  Extracted from main unchanged so it can be run once per
// population: the four values a population SELECTS FOR ITSELF -- window size,
// overlap fraction, LOD cutoff and size-class boundaries -- are locals here,
// initialised from the command line.  When the user gave them explicitly the
// local simply keeps that value, which is how an explicit --winsize,
// --overlap-frac, --lod-cutoff or --size-bounds applies to every population.
//
// Ownership: this function releases the PER-POPULATION data it is handed
// (genotypes, frequencies, likelihoods, genotype frequencies) and everything
// it derives from them.  The caller keeps what is SHARED across populations:
// the map, the centromere table, the full individual metadata and the
// command line.
//----------------------------------------------------------------------------
static const int POP_OK    = 0;   //analysed and written
static const int POP_ERROR = 2;   //stop the run, exit 2
static const int POP_DONE  = 3;   //the run finished inside (--winsize-multi)

struct PopResult
{
    int            status;
    int            winsize;
    double         overlapFrac;
    double         lodCutoff;
    vector<double> boundSizes;
    //What the shared sex chromosome was actually called at.  Equal to the
    //autosomal values unless --sexchr-winsize/--sexchr-lod-cutoff changed
    //them, and only meaningful when the run had a shared sex chromosome in
    //it: recorded so the parameter file says what was applied where, rather
    //than leaving a reader to work it out from which flags were set.
    bool           haveSexChr;
    int            sexWinsize;
    double         sexLodCutoff;
};

static PopResult analyzePopulation(const GarlicOptions &opt,
                                   param_t *params,
                                   const string &popLabel,
                                   vector< HapData * > *hapDataByChr,
                                   vector< FreqData * > *freqDataByChr,
                                   vector< MapData * > *mapDataByChr,
                                   vector< GenoLikeData * > *GLDataByChr,
                                   vector< GenoFreqData * > *genoFreqDataByChr,
                                   IndData *indData,
                                   centromere *centro,
                                   bool USE_GL,
                                   double variantDensity,
                                   const vector<ChrRole> *chrRole,
                                   const ExcludedRegions *parRegions,
                                   const ChrLengths *chrLengths)
{
    PopResult res;
    res.status = POP_OK;

    //Empty for a single-population run, so its messages read exactly as they
    //always did.  With several, every failure has to say WHICH population it
    //is talking about -- otherwise the user gets the existing advice with no
    //way to know where to apply it.
    const string where = popLabel.empty() ? string("") : (" [population " + popLabel + "]");

    //Read-only for this population.
    const string &outfile          = opt.outfile;
    const bool   &WEIGHTED         = opt.WEIGHTED;
    const bool   &CM               = opt.CM;
    //Not const: selectWinsizeFromList and exploreWinsizes take it by
    //non-const pointer/reference.  A copy, so one population cannot disturb
    //the next.
    vector<int> multiWinsizes = opt.multiWinsizes;
    const bool   &WINSIZE_EXPLORE  = opt.WINSIZE_EXPLORE;
    const bool   &AUTO_WINSIZE     = opt.AUTO_WINSIZE;
    const int    &AUTO_WINSIZE_STEP = opt.AUTO_WINSIZE_STEP;
    const bool   &AUTO_CUTOFF      = opt.AUTO_CUTOFF;
    const bool   &AUTO_BOUNDS      = opt.AUTO_BOUNDS;
    const int    &numThreads       = opt.numThreads;
    const double &error            = opt.error;
    const int    &MAX_GAP          = opt.MAX_GAP;
    const bool   &AUTO_OVERLAP_FRAC = opt.AUTO_OVERLAP_FRAC;
    const double &mu               = opt.mu;
    const int    &M                = opt.M;
    const int    &NCLUST           = opt.NCLUST;
    const int    &KDE_SUBSAMPLE    = opt.KDE_SUBSAMPLE;
    const int    &LD_SUBSAMPLE     = opt.LD_SUBSAMPLE;
    const bool   &RAW_LOD          = opt.RAW_LOD;
    const bool   &PHASED           = opt.PHASED;
    const int    &KDE_THIN_STEP    = opt.KDE_THIN_STEP;
    const int    &MAX_WINSIZE      = opt.MAX_WINSIZE;
    const double &SEXCHR_LOD_CUTOFF = opt.SEXCHR_LOD_CUTOFF;
    const bool   &SEXCHR_CUTOFF_SET = opt.SEXCHR_CUTOFF_SET;
    const int    &SEXCHR_WINSIZE    = opt.SEXCHR_WINSIZE;

    bool haveSexChr = false;
    if (chrRole != NULL)
        for (unsigned int c = 0; c < chrRole->size(); c++)
            if (chrRole->at(c) == CHR_SEX_SHARED) haveSexChr = true;

    //Chosen per population; the command-line value is the starting point, and
    //stays untouched when it was given explicitly.
    int            winsize      = opt.winsize;
    double         OVERLAP_FRAC = opt.OVERLAP_FRAC;
    double         LOD_CUTOFF   = opt.LOD_CUTOFF;
    vector<double> boundSizes   = opt.boundSizes;

    vector< WinData * > *winDataByChr = NULL;
    vector< LDData * >  *ldDataByChr  = NULL;
    KDEResult           *kdeResult    = NULL;

//++++++++++Pipeline begins++++++++++
    if (WINSIZE_EXPLORE && AUTO_WINSIZE && !WEIGHTED)
    {
        kdeResult = selectWinsizeFromList(hapDataByChr, freqDataByChr, mapDataByChr,
                                          indData, centro, &multiWinsizes, winsize, error,
                                          GLDataByChr, USE_GL,
                                          MAX_GAP, KDE_SUBSAMPLE, outfile, WEIGHTED, genoFreqDataByChr, PHASED, KDE_THIN_STEP,
                                          chrRole);
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
                        MAX_GAP, KDE_SUBSAMPLE, outfile, WEIGHTED, M, mu, numThreads, PHASED, KDE_THIN_STEP, LD_SUBSAMPLE,
                        chrRole);

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
        //Only this population's data: the map, the metadata, the centromere
        //table and the command line belong to the caller.
        releaseHapData(hapDataByChr);
        if (freqDataByChr != NULL) releaseFreqData(freqDataByChr);
        if (USE_GL) releaseGLData(GLDataByChr);
        if (genoFreqDataByChr != NULL) releaseGenoFreq(genoFreqDataByChr);
        res.status = POP_DONE;
        return res;
    }
    else if (AUTO_WINSIZE)
    {
        if(!WEIGHTED){
            try{
                kdeResult = selectWinsize(hapDataByChr, freqDataByChr, mapDataByChr,
                                          indData, centro, winsize, AUTO_WINSIZE_STEP, error,
                                          GLDataByChr, USE_GL,
                                          MAX_GAP, KDE_SUBSAMPLE, outfile, WEIGHTED, genoFreqDataByChr, PHASED, KDE_THIN_STEP,
                                          MAX_WINSIZE, chrRole);
            }
            catch (...){
                logCurrentException("automatic window size selection");
                res.status = POP_ERROR; return res;
            }
        }
        else{
            winsize = selectWinsizeWeighted(variantDensity);
        }
        
        LOG.log("Selected window size:", winsize);
    }

    cout << "Window size: " << winsize << endl;
    //Said once, where the autosomal size is said, because a reader who sees
    //only one number will assume it applied to everything.
    if (SEXCHR_WINSIZE > 0 && haveSexChr)
        LOG.log("Window size on the shared sex chromosome:", SEXCHR_WINSIZE);

    if(AUTO_OVERLAP_FRAC){
        OVERLAP_FRAC = selectOverlapFrac(variantDensity, winsize);
        LOG.log("Selected overlap fraction:", OVERLAP_FRAC);
    }


    if(WEIGHTED){
        cerr << "Calculating LD matrix.\n";
        ldDataByChr = calcLDData(hapDataByChr, freqDataByChr, mapDataByChr, genoFreqDataByChr, centro, winsize, MAX_GAP, PHASED, numThreads, LD_SUBSAMPLE,
                                 SEXCHR_WINSIZE, chrRole, &(indData->zygo));
        if(!PHASED) releaseGenoFreq(genoFreqDataByChr);
        winDataByChr = calcwLODWindows(hapDataByChr, freqDataByChr, mapDataByChr,
                                       GLDataByChr, ldDataByChr,
                                       centro, winsize, error,
                                       MAX_GAP, USE_GL, M, mu, numThreads,
                                       SEXCHR_WINSIZE, chrRole);
        releaseLDData(ldDataByChr);
    }
    else{
        winDataByChr = calcLODWindows(hapDataByChr, freqDataByChr, mapDataByChr,
                                      GLDataByChr,
                                      centro, winsize, error,
                                      MAX_GAP, USE_GL,
                                      SEXCHR_WINSIZE, chrRole);
    }
    releaseHapData(hapDataByChr);
    if (freqDataByChr != NULL) releaseFreqData(freqDataByChr);
    if (USE_GL) releaseGLData(GLDataByChr);

    //A heterogametic individual's windows on the shared sex chromosome are
    //sums of lod() over genotypes it does not have, which is exactly 0 rather
    //than MISSING -- so without this they would be a spike at zero in the
    //density and would be CALLED under any negative cutoff.  Masked before
    //--raw-lod writes, so that file says missing rather than zero too.
    {
        long long nMasked = maskIneligibleWindows(winDataByChr, indData, chrRole);
        if (nMasked > 0) LOG.log("Windows excluded as not diploid in the individual:", nMasked);
    }

    if (RAW_LOD){
        //Output raw windows
        try { writeWinData(winDataByChr, indData, mapDataByChr, outfile); }
        catch (...) { logCurrentException("writing the raw LOD windows"); res.status = POP_ERROR; return res; }
    }

    if (AUTO_CUTOFF){
        //if ((!AUTO_WINSIZE && !WINSIZE_EXPLORE) || (AUTO_WINSIZE && WEIGHTED) )
        bool cutoffOK = true;
        if(kdeResult == NULL)
        {
            LOD_CUTOFF = selectLODCutoff(winDataByChr, indData, KDE_SUBSAMPLE, makeKDEFilename(outfile, winsize), (KDE_THIN_STEP > 0 ? KDE_THIN_STEP : winsize), winsize, cutoffOK,
                                         chrRole, CHR_AUTOSOME);
        }
        else LOD_CUTOFF = selectLODCutoff(kdeResult, winsize, cutoffOK);

        //A failed KDE used to return -1, which main then used as the cutoff and
        //happily wrote a complete, plausible-looking .roh.bed from it.
        if (!cutoffOK)
        {
            LOG.err("ERROR: Could not select a LOD score cutoff automatically" + where + ". Stopping.");
            res.status = POP_ERROR; return res;
        }

        LOG.log("Selected LOD score cutoff:", LOD_CUTOFF);
        //Estimated on the autosomes and applied to the sex chromosome, which
        //is what Cotter et al. (2024) did and for the same reason.  Report
        //what the sex chromosome alone would have given, so the assumption is
        //visible rather than implicit.
        const int sexWin = (SEXCHR_WINSIZE > 0) ? SEXCHR_WINSIZE : winsize;
        reportSexChrLODCutoff(winDataByChr, indData, chrRole,
                              (KDE_THIN_STEP > 0 ? KDE_THIN_STEP : sexWin), sexWin,
                              SEXCHR_CUTOFF_SET ? SEXCHR_LOD_CUTOFF : LOD_CUTOFF,
                              SEXCHR_CUTOFF_SET);
    }
    else cout << "User defined LOD score cutoff: " << LOD_CUTOFF << "\n";
    if (SEXCHR_CUTOFF_SET && haveSexChr)
        LOG.log("User defined LOD score cutoff on the shared sex chromosome:", SEXCHR_LOD_CUTOFF);

    cout << "Assembling ROH windows\n";
    //Assemble ROH for each individual in each pop
    ROHLength *rohLength;
    vector< ROHData * > *rohDataByInd = assembleROHWindows(winDataByChr, mapDataByChr, indData,
                                        centro, LOD_CUTOFF, &rohLength, winsize, MAX_GAP, OVERLAP_FRAC, CM,
                                        chrRole, SEXCHR_LOD_CUTOFF, SEXCHR_CUTOFF_SET, SEXCHR_WINSIZE);

    releaseWinData(winDataByChr);
    
    if (AUTO_BOUNDS){
        cout << "Fitting " << NCLUST << "-component GMM for size classification\n";
        //A GMM with more components than observations has no solution; the GMM
        //used to throw an int that nothing caught, so the process died with
        //SIGABRT after having already done all the work.  Report it instead.
        if (rohLength == NULL || rohLength->size < NCLUST)
        {
            LOG.err("ERROR:" + where + " Cannot fit a", NCLUST, false);
            LOG.err("-component GMM to", int(rohLength == NULL ? 0 : rohLength->size), false);
            LOG.err(" ROH.");
            LOG.err("\tNo (or too few) ROH were called, so size classes cannot be chosen automatically.");
            LOG.err("\tCheck --lod-cutoff, or pass --size-bounds to set the class boundaries yourself.");
            //rohDataByInd and rohLength already exist here; with a loop over
            //populations an error return that skipped them would leak one
            //population's tracts per failure.
            releaseROHLength(rohLength); releaseROHData(rohDataByInd);
            res.status = POP_ERROR; return res;
        }
        try {
            boundSizes = selectSizeClasses(rohLength, NCLUST);
        }
        catch (...) {
            logCurrentException("GMM size-class fitting");
            LOG.err("\tPass --size-bounds to set the ROH size class boundaries explicitly.");
            //rohDataByInd and rohLength already exist here; with a loop over
            //populations an error return that skipped them would leak one
            //population's tracts per failure.
            releaseROHLength(rohLength); releaseROHData(rohDataByInd);
            res.status = POP_ERROR; return res;
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
                  indData, centro, CM, popLabel, opt.POOLED, chrRole, parRegions,
                  opt.FROH_DENOM, chrLengths);
        }
    }
    catch (...) { logCurrentException("writing the ROH calls"); writeStatus = 2; }
    //Per-population state, released here rather than by the caller: with more
    //than one population these would otherwise accumulate across the run.
    releaseROHLength(rohLength);
    releaseROHData(rohDataByInd);
    if (kdeResult != NULL) releaseKDEResult(kdeResult);

    res.winsize     = winsize;
    res.overlapFrac = OVERLAP_FRAC;
    res.lodCutoff   = LOD_CUTOFF;
    res.boundSizes  = boundSizes;
    res.haveSexChr  = haveSexChr;
    res.sexWinsize  = (SEXCHR_WINSIZE > 0) ? SEXCHR_WINSIZE : winsize;
    res.sexLodCutoff = SEXCHR_CUTOFF_SET ? SEXCHR_LOD_CUTOFF : LOD_CUTOFF;
    if (writeStatus != 0) res.status = POP_ERROR;
    return res;
}


int main(int argc, char *argv[])
{
//++++++++++CLI handling++++++++++
    int cliStatus = PARAM_OK;
    param_t *params = getCLI(argc, argv, cliStatus);
    //0 success, 1 usage error, 2 runtime error.
    if (params == NULL) return (cliStatus == PARAM_HELP) ? 0 : 1;

    //Whatever a VCF's ##contig lines state.  Empty for TPED input and for a
    //VCF that declares no lengths; the lowest-precedence source of the
    //lengths --froh-denominator assembly needs, and a cross-check on --build
    //whether or not that mode is in use.
    ChrLengths headerLengths;

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
    //The options the ANALYSIS uses are read inside analyzePopulation, from
    //its own GarlicOptions reference; main keeps only what it needs to load
    //the data and write the run record.
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
    bool &AUTO_WINSIZE = opt.AUTO_WINSIZE;
    bool &AUTO_CUTOFF = opt.AUTO_CUTOFF;
    bool &AUTO_BOUNDS = opt.AUTO_BOUNDS;
    bool &AUTO_OVERLAP_FRAC = opt.AUTO_OVERLAP_FRAC;
    bool &PHASED = opt.PHASED;
    int &KDE_THIN_STEP = opt.KDE_THIN_STEP;

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
    //Filled once --pop has had its say; see enumeratePopulations.  Declared
    //here rather than at the point of use because that is inside the reading
    //try block, and the parameter record below needs it.
    vector< pair<string, int> > populations;
    vector< HapData * > *hapDataByChr = NULL;
    vector< FreqData * > *freqDataByChr = NULL;
    vector< GenoFreqData * > *genoFreqDataByChr = NULL;
    vector< GenoLikeData * > *GLDataByChr = NULL;
    vector< GenMapScaffold *> *scaffoldMapByChr = NULL;
    //One frequency set per population, when --freq-file supplied them.
    vector< vector< FreqData * >* > *fileFreq = NULL;
    //The same, when garlic computed them itself.
    vector< vector< FreqData * >* > popFreq;
    bool USE_GL = false;
    //Parsed before the read, not after: a haploid genotype is an error on an
    //autosome and the normal encoding of a hemizygous call on a sex
    //chromosome, and the VCF reader has to tell them apart while it is still
    //reading.  These are the flags alone; the role table proper is built from
    //the data further down, and it is what the refusals are based on.
    bool AUTOSOMES_ONLY = params->getBoolFlag(ARG_AUTOSOMES_ONLY);
    string sexSystemArg = params->getStringFlag(ARG_SEX_SYSTEM);
    int sexSystem = SEX_SYSTEM_UNSET;
    if (sexSystemArg.compare(DEFAULT_SEX_SYSTEM) != 0)
    {
        sexSystem = parseSexSystem(sexSystemArg);
        if (sexSystem < 0)
        {
            LOG.err("ERROR:", ARG_SEX_SYSTEM, false);
            LOG.err(" must be xy, zw or none, not", sexSystemArg, false);
            LOG.err(".");
            return 1;
        }
    }
    vector<string> sexChrNames, degenerateChrNames, haploidChrNames;
    if (params->isFlagSet(ARG_SEX_CHR)) sexChrNames = params->getStringListFlag(ARG_SEX_CHR);
    if (params->isFlagSet(ARG_SEX_CHR_DEGENERATE)) degenerateChrNames = params->getStringListFlag(ARG_SEX_CHR_DEGENERATE);
    if (params->isFlagSet(ARG_HAPLOID_CHR)) haploidChrNames = params->getStringListFlag(ARG_HAPLOID_CHR);
    map<string, ChrRole> declaredRoles;
    if (declaredRolesByKey(declaredRoles, sexSystem, sexChrNames, degenerateChrNames,
                           haploidChrNames, isHumanBuild(BUILD)) < 0) return 1;

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
                        USE_GL ? GL_TYPE : string("none"), sampleIDs, &declaredRoles,
                        &headerLengths);

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

        //Recorded before anything acts on it.  garlic still pools every
        //individual for allele frequencies at this commit -- the warning in
        //checkIndData still applies and is still emitted -- but which
        //populations are present, and in what order, is now on the record.
        populations = enumeratePopulations(indData);
        {
            stringstream ps;
            ps << "Populations found: " << populations.size() << " (";
            for (unsigned int i = 0; i < populations.size(); i++)
            {
                if (i) ps << ", ";
                ps << populations[i].first << ": " << populations[i].second;
            }
            ps << ")";
            LOG.log(ps.str());
        }
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

    //--pool-populations restores the pre-per-population behaviour: one
    //analysis over everyone, pooled frequencies, unlabelled output names.
    //Implemented by making the loop take the single-population path rather
    //than by rewriting the population list, so the run record still reports
    //which populations were actually present.
    bool POOL = params->getBoolFlag(ARG_POOL);
    opt.POOLED = POOL;
    if (POOL && populations.size() > 1)
    {
        LOG.log("--pool-populations: analysing all", indData->nind, false);
        LOG.log(" individuals as one population.");
        LOG.err("WARNING: --pool-populations pools allele frequencies across", int(populations.size()), false);
        LOG.err(" populations,");
        LOG.err("\twhich inflates heterozygosity relative to any one of them and biases");
        LOG.err("\tthe LOD scores for all of them.");
    }
    bool singlePop = (populations.size() <= 1) || POOL;

    //A population's allele frequencies are estimated from that population
    //alone, so a small one estimates them coarsely: with n individuals the
    //only attainable values are multiples of 1/(2n).  Worth saying, because
    //the split is silent otherwise.
    if (!singlePop)
    {
        int smallest = indData->nind;
        for (unsigned int k = 0; k < populations.size(); k++)
            if (populations[k].second < smallest) smallest = populations[k].second;
        if (smallest < 5)
        {
            LOG.err("WARNING: the smallest population has", smallest, false);
            LOG.err(" individuals.");
            LOG.err("\tAllele frequencies are estimated within each population, so they");
            LOG.err("\tresolve only to multiples of 1/(2n) there.");
        }
        //Column 1 of a TFAM is garlic's population label, but PLINK calls it
        //the FAMILY ID and it is often one family -- sometimes one individual
        //-- per value.  Averaging under two individuals per population is a
        //much better sign of that than of real populations.
        if (populations.size() * 2 > (unsigned int)indData->nind)
        {
            LOG.err("WARNING:", int(populations.size()), false);
            LOG.err(" populations for", indData->nind, false);
            LOG.err(" individuals.");
            LOG.err("\tgarlic takes the population label from column 1 of the TFAM, which");
            LOG.err("\tPLINK calls the family ID; if those are family or sample identifiers");
            LOG.err("\trather than populations, this run is analysing each one separately.");
            LOG.err("\tUse --pop to supply the population labels.");
        }
    }

    //Which individuals belong to each population, in file order.
    vector< vector<int> > popIndex(populations.size());
    for (int i = 0; i < indData->nind; i++)
        for (unsigned int k = 0; k < populations.size(); k++)
            if (indData->pop[i].compare(populations[k].first) == 0)
            { popIndex[k].push_back(i); break; }

    //Chromosome names are matched case-insensitively and ignoring a "chr"
    //prefix from here on, so two names that differ only that way would be
    //indistinguishable.  Checked against the data as read, before --chr can
    //hide one of them.
    if (checkChrKeyCollisions(mapDataByChr) < 0) return 1;

    //---- what each chromosome is ------------------------------------------
    //
    //garlic does not decide this from a name.  --sex-system supplies the one
    //bit the data cannot (which sex code is heterogametic), the conventional
    //names follow from it, and anything that looks sex-linked but has not
    //been accounted for stops the run.  A hemizygous genotype is written as a
    //homozygous call, so the alternative to stopping is a plausible FROH that
    //is wrong, and a warning does not prevent that number being published.
    SexModel sexModel;
    if (buildSexModel(sexModel, mapDataByChr, indData, sexSystem,
                      sexChrNames, degenerateChrNames, haploidChrNames,
                      isHumanBuild(BUILD), AUTOSOMES_ONLY) < 0) return 1;

    if (sexSystem != SEX_SYSTEM_UNSET)
        LOG.log("Sex determination system:", sexSystemName(sexSystem));

    //Every role other than an autosome is dropped here.  A degenerate sex
    //chromosome is hemizygous in one sex and absent in the other, a haploid
    //chromosome is hemizygous in everyone, and a pseudoautosomal region is
    //not called by choice -- so none of them can contribute a run of
    //homozygosity to anybody, and keeping them would only feed the LOD score
    //density windows that are not evidence of anything.
    //
    //The SHARED sex chromosome is kept unless --autosomes-only asked for it
    //to go: it is diploid in the homogametic sex, which can be autozygous on
    //it, and the individuals who cannot are handled per individual rather
    //than by dropping the chromosome for everybody.
    {
        vector<string> keepNames;
        int nDropped[5] = {0, 0, 0, 0, 0};
        for (unsigned int chr = 0; chr < mapDataByChr->size(); chr++)
        {
            ChrRole r = sexModel.role[chr];
            bool keep = (r == CHR_AUTOSOME) || (r == CHR_SEX_SHARED && !AUTOSOMES_ONLY);
            if (keep) { keepNames.push_back(mapDataByChr->at(chr)->chr); continue; }
            nDropped[int(r)]++;
            LOG.log("Dropping", mapDataByChr->at(chr)->chr, false);
            if (r == CHR_SEX_SHARED)
                LOG.log(": shared sex chromosome.");
            else if (r == CHR_SEX_DEGENERATE)
                LOG.log(": carried only by the heterogametic sex, so no individual can be autozygous on it.");
            else if (r == CHR_HAPLOID)
                LOG.log(": haploid in every individual.");
            else
                LOG.log(": pseudoautosomal under PLINK's human coding; see --par.");
        }

        int ndrop = nDropped[1] + nDropped[2] + nDropped[3] + nDropped[4];
        if (ndrop > 0)
        {
            if (keepNames.empty())
            {
                LOG.err("ERROR: nothing left to analyse: no chromosome in the data can carry a run");
                LOG.err("\tof homozygosity in any individual.");
                return 1;
            }
            int nkept = filterChromosomes(keepNames, &mapDataByChr, &hapDataByChr,
                                          &freqDataByChr, &GLDataByChr, USE_GL);
            if (nkept < 0) return 1;
            LOG.log("Chromosomes dropped:", ndrop, false);
            LOG.log(", analysed:", nkept);
            //Rebuilt from the keyed copy rather than re-derived: filtering
            //moves the positional vector, and a role must not be able to come
            //out differently after a filter than it did before one.
            sexModel.rebuild(mapDataByChr);
        }
        else if (AUTOSOMES_ONLY)
        {
            LOG.log("--autosomes-only: every chromosome in the data is an autosome, nothing to drop.");
        }
    }

    //---- pseudoautosomal regions -------------------------------------------
    //
    //Dropped before anything reads the genotypes: they are diploid in both
    //sexes, so a heterozygous call there is real rather than impossible, and
    //leaving them in would both inflate the sex check's count of impossible
    //calls and bias the allele frequency of the heterogametic sex -- dropping
    //heterozygotes while keeping homozygotes is not a random thinning.
    ExcludedRegions par;
    {
        bool any = false, suppressed = false, fromBuild = false;
        if (params->isFlagSet(ARG_PAR))
        {
            vector<string> specs = params->getStringListFlag(ARG_PAR);
            //"--par none" is how a run on a human build says it does not want
            //the built-in regions: the alternative would be another flag whose
            //only job is to turn one table off.
            if (specs.size() == 1 && specs[0].compare("none") == 0) suppressed = true;
            else
            {
                if (parsePARSpecs(specs, par) < 0) return 1;
                any = true;
            }
        }
        string parFile = params->getStringFlag(ARG_PAR_FILE);
        if (parFile.compare(DEFAULT_PAR_FILE) != 0)
        {
            if (readPARFile(parFile, par) < 0) return 1;
            any = true;
        }

        //The assembly defines these, so --build supplies them.  Only when the
        //user has named none of their own: a run that says --par means the
        //regions it names, not those plus a table it did not ask for.
        if (!any && !suppressed && sexModel.anyOfRole(CHR_SEX_SHARED) && isHumanBuild(BUILD))
        {
            //Only for the human X, and filed under the name this data set
            //calls it -- chrX in one file and 23 in the next.  A shared sex
            //chromosome under a human build that is not the X is someone
            //else's chromosome in human coordinates, and guessing there would
            //be worse than doing nothing.
            string sharedName;
            for (unsigned int chr = 0; chr < mapDataByChr->size(); chr++)
            {
                if (sexModel.role[chr] != CHR_SEX_SHARED) continue;
                string key = canonChrKey(mapDataByChr->at(chr)->chr);
                if (key.compare("x") == 0 || key.compare("23") == 0)
                    sharedName = mapDataByChr->at(chr)->chr;
            }
            vector<Interval> builtin;
            if (!sharedName.empty() && builtinPAR(BUILD, builtin))
            {
                for (unsigned int i = 0; i < builtin.size(); i++)
                    par.add(sharedName, builtin[i].start, builtin[i].end);
                any = true; fromBuild = true;
                LOG.log("Pseudoautosomal regions for", BUILD, false);
                LOG.log(" applied to", sharedName, false);
                LOG.log("; see centromeres/par_regions.txt. Pass --par none to keep them.");
            }
        }

        if (any)
        {
            if (par.finalise() < 0) return 1;
            if (!sexModel.anyOfRole(CHR_SEX_SHARED))
            {
                LOG.err("ERROR: there is no shared sex chromosome to hold a pseudoautosomal region.");
                LOG.err("\tSee --sex-system and --sex-chr.");
                return 1;
            }
            if (dropExcludedSites(&mapDataByChr, &hapDataByChr, &freqDataByChr,
                                  &GLDataByChr, par, sexModel, USE_GL, PHASED,
                                  fromBuild) < 0) return 1;
        }
    }

    //---- who is diploid on the shared sex chromosome -----------------------
    //
    //Runs on every dataset that has one, because it is one pass over that
    //chromosome and it is the only thing that can check a declared sex, infer
    //an unrecorded one, and find the genotypes that the answer makes
    //impossible.  It also does the recode, which has to happen before allele
    //frequencies are touched: a hemizygous call becomes a half call, which
    //contributes its one observed allele and no genotype, and that is what
    //makes the frequency come out as "the heterogametic sex contributes one
    //allele and the homogametic two" with no change to the frequency code.
    if (sexModel.anyOfRole(CHR_SEX_SHARED))
    {
        double HET_RATE_LO = 0.02, HET_RATE_HI = 0.10;
        if (params->isFlagSet(ARG_HET_RATE_BOUNDS))
        {
            vector<double> b = params->getDoubleListFlag(ARG_HET_RATE_BOUNDS);
            if (b.size() != 2 || b[0] < 0 || b[1] > 1 || b[0] > b[1])
            {
                LOG.err("ERROR:", ARG_HET_RATE_BOUNDS, false);
                LOG.err(" takes two rates, <lo> <hi>, with 0 <= lo <= hi <= 1.");
                return 1;
            }
            HET_RATE_LO = b[0]; HET_RATE_HI = b[1];
        }

        if (runSexCheck(hapDataByChr, mapDataByChr, indData, sexModel,
                        HET_RATE_LO, HET_RATE_HI, outfile, false) < 0) return 1;

        //Only the chromosomes whose genotypes just changed.  The readers
        //compute frequencies while reading, before anything knows what a
        //chromosome is, so the sex chromosome's are wrong and every autosome's
        //is right -- and --resample draws random numbers, so recomputing an
        //autosome would not even be a no-op.
        recomputeFreqForRole(hapDataByChr, freqDataByChr, sexModel, CHR_SEX_SHARED, nresample);

        if (!AUTO_FREQ)
        {
            LOG.err("WARNING: --freq-file supplies the allele frequencies of the shared sex");
            LOG.err("\tchromosome as well. They must already be computed with the heterogametic");
            LOG.err("\tsex contributing one allele and the homogametic sex two; garlic cannot");
            LOG.err("\ttell from the file whether they were.");
        }
    }


//++++++++++Allele frequencies++++++++++
    if (AUTO_FREQ)
    {
        //cout << "Calculating allele frequencies\n";
        //freqDataByChr = calcFreqData2(hapDataByChr, nresample);

        string freqOutfile = outfile;
        freqOutfile += ".freq";
        if (singlePop)
        {
            writeFreqData(freqOutfile, freqDataByChr, mapDataByChr, indData);
        }
        else
        {
            //One column per population, which is the format --freq-file reads,
            //so a run's own frequency file reproduces that run.  A single
            //pooled column would describe an analysis that did not happen.
            //Computed here rather than in the loop because the file is written
            //against the UNFILTERED map, and filtering is per population.
            vector<string> popNames;
            for (unsigned int k = 0; k < populations.size(); k++)
            {
                popNames.push_back(populations[k].first);
                popFreq.push_back(calcFreqDataForIndices(hapDataByChr, popIndex[k], nresample));
            }
            try { writeFreqDataWide(freqOutfile, popFreq, popNames, mapDataByChr); }
            catch (...) { logCurrentException("writing the allele frequency file"); return 2; }
        }
    }
    else //(!AUTO_FREQ)
    {
        cout << "Loading user provided allele frequencies from " << freqfile << "\n";
        //One request, naming the populations this run will analyse.  With a
        //single population that is an empty list, which asks for one pooled
        //set and accepts the one-column file every earlier version wrote.
        //With several it asks for a column each, by name, so the file's column
        //ORDER does not have to match the order they appear in the TFAM.
        try {
            vector<string> want;
            if (populations.size() > 1 && !params->getBoolFlag(ARG_POOL))
                for (unsigned int k = 0; k < populations.size(); k++)
                    want.push_back(populations[k].first);
            vector< vector< FreqData * >* > *sets =
                readFreqData(freqfile, mapDataByChr, want);
            if (!want.empty()) fileFreq = sets;
            else { freqDataByChr = sets->at(0); sets->clear(); delete sets; }
        }
        catch (...)
        {
            //Not new, but newly reachable: a wide file adds ways to fail here
            //(a missing column, a population named in the TFAM and not in the
            //file), and the genotypes are already loaded by this point.
            logCurrentException("reading the allele frequency file");
            releaseHapData(hapDataByChr);
            if (freqDataByChr != NULL) releaseFreqData(freqDataByChr);
            if (USE_GL) releaseGLData(GLDataByChr);
            if (genoFreqDataByChr != NULL) releaseGenoFreq(genoFreqDataByChr);
            if (scaffoldMapByChr != NULL) releaseGenMapScaffold(scaffoldMapByChr);
            releaseMapData(mapDataByChr);
            releaseIndData(indData);
            delete centro; delete params; freeRNG();
            return 2;
        }
    }

    //Site filtering happens PER POPULATION, inside the loop below, because
    //each population filters on its own allele frequencies.  A site
    //monomorphic in one population can be polymorphic in another, and a
    //population analysed alongside others has to see the same locus set it
    //would see alone -- which is the difference between "its frequencies are
    //its own" and "its analysis is its own".  The scaffold map therefore has
    //to stay alive until the loop is done.

    //chrCoordList->clear();
    //delete chrCoordList;

    vector<string> keepChr = params->getStringListFlag(ARG_CHR);
    if (params->isFlagSet(ARG_CHR))
    {
        int nkept = filterChromosomes(keepChr, &mapDataByChr, &hapDataByChr,
                                      &freqDataByChr, &GLDataByChr, USE_GL);
        if (nkept < 0) return 1;
        sexModel.rebuild(mapDataByChr);
        LOG.loga("Restricted to chromosomes:", &(keepChr[0]), int(keepChr.size()));
        LOG.log("Chromosomes analysed:", nkept);
        if (!PHASED && WEIGHTED)
        {
            releaseGenoFreq(genoFreqDataByChr);
            genoFreqDataByChr = calculateGenoFreq(hapDataByChr);
        }
    }

    //---- chromosome lengths ------------------------------------------------
    //
    //Resolved here because it needs three things the command line alone does
    //not have: which chromosomes survived filtering, what a VCF's header said,
    //and where each chromosome's last marker lies.
    //
    //Precedence: --chr-lengths, then --build, then the VCF header.  The user's
    //own file is the most specific statement.  --build beats the header
    //because it already supplies the centromeres and the pseudoautosomal
    //regions, and taking lengths from a different source than those would be
    //incoherent.  The header is the file's own claim and is the weakest, but
    //it is worth having: it is the only length a VCF from a non-human assembly
    //carries, and it checks --build for free.
    ChrLengths chrLengths;
    {
        ChrLengths buildLengths;
        const bool haveBuild = builtinChrLengths(BUILD, buildLengths);
        const bool haveFile  = (opt.CHR_LENGTHS_FILE.compare(DEFAULT_CHR_LENGTHS) != 0);

        if (haveFile)
        {
            if (readChrLengthsFile(opt.CHR_LENGTHS_FILE, chrLengths) < 0) return 1;
        }
        else if (haveBuild) chrLengths = buildLengths;
        else if (!headerLengths.empty())
        {
            chrLengths = headerLengths;
            chrLengths.setSource("VCF ##contig header");
        }
        if (!chrLengths.empty()) LOG.log("Chromosome lengths:", chrLengths.source());

        //Two assertions about the same thing.  A hg19 VCF analysed under
        //--build hg38 is a real and otherwise undetectable mistake, and the
        //header is the only witness to it.  Compared only where both speak,
        //and only for chromosomes the run actually uses -- a decoy or an alt
        //contig the build has never heard of is not a disagreement.
        if (haveBuild && !headerLengths.empty())
        {
            for (unsigned int chr = 0; chr < mapDataByChr->size(); chr++)
            {
                const string &name = mapDataByChr->at(chr)->chr;
                pos_t fromHeader = headerLengths.get(name);
                pos_t fromBuild  = buildLengths.get(name);
                if (fromHeader == 0 || fromBuild == 0 || fromHeader == fromBuild) continue;
                ostringstream ss;
                ss << "ERROR: the VCF header gives " << name << " a length of " << (long long)fromHeader
                   << ", and " << BUILD << " has " << (long long)fromBuild << ".";
                LOG.err(ss.str());
                LOG.err("\tOne of them is not this assembly.");
                return 1;
            }
        }

        //A marker past the end of its chromosome means the lengths and the
        //coordinates are from different assemblies.  Under assembly that makes
        //the denominator wrong, so it stops the run; under analyzed the
        //lengths feed nothing, so it is a warning -- but still worth saying,
        //because the centromere and pseudoautosomal coordinates --build
        //supplied are from that same wrong assembly.
        if (!chrLengths.empty())
        {
            string over;
            for (unsigned int chr = 0; chr < mapDataByChr->size(); chr++)
            {
                MapData *md = mapDataByChr->at(chr);
                if (md->nloci < 1) continue;
                pos_t len = chrLengths.get(md->chr);
                if (len == 0 || md->physicalPos[md->nloci - 1] <= len) continue;
                if (!over.empty()) over += ", ";
                ostringstream ss;
                ss << md->chr << " (last marker " << (long long)md->physicalPos[md->nloci - 1]
                   << " > " << (long long)len << ")";
                over += ss.str();
            }
            if (!over.empty())
            {
                if (opt.FROH_DENOM == FROH_ASSEMBLY)
                {
                    LOG.err("ERROR: markers lie past the end of the chromosome according to");
                    LOG.err("\t" + chrLengths.source() + ":", over);
                    LOG.err("\tThose lengths are the denominator, so this run would divide by the");
                    LOG.err("\twrong numbers. Check --build, or give --chr-lengths.");
                    return 1;
                }
                LOG.err("WARNING: markers lie past the end of the chromosome according to");
                LOG.err("\t" + chrLengths.source() + ":", over);
                LOG.err("\tThe coordinates and that assembly do not agree, which also means the");
                LOG.err("\tcentromere and pseudoautosomal positions in use are the wrong ones.");
            }
        }

        //Under assembly the lengths are load-bearing, so every analysed
        //chromosome must have one.  Falling back to a chromosome's marker span
        //for the ones that are missing would produce a denominator that is
        //neither convention, and a warning about it does not survive into a
        //figure.
        if (opt.FROH_DENOM == FROH_ASSEMBLY)
        {
            string missing;
            for (unsigned int chr = 0; chr < mapDataByChr->size(); chr++)
            {
                MapData *md = mapDataByChr->at(chr);
                if (md->nloci < 2) continue;   //not in any denominator anyway
                if (sexModel.role[chr] != CHR_AUTOSOME &&
                    sexModel.role[chr] != CHR_SEX_SHARED) continue;
                if (chrLengths.get(md->chr) > 0) continue;
                if (!missing.empty()) missing += ", ";
                missing += md->chr;
            }
            if (!missing.empty())
            {
                LOG.err("ERROR:", ARG_FROH_DENOM, false);
                LOG.err(" assembly needs the length of every analysed chromosome, and");
                if (chrLengths.empty()) LOG.err("\tnone is available. Supply them with one of:");
                else
                {
                    LOG.err("\t" + chrLengths.source() + " does not give one for:", missing);
                    LOG.err("\tSupply them with one of:");
                }
                LOG.err(string("\t  ") + ARG_CHR_LENGTHS + " <file>   two columns, chromosome and length, or a .fai");
                LOG.err(string("\t  ") + ARG_BUILD + " <assembly>    for a human assembly garlic knows");
                LOG.err("\t  a VCF whose header carries ##contig=<ID=...,length=...>");
                return 1;
            }
        }
    }

    //Checked here rather than in the option parser: whether the run has a
    //shared sex chromosome depends on the data and on --autosomes-only and
    //--chr, none of which the command line alone can answer.  Silently
    //ignoring the flags would leave the user believing the sex chromosome was
    //called at the value they gave.
    if ((opt.SEXCHR_CUTOFF_SET || opt.SEXCHR_WINSIZE > 0) &&
        !sexModel.anyOfRole(CHR_SEX_SHARED))
    {
        LOG.err("ERROR:", opt.SEXCHR_CUTOFF_SET ? ARG_SEXCHR_LOD_CUTOFF : ARG_SEXCHR_WINSIZE, false);
        LOG.err(" was given, but no shared sex chromosome is being analysed.");
        LOG.err("\tSee --sex-system and --sex-chr; --autosomes-only and --chr can also");
        LOG.err("\tremove it from the run.");
        return 1;
    }

    //++++++++++Pipeline begins++++++++++
    //Each population is analysed on its own: its own allele frequencies, its
    //own window size, its own LOD cutoff and its own size classes.  Anything
    //the user gave explicitly is not re-derived, so --winsize, --lod-cutoff,
    //--overlap-frac and --size-bounds apply to every population -- that falls
    //out of analyzePopulation initialising its locals from the option struct.
    //
    //Sequentially, not together: the per-population genotype matrix is the
    //largest structure in the run, and holding one at a time bounds peak
    //memory by the largest population rather than by their sum.  The inner
    //stages are already threaded, so nothing is left idle.

    int writeStatus = 0;
    unsigned int npass = POOL ? 1 : (unsigned int)populations.size();
    for (unsigned int k = 0; k < npass && writeStatus == 0; k++)
    {
        const string &popName = populations[k].first;

        GarlicOptions popOpt = opt;
        //Naming: unchanged for one population, so every existing invocation
        //writes exactly the files it always did.  With more than one, the
        //label goes in the middle -- <out>.<POP>.roh.bed -- because appending
        //it after the extension would break every tool that dispatches on it.
        if (!singlePop)
        {
            popOpt.outfile = opt.outfile + "." + popName;
            LOG.log("");
            LOG.log("Population:", popName, false);
            LOG.log(" (", int(popIndex[k].size()), false);
            LOG.log(" individuals)");

            //A fresh stream per population, derived from the run's seed, so a
            //population's result does not depend on how many populations
            //preceded it in the file.  --seed still reproduces the whole run.
            initRNG(opt.SEED + (unsigned long int)(k));
        }

        vector< HapData * >      *pHap = hapDataByChr;
        vector< FreqData * >     *pFreq = freqDataByChr;
        vector< GenoLikeData * > *pGL = GLDataByChr;
        vector< GenoFreqData * > *pGF = genoFreqDataByChr;
        vector< MapData * >      *pMap = mapDataByChr;
        IndData                  *pInd = indData;

        if (!singlePop)
        {
            //Gathered from the global matrix, and the frequencies computed
            //from the same indices -- not from the gathered copy -- so the
            //two cannot disagree about who is in the population.
            subsetDataByIndex(hapDataByChr, GLDataByChr, indData, popIndex[k],
                              &pHap, &pGL, &pInd, USE_GL, PHASED);
            if (fileFreq != NULL)
            {
                //Handed over, not copied: analyzePopulation releases pFreq, so
                //the slot is detached to keep the two from freeing it twice.
                pFreq = fileFreq->at(k);
                fileFreq->at(k) = NULL;
            }
            else
            {
                //Computed and written above; handed over the same way, so the
                //frequencies in <out>.freq.gz are exactly the ones used.
                pFreq = popFreq[k];
                popFreq[k] = NULL;
            }
            pMap = cloneMapData(mapDataByChr);
        }

        //Filter on THIS population's frequencies.  A site monomorphic here is
        //uninformative here, whatever it looks like in another population, and
        //dropping it is what makes a population's windows the ones it would
        //have on its own.  With one population this is the filtering that used
        //to run once before the loop, on the same data, in the same order.
        {
            int before = 0;
            for (unsigned int c = 0; c < pMap->size(); c++) before += pMap->at(c)->nloci;
            int after;
            if (WEIGHTED || CM)
            {
                after = filterMonomorphicAndOOBSites(&pMap, &pHap, &pFreq, &pGL,
                                                     scaffoldMapByChr, USE_GL, PHASED);
                LOG.log("Monomorphic or out of bounds loci filtered:", before - after);
                int numInterpolated = interpolateGeneticmap(&pMap, scaffoldMapByChr);
                LOG.log("Number of genetic map locations interpolated:", numInterpolated);
            }
            else
            {
                after = filterMonomorphicSites(&pMap, &pHap, &pFreq, &pGL, USE_GL, PHASED);
                LOG.log("Monomorphic loci filtered:", before - after);
            }
            LOG.log("Total loci used for analysis:", after);
            numLoci = after;

            //After filtering, and from this population's own map: the density
            //that drives --auto-winsize and --auto-overlap-frac is a property
            //of the sites this population actually uses.
            if ((AUTO_WINSIZE && WEIGHTED) || AUTO_OVERLAP_FRAC)
                variantDensity = calcDensity(numLoci, pMap, centro);

            //Genotype frequencies are derived from the filtered genotypes, so
            //they follow the filtering rather than precede it.
            if (!PHASED && WEIGHTED)
            {
                if (pGF != NULL && pGF != genoFreqDataByChr) releaseGenoFreq(pGF);
                pGF = calculateGenoFreq(pHap);
            }

            //filterSites does not prune in place: it builds new vectors and
            //releases the ones it was given.  With one population those ARE
            //main's, so main's pointers have to follow, or the release at the
            //end of the run would free memory that is already gone.
            if (singlePop)
            {
                mapDataByChr = pMap; hapDataByChr = pHap;
                freqDataByChr = pFreq; GLDataByChr = pGL;
                genoFreqDataByChr = pGF;
            }
        }

        PopResult pop = analyzePopulation(popOpt, params,
                                          singlePop ? string("") : popName,
                                          pHap, pFreq, pMap,
                                          pGL, pGF,
                                          pInd, centro, USE_GL, variantDensity,
                                          &(sexModel.role), &par, &chrLengths);
        //analyzePopulation has released pHap/pFreq/pGL/pGF by now; the map and
        //the per-population IndData are the caller's.
        if (!singlePop) { releaseIndData(pInd); releaseMapData(pMap); }

        if (pop.status == POP_DONE)
        {
            //--winsize-multi finished inside.  It is a diagnostic over window
            //sizes, so it ends the run rather than continuing to the next
            //population.
            if (!singlePop)
            {
                releaseHapData(hapDataByChr); if (freqDataByChr != NULL) releaseFreqData(freqDataByChr);
                if (USE_GL) releaseGLData(GLDataByChr);
                if (genoFreqDataByChr != NULL) releaseGenoFreq(genoFreqDataByChr);
            }
            releaseMapData(mapDataByChr); releaseIndData(indData);
            delete centro; delete params; freeRNG();
            return 0;
        }
        if (pop.status == POP_ERROR) { writeStatus = 2; break; }

        //Machine-readable record of what this run actually did.  The
        //auto-selected values were previously only prose lines in the .log,
        //which is what the documented workflow for reusing a cutoff asks users
        //to parse by hand.  One per population, named like its other outputs.
        {
            vector< pair<string,string> > resolved;
            ostringstream v;
            v << opt.SEED;                    resolved.push_back(make_pair("seed", v.str()));
            if (POOL) resolved.push_back(make_pair("populations_pooled", "true"));
            if (!singlePop)
            {
                //Quoted here: writeParamsJSON emits a resolved value verbatim,
                //so a bare label produced invalid JSON -- caught by --load-params
                //failing to parse its own output.
                v.str(""); v << "\"" << popName << "\"";
                resolved.push_back(make_pair("population", v.str()));
                v.str(""); v << (opt.SEED + (unsigned long int)(k));
                resolved.push_back(make_pair("population_seed", v.str()));
            }
            //From the RESULT, not from opt: these four are selected per
            //population and the option struct still holds whatever the command
            //line said.  Reading opt here would have recorded the starting
            //window size under --auto-winsize, and the sentinel cutoff under
            //the automatic one.
            v.str(""); v << pop.winsize;      resolved.push_back(make_pair("winsize", v.str()));
            v.str(""); v << pop.overlapFrac;  resolved.push_back(make_pair("overlap_frac", v.str()));
            v.str(""); v << pop.lodCutoff;    resolved.push_back(make_pair("lod_cutoff", v.str()));
            //Only when the run had one: on an autosome-only run these would
            //be the autosomal numbers wearing another name, which reads as a
            //sex chromosome having been analysed.
            if (pop.haveSexChr)
            {
                v.str(""); v << pop.sexWinsize;   resolved.push_back(make_pair("sexchr_winsize", v.str()));
                v.str(""); v << pop.sexLodCutoff; resolved.push_back(make_pair("sexchr_lod_cutoff", v.str()));
            }
            //Recorded whether or not the flag was given: the two conventions
            //differ by several percent, so which one produced a column of
            //FROH is not something a reader should have to infer from which
            //flags happen to be in the "set" list.  The source goes with it
            //when there is one -- under analyzed it is what the span check
            //used, under assembly it is the denominator itself.
            v.str(""); v << "\"" << (opt.FROH_DENOM == FROH_ASSEMBLY ? "assembly" : "analyzed") << "\"";
            resolved.push_back(make_pair("froh_denominator", v.str()));
            if (!chrLengths.empty())
            {
                v.str(""); v << "\"" << chrLengths.source() << "\"";
                resolved.push_back(make_pair("chromosome_lengths", v.str()));
            }
            v.str(""); v << KDE_THIN_STEP;    resolved.push_back(make_pair("kde_thin_step", v.str()));
            v.str(""); v << mapDataByChr->size(); resolved.push_back(make_pair("chromosomes_analysed", v.str()));
            v.str(""); v << "[";
            for (unsigned int i = 0; i < populations.size(); i++)
            {
                if (i) v << ", ";
                v << "\"" << populations[i].first << "\"";
            }
            v << "]";                         resolved.push_back(make_pair("populations", v.str()));
            v.str(""); v << "[";
            for (unsigned int i = 0; i < pop.boundSizes.size(); i++) { if (i) v << ", "; v << pop.boundSizes[i]; }
            v << "]";                         resolved.push_back(make_pair("size_bounds", v.str()));
            v.str(""); v << (AUTO_CUTOFF ? "true" : "false");  resolved.push_back(make_pair("cutoff_was_automatic", v.str()));
            v.str(""); v << (AUTO_BOUNDS ? "true" : "false");  resolved.push_back(make_pair("bounds_were_automatic", v.str()));
            try { writeParamsJSON(popOpt.outfile + ".params.json", params, resolved); }
            catch (...) { logCurrentException("writing the parameter record"); writeStatus = 2; }
        }
    }

    //The scaffold is read by every population's filtering and interpolation,
    //so it is released here rather than after the first one.
    if (scaffoldMapByChr != NULL) releaseGenMapScaffold(scaffoldMapByChr);
    for (unsigned int i = 0; i < popFreq.size(); i++)
        if (popFreq[i] != NULL) releaseFreqData(popFreq[i]);
    popFreq.clear();
    if (fileFreq != NULL)
    {
        //Whatever the loop did not consume -- everything, if it stopped early.
        for (unsigned int i = 0; i < fileFreq->size(); i++)
            if (fileFreq->at(i) != NULL) releaseFreqData(fileFreq->at(i));
        fileFreq->clear(); delete fileFreq;
    }

    //The global data outlives the loop only when populations were copied out
    //of it; with one population analyzePopulation was handed the originals.
    if (!singlePop)
    {
        releaseHapData(hapDataByChr);
        if (freqDataByChr != NULL) releaseFreqData(freqDataByChr);
        if (USE_GL) releaseGLData(GLDataByChr);
        if (genoFreqDataByChr != NULL) releaseGenoFreq(genoFreqDataByChr);
    }

    //centro is read by writeFROH; it used to be deleted before the writers ran.
    delete centro;

    releaseIndData(indData);
    //rohLength, rohDataByInd and kdeResult are per-population and released by
    //analyzePopulation, which is what keeps them from accumulating once the
    //loop runs more than once.
    releaseMapData(mapDataByChr);
    delete params;
    freeRNG();
    cout << "Finished.\n";

    //The PTW32_STATIC_LIB process attach/detach calls that used to bracket
    //main are gone: they are pthreads-win32 specific, and std::thread needs no
    //per-process initialisation on any platform.

    return writeStatus;
}
