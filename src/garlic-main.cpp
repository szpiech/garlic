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
};

static PopResult analyzePopulation(const GarlicOptions &opt,
                                   param_t *params,
                                   vector< HapData * > *hapDataByChr,
                                   vector< FreqData * > *freqDataByChr,
                                   vector< MapData * > *mapDataByChr,
                                   vector< GenoLikeData * > *GLDataByChr,
                                   vector< GenoFreqData * > *genoFreqDataByChr,
                                   IndData *indData,
                                   centromere *centro,
                                   bool USE_GL,
                                   double variantDensity)
{
    PopResult res;
    res.status = POP_OK;

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
                                          MAX_WINSIZE);
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
    if (freqDataByChr != NULL) releaseFreqData(freqDataByChr);
    if (USE_GL) releaseGLData(GLDataByChr);

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
            LOD_CUTOFF = selectLODCutoff(winDataByChr, indData, KDE_SUBSAMPLE, makeKDEFilename(outfile, winsize), (KDE_THIN_STEP > 0 ? KDE_THIN_STEP : winsize), winsize, cutoffOK);
        }
        else LOD_CUTOFF = selectLODCutoff(kdeResult, winsize, cutoffOK);

        //A failed KDE used to return -1, which main then used as the cutoff and
        //happily wrote a complete, plausible-looking .roh.bed from it.
        if (!cutoffOK)
        {
            LOG.err("ERROR: Could not select a LOD score cutoff automatically. Stopping.");
            res.status = POP_ERROR; return res;
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
                      indData->pop, centro, CM);
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
    //--pool-populations restores the pre-per-population behaviour: one
    //analysis over everyone, pooled frequencies, unlabelled output names.
    //Implemented by making the loop take the single-population path rather
    //than by rewriting the population list, so the run record still reports
    //which populations were actually present.
    bool POOL = params->getBoolFlag(ARG_POOL);
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
            else pFreq = calcFreqDataForIndices(hapDataByChr, popIndex[k], nresample);
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
                                          pHap, pFreq, pMap,
                                          pGL, pGF,
                                          pInd, centro, USE_GL, variantDensity);
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
