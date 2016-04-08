#include "garlic-errlog.h"
#include "garlic-cli.h"
#include <iostream>
#include <cstdio>
#include <fstream>
#include <cmath>
//#include <pthread.h>
#include "garlic-data.h"
#include "garlic-roh.h"
#include "garlic-kde.h"
#include "param_t.h"
#include "garlic-centromeres.h"

using namespace std;

string getCommandLineString(int argc, char *argv[])
{
    string str = argv[0];
    for (int i = 1; i < argc; i++) {
        str += " " + string(argv[i]);
    }
    return str;
}

int main(int argc, char *argv[])
{
//++++++++++CLI handling++++++++++
    param_t *params = getCLI(argc, argv);
    if (params == NULL) return 0;

    string outfile = params->getStringFlag(ARG_OUTFILE);
    LOG.init(outfile);
    LOG.log(getCommandLineString(argc, argv));
    LOG.log("Output file basename:", outfile);

    bool argerr = false;

    string tpedfile = params->getStringFlag(ARG_TPED);
    string tfamfile = params->getStringFlag(ARG_TFAM);
    argerr = checkRequiredFiles(tpedfile, tfamfile);
    LOG.log("TPED file:", tpedfile);
    char TPED_MISSING = params->getCharFlag(ARG_TPED_MISSING);
    LOG.log("TPED missing data code:", TPED_MISSING);
    LOG.log("TFAM file:", tfamfile);

    string BUILD = params->getStringFlag(ARG_BUILD);
    argerr = checkBuild(BUILD);
    LOG.log("Genome build:", BUILD);

    string centromereFile = params->getStringFlag(ARG_CENTROMERE_FILE);
    LOG.log("User defined centromere file:", centromereFile);

    int nresample = params->getIntFlag(ARG_RESAMPLE);
    string freqfile = params->getStringFlag(ARG_FREQ_FILE);
    bool FREQ_ONLY = params->getBoolFlag(ARG_FREQ_ONLY);
    bool AUTO_FREQ = true;
    argerr = checkAutoFreq(freqfile, FREQ_ONLY, AUTO_FREQ);
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
    argerr = checkMultiWinsizes(multiWinsizes, WINSIZE_EXPLORE);
    LOG.log("Explore window sizes:", WINSIZE_EXPLORE);
    if (WINSIZE_EXPLORE) LOG.logv("User defined window sizes:", multiWinsizes);

    bool AUTO_WINSIZE = params->getBoolFlag(ARG_AUTO_WINSIZE);
    argerr = checkAutoWinsize(WINSIZE_EXPLORE, AUTO_WINSIZE);
    LOG.log("Automatic window size:", AUTO_WINSIZE);

    int winsize = params->getIntFlag(ARG_WINSIZE);
    argerr = checkWinsize(winsize);
    if (!WINSIZE_EXPLORE && !AUTO_WINSIZE) LOG.log("User defined window size:", winsize);

    double LOD_CUTOFF = params->getDoubleFlag(ARG_LOD_CUTOFF);
    bool AUTO_CUTOFF = true;
    argerr = checkAutoCutoff(LOD_CUTOFF, AUTO_CUTOFF);
    LOG.log("Choose LOD score cutoff automatically:", AUTO_CUTOFF);
    if (!AUTO_CUTOFF) LOG.log("User defined LOD score cutoff:", LOD_CUTOFF);

    vector<double> boundSizes = params->getDoubleListFlag(ARG_BOUND_SIZE);
    bool AUTO_BOUNDS = true;
    argerr = checkBoundSizes(boundSizes, AUTO_BOUNDS);
    LOG.log("Choose ROH class thresholds automatically:", AUTO_BOUNDS);
    if (!AUTO_BOUNDS) LOG.logv("User defined ROH class thresholds:", boundSizes);

    int numThreads = params->getIntFlag(ARG_THREADS);
    argerr = checkThreads(numThreads);
    LOG.log("Threads:", numThreads);

    double error = params->getDoubleFlag(ARG_ERROR);
    argerr = checkError(error);
    LOG.log("Genotyping error:", error);

    int MAX_GAP = params->getIntFlag(ARG_MAX_GAP);
    argerr = checkMaxGap(MAX_GAP);
    LOG.log("Max gap:", MAX_GAP);

    if (argerr) return -1;

    int KDE_SUBSAMPLE = params->getIntFlag(ARG_KDE_SUBSAMPLE);
    if (KDE_SUBSAMPLE <= 0) LOG.log("# of rand individuals for KDE: ALL");
    else LOG.log("# of rand individuals for KDE:", KDE_SUBSAMPLE);

    bool RAW_LOD = params->getBoolFlag(ARG_RAW_LOD);
    LOG.log("Output raw LOD scores:", RAW_LOD);

    //double AUTO_WINSIZE_THRESHOLD = 0.5;

//++++++++++Datafile reading++++++++++
    centromere *centro;
    centro = new centromere(BUILD, centromereFile, DEFAULT_CENTROMERE_FILE);

    int numLoci, numInd, numCols;
    vector< int_pair_t > *chrCoordList;
    vector< MapData * > *mapDataByChr;
    string popName;
    IndData *indData;
    vector< HapData * > *hapDataByChr;
    vector< FreqData * > *freqDataByChr;
    vector< WinData * > *winDataByChr;
    KDEResult *kdeResult;
    try
    {
        chrCoordList = scanTPEDMapData(tpedfile, numLoci, numCols);
        mapDataByChr = readTPEDMapData(tpedfile, numCols, chrCoordList, TPED_MISSING);

        LOG.log("Total loci:", numLoci);

        scanIndData3(tfamfile, numInd, popName);
        indData = readIndData3(tfamfile, numInd);

        LOG.log("Population:", popName);
        LOG.log("Total diploid individuals:", numInd);

        hapDataByChr = readTPEDHapData3(tpedfile, numLoci, numInd, TPED_MISSING, mapDataByChr);
    }
    catch (...) { return 1; }

//++++++++++Allele frequencies++++++++++
    if (AUTO_FREQ)
    {
        cout << "Calculating allele frequencies\n";
        freqDataByChr = calcFreqData2(hapDataByChr, nresample);

        string freqOutfile = outfile;
        freqOutfile += ".freq";
        writeFreqData(freqOutfile, popName, freqDataByChr, mapDataByChr, indData);
    }
    else //(!AUTO_FREQ)
    {
        cout << "Loading user provided allele frequencies from " << freqfile << "\n";
        try { freqDataByChr = readFreqData(freqfile, popName, chrCoordList, mapDataByChr); }
        catch (...) { return -1; }
    }
    if (FREQ_ONLY) return 0;

    chrCoordList->clear();
    delete chrCoordList;


//Filter data based on frequency data.
//Remove all monomorphic sites.
//If a frequency file is provided that reports
//a frequency in (0,1) the site will be retained
//even if it appears monomorphic in the sample.
    int newLoci = filterMonomorphicSites(&mapDataByChr, &hapDataByChr, &freqDataByChr);

    LOG.log("Monomorphic loci filtered:", numLoci - newLoci);
    LOG.log("Total loci used for analysis:", newLoci);

    numLoci = newLoci;


//++++++++++Pipeline begins++++++++++
    if (WINSIZE_EXPLORE && AUTO_WINSIZE)
    {
        kdeResult = selectWinsizeFromList(hapDataByChr, freqDataByChr, mapDataByChr,
                                          indData, centro, &multiWinsizes, winsize, error,
                                          MAX_GAP, KDE_SUBSAMPLE, outfile);
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
                        MAX_GAP, KDE_SUBSAMPLE, outfile);

        return 0;
    }
    else if (AUTO_WINSIZE)
    {
        try
        {
            kdeResult = selectWinsize(hapDataByChr, freqDataByChr, mapDataByChr,
                                      indData, centro, winsize, error,
                                      MAX_GAP, KDE_SUBSAMPLE, outfile);
        }
        catch (...)
        {
            return 1;
        }
        /*
        kdeResult = automaticallyChooseWindowSize(hapDataByChr, freqDataByChr, mapDataByChr,
                    indData, centro, winsize, error, MAX_GAP, KDE_SUBSAMPLE, numThreads,
                    WINSIZE_EXPLORE, AUTO_WINSIZE_THRESHOLD, outfile);
                    */
        LOG.log("Selected window size:", winsize);
    }

    cout << "Window size: " << winsize << endl;

    winDataByChr = calcLODWindows(hapDataByChr, freqDataByChr, mapDataByChr,
                                  indData, centro, winsize, error, MAX_GAP);

    releaseHapData(hapDataByChr);
    releaseFreqData(freqDataByChr);

    if (RAW_LOD)
    {
        //Output raw windows
        try { writeWinData(winDataByChr, indData, mapDataByChr, outfile); }
        catch (...) { return -1; }
    }

    if (AUTO_CUTOFF)
    {
        if (!AUTO_WINSIZE && !WINSIZE_EXPLORE)
        {
            LOD_CUTOFF = selectLODCutoff(winDataByChr, indData, KDE_SUBSAMPLE, makeKDEFilename(outfile, winsize));
        }
        else
        {
            LOD_CUTOFF = selectLODCutoff(kdeResult);
        }
        LOG.log("Selected LOD score cutoff:", LOD_CUTOFF);
        //cout << "Selected LOD score cutoff: " << LOD_CUTOFF << endl;
    }
    else
    {
        cout << "User defined LOD score cutoff: " << LOD_CUTOFF << "\n";
    }

    cout << "Assembling ROH windows\n";
    //Assemble ROH for each individual in each pop
    ROHLength *rohLength;
    vector< ROHData * > *rohDataByInd = assembleROHWindows(winDataByChr, mapDataByChr, indData,
                                        centro, LOD_CUTOFF, &rohLength, winsize, MAX_GAP);

    releaseWinData(winDataByChr);
    delete centro;
    int_pair_t bounds;

    if (AUTO_BOUNDS)
    {
        cout << "Fitting 3-component GMM for size classification\n";
        bounds = selectSizeClasses(rohLength);
        LOG.log("Selected ROH size boundaries ( A/B, B/C ) = (", bounds.first, false);
        LOG.log(",", bounds.second, false);
        LOG.log(" )");
    }
    else
    {
        bounds.first = boundSizes[0];
        bounds.second = boundSizes[1];
        cout << "User provided size boundaries.\n";
    }

    cout << "ROH size boundaries ( A/B, B/C ) = ( " << bounds.first << ", " << bounds.second << " )\n";

    //Output ROH calls to file, one for each individual
    //includes A/B/C size classifications
    cout << "Writing ROH tracts.\n";
    writeROHData(makeROHFilename(outfile), rohDataByInd, mapDataByChr, bounds, popName, VERSION);

    //releaseIndData(indData);
    //releaseROHLength(rohLength);
    //releaseROHData(rohDataByInd);
    //releaseKDEResult(kdeResult);
    //releaseMapData(mapDataByChr);
    //delete params;
    cout << "Finished.\n";
    return 0;
}
