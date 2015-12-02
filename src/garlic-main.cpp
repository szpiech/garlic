#include <iostream>
#include <cstdio>
#include <fstream>
#include <cmath>
#include <pthread.h>
#include "garlic-cli.h"
#include "garlic-data.h"
#include "garlic-roh.h"
#include "garlic-kde.h"
#include "param_t.h"
#include "gmm.h"
#include "gsl/gsl_statistics.h"
#include "gsl/gsl_sort.h"
#include "BoundFinder.h"
#include "garlic-centromeres.h"

using namespace std;

void exploreWinsizes(vector< HapData * > *hapDataByChr,
                     vector< FreqData * > *freqDataByChr,
                     vector< MapData * > *mapDataByChr,
                     IndData *indData, centromere *centro,
                     vector<int> &multiWinsizes, double error,
                     int MAX_GAP, int KDE_SUBSAMPLE, string outfile);

int selectWinsize(vector< HapData * > *hapDataByChr,
                  vector< FreqData * > *freqDataByChr,
                  vector< MapData * > *mapDataByChr,
                  IndData *indData, centromere *centro,
                  int winsize, double error,
                  int MAX_GAP, int KDE_SUBSAMPLE);

double selectLODCutoff(vector< WinData * > *winDataByChr, int KDE_SUBSAMPLE, string outfile);

int_pair_t selectSizeClasses(ROHLength *rohLength);



int main(int argc, char *argv[])
{
    param_t *params = getCLI(argc, argv);

    bool argerr = false;

    string BUILD = params->getStringFlag(ARG_BUILD);
    argerr = checkBuild(BUILD);

    vector<int> multiWinsizes = params->getIntListFlag(ARG_WINSIZE_MULTI);
    bool WINSIZE_EXPLORE = false;
    argerr = checkMultiWinsizes(multiWinsizes, WINSIZE_EXPLORE);

    bool AUTO_WINSIZE = params->getBoolFlag(ARG_AUTO_WINSIZE);
    argerr = checkAutoWinsize(WINSIZE_EXPLORE, AUTO_WINSIZE);

    string freqfile = params->getStringFlag(ARG_FREQ_FILE);
    bool FREQ_ONLY = params->getBoolFlag(ARG_FREQ_ONLY);
    bool AUTO_FREQ = true;
    argerr = checkAutoFreq(freqfile, FREQ_ONLY, AUTO_FREQ);

    double LOD_CUTOFF = params->getDoubleFlag(ARG_LOD_CUTOFF);
    bool AUTO_CUTOFF = true;
    argerr = checkAutoCutoff(LOD_CUTOFF, AUTO_CUTOFF);

    vector<double> boundSizes = params->getDoubleListFlag(ARG_BOUND_SIZE);
    bool AUTO_BOUNDS = true;
    argerr = checkBoundSizes(boundSizes, AUTO_BOUNDS);

    string tpedfile = params->getStringFlag(ARG_TPED);
    string tfamfile = params->getStringFlag(ARG_TFAM);
    argerr = checkRequiredFiles(tpedfile, tfamfile);

    int numThreads = params->getIntFlag(ARG_THREADS);
    argerr = checkThreads(numThreads);

    double error = params->getDoubleFlag(ARG_ERROR);
    argerr = checkError(error);

    int winsize = params->getIntFlag(ARG_WINSIZE);
    argerr = checkWinsize(winsize);

    int MAX_GAP = params->getIntFlag(ARG_MAX_GAP);
    argerr = checkMaxGap(MAX_GAP);

    if(argerr) return -1;

    int KDE_SUBSAMPLE = params->getIntFlag(ARG_KDE_SUBSAMPLE);
    int nresample = params->getIntFlag(ARG_RESAMPLE);
    string outfile = params->getStringFlag(ARG_OUTFILE);    
    bool RAW_LOD = params->getBoolFlag(ARG_RAW_LOD);
    char TPED_MISSING = params->getCharFlag(ARG_TPED_MISSING);
    string centromereFile = params->getStringFlag(ARG_CENTROMERE_FILE);

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

    try
    {
        chrCoordList = scanTPEDMapData(tpedfile, numLoci, numCols);
        mapDataByChr = readTPEDMapData(tpedfile, numCols, chrCoordList, TPED_MISSING);

        scanIndData3(tfamfile, numInd, popName);
        indData = readIndData3(tfamfile, numInd);

        hapDataByChr = readTPEDHapData3(tpedfile, numLoci, numInd, TPED_MISSING, mapDataByChr);
    }
    catch (...) { return 1; }

    if (AUTO_FREQ)
    {
        cerr << "Calculating allele frequencies...\n";
        freqDataByChr = calcFreqData2(hapDataByChr, nresample);

        string freqOutfile = outfile;
        freqOutfile += ".freq";
        writeFreqData(freqOutfile, popName, freqDataByChr, mapDataByChr, indData);
    }
    else //(!AUTO_FREQ)
    {
        cerr << "Loading user provided allele frequencies from " << freqfile << "...\n";
        try { freqDataByChr = readFreqData(freqfile, popName, chrCoordList, mapDataByChr); }
        catch (...) { return -1; }
    }
    if (FREQ_ONLY) return 0;

    chrCoordList->clear();
    delete chrCoordList;

    //Calcuate LOD scores for a range of window sizes
    //Calculate KDE for these LOD scores
    //Output results to file and exit
    if (WINSIZE_EXPLORE)
    {
        exploreWinsizes(hapDataByChr, freqDataByChr, mapDataByChr,
                        indData, centro, multiWinsizes, error,
                        MAX_GAP, KDE_SUBSAMPLE, outfile);
        return 0;
    }

    if (AUTO_WINSIZE)
    {
        winsize = selectWinsize(hapDataByChr, freqDataByChr, mapDataByChr,
                                indData, centro, winsize, error,
                                MAX_GAP, KDE_SUBSAMPLE);
    }

    cerr << "Window size: " << winsize << endl;

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
        LOD_CUTOFF = selectLODCutoff(winDataByChr, KDE_SUBSAMPLE, makeKDEFilename(outfile, winsize));

        //This can be shifted to a log file
        string lodCutoffOutfile = outfile;
        lodCutoffOutfile += ".lod.cutoff";

        ofstream fout;
        fout.open(lodCutoffOutfile.c_str());
        if (fout.fail())
        {
            cerr << "ERROR: Could not open " << lodCutoffOutfile << " for writing.\n";
            return -1;
        }

        fout << fixed;

        string popName = indData->pop;

        fout << popName << " " << LOD_CUTOFF << "\n";
        cerr << popName << " LOD score cutoff: " << LOD_CUTOFF << "\n";

        cerr << "Wrote " << lodCutoffOutfile << "\n";
        fout.close();
    }
    else
    {
        cerr << "User provided LOD score cutoff: " << LOD_CUTOFF << "\n";
    }

    cerr << "Assembling ROH windows...\n";
    //Assemble ROH for each individual in each pop
    ROHLength *rohLength;
    vector< ROHData * > *rohDataByInd = assembleROHWindows(winDataByChr, mapDataByChr, indData,
                                        centro, LOD_CUTOFF, &rohLength, winsize, MAX_GAP);

    //GMM set up for size classifications
    //There might be some benefit in doing this across a range of ngaussians and choosing the classification
    //That has highest BIC

    int_pair_t bounds;

    if (AUTO_BOUNDS)
    {
        cerr << "Fitting 3-component GMM for size classification.\n";

        bounds = selectSizeClasses(rohLength);

        //This can be shifted to a log file
        string sizeBoundaryOutfile = outfile;
        sizeBoundaryOutfile += ".size.bounds";
        ofstream fout;
        fout.open(sizeBoundaryOutfile.c_str());
        if (fout.fail())
        {
            cerr << "ERROR: Could not open " << sizeBoundaryOutfile << " for writing.\n";
            return -1;
        }

        fout << fixed;
        fout << rohLength->pop << " " << bounds.first << " " << bounds.second << endl;
        cerr << "Wrote " << sizeBoundaryOutfile << endl;
        fout.close();
    }
    else
    {
        bounds.first = boundSizes[0];
        bounds.second = boundSizes[1];

        cerr << "User provided size boundaries.\n";
    }

    cerr << "Population " << rohLength->pop << " A/B: " << bounds.first << " B/C: " << bounds.second << endl;

    //Output ROH calls to file, one for each individual
    //includes A/B/C size classifications
    cerr << "Writing ROH tracts...\n";
    writeROHData(makeROHFilename(outfile), rohDataByInd, mapDataByChr, bounds, popName, VERSION);
    return 0;
}

void exploreWinsizes(vector< HapData * > *hapDataByChr,
                     vector< FreqData * > *freqDataByChr,
                     vector< MapData * > *mapDataByChr,
                     IndData *indData,
                     centromere *centro,
                     vector<int> &multiWinsizes,
                     double error,
                     int MAX_GAP, int KDE_SUBSAMPLE, string outfile)
{
    vector< WinData * > *winDataByChr;
    //This could be paralellized
    for (int i = 0; i < multiWinsizes.size(); i++)
    {
        //Since we are not proceeding with the full ROH calling
        //procedure, we only calculate LOD scores for a random
        //subset of individuals as governed by --kde-subsample
        if (KDE_SUBSAMPLE <= 0)
        {
            winDataByChr = calcLODWindows(hapDataByChr, freqDataByChr, mapDataByChr,
                                          indData, centro, multiWinsizes[i],
                                          error, MAX_GAP);
        }
        else
        {
            vector< HapData * > *subsetHapDataByChr;
            IndData *subsetIndData;
            subsetData(hapDataByChr, indData, &subsetHapDataByChr, &subsetIndData, KDE_SUBSAMPLE);
            winDataByChr = calcLODWindows(subsetHapDataByChr, freqDataByChr, mapDataByChr,
                                          subsetIndData, centro, multiWinsizes[i],
                                          error, MAX_GAP);
            releaseHapData(subsetHapDataByChr);
            releaseIndData(subsetIndData);
        }

        //Format the LOD window data into a single array per pop with no missing data
        //Prepped for KDE
        DoubleData *rawWinData = convertWinData2DoubleData(winDataByChr);

        //Compute KDE of LOD score distribution
        cerr << "Estimating distribution of raw LOD score windows:\n";
        KDEResult *kdeResult = computeKDE(rawWinData->data, rawWinData->size);
        releaseDoubleData(rawWinData);

        //Output kde points
        try { writeKDEResult(kdeResult, makeKDEFilename(outfile, multiWinsizes[i])); }
        catch (...) { throw 0; }

        releaseWinData(winDataByChr);
    }
    return;
}

int selectWinsize(vector< HapData * > *hapDataByChr,
                  vector< FreqData * > *freqDataByChr,
                  vector< MapData * > *mapDataByChr,
                  IndData *indData, centromere *centro,
                  int winsize, double error,
                  int MAX_GAP, int KDE_SUBSAMPLE)
{
    vector< WinData * > *winDataByChr;
    double AUTO_WINSIZE_THRESHOLD = 0.05;
    int winsizeQuery = winsize - 10;

    vector< HapData * > *subsetHapDataByChr;
    IndData *subsetIndData;
    if (KDE_SUBSAMPLE > 0) {
        subsetData(hapDataByChr, indData, &subsetHapDataByChr, &subsetIndData, KDE_SUBSAMPLE);
    }

    cerr << "Searching for acceptable window size:\n\
        winsize\tsummse\n";

    double prev = 1, curr = 1;
    bool finished = false;
    while (!finished)
    {
        //Since we are not proceeding with the full ROH calling
        //procedure, we only calculate LOD scores for a random
        //subset of individuals as governed by --kde-subsample
        if (KDE_SUBSAMPLE <= 0)
        {
            winDataByChr = calcLODWindows(hapDataByChr,
                                          freqDataByChr,
                                          mapDataByChr,
                                          indData,
                                          centro,
                                          winsizeQuery,
                                          error,
                                          MAX_GAP);

        }
        else
        {
            winDataByChr = calcLODWindows(subsetHapDataByChr,
                                          freqDataByChr,
                                          mapDataByChr,
                                          subsetIndData,
                                          centro,
                                          winsizeQuery,
                                          error,
                                          MAX_GAP);
        }

        //Format the LOD window data into a single array per pop with no missing data
        //Prepped for KDE
        DoubleData *rawWinData = convertWinData2DoubleData(winDataByChr);

        //Compute KDE of LOD score distribution
        cerr << "Estimating distribution of raw LOD score windows:\n";
        KDEResult *kdeResult = computeKDE(rawWinData->data, rawWinData->size);
        releaseDoubleData(rawWinData);

        curr = calculateWiggle(kdeResult);
        if (winsizeQuery >= winsize) cerr << winsizeQuery << "\t" << (prev - curr) / prev << "\n";
        if (winsizeQuery >= winsize && (prev - curr) / prev <= AUTO_WINSIZE_THRESHOLD) finished = true;
        else {
            winsizeQuery += 10;
            prev = curr;
        }
        releaseWinData(winDataByChr);
    }

    if (KDE_SUBSAMPLE > 0) {
        releaseHapData(subsetHapDataByChr);
        releaseIndData(subsetIndData);
    }
    return winsizeQuery;
}

double selectLODCutoff(vector< WinData * > *winDataByChr, int KDE_SUBSAMPLE, string kdeoutfile)
{
    //Format the LOD window data into a single array per pop with no missing data
    //Prepped for KDE
    DoubleData *rawWinData;
    double LOD_CUTOFF;

    if (KDE_SUBSAMPLE <= 0) rawWinData = convertWinData2DoubleData(winDataByChr);
    else rawWinData = convertSubsetWinData2DoubleData(winDataByChr, KDE_SUBSAMPLE);

    //Compute KDE of LOD score distribution
    cerr << "Estimating distribution of raw LOD score windows:\n";
    KDEResult *kdeResult = computeKDE(rawWinData->data, rawWinData->size);
    releaseDoubleData(rawWinData);

    //Output kde points
    try { writeKDEResult(kdeResult, kdeoutfile); }
    catch (...) { return -1; }

    try { LOD_CUTOFF = get_min_btw_modes(kdeResult->x, kdeResult->y, 512); }
    catch (...)
    {
        cerr << "ERROR: Failed to find the minimum between modes in the LOD score density.\n";
        cerr << "\tResults from density estimation have been written to file for inspection.\n";
        cerr << "\tA cutoff can be manually specified on the command line with " << ARG_LOD_CUTOFF << ".\n";
        return -1;
    }

    releaseKDEResult(kdeResult);
    return LOD_CUTOFF;
}

int_pair_t selectSizeClasses(ROHLength *rohLength)
{
    int_pair_t bounds;
    size_t *sortIndex;

    int ngaussians = 3;
    size_t maxIter = 1000;
    double tolerance = 1e-8;
    double * W;
    double * Mu;
    double * Sigma;

    W = new double[ngaussians];
    Mu = new double[ngaussians];
    Sigma = new double[ngaussians];
    sortIndex = new size_t[ngaussians];

    //calculate mean and var for the population size distribution to use for initial guess
    double var = gsl_stats_variance(rohLength->length, 1, rohLength->size);
    double mu = gsl_stats_mean(rohLength->length, 1, rohLength->size);
    for (int n = 0; n < ngaussians; n++)
    {
        W[n] = 1.0 / double(ngaussians);
        Mu[n] = mu * double(n + 1) / double(ngaussians + 1);
        Sigma[n] = var * (n + 1) / double(ngaussians);
    }

    GMM gmm(ngaussians, W, Mu, Sigma, maxIter, tolerance, false);

    gmm.estimate(rohLength->length, rohLength->size);

    for (int n = 0; n < ngaussians; n++)
    {
        W[n] = gmm.getMixCoefficient(n);
        Mu[n] = gmm.getMean(n);
        Sigma[n] = gmm.getVar(n);
        sortIndex[n] = n;
    }

    gsl_sort_index(sortIndex, Mu, 1, ngaussians);

    cerr << rohLength->pop << " ["
         << W[sortIndex[0]] << " ," << W[sortIndex[1]] << " ," << W[sortIndex[2]] << "] ["
         << Mu[sortIndex[0]] << " ," << Mu[sortIndex[1]] << " ," << Mu[sortIndex[2]] << "] ["
         << Sigma[sortIndex[0]] << " ," << Sigma[sortIndex[1]] << " ," << Sigma[sortIndex[2]] << "]\n";

    //Find boundaries, there are ngaussians-1 of them, but for the moment this is defined to be 2
    //This finds the 'first' root of the difference between two gaussians
    BoundFinder SM(Mu[sortIndex[0]], Sigma[sortIndex[0]], W[sortIndex[0]], Mu[sortIndex[1]], Sigma[sortIndex[1]], W[sortIndex[1]], 1000, 1e-4, false);
    bounds.first = SM.findBoundary();
    BoundFinder ML(Mu[sortIndex[1]], Sigma[sortIndex[1]], W[sortIndex[1]], Mu[sortIndex[2]], Sigma[sortIndex[2]], W[sortIndex[2]], 1000, 1e-4, false);
    bounds.second = ML.findBoundary();

    delete [] W;
    delete [] Mu;
    delete [] Sigma;

    return bounds;
}