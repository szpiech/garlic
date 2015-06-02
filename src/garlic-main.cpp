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

using namespace std;

int main(int argc, char *argv[])
{
    param_t params;

    params.addFlag(ARG_OUTFILE, DEFAULT_OUTFILE, "", HELP_OUTFILE);
    params.addFlag(ARG_THREADS, DEFAULT_THREADS, "", HELP_THREADS);
    params.addFlag(ARG_ERROR, DEFAULT_ERROR, "", HELP_ERROR);
    params.addFlag(ARG_WINSIZE, DEFAULT_WINSIZE, "", HELP_WINSIZE);
    params.addFlag(ARG_POINTS, DEFAULT_POINTS, "", HELP_POINTS);
    params.addFlag(ARG_BW, DEFAULT_BW, "", HELP_BW);
    params.addFlag(ARG_MAX_GAP, DEFAULT_MAX_GAP, "", HELP_MAX_GAP);
    params.addFlag(ARG_RESAMPLE, DEFAULT_RESAMPLE, "", HELP_RESAMPLE);
    params.addFlag(ARG_TPED, DEFAULT_TPED, "", HELP_TPED);
    params.addFlag(ARG_TFAM, DEFAULT_TFAM, "", HELP_TFAM);
    params.addFlag(ARG_RAW_LOD, DEFAULT_RAW_LOD, "", HELP_RAW_LOD);
    params.addListFlag(ARG_BOUND_SIZE, DEFAULT_BOUND_SIZE, "", HELP_BOUND_SIZE);
    params.addFlag(ARG_LOD_CUTOFF, DEFAULT_LOD_CUTOFF, "", HELP_LOD_CUTOFF);
    params.addFlag(ARG_LOD_CUTOFF_FILE, DEFAULT_LOD_CUTOFF_FILE, "", HELP_LOD_CUTOFF_FILE);
    params.addFlag(ARG_BOUND_SIZE_FILE, DEFAULT_BOUND_SIZE_FILE, "", HELP_BOUND_SIZE_FILE);
    params.addFlag(ARG_TPED_MISSING, DEFAULT_TPED_MISSING, "", HELP_TPED_MISSING);
    params.addFlag(ARG_FREQ_FILE, DEFAULT_FREQ_FILE, "", HELP_FREQ_FILE);
    params.addFlag(ARG_FREQ_ONLY, DEFAULT_FREQ_ONLY, "", HELP_FREQ_ONLY);
    params.addListFlag(ARG_WINSIZE_MULTI, DEFAULT_WINSIZE_MULTI, "", HELP_WINSIZE_MULTI);
    params.addFlag(ARG_POP_SPLIT, DEFAULT_POP_SPLIT , "", HELP_POP_SPLIT);
    params.addFlag(ARG_KDE_SUBSAMPLE, DEFAULT_KDE_SUBSAMPLE , "", HELP_KDE_SUBSAMPLE);
    params.addFlag(ARG_AUTO_WINSIZE, DEFAULT_AUTO_WINSIZE, "", HELP_AUTO_WINSIZE);

    try
    {
        params.parseCommandLine(argc, argv);
    }
    catch (...)
    {
        return -1;
    }

    int KDE_SUBSAMPLE = params.getIntFlag(ARG_KDE_SUBSAMPLE);
    int POP_SPLIT = params.getBoolFlag(ARG_POP_SPLIT);
    int winsize = params.getIntFlag(ARG_WINSIZE);
    int* winsizeByPop = NULL;
    string outfile = params.getStringFlag(ARG_OUTFILE);
    string tpedfile = params.getStringFlag(ARG_TPED);
    string tfamfile = params.getStringFlag(ARG_TFAM);
    int numThreads = params.getIntFlag(ARG_THREADS);
    double error = params.getDoubleFlag(ARG_ERROR);
    int MAX_GAP = params.getIntFlag(ARG_MAX_GAP);
    int nresample = params.getIntFlag(ARG_RESAMPLE);
    bool TPED = false;
    bool RAW_LOD = params.getBoolFlag(ARG_RAW_LOD);
    vector<double> boundSizes = params.getDoubleListFlag(ARG_BOUND_SIZE);
    double LOD_CUTOFF = params.getDoubleFlag(ARG_LOD_CUTOFF);
    string lodCutoffFile = params.getStringFlag(ARG_LOD_CUTOFF_FILE);
    string boundSizeFile = params.getStringFlag(ARG_BOUND_SIZE_FILE);
    bool AUTO_BOUNDS = true;
    bool AUTO_CUTOFF = true;
    char TPED_MISSING = params.getCharFlag(ARG_TPED_MISSING);
    string freqfile = params.getStringFlag(ARG_FREQ_FILE);
    bool AUTO_FREQ = true;
    bool FREQ_ONLY = params.getBoolFlag(ARG_FREQ_ONLY);
    vector<int> multiWinsizes = params.getIntListFlag(ARG_WINSIZE_MULTI);
    bool WINSIZE_EXPLORE = false;
    bool AUTO_WINSIZE = params.getBoolFlag(ARG_AUTO_WINSIZE);

    double AUTO_WINSIZE_THRESHOLD = 0.05;

    if (multiWinsizes[0] != DEFAULT_WINSIZE_MULTI)
    {
        for (int i = 0; i < multiWinsizes.size(); i++)
        {
            if (multiWinsizes[i] <= 0)
            {
                cerr << "ERROR: SNP window sizes must be > 1.\n";
                return -1;
            }
        }
        WINSIZE_EXPLORE = true;
    }

    if (freqfile.compare(DEFAULT_FREQ_FILE) != 0)
    {
        AUTO_FREQ = false;
        if (FREQ_ONLY)
        {
            cerr << "ERROR: Specifying both " << ARG_FREQ_ONLY << " and " << ARG_FREQ_FILE << " accomplishes nothing.\n";
            return -1;
        }
    }

    //Check if both AUTO_WINSIZE and WINSIZE_EXPLORE are set
    //If so, exit with error.
    if (WINSIZE_EXPLORE && AUTO_WINSIZE)
    {
        cerr << "ERROR: Must set only one of " << ARG_WINSIZE_MULTI << " and " << ARG_AUTO_WINSIZE << ".\n";
        return -1;
    }

    //Check if both LOD_CUTOFF and LOD_CUTOFF_FILE defined
    //and error if so
    if (lodCutoffFile.compare(DEFAULT_LOD_CUTOFF_FILE) != 0 && LOD_CUTOFF != DEFAULT_LOD_CUTOFF)
    {
        cerr << "ERROR: At most, only one of " << ARG_LOD_CUTOFF << " and " << ARG_LOD_CUTOFF_FILE << " should be specified.\n";
        return -1;
    }
    else if (LOD_CUTOFF != DEFAULT_LOD_CUTOFF || lodCutoffFile.compare(DEFAULT_LOD_CUTOFF_FILE) != 0)
    {
        AUTO_CUTOFF = false;
    }

    if (boundSizes[0] != DEFAULT_BOUND_SIZE && boundSizeFile.compare(DEFAULT_BOUND_SIZE_FILE) != 0)
    {
        cerr << "ERROR: At most, only one of " << ARG_BOUND_SIZE << " and " << ARG_BOUND_SIZE_FILE << " should be specified.\n";
        return -1;
    }
    else if (boundSizes[0] != DEFAULT_BOUND_SIZE || boundSizeFile.compare(DEFAULT_BOUND_SIZE_FILE) != 0)
    {
        AUTO_BOUNDS = false;
    }

    if (boundSizes.size() == 2)
    {
        double tmp;
        AUTO_BOUNDS = false;
        if (boundSizes[0] <= 0 || boundSizes[1] <= 0)
        {
            cerr << "ERROR: User provided size boundaries must be positive.\n";
            return -1;
        }
        else if (boundSizes[0] > boundSizes[1])
        {
            tmp = boundSizes[0];
            boundSizes[0] = boundSizes[1];
            boundSizes[1] = tmp;
        }
        else if (boundSizes[0] == boundSizes[1])
        {
            cerr << "ERROR: Size boundaries must be different.\n";
            return -1;
        }
    }
    else if (boundSizes.size() > 2)
    {
        cerr << "ERROR: Must provide exactly two boundaries.\n";
        return -1;
    }


    if (tpedfile.compare(DEFAULT_TPED) != 0 && tfamfile.compare(DEFAULT_TFAM) != 0) TPED = true;
    else
    {
        cerr << "ERROR: Must provide both a tped and tfam file.\n";
        return -1;
    }

    /*
    if (TPED == false)
    {
        if (mapfile.compare(DEFAULT_MAPFILE) == 0 ||
                hapfile.compare(DEFAULT_HAPFILE) == 0 ||
                indfile.compare(DEFAULT_INDFILE) == 0)
        {
            cerr << "ERROR: Must specify map/hap/ind files or tped/tfam files.\n";
            return 1;
        }
        if (tpedfile.compare(DEFAULT_TPED) != 0 || tfamfile.compare(DEFAULT_TFAM) != 0)
        {
            cerr << "ERROR: Must specify map/hap/ind files or tped/tfam files, but not both.\n";
            return 1;
        }
    }
    else
    {
        if (mapfile.compare(DEFAULT_MAPFILE) != 0 ||
                hapfile.compare(DEFAULT_HAPFILE) != 0 ||
                indfile.compare(DEFAULT_INDFILE) != 0)
        {
            cerr << "ERROR: Must specify map/hap/ind files or tped/tfam files, but not both.\n";
            return 1;
        }
    }
    */
    if (numThreads <= 0)
    {
        cerr << "ERROR: Number of threads must be > 0.\n";
        return 1;
    }

    if (error <= 0 || error >= 1)
    {
        cerr << "ERROR: Genotype error rate must be > 0 and < 1.\n";
        return 1;
    }

    if (winsize <= 1)
    {
        cerr << "ERROR: SNP window size must be > 1.\n";
        return 1;
    }

    int numLoci, numInd;
    vector< int_pair_t > *chrCoordList;
    vector< MapData * > *mapDataByChr;

    //vector< int_pair_t > *indCoordList;

    string *indList;
    map<string, string> ind2pop;
    map<string, int> pop2size;
    map<string, int> pop2index;
    vector< IndData * > *indDataByPop;

    vector< vector< HapData * >* > *hapDataByPopByChr;
    vector< vector< FreqData * >* > *freqDataByPopByChr;

    vector< vector< WinData * >* > *winDataByPopByChr;

    map<string, double> pop2lodcutoff;
    map<string, double> pop2SMbound;
    map<string, double> pop2MLbound;
    string *oneAllele;
    try
    {
        if (TPED)
        {
            chrCoordList = scanTPEDMapData(tpedfile, numLoci);
            mapDataByChr = readTPEDMapData(tpedfile, chrCoordList, TPED_MISSING);

            scanIndData2(tfamfile, numInd, ind2pop, pop2size);
            indList = new string[numInd];
            indDataByPop = readIndData2(tfamfile, numInd, ind2pop, pop2size, indList, pop2index);

            //Read user pprovided freq data here
            if (!AUTO_FREQ)
            {
                cerr << "Loading user provided allele frequencies from " << freqfile << "...\n";
                freqDataByPopByChr = readFreqData(freqfile, chrCoordList, mapDataByChr, pop2index);
            }

            hapDataByPopByChr = readTPEDHapData2(tpedfile, numLoci, numInd, chrCoordList, indList,
                                                 ind2pop, pop2size, pop2index, TPED_MISSING, mapDataByChr);
        }
        //load LOD score cutoff if specified
        if (LOD_CUTOFF != DEFAULT_LOD_CUTOFF)
        {
            for (map<string, int>::iterator it = pop2size.begin(); it != pop2size.end(); it++)
            {
                pop2lodcutoff[it->first] = LOD_CUTOFF;
            }
        }
        else if (lodCutoffFile.compare(DEFAULT_LOD_CUTOFF_FILE) != 0)
        {
            pop2lodcutoff = readLODCutoff(lodCutoffFile, pop2size);
        }

        //load size bounds if specified
        if (boundSizes[0] != DEFAULT_BOUND_SIZE)
        {
            if (boundSizes.size() == 2)
            {
                double tmp;
                if (boundSizes[0] <= 0 || boundSizes[1] <= 0)
                {
                    cerr << "ERROR: User provided size boundaries must be positive.\n";
                    return -1;
                }
                else if (boundSizes[0] > boundSizes[1])
                {
                    tmp = boundSizes[0];
                    boundSizes[0] = boundSizes[1];
                    boundSizes[1] = tmp;
                }
                else if (boundSizes[0] == boundSizes[1])
                {
                    cerr << "ERROR: Size boundaries must be different.\n";
                    return -1;
                }

                for (map<string, int>::iterator it = pop2size.begin(); it != pop2size.end(); it++)
                {
                    pop2SMbound[it->first] = boundSizes[0];
                    pop2MLbound[it->first] = boundSizes[1];
                }
            }
            else if (boundSizes.size() > 2)
            {
                cerr << "ERROR: Must provide exactly two boundaries.\n";
                return -1;
            }
        }
        else if (boundSizeFile.compare(DEFAULT_BOUND_SIZE_FILE) != 0)
        {
            readBoundSizes(boundSizeFile, pop2SMbound, pop2MLbound, pop2size);
        }

        if (AUTO_FREQ)
        {
            cerr << "Calculating allele frequencies...\n";
            freqDataByPopByChr = calcFreqData(hapDataByPopByChr, nresample);

            string freqOutfile = outfile;
            freqOutfile += ".freq";
            writeFreqData(freqOutfile, freqDataByPopByChr, mapDataByChr, indDataByPop);
        }


        if (POP_SPLIT)
        {
            writeTFAMDataByPop(outfile, indDataByPop, pop2index);
            writeTPEDDataByPop(outfile, hapDataByPopByChr, mapDataByChr, pop2index);
        }
        if (FREQ_ONLY || POP_SPLIT) return 0;
    }
    catch (...)
    {
        return 1;
    }

    int numPop = pop2size.size();
    chrCoordList->clear();
    delete chrCoordList;

    //Initialize winsizeByPop with the value given by --winsize
    //If --auto-winsize is set, will be dealt with below
    winsizeByPop = new int [numPop];
    for (int i = 0; i < numPop; i++)
    {
        winsizeByPop[i] = winsize;
    }

    //Calcuate LOD scores for a range of window sizes
    //Calculate KDE for these LOD scores
    //Output results to file and exit
    if (WINSIZE_EXPLORE)
    {
        for (int i = 0; i < multiWinsizes.size(); i++)
        {
            char winStr[10];
            sprintf(winStr, "%d", multiWinsizes[i]);
            string kdeoutfile = outfile;
            kdeoutfile += ".";
            kdeoutfile += winStr;

            for (int j = 0; j < numPop; j++)
            {
                winsizeByPop[j] = multiWinsizes[i];
            }

            //Since we are not proceeding with the full ROH calling
            //procedure, we only calculate LOD scores for a random
            //subset of individuals as governed by --kde-subsample
            if (KDE_SUBSAMPLE <= 0)
            {
                winDataByPopByChr = calcLODWindows(hapDataByPopByChr,
                                                   freqDataByPopByChr,
                                                   mapDataByChr,
                                                   indDataByPop,
                                                   winsizeByPop,
                                                   error,
                                                   MAX_GAP,
                                                   numThreads);
            }
            else
            {
                vector< vector< HapData * >* > *subsetHapDataByPopByChr;
                vector< IndData * > *subsetIndDataByPop;
                subsetData(hapDataByPopByChr, indDataByPop, &subsetHapDataByPopByChr, &subsetIndDataByPop, KDE_SUBSAMPLE);
                winDataByPopByChr = calcLODWindows(subsetHapDataByPopByChr,
                                                   freqDataByPopByChr,
                                                   mapDataByChr,
                                                   subsetIndDataByPop,
                                                   winsizeByPop,
                                                   error,
                                                   MAX_GAP,
                                                   numThreads);
                releaseHapData(subsetHapDataByPopByChr);
                releaseIndData(subsetIndDataByPop);
            }


            //Format the LOD window data into a single array per pop with no missing data
            //Prepped for KDE
            vector < DoubleData * > *rawWinDataByPop;
            rawWinDataByPop = convertWinData2DoubleData(winDataByPopByChr);

            //Compute KDE of LOD score distribution
            cerr << "Estimating distribution of raw LOD score windows:\n";
            vector < KDEResult * > *kdeResultByPop = computeKDE(rawWinDataByPop, indDataByPop, numThreads);
            releaseDoubleData(rawWinDataByPop);

            //Output kde points
            try
            {
                writeKDEResult(kdeResultByPop, indDataByPop, kdeoutfile);
            }
            catch (...)
            {
                return -1;
            }

            releaseWinData(winDataByPopByChr);
        }

        return 0;
    }
    else if (AUTO_WINSIZE)
    {
        for (int i = 0; i < numPop; i++)
        {
            winsizeByPop[i] -= 10;
        }
        vector< vector< HapData * >* > *subsetHapDataByPopByChr;
        vector< IndData * > *subsetIndDataByPop;
        if (KDE_SUBSAMPLE > 0) {
            subsetData(hapDataByPopByChr, indDataByPop, &subsetHapDataByPopByChr, &subsetIndDataByPop, KDE_SUBSAMPLE);
        }
        /*
        for (int pop = 0; pop < numPop; pop++) {
            cerr << subsetIndDataByPop->at(pop)->pop << "\n";
            for (int chr = 0; chr < hapDataByPopByChr->at(pop)->size(); chr++) {
                cerr << mapDataByChr->at(chr)->chr << endl;
                for (int locus = 0; locus < 3q0; locus++) {
                    cerr << locus << " ";
                    for (int ind = 0; ind < subsetIndDataByPop->at(pop)->nind; ind++) {
                        cerr << subsetHapDataByPopByChr->at(pop)->at(chr)->data[locus][ind] << " ";
                    }
                    cerr << endl;
                }

            }
        }
        */
        cerr << "Searching for acceptable window size:\n\
        pop\twinsize\tsummse\n";
        for (int pop = 0; pop < numPop; pop++)
        {
            double prev = 1, curr = 1;
            bool finished = false;
            while (!finished)
            {
                //Since we are not proceeding with the full ROH calling
                //procedure, we only calculate LOD scores for a random
                //subset of individuals as governed by --kde-subsample
                if (KDE_SUBSAMPLE <= 0)
                {
                    winDataByPopByChr = calcLODWindowsSinglePop(hapDataByPopByChr,
                                        freqDataByPopByChr,
                                        mapDataByChr,
                                        indDataByPop,
                                        winsizeByPop,
                                        error,
                                        MAX_GAP,
                                        numThreads,
                                        pop);

                }
                else
                {
                    winDataByPopByChr = calcLODWindowsSinglePop(subsetHapDataByPopByChr,
                                        freqDataByPopByChr,
                                        mapDataByChr,
                                        subsetIndDataByPop,
                                        winsizeByPop,
                                        error,
                                        MAX_GAP,
                                        numThreads,
                                        pop);
                }

                //Format the LOD window data into a single array per pop with no missing data
                //Prepped for KDE
                vector < DoubleData * > *rawWinDataByPop;
                rawWinDataByPop = convertWinData2DoubleData(winDataByPopByChr);

                //Compute KDE of LOD score distribution
                cerr << "Estimating distribution of raw LOD score windows:\n";
                KDEResult *kdeResult = computeKDE(rawWinDataByPop->at(0)->data, rawWinDataByPop->at(0)->size);
                releaseDoubleData(rawWinDataByPop);

                curr = calculateWiggle(kdeResult);
                if(winsizeByPop[pop] >= winsize) cerr << indDataByPop->at(pop)->pop << "\t" << winsizeByPop[pop] << "\t" << (prev - curr) / prev << "\n";
                if (winsizeByPop[pop] >= winsize && (prev - curr) / prev <= AUTO_WINSIZE_THRESHOLD) finished = true;
                else {
                    winsizeByPop[pop] += 10;
                    prev = curr;
                }
                releaseWinData(winDataByPopByChr);
            }
        }
        if (KDE_SUBSAMPLE > 0) {
            releaseHapData(subsetHapDataByPopByChr);
            releaseIndData(subsetIndDataByPop);
        }
    }
    cerr << "\nFinal window size:\n\
    pop\twinsize\n";
    for (int i = 0; i < numPop; i++) {
        cerr << indDataByPop->at(i)->pop << "\t" << winsizeByPop[i] << endl;
    }

    return 0;

    winDataByPopByChr = calcLODWindows(hapDataByPopByChr,
                                       freqDataByPopByChr,
                                       mapDataByChr,
                                       indDataByPop,
                                       winsizeByPop,
                                       error,
                                       MAX_GAP,
                                       numThreads);


    releaseHapData(hapDataByPopByChr);
    releaseFreqData(freqDataByPopByChr);


    if (RAW_LOD)
    {
        //Output raw windows
        try
        {
            writeWinData(winDataByPopByChr, indDataByPop, mapDataByChr, outfile);
        }
        catch (...)
        {
            return -1;
        }
    }

    double *lodScoreCutoffByPop = new double[numPop];
    if (AUTO_CUTOFF)
    {
        char winStr[10];
        sprintf(winStr, "%d", winsize);
        string kdeoutfile = outfile;
        kdeoutfile += ".";
        kdeoutfile += winStr;
        //Format the LOD window data into a single array per pop with no missing data
        //Prepped for KDE
        vector < DoubleData * > *rawWinDataByPop;

        if (KDE_SUBSAMPLE <= 0)
        {
            rawWinDataByPop = convertWinData2DoubleData(winDataByPopByChr);
        }
        else
        {
            rawWinDataByPop = convertSubsetWinData2DoubleData(winDataByPopByChr, KDE_SUBSAMPLE);
        }

        //Compute KDE of LOD score distribution
        cerr << "Estimating distribution of raw LOD score windows:\n";
        vector < KDEResult * > *kdeResultByPop = computeKDE(rawWinDataByPop, indDataByPop, numThreads);
        releaseDoubleData(rawWinDataByPop);

        //Output kde points
        try
        {
            writeKDEResult(kdeResultByPop, indDataByPop, kdeoutfile);
        }
        catch (...)
        {
            return -1;
        }

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

        //For each population, find the LOD score cutoff
        for (int pop = 0; pop < numPop; pop++)
        {
            try
            {
                lodScoreCutoffByPop[pop] = get_min_btw_modes(kdeResultByPop->at(pop)->x, kdeResultByPop->at(pop)->y, 512);
            }
            catch (...)
            {
                cerr << "ERROR: Failed to find the minimum between modes in the LOD score density.\n";
                cerr << "\tResults from density estimation have been written to file for inspection.\n";
                cerr << "\tA cutoff can be manually specified on the command line with " << ARG_LOD_CUTOFF << endl;
                cerr << "\tor " << ARG_LOD_CUTOFF_FILE << endl;
                return -1;
            }

            string popName = indDataByPop->at(pop)->pop;

            fout << popName << " " << lodScoreCutoffByPop[pop] << "\n";
            cerr << popName << " LOD score cutoff: " << lodScoreCutoffByPop[pop] << "\n";
        }
        cerr << "Wrote " << lodCutoffOutfile << "\n";
        fout.close();


        releaseKDEResult(kdeResultByPop);
    }
    else
    {
        for (map<string, double>::iterator it = pop2lodcutoff.begin(); it != pop2lodcutoff.end(); it++)
        {
            lodScoreCutoffByPop[pop2index[it->first]] = it->second;
            cerr << "User provided LOD score cutoff: " << it->first << " " << lodScoreCutoffByPop[pop2index[it->first]] << "\n";
        }
    }

    cerr << "Assembling ROH windows...\n";
    //Assemble ROH for each individual in each pop
    vector< ROHLength * > *rohLengthByPop = new vector< ROHLength * >;
    vector< vector< ROHData * >* > *rohDataByPopByInd = assembleROHWindows(winDataByPopByChr,
            mapDataByChr,
            indDataByPop,
            lodScoreCutoffByPop,
            &rohLengthByPop,
            winsize);

    //GMM set up for size classifications
    //There might be some benefit in doing this across a range of ngaussians and choosing the classification
    //That has highest BIC
    int ngaussians = 3;
    size_t maxIter = 1000;
    double tolerance = 1e-8;
    double *W;
    double *Mu;
    double *Sigma;
    double *shortMedBound = new double[rohLengthByPop->size()];
    double *medLongBound = new double[rohLengthByPop->size()];
    size_t *sortIndex;

    if (AUTO_BOUNDS)
    {
        cerr << "Fitting 3-component GMM for size classification.\n";
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

        W = new double[ngaussians];
        Mu = new double[ngaussians];
        Sigma = new double[ngaussians];
        sortIndex = new size_t[ngaussians];

        for (int pop = 0; pop < rohLengthByPop->size(); pop++)
        {
            //calculate mean and var for the population size distribution to use for initial guess
            double var = gsl_stats_variance(rohLengthByPop->at(pop)->length, 1, rohLengthByPop->at(pop)->size);
            double mu = gsl_stats_mean(rohLengthByPop->at(pop)->length, 1, rohLengthByPop->at(pop)->size);
            for (int n = 0; n < ngaussians; n++)
            {
                W[n] = 1.0 / double(ngaussians);
                Mu[n] = mu * double(n + 1) / double(ngaussians + 1);
                //gsl_stats_quantile_from_sorted_data(rohLengthByPop->at(pop)->length, 1, rohLengthByPop->at(pop)->size, double(n+1)/double(ngaussians+1));
                Sigma[n] = var * (n + 1) / double(ngaussians);
            }

            GMM gmm(ngaussians, W, Mu, Sigma, maxIter, tolerance, false);

            gmm.estimate(rohLengthByPop->at(pop)->length, rohLengthByPop->at(pop)->size);

            for (int n = 0; n < ngaussians; n++)
            {
                W[n] = gmm.getMixCoefficient(n);
                Mu[n] = gmm.getMean(n);
                Sigma[n] = gmm.getVar(n);
                sortIndex[n] = n;
            }

            gsl_sort_index(sortIndex, Mu, 1, ngaussians);

            cerr << rohLengthByPop->at(pop)->pop << " ["
                 << W[sortIndex[0]] << " ," << W[sortIndex[1]] << " ," << W[sortIndex[2]] << "] ["
                 << Mu[sortIndex[0]] << " ," << Mu[sortIndex[1]] << " ," << Mu[sortIndex[2]] << "] ["
                 << Sigma[sortIndex[0]] << " ," << Sigma[sortIndex[1]] << " ," << Sigma[sortIndex[2]] << "]\n";

            //Find boundaries, there are ngaussians-1 of them, but for the moment this is defined to be 2
            //This finds the 'first' root of the difference between two gaussians
            BoundFinder SM(Mu[sortIndex[0]], Sigma[sortIndex[0]], W[sortIndex[0]], Mu[sortIndex[1]], Sigma[sortIndex[1]], W[sortIndex[1]], 1000, 1e-4, false);
            shortMedBound[pop] = SM.findBoundary();
            BoundFinder ML(Mu[sortIndex[1]], Sigma[sortIndex[1]], W[sortIndex[1]], Mu[sortIndex[2]], Sigma[sortIndex[2]], W[sortIndex[2]], 1000, 1e-4, false);
            medLongBound[pop] = ML.findBoundary();

            cerr << rohLengthByPop->at(pop)->pop << " A/B: " << shortMedBound[pop] << " B/C: " << medLongBound[pop] << endl;
            fout << rohLengthByPop->at(pop)->pop << " " << shortMedBound[pop] << " " << medLongBound[pop] << endl;
        }
        cerr << "Wrote " << sizeBoundaryOutfile << endl;
        fout.close();
    }
    else
    {
        for (map<string, double>::iterator it = pop2SMbound.begin(); it != pop2SMbound.end(); it++)
        {
            shortMedBound[pop2index[it->first]] = pop2SMbound[it->first];
            medLongBound[pop2index[it->first]] = pop2MLbound[it->first];
            cerr << "User provided size boundaries. " << it->first
                 << " A/B: " << shortMedBound[pop2index[it->first]]
                 << " B/C: " << medLongBound[pop2index[it->first]] << endl;
        }
    }

    //Output ROH calls to file, one for each individual
    //includes A/B/C size classifications
    cerr << "Writing ROH tracts...\n";
    writeROHData(outfile, rohDataByPopByInd, mapDataByChr, shortMedBound, medLongBound, ind2pop);

    return 0;

}

