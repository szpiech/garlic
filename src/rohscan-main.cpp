#include <iostream>
#include <cstdio>
#include <fstream>
#include <cmath>
#include <pthread.h>
#include "rohscan-data.h"
#include "rohscan-roh.h"
#include "rohscan-kde.h";
#include "param_t.h"
#include "gmm.h"
#include "gsl/gsl_statistics.h"
#include "gsl/gsl_sort.h"
#include "BoundFinder.h"

using namespace std;

const string ARG_MAPFILE = "--map";
const string DEFAULT_MAPFILE = "__mapfile";
const string HELP_MAPFILE = "A mapfile with one row per variant site.\n\
\tFormatted <chr#> <locusID> <genetic pos> <physical pos>";

const string ARG_HAPFILE = "--hap";
const string DEFAULT_HAPFILE = "__hapfile";
const string HELP_HAPFILE = "A hapfile with one row per individual,\n\
\tand one column per variant.\n\
\tVariants should be coded 0/1/-9.";

const string ARG_INDFILE = "--ind";
const string DEFAULT_INDFILE = "__indfile";
const string HELP_INDFILE = "An indfile containing population and individual IDs.\n\
\tOne row per individual, formatted <popID> <indID>";

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

const string ARG_POINTS = "--kde-points";
const int DEFAULT_POINTS = 512;
const string HELP_POINTS = "The number of equally spaced points at which to do KDE.";

const string ARG_BW = "--kde-bw";
const double DEFAULT_BW = -1;
const string HELP_BW = "Manually set the bandwidth for the KDE of lod scores.\n\
\tBy default, the nrd0 rule of thumb is used.";

const string ARG_MAX_GAP = "--max-gap";
const int DEFAULT_MAX_GAP = 200000;
const string HELP_MAX_GAP = "A LOD score window is not calculated if the gap (in bps)\n\
\tbetween two loci is greater than this value.";

const string ARG_RESAMPLE = "--resample";
const int DEFAULT_RESAMPLE = 0;
const string HELP_RESAMPLE = "Number of resamples for estimating allele frequencies.\n\
\tWhen set to 0 (default), rohscan will use allele\n\
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

const string ARG_LOD_CUTOFF_FILE = "--lod-cutoff-file";
const string DEFAULT_LOD_CUTOFF_FILE = "_none";
const string HELP_LOD_CUTOFF_FILE = "For LOD based ROH calling, specify a file with LOD score cutoffs\n\
\tabove which ROH are called for each population.\n\
\tFile format is <pop ID> <cutoff>.\n\
\tBy default, these cutoffs are chosen automatically per population with KDE.";

const string ARG_BOUND_SIZE = "--size-bounds";
const double DEFAULT_BOUND_SIZE = -1;
const string HELP_BOUND_SIZE = "Specify the short/medium and medium/long\n\
\tROH boundaries.  By default, this is chosen automatically\n\
\twith a 3-component GMM.  Must provide 2 numbers.";

const string ARG_BOUND_SIZE_FILE = "--size-bounds-file";
const string DEFAULT_BOUND_SIZE_FILE = "_none";
const string HELP_BOUND_SIZE_FILE = "A file speficying the short/medium and medium/long\n\
\tROH boundaries per population.\n\
\tFile format <pop ID> <short/medium boundary> <medium/long boundary>\n\
\tBy default, this is chosen automatically\n\
\twith a 3-component GMM.  Must provide 2 numbers.";

const string ARG_TPED_MISSING = "--tped-missing";
const string DEFAULT_TPED_MISSING = "0";
const string HELP_TPED_MISSING = "Missing data code for TPED files.";

int main(int argc, char *argv[])
{

    param_t params;

    params.addFlag(ARG_MAPFILE, DEFAULT_MAPFILE, "", HELP_MAPFILE);
    params.addFlag(ARG_HAPFILE, DEFAULT_HAPFILE, "", HELP_HAPFILE);
    params.addFlag(ARG_INDFILE, DEFAULT_INDFILE, "", HELP_INDFILE);
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

    try
    {
        params.parseCommandLine(argc, argv);
    }
    catch (...)
    {
        return -1;
    }


    int winsize = params.getIntFlag(ARG_WINSIZE);
    string mapfile = params.getStringFlag(ARG_MAPFILE);
    string hapfile = params.getStringFlag(ARG_HAPFILE);
    string indfile = params.getStringFlag(ARG_INDFILE);
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
    string TPED_MISSING = params.getStringFlag(ARG_TPED_MISSING);

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
    try
    {
        if (TPED)
        {
            chrCoordList = scanTPEDMapData(tpedfile, numLoci);
            mapDataByChr = readTPEDMapData(tpedfile, chrCoordList);

            scanIndData2(tfamfile, numInd, ind2pop, pop2size);
            indList = new string[numInd];
            indDataByPop = readIndData2(tfamfile, numInd, ind2pop, pop2size, indList, pop2index);

            hapDataByPopByChr = readTPEDHapData2(tpedfile, numLoci, numInd, chrCoordList, indList, ind2pop, pop2size, pop2index, TPED_MISSING);
        }
        else
        {
            chrCoordList = scanMapData(mapfile, numLoci);
            mapDataByChr = readMapData(mapfile, chrCoordList);

            scanIndData2(indfile, numInd, ind2pop, pop2size);
            indList = new string[numInd];
            indDataByPop = readIndData2(indfile, numInd, ind2pop, pop2size, indList, pop2index);

            hapDataByPopByChr = readHapData2(hapfile, numLoci, numInd, chrCoordList, indList, ind2pop, pop2size, pop2index);
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

        freqDataByPopByChr = calcFreqData(hapDataByPopByChr, nresample);

        winDataByPopByChr = initWinData(mapDataByChr, indDataByPop);
    }
    catch (...)
    {
        return 1;
    }


    int numChr = chrCoordList->size();
    int numPop = pop2size.size();
    chrCoordList->clear();
    //indCoordList->clear();
    delete chrCoordList;
    // delete indCoordList;

    //Create a vector of pop/chr pairs
    //These will be distributed across threads for LOD score calculation
    vector<int_pair_t> *popChrPairs = new vector<int_pair_t>;
    int_pair_t pair;
    for (int pop = 0; pop < numPop; pop++)
    {
        for (int chr = 0; chr < numChr; chr++)
        {
            pair.first = pop;
            pair.second = chr;
            popChrPairs->push_back(pair);
        }
    }

    cerr << "There are " << popChrPairs->size() << " pop/chr combinations to compute.\n";

    int numThreadsLODcalc = numThreads;

    if (popChrPairs->size() < numThreads)
    {
        numThreadsLODcalc = popChrPairs->size();
        cerr << "WARNING: there are fewer pop/chr pairs than threads requested.  Running with "
             << numThreads << " threads instead.\n";
    }

    //Partition pop/chr pairs amongst the specified threads
    unsigned long int *NUM_PER_THREAD = new unsigned long int[numThreadsLODcalc];
    unsigned long int div = popChrPairs->size() / numThreadsLODcalc;
    for (int i = 0; i < numThreadsLODcalc; i++) NUM_PER_THREAD[i] = div;
    for (int i = 0; i < (popChrPairs->size()) % numThreadsLODcalc; i++) NUM_PER_THREAD[i]++;


    work_order_t *order;
    pthread_t *peer = new pthread_t[numThreadsLODcalc];
    int prev_index = 0;
    for (int i = 0; i < numThreadsLODcalc; i++)
    {
        order = new work_order_t;
        order->first_index = prev_index;
        order->last_index = prev_index + NUM_PER_THREAD[i];
        prev_index += NUM_PER_THREAD[i];

        order->winsize = winsize;
        order->error = error;
        order->MAX_GAP = MAX_GAP;

        order->popChrPairs = popChrPairs;
        order->indDataByPop = indDataByPop;
        order->mapDataByChr = mapDataByChr;
        order->hapDataByPopByChr = hapDataByPopByChr;
        order->freqDataByPopByChr = freqDataByPopByChr;
        order->winDataByPopByChr = winDataByPopByChr;

        order->id = i;
        pthread_create(&(peer[i]),
                       NULL,
                       (void *(*)(void *))scan,
                       (void *)order);

    }

    for (int i = 0; i < numThreadsLODcalc; i++)
    {
        pthread_join(peer[i], NULL);
    }

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
        //Format the LOD window data into a single array per pop with no missing data
        //Prepped for KDE
        vector < DoubleData * > *rawWinDataByPop = convertWinData2DoubleData(winDataByPopByChr);

        //Compute KDE of LOD score distribution
        cerr << "Estimating distribution of raw LOD score windows:\n";
        vector < KDEResult * > *kdeResultByPop = computeKDE(rawWinDataByPop, indDataByPop, numThreads);
        releaseDoubleData(rawWinDataByPop);

        //Output kde points
        try
        {
            writeKDEResult(kdeResultByPop, indDataByPop, outfile);
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

    cerr << "Begin ROH window assembly.\n";
    //Assemble ROH for each individual in each pop
    vector< ROHLength * > *rohLengthByPop = new vector< ROHLength * >;
    vector< vector< ROHData * >* > *rohDataByPopByInd = assembleROHWindows(winDataByPopByChr,
            mapDataByChr,
            indDataByPop,
            lodScoreCutoffByPop,
            &rohLengthByPop,
            winsize);

    cerr << "Complete.\n";

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
    //Could be modified to allow for arbitrary number of size classifications
    for (int pop = 0; pop < rohDataByPopByInd->size(); pop++)
    {
        vector< ROHData * > *rohDataByInd = rohDataByPopByInd->at(pop);
        for (int ind = 0; ind < rohDataByInd->size(); ind++)
        {
            ROHData *rohData = rohDataByInd->at(ind);
            cerr << "Writing ROH tracts for " << rohData->indID << endl;
            string rohOutfile = outfile;
            rohOutfile += ".";
            rohOutfile += rohData->indID;
            rohOutfile += ".roh";
            ofstream out;
            out.open(rohOutfile.c_str());
            if (out.fail())
            {
                cerr << "ERROR: Failed to open " << rohOutfile << " for writing.\n";
                return -1;
            }

            for (int roh = 0; roh < rohData->chr.size(); roh++)
            {
                int size = (rohData->stop[roh] - rohData->start[roh]);
                char sizeClass = 'C';
                if (size < shortMedBound[pop]) sizeClass = 'A';
                else if (size < medLongBound[pop]) sizeClass = 'B';
                out << rohData->chr[roh] << " " << rohData->start[roh] << " " << rohData->stop[roh] << " " << size << " " << sizeClass << endl;
            }

            out.close();
        }
    }


    return 0;
}

