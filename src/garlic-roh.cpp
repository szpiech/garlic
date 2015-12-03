#include "garlic-roh.h"

pthread_mutex_t cerr_mutex = PTHREAD_MUTEX_INITIALIZER;

bool inGap(int qStart, int qEnd, int targetStart, int targetEnd) {
    return ( (targetStart <= qStart && targetEnd >= qStart) ||
             (targetStart <= qEnd && targetEnd >= qEnd) ||
             (targetStart >= qStart && targetEnd <= qEnd) );
}

void calcLOD(IndData *indData, MapData *mapData,
             HapData *hapData, FreqData *freqData,
             WinData *winData, centromere *centro,
             int winsize, double error, int MAX_GAP)
{
    short **data = hapData->data;
    int nloci = hapData->nloci;
    int nind = hapData->nind;
    int *physicalPos = mapData->physicalPos;
    //double *geneticPos = mapData->geneticPos;
    //string *locusName = mapData->locusName;
    double *freq = freqData->freq;
    int start = 0;
    int stop = mapData->nloci;
    double **win = winData->data;
    //int nmiss = 0;

    int cStart = centro->centromereStart(mapData->chr);
    int cEnd = centro->centromereEnd(mapData->chr);

    //Check if the last window would overshoot the last locus in the data
    if (nloci - stop < winsize) stop = nloci - winsize + 1;

    //For each individual
    for (int ind = 0; ind < nind; ind++)
    {
        //starting locus of the window
        for (int locus = start; locus < stop; locus++)
        {
            win[ind][locus] = 0;

            //First window?  If so we have to calcualte the whole thing
            if (locus == start)
            {
                int prevI = locus;
                for (int i = locus; i < locus + winsize; i++)
                {
                    if (physicalPos[i] - physicalPos[prevI] > MAX_GAP ||
                            inGap(physicalPos[prevI], physicalPos[i], cStart, cEnd))
                    {
                        win[ind][locus] = MISSING;
                        //nmiss++;
                        locus = prevI;
                        break;
                    }
                    win[ind][locus] += lod(data[i][ind], freq[i], error);
                    prevI = i;
                }

                //if(skip) continue;

            }
            else //Otherwise, we can just subtract the locus that falls off and add the new one
            {
                //But first we have to check if the previous window was MISSING
                if (win[ind][locus - 1] != MISSING)
                {
                    //If the gap to the next locus is > MAX_GAP then make the window MISSING
                    if (physicalPos[locus + winsize - 1] - physicalPos[locus + winsize - 2] > MAX_GAP ||
                            inGap(physicalPos[locus + winsize - 2], physicalPos[locus + winsize - 1], cStart, cEnd))
                    {
                        win[ind][locus] = MISSING;
                        //nmiss++;
                        locus = locus + winsize - 2;
                    }
                    else
                    {
                        win[ind][locus] = win[ind][locus - 1] -
                                          lod(data[locus - 1][ind], freq[locus - 1], error) +
                                          lod(data[locus + winsize - 1][ind], freq[locus + winsize - 1], error);
                    }
                }
                else
                {
                    int prevI = locus;
                    for (int i = locus; i < locus + winsize; i++)
                    {
                        if (physicalPos[i] - physicalPos[prevI] > MAX_GAP ||
                                inGap(physicalPos[prevI], physicalPos[i], cStart, cEnd))
                        {
                            win[ind][locus] = MISSING;
                            //nmiss++;
                            locus = prevI;
                            break;
                        }
                        win[ind][locus] += lod(data[i][ind], freq[i], error);
                        prevI = i;
                    }

                    //if(skip) continue;
                }
            }
        }
    }

    //winData->nmiss = nmiss;

    return;
}

vector< WinData * > *calcLODWindows(vector< HapData * > *hapDataByChr,
                                    vector< FreqData * > *freqDataByChr,
                                    vector< MapData * > *mapDataByChr,
                                    IndData *indData,
                                    centromere *centro,
                                    int winsize, double error, int MAX_GAP)
{
    vector< WinData * > *winDataByChr = initWinData(mapDataByChr, indData);

    for (unsigned int chr = 0; chr < winDataByChr->size(); chr++) {
        calcLOD(indData, mapDataByChr->at(chr),
                hapDataByChr->at(chr), freqDataByChr->at(chr),
                winDataByChr->at(chr), centro,
                winsize, error, MAX_GAP);
    }
    return winDataByChr;
}


/*
 * Genotype is 0/1/2 counting the number of alternate alleles
 *
 */
double lod(const short &genotype, const double &freq, const double &error)
{

    double autozygous, nonAutozygous;
    if (freq == 0 || freq == 1)
    {
        autozygous = 1;
        nonAutozygous = 1;
    }
    else if (genotype == 0)
    {
        nonAutozygous = (1 - freq) * (1 - freq);
        autozygous = (1 - error) * (1 - freq) + error * nonAutozygous;
    }
    else if (genotype == 1)
    {
        nonAutozygous = 2 * (freq) * (1 - freq);
        autozygous = error * nonAutozygous;
    }
    else if (genotype == 2)
    {
        nonAutozygous = (freq) * (freq);
        autozygous = (1 - error) * (freq) + error * nonAutozygous;
    }
    else
    {
        autozygous = 1;
        nonAutozygous = 1;
    }

    return log10(autozygous / nonAutozygous);
}

vector< ROHData * > *initROHData(IndData *indData)
{
    vector< ROHData * > *rohDataByInd = new vector< ROHData * >;
    for (int ind = 0; ind < indData->nind; ind++)
    {
        ROHData *rohData = new ROHData;
        rohDataByInd->push_back(rohData);
    }
    return rohDataByInd;
}

vector< ROHData * > *assembleROHWindows(vector< WinData * > *winDataByChr,
                                        vector< MapData * > *mapDataByChr,
                                        IndData *indData,
                                        centromere *centro,
                                        double lodScoreCutoff,
                                        ROHLength **rohLength,
                                        int winSize,
                                        int MAX_GAP)
{
    vector<int> lengths;
    vector< ROHData * > *rohDataByInd = initROHData(indData);

    for (int ind = 0; ind < indData->nind; ind++)
    {
        ROHData *rohData = rohDataByInd->at(ind);
        rohData->indID = indData->indID[ind];

        for (unsigned int chr = 0; chr < winDataByChr->size(); chr++)
        {
            WinData *winData = winDataByChr->at(chr);
            MapData *mapData = mapDataByChr->at(chr);

            int cStart = centro->centromereStart(mapData->chr);
            int cEnd = centro->centromereEnd(mapData->chr);

            //translation of the perl script here###Updated to match trevor's algorithm
            //int winStart = -1;
            //int winStop = -1;
            bool *inWin = new bool[mapData->nloci];
            for (int w = 0; w < mapData->nloci; w++) inWin[w] = false;
            for (int w = 0; w < winData->nloci; w++)
            {
                if (winData->data[ind][w] >= lodScoreCutoff)
                {
                    //There is a faster way to do this but I'm lazy right now
                    //for consecutive w, we obiously don't have to assign true again.
                    for (int i = 0; i < winSize; i++) inWin[w + i] = true;
                }
            }

            int winStart = -1;
            int winStop = -1;
            for (int w = 0; w < mapData->nloci; w++)
            {
                //No window being extended and the snp is in ROH
                //Start the window
                if (winStart < 0 && inWin[w])
                {
                    winStart = mapData->physicalPos[w];
                }
                //Window being extended and snp is not in ROH
                //end the window at w-1
                //reset winStart to -1
                else if (inWin[w] && (mapData->physicalPos[w] - mapData->physicalPos[w - 1] > MAX_GAP ||
                                      inGap(mapData->physicalPos[w - 1], mapData->physicalPos[w], cStart, cEnd)) ) {
                    winStop = mapData->physicalPos[w - 1];
                    int size = winStop - winStart + 1;
                    lengths.push_back(size);
                    rohData->chr.push_back(chr);
                    rohData->start.push_back(winStart);
                    rohData->stop.push_back(winStop);
                    winStop = -1;
                    winStart = mapData->physicalPos[w];
                }
                else if (winStart > 0 && !inWin[w])
                {
                    winStop = mapData->physicalPos[w - 1];
                    int size = winStop - winStart + 1;
                    lengths.push_back(size);
                    rohData->chr.push_back(chr);
                    rohData->start.push_back(winStart);
                    rohData->stop.push_back(winStop);
                    winStart = -1;
                    winStop = -1;
                }
            }

            delete [] inWin;
        }
    }

    ROHLength *rohLengths = initROHLength(lengths.size(), indData->pop);
    for (unsigned int i = 0; i < lengths.size(); i++)
    {
        rohLengths->length[i] = lengths[i];
    }
    (*rohLength) = rohLengths;

    return rohDataByInd;
}

ROHLength *initROHLength(int size, string pop)
{
    ROHLength *rohLength = new ROHLength;
    rohLength->pop = pop;
    rohLength->length = new double[size];
    rohLength->size = size;
    return rohLength;
}

void releaseROHLength(ROHLength *rohLength)
{
    delete [] rohLength->length;
    delete rohLength;
    return;
}

void releaseROHLength(vector< ROHLength * > *rohLengthByPop)
{
    for (unsigned int pop = 0; pop < rohLengthByPop->size(); pop++)
    {
        releaseROHLength(rohLengthByPop->at(pop));
    }
    delete rohLengthByPop;
    return;
}

void writeROHData(string outfile,
                  vector< ROHData * > *rohDataByInd,
                  vector< MapData * > *mapDataByChr,
                  int_pair_t bounds,
                  string popName,
                  string version)
{
    ofstream out;
    out.open(outfile.c_str());
    if (out.fail())
    {
        cerr << "ERROR: Failed to open " << outfile << " for writing.\n";
        LOG.err("ERROR: Failed to open", outfile);
        throw 0;
    }

    for (unsigned int ind = 0; ind < rohDataByInd->size(); ind++)
    {
        ROHData *rohData = rohDataByInd->at(ind);
        out << "track name=\"Ind: " + rohData->indID + " Pop:" + popName +
            " ROH\" description=\"Ind: " + rohData->indID + " Pop:" + popName +
            " ROH from GARLIC v" + version + "\" visibility=2 itemRgb=\"On\"\n";

        for (unsigned int roh = 0; roh < rohData->chr.size(); roh++)
        {
            int size = (rohData->stop[roh] - rohData->start[roh]);
            char sizeClass = 'C';
            string color = "0,76,153";
            if (size < bounds.first) {
                sizeClass = 'A';
                color = "153,0,0";
            }
            else if (size < bounds.second) {
                sizeClass = 'B';
                color = "0,153,0";
            }
            string chr = mapDataByChr->at(rohData->chr[roh])->chr;
            if (chr[0] != 'c' && chr[0] != 'C') chr = "chr" + chr;
            out << chr << "\t" << rohData->start[roh] << "\t" << rohData->stop[roh]
                << "\t" << sizeClass << "\t" << size << "\t.\t0\t0\t" << color << endl;
        }
    }
    LOG.log("ROH calls:", outfile);
    out.close();
    return;
}

string makeROHFilename(string outfile)
{
    outfile += ".roh.bed";
    return outfile;
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
        LOG.err("ERROR: Failed to find the minimum between modes in the LOD score density.");
        LOG.err("\tResults from density estimation have been written to file for inspection.");
        LOG.err("\tA cutoff can be manually specified on the command line with", ARG_LOD_CUTOFF);
        return -1;
    }

    releaseKDEResult(kdeResult);
    return LOD_CUTOFF;
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
    for (unsigned int i = 0; i < multiWinsizes.size(); i++)
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

    cerr << "Gaussian class A ( mixture, mean, std ) = ( "
         << W[sortIndex[0]] << ", " << Mu[sortIndex[0]] << ", " << Sigma[sortIndex[0]] << " )\n";
    LOG.log("Gaussian class A ( mixture, mean, std ) = (", W[sortIndex[0]], false);
    LOG.log(",", Mu[sortIndex[0]], false);
    LOG.log(",", Sigma[sortIndex[0]], false);
    LOG.log(" )");

    cerr << "Gaussian class B ( mixture, mean, std ) = ( "
         << W[sortIndex[1]] << ", " << Mu[sortIndex[1]] << ", " << Sigma[sortIndex[1]] << " )\n";
    LOG.log("Gaussian class B ( mixture, mean, std ) = (", W[sortIndex[1]], false);
    LOG.log(",", Mu[sortIndex[1]], false);
    LOG.log(",", Sigma[sortIndex[1]], false);
    LOG.log(" )");

    cerr << "Gaussian class C ( mixture, mean, std ) = ( "
         << W[sortIndex[2]] << ", " << Mu[sortIndex[2]] << ", " << Sigma[sortIndex[2]] << " )\n";
    LOG.log("Gaussian class C ( mixture, mean, std ) = (", W[sortIndex[2]], false);
    LOG.log(",", Mu[sortIndex[2]], false);
    LOG.log(",", Sigma[sortIndex[2]], false);
    LOG.log(" )");

    //Find boundaries, there are ngaussians-1 of them, but for the moment this is defined to be 2
    //This finds the 'first' root of the difference between two gaussians
    BoundFinder SM(Mu[sortIndex[0]],
                   Sigma[sortIndex[0]],
                   W[sortIndex[0]],
                   Mu[sortIndex[1]],
                   Sigma[sortIndex[1]],
                   W[sortIndex[1]],
                   1000, 1e-4, false);
    bounds.first = SM.findBoundary();

    BoundFinder ML(Mu[sortIndex[1]],
                   Sigma[sortIndex[1]],
                   W[sortIndex[1]],
                   Mu[sortIndex[2]],
                   Sigma[sortIndex[2]],
                   W[sortIndex[2]],
                   1000, 1e-4, false);
    bounds.second = ML.findBoundary();

    delete [] W;
    delete [] Mu;
    delete [] Sigma;

    return bounds;
}


/*
vector< vector< WinData * >* > *calcLODWindowsSinglePop(vector< vector< HapData * >* > *hapDataByPopByChr,
        vector< vector< FreqData * >* > *freqDataByPopByChr,
        vector< MapData * > *mapDataByChr,
        vector< IndData * > *indDataByPop,
        centromere *centro,
        int* winsize, double error, int MAX_GAP, int numThreads, int pop)
{

    int numChr = mapDataByChr->size();
    int numPop = 1;

    vector< vector< WinData * >* > *winDataByPopByChr = initWinData(mapDataByChr, indDataByPop, pop);

    //Create a vector of pop/chr pairs
    //These will be distributed across threads for LOD score calculation
    vector<int_pair_t> *popChrPairs = new vector<int_pair_t>;
    int_pair_t pair;

    for (int chr = 0; chr < numChr; chr++)
    {
        pair.first = pop;
        pair.second = chr;
        popChrPairs->push_back(pair);
    }

    //cerr << "There are " << popChrPairs->size() << " pop/chr combinations to compute.\n";

    int numThreadsLODcalc = numThreads;

    if (popChrPairs->size() < numThreads)
    {
        numThreadsLODcalc = popChrPairs->size();
        //cerr << "WARNING: there are fewer pop/chr pairs than threads requested.  Running with "
        //     << numThreads << " threads instead.\n";
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
        order->centro = centro;

        order->id = i;
        pthread_create(&(peer[i]),
                       NULL,
                       (void *(*)(void *))scanSinglePop,
                       (void *)order);

    }

    for (int i = 0; i < numThreadsLODcalc; i++)
    {
        pthread_join(peer[i], NULL);
    }
    return winDataByPopByChr;
}
*/


/*
void scan(void *order)
{
    work_order_t *w = (work_order_t *)order;
    int id = w->id;
    int first_index = w->first_index;
    int last_index = w->last_index;
    vector<int_pair_t> *popChrPairs = w->popChrPairs;
    vector< IndData * > *indDataByPop = w->indDataByPop;
    vector< MapData * > *mapDataByChr = w->mapDataByChr;
    vector< vector< HapData * >* > *hapDataByPopByChr = w->hapDataByPopByChr;
    vector< vector< FreqData * >* > *freqDataByPopByChr = w->freqDataByPopByChr;
    vector< vector< WinData * >* > *winDataByPopByChr = w->winDataByPopByChr;
    int* winsize = w->winsize;
    double error = w->error;
    int MAX_GAP = w->MAX_GAP;
    centromere *centro = w->centro;

    //pthread_mutex_lock(&cerr_mutex);
    //cerr << "Thread " << id << ":\n";
    //pthread_mutex_unlock(&cerr_mutex);

    for (int i = first_index; i < last_index; i++)
    {
        int pop = popChrPairs->at(i).first;
        int chr = popChrPairs->at(i).second;

        IndData *indData = indDataByPop->at(pop);
        MapData *mapData = mapDataByChr->at(chr);
        HapData *hapData = hapDataByPopByChr->at(pop)->at(chr);
        FreqData *freqData = freqDataByPopByChr->at(pop)->at(chr);
        WinData *winData = winDataByPopByChr->at(pop)->at(chr);

        calcLOD(indData, mapData, hapData, freqData, winData, centro, winsize[pop], error, MAX_GAP);

        pthread_mutex_lock(&cerr_mutex);
        cerr << indDataByPop->at(pop)->pop << " chromosome " << mapDataByChr->at(chr)->chr << " LOD windows finished.\n";
        pthread_mutex_unlock(&cerr_mutex);
    }

    return;
}
*/
/*
//ugly hack
void scanSinglePop(void *order)
{
    work_order_t *w = (work_order_t *)order;
    int id = w->id;
    int first_index = w->first_index;
    int last_index = w->last_index;
    vector<int_pair_t> *popChrPairs = w->popChrPairs;
    vector< IndData * > *indDataByPop = w->indDataByPop;
    vector< MapData * > *mapDataByChr = w->mapDataByChr;
    vector< vector< HapData * >* > *hapDataByPopByChr = w->hapDataByPopByChr;
    vector< vector< FreqData * >* > *freqDataByPopByChr = w->freqDataByPopByChr;
    vector< vector< WinData * >* > *winDataByPopByChr = w->winDataByPopByChr;
    int* winsize = w->winsize;
    double error = w->error;
    int MAX_GAP = w->MAX_GAP;
    centromere *centro = w->centro;

    //pthread_mutex_lock(&cerr_mutex);
    //cerr << "Thread " << id << ":\n";
    //pthread_mutex_unlock(&cerr_mutex);

    for (int i = first_index; i < last_index; i++)
    {
        int pop = popChrPairs->at(i).first;
        int chr = popChrPairs->at(i).second;

        IndData *indData = indDataByPop->at(pop);
        MapData *mapData = mapDataByChr->at(chr);
        HapData *hapData = hapDataByPopByChr->at(pop)->at(chr);
        FreqData *freqData = freqDataByPopByChr->at(pop)->at(chr);
        WinData *winData = winDataByPopByChr->at(0)->at(chr);

        calcLOD(indData, mapData, hapData, freqData, winData, centro, winsize[pop], error, MAX_GAP);

        pthread_mutex_lock(&cerr_mutex);
        cerr << indDataByPop->at(pop)->pop << " chromosome " << mapDataByChr->at(chr)->chr << " LOD windows finished.\n";
        pthread_mutex_unlock(&cerr_mutex);
    }

    return;
}
*/

