#include "garlic-roh.h"

//pthread_mutex_t io_mutex = PTHREAD_MUTEX_INITIALIZER;

bool inGap(int qStart, int qEnd, int targetStart, int targetEnd)
{
    return ( (targetStart <= qStart && targetEnd >= qStart) ||
             (targetStart <= qEnd && targetEnd >= qEnd) ||
             (targetStart >= qStart && targetEnd <= qEnd) );
}

void calcLOD(IndData *indData, MapData *mapData,
             HapData *hapData, FreqData *freqData,
             GenoLikeData *GLData,
             WinData *winData, centromere *centro,
             int winsize, double error, int MAX_GAP, bool USE_GL)
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
                    if (USE_GL) error = GLData->data[i][ind];
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
                        if (USE_GL) {
                            win[ind][locus] = win[ind][locus - 1] -
                                              lod(data[locus - 1][ind], freq[locus - 1], GLData->data[locus - 1][ind]) +
                                              lod(data[locus + winsize - 1][ind], freq[locus + winsize - 1], GLData->data[locus + winsize - 1][ind]);

                        }
                        else {
                            win[ind][locus] = win[ind][locus - 1] -
                                              lod(data[locus - 1][ind], freq[locus - 1], error) +
                                              lod(data[locus + winsize - 1][ind], freq[locus + winsize - 1], error);
                        }
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
                        if (USE_GL) error = GLData->data[i][ind];
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

double nomut(double M, double mu, double interval) {
    return exp(-2.0 * M * mu * interval);
}

double norec(double M, double interval) {
    return nomut(M, 1, interval);
}

double ld(HapData *hapData, GenoFreqData *genoFreqData, int site, int start, int end) {
    double sum = 0;
    for (int i = start; i <= end; i++) {
        if (i != site) sum += hr2(hapData, genoFreqData, i, site);
        else sum += 1;
    }
    //sum *= 2.0;
    return (1.0 / sum);
}

double hr2(HapData *hapData, GenoFreqData *genoFreqData, int i, int j) {
    double HA = genoFreqData->homFreq[i];
    double HB = genoFreqData->homFreq[j];
    double HAB = 0;
    double total = 0;
    for (int ind = 0; ind < hapData->nind; ind++) {
        if (hapData->data[i][ind] != -9 && hapData->data[j][ind] != -9) {
            total++;
            if (hapData->data[i][ind] != 1 && hapData->data[j][ind] != 1) {
                HAB++;
            }
        }
    }
    HAB /= total;
    double H = HAB - HA * HB;
    return H * H / (HA * (1 - HA) * HB * (1 - HB));
}

void calcwLOD(IndData *indData, MapData *mapData,
              HapData *hapData, FreqData *freqData,
              GenoLikeData *GLData,
              GenoFreqData *genoFreqData,
              WinData *winData, centromere *centro,
              int winsize, double error, int MAX_GAP, bool USE_GL)
{
    double mu = 1e-9;
    int M = 7;
    short **data = hapData->data;
    int nloci = hapData->nloci;
    int nind = hapData->nind;
    int *physicalPos = mapData->physicalPos;
    double *geneticPos = mapData->geneticPos;
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
                    if (USE_GL) error = GLData->data[i][ind];
                    double physInterval = ( i > 0 ) ? (physicalPos[i] - physicalPos[i - 1]) : physicalPos[i];
                    double geneInterval = ( i > 0 ) ? (geneticPos[i] - geneticPos[i - 1]) : geneticPos[i];
                    win[ind][locus] += lod(data[i][ind], freq[i], error) * nomut(M, mu, physInterval) * norec(M, geneInterval) * ld(hapData, genoFreqData, i, locus, locus + winsize - 1);
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
                        double physIntervalPrev = ( locus - 1 > 0 ) ? (physicalPos[locus - 1] - physicalPos[locus - 2]) : physicalPos[locus - 1];
                        double geneIntervalPrev = ( locus - 1 > 0 ) ? (geneticPos[locus - 1] - geneticPos[locus - 2]) : geneticPos[locus - 1];
                        double physInterval = ( locus + winsize - 1 > 0 ) ? (physicalPos[locus + winsize - 1] - physicalPos[locus + winsize - 2]) : physicalPos[locus + winsize - 1];
                        double geneInterval = ( locus + winsize - 1 > 0 ) ? (geneticPos[locus + winsize - 1] - geneticPos[locus + winsize - 2]) : geneticPos[locus + winsize - 1];

                        if (USE_GL) {
                            win[ind][locus] = win[ind][locus - 1] -
                                              (lod(data[locus - 1][ind], freq[locus - 1], GLData->data[locus - 1][ind]) *
                                               nomut(M, mu, physIntervalPrev) *
                                               norec(M, geneIntervalPrev) *
                                               ld(hapData, genoFreqData, locus - 1, locus - 1, locus + winsize - 2)) +
                                              (lod(data[locus + winsize - 1][ind], freq[locus + winsize - 1], GLData->data[locus + winsize - 1][ind]) *
                                               nomut(M, mu, physInterval) *
                                               norec(M, geneInterval) *
                                               ld(hapData, genoFreqData, locus + winsize - 1, locus, locus + winsize - 1));

                        }
                        else {
                            win[ind][locus] = win[ind][locus - 1] -
                                              (lod(data[locus - 1][ind], freq[locus - 1], error) *
                                               nomut(M, mu, physIntervalPrev) *
                                               norec(M, geneIntervalPrev) *
                                               ld(hapData, genoFreqData, locus - 1, locus - 1, locus + winsize - 2) ) +
                                              (lod(data[locus + winsize - 1][ind], freq[locus + winsize - 1], error) *
                                               nomut(M, mu, physInterval) *
                                               norec(M, geneInterval) *
                                               ld(hapData, genoFreqData, locus + winsize - 1, locus, locus + winsize - 1));
                        }
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
                        if (USE_GL) error = GLData->data[i][ind];
                        //win[ind][locus] += lod(data[i][ind], freq[i], error);

                        double physInterval = ( i > 0 ) ? (physicalPos[i] - physicalPos[i - 1]) : physicalPos[i];
                        double geneInterval = ( i > 0 ) ? (geneticPos[i] - geneticPos[i - 1]) : geneticPos[i];
                        win[ind][locus] += lod(data[i][ind], freq[i], error) * nomut(M, mu, physInterval) * norec(M, geneInterval) * ld(hapData, genoFreqData, i, locus, locus + winsize - 1);

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
                                    vector< GenoLikeData * > *GLDataByChr,
                                    IndData *indData,
                                    centromere *centro,
                                    int winsize, double error, int MAX_GAP, bool USE_GL)
{
    vector< WinData * > *winDataByChr = initWinData(mapDataByChr, indData);

    for (unsigned int chr = 0; chr < winDataByChr->size(); chr++)
    {
        calcLOD(indData, mapDataByChr->at(chr),
                hapDataByChr->at(chr), freqDataByChr->at(chr),
                GLDataByChr->at(chr),
                winDataByChr->at(chr), centro,
                winsize, error, MAX_GAP, USE_GL);
    }
    return winDataByChr;
}

vector< WinData * > *calcwLODWindows(vector< HapData * > *hapDataByChr,
                                     vector< FreqData * > *freqDataByChr,
                                     vector< MapData * > *mapDataByChr,
                                     vector< GenoLikeData * > *GLDataByChr,
                                     vector< GenoFreqData * > *genoFreqDataByChr,
                                     IndData *indData,
                                     centromere *centro,
                                     int winsize, double error, int MAX_GAP, bool USE_GL)
{
    vector< WinData * > *winDataByChr = initWinData(mapDataByChr, indData);

    for (unsigned int chr = 0; chr < winDataByChr->size(); chr++)
    {
        calcwLOD(indData, mapDataByChr->at(chr),
                 hapDataByChr->at(chr), freqDataByChr->at(chr),
                 GLDataByChr->at(chr),
                 genoFreqDataByChr->at(chr),
                 winDataByChr->at(chr), centro,
                 winsize, error, MAX_GAP, USE_GL);
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

void releaseROHData(vector< ROHData * > *rohDataByInd)
{
    for (unsigned int ind = 0; ind < rohDataByInd->size(); ind++)
    {
        delete rohDataByInd->at(ind);
    }
    rohDataByInd->clear();
    delete rohDataByInd;
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

double selectLODCutoff(KDEResult *kdeResult)
{
    double LOD_CUTOFF;
    try { LOD_CUTOFF = get_min_btw_modes(kdeResult->x, kdeResult->y, 512); }
    catch (...)
    {
        LOG.err("ERROR: Failed to find the minimum between modes in the LOD score density.");
        LOG.err("\tResults from density estimation have been written to file for inspection.");
        LOG.err("\tA cutoff can be manually specified on the command line with", ARG_LOD_CUTOFF);
        return -1;
    }
    return LOD_CUTOFF;
}


double selectLODCutoff(vector< WinData * > *winDataByChr, IndData *indData, int KDE_SUBSAMPLE, string kdeoutfile)
{
    //Format the LOD window data into a single array per pop with no missing data
    //Prepped for KDE
    DoubleData *rawWinData;
    double LOD_CUTOFF;

    if (KDE_SUBSAMPLE <= 0) rawWinData = convertWinData2DoubleData(winDataByChr);
    else rawWinData = convertSubsetWinData2DoubleData(winDataByChr, indData, KDE_SUBSAMPLE);

    //Compute KDE of LOD score distribution
    cout << "Estimating distribution of raw LOD score windows:\n";
    KDEResult *kdeResult = computeKDE(rawWinData->data, rawWinData->size);
    releaseDoubleData(rawWinData);

    //Output kde points
    try { writeKDEResult(kdeResult, kdeoutfile); }
    catch (...) { return -1; }

    try { LOD_CUTOFF = get_min_btw_modes(kdeResult->x, kdeResult->y, 512); }
    catch (...)
    {
        LOG.err("ERROR: Failed to find the minimum between modes in the LOD score density.");
        LOG.err("\tResults from density estimation have been written to file for inspection.");
        LOG.err("\tA cutoff can be manually specified on the command line with", ARG_LOD_CUTOFF);
        return -1;
    }

    releaseKDEResult(kdeResult);
    return LOD_CUTOFF;
}

/*
void winsizeExplore(vector< HapData * > *hapDataByChr, vector< FreqData * > *freqDataByChr,
                    vector< MapData * > *mapDataByChr, IndData *indData,
                    centromere *centro, vector<int> *multiWinsizes,
                    KDEWinsizeReport *winsizeReport, double error, int MAX_GAP,
                    int KDE_SUBSAMPLE, int numThreads, bool WINSIZE_EXPLORE, string outfile)
{

}
*/

/*
KDEResult *automaticallyChooseWindowSizeFromCustomList();

KDEResult *automaticallyChooseWindowSize(vector< HapData * > *hapDataByChr, vector< FreqData * > *freqDataByChr,
        vector< MapData * > *mapDataByChr, IndData *indData,
        centromere *centro, int &winsize, double error, int MAX_GAP,
        int KDE_SUBSAMPLE, int numThreads, bool WINSIZE_EXPLORE, double AUTO_WINSIZE_THRESHOLD, string outfile)
{
    KDEWinsizeReport *winsizeReport;
    KDEResult *kdeResult;
    int selectedWinsize = -1;
    int stepSize = 10;
    int lastWinsize = winsize - stepSize;
    vector<int> *multiWinsizes;
    bool finished = false;
    do
    {
        multiWinsizes = getWinsizeList(lastWinsize, stepSize, numThreads);
        cerr << "before.\n";
        winsizeReport = calculateLODOverWinsizeRange(hapDataByChr, freqDataByChr,
                        mapDataByChr, indData, centro, multiWinsizes, error, MAX_GAP,
                        KDE_SUBSAMPLE, numThreads, WINSIZE_EXPLORE, outfile);
        cerr << "end\n";
        lastWinsize = multiWinsizes->at(multiWinsizes->size() - 1);
        multiWinsizes->clear();
        delete multiWinsizes;

        selectedWinsize = selectWinsize(winsizeReport, AUTO_WINSIZE_THRESHOLD);
        if (selectedWinsize > 0)
        {
            winsize = selectedWinsize;
            kdeResult = cloneKDEResult(winsizeReport->kdeResultByWinsize->at(selectedWinsize));
            finished = true;
        }
        releaseKDEWinsizeReport(winsizeReport);

    } while (!finished);

    return kdeResult;
}

vector<int> *getWinsizeList(int lastWinsize, int stepSize, int numThreads)
{
    vector<int> *winsizeList = new vector<int>;
    for (int i = 1; i <= numThreads; i++)
    {
        winsizeList->push_back(lastWinsize + (i * stepSize));
        //cerr << lastWinsize + (i * stepSize) << " ";
    }
    //cerr << endl;
    return winsizeList;
}
*/
/*
int selectWinsize(KDEWinsizeReport *winsizeReport, double AUTO_WINSIZE_THRESHOLD)
{
    map<int, double>::iterator it;
    int winsize = -1;
    double diff = numeric_limits<double>::max();
    for (it = winsizeReport->win2mse->begin(); it != winsizeReport->win2mse->end(); it++)
    {
        if (it->second < AUTO_WINSIZE_THRESHOLD && (AUTO_WINSIZE_THRESHOLD - it->second) < diff)
        {
            winsize = it->first;
            diff = AUTO_WINSIZE_THRESHOLD - it->second;
        }
    }
    return winsize;
}
*/
/*
KDEWinsizeReport *calculateLODOverWinsizeRange(vector< HapData * > *hapDataByChr, vector< FreqData * > *freqDataByChr,
        vector< MapData * > *mapDataByChr, IndData *indData,
        centromere *centro, vector<int> *multiWinsizes, double error, int MAX_GAP,
        int KDE_SUBSAMPLE, int numThreads, bool WINSIZE_EXPLORE, string outfile)
{
    vector< HapData * > *hapDataByChrToCalc;
    IndData *indDataToCalc;

    if (KDE_SUBSAMPLE > 0) subsetData(hapDataByChr, indData, &hapDataByChrToCalc, &indDataToCalc, KDE_SUBSAMPLE);
    else
    {
        hapDataByChrToCalc = hapDataByChr;
        indDataToCalc = indData;
    }

    KDEWinsizeReport *winsizeReport = initKDEWinsizeReport();
    work_order_t *order;
    vector< work_order_t * > orderList;
    pthread_t *peer = new pthread_t[numThreads];
    for (int i = 0; i < numThreads; i++)
    {
        cerr << "set up thread " << i << endl;
        order = new work_order_t;
        order->multiWinsizes = multiWinsizes;
        order->error = error;
        order->MAX_GAP = MAX_GAP;

        order->indData = indDataToCalc;
        order->mapDataByChr = mapDataByChr;
        order->hapDataByChr = hapDataByChrToCalc;
        order->freqDataByChr = freqDataByChr;
        order->winsizeReport = winsizeReport;
        order->centro = centro;
        order->numThreads = numThreads;
        order->WINSIZE_EXPLORE = WINSIZE_EXPLORE;
        order->winsizeReport = winsizeReport;
        order->outfile = outfile;
        order->id = i;
        orderList.push_back(order);
        pthread_create(&(peer[i]),
                       NULL,
                       (void *(*)(void *))compute,
                       (void *)order);
        order = NULL;
    }

    for (int i = 0; i < numThreads; i++)
    {
        pthread_join(peer[i], NULL);
        delete orderList[i];
        cerr << "ended " << i << endl;
    }

    delete [] peer;

    if (KDE_SUBSAMPLE > 0) {
        releaseHapData(hapDataByChrToCalc);
        releaseIndData(indDataToCalc);
    }

    hapDataByChrToCalc = NULL;
    indDataToCalc = NULL;

    return winsizeReport;
}
*/
/*
void compute(void *order)
{
    work_order_t *w = (work_order_t *)order;
    int id = w->id;
    IndData *indData = w->indData;
    vector< MapData * > *mapDataByChr = w->mapDataByChr;
    vector< HapData * > *hapDataByChr = w->hapDataByChr;
    vector< FreqData * > *freqDataByChr = w->freqDataByChr;
    KDEWinsizeReport *winsizeReport = w->winsizeReport;
    vector<int> *multiWinsizes = w->multiWinsizes;
    double error = w->error;
    int MAX_GAP = w->MAX_GAP;
    centromere *centro = w->centro;
    int numThreads = w->numThreads;
    bool WINSIZE_EXPLORE = w->WINSIZE_EXPLORE;
    string outfile = w->outfile;

    pthread_mutex_lock(&io_mutex);
    cerr << "thread " << id << endl;
    pthread_mutex_unlock(&io_mutex);

    vector< WinData * > *winDataByChr;
    for (unsigned int i = id; i < multiWinsizes->size(); i += numThreads)
    {

        winDataByChr = calcLODWindows(hapDataByChr, freqDataByChr, mapDataByChr,
                                      indData, centro, multiWinsizes->at(i),
                                      error, MAX_GAP);

        pthread_mutex_lock(&io_mutex);
        cerr << "thread " << id << " LOD windows completed." << endl;
        pthread_mutex_unlock(&io_mutex);

        DoubleData *rawWinData = convertWinData2DoubleData(winDataByChr);
        releaseWinData(winDataByChr);

        pthread_mutex_lock(&io_mutex);
        KDEResult *kdeResult = computeKDE(rawWinData->data, rawWinData->size);
        releaseDoubleData(rawWinData);
        pthread_mutex_unlock(&io_mutex);

        pthread_mutex_lock(&io_mutex);
        winsizeReport->kdeResultByWinsize->operator[](multiWinsizes->at(i)) = kdeResult;
        winsizeReport->win2mse->operator[](multiWinsizes->at(i)) = calculateWiggle(kdeResult);
        try { if (WINSIZE_EXPLORE) writeKDEResult(kdeResult, makeKDEFilename(outfile, multiWinsizes->at(i))); }
        catch (...) { throw 0; }
        pthread_mutex_unlock(&io_mutex);
    }

    return;
}
*/
void exploreWinsizes(vector< HapData * > *hapDataByChr,
                     vector< FreqData * > *freqDataByChr,
                     vector< MapData * > *mapDataByChr,
                     IndData *indData,
                     centromere *centro,
                     vector<int> &multiWinsizes,
                     double error,
                     vector< GenoLikeData * > *GLDataByChr, bool USE_GL,
                     int MAX_GAP, int KDE_SUBSAMPLE, string outfile,
                     bool WEIGHTED, vector< GenoFreqData * > *genoFreqDataByChr)
{
    vector< WinData * > *winDataByChr;
    vector< HapData * > *hapDataByChrToCalc;
    vector< GenoLikeData * > *GLDataByChrToCalc;
    IndData *indDataToCalc;

    if (KDE_SUBSAMPLE > 0) subsetData(hapDataByChr, GLDataByChr, indData, &hapDataByChrToCalc, &GLDataByChrToCalc, &indDataToCalc, KDE_SUBSAMPLE, USE_GL);
    else
    {
        hapDataByChrToCalc = hapDataByChr;
        if (USE_GL) GLDataByChrToCalc = GLDataByChr;
        indDataToCalc = indData;
    }


    for (unsigned int i = 0; i < multiWinsizes.size(); i++)
    {
        if (WEIGHTED) {
            winDataByChr = calcwLODWindows(hapDataByChrToCalc, freqDataByChr, mapDataByChr,
                                           GLDataByChrToCalc,
                                           genoFreqDataByChr,
                                           indDataToCalc, centro, multiWinsizes[i],
                                           error, MAX_GAP, USE_GL);
        }
        else {
            winDataByChr = calcLODWindows(hapDataByChrToCalc, freqDataByChr, mapDataByChr,
                                          GLDataByChrToCalc,
                                          indDataToCalc, centro, multiWinsizes[i],
                                          error, MAX_GAP, USE_GL);
        }
        DoubleData *rawWinData = convertWinData2DoubleData(winDataByChr);
        releaseWinData(winDataByChr);

        KDEResult *kdeResult = computeKDE(rawWinData->data, rawWinData->size);
        releaseDoubleData(rawWinData);

        try { writeKDEResult(kdeResult, makeKDEFilename(outfile, multiWinsizes[i])); }
        catch (...) { throw 0; }
    }

    if (KDE_SUBSAMPLE > 0) {
        releaseHapData(hapDataByChrToCalc);
        releaseIndData(indDataToCalc);
        if (USE_GL) releaseGLData(GLDataByChrToCalc);
    }

    hapDataByChrToCalc = NULL;
    indDataToCalc = NULL;
    GLDataByChrToCalc = NULL;
    return;
}


KDEResult *selectWinsize(vector< HapData * > *hapDataByChr,
                         vector< FreqData * > *freqDataByChr,
                         vector< MapData * > *mapDataByChr,
                         IndData *indData, centromere *centro,
                         int &winsize, int step, double error,
                         vector< GenoLikeData * > *GLDataByChr, bool USE_GL,
                         int MAX_GAP, int KDE_SUBSAMPLE, string outfile,
                         bool WEIGHTED, vector< GenoFreqData * > *genoFreqDataByChr)
{
    double AUTO_WINSIZE_THRESHOLD = 0.5;
    vector< WinData * > *winDataByChr = NULL;
    vector< HapData * > *hapDataByChrToCalc = NULL;
    vector< GenoLikeData * > *GLDataByChrToCalc = NULL;
    IndData *indDataToCalc = NULL;
    KDEResult *selectedKDEResult = NULL;

    //subset of individuals as given by --kde-subsample
    if (KDE_SUBSAMPLE > 0) subsetData(hapDataByChr, GLDataByChr, indData, &hapDataByChrToCalc, &GLDataByChrToCalc, &indDataToCalc, KDE_SUBSAMPLE, USE_GL);
    else
    {
        hapDataByChrToCalc = hapDataByChr;
        if (USE_GL) GLDataByChrToCalc = GLDataByChr;
        indDataToCalc = indData;
    }


    LOG.log("Searching for acceptable window size, smoothness threshold:", AUTO_WINSIZE_THRESHOLD);
    LOG.log("winsize\tsmoothness");

    int winsizeQuery = winsize;
    double mse;
    bool finished = false;
    while (!finished)
    {
        if (WEIGHTED) {
            winDataByChr = calcwLODWindows(hapDataByChrToCalc, freqDataByChr, mapDataByChr,
                                           GLDataByChrToCalc,
                                           genoFreqDataByChr,
                                           indDataToCalc, centro, winsizeQuery,
                                           error, MAX_GAP, USE_GL);
        }
        else {
            winDataByChr = calcLODWindows(hapDataByChrToCalc, freqDataByChr, mapDataByChr,
                                          GLDataByChrToCalc,
                                          indDataToCalc, centro, winsizeQuery,
                                          error, MAX_GAP, USE_GL);
        }
        DoubleData *rawWinData = convertWinData2DoubleData(winDataByChr);
        releaseWinData(winDataByChr);

        KDEResult *kdeResult = computeKDE(rawWinData->data, rawWinData->size);
        releaseDoubleData(rawWinData);

        mse = calculateWiggle(kdeResult);
        LOG.logn(winsizeQuery);
        LOG.logn("\t");
        LOG.log(mse);

        if (mse <= AUTO_WINSIZE_THRESHOLD)
        {
            finished = true;
            selectedKDEResult = cloneKDEResult(kdeResult);
            winsize = winsizeQuery;
            try { writeKDEResult(selectedKDEResult, makeKDEFilename(outfile, winsize)); }
            catch (...) { throw 0; }
        }
        else winsizeQuery += step;
        releaseKDEResult(kdeResult);
    }

    if (KDE_SUBSAMPLE > 0) {
        releaseHapData(hapDataByChrToCalc);
        releaseIndData(indDataToCalc);
        if (USE_GL) releaseGLData(GLDataByChrToCalc);
    }

    hapDataByChrToCalc = NULL;
    indDataToCalc = NULL;
    GLDataByChrToCalc = NULL;

    return selectedKDEResult;
}

KDEResult *selectWinsizeFromList(vector< HapData * > *hapDataByChr,
                                 vector< FreqData * > *freqDataByChr,
                                 vector< MapData * > *mapDataByChr,
                                 IndData *indData, centromere *centro,
                                 vector<int> *multiWinsizes, int &winsize, double error,
                                 vector< GenoLikeData * > *GLDataByChr, bool USE_GL,
                                 int MAX_GAP, int KDE_SUBSAMPLE, string outfile,
                                 bool WEIGHTED, vector< GenoFreqData * > *genoFreqDataByChr)
{
    double AUTO_WINSIZE_THRESHOLD = 0.5;
    vector< WinData * > *winDataByChr = NULL;
    vector< HapData * > *hapDataByChrToCalc = NULL;
    vector< GenoLikeData * > *GLDataByChrToCalc = NULL;
    IndData *indDataToCalc = NULL;
    KDEResult *selectedKDEResult = NULL;

    //subset of individuals as given by --kde-subsample
    if (KDE_SUBSAMPLE > 0) subsetData(hapDataByChr, GLDataByChr, indData, &hapDataByChrToCalc, &GLDataByChrToCalc, &indDataToCalc, KDE_SUBSAMPLE, USE_GL);
    else
    {
        hapDataByChrToCalc = hapDataByChr;
        if (USE_GL) GLDataByChrToCalc = GLDataByChr;
        indDataToCalc = indData;
    }

    LOG.log("Searching for acceptable window size, smoothness threshold:", AUTO_WINSIZE_THRESHOLD);
    LOG.log("winsize\tsmoothness");

    double mse;
    for (unsigned int i = 0; i < multiWinsizes->size(); i++)
    {
        if (WEIGHTED) {
            winDataByChr = calcwLODWindows(hapDataByChrToCalc, freqDataByChr, mapDataByChr,
                                           GLDataByChrToCalc,
                                           genoFreqDataByChr,
                                           indDataToCalc, centro, multiWinsizes->at(i),
                                           error, MAX_GAP, USE_GL);
        }
        else {
            winDataByChr = calcLODWindows(hapDataByChrToCalc, freqDataByChr, mapDataByChr,
                                          GLDataByChrToCalc,
                                          indDataToCalc, centro, multiWinsizes->at(i),
                                          error, MAX_GAP, USE_GL);
        }
        DoubleData *rawWinData = convertWinData2DoubleData(winDataByChr);
        releaseWinData(winDataByChr);

        KDEResult *kdeResult = computeKDE(rawWinData->data, rawWinData->size);
        releaseDoubleData(rawWinData);

        mse = calculateWiggle(kdeResult);
        LOG.logn(multiWinsizes->at(i));
        LOG.logn("\t");
        LOG.log(mse);

        if (mse <= AUTO_WINSIZE_THRESHOLD || i == multiWinsizes->size() - 1)
        {
            selectedKDEResult = cloneKDEResult(kdeResult);
            winsize = multiWinsizes->at(i);
            try { writeKDEResult(selectedKDEResult, makeKDEFilename(outfile, winsize)); }
            catch (...) { throw 0; }
            releaseKDEResult(kdeResult);
            break;
        }
        releaseKDEResult(kdeResult);
    }

    if (KDE_SUBSAMPLE > 0) {
        releaseHapData(hapDataByChrToCalc);
        releaseIndData(indDataToCalc);
        if (USE_GL) releaseGLData(GLDataByChrToCalc);
    }

    hapDataByChrToCalc = NULL;
    indDataToCalc = NULL;
    GLDataByChrToCalc = NULL;

    return selectedKDEResult;
}

int_pair_t selectSizeClasses(ROHLength *rohLength)
{
    int_pair_t bounds;
    size_t *sortIndex;

    int ngaussians = 3;
    size_t maxIter = 1000;
    double tolerance = 1e-5;
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

    GMM gmm(ngaussians, W, Mu, Sigma, maxIter, tolerance, true);

    gmm.estimate(rohLength->length, rohLength->size);

    for (int n = 0; n < ngaussians; n++)
    {
        W[n] = gmm.getMixCoefficient(n);
        Mu[n] = gmm.getMean(n);
        Sigma[n] = gmm.getVar(n);
        sortIndex[n] = n;
    }

    gsl_sort_index(sortIndex, Mu, 1, ngaussians);

    LOG.log("Gaussian class A ( mixture, mean, std ) = (", W[sortIndex[0]], false);
    LOG.log(",", Mu[sortIndex[0]], false);
    LOG.log(",", Sigma[sortIndex[0]], false);
    LOG.log(" )");

    LOG.log("Gaussian class B ( mixture, mean, std ) = (", W[sortIndex[1]], false);
    LOG.log(",", Mu[sortIndex[1]], false);
    LOG.log(",", Sigma[sortIndex[1]], false);
    LOG.log(" )");

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
    delete [] sortIndex;
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

