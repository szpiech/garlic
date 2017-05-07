#include "garlic-roh.h"

int selectWinsizeWeighted(double density){
   /*
    Window size = 8.3235*log(SNVdensity)+138.0521
    Overlap (%) = 6.375*log(SNVdensity)+63.888;*/ 
    int size = int(8.3235*log(density)+138.0521 + 0.5);
    return (size >= 10 ? size : 10);
}

bool inGap(int qStart, int qEnd, int targetStart, int targetEnd)
{
    return ( (targetStart <= qStart && targetEnd >= qStart) ||
             (targetStart <= qEnd && targetEnd >= qEnd) ||
             (targetStart >= qStart && targetEnd <= qEnd) );
}

void calcLOD(MapData *mapData,
             HapData *hapData, FreqData *freqData,
             GenoLikeData *GLData,
             WinData *winData, centromere *centro,
             int winsize, double error, int MAX_GAP, bool USE_GL)
{
    short **data = hapData->data;
    int nloci = hapData->nloci;
    int nind = hapData->nind;
    double *physicalPos = mapData->physicalPos;
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



void calcwLOD(MapData *mapData,
              HapData *hapData,
              FreqData *freqData,
              GenoLikeData *GLData,
              LDData *LD,
              WinData *winData, centromere *centro,
              int winsize, double error, int MAX_GAP, bool USE_GL, double mu, int M, int numThreads)
{
    unsigned int *NUM_PER_THREAD = make_thread_partition(numThreads, hapData->nloci);

    WLOD_work_order_t *order;
    pthread_t *peer = new pthread_t[numThreads];
    vector< WLOD_work_order_t * > orders;
    unsigned int previous = 0;
    for (int i = 0; i < numThreads; i++)
    {
        order = new WLOD_work_order_t;
        order->mapData = mapData;
        order->hapData = hapData;
        order->freqData = freqData;
        order->GLData = GLData;
        order->LD = LD;
        order->winData = winData;
        order->cStart = centro->centromereStart(mapData->chr);
        order->cEnd = centro->centromereEnd(mapData->chr);
        order->winsize = winsize;
        order->error = error;
        order->MAX_GAP = MAX_GAP;
        order->USE_GL = USE_GL;
        order->mu = mu;
        order->M = M;
        order->start = previous;
        previous += NUM_PER_THREAD[i];
        order->stop = previous;

        pthread_create(&(peer[i]),
                       NULL,
                       (void *(*)(void *))parallelwLOD,
                       (void *)order);
        orders.push_back(order);
    }

    for (int i = 0; i < numThreads; i++){
        pthread_join(peer[i], NULL);
        delete orders[i];
    }

    orders.clear();
    delete [] NUM_PER_THREAD;
    delete [] peer;
    return;
}

void parallelwLOD(void *order){
    WLOD_work_order_t *p = (WLOD_work_order_t *)order;
    GenoLikeData *GLData = p->GLData;
    LDData *LD = p->LD;
    int winsize = p->winsize;
    double error = p->error;
    int MAX_GAP = p->MAX_GAP;
    bool USE_GL = p->USE_GL;
    double mu = p->mu;
    int M = p->M;
    int start = p->start;
    int stop = p->stop;

    short **data = p->hapData->data;
    int nloci = p->hapData->nloci;
    int nind = p->hapData->nind;
    double *physicalPos = p->mapData->physicalPos;
    double *geneticPos = p->mapData->geneticPos;
    double *freq = p->freqData->freq;
    double **win = p->winData->data;

    int cStart = p->cStart;
    int cEnd = p->cEnd;

    //Check if the last window would overshoot the last locus in the data
    if (nloci - stop < winsize) stop = nloci - winsize + 1;

    int size = stop-start+winsize+1;
    double *score;

    //cerr << "XXX " << start << " " << stop << endl;

  //  cerr << "wLOD " << p->mapData->chr << endl;

    //For each individual
    for (int ind = 0; ind < nind; ind++){
        score = new double[size];
        for (int locus = start; locus < ((stop+winsize+1 > nloci) ? nloci : stop+winsize+1); locus++){
            if (USE_GL) error = GLData->data[locus][ind];
            double physInterval = ( locus > 0 ) ? (physicalPos[locus] - physicalPos[locus - 1]) : physicalPos[locus];
            double geneInterval = ( locus > 0 ) ? (geneticPos[locus] - geneticPos[locus - 1]) : geneticPos[locus];
            //if(p->mapData->chr.compare("chr21") == 0 && locus > 0) cerr << nloci << " " << locus << " " <<  geneticPos[locus] << " - " << geneticPos[locus - 1] << " " << geneInterval << " " << physicalPos[locus] << endl;
            score[locus-start] = lod(data[locus][ind], freq[locus], error) * nomut(M, mu, physInterval) * norec(M, geneInterval);
        }

        //starting locus of the window
        for (int locus = start; locus < stop; locus++)
        {
            win[ind][locus] = 0;

            //First window?  If so we have to calcualte the whole thing
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

                win[ind][locus] += score[i-start] * (1.0 / LD->LD[locus][i-locus]);
                prevI = i;
            }
        }

        delete [] score;
    }
}

vector< WinData * > *calcLODWindows(vector< HapData * > *hapDataByChr,
                                    vector< FreqData * > *freqDataByChr,
                                    vector< MapData * > *mapDataByChr,
                                    vector< GenoLikeData * > *GLDataByChr,
                                    centromere *centro,
                                    int winsize, double error, int MAX_GAP, bool USE_GL)
{
    vector< WinData * > *winDataByChr = initWinData(mapDataByChr, hapDataByChr->at(0)->nind);

    for (unsigned int chr = 0; chr < winDataByChr->size(); chr++)
    {
        if(USE_GL){
            calcLOD(mapDataByChr->at(chr),
                    hapDataByChr->at(chr), freqDataByChr->at(chr),
                    GLDataByChr->at(chr),
                    winDataByChr->at(chr), centro,
                    winsize, error, MAX_GAP, USE_GL);
        }
        else{
            calcLOD(mapDataByChr->at(chr),
                    hapDataByChr->at(chr), freqDataByChr->at(chr),
                    NULL,
                    winDataByChr->at(chr), centro,
                    winsize, error, MAX_GAP, USE_GL);
        }
    }
    return winDataByChr;
}

vector< WinData * > *calcwLODWindows(vector< HapData * > *hapDataByChr,
                                     vector< FreqData * > *freqDataByChr,
                                     vector< MapData * > *mapDataByChr,
                                     vector< GenoLikeData * > *GLDataByChr,
                                     vector< LDData * > *ldDataByChr,
                                     centromere *centro,
                                     int winsize, double error, int MAX_GAP, bool USE_GL, 
                                     int M, double mu, int numThreads)
{
    vector< WinData * > *winDataByChr = initWinData(mapDataByChr, hapDataByChr->at(0)->nind);

    for (unsigned int chr = 0; chr < winDataByChr->size(); chr++)
    {
        if(USE_GL){
            calcwLOD(mapDataByChr->at(chr),
                     hapDataByChr->at(chr),
                     freqDataByChr->at(chr),
                     GLDataByChr->at(chr),
                     ldDataByChr->at(chr),
                     winDataByChr->at(chr), centro,
                     winsize, error, MAX_GAP, USE_GL, mu, M, numThreads);
        }
        else{
            calcwLOD(mapDataByChr->at(chr),
                     hapDataByChr->at(chr),
                     freqDataByChr->at(chr),
                     NULL,
                     ldDataByChr->at(chr),
                     winDataByChr->at(chr), centro,
                     winsize, error, MAX_GAP, USE_GL, mu, M, numThreads);   
        }
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
                                        int MAX_GAP,
                                        double OVERLAP_FRAC, bool CM)
{
    vector<double> lengths;
    vector< ROHData * > *rohDataByInd = initROHData(indData);

    double OVERLAP_THRESHOLD = OVERLAP_FRAC * winSize;
    OVERLAP_THRESHOLD = (OVERLAP_THRESHOLD >= 1) ? OVERLAP_THRESHOLD : 1;

    for (int ind = 0; ind < indData->nind; ind++)
    {
        ROHData *rohData = rohDataByInd->at(ind);
        rohData->indID = indData->indID[ind];

        for (unsigned int chr = 0; chr < winDataByChr->size(); chr++)
        {
            WinData *winData = winDataByChr->at(chr);
            MapData *mapData = mapDataByChr->at(chr);
            double *pos;

            if(CM) pos = mapData->geneticPos;
            else pos = mapData->physicalPos;

            int cStart = centro->centromereStart(mapData->chr);
            int cEnd = centro->centromereEnd(mapData->chr);

            //translation of the perl script here###Updated to match trevor's algorithm
            //int winStart = -1;
            //int winStop = -1;
            short *inWin = new short[mapData->nloci];
            for (int w = 0; w < mapData->nloci; w++) inWin[w] = 0;
            for (int w = 0; w < winData->nloci; w++)
            {
                if (winData->data[ind][w] >= lodScoreCutoff)
                {
                    for (int i = 0; i < winSize; i++) inWin[w + i]++;
                }
            }

            double winStart = -1;
            int winStartIndex = -1;
            double winStop = -1;
            int winStopIndex = -1;
            for (int w = 0; w < mapData->nloci; w++)
            {
                //No window being extended and the snp is in ROH
                //Start the window
                if (winStart < 0 && inWin[w] >= OVERLAP_THRESHOLD)
                {
                    //winStart = mapData->physicalPos[w];
                    winStart = pos[w];
                    winStartIndex = w;
                }
                else if (inWin[w] >= OVERLAP_THRESHOLD && (mapData->physicalPos[w] - mapData->physicalPos[w - 1] > MAX_GAP ||
                         inGap(mapData->physicalPos[w - 1], mapData->physicalPos[w], cStart, cEnd)) ) {
                    //winStop = mapData->physicalPos[w - 1];
                    winStop = pos[w - 1];
                    winStopIndex = w - 1;
                    if(winStopIndex - winStartIndex + 1 >= OVERLAP_THRESHOLD){
                        double size = winStop - winStart + (!CM ? 1 : 0);
                        lengths.push_back(size);
                        rohData->chr.push_back(chr);
                        rohData->start.push_back(winStart);
                        rohData->stop.push_back(winStop);
                    }
                    winStop = -1;
                    winStopIndex = -1;
                    winStart = pos[w];
                    winStartIndex = w;
                }
                else if (winStart > 0 && ! (inWin[w] >= OVERLAP_THRESHOLD) )
                {
                    //winStop = mapData->physicalPos[w - 1];
                    winStop = pos[w - 1];
                    winStopIndex = w - 1;
                    if(winStopIndex - winStartIndex + 1 >= OVERLAP_THRESHOLD){
                        double size = winStop - winStart + (!CM ? 1 : 0);
                        lengths.push_back(size);
                        rohData->chr.push_back(chr);
                        rohData->start.push_back(winStart);
                        rohData->stop.push_back(winStop);
                    }
                    winStart = -1;
                    winStartIndex = -1;
                    winStop = -1;
                    winStopIndex = -1;

                }
                else if (winStart > 0 && w + 1 >= mapData->nloci)
                {
                    //winStop = mapData->physicalPos[w];
                    winStop = pos[w];
                    winStopIndex = w;
                    if(winStopIndex - winStartIndex + 1 >= OVERLAP_THRESHOLD){
                        double size = winStop - winStart + (!CM ? 1 : 0);
                        lengths.push_back(size);
                        rohData->chr.push_back(chr);
                        rohData->start.push_back(winStart);
                        rohData->stop.push_back(winStop);
                    }
                    winStart = -1;
                    winStartIndex = -1;
                    winStop = -1;
                    winStopIndex = -1;
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
                  vector< double > bounds,
                  string popName,
                  string version, bool CM)
{
    string colors[9];
    colors[0] = "228,26,28";
    colors[1] = "77,175,74";
    colors[2] = "55,126,184";
    colors[3] = "152,78,163";
    colors[4] = "255,127,0";
    colors[5] = "255,255,51";
    colors[6] = "166,86,40";
    colors[7] = "247,129,191";
    colors[8] = "153,153,153";
   
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
            double size = (rohData->stop[roh] - rohData->start[roh]);
            char sc = 'A';
            char sizeClass = 'X';
            string color = "X";

            unsigned int i = 0;
            for(i = 0; i < bounds.size(); i++){
                if (size < bounds[i]) {
                    sizeClass = sc;
                    color = colors[(i <= 8) ? i : 8];
                    break;
                }
                sc++;
            }

            if(color.compare("X") == 0){
                sizeClass = sc;
                color = colors[(i <= 8) ? i : 8];
            }

            string chr = mapDataByChr->at(rohData->chr[roh])->chr;
            if (chr[0] != 'c' && chr[0] != 'C') chr = "chr" + chr;
            if(CM){
                out << chr << "\t" << rohData->start[roh] << "\t" << rohData->stop[roh]
                    << "\t" << sizeClass << "\t" << size << "\t.\t0\t0\t" << color << endl;
            }
            else{
                out << chr << "\t" << int(rohData->start[roh]) << "\t" << int(rohData->stop[roh])
                    << "\t" << sizeClass << "\t" << int(size) << "\t.\t0\t0\t" << color << endl;
            }
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


double selectLODCutoff(vector< WinData * > *winDataByChr, IndData *indData, int KDE_SUBSAMPLE, string kdeoutfile, int step)
{
    //Format the LOD window data into a single array per pop with no missing data
    //Prepped for KDE
    DoubleData *rawWinData;
    double LOD_CUTOFF;

    if (KDE_SUBSAMPLE <= 0) rawWinData = convertWinData2DoubleData(winDataByChr, step);
    else rawWinData = convertSubsetWinData2DoubleData(winDataByChr, indData, KDE_SUBSAMPLE, step);

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
                     vector< GenoLikeData * > *GLDataByChr,
                     vector< GenoFreqData * > *genoFreqDataByChr, bool USE_GL,
                     int MAX_GAP, int KDE_SUBSAMPLE, string outfile,
                     bool WEIGHTED, int M, double mu, int numThreads, bool PHASED, bool THIN)
{
    vector< WinData * > *winDataByChr;
    vector< HapData * > *hapDataByChrToCalc;
    vector< GenoLikeData * > *GLDataByChrToCalc;
    IndData *indDataToCalc;

    if (KDE_SUBSAMPLE > 0) subsetData(hapDataByChr, GLDataByChr, indData, &hapDataByChrToCalc, &GLDataByChrToCalc, &indDataToCalc, KDE_SUBSAMPLE, USE_GL, PHASED);
    else
    {
        hapDataByChrToCalc = hapDataByChr;
        if (USE_GL) GLDataByChrToCalc = GLDataByChr;
        indDataToCalc = indData;
    }

    vector< LDData * > *ldDataByChr;

    for (unsigned int i = 0; i < multiWinsizes.size(); i++)
    {
        if (WEIGHTED) {
            ldDataByChr = calcLDData(hapDataByChr, freqDataByChr, mapDataByChr, genoFreqDataByChr, centro, multiWinsizes[i], MAX_GAP, PHASED, numThreads);
            winDataByChr = calcwLODWindows(hapDataByChrToCalc, freqDataByChr, mapDataByChr,
                                           GLDataByChrToCalc,
                                           ldDataByChr,
                                           centro, multiWinsizes[i],
                                           error, MAX_GAP, USE_GL, M, mu, numThreads);
            releaseLDData(ldDataByChr);
        }
        else {
            winDataByChr = calcLODWindows(hapDataByChrToCalc, freqDataByChr, mapDataByChr,
                                          GLDataByChrToCalc,
                                          centro, multiWinsizes[i],
                                          error, MAX_GAP, USE_GL);
        }
        DoubleData *rawWinData = convertWinData2DoubleData(winDataByChr, (THIN ? multiWinsizes[i] : 1));
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
                         bool WEIGHTED, vector< GenoFreqData * > *genoFreqDataByChr, bool PHASED, bool THIN)
{
    double AUTO_WINSIZE_THRESHOLD = 0.5;
    vector< WinData * > *winDataByChr = NULL;
    vector< HapData * > *hapDataByChrToCalc = NULL;
    vector< GenoLikeData * > *GLDataByChrToCalc = NULL;
    IndData *indDataToCalc = NULL;
    KDEResult *selectedKDEResult = NULL;

    //subset of individuals as given by --kde-subsample
    if (KDE_SUBSAMPLE > 0) subsetData(hapDataByChr, GLDataByChr, indData, &hapDataByChrToCalc, &GLDataByChrToCalc, &indDataToCalc, KDE_SUBSAMPLE, USE_GL, PHASED);
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
            /*
            winDataByChr = calcwLODWindows(hapDataByChrToCalc, freqDataByChr, mapDataByChr,
                                           GLDataByChrToCalc,
                                           genoFreqDataByChr,
                                           centro, winsizeQuery,
                                           error, MAX_GAP, USE_GL);
            */
            cerr << "Not currently supported.\n";
            throw 0;
        }
        else {
            winDataByChr = calcLODWindows(hapDataByChrToCalc, freqDataByChr, mapDataByChr,
                                          GLDataByChrToCalc,
                                          centro, winsizeQuery,
                                          error, MAX_GAP, USE_GL);
        }
        DoubleData *rawWinData = convertWinData2DoubleData(winDataByChr, (THIN ? winsizeQuery : 1));
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
                                 bool WEIGHTED, vector< GenoFreqData * > *genoFreqDataByChr, bool PHASED, bool THIN)
{
    double AUTO_WINSIZE_THRESHOLD = 0.5;
    vector< WinData * > *winDataByChr = NULL;
    vector< HapData * > *hapDataByChrToCalc = NULL;
    vector< GenoLikeData * > *GLDataByChrToCalc = NULL;
    IndData *indDataToCalc = NULL;
    KDEResult *selectedKDEResult = NULL;

    //subset of individuals as given by --kde-subsample
    if (KDE_SUBSAMPLE > 0) subsetData(hapDataByChr, GLDataByChr, indData, &hapDataByChrToCalc, &GLDataByChrToCalc, &indDataToCalc, KDE_SUBSAMPLE, USE_GL, PHASED);
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
            /*
            winDataByChr = calcwLODWindows(hapDataByChrToCalc, freqDataByChr, mapDataByChr,
                                           GLDataByChrToCalc,
                                           genoFreqDataByChr,
                                           centro, multiWinsizes->at(i),
                                           error, MAX_GAP, USE_GL);
            */
            cerr << "Not currently supported.\n";
            throw 0;
        }
        else {
            winDataByChr = calcLODWindows(hapDataByChrToCalc, freqDataByChr, mapDataByChr,
                                          GLDataByChrToCalc,
                                          centro, multiWinsizes->at(i),
                                          error, MAX_GAP, USE_GL);
        }
        DoubleData *rawWinData = convertWinData2DoubleData(winDataByChr, (THIN ? multiWinsizes->at(i) : 1));
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

vector<double> selectSizeClasses(ROHLength *rohLength, int NCLUST)
{
    vector<double> bounds;
    size_t *sortIndex;

    int ngaussians = NCLUST;
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
    char sizeClass = 'A';

    for(int i = 0; i < ngaussians; i++){
        LOG.log("Gaussian class", sizeClass, false); 
        LOG.log(" ( mixture, mean, std ) = (", W[sortIndex[i]], false);
        LOG.log(",", Mu[sortIndex[i]], false);
        LOG.log(",", Sigma[sortIndex[i]], false);
        LOG.log(" )");
        sizeClass++;
    }
    //Find boundaries, there are ngaussians-1 of them, but for the moment this is defined to be 2
    //This finds the 'first' root of the difference between two gaussians
    for(int i = 1; i < ngaussians; i++){
        BoundFinder BF(Mu[sortIndex[i-1]],
                       Sigma[sortIndex[i-1]],
                       W[sortIndex[i-1]],
                       Mu[sortIndex[i]],
                       Sigma[sortIndex[i]],
                       W[sortIndex[i]],
                       1000, 1e-4, false);
        bounds.push_back(BF.findBoundary());
    }

    delete [] W;
    delete [] Mu;
    delete [] Sigma;
    delete [] sortIndex;
    return bounds;
}




