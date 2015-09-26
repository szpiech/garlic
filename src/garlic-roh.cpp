#include "garlic-roh.h"

pthread_mutex_t cerr_mutex = PTHREAD_MUTEX_INITIALIZER;

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

        calcLOD(indData, mapData, hapData, freqData, winData, winsize[pop], error, MAX_GAP);

        pthread_mutex_lock(&cerr_mutex);
        cerr << indDataByPop->at(pop)->pop << " chromosome " << mapDataByChr->at(chr)->chr << " LOD windows finished.\n";
        pthread_mutex_unlock(&cerr_mutex);
    }

    return;
}

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

        calcLOD(indData, mapData, hapData, freqData, winData, winsize[pop], error, MAX_GAP);

        pthread_mutex_lock(&cerr_mutex);
        cerr << indDataByPop->at(pop)->pop << " chromosome " << mapDataByChr->at(chr)->chr << " LOD windows finished.\n";
        pthread_mutex_unlock(&cerr_mutex);
    }

    return;
}


void calcLOD(IndData *indData, MapData *mapData,
             HapData *hapData, FreqData *freqData,
             WinData *winData, int winsize, double error,
             int MAX_GAP)
{
    short **data = hapData->data;
    int nloci = hapData->nloci;
    int nind = hapData->nind;
    int *physicalPos = mapData->physicalPos;
    double *geneticPos = mapData->geneticPos;
    string *locusName = mapData->locusName;
    double *freq = freqData->freq;
    int start = 0;
    int stop = mapData->nloci;
    double **win = winData->data;
    //int nmiss = 0;

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
                    if (physicalPos[i] - physicalPos[prevI] > MAX_GAP)
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
                    if (physicalPos[locus + winsize - 1] - physicalPos[locus + winsize - 2] > MAX_GAP)
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
                        if (physicalPos[i] - physicalPos[prevI] > MAX_GAP)
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

vector< vector< WinData * >* > *calcLODWindows(vector< vector< HapData * >* > *hapDataByPopByChr,
        vector< vector< FreqData * >* > *freqDataByPopByChr,
        vector< MapData * > *mapDataByChr,
        vector< IndData * > *indDataByPop,
        int* winsize, double error, int MAX_GAP, int numThreads)
{

    int numChr = mapDataByChr->size();
    int numPop = indDataByPop->size();

    vector< vector< WinData * >* > *winDataByPopByChr = initWinData(mapDataByChr, indDataByPop);

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
    return winDataByPopByChr;
}

vector< vector< WinData * >* > *calcLODWindowsSinglePop(vector< vector< HapData * >* > *hapDataByPopByChr,
        vector< vector< FreqData * >* > *freqDataByPopByChr,
        vector< MapData * > *mapDataByChr,
        vector< IndData * > *indDataByPop,
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


vector< vector< ROHData * >* > *initROHData(vector< IndData * > *indDataByPop)
{
    vector< vector< ROHData * >* > *rohDataByPopByInd = new vector< vector< ROHData * >* >;

    for (int pop = 0; pop < indDataByPop->size(); pop++)
    {
        vector< ROHData * > *rohDataByInd = new vector< ROHData * >;
        for (int ind = 0; ind < indDataByPop->at(pop)->nind; ind++)
        {
            ROHData *rohData = new ROHData;
            rohDataByInd->push_back(rohData);
        }
        rohDataByPopByInd->push_back(rohDataByInd);
    }

    return rohDataByPopByInd;
}

vector< vector< ROHData * >* > *assembleROHWindows(vector< vector< WinData * >* > *winDataByPopByChr,
        vector< MapData * > *mapDataByChr,
        vector< IndData * > *indDataByPop,
        double *lodScoreCutoffByPop,
        vector< ROHLength * > **rohLengthByPop,
        int winSize,
        int MAX_GAP)
{

    //rohLengthByPop = new vector< ROHLength* >;
    vector< vector< ROHData * >* > *rohDataByPopByInd = initROHData(indDataByPop);

    for (int pop = 0; pop < indDataByPop->size(); pop++)
    {
        vector<int> lengths;
        vector< ROHData * > *rohDataByInd = rohDataByPopByInd->at(pop);
        vector< WinData * > *winDataByChr = winDataByPopByChr->at(pop);
        IndData *indData = indDataByPop->at(pop);
        double lodScoreCutoff = lodScoreCutoffByPop[pop];

        for (int ind = 0; ind < indData->nind; ind++)
        {
            ROHData *rohData = rohDataByInd->at(ind);
            rohData->indID = indData->indID[ind];
            //cerr << "Assembling ROH for individual " << rohData->indID << endl;

            for (int chr = 0; chr < winDataByChr->size(); chr++)
            {
                WinData *winData = winDataByChr->at(chr);
                MapData *mapData = mapDataByChr->at(chr);

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
                    //Window being extended and snp is in ROH
                    //else if(winStart > 0 && inWin[w])
                    //  {
                    //    continue;
                    //  }
                    //Window being extended and snp is not in ROH
                    //end the window at w-1
                    //reset winStart to -1
                    else if (inWin[w] && mapData->physicalPos[w] - mapData->physicalPos[w - 1] > MAX_GAP){
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

        ROHLength *rohLength = initROHLength(lengths.size(), indData->pop);
        for (int i = 0; i < lengths.size(); i++)
        {
            rohLength->length[i] = lengths[i];
        }
        (*rohLengthByPop)->push_back(rohLength);

    }

    return rohDataByPopByInd;
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
    for (int pop = 0; pop < rohLengthByPop->size(); pop++)
    {
        releaseROHLength(rohLengthByPop->at(pop));
    }
    delete rohLengthByPop;
    return;
}

void writeROHData(string outfile,
                  vector< vector< ROHData * >* > *rohDataByPopByInd,
                  vector< MapData * > *mapDataByChr,
                  double *shortMedBound,
                  double *medLongBound,
                  map<string, string> &ind2pop,
                  string version)
{

    string rohOutfile = outfile;
    rohOutfile += ".roh.bed";
    ofstream out;
    out.open(rohOutfile.c_str());
    if (out.fail())
    {
        cerr << "ERROR: Failed to open " << rohOutfile << " for writing.\n";
        throw 0;
    }

    for (int pop = 0; pop < rohDataByPopByInd->size(); pop++)
    {
        vector< ROHData * > *rohDataByInd = rohDataByPopByInd->at(pop);
        for (int ind = 0; ind < rohDataByInd->size(); ind++)
        {
            ROHData *rohData = rohDataByInd->at(ind);
            string popStr = ind2pop[rohData->indID];
            out << "track name=\"Ind: " + rohData->indID + " Pop:" + popStr + 
                " ROH\" description=\"Ind: " + rohData->indID + " Pop:" + popStr + 
                " ROH from GARLIC v" + version + "\" visibility=2 itemRgb=\"On\"\n";

            for (int roh = 0; roh < rohData->chr.size(); roh++)
            {
                int size = (rohData->stop[roh] - rohData->start[roh]);
                char sizeClass = 'C';
                string color = "0,76,153";
                if (size < shortMedBound[pop]){
                    sizeClass = 'A';
                    color = "153,0,0";
                }
                else if (size < medLongBound[pop]){
                    sizeClass = 'B';
                    color = "0,153,0";
                }
                string chr = mapDataByChr->at(rohData->chr[roh])->chr;
                if(chr[0] != 'c' && chr[0] != 'C') chr = "chr" + chr;
                out << chr << "\t" << rohData->start[roh] << "\t" << rohData->stop[roh] 
                    << "\t" << sizeClass << "\t" << size << "\t.\t0\t0\t" << color << endl;
            }
        }
    }
    out.close();
    return;
}