#include "garlic-roh.h"
#include <iomanip>
#include <sstream>

static double AUTO_WINSIZE_THRESHOLD_G = 0.50;
static double AUTO_WINSIZE_SLOPE = 8.3235;
static double AUTO_WINSIZE_INTERCEPT = 138.0521;
static int GMM_MAX_ITER = 1000;
static double GMM_TOL = 1e-5;

void setAutoWinsizeThreshold(double t) { AUTO_WINSIZE_THRESHOLD_G = t; }
void setAutoWinsizeCoef(double slope, double intercept) { AUTO_WINSIZE_SLOPE = slope; AUTO_WINSIZE_INTERCEPT = intercept; }
void setGMMParams(int maxIter, double tol) { GMM_MAX_ITER = maxIter; GMM_TOL = tol; }

int selectWinsizeWeighted(double density){
   /*
    Window size = 8.3235*log(SNVdensity)+138.0521
    Overlap (%) = 6.375*log(SNVdensity)+63.888;*/ 
    int size = int(AUTO_WINSIZE_SLOPE*log(density)+AUTO_WINSIZE_INTERCEPT + 0.5);
    return (size >= 10 ? size : 10);
}

bool inGap(pos_t qStart, pos_t qEnd, pos_t targetStart, pos_t targetEnd)
{
    return ( (targetStart <= qStart && targetEnd >= qStart) ||
             (targetStart <= qEnd && targetEnd >= qEnd) ||
             (targetStart >= qStart && targetEnd <= qEnd) );
}

static int LOD_NUM_THREADS = 1;

void setLODThreads(int n) { LOD_NUM_THREADS = (n > 0 ? n : 1); }

struct LOD_work_order_t
{
    MapData *mapData;
    HapData *hapData;
    FreqData *freqData;
    GenoLikeData *GLData;
    WinData *winData;
    int winsize;
    double error;
    int MAX_GAP;
    bool USE_GL;
    int cStart;
    int cEnd;
    const double *lut;
    Bar *bar;
    int indStart;
    int indStop;
};

//Per-individual LOD windows.  Each individual writes only win[ind][*], so
//there is nothing shared between workers except the (mutex-guarded) progress
//bar; `error` is taken by value so the USE_GL path can overwrite it locally.
static void parallelLOD(LOD_work_order_t *p)
{
    Matrix<geno_t> &data = p->hapData->data;
    const int nloci = p->hapData->nloci;
    const pos_t *physicalPos = p->mapData->physicalPos.data();
    const double *freq = p->freqData->freq.data();
    Matrix<double> &win = p->winData->data;
    GenoLikeData *GLData = p->GLData;
    const double *lut = p->lut;
    const int winsize = p->winsize;
    const int MAX_GAP = p->MAX_GAP;
    const bool USE_GL = p->USE_GL;
    const int cStart = p->cStart;
    const int cEnd = p->cEnd;

    int start = 0;
    int stop = p->mapData->nloci;
    if (nloci - stop < winsize) stop = nloci - winsize + 1;

    //lod() depends on the genotype only through {0,1,2,missing}, so with a
    //scalar error rate there are exactly four possible values per locus.
    //LODV reads them from the precomputed table; with --tgls the error rate
    //varies per genotype and the table does not apply.
    #define LODV(I, IND) ( lut ? lut[4 * (I) + ((data[(I)][(IND)] >= 0 && data[(I)][(IND)] <= 2) ? data[(I)][(IND)] : 3)] \
                               : lod(data[(I)][(IND)], freq[(I)], error) )

    for (int ind = p->indStart; ind < p->indStop; ind++)
    {
        double error = p->error;
        advanceBar(*(p->bar), 1);
        for (int locus = start; locus < stop; locus++)
        {
            win[ind][locus] = 0;

            if (locus == start)
            {
                int prevI = locus;
                for (int i = locus; i < locus + winsize; i++)
                {
                    if (physicalPos[i] - physicalPos[prevI] > MAX_GAP ||
                            inGap(physicalPos[prevI], physicalPos[i], cStart, cEnd))
                    {
                        win[ind][locus] = MISSING;
                        locus = prevI;
                        break;
                    }
                    if (USE_GL) error = GLData->data[i][ind];
                    win[ind][locus] += LODV(i, ind);
                    prevI = i;
                }
            }
            else
            {
                if (win[ind][locus - 1] != MISSING)
                {
                    if (physicalPos[locus + winsize - 1] - physicalPos[locus + winsize - 2] > MAX_GAP ||
                            inGap(physicalPos[locus + winsize - 2], physicalPos[locus + winsize - 1], cStart, cEnd))
                    {
                        win[ind][locus] = MISSING;
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
                                              LODV(locus - 1, ind) +
                                              LODV(locus + winsize - 1, ind);
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
                            locus = prevI;
                            break;
                        }
                        if (USE_GL) error = GLData->data[i][ind];
                        win[ind][locus] += LODV(i, ind);
                        prevI = i;
                    }
                }
            }
        }
    }
    #undef LODV
    return;
}

void calcLOD(MapData *mapData,
             HapData *hapData, FreqData *freqData,
             GenoLikeData *GLData,
             WinData *winData, centromere *centro,
             int winsize, double error, int MAX_GAP, bool USE_GL)
{
    const int nloci = hapData->nloci;
    const int nind = hapData->nind;

    Bar bar;
    barInit(bar, hapData->nind, 100);

    //Four LOD values per locus, computed once instead of once per genotype:
    //this replaces nind * nloci log10() calls with 4 * nloci.  Not applicable
    //when --tgls supplies a per-genotype error rate.
    double *lut = NULL;
    vector<double> lutBuf;
    if (!USE_GL)
    {
        lutBuf.assign(4 * size_t(nloci), 0.0);
        lut = lutBuf.data();
        const double *freq = freqData->freq.data();
        for (int i = 0; i < nloci; i++)
        {
            lut[4 * size_t(i) + 0] = lod(0, freq[i], error);
            lut[4 * size_t(i) + 1] = lod(1, freq[i], error);
            lut[4 * size_t(i) + 2] = lod(2, freq[i], error);
            lut[4 * size_t(i) + 3] = lod(-9, freq[i], error);
        }
    }

    LOD_work_order_t proto;
    proto.mapData = mapData; proto.hapData = hapData; proto.freqData = freqData;
    proto.GLData = GLData; proto.winData = winData; proto.winsize = winsize;
    proto.error = error; proto.MAX_GAP = MAX_GAP; proto.USE_GL = USE_GL;
    proto.cStart = centro->centromereStart(mapData->chr);
    proto.cEnd = centro->centromereEnd(mapData->chr);
    proto.lut = lut; proto.bar = &bar;
    proto.indStart = 0; proto.indStop = nind;

    int nt = LOD_NUM_THREADS;
    if (nt > nind) nt = nind;
    if (nt < 1) nt = 1;

    if (nt == 1)
    {
        parallelLOD(&proto);
    }
    else
    {
        //orders is sized and filled before any thread starts; a thread holds
        //&orders[i], so it must not be grown afterwards.
        vector<LOD_work_order_t> orders(nt);
        int per = nind / nt;
        int extra = nind % nt;
        int at = 0;
        for (int i = 0; i < nt; i++)
        {
            orders[i] = proto;
            orders[i].indStart = at;
            at += per + (i < extra ? 1 : 0);
            orders[i].indStop = at;
        }
        vector<std::thread> peer;
        peer.reserve(nt);
        for (int i = 0; i < nt; i++) peer.emplace_back(parallelLOD, &(orders[i]));
        for (int i = 0; i < nt; i++) peer[i].join();
    }

    finalize(bar);
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
    vector<unsigned int> NUM_PER_THREAD = make_thread_partition(numThreads, hapData->nloci);

    Bar bar;
    barInit(bar,hapData->nind,100);

    WLOD_work_order_t *order;
    vector<std::thread> peer;
    peer.reserve(numThreads);
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
        order->bar = &bar;
        previous += NUM_PER_THREAD[i];
        order->stop = previous;
        order->numThreads = numThreads;

        peer.emplace_back(parallelwLOD, order);
        orders.push_back(order);
    }

    for (int i = 0; i < numThreads; i++){
        peer[i].join();
        delete orders[i];
    }

    finalize(bar);

    orders.clear();
    return;
}

//Was launched as (void *(*)(void *))parallelwLOD, a cast between incompatible
//function pointer types -- it returns void, not void *.  Calling through such
//a cast is undefined behaviour; std::thread calls it with its real signature.
void parallelwLOD(WLOD_work_order_t *p){
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
    Bar *bar = p->bar;
    int numThreads = p->numThreads;

    Matrix<geno_t> &data = p->hapData->data;
    int nloci = p->hapData->nloci;
    int nind = p->hapData->nind;
    pos_t *physicalPos = p->mapData->physicalPos.data();
    double *geneticPos = p->mapData->geneticPos.data();
    double *freq = p->freqData->freq.data();
    Matrix<double> &win = p->winData->data;

    int cStart = p->cStart;
    int cEnd = p->cEnd;

    //Check if the last window would overshoot the last locus in the data
    if (nloci - stop < winsize) stop = nloci - winsize + 1;

    int size = stop-start+winsize+1;
    double *score;
    vector<double> scoreBuf;

    //cerr << "XXX " << start << " " << stop << endl;

  //  cerr << "wLOD " << p->mapData->chr << endl;

    //For each individual
    for (int ind = 0; ind < nind; ind++){
        advanceBar(*bar,1.0/double(numThreads));
        scoreBuf.resize(size);
        score = scoreBuf.data();
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

    }
}

vector< WinData * > *calcLODWindows(vector< HapData * > *hapDataByChr,
                                    vector< FreqData * > *freqDataByChr,
                                    vector< MapData * > *mapDataByChr,
                                    vector< GenoLikeData * > *GLDataByChr,
                                    centromere *centro,
                                    int winsize, double error, int MAX_GAP, bool USE_GL,
                                    int sexWinsize, const vector<ChrRole> *role)
{
    if (!LOG.isQuiet())
    {
        //Named separately when they differ: this line is the record of what
        //the windows in this run are, and one number would not be it.
        cerr << "Calculating LOD scores with winsize " << winsize;
        if (sexWinsize > 0 && role != NULL)
        {
            for (unsigned int c = 0; c < role->size(); c++)
                if (role->at(c) == CHR_SEX_SHARED)
                { cerr << " (" << sexWinsize << " on the shared sex chromosome)"; break; }
        }
        cerr << ".\n";
    }

    vector< WinData * > *winDataByChr = initWinData(mapDataByChr, hapDataByChr->at(0)->nind);

    for (unsigned int chr = 0; chr < winDataByChr->size(); chr++)
    {
        if (progressEnabled()) cerr << mapDataByChr->at(chr)->chr << "    ";
        const int w = winsizeForChr(winsize, sexWinsize, role, chr);
        if(USE_GL){
            calcLOD(mapDataByChr->at(chr),
                    hapDataByChr->at(chr), freqDataByChr->at(chr),
                    GLDataByChr->at(chr),
                    winDataByChr->at(chr), centro,
                    w, error, MAX_GAP, USE_GL);
        }
        else{
            calcLOD(mapDataByChr->at(chr),
                    hapDataByChr->at(chr), freqDataByChr->at(chr),
                    NULL,
                    winDataByChr->at(chr), centro,
                    w, error, MAX_GAP, USE_GL);
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
                                     int M, double mu, int numThreads,
                                     int sexWinsize, const vector<ChrRole> *role)
{
    if (!LOG.isQuiet())
    {
        //Named separately when they differ: this line is the record of what
        //the windows in this run are, and one number would not be it.
        cerr << "Calculating LOD scores with winsize " << winsize;
        if (sexWinsize > 0 && role != NULL)
        {
            for (unsigned int c = 0; c < role->size(); c++)
                if (role->at(c) == CHR_SEX_SHARED)
                { cerr << " (" << sexWinsize << " on the shared sex chromosome)"; break; }
        }
        cerr << ".\n";
    }

    vector< WinData * > *winDataByChr = initWinData(mapDataByChr, hapDataByChr->at(0)->nind);

    for (unsigned int chr = 0; chr < winDataByChr->size(); chr++)
    {
        if (progressEnabled()) cerr << mapDataByChr->at(chr)->chr << "    ";
        const int w = winsizeForChr(winsize, sexWinsize, role, chr);
        if(USE_GL){
            calcwLOD(mapDataByChr->at(chr),
                     hapDataByChr->at(chr),
                     freqDataByChr->at(chr),
                     GLDataByChr->at(chr),
                     ldDataByChr->at(chr),
                     winDataByChr->at(chr), centro,
                     w, error, MAX_GAP, USE_GL, mu, M, numThreads);
        }
        else{
            calcwLOD(mapDataByChr->at(chr),
                     hapDataByChr->at(chr),
                     freqDataByChr->at(chr),
                     NULL,
                     ldDataByChr->at(chr),
                     winDataByChr->at(chr), centro,
                     w, error, MAX_GAP, USE_GL, mu, M, numThreads);   
        }
    }
    return winDataByChr;
}



/*
 * Genotype is 0/1/2 counting the number of alternate alleles
 *
 */
double lod(const int genotype, const double &freq, const double &error)
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

struct ROH_work_order_t
{
    vector< WinData * > *winDataByChr;
    vector< MapData * > *mapDataByChr;
    IndData *indData;
    centromere *centro;
    double lodScoreCutoff;
    int winSize;
    int MAX_GAP;
    double OVERLAP_FRAC;
    bool CM;
    const vector<ChrRole> *role;
    double sexCutoff;
    bool haveSexCutoff;
    int sexWinsize;
    vector< ROHData * > *rohDataByInd;
    //Per-thread so there is no shared push_back; concatenated in thread order
    //below, which reproduces the serial order exactly because each thread owns
    //a contiguous block of individuals.
    vector<double> lengths;
    int indStart;
    int indStop;
};

static void parallelAssembleROH(ROH_work_order_t *p)
{
    vector< WinData * > *winDataByChr = p->winDataByChr;
    vector< MapData * > *mapDataByChr = p->mapDataByChr;
    IndData *indData = p->indData;
    centromere *centro = p->centro;
    const int winSizeRun = p->winSize;
    const int MAX_GAP = p->MAX_GAP;
    const double OVERLAP_FRAC = p->OVERLAP_FRAC;
    const bool CM = p->CM;
    const vector<ChrRole> *role = p->role;
    vector< ROHData * > *rohDataByInd = p->rohDataByInd;
    vector<double> &lengths = p->lengths;

    for (int ind = p->indStart; ind < p->indStop; ind++)
    {
        ROHData *rohData = rohDataByInd->at(ind);
        rohData->indID = indData->indID[ind];

        for (unsigned int chr = 0; chr < winDataByChr->size(); chr++)
        {
            //Every tract this chromosome yields is written; only an autosome's
            //feeds the size-class GMM.  With role == NULL -- every caller
            //before roles existed, and every run with no sex chromosome --
            //this is true everywhere and the behaviour is unchanged.
            const bool sizeClassChr = (role == NULL || chr >= role->size() ||
                                       role->at(chr) == CHR_AUTOSOME);

            //Both of these belong to the chromosome, not to the run: a cutoff
            //is a statement about windows of a given size, so a chromosome
            //given its own size needs its own cutoff to go with it.
            const int winSize = winsizeForChr(winSizeRun, p->sexWinsize, role, chr);
            double OVERLAP_THRESHOLD = OVERLAP_FRAC * winSize;
            OVERLAP_THRESHOLD = (OVERLAP_THRESHOLD >= 1) ? OVERLAP_THRESHOLD : 1;
            OVERLAP_THRESHOLD = (OVERLAP_THRESHOLD <= winSize) ? OVERLAP_THRESHOLD : winSize;
            const double lodScoreCutoff =
                (p->haveSexCutoff && role != NULL && chr < role->size() &&
                 role->at(chr) == CHR_SEX_SHARED) ? p->sexCutoff : p->lodScoreCutoff;

            WinData *winData = winDataByChr->at(chr);
            MapData *mapData = mapDataByChr->at(chr);
            pos_t *pos;
            double *gpos;
            gpos = mapData->geneticPos.data();
            pos = mapData->physicalPos.data();

            pos_t cStart = centro->centromereStart(mapData->chr);
            pos_t cEnd = centro->centromereEnd(mapData->chr);

            //translation of the perl script here###Updated to match trevor's algorithm
            //int winStart = -1;
            //int winStop = -1;
            //Difference array + prefix sum instead of incrementing winSize
            //entries per passing window: O(1) per window rather than
            //O(winSize), and it removes an out-of-bounds write.
            //
            //winData->nloci == mapData->nloci (initWinData sizes it from the
            //map), so the old inner loop reached inWin[nloci + winSize - 2],
            //up to winSize-1 elements past the end.  It was masked only
            //because the tail windows hold MISSING (-9999) and so failed the
            //cutoff test -- passing --lod-cutoff below -9999 wrote past the end.
            vector<int> inWinDiff(mapData->nloci + 1);
            for (int w = 0; w <= mapData->nloci; w++) inWinDiff[w] = 0;
            for (int w = 0; w < winData->nloci; w++)
            {
                if (winData->data[ind][w] >= lodScoreCutoff)
                {
                    int lo = w;
                    int hi = w + winSize;
                    if (hi > mapData->nloci) hi = mapData->nloci;
                    if (lo < mapData->nloci) { inWinDiff[lo]++; inWinDiff[hi]--; }
                }
            }
            vector<short> inWin(mapData->nloci);
            int running = 0;
            for (int w = 0; w < mapData->nloci; w++)
            {
                running += inWinDiff[w];
                inWin[w] = short(running);
            }

            double gwinStart = -1;
            double gwinStop = -1;
            pos_t winStart = -1;
            int winStartIndex = -1;
            pos_t winStop = -1;
            int winStopIndex = -1;
            for (int w = 0; w < mapData->nloci; w++)
            {
                //No window being extended and the snp is in ROH
                //Start the window
                if (winStart < 0 && inWin[w] >= OVERLAP_THRESHOLD)
                {
                    gwinStart = gpos[w];
                    winStart = pos[w];
                    winStartIndex = w;
                }
                else if (inWin[w] >= OVERLAP_THRESHOLD && (mapData->physicalPos[w] - mapData->physicalPos[w - 1] > MAX_GAP ||
                         inGap(mapData->physicalPos[w - 1], mapData->physicalPos[w], cStart, cEnd)) ) {
                    gwinStop = gpos[w - 1];
                    winStop = pos[w - 1];
                    winStopIndex = w - 1;
                    if(winStopIndex - winStartIndex + 1 >= OVERLAP_THRESHOLD){
                        double size = CM ? gwinStop - gwinStart : winStop - winStart + 1;
                        if (sizeClassChr) lengths.push_back(size);
                        rohData->length.push_back(size);
                        rohData->chr.push_back(chr);
                        rohData->start.push_back(winStart);
                        rohData->stop.push_back(winStop);
                    }
                    gwinStop = -1;
                    winStop = -1;
                    winStopIndex = -1;
                    gwinStart = gpos[w];
                    winStart = pos[w];
                    winStartIndex = w;
                }
                else if (winStart >= 0 && ! (inWin[w] >= OVERLAP_THRESHOLD) )
                {
                    gwinStop = gpos[w - 1];
                    winStop = pos[w - 1];
                    winStopIndex = w - 1;
                    if(winStopIndex - winStartIndex + 1 >= OVERLAP_THRESHOLD){
                        double size = CM ? gwinStop - gwinStart : winStop - winStart + 1;
                        if (sizeClassChr) lengths.push_back(size);
                        rohData->length.push_back(size);
                        rohData->chr.push_back(chr);
                        rohData->start.push_back(winStart);
                        rohData->stop.push_back(winStop);
                    }
                    gwinStart = -1;
                    winStart = -1;
                    winStartIndex = -1;
                    gwinStop = -1;
                    winStop = -1;
                    winStopIndex = -1;
                }
                else if (winStart > 0 && w + 1 >= mapData->nloci)
                {
                    gwinStop = gpos[w];
                    winStop = pos[w];
                    winStopIndex = w;
                    if(winStopIndex - winStartIndex + 1 >= OVERLAP_THRESHOLD){
                        double size = CM ? gwinStop - gwinStart : winStop - winStart + 1;
                        if (sizeClassChr) lengths.push_back(size);
                        rohData->length.push_back(size);
                        rohData->chr.push_back(chr);
                        rohData->start.push_back(winStart);
                        rohData->stop.push_back(winStop);
                    }
                    gwinStart = -1;
                    winStart = -1;
                    winStartIndex = -1;
                    gwinStop = -1;
                    winStop = -1;
                    winStopIndex = -1;
                }
            }

        }
    }
    return;
}

vector< ROHData * > *assembleROHWindows(vector< WinData * > *winDataByChr,
                                        vector< MapData * > *mapDataByChr,
                                        IndData *indData,
                                        centromere *centro,
                                        double lodScoreCutoff,
                                        ROHLength **rohLength,
                                        int winSize,
                                        int MAX_GAP,
                                        double OVERLAP_FRAC, bool CM,
                                        const vector<ChrRole> *role,
                                        double sexCutoff, bool haveSexCutoff,
                                        int sexWinsize)
{
    vector< ROHData * > *rohDataByInd = initROHData(indData);

    //Individuals are independent: each writes only its own ROHData.
    int nt = LOD_NUM_THREADS;
    if (nt > indData->nind) nt = indData->nind;
    if (nt < 1) nt = 1;

    vector<ROH_work_order_t> orders(nt);
    int per = indData->nind / nt;
    int extra = indData->nind % nt;
    int at = 0;
    for (int i = 0; i < nt; i++)
    {
        orders[i].winDataByChr = winDataByChr;
        orders[i].mapDataByChr = mapDataByChr;
        orders[i].indData = indData;
        orders[i].centro = centro;
        orders[i].lodScoreCutoff = lodScoreCutoff;
        orders[i].winSize = winSize;
        orders[i].MAX_GAP = MAX_GAP;
        orders[i].OVERLAP_FRAC = OVERLAP_FRAC;
        orders[i].sexCutoff = sexCutoff;
        orders[i].haveSexCutoff = haveSexCutoff;
        orders[i].sexWinsize = sexWinsize;
        orders[i].CM = CM;
        orders[i].role = role;
        orders[i].rohDataByInd = rohDataByInd;
        orders[i].indStart = at;
        at += per + (i < extra ? 1 : 0);
        orders[i].indStop = at;
    }

    if (nt == 1)
    {
        parallelAssembleROH(&(orders[0]));
    }
    else
    {
        vector<std::thread> peer;
        peer.reserve(nt);
        for (int i = 0; i < nt; i++) peer.emplace_back(parallelAssembleROH, &(orders[i]));
        for (int i = 0; i < nt; i++) peer[i].join();
    }

    vector<double> lengths;
    for (int i = 0; i < nt; i++)
        lengths.insert(lengths.end(), orders[i].lengths.begin(), orders[i].lengths.end());

    ROHLength *rohLengths = initROHLength(lengths.size());
    for (unsigned int i = 0; i < lengths.size(); i++)
    {
        rohLengths->length[i] = lengths[i];
    }
    (*rohLength) = rohLengths;

    return rohDataByInd;
}

ROHLength *initROHLength(int size)
{
    ROHLength *rohLength = new ROHLength;
    //rohLength->pop = pop;
    rohLength->length.resize(size);
    rohLength->size = size;
    return rohLength;
}

void releaseROHLength(ROHLength *rohLength)
{
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

//Default ostream precision is 6 significant digits, which turns a 67,657,700 bp
//total into "6.76577e+07".  Physical lengths are whole base pairs; genetic
//lengths need decimals.
static string fmtLength(double v, bool CM)
{
    ostringstream ss;
    if (CM) ss << fixed << setprecision(6) << v;
    else    ss << (long long)(v + 0.5);
    return ss.str();
}

//Spreadsheet-style labels so more than 26 size classes stay distinguishable:
//A..Z, then AA, AB, ...  The letter used to be incremented past 'Z' into
//punctuation.
string sizeClassLabel(int k)
{
    string s;
    int n = k;
    do { s.insert(s.begin(), char('A' + (n % 26))); n = n / 26 - 1; } while (n >= 0);
    return s;
}

//The nine ColorBrewer entries are kept for the first nine classes so existing
//output is unchanged; beyond that, walk the hue circle.
vector<string> makeClassColors(int nclass)
{
    const char *base[9] = {"228,26,28", "77,175,74", "55,126,184", "152,78,163",
                           "255,127,0", "255,255,51", "166,86,40", "247,129,191",
                           "153,153,153"};
    vector<string> out;
    for (int i = 0; i < nclass; i++)
    {
        if (i < 9) { out.push_back(string(base[i])); continue; }
        double h = fmod(0.13 + 0.618033988749895 * double(i - 9), 1.0) * 6.0;
        int seg = int(h);
        double frac = h - double(seg);
        double v = 210, p = 60;
        double q = v - (v - p) * frac, t = p + (v - p) * frac;
        double rr, gg, bb;
        if      (seg == 0) { rr = v; gg = t; bb = p; }
        else if (seg == 1) { rr = q; gg = v; bb = p; }
        else if (seg == 2) { rr = p; gg = v; bb = t; }
        else if (seg == 3) { rr = p; gg = q; bb = v; }
        else if (seg == 4) { rr = t; gg = p; bb = v; }
        else               { rr = v; gg = p; bb = q; }
        ostringstream ss;
        ss << int(rr) << "," << int(gg) << "," << int(bb);
        out.push_back(ss.str());
    }
    return out;
}

//Class index for a ROH of this size, using the same rule as writeROHData.
static int rohSizeClassIndex(double size, vector<double> &bounds)
{
    for (unsigned int i = 0; i < bounds.size(); i++)
        if (size < bounds[i]) return int(i);
    return int(bounds.size());
}

void writeFROH(string outfile,
               vector< ROHData * > *rohDataByInd,
               vector< MapData * > *mapDataByChr,
               vector< double > bounds,
               IndData *indData,
               centromere *centro,
               bool CM,
               const string &popLabel,
               bool pooled,
               const vector<ChrRole> *role,
               const ExcludedRegions *excluded,
               int denomMode,
               const ChrLengths *chrLengths)
{
    const int nclass = int(bounds.size()) + 1;   //A..  plus the ALL column set

    //Two denominators, not one.  An individual that cannot carry a run of
    //homozygosity on the sex chromosome must not have that chromosome in the
    //denominator of anything -- it would deflate its FROH by the sex
    //chromosome's share of the analysed span -- and autosomal and
    //sex-chromosomal FROH are not summable in any case: the X coalesces
    //faster than the autosomes and carries more ROH at the same level of
    //consanguinity, so the two are quantities to compare, not to pool.
    //
    //Because eligibility is per individual but the SPAN is not, this is two
    //constants and a rule about who gets the second one, rather than a
    //denominator per individual.
    double denomAuto = 0, denomSex = 0;
    bool haveSexChr = false, haveExcluded = false;
    string sexChrNames, excludedList;
    for (unsigned int chr = 0; chr < mapDataByChr->size(); chr++)
    {
        MapData *md = mapDataByChr->at(chr);
        ChrRole r = (role != NULL && chr < role->size()) ? role->at(chr) : CHR_AUTOSOME;
        if (md->nloci < 2) continue;

        pos_t first = md->physicalPos[0];
        pos_t last  = md->physicalPos[md->nloci - 1];
        const vector<Interval> *iv = (excluded != NULL) ? excluded->get(md->chr) : NULL;

        double span;
        if (denomMode == FROH_ASSEMBLY)
        {
            //The RAW length of the chromosome.  Nothing is subtracted -- not
            //the centromere, not an excluded region, not the part beyond the
            //outermost marker -- because the point of this mode is a number
            //reproducible from an external table and comparable with a
            //published one, and "the assembly minus whatever this run
            //happened to exclude" is neither.
            //
            //main has already refused the run if any analysed chromosome has
            //no length here, so the lookup cannot come back empty.
            span = double(chrLengths->get(md->chr));
        }
        else if (CM)
        {
            span = md->geneticPos[md->nloci - 1] - md->geneticPos[0];
            //The centromere is NOT subtracted in cM: a genetic map already
            //spans about 0 cM across it, and taking it off again would
            //double-count.  That reasoning does not carry to a
            //pseudoautosomal region, which carries an obligate crossover and
            //is the most recombinationally active part of the genome per base
            //-- so an interior one overstates a cM denominator by more than a
            //centromere would.  Subtracted here as the genetic distance
            //between the loci that flank it, which needs no interpolation
            //because both of them are in the map.
            for (unsigned int k = 0; iv != NULL && k < iv->size(); k++)
            {
                if (iv->at(k).start <= first || iv->at(k).end >= last) continue;
                int after = -1;
                for (int locus = 0; locus < md->nloci; locus++)
                    if (md->physicalPos[locus] > iv->at(k).end) { after = locus; break; }
                if (after > 0) span -= (md->geneticPos[after] - md->geneticPos[after - 1]);
            }
        }
        else
        {
            double lo = md->physicalPos[0];
            double hi = md->physicalPos[md->nloci - 1];
            span = hi - lo;
            double gs = centro->centromereStart(md->chr);
            double ge = centro->centromereEnd(md->chr);
            if (ge > gs)
            {
                double ols = (gs > lo ? gs : lo);
                double ole = (ge < hi ? ge : hi);
                if (ole > ols) span -= (ole - ols);
            }
            //Merged, so overlapping exclusions are not subtracted twice, and
            //clipped, so a terminal one contributes nothing.
            if (iv != NULL) span -= double(excluded->overlap(md->chr, first, last));
        }

        for (unsigned int k = 0; iv != NULL && k < iv->size(); k++)
        {
            haveExcluded = true;
            if (!excludedList.empty()) excludedList += ",";
            ostringstream ss;
            ss << md->chr << ":" << (long long)iv->at(k).start << "-" << (long long)iv->at(k).end;
            excludedList += ss.str();
        }

        if (r == CHR_SEX_SHARED)
        {
            denomSex += span;
            haveSexChr = true;
            if (!sexChrNames.empty()) sexChrNames += ",";
            sexChrNames += md->chr;
        }
        else denomAuto += span;
    }

    ofstream out;
    //ios::binary so the bytes written are the bytes on disk.  Without it the
//Windows CRT translates every \n into \r\n, and garlic's text outputs came
//out CRLF there while being otherwise byte-identical -- confirmed by taking a
//Unix .roh.bed, converting it, and reproducing the exact md5 CI reported.
//A BED file should be LF, and a result that differs by platform is a result
//users cannot checksum against each other.  No effect on POSIX.
    out.open(outfile.c_str(), ios::binary);
    if (out.fail())
    {
        LOG.err("ERROR: Failed to open", outfile);
        throw 0;
    }

    out << "## garlic FROH\n";
    //Whether the individuals below were analysed TOGETHER.  Only the pooled
    //case needs saying: in a per-population table the pop column is already
    //constant and names the population, so a "## population" line would
    //repeat it -- and would make a population's table differ depending on
    //whether it was run alongside others, which is exactly the equivalence
    //the suite asserts.
    if (pooled) out << "## populations_pooled\ttrue\n";
    out << "## units\t" << (CM ? "cM" : "bp") << "\n";
    out << "## size_class_boundaries";
    for (unsigned int i = 0; i < bounds.size(); i++) out << "\t" << bounds[i];
    out << "\n";
    if (haveExcluded)
        out << "## excluded_regions\t" << excludedList << "\n";
    //Named, because the two conventions differ by several percent and a
    //denominator is not self-describing.  The source goes with it: "assembly"
    //is not a number until you know which assembly.
    out << "## denominator_mode\t"
        << (denomMode == FROH_ASSEMBLY ? "assembly" : "analyzed") << "\n";
    if (denomMode == FROH_ASSEMBLY && chrLengths != NULL)
        out << "## chromosome_lengths\t" << chrLengths->source() << "\n";
    out << "## autosome_denominator\t" << (long long)(denomAuto + 0.5) << "\t";
    out << (denomMode == FROH_ASSEMBLY
               ? "sum over analysed autosomes of the full chromosome length, nothing subtracted"
               : (CM ? "sum over analysed autosomes of (last - first genetic position)"
                     : "sum over analysed autosomes of (last - first physical position), minus the assembly gap inside that span"))
        << "\n";
    //True in both modes and easy to misread in one: a chromosome carrying a
    //single locus is not analysed and is in no denominator, so under assembly
    //this is not the whole assembly either.
    out << "## note\ta chromosome with fewer than two analysed loci is in neither denominator\n";
    if (haveSexChr)
    {
        out << "## sex_chromosomes\t" << sexChrNames << "\n";
        out << "## sexchr_denominator\t" << (long long)(denomSex + 0.5) << "\t"
            << "the same over the shared sex chromosome; applies only to individuals diploid there\n";
        if (denomMode == FROH_ASSEMBLY && haveExcluded)
            out << "## sexchr_denominator_note\tincludes the excluded regions above, which assembly does not subtract\n";
        out << "## na\tNA where an individual cannot carry a run of homozygosity on the sex chromosome\n";
    }

    //One row per individual, metrics across columns: a row per size class was
    //a shape nothing reads without reshaping it first, and it cannot carry
    //two regions at all without repeating every individual four times.
    out << "ind\tpop\tzygosity";
    for (int pass = 0; pass < (haveSexChr ? 2 : 1); pass++)
    {
        string pre = (pass == 0 ? "auto_" : "sexchr_");
        for (int k = 0; k <= nclass; k++)
        {
            string cls = (k < nclass ? sizeClassLabel(k) : string("ALL"));
            out << "\t" << pre << cls << "_n"
                << "\t" << pre << cls << "_len"
                << "\t" << pre << cls << "_froh";
        }
    }
    out << "\n";

    for (unsigned int ind = 0; ind < rohDataByInd->size(); ind++)
    {
        ROHData *rohData = rohDataByInd->at(ind);
        //[region][class]; region 0 autosomes, 1 the shared sex chromosome.
        vector< vector<int> >    n(2, vector<int>(nclass + 1, 0));
        vector< vector<double> > tot(2, vector<double>(nclass + 1, 0.0));

        for (unsigned int roh = 0; roh < rohData->length.size(); roh++)
        {
            double size = rohData->length[roh];
            int c = rohData->chr[roh];
            ChrRole r = (role != NULL && c >= 0 && (unsigned int)c < role->size())
                        ? role->at(c) : CHR_AUTOSOME;
            int region = (r == CHR_SEX_SHARED) ? 1 : 0;
            int k = rohSizeClassIndex(size, bounds);
            n[region][k]++;            tot[region][k] += size;
            n[region][nclass]++;       tot[region][nclass] += size;
        }

        int zyg = (ind < indData->zygo.size()) ? indData->zygo[ind] : ZYG_UNKNOWN;
        const char *zname = (zyg == ZYG_HOMOGAMETIC ? "homogametic"
                          : (zyg == ZYG_HETEROGAMETIC ? "heterogametic" : "unknown"));

        out << rohData->indID << "\t" << indData->pop[ind] << "\t"
            << (haveSexChr ? zname : "NA");

        for (int region = 0; region < (haveSexChr ? 2 : 1); region++)
        {
            //NA rather than 0 for an individual that cannot be autozygous
            //there: zero is a measurement and this is not one.
            bool eligible = (region == 0) || eligibleForCalling(CHR_SEX_SHARED, zyg);
            double denom = (region == 0) ? denomAuto : denomSex;
            for (int k = 0; k <= nclass; k++)
            {
                if (!eligible) { out << "\tNA\tNA\tNA"; continue; }
                out << "\t" << n[region][k]
                    << "\t" << fmtLength(tot[region][k], CM)
                    << "\t" << fixed << setprecision(8)
                    << (denom > 0 ? tot[region][k] / denom : 0.0) << defaultfloat;
            }
        }
        out << "\n";
    }

    out.close();
    LOG.log("FROH table:", outfile);
    return;
}

void writeROHData(string outfile,
                  vector< ROHData * > *rohDataByInd,
                  vector< MapData * > *mapDataByChr,
                  vector< double > bounds,
                  const vector<string> &pop,
                  string version, bool CM)
{
    //--nclust is unbounded, but the palette had nine entries and the index was
    //clamped at 8, so every class past the ninth was drawn in the same grey.
    const int nclass = int(bounds.size()) + 1;
    vector<string> colors = makeClassColors(nclass);
   
    ofstream out;
    //ios::binary so the bytes written are the bytes on disk.  Without it the
//Windows CRT translates every \n into \r\n, and garlic's text outputs came
//out CRLF there while being otherwise byte-identical -- confirmed by taking a
//Unix .roh.bed, converting it, and reproducing the exact md5 CI reported.
//A BED file should be LF, and a result that differs by platform is a result
//users cannot checksum against each other.  No effect on POSIX.
    out.open(outfile.c_str(), ios::binary);
    if (out.fail())
    {
        LOG.err("ERROR: Failed to open", outfile);
        throw 0;
    }

    string popName;

    for (unsigned int ind = 0; ind < rohDataByInd->size(); ind++)
    {
        ROHData *rohData = rohDataByInd->at(ind);
        popName = pop[ind];
        out << "track name=\"Ind: " + rohData->indID + " Pop:" + popName +
            " ROH\" description=\"Ind: " + rohData->indID + " Pop:" + popName +
            " ROH from GARLIC v" + version + "\" visibility=2 itemRgb=\"On\"\n";

        for (unsigned int roh = 0; roh < rohData->chr.size(); roh++)
        {
            double size = rohData->length[roh];
            char sc = 'A';
            char sizeClass = 'X';
            string color = "X";

            int k = rohSizeClassIndex(size, bounds);
            string sizeClassLab = sizeClassLabel(k);
            color = colors[k];
            (void)sc; (void)sizeClass;

            string chr = mapDataByChr->at(rohData->chr[roh])->chr;
            if (chr[0] != 'c' && chr[0] != 'C') chr = "chr" + chr;
            //BED chromStart is 0-based and chromEnd is exclusive, but
            //rohData->start is the 1-based physical position of the ROH's
            //first variant, and it used to be written verbatim.  So
            //chromEnd - chromStart came out as length-1 for every tract and a
            //browser drew each ROH one base short and one base right, while
            //column 5 carried the true length -- the file disagreed with
            //itself about its own intervals.  Emitting start-1 makes
            //chromEnd - chromStart == length.
            //
            //start and stop are physical positions in both bp and --cm mode
            //(only column 5 changes to a genetic length), so the shift applies
            //to both.  Clamped at 0 in case an input carries position 0, which
            //has no 0-based representation.
            pos_t bedStart = rohData->start[roh] - 1;
            if (bedStart < 0) bedStart = 0;
            if(CM){
                out << chr << "\t" << bedStart << "\t" << rohData->stop[roh]
                    << "\t" << sizeClassLab << "\t" << size << "\t.\t0\t0\t" << color << endl;
            }
            else{
                out << chr << "\t" << bedStart << "\t" << rohData->stop[roh]
                    << "\t" << sizeClassLab << "\t" << int(size) << "\t.\t0\t0\t" << color << endl;
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

double selectLODCutoff(KDEResult *kdeResult, int wsize, bool &ok)
{
    double LOD_CUTOFF;
    ok = true;
    try { LOD_CUTOFF = get_min_btw_modes(kdeResult->x.data(), kdeResult->y.data(), kdeResult->size, wsize); }
    catch (...)
    {
        logCurrentException("locating the minimum between LOD score modes");
        LOG.err("ERROR: Failed to find the minimum between modes in the LOD score density.");
        LOG.err("\tResults from density estimation have been written to file for inspection.");
        LOG.err("\tA cutoff can be manually specified on the command line with", ARG_LOD_CUTOFF);
        ok = false;
        return -1;
    }
    return LOD_CUTOFF;
}


long long maskIneligibleWindows(vector< WinData * > *winDataByChr, IndData *indData,
                                const vector<ChrRole> *role)
{
    if (role == NULL) return 0;
    long long n = 0;
    for (unsigned int chr = 0; chr < winDataByChr->size(); chr++)
    {
        if (chr >= role->size() || role->at(chr) == CHR_AUTOSOME) continue;
        WinData *win = winDataByChr->at(chr);
        for (int ind = 0; ind < win->nind; ind++)
        {
            if (eligibleForCalling(role->at(chr), indData->zygo[ind])) continue;
            for (int locus = 0; locus < win->nloci; locus++)
            {
                if (win->data[ind][locus] == MISSING) continue;
                win->data[ind][locus] = MISSING;
                n++;
            }
        }
    }
    return n;
}

bool reportSexChrLODCutoff(vector< WinData * > *winDataByChr, IndData *indData,
                           const vector<ChrRole> *role, int step, int wsize,
                           double appliedCutoff, bool cutoffIsItsOwn)
{
    if (role == NULL) return false;
    bool haveSexChr = false;
    for (unsigned int chr = 0; chr < role->size(); chr++)
        if (role->at(chr) == CHR_SEX_SHARED) haveSexChr = true;
    if (!haveSexChr) return false;

    //Subsampling is deliberately not applied: there is little enough of this
    //chromosome as it is, and this costs one KDE on a fraction of the genome.
    DoubleData *raw = convertWinData2DoubleData(winDataByChr, step, role, CHR_SEX_SHARED);
    if (raw == NULL || raw->size < 2) { if (raw) releaseDoubleData(raw); return false; }

    bool found = false;
    double cutoff = 0;
    try
    {
        KDEResult *kde = computeKDE(raw->data.data(), raw->size);
        try { cutoff = get_min_btw_modes(kde->x.data(), kde->y.data(), kde->size, wsize); found = true; }
        catch (...) { found = false; }
        releaseKDEResult(kde);
    }
    catch (...) { found = false; }
    releaseDoubleData(raw);

    if (!found)
    {
        //Expected, not a failure: one chromosome's windows from part of the
        //cohort often have no second mode to find.  That is the reason the
        //autosomal cutoff is the one used.
        LOG.log("Sex chromosome: no separate LOD score cutoff could be estimated from its own");
        if (cutoffIsItsOwn) LOG.log("\twindows; the one given for it is the one applied.");
        else LOG.log("\twindows; the autosomal cutoff is the one applied.");
        return false;
    }

    LOG.log("Sex chromosome: its own windows would give a LOD score cutoff of", cutoff, false);
    if (cutoffIsItsOwn) LOG.log("; the cutoff given for it,", appliedCutoff, false);
    else LOG.log("; the autosomal cutoff", appliedCutoff, false);
    LOG.log(" is the one applied.");
    return true;
}

double selectLODCutoff(vector< WinData * > *winDataByChr, IndData *indData, int KDE_SUBSAMPLE, string kdeoutfile, int step, int wsize, bool &ok,
                       const vector<ChrRole> *role, ChrRole keep)
{
    ok = true;
    //Format the LOD window data into a single array per pop with no missing data
    //Prepped for KDE
    DoubleData *rawWinData;
    double LOD_CUTOFF;

    if (KDE_SUBSAMPLE <= 0) rawWinData = convertWinData2DoubleData(winDataByChr, step, role, keep);
    else rawWinData = convertSubsetWinData2DoubleData(winDataByChr, indData, KDE_SUBSAMPLE, step, role, keep);

    //Compute KDE of LOD score distribution
    if (!LOG.isQuiet()) cerr << "Estimating distribution of raw LOD score windows:\n";
    KDEResult *kdeResult = computeKDE(rawWinData->data.data(), rawWinData->size);
    releaseDoubleData(rawWinData);

    //Output kde points
    try { writeKDEResult(kdeResult, kdeoutfile); }
    catch (...) { logCurrentException("writing the KDE"); ok = false; return -1; }

    try { LOD_CUTOFF = get_min_btw_modes(kdeResult->x.data(), kdeResult->y.data(), kdeResult->size, wsize); }
    catch (...)
    {
        logCurrentException("locating the minimum between LOD score modes");
        LOG.err("ERROR: Failed to find the minimum between modes in the LOD score density.");
        LOG.err("\tResults from density estimation have been written to file for inspection.");
        LOG.err("\tA cutoff can be manually specified on the command line with", ARG_LOD_CUTOFF);
        ok = false;
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
                     bool WEIGHTED, int M, double mu, int numThreads, bool PHASED, int thinStep, int LD_SUBSAMPLE,
                     const vector<ChrRole> *role)
{
    //--winsize-multi values come straight from the command line and were never
    //checked against the data; a size >= the shortest chromosome indexes the
    //window loops out of bounds (and leaves figtree with N = 0 points).
    {
        int dataBound = mapDataByChr->at(0)->nloci;
        for (unsigned int c = 1; c < mapDataByChr->size(); c++)
            if (mapDataByChr->at(c)->nloci < dataBound) dataBound = mapDataByChr->at(c)->nloci;
        for (unsigned int i = 0; i < multiWinsizes.size(); i++) {
            if (multiWinsizes[i] >= dataBound) {
                LOG.err("ERROR: --winsize-multi value", multiWinsizes[i], false);
                LOG.err(" is >= the number of loci on the shortest chromosome (", dataBound, false);
                LOG.err(").");
                throw 0;
            }
        }
    }

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
            ldDataByChr = calcLDData(hapDataByChr, freqDataByChr, mapDataByChr, genoFreqDataByChr, centro, multiWinsizes[i], MAX_GAP, PHASED, numThreads, LD_SUBSAMPLE);
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
        DoubleData *rawWinData = convertWinData2DoubleData(winDataByChr, (thinStep > 0 ? thinStep : multiWinsizes[i]), role, CHR_AUTOSOME);
        releaseWinData(winDataByChr);

        KDEResult *kdeResult = computeKDE(rawWinData->data.data(), rawWinData->size);
        releaseDoubleData(rawWinData);

        try { writeKDEResult(kdeResult, makeKDEFilename(outfile, multiWinsizes[i])); }
        catch (...) { releaseKDEResult(kdeResult); throw; }   //rethrow, not throw 0: keep the type
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
    return;
}


KDEResult *selectWinsize(vector< HapData * > *hapDataByChr,
                         vector< FreqData * > *freqDataByChr,
                         vector< MapData * > *mapDataByChr,
                         IndData *indData, centromere *centro,
                         int &winsize, int step, double error,
                         vector< GenoLikeData * > *GLDataByChr, bool USE_GL,
                         int MAX_GAP, int KDE_SUBSAMPLE, string outfile,
                         bool WEIGHTED, vector< GenoFreqData * > *genoFreqDataByChr, bool PHASED, int thinStep,
                         int MAX_WINSIZE,
                         const vector<ChrRole> *role)
{
    double AUTO_WINSIZE_THRESHOLD = AUTO_WINSIZE_THRESHOLD_G;
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

    //The search used to be unbounded: if the smoothness criterion was never
    //met, winsizeQuery grew until it exceeded the number of loci and the
    //window loops indexed out of bounds.  Bound it by --max-winsize and by the
    //shortest chromosome, and fail with something actionable.
    int dataBound = mapDataByChr->at(0)->nloci;
    for (unsigned int c = 1; c < mapDataByChr->size(); c++)
        if (mapDataByChr->at(c)->nloci < dataBound) dataBound = mapDataByChr->at(c)->nloci;

    int winsizeQuery = winsize;
    double mse;
    bool finished = false;
    while (!finished)
    {
        if (winsizeQuery > MAX_WINSIZE || winsizeQuery >= dataBound)
        {
            LOG.err("ERROR: Automatic window size search reached", winsizeQuery, false);
            LOG.err(" without meeting the smoothness threshold", AUTO_WINSIZE_THRESHOLD);
            if (winsizeQuery >= dataBound)
                LOG.err("\tThe shortest chromosome has only", dataBound, false);
            if (winsizeQuery >= dataBound) LOG.err(" loci.");
            LOG.err("\tRaise --max-winsize, set --winsize explicitly, or use --winsize-multi.");
            throw 0;
        }
        if (WEIGHTED) {
            /*
            winDataByChr = calcwLODWindows(hapDataByChrToCalc, freqDataByChr, mapDataByChr,
                                           GLDataByChrToCalc,
                                           genoFreqDataByChr,
                                           centro, winsizeQuery,
                                           error, MAX_GAP, USE_GL);
            */
            if (!LOG.isQuiet()) cerr << "Not currently supported.\n";
            throw 0;
        }
        else {
            winDataByChr = calcLODWindows(hapDataByChrToCalc, freqDataByChr, mapDataByChr,
                                          GLDataByChrToCalc,
                                          centro, winsizeQuery,
                                          error, MAX_GAP, USE_GL);
        }
        DoubleData *rawWinData = convertWinData2DoubleData(winDataByChr, (thinStep > 0 ? thinStep : winsizeQuery), role, CHR_AUTOSOME);
        releaseWinData(winDataByChr);

        KDEResult *kdeResult = computeKDE(rawWinData->data.data(), rawWinData->size);
        releaseDoubleData(rawWinData);

        mse = calculateWiggle(kdeResult);
        LOG.log("",winsizeQuery,false);
        LOG.log("\t",mse);

        if (mse <= AUTO_WINSIZE_THRESHOLD)
        {
            finished = true;
            selectedKDEResult = cloneKDEResult(kdeResult);
            winsize = winsizeQuery;
            try { writeKDEResult(selectedKDEResult, makeKDEFilename(outfile, winsize)); }
            catch (...) { throw; }   //rethrow, not throw 0: keep the type
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
                                 bool WEIGHTED, vector< GenoFreqData * > *genoFreqDataByChr, bool PHASED, int thinStep,
                                 const vector<ChrRole> *role)
{
    //--winsize-multi values are taken verbatim from the command line and were
    //never checked against the data; a size >= the shortest chromosome indexes
    //the window loops out of bounds.
    {
        int dataBound = mapDataByChr->at(0)->nloci;
        for (unsigned int c = 1; c < mapDataByChr->size(); c++)
            if (mapDataByChr->at(c)->nloci < dataBound) dataBound = mapDataByChr->at(c)->nloci;
        for (unsigned int i = 0; i < multiWinsizes->size(); i++) {
            if (multiWinsizes->at(i) >= dataBound) {
                LOG.err("ERROR: --winsize-multi value", multiWinsizes->at(i), false);
                LOG.err(" is >= the number of loci on the shortest chromosome (", dataBound, false);
                LOG.err(").");
                throw 0;
            }
        }
    }
    double AUTO_WINSIZE_THRESHOLD = AUTO_WINSIZE_THRESHOLD_G;
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
            if (!LOG.isQuiet()) cerr << "Not currently supported.\n";
            throw 0;
        }
        else {
            winDataByChr = calcLODWindows(hapDataByChrToCalc, freqDataByChr, mapDataByChr,
                                          GLDataByChrToCalc,
                                          centro, multiWinsizes->at(i),
                                          error, MAX_GAP, USE_GL);
        }
        DoubleData *rawWinData = convertWinData2DoubleData(winDataByChr, (thinStep > 0 ? thinStep : multiWinsizes->at(i)), role, CHR_AUTOSOME);
        releaseWinData(winDataByChr);

        KDEResult *kdeResult = computeKDE(rawWinData->data.data(), rawWinData->size);
        releaseDoubleData(rawWinData);

        mse = calculateWiggle(kdeResult);
        LOG.log("",multiWinsizes->at(i),false);
        LOG.log("\t",mse);

        if (mse <= AUTO_WINSIZE_THRESHOLD || i == multiWinsizes->size() - 1)
        {
            selectedKDEResult = cloneKDEResult(kdeResult);
            winsize = multiWinsizes->at(i);
            try { writeKDEResult(selectedKDEResult, makeKDEFilename(outfile, winsize)); }
            catch (...) { throw; }   //rethrow, not throw 0: keep the type
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

    int ngaussians = NCLUST;
    size_t maxIter = size_t(GMM_MAX_ITER);
    double tolerance = GMM_TOL;
    //Vectors, not new[]: gmm.estimate() below can throw (a collapsed component
    //trips garlicLogChecked), and the raw version leaked all four on that path.
    vector<double> W(ngaussians);
    vector<double> Mu(ngaussians);
    vector<double> Sigma(ngaussians);
    vector<size_t> sortIndex(ngaussians);

    //calculate mean and var for the population size distribution to use for initial guess
    double var = garlicVariance(rohLength->length.data(), rohLength->size);
    double mu = garlicMean(rohLength->length.data(), rohLength->size);
    for (int n = 0; n < ngaussians; n++)
    {
        W[n] = 1.0 / double(ngaussians);
        Mu[n] = mu * double(n + 1) / double(ngaussians + 1);
        Sigma[n] = var * (n + 1) / double(ngaussians);
    }

    //verbose = !isQuiet: the EM iteration counter is progress output, and it
    //went to stderr even under --quiet.  GMM already gates it on its own
    //`verbose` member; the caller was passing an unconditional true.
    GMM gmm(ngaussians, W.data(), Mu.data(), Sigma.data(), maxIter, tolerance, !LOG.isQuiet(), true);

    gmm.estimate(rohLength->length.data(), rohLength->size);

    for (int n = 0; n < ngaussians; n++)
    {
        W[n] = gmm.getMixCoefficient(n);
        Mu[n] = gmm.getMean(n);
        Sigma[n] = gmm.getVar(n);
        sortIndex[n] = n;
    }

    garlicSortIndex(sortIndex.data(), Mu.data(), ngaussians);
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

    return bounds;
}




