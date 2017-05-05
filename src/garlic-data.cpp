#include "garlic-data.h"

void loadTPEDData(string tpedfile, int &numLoci, int &numInd,
                                   vector< HapData * > **hapDataByChr,
                                   vector< MapData * > **mapDataByChr,
                                   vector< FreqData * > **freqDataByChr,
                                   char TPED_MISSING, int nresample, bool PHASED, bool AUTO_FREQ){
    const gsl_rng_type *T;
    gsl_rng *r;
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r, time(NULL));

    igzstream fin;
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        LOG.err("ERROR: Failed to open", filename);
        throw 0;
    }

    string line;
    stringstream ss;
    char oneAllele = TPED_MISSING;
    int currChrLoci = 0;
    int ncols = 0;
    int nalleles = 0;
    int total = 0;
    string chr, locusName, prevChr;
    string emptyChr = "_nochr";
    double gpos, ppos;

    vector<double> geneticPos, physicalPos;
    vector<string> locusNames;
    vector<char> allele;
    
    char alleleStr1, alleleStr2;

    vector< double > freq;
    vector< short * > hap;
    vector< bool * > fc;
    short *data;
    bool *firstCopy;

    while(getline(fin,line)){
        numLoci++;
        ncols = countFields(line) - 4;
        numInd = ncols/2;

        ss.str(line);

        ss >> chr;

        if (prevChr.compare(emptyChr) == 0 && numLoci == 1) prevChr = chr;

        if (chr.compare(prevChr) != 0){
            
            LOG.logn("Chromosome",checkChrName(chr));
            LOG.logn(":",currChrLoci);
            LOG.log(" sites.");

            (*mapDataByChr)->push_back(initMapData(geneticPos,physicalPos,locusNames, allele, currChrLoci, checkChrName(chr)));
            geneticPos.clear();
            physicalPos.clear();
            allele.clear();
            locusNames.clear();

            (*hapDataByChr)->push_back(initHapData(hap, fc, currChrLoci, numInd, PHASED));
            hap.clear();
            fc.clear();

            (*freqDataByChr)->push_back(initFreqData(freq, nloci));
            freq.clear();

            prevChr = chr;
            currChrLoci = 0;
        }

        currChrLoci++;
        
        ss >> locusName;
        locusNames.push_back(locusName);
        ss >> gpos;
        geneticPos.push_back(gpos);
        ss >> ppos;
        physicalPos.push_back(ppos);

        //alleles remain in ss
        nalleles = 0;
        total = 0;
        data = new short[numInd];
        if(PHASED) firstCopy = new bool[numInd];

        for(int i = 0; i < numInd; i++){
            data[i] = 0;
            ss >> alleleStr1 >> alleleStr2;
            if(oneAllele == TPED_MISSING && alleleStr1 != TPED_MISSING) oneAllele = alleleStr1;
            if(oneAllele == TPED_MISSING && alleleStr2 != TPED_MISSING) oneAllele = alleleStr2;
            if (alleleStr1 == TPED_MISSING) data[i] += -9;
            else if (alleleStr1 == oneAllele){
                data[i] += 1;
                nalleles++;
                total++;
            }
            else total++;
            if (alleleStr2 == TPED_MISSING) data[i] += -9;
            else if (alleleStr2 == oneAllele){
                data[i] += 1;
                nalleles++;
                total++;
            }
            else total++;
            if (data[i] < 0) data[i] = -9;
            if(PHASED) firstCopy[i] = (alleleStr1 == oneAllele);
        }

        allele.push_back(oneAllele);
        hap.push_back(data);
        data = NULL;
        if(PHASED){
            fc.push_back(firstCopy);
            firstCopy = NULL;
        }

        if(AUTO_FREQ){
            double freqtmp = double(nalleles)/double(total);
            if (nresample > 0){
                int count = 0;
                for (int i = 0; i < nresample; i++){
                    if (gsl_rng_uniform(r) <= freqtmp) count++;
                }
                freqtmp = count / nresample;
            }
            freq.push_back(freqtmp);
        }

        ss.clear();
    }

    LOG.logn("Chromosome",checkChrName(chr));
    LOG.logn(":",currChrLoci);
    LOG.log(" sites.");

    (*mapDataByChr)->push_back(initMapData(geneticPos,physicalPos,locusNames, allele, currChrLoci, checkChrName(chr)));
    geneticPos.clear();
    physicalPos.clear();
    allele.clear();
    locusNames.clear();

    (*hapDataByChr)->push_back(initHapData(hap, fc, currChrLoci, numInd, PHASED));
    hap.clear();
    fc.clear();

    (*freqDataByChr)->push_back(initFreqData(freq, nloci));
    freq.clear();

    return;
}

FreqData *initFreqData(const vector<double> &freq, int nloci){
    FreqData *freqData = new FreqData;
    freqData->freq = new double[nloci];

    for(int i = 0; i < nloci; i++){
        freqData->freq[i] = freq[i];
    }

    return freqData;
}

HapData *initHapData(const vector< short * > &hap,
                     const vector< bool * > &fc,
                     int nloci, int nind, bool PHASED){
    HapData *hapData = new HapData;
    hapData->nind = nind;
    hapData->nloci = nloci;
    hapData->data = new short*[nloci];
    if(PHASED) hapData->firstCopy = new bool*[nloci];

    for(int i = 0; i < nloci; i++){
        hapData->data[i] = hap[i];
        if(PHASED) hapData->firstCopy[i] = fc[i];
    }

    return hapData;
}

MapData *initMapData(const vector<double> &geneticPos,
                     const vector<double> &physicalPos,
                     const vector<string> &locusNames,
                     const vector<char> &allele,
                     int nloci, string chr){

    MapData *mapData = new MapData;
    mapData->physicalPos = new double[nloci];
    mapData->geneticPos = new double[nloci];
    mapData->locusName = new string[nloci];
    mapData->allele = new char[nloci];
    mapData->nloci = nloci;
    mapData->chr = chr;

    for(int i = 0; i < nloci; i++){
        mapData->physicalPos[i] = physicalPos[i];
        mapData->geneticPos[i] = geneticPos[i];
        mapData->locusName[i] = locusNames[i];
        mapData->allele[i] = allele[i];
    }

    return mapData;
}

void freqOnly(string filename, string outfile, int nresample, char TPED_MISSING){
    
    const gsl_rng_type *T;
    gsl_rng *r;
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r, time(NULL));

    string freqoutfile = outfile + ".freq.gz";

    ogzstream fout;
    fout.open(freqoutfile.c_str());
    if (fout.fail())
    {
        cerr << "ERROR: Failed to open " << freqoutfile << " for writing.\n";
        LOG.err("ERROR: Failed to open", freqoutfile);
        throw 0;
    }

    fout << "CHR\tSNP\tPOS\tALLELE\tFREQ\n";

    igzstream fin;
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        LOG.err("ERROR: Failed to open", filename);
        throw 0;
    }

    string junk;
    char oneAllele;
    int count;
    string line;
    int nloci = 0;
    int ncols;
    string chr, locusName;
    double gpos, ppos;
    double nalleles = 0;
    double total = 0;
    stringstream ss;
    while(getline(fin,line)){
        nloci++;
        ncols = countFields(line);

        ss.str(line);
        ss >> chr;
        ss >> locusName;
        ss >> gpos;
        ss >> ppos;

        oneAllele = TPED_MISSING;
        total = 0;
        nalleles = 0;
        for(count = 0; count < ncols-4; count++){
            ss >> junk;
            if(junk[0] != TPED_MISSING){
                total++;
                if(oneAllele == TPED_MISSING) oneAllele = junk.c_str()[0]; 
                if(junk[0] == oneAllele) nalleles++;
            }
        }

        double freq = nalleles/total;
        if (nresample > 0){
            count = 0;
            for (int i = 0; i < nresample; i++){
                if (gsl_rng_uniform(r) <= freq) count++;
            }
            freq = count / nresample;
        }
        ss.clear();
        
        fout << checkChrName(chr) << "\t" << locusName << "\t" << int(ppos) << "\t" << oneAllele << "\t" << freq << endl;
    }

    gsl_rng_free(r);

    fin.close();
    fout.close();
}


double calcDensity(int numLoci, vector< MapData * > *mapDataByChr, centromere *centro){
    double density = numLoci;
    double length = 0;
    for(unsigned int i = 0; i < mapDataByChr->size(); i++){
        char str[3];
        sprintf(str,"%d",i+1);
        string chrstr = checkChrName(string(str));
        int nloci = mapDataByChr->at(i)->nloci;
        length += mapDataByChr->at(i)->physicalPos[nloci-1] - mapDataByChr->at(i)->physicalPos[0] + 1 - (centro->centromereEnd(chrstr) - centro->centromereStart(chrstr));
    }
    return density/length;
}

vector< LDData * > *calcLDData(vector< HapData * > *hapDataByChr, 
                               vector< FreqData * > *freqDataByChr,
                               vector< MapData * > *mapDataByChr,
                               vector< GenoFreqData * > *genoFreqDataByChr,
                               centromere *centro,
                               int winsize,
                               int MAX_GAP,
                               bool PHASED,
                               int numThreads){
    vector< LDData * > *ldDataByChr = new vector< LDData * >;
    for(unsigned int chr = 0; chr < hapDataByChr->size(); chr++){

        if(!PHASED) ldDataByChr->push_back(calcHR2LD(hapDataByChr->at(chr), genoFreqDataByChr->at(chr), winsize, numThreads));
        else ldDataByChr->push_back(calcR2LD(hapDataByChr->at(chr), freqDataByChr->at(chr), winsize, numThreads));
    }
    return ldDataByChr;
}

LDData *calcHR2LD(HapData *hapData, GenoFreqData *genoFreqData, int winsize, int numThreads){
    
    LDData *LD = initLDData(hapData->nloci, winsize);

    unsigned int *NUM_PER_THREAD = make_thread_partition(numThreads, hapData->nloci);

    HR2_work_order_t *order;
    pthread_t *peer = new pthread_t[numThreads];
    vector< HR2_work_order_t * > orders;
    unsigned int previous = 0;
    for (int i = 0; i < numThreads; i++)
    {
        order = new HR2_work_order_t;
        order->hapData = hapData;
        order->genoFreqData = genoFreqData;
        order->LD = LD;
        order->winsize = winsize;
        //order->id = i;
        order->start = previous;
        previous += NUM_PER_THREAD[i];
        order->stop = previous;

        pthread_create(&(peer[i]),
                       NULL,
                       (void *(*)(void *))parallelHR2,
                       (void *)order);
        orders.push_back(order);
    }

    for (int i = 0; i < numThreads; i++){
        pthread_join(peer[i], NULL);
        delete orders[i];
    }

    orders.clear();

    delete [] peer;

    return LD;
}

LDData *calcR2LD(HapData *hapData, FreqData *freqData, int winsize, int numThreads){
    
    LDData *LD = initLDData(hapData->nloci, winsize);

    unsigned int *NUM_PER_THREAD = make_thread_partition(numThreads, hapData->nloci);

    R2_work_order_t *order;
    pthread_t *peer = new pthread_t[numThreads];
    vector< R2_work_order_t * > orders;
    unsigned int previous = 0;
    for (int i = 0; i < numThreads; i++)
    {
        order = new R2_work_order_t;
        order->hapData = hapData;
        order->freqData = freqData;
        order->LD = LD;
        order->winsize = winsize;
        //order->id = i;
        order->start = previous;
        previous += NUM_PER_THREAD[i];
        order->stop = previous;

        pthread_create(&(peer[i]),
                       NULL,
                       (void *(*)(void *))parallelR2,
                       (void *)order);
        orders.push_back(order);
    }

    for (int i = 0; i < numThreads; i++){
        pthread_join(peer[i], NULL);
        delete orders[i];
    }

    orders.clear();

    delete [] peer;

    return LD;
}


void parallelHR2(void *order){
    HR2_work_order_t *p = (HR2_work_order_t *)order;
    LDData *LD = p->LD;
    HapData *hapData = p->hapData;
    GenoFreqData *genoFreqData = p->genoFreqData;
    int winsize = p->winsize;
    int start = p->start;
    int stop = p->stop;

    if (hapData->nloci - stop < winsize) stop = hapData->nloci - winsize + 1;

    for (int locus = start; locus < stop; locus++){
        for (int i = locus; i < locus + winsize; i++){
            ldHR2(LD, hapData, genoFreqData, i, locus, locus + winsize - 1);
        }
    }
}

void parallelR2(void *order){
    R2_work_order_t *p = (R2_work_order_t *)order;
    LDData *LD = p->LD;
    HapData *hapData = p->hapData;
    FreqData *freqData = p->freqData;
    int winsize = p->winsize;
    int start = p->start;
    int stop = p->stop;

    if (hapData->nloci - stop < winsize) stop = hapData->nloci - winsize + 1;

    for (int locus = start; locus < stop; locus++){
        for (int i = locus; i < locus + winsize; i++){
            ldR2(LD, hapData, freqData, i, locus, locus + winsize - 1);
        }
    }
}


void ldHR2(LDData *LD, HapData *hapData, GenoFreqData *genoFreqData, int site, int start, int end) {
    for (int i = start; i <= end; i++) {
        if (i != site) LD->LD[start][site-start] += hr2(hapData, genoFreqData, i, site);
        else LD->LD[start][site-start] += 1;
    }
    return;
}

void ldR2(LDData *LD, HapData *hapData, FreqData *freqData, int site, int start, int end) {
    for (int i = start; i <= end; i++) {
        if (i != site) LD->LD[start][site-start] += r2(hapData, freqData, i, site);
        else LD->LD[start][site-start] += 1;
    }
    return;
}


unsigned int *make_thread_partition(int &numThreads, int ncols) {
    if (numThreads > ncols) numThreads = ncols;
    unsigned int *NUM_PER_THREAD = new unsigned int[numThreads];
    unsigned int div = ncols / numThreads;

    for (int i = 0; i < numThreads; i++)
    {
        NUM_PER_THREAD[i] = 0;
        NUM_PER_THREAD[i] += div;
    }

    for (int i = 0; i < ncols % numThreads; i++)
    {
        NUM_PER_THREAD[i]++;
    }

    return NUM_PER_THREAD;
}


double hr2(HapData *hapData, GenoFreqData *genoFreqData, int i, int j) {
    double HA = genoFreqData->homFreq[i];
    double HB = genoFreqData->homFreq[j];
    if(HA > 0 && HA < 1 && HB > 0 && HB < 1){
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
        double HR2 = H * H / (HA * (1 - HA) * HB * (1 - HB));
        if(HR2 > 1) return 1;
        else return HR2;
    }
    else{
        return 0;
    }
}

double r2(HapData *hapData, FreqData *freqData, int i, int j) {
    
    double pi = freqData->freq[i];
    double pj = freqData->freq[j];

    if(pi > 0 && pi < 1 && pj > 0 && pj < 1){
        double x11 = 0;
        double total = 0;
        for (int ind = 0; ind < hapData->nind; ind++) {
            if (hapData->data[i][ind] != -9 && hapData->data[j][ind] != -9) {
                total+=2;
                if (hapData->data[i][ind] == 2 && hapData->data[j][ind] == 2) x11+=2;
                else if (hapData->data[i][ind] == 1 && hapData->data[j][ind] == 2) x11++;
                else if (hapData->data[i][ind] == 2 && hapData->data[j][ind] == 1) x11++;
                else if (hapData->data[i][ind] == 1 && 
                         hapData->data[j][ind] == 1 && 
                         hapData->firstCopy[j][ind] == hapData->firstCopy[i][ind]){
                    x11++;
                }
            }
        }
        x11 /= total;
        double D = x11 - pi * pj;
        double R2 = D * D / (pi * (1 - pi) * pj * (1 - pj));
        if(R2 > 1) return 1;
        else return R2;
    }
    else{
        return 0;
    }
}

LDData *initLDData(int nloci, int winsize){
    LDData *data = new LDData;
    data->nloci = nloci;
    data->winsize = winsize;
    data->LD = new double*[nloci];
    for(int i = 0; i < nloci; i++){
        data->LD[i] = new double[winsize];
        for(int j = 0; j < winsize; j++){
            data->LD[i][j] = 0;
        }
    }
    return data;
}
void releaseLDData(LDData *data){
    for(int i = 0; i < data->nloci; i++){
        delete [] data->LD[i];
    }
    delete [] data->LD;
    delete data;
}
void releaseLDData(vector< LDData * > *ldDataByChr){
    for (unsigned int chr = 0; chr < ldDataByChr->size(); chr++){
        releaseLDData(ldDataByChr->at(chr));
    }
    ldDataByChr->clear();
    delete ldDataByChr;
    return;
}

vector< GenoFreqData * > *calculateGenoFreq(vector <HapData *> *hapDataByChr){
    vector< GenoFreqData * > *genoFreqDataByChr = new vector< GenoFreqData * >;
    for(unsigned int chr = 0; chr < hapDataByChr->size(); chr++){
        genoFreqDataByChr->push_back(calculateGenoFreq(hapDataByChr->at(chr)));
    }
    return genoFreqDataByChr;
}

GenoFreqData *calculateGenoFreq(HapData *hapData){
    GenoFreqData *genoFreqData = initGenoFreq(hapData->nloci);
    double total, freqHom;

    for (int locus = 0; locus < hapData->nloci; locus++)
    {
        total = 0;
        freqHom = 0;
        for (int ind = 0; ind < hapData->nind; ind++)
        {
            if (hapData->data[locus][ind] != -9)
            {
                if(hapData->data[locus][ind] == 2 || hapData->data[locus][ind] == 0) freqHom++;
                total++;
            }
        }
        freqHom /= total;
        genoFreqData->homFreq[locus] = freqHom;
    }
    return genoFreqData;
}

GenoFreqData *initGenoFreq(int nloci){
    GenoFreqData *data = new GenoFreqData;
    data->homFreq = new double[nloci];
    data->nloci = nloci;
    return data;
}

void releaseGenoFreq(GenoFreqData *genoFreqData){
    delete [] genoFreqData->homFreq;
    delete genoFreqData;
    return;
}

void releaseGenoFreq(vector< GenoFreqData * > *genoFreqDataByChr)
{
    for (unsigned int chr = 0; chr < genoFreqDataByChr->size(); chr++)
    {
        releaseGenoFreq(genoFreqDataByChr->at(chr));
    }
    genoFreqDataByChr->clear();
    delete genoFreqDataByChr;
    return;
}

int interpolateGeneticmap(vector< MapData * > **mapDataByChr, vector< GenMapScaffold * > *scaffoldMapByChr) {
    int numInterpolated = 0;
    for (unsigned int i = 0; i < (*mapDataByChr)->size(); i++) {
        numInterpolated += interpolateGeneticmap((*mapDataByChr)->at(i), scaffoldMapByChr->at(i));
    }
    return numInterpolated;
}

int interpolateGeneticmap(MapData *mapData, GenMapScaffold *scaffold) {
    int numInterpolated = 0;
    for (int i = 0; i < mapData->nloci; i++) {
        mapData->geneticPos[i] = getMapInfo(mapData->physicalPos[i], scaffold, numInterpolated);
    }
    return numInterpolated;
}

double getMapInfo(int queryPos, GenMapScaffold *scaffold, int &count) {
    if (queryPos < scaffold->physicalPos[0])
    {
        LOG.err("ERROR: Sites outside of map scaffold should have been filtered out by GARLIC.");
        throw 0;
    }
    else if (queryPos > scaffold->physicalPos[scaffold->nloci - 1])
    {
        LOG.err("ERROR: Sites outside of map scaffold should have been filtered out by GARLIC.");
        throw 0;
    }
    else if (scaffold->ppos2index.count(queryPos) > 0)
    {
        int index = scaffold->ppos2index[queryPos];
        return scaffold->geneticPos[index];
    }
    else
    {
        int startIndex;
        int endIndex;
        for (/*scaffold->currentIndex*/; scaffold->currentIndex < scaffold->nloci - 1; scaffold->currentIndex++)
        {
            if (queryPos > scaffold->physicalPos[scaffold->currentIndex] && queryPos < scaffold->physicalPos[scaffold->currentIndex + 1])
            {
                startIndex = scaffold->currentIndex;
                endIndex = scaffold->currentIndex + 1;
                break;
            }
        }
        count++;
        return interpolate(scaffold->physicalPos[startIndex], scaffold->geneticPos[startIndex],
                           scaffold->physicalPos[endIndex], scaffold->geneticPos[endIndex],
                           queryPos);
    }
}

double interpolate(double x0, double y0, double x1, double y1, double query)
{
    return ( ( (y1 - y0) / (x1 - x0) ) * query + ( y0 - ((y1 - y0) / (x1 - x0)) * x0 ) );
}


vector< GenMapScaffold *> *loadMapScaffold(string mapfile, centromere *centro) {
    igzstream fin;
    cerr << "Opening " << mapfile << "...\n";
    fin.open(mapfile.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << mapfile << " for reading.\n";
        throw 0;
    }

    string line;
    vector<int> nlociByChr;
    int nloci = 0;
    int currentLoci = 0;
    int num_cols = 4;
    int current_cols = 0;
    stringstream ss;
    string currentChr;
    string chrstr;
    while (getline(fin, line))
    {
        nloci++;
        current_cols = countFields(line);
        if (current_cols != num_cols)
        {
            cerr << "ERROR: line " << nloci << " of " << mapfile << " has " << current_cols
                 << ", but expected " << num_cols << ".\n";
            throw 0;
        }

        ss.str(line);
        ss >> chrstr;

        if (nloci == 1) {
            currentChr = chrstr;
        }

        if (currentChr.compare(chrstr) == 0) {
            currentLoci++;
        }
        else {
            nlociByChr.push_back(currentLoci);
            currentLoci = 1;
            currentChr = chrstr;
        }
        ss.clear();
    }

    nlociByChr.push_back(currentLoci);

    fin.close();
    fin.clear();
    fin.open(mapfile.c_str());

    cerr << "Loading genetic map scaffold for " << nloci << " loci.\n";

    string locusName;
    string chrname;
    vector< GenMapScaffold * > *scaffoldMapByChr = new vector< GenMapScaffold * >;

    for (unsigned int chr = 0; chr < nlociByChr.size(); chr++) {
        GenMapScaffold *scaffold = initGenMapScaffold(nlociByChr.at(chr));
        for (int locus = 0; locus < nlociByChr.at(chr); locus++)
        {
            getline(fin,line);
            ss.str(line);
            //fin >> scaffold->chr;
            ss >> chrname;
            scaffold->chr = checkChrName(chrname);
            ss >> locusName;
            ss >> scaffold->geneticPos[locus];
            ss >> scaffold->physicalPos[locus];
            scaffold->ppos2index[scaffold->physicalPos[locus]] = locus;
            ss.clear();
        }
        scaffold->centroStart = centro->centromereStart(scaffold->chr);
        scaffold->centroEnd = centro->centromereEnd(scaffold->chr);
        scaffoldMapByChr->push_back(scaffold);
    }

    fin.close();

    return scaffoldMapByChr;
}

GenMapScaffold *initGenMapScaffold(int nloci) {
    GenMapScaffold *scaffoldMap = new GenMapScaffold;
    scaffoldMap->physicalPos = new int[nloci];
    scaffoldMap->geneticPos = new double[nloci];
    scaffoldMap->nloci = nloci;
    scaffoldMap->currentIndex = 0;
    return scaffoldMap;
}

void releaseGenMapScaffold(GenMapScaffold *scaffoldMap) {
    delete [] scaffoldMap->physicalPos;
    delete [] scaffoldMap->geneticPos;
    scaffoldMap->ppos2index.clear();
    delete scaffoldMap;
    return;
}

void releaseGenMapScaffold(vector< GenMapScaffold * > *scaffoldMapByChr) {
    for (unsigned int i = 0; i < scaffoldMapByChr->size(); i++) {
        releaseGenMapScaffold(scaffoldMapByChr->at(i));
    }
    delete scaffoldMapByChr;
    return;
}

int filterMonomorphicSites(vector< MapData * > **mapDataByChr,
                           vector< HapData * > **hapDataByChr,
                           vector< FreqData * > **freqDataByChr,
                           vector< GenoLikeData * > **GLDataByChr,
                           bool USE_GL, bool PHASED)
{
    vector< MapData * > *mapDataByChr2 = new vector< MapData * >;
    vector< HapData * > *hapDataByChr2 = new vector< HapData * >;
    vector< FreqData * > *freqDataByChr2 = new vector< FreqData * >;
    vector< GenoLikeData * > *GLDataByChr2;

    if(USE_GL){
        GLDataByChr2 = new vector< GenoLikeData * >;
    }

    int numLoci = 0;
    for (unsigned int i = 0; i < (*mapDataByChr)->size(); i++) {
        int newLoci = 0;
        MapData *mapData2 = filterMonomorphicSites((*mapDataByChr)->at(i), (*freqDataByChr)->at(i), newLoci);
        HapData *hapData2 = filterMonomorphicSites((*hapDataByChr)->at(i), (*freqDataByChr)->at(i), newLoci, PHASED);
        GenoLikeData *GLData2;
        if(USE_GL){
            GLData2 = filterMonomorphicSites((*GLDataByChr)->at(i), (*freqDataByChr)->at(i), newLoci);
        }
        FreqData *freqData2 = filterMonomorphicSites((*freqDataByChr)->at(i), newLoci);
        mapDataByChr2->push_back(mapData2);
        hapDataByChr2->push_back(hapData2);
        if(USE_GL) GLDataByChr2->push_back(GLData2);
        freqDataByChr2->push_back(freqData2);
        numLoci += newLoci;
    }

    releaseMapData(*mapDataByChr);
    releaseHapData(*hapDataByChr);
    releaseFreqData(*freqDataByChr);
    if(USE_GL) releaseGLData(*GLDataByChr);

    *mapDataByChr = mapDataByChr2;
    *hapDataByChr = hapDataByChr2;
    if(USE_GL) *GLDataByChr = GLDataByChr2;
    *freqDataByChr = freqDataByChr2;

    return numLoci;
}

int filterMonomorphicAndOOBSites(vector< MapData * > **mapDataByChr,
                                 vector< HapData * > **hapDataByChr,
                                 vector< FreqData * > **freqDataByChr,
                                 vector< GenoLikeData * > **GLDataByChr,
                                 vector< GenMapScaffold * > *scaffoldMapByChr,
                                 bool USE_GL, bool PHASED)
{
    vector< MapData * > *mapDataByChr2 = new vector< MapData * >;
    vector< HapData * > *hapDataByChr2 = new vector< HapData * >;
    vector< FreqData * > *freqDataByChr2 = new vector< FreqData * >;
    vector< GenoLikeData * > *GLDataByChr2;

    if(USE_GL){
        GLDataByChr2 = new vector< GenoLikeData * >;
    }

    int numLoci = 0;
    for (unsigned int i = 0; i < (*mapDataByChr)->size(); i++) {
        int newLoci = 0;
        MapData *mapData2 = filterMonomorphicAndOOBSites((*mapDataByChr)->at(i), (*freqDataByChr)->at(i), scaffoldMapByChr->at(i), newLoci);
        HapData *hapData2 = filterMonomorphicAndOOBSites((*hapDataByChr)->at(i), (*mapDataByChr)->at(i), (*freqDataByChr)->at(i), scaffoldMapByChr->at(i), newLoci, PHASED);
        GenoLikeData *GLData2; 
        if(USE_GL){
            GLData2 = filterMonomorphicAndOOBSites((*GLDataByChr)->at(i), (*mapDataByChr)->at(i), (*freqDataByChr)->at(i), scaffoldMapByChr->at(i), newLoci);
        }
        FreqData *freqData2 = filterMonomorphicAndOOBSites((*freqDataByChr)->at(i), (*mapDataByChr)->at(i), scaffoldMapByChr->at(i), newLoci);
        mapDataByChr2->push_back(mapData2);
        hapDataByChr2->push_back(hapData2);
        if(USE_GL) GLDataByChr2->push_back(GLData2);
        freqDataByChr2->push_back(freqData2);
        numLoci += newLoci;
    }

    releaseMapData(*mapDataByChr);
    releaseHapData(*hapDataByChr);
    releaseFreqData(*freqDataByChr);
    if(USE_GL) releaseGLData(*GLDataByChr);

    *mapDataByChr = mapDataByChr2;
    *hapDataByChr = hapDataByChr2;
    if(USE_GL) *GLDataByChr = GLDataByChr2;
    *freqDataByChr = freqDataByChr2;

    return numLoci;

}

MapData *filterMonomorphicSites(MapData *mapData, FreqData *freqData, int &newLoci)
{
    if (newLoci <= 0) {
        newLoci = 0;
        for (int i = 0; i < freqData->nloci; i++)
            if (freqData->freq[i] > 0 && freqData->freq[i] < 1)
                newLoci++;
    }

    MapData *mapData2 = initMapData(newLoci);
    mapData2->chr = mapData->chr;
    int index = 0;
    for (int i = 0; i < freqData->nloci; i++)
    {
        if (freqData->freq[i] > 0 && freqData->freq[i] < 1)
        {
            mapData2->physicalPos[index] = mapData->physicalPos[i];
            mapData2->geneticPos[index] = mapData->geneticPos[i];
            mapData2->locusName[index] = mapData->physicalPos[i];
            mapData2->allele[index] = mapData->allele[i];
            index++;
        }
    }

    return mapData2;
}

HapData *filterMonomorphicSites(HapData *hapData, FreqData *freqData, int &newLoci, bool PHASED)
{
    if (newLoci <= 0) {
        newLoci = 0;
        for (int i = 0; i < freqData->nloci; i++)
            if (freqData->freq[i] > 0 && freqData->freq[i] < 1)
                newLoci++;
    }

    HapData *hapData2 = initHapData(hapData->nind, newLoci, PHASED);
    int index = 0;
    for (int i = 0; i < freqData->nloci; i++)
    {
        if (freqData->freq[i] > 0 && freqData->freq[i] < 1)
        {
            for (int j = 0; j < hapData->nind; j++)
            {
                hapData2->data[index][j] = hapData->data[i][j];
                if(PHASED) hapData2->firstCopy[index][j] = hapData->firstCopy[i][j];
            }
            index++;
        }
    }

    return hapData2;
}

GenoLikeData *filterMonomorphicSites(GenoLikeData *GLData, FreqData *freqData, int &newLoci)
{
    if (newLoci <= 0) {
        newLoci = 0;
        for (int i = 0; i < freqData->nloci; i++)
            if (freqData->freq[i] > 0 && freqData->freq[i] < 1)
                newLoci++;
    }

    GenoLikeData *GLData2 = initGLData(GLData->nind, newLoci);
    int index = 0;
    for (int i = 0; i < freqData->nloci; i++)
    {
        if (freqData->freq[i] > 0 && freqData->freq[i] < 1)
        {
            for (int j = 0; j < GLData->nind; j++)
            {
                GLData2->data[index][j] = GLData->data[i][j];
            }
            index++;
        }
    }

    return GLData2;
}

FreqData *filterMonomorphicSites(FreqData *freqData, int &newLoci)
{
    if (newLoci <= 0) {
        newLoci = 0;
        for (int i = 0; i < freqData->nloci; i++)
            if (freqData->freq[i] > 0 && freqData->freq[i] < 1)
                newLoci++;
    }

    FreqData *freqData2 = initFreqData(newLoci);
    int index = 0;
    for (int i = 0; i < freqData->nloci; i++)
    {
        if (freqData->freq[i] > 0 && freqData->freq[i] < 1)
        {
            freqData2->freq[index] = freqData->freq[i];
            index++;
        }
    }

    return freqData2;
}

MapData *filterMonomorphicAndOOBSites(MapData *mapData, FreqData *freqData, GenMapScaffold *scaffold, int &newLoci) {
    if (newLoci <= 0) {
        newLoci = 0;
        for (int i = 0; i < freqData->nloci; i++) {
            if ((freqData->freq[i] > 0 && freqData->freq[i] < 1) &&
                    !(mapData->physicalPos[i] < scaffold->physicalPos[0]) &&
                    !(mapData->physicalPos[i] > scaffold->physicalPos[scaffold->nloci - 1]) &&
                    !(mapData->physicalPos[i] > scaffold->centroStart && mapData->physicalPos[i] < scaffold->centroEnd)) {
                newLoci++;
            }
        }
    }

    MapData *mapData2 = initMapData(newLoci);
    mapData2->chr = mapData->chr;
    int index = 0;
    for (int i = 0; i < freqData->nloci; i++)
    {
        if ((freqData->freq[i] > 0 && freqData->freq[i] < 1) &&
                !(mapData->physicalPos[i] < scaffold->physicalPos[0]) &&
                !(mapData->physicalPos[i] > scaffold->physicalPos[scaffold->nloci - 1]) &&
                !(mapData->physicalPos[i] > scaffold->centroStart && mapData->physicalPos[i] < scaffold->centroEnd))
        {
            mapData2->physicalPos[index] = mapData->physicalPos[i];
            mapData2->geneticPos[index] = mapData->geneticPos[i];
            mapData2->locusName[index] = mapData->physicalPos[i];
            mapData2->allele[index] = mapData->allele[i];
            index++;
        }
    }

    return mapData2;
}

HapData *filterMonomorphicAndOOBSites(HapData *hapData, MapData *mapData, FreqData *freqData, GenMapScaffold *scaffold, int &newLoci, bool PHASED) {
    if (newLoci <= 0) {
        newLoci = 0;
        for (int i = 0; i < freqData->nloci; i++) {
            if ((freqData->freq[i] > 0 && freqData->freq[i] < 1) &&
                    !(mapData->physicalPos[i] < scaffold->physicalPos[0]) &&
                    !(mapData->physicalPos[i] > scaffold->physicalPos[scaffold->nloci - 1]) &&
                    !(mapData->physicalPos[i] > scaffold->centroStart && mapData->physicalPos[i] < scaffold->centroEnd)) {
                newLoci++;
            }
        }
    }

    HapData *hapData2 = initHapData(hapData->nind, newLoci, PHASED);
    int index = 0;
    for (int i = 0; i < freqData->nloci; i++)
    {
        if ((freqData->freq[i] > 0 && freqData->freq[i] < 1) &&
                !(mapData->physicalPos[i] < scaffold->physicalPos[0]) &&
                !(mapData->physicalPos[i] > scaffold->physicalPos[scaffold->nloci - 1]) &&
                !(mapData->physicalPos[i] > scaffold->centroStart && mapData->physicalPos[i] < scaffold->centroEnd))
        {
            for (int j = 0; j < hapData->nind; j++)
            {
                hapData2->data[index][j] = hapData->data[i][j];
                if(PHASED) hapData2->firstCopy[index][j] = hapData->firstCopy[i][j];
            }
            index++;
        }
    }

    return hapData2;
}

GenoLikeData *filterMonomorphicAndOOBSites(GenoLikeData *GLData, MapData *mapData, FreqData *freqData, GenMapScaffold *scaffold, int &newLoci) {
    if (newLoci <= 0) {
        newLoci = 0;
        for (int i = 0; i < freqData->nloci; i++) {
            if ((freqData->freq[i] > 0 && freqData->freq[i] < 1) &&
                    !(mapData->physicalPos[i] < scaffold->physicalPos[0]) &&
                    !(mapData->physicalPos[i] > scaffold->physicalPos[scaffold->nloci - 1]) &&
                    !(mapData->physicalPos[i] > scaffold->centroStart && mapData->physicalPos[i] < scaffold->centroEnd)) {
                newLoci++;
            }
        }
    }

    GenoLikeData *GLData2 = initGLData(GLData->nind, newLoci);
    int index = 0;
    for (int i = 0; i < freqData->nloci; i++)
    {
        if ((freqData->freq[i] > 0 && freqData->freq[i] < 1) &&
                !(mapData->physicalPos[i] < scaffold->physicalPos[0]) &&
                !(mapData->physicalPos[i] > scaffold->physicalPos[scaffold->nloci - 1]) &&
                !(mapData->physicalPos[i] > scaffold->centroStart && mapData->physicalPos[i] < scaffold->centroEnd))
        {
            for (int j = 0; j < GLData->nind; j++)
            {
                GLData2->data[index][j] = GLData->data[i][j];
            }
            index++;
        }
    }

    return GLData2;
}

FreqData *filterMonomorphicAndOOBSites(FreqData *freqData, MapData *mapData, GenMapScaffold *scaffold, int &newLoci) {
    if (newLoci <= 0) {
        newLoci = 0;
        for (int i = 0; i < freqData->nloci; i++) {
            if ((freqData->freq[i] > 0 && freqData->freq[i] < 1) &&
                    !(mapData->physicalPos[i] < scaffold->physicalPos[0]) &&
                    !(mapData->physicalPos[i] > scaffold->physicalPos[scaffold->nloci - 1]) &&
                    !(mapData->physicalPos[i] > scaffold->centroStart && mapData->physicalPos[i] < scaffold->centroEnd)) {
                newLoci++;
            }
        }
    }

    FreqData *freqData2 = initFreqData(newLoci);
    int index = 0;
    for (int i = 0; i < freqData->nloci; i++)
    {
        if ((freqData->freq[i] > 0 && freqData->freq[i] < 1) &&
                !(mapData->physicalPos[i] < scaffold->physicalPos[0]) &&
                !(mapData->physicalPos[i] > scaffold->physicalPos[scaffold->nloci - 1]) &&
                !(mapData->physicalPos[i] > scaffold->centroStart && mapData->physicalPos[i] < scaffold->centroEnd))
        {
            freqData2->freq[index] = freqData->freq[i];
            index++;
        }
    }

    return freqData2;
}

bool goodDouble(string str)
{
    string::iterator it;
    //int dashCount = 0;
    int decimalCount = 0;
    for (it = str.begin(); it != str.end(); it++)
    {
        if (!isdigit(*it) && *it != '.' && *it != '-') return 0;
        if (*it == '.') decimalCount++;
        if (*it == '-' && it != str.begin()) return 0;
        if (/*dashCount > 1 || */decimalCount > 1) return 0;
    }
    return 1;
}

FreqData *calcFreqData(HapData *hapData, int nresample, const gsl_rng *r)
{
    FreqData *freqData = initFreqData(hapData->nloci);
    double total, freq, count;

    for (int locus = 0; locus < hapData->nloci; locus++)
    {
        total = 0;
        count = 0;
        for (int ind = 0; ind < hapData->nind; ind++)
        {
            if (hapData->data[locus][ind] != -9)
            {
                count += hapData->data[locus][ind];
                total += 2;
            }
        }
        freq = count / total;
        if (nresample == 0) freqData->freq[locus] = freq;
        else
        {
            count = 0;
            for (int i = 0; i < nresample; i++)
            {
                if (gsl_rng_uniform(r) <= freq) count++;
            }
            freqData->freq[locus] = count / nresample;
        }

    }
    return freqData;
}

vector< FreqData * > *calcFreqData2(vector< HapData * > *hapDataByChr, int nresample)
{
    const gsl_rng_type *T;
    gsl_rng *r;
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r, time(NULL));

    vector< FreqData * > *freqDataByChr = new vector< FreqData * >;
    for (unsigned int chr = 0; chr < hapDataByChr->size(); chr++)
    {
        FreqData *data = calcFreqData(hapDataByChr->at(chr), nresample, r);
        freqDataByChr->push_back(data);
    }

    gsl_rng_free(r);

    return freqDataByChr;
}

//allocates the arrays and populates them with MISSING
FreqData *initFreqData(int nloci)
{
    if (nloci < 1)
    {
        LOG.err("ERROR: number of loci must be positive: ", nloci);
        throw 0;
    }

    FreqData *data = new FreqData;
    data->nloci = nloci;
    data->freq = new double[nloci];
    //data->allele = new string[nloci];

    for (int locus = 0; locus < nloci; locus++)
    {
        data->freq[locus] = MISSING;
        //data->allele[locus] = " ";
    }

    return data;
}

void releaseFreqData(FreqData *data)
{
    if (data == NULL) return;
    data->nloci = -9;
    delete [] data->freq;
    //delete [] data->allele;
    delete data;
    data = NULL;
    return;
}

void releaseFreqData(vector< FreqData * > *freqDataByChr)
{
    for (unsigned int chr = 0; chr < freqDataByChr->size(); chr++)
    {
        releaseFreqData(freqDataByChr->at(chr));
    }
    freqDataByChr->clear();
    delete freqDataByChr;
    return;
}

void writeFreqData(string freqOutfile, string popName,
                   vector< FreqData * > *freqDataByChr,
                   vector< MapData * > *mapDataByChr,
                   IndData *indData)
{
    freqOutfile += ".gz";
    ogzstream fout;
    fout.open(freqOutfile.c_str());

    if (fout.fail())
    {
        //cerr << "ERROR: Failed to open " << freqOutfile << " for writing.\n";
        LOG.err("ERROR: Failed to open", freqOutfile);
        throw 0;
    }

    fout << "CHR\tSNP\tPOS\tALLELE\tFREQ\n";

    for (unsigned int chr = 0; chr < mapDataByChr->size(); chr++)
    {
        for (int locus = 0; locus < mapDataByChr->at(chr)->nloci; locus++)
        {
            fout << mapDataByChr->at(chr)->chr << "\t"
                 << mapDataByChr->at(chr)->locusName[locus] << "\t"
                 << int(mapDataByChr->at(chr)->physicalPos[locus]) << "\t"
                 << mapDataByChr->at(chr)->allele[locus] << "\t"
                 << freqDataByChr->at(chr)->freq[locus] << "\n";
        }
    }
    cout << "Wrote allele frequency data to " << freqOutfile << endl;
    fout.close();
    return;
}

vector< FreqData * > *readFreqData(string freqfile, string popName,
                                   vector< MapData * > *mapDataByChr)
{
    //scan file for format integrity
    int minCols = 5;
    int expectedRows = 1;

    for (unsigned int chr = 0; chr < mapDataByChr->size(); chr++)
    {
        expectedRows += mapDataByChr->at(chr)->nloci;
    }

    igzstream fin;
    fin.open(freqfile.c_str());
    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << freqfile << " for reading.\n";
        LOG.err("ERROR: Failed to open", freqfile);
        throw 0;
    }
    cerr << "Reading " << freqfile << "\n";
    string line;
    int currentRows = 0;
    int currentCols = 0;
    int previousCols = -1;
    while (getline(fin, line))
    {
        currentRows++;
        currentCols = countFields(line);
        if (currentCols < minCols)
        {
            cerr << "ERROR: Found " << currentCols << " in " << freqfile
                 << " on line " << currentRows << " but expected at least "
                 << minCols << ".\n";
            LOG.err("ERROR: Found", currentCols, false);
            LOG.err(" in", freqfile, false);
            LOG.err(" on line", currentRows, false);
            LOG.err(" but expected at least", minCols);
            throw 0;
        }
        if (currentCols != previousCols && previousCols != -1)
        {
            cerr << "ERROR: " << freqfile << " has differing number of columns across rows.\n";
            LOG.err("ERROR: Differing number of columns across rows found in", freqfile);
            throw 0;
        }
        previousCols = currentCols;
    }

    if (currentRows != expectedRows)
    {
        cerr << "ERROR: " << freqfile << " has " << currentRows << " rows but expected "
             << expectedRows << ".\n";
        LOG.err("ERROR:", freqfile, false);
        LOG.err(" has", currentRows, false);
        LOG.err(" rows but expected", expectedRows);
        throw 0;
    }

    fin.close();
    fin.clear();

    //allocate
    vector< FreqData * > *freqDataByChr = new vector< FreqData * >;
    for (unsigned int chr = 0; chr < mapDataByChr->size(); chr++)
    {
        FreqData *data = initFreqData(mapDataByChr->at(chr)->nloci);
        freqDataByChr->push_back(data);
    }

    fin.open(freqfile.c_str());
    if (fin.fail())
    {
        LOG.err("ERROR: Failed to open", freqfile);
        throw 0;
    }

    string header, junk;
    //int popsFound = 0;
    getline(fin, header);
    //int headerSize = countFields(header) - 2;
    stringstream ss;
    //ss.str(header);
    //ss >> junk >> junk >> junk >> junk;
    //int popLocation = -1;
    //string popStr;
    //for (int i = 0; i < headerSize; i++)
    //{
    //    ss >> popStr;
    //    if (popStr.compare(popName) == 0) popLocation = i;
    //}

    //if (popLocation < 0) {
    //    cerr << "ERROR: Could not find " << popName << " in " << freqfile << endl;
    //    LOG.err("ERROR: Could not find", popName, false);
    //    LOG.err(" in", freqfile);
    //    throw 0;
    //}

    string locusID;//, allele;
    char allele;
    string chromosome;
    double position;
    for (unsigned int chr = 0; chr < mapDataByChr->size(); chr++)
    {
        for (int locus = 0; locus < mapDataByChr->at(chr)->nloci; locus++)
        {
            fin >> chromosome >> locusID >> position >> allele >> freqDataByChr->at(chr)->freq[locus];
            if (mapDataByChr->at(chr)->locusName[locus].compare(locusID) != 0){
                LOG.err("ERROR: Loci appear out of order in", freqfile);
                throw 0;
            }
            else mapDataByChr->at(chr)->allele[locus] = allele;
        }
    }

    fin.close();

    return freqDataByChr;
}

//allocates the arrays and populates them with MISSING or "--" depending on type
MapData *initMapData(int nloci)
{
    if (nloci < 1)
    {
        cerr << "ERROR: number of loci (" << nloci << ") must be positive.\n";
        LOG.err("ERROR: number of loci must be positive:", nloci);
        throw 0;
    }

    MapData *data = new MapData;
    data->nloci = nloci;
    data->locusName = new string[nloci];
    data->physicalPos = new double[nloci];
    data->geneticPos = new double[nloci];
    data->allele = new char[nloci];
    //data->allele0 = new char[nloci];
    data->chr = "--";

    for (int locus = 0; locus < nloci; locus++)
    {
        data->locusName[locus] = "--";
        data->physicalPos[locus] = MISSING;
        data->geneticPos[locus] = MISSING;
        data->allele[locus] = '-';
        //data->allele0[locus] = '-';
    }

    return data;
}

void releaseMapData(MapData *data)
{
    if (data == NULL) return;
    data->nloci = -9;
    delete [] data->locusName;
    delete [] data->physicalPos;
    delete [] data->geneticPos;
    delete [] data->allele;
    //delete [] data->allele0;
    delete data;
    data = NULL;
    return;
}

void releaseMapData(vector< MapData * > *mapDataByChr)
{
    for (unsigned int i = 0; i < mapDataByChr->size(); i++)
    {
        releaseMapData(mapDataByChr->at(i));
    }
    mapDataByChr->clear();
    delete mapDataByChr;
    return;
}

void releaseIndData(vector< IndData * > *indDataByPop)
{
    for (unsigned int i = 0; i < indDataByPop->size(); i++)
    {
        releaseIndData(indDataByPop->at(i));
    }
    indDataByPop->clear();
    delete indDataByPop;
    return;
}

void releaseIndData(IndData *data)
{
    delete [] data->indID;
    delete data;
    return;
}

vector< HapData * > *readTPEDHapData3(string filename,
                                      int expectedLoci,
                                      int expectedInd,
                                      char TPED_MISSING,
                                      vector< MapData * > *mapDataByChr, bool PHASED)
{
    int expectedHaps = 2 * expectedInd;
    igzstream fin;
    //cerr << "Checking " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        LOG.err("ERROR: Failed to open", filename);
        throw 0;
    }

    cout << "Loading genotypes from " << filename << "\n";

    //int fileStart = fin.tellg();
    string line;
    int nhaps = -1;
    int nloci = 0;
    while (getline(fin, line))
    {
        nloci++;
        nhaps = countFields(line);
        //cout << "nhaps: " << current_nhaps << endl;
        if (nhaps != expectedHaps + 4)
        {
            cerr << "ERROR: line " << nloci << " of " << filename << " has " << nhaps
                 << " columns, but expected " << expectedHaps + 4 << ".\n";
            LOG.err("ERROR: line", nloci, false);
            LOG.err(" of", filename, false);
            LOG.err(" has", nhaps, false);
            LOG.err(" columns, but expected", expectedHaps);
            throw 0;
        }
    }
    if (nloci != expectedLoci)
    {
        cerr << "ERROR: " << filename << " has " << nloci
             << " loci, but expected " << expectedLoci << ".\n";
        LOG.err("ERROR:", filename, false);
        LOG.err(" has", nloci, false);
        LOG.err(" loci, but expected", expectedLoci);
        throw 0;
    }

    fin.close();
    fin.clear();

    fin.open(filename.c_str());

    if (fin.fail())
    {
        LOG.err("ERROR: Failed to open", filename);
        throw 0;
    }

    vector< HapData * > *hapDataByChr = new vector< HapData * >;
    for (unsigned int chr = 0; chr < mapDataByChr->size(); chr++)
    {
        HapData *data = initHapData(expectedInd, mapDataByChr->at(chr)->nloci, PHASED);
        hapDataByChr->push_back(data);
    }

    //string alleleStr1, alleleStr2;
    char alleleStr1, alleleStr2;
    string junk;//, oneAllele;
    char oneAllele;
    short allele;
    //For each chromosome
    for (unsigned int chr = 0; chr < mapDataByChr->size(); chr++)
    {
        //For each locus on the chromosome
        for (int locus = 0; locus < mapDataByChr->at(chr)->nloci; locus++)
        {
            oneAllele = mapDataByChr->at(chr)->allele[locus];
            fin >> junk;
            fin >> junk;
            fin >> junk;
            fin >> junk;
            for (int ind = 0; ind < expectedInd; ind++)
            {
                fin >> alleleStr1 >> alleleStr2;
                allele = 0;

                if (alleleStr1 == TPED_MISSING) allele += -9;
                else if (alleleStr1 == oneAllele) allele += 1;

                if (alleleStr2 == TPED_MISSING) allele += -9;
                else if (alleleStr2 == oneAllele) allele += 1;

                if (allele < 0) hapDataByChr->at(chr)->data[locus][ind] = -9;
                else hapDataByChr->at(chr)->data[locus][ind] = allele;

                if(PHASED) hapDataByChr->at(chr)->firstCopy[locus][ind] = (alleleStr1 == oneAllele);
            }
        }
    }

    fin.close();
    return hapDataByChr;
}

vector< GenoLikeData * > *readTGLSData(string filename,
                                       int expectedLoci,
                                       int expectedInd,
                                       vector< MapData * > *mapDataByChr,
                                       string GL_TYPE)
{
    //int expectedHaps = 2 * expectedInd;
    igzstream fin;
    //cerr << "Checking " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        //cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        LOG.err("ERROR: Failed to open", filename);
        throw 0;
    }

    cout << "Loading genotype likelihoods from " << filename << "\n";

    //int fileStart = fin.tellg();
    string line;
    int nhaps = -1;
    int nloci = 0;
    while (getline(fin, line))
    {
        nloci++;
        nhaps = countFields(line);
        //cout << "nhaps: " << current_nhaps << endl;
        if (nhaps != expectedInd + 4)
        {
            /*
            cerr << "ERROR: line " << nloci << " of " << filename << " has " << nhaps
                 << " columns, but expected " << expectedInd + 4 << ".\n";
            */
            LOG.err("ERROR: line", nloci, false);
            LOG.err(" of", filename, false);
            LOG.err(" has", nhaps, false);
            LOG.err(" columns, but expected", expectedInd);
            throw 0;
        }
    }
    if (nloci != expectedLoci)
    {
        /*
        cerr << "ERROR: " << filename << " has " << nloci
             << " loci, but expected " << expectedLoci << ".\n";
        */
        LOG.err("ERROR:", filename, false);
        LOG.err(" has", nloci, false);
        LOG.err(" loci, but expected", expectedLoci);
        throw 0;
    }

    fin.close();
    fin.clear();

    fin.open(filename.c_str());

    if (fin.fail())
    {
        //cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        LOG.err("ERROR: Failed to open", filename);
        throw 0;
    }

    vector< GenoLikeData * > *GLDataByChr = new vector< GenoLikeData * >;
    for (unsigned int chr = 0; chr < mapDataByChr->size(); chr++)
    {
        GenoLikeData *data = initGLData(expectedInd, mapDataByChr->at(chr)->nloci);
        GLDataByChr->push_back(data);
    }

    //string alleleStr1, alleleStr2;
    string junk;//, oneAllele;
    double gl;
    //For each chromosome
    for (unsigned int chr = 0; chr < mapDataByChr->size(); chr++)
    {
        //For each locus on the chromosome
        for (int locus = 0; locus < mapDataByChr->at(chr)->nloci; locus++)
        {
            fin >> junk;
            fin >> junk;
            fin >> junk;
            fin >> junk;
            for (int ind = 0; ind < expectedInd; ind++)
            {
                fin >> gl;
                if (GL_TYPE.compare("GL") == 0) {
                }
                else if (GL_TYPE.compare("PL") == 0) {
                    gl = -gl / 10;
                }
                else {
                    LOG.err("ERROR: This error should never get triggered. Something really bad happened.");
                }
                gl = ( gl > -10 ) ? gl : -10;
                GLDataByChr->at(chr)->data[locus][ind] = 1 - pow(10,gl);

            }
        }
    }

    fin.close();
    return GLDataByChr;
}

void releaseHapData(vector< HapData * > *hapDataByChr)
{
    for (unsigned int chr = 0; chr < hapDataByChr->size(); chr++)
    {
        releaseHapData(hapDataByChr->at(chr));
    }
    hapDataByChr->clear();
    delete hapDataByChr;
}

string getPost(int num)
{
    string post;
    if (num == 1) post = "st";
    else if (num == 2) post = "nd";
    else if (num == 3) post = "rd";
    else post = "th";
    return post;
}

WinData *initWinData(unsigned int nind, unsigned int nloci)
{
    if (nind < 1 || nloci < 1)
    {
        cerr << "ERROR: Can't allocate WinData object.  Number of individuals (" << nind
             << ") and number of loci (" << nloci
             << ") must be positive.\n";
        LOG.err("ERROR: Can't allocate WinData object.  Number of individuals (", int(nind), false);
        LOG.err(" ) and number of loci (", int(nloci), false);
        LOG.err(" ) must be positive.");

        throw 0;
    }

    WinData *data = new WinData;
    data->nind = nind;
    data->nloci = nloci;
    //data->nmiss = 0;

    data->data = new double*[nind];
    for (unsigned int i = 0; i < nind; i++)
    {
        data->data[i] = new double[nloci];
        for (unsigned int j = 0; j < nloci; j++)
        {
            data->data[i][j] = MISSING;
        }
    }

    return data;
}

void releaseWinData(WinData *data)
{
    if (data == NULL) return;
    for (int i = 0; i < data->nind; i++)
    {
        delete [] data->data[i];
    }

    delete [] data->data;

    data->data = NULL;
    data->nind = -9;
    data->nloci = -9;
    //data->nmiss = -9;
    delete data;
    data = NULL;
    return;
}

void releaseWinData(vector< WinData * > *winDataByChr)
{
    for (unsigned int chr = 0; chr < winDataByChr->size(); chr++)
    {
        releaseWinData(winDataByChr->at(chr));
    }
    winDataByChr->clear();
    delete winDataByChr;
}

vector< vector< WinData * >* > *initWinData(vector< MapData * > *mapDataByChr,
        vector< IndData * > *indDataByPop)
{
    vector< vector< WinData * >* > *winDataByPopByChr = new vector< vector< WinData * >* >;

    for (unsigned int pop = 0; pop < indDataByPop->size(); pop++)
    {
        int nind = indDataByPop->at(pop)->nind;
        vector< WinData * > *winDataByChr = new vector< WinData * >;
        for (unsigned int chr = 0; chr < mapDataByChr->size(); chr++)
        {
            int nloci = mapDataByChr->at(chr)->nloci;
            WinData *data = initWinData(nind, nloci);
            winDataByChr->push_back(data);
        }
        winDataByPopByChr->push_back(winDataByChr);
    }

    return winDataByPopByChr;
}

vector< WinData * > *initWinData(vector< MapData * > *mapDataByChr, int nind)
{
    //int nind = indData->nind;
    vector< WinData * > *winDataByChr = new vector< WinData * >;
    for (unsigned int chr = 0; chr < mapDataByChr->size(); chr++)
    {
        int nloci = mapDataByChr->at(chr)->nloci;
        WinData *data = initWinData(nind, nloci);
        winDataByChr->push_back(data);
    }

    return winDataByChr;
}

void writeWinData(vector< WinData * > *winDataByChr,
                  IndData *indData,
                  vector< MapData * > *mapDataByChr,
                  string outfile)
{
    ogzstream fout;
    int numChr = mapDataByChr->size();
    string popName = indData->pop;

    for (int chr = 0; chr < numChr; chr++)
    {
        string rawWinOutfile = outfile;
        rawWinOutfile += ".";
        rawWinOutfile += popName;
        rawWinOutfile += ".";
        rawWinOutfile += mapDataByChr->at(chr)->chr;
        rawWinOutfile += ".raw.lod.windows.gz";

        fout.open(rawWinOutfile.c_str());
        if (fout.fail())
        {
            cerr << "ERROR: Failed to open " << rawWinOutfile << " for writing.\n";
            LOG.err("ERROR: Failed to open", rawWinOutfile);
            throw - 1;
        }

        WinData *winData = winDataByChr->at(chr);

        for (int ind = 0; ind < winData->nind; ind++)
        {
            for (int locus = 0; locus < winData->nloci; locus++)
            {
                if (winData->data[ind][locus] == MISSING) fout << "NA";
                else fout << winData->data[ind][locus];
                if (locus < winData->nloci - 1) fout << " ";
            }
            fout << endl;
        }
        cerr << "Wrote " << rawWinOutfile << "\n";
        fout.close();
    }

    return;
}

HapData *initHapData(unsigned int nind, unsigned int nloci, bool PHASED)
{
    if (nind < 1 || nloci < 1)
    {
        LOG.err("ERROR: Can not allocate HapData object. Number of haplotypes (", int(nind), false);
        LOG.err(" ) and number of loci (", int(nloci), false);
        LOG.err(" ) must be positive.");
        throw 0;
    }

    HapData *data = new HapData;
    data->nind = nind;
    data->nloci = nloci;

    data->data = new short*[nloci];
    if(PHASED) data->firstCopy = new bool*[nloci];
    else data->firstCopy = NULL;
    for (unsigned int i = 0; i < nloci; i++)
    {
        data->data[i] = new short[nind];
        if(PHASED) data->firstCopy[i] = new bool[nind];
        for (unsigned int j = 0; j < nind; j++)
        {
            data->data[i][j] = MISSING;
            if(PHASED) data->firstCopy[i][j] = 0;
        }
    }

    return data;
}

void releaseHapData(HapData *data)
{
    if (data == NULL) return;

    for (int i = 0; i < data->nloci; i++)
    {
        if(data->firstCopy != NULL) delete [] data->firstCopy[i];
        delete [] data->data[i];
    }

    delete [] data->data;

    if(data->firstCopy != NULL) delete [] data->firstCopy;

    data->data = NULL;
    data->firstCopy = NULL;
    data->nind = -9;
    data->nloci = -9;
    delete data;
    data = NULL;
    return;
}

GenoLikeData *initGLData(unsigned int nind, unsigned int nloci) {
    if (nind < 1 || nloci < 1)
    {
        LOG.err("ERROR: Can not allocate GenoLikeData object. Number of individuals (", int(nind), false);
        LOG.err(" ) and number of loci (", int(nloci), false);
        LOG.err(" ) must be positive.");
        throw 0;
    }

    GenoLikeData *data = new GenoLikeData;
    data->nind = nind;
    data->nloci = nloci;

    data->data = new double*[nloci];
    for (unsigned int i = 0; i < nloci; i++)
    {
        data->data[i] = new double[nind];
        for (unsigned int j = 0; j < nind; j++)
        {
            data->data[i][j] = 1;
        }
    }

    return data;
}

void releaseGLData(GenoLikeData *data) {
    if (data == NULL) return;
    for (int i = 0; i < data->nloci; i++)
    {
        delete [] data->data[i];
    }

    delete [] data->data;

    data->data = NULL;
    data->nind = -9;
    data->nloci = -9;
    delete data;
    data = NULL;
    return;
}

void releaseGLData(vector< GenoLikeData * > *GLDataByChr) {
    for (unsigned int chr = 0; chr < GLDataByChr->size(); chr++)
    {
        releaseGLData(GLDataByChr->at(chr));
    }
    GLDataByChr->clear();
    delete GLDataByChr;
}

int countFields(const string &str)
{
    string::const_iterator it;
    int result;
    int numFields = 0;
    int seenChar = 0;
    for (it = str.begin() ; it < str.end(); it++)
    {
        result = isspace(*it);
        if (result == 0 && seenChar == 0)
        {
            numFields++;
            seenChar = 1;
        }
        else if (result != 0)
        {
            seenChar = 0;
        }
    }
    return numFields;
}

vector< int_pair_t > *scanTPEDMapData(string filename, int &numLoci, int &numCols)
{
    igzstream fin;
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        LOG.err("ERROR: Failed to open", filename);
        throw 0;
    }

    cout << "Reading " << filename << "\n";

    vector< int_pair_t > *chrStartStop = new vector< int_pair_t >;
    stringstream ss;
    string line;
    int nloci = 0;
    int index;
    int prevCols = -1;
    int currCols = -1;
    string emptyChr = "_nochr";
    string prevChr = emptyChr;
    string currChr = prevChr;
    int_pair_t currChrCoordinates;
    while (getline(fin, line))
    {
        nloci++;
        index = nloci - 1;
        currCols = countFields(line);
        if (currCols != prevCols && nloci > 1)
        {
            cerr << "ERROR: line " << nloci << " of " << filename
                 << " has different number of columns from earlier lines.\n";
            LOG.err("ERROR: line", nloci, false);
            LOG.err(" has different number of columns from earlier lines in ", filename);
            throw 0;
        }

        ss.str(line);
        ss >> currChr;
        if (prevChr.compare(emptyChr) == 0 && index == 0)
        {
            prevChr = currChr;
            currChrCoordinates.first = index;
        }

        if (currChr.compare(prevChr) != 0)
        {
            currChrCoordinates.second = index - 1;
            chrStartStop->push_back(currChrCoordinates);
            currChrCoordinates.first = index;
            prevChr = currChr;
        }
        ss.clear();
        prevCols = currCols;
    }

    fin.close();

    numLoci = nloci;

    currChrCoordinates.second = index;
    chrStartStop->push_back(currChrCoordinates);

    numCols = currCols;

    return chrStartStop;
}

string lc(string str) {
    char c[2] = {' ', '\0'};
    for (unsigned int i = 0; i < str.size(); i++) {
        c[0] = char(tolower(str[i]));
        str.replace(i, 1, c);
    }
    return str;
}

string checkChrName(string chr) {
    if (chr[0] != 'c') {
        chr = "chr" + chr;
    }
    return chr;
}

vector< MapData * > *readTPEDMapData(string filename, int numCols, vector< int_pair_t > *chrCoordList, char TPED_MISSING)
{
    vector< MapData * > *mapDataByChr = new vector< MapData * >;

    igzstream fin;
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        LOG.err("ERROR: Failed to open", filename);
        throw 0;
    }

    cout << "Loading map from " << filename << "\n";

    string junk;
    char oneAllele;
    int count = 4;
    //For each chromosome
    for (unsigned int i = 0; i < chrCoordList->size(); i++)
    {
        int size = chrCoordList->at(i).second - chrCoordList->at(i).first + 1;
        MapData *data = initMapData(size);
        for (int locus = 0; locus < size; locus++)
        {
            fin >> data->chr;
            fin >> data->locusName[locus];
            fin >> data->geneticPos[locus];
            fin >> data->physicalPos[locus];

            //assign 'one' allele coding
            oneAllele = TPED_MISSING;
            count = 4;
            while (oneAllele == TPED_MISSING && count < numCols) {
                fin >> junk;
                if (junk[0] != TPED_MISSING) oneAllele = junk.c_str()[0];
                count++;
            }

            if (count >= numCols && oneAllele == TPED_MISSING) {
                cerr << "ERROR: locus " << data->locusName[locus] << " appears to have no data.\n";
                LOG.err("ERROR: Locus appears to have no data:", data->locusName[locus]);
                throw 0;
            }

            data->allele[locus] = oneAllele;
            if (count < numCols) getline(fin, junk);
        }
        data->chr = lc(data->chr);
        data->chr = checkChrName(data->chr);
        cout << size << " loci on chromosome " << data->chr << endl;
        mapDataByChr->push_back(data);
    }

    fin.close();

    return mapDataByChr;
}


void scanIndData3(string filename, int &numInd, string &popName) {
    igzstream fin;
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        LOG.err("ERROR: Failed to open", filename);
        throw 0;
    }

    cout << "Reading " << filename << "\n";

    map<string, int> indList;

    string line;
    int nind = 0;
    int min_cols = 2;
    int current_cols = 0;
    string pop, ind;
    stringstream ss;
    while (getline(fin, line))
    {
        nind++;
        current_cols = countFields(line);
        if (current_cols < min_cols)
        {
            cerr << "ERROR: line " << nind << " of " << filename << " has " << current_cols
                 << ", but expected at least " << min_cols << ".\n";
            LOG.err("ERROR: Line", nind, false);
            LOG.err(" of", filename, false);
            LOG.err(" has", current_cols, false);
            LOG.err(", but expected at least", min_cols);
            throw 0;
        }
        //stringstream ss;
        ss.str(line);
        ss >> pop >> ind;
        if (indList.count(ind) > 0)
        {
            cerr << "ERROR: Found duplicate individual ID (" << ind << ") in " << filename << endl;
            LOG.err("ERROR: Found duplicate individual ID ( ", ind, false);
            LOG.err(" ) in", filename);
            throw 0;
        }
        else indList[ind] = 1;

        //Single population only, check to see if more than
        if (nind == 1) {
            popName = pop;
        }
        else if (pop.compare(popName) != 0) {
            cerr << "ERROR: Found multiple population IDs (" << pop << ", " << popName << ") in " << filename << endl;
            cerr << "\tGARLIC must be given only a single population at a time.\n";
            LOG.err("ERROR: Found multiple population IDs ( ", pop, false);
            LOG.err(",", popName, false);
            LOG.err(" ) in", filename);
            throw 0;
        }
        ss.clear();
    }

    fin.close();

    numInd = nind;

    return;
}


IndData *readIndData3(string filename, int numInd)
{
    igzstream fin;
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        LOG.err("ERROR: Failed to open", filename);
        throw 0;
    }

    cout << "Loading individual IDs\n";

    IndData *indData = initIndData(numInd);

    string line, pop, ind;
    stringstream ss;
    for (int i = 0; i < numInd; i++)
    {
        getline(fin, line);
        ss.str(line);
        ss >> pop >> ind;
        indData->indID[i] = ind;
        ss.clear();
    }
    fin.close();

    indData->pop = pop;
    return indData;
}

IndData *initIndData(int nind)
{
    if (nind < 1)
    {
        cerr << "ERROR: number of individuals (" << nind << ") must be positive.\n";
        LOG.err("ERROR: Number of individuals must be positive:", nind);
        throw 0;
    }

    IndData *data = new IndData;
    data->nind = nind;
    data->indID = new string[nind];

    for (int ind = 0; ind < nind; ind++)
    {
        data->indID[ind] = "--";
    }

    return data;
}

DoubleData *initDoubleData(int n)
{
    DoubleData *data = new DoubleData;

    data->size = n;
    data->data = new double[n];

    return data;
}

DoubleData *convertWinData2DoubleData(vector< WinData * > *winDataByChr, int step)
{
    //int nmiss = 0;
    //int ncols = 0;
    //int nrows = 0;
    int size = 0;
    DoubleData *rawWinData;
    for (unsigned int chr = 0; chr < winDataByChr->size(); chr++)
    {
        for (int ind = 0; ind < winDataByChr->at(chr)->nind; ind++)
        {
            for (int locus = 0; locus < winDataByChr->at(chr)->nloci; locus+=step)
            {
                //if (winDataByChr->at(chr)->data[ind][locus] == MISSING) nmiss++;
                if (winDataByChr->at(chr)->data[ind][locus] != MISSING) size++;
            }
        }

        //ncols += winDataByChr->at(chr)->nloci;
        //nrows = winDataByChr->at(chr)->nind;
    }
    //rawWinData = initDoubleData(ncols * nrows - nmiss);
    rawWinData = initDoubleData(size);
    
    int i = 0;
    for (unsigned int chr = 0; chr < winDataByChr->size(); chr++)
    {
        for (int ind = 0; ind < winDataByChr->at(chr)->nind; ind++)
        {
            for (int locus = 0; locus < winDataByChr->at(chr)->nloci; locus+=step)
            {
                if (winDataByChr->at(chr)->data[ind][locus] != MISSING)
                {
                    rawWinData->data[i] = winDataByChr->at(chr)->data[ind][locus];
                    i++;
                }
            }
        }
    }

    return rawWinData;
}

DoubleData *convertSubsetWinData2DoubleData(vector< WinData * > *winDataByChr, IndData *indData, int subsample, int step)
{
    const gsl_rng_type *T;
    gsl_rng *r;
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r, time(NULL));

    //to hold the indicies of the randomly selected individuals
    int nind = winDataByChr->at(0)->nind;
    int *randInd;
    if (subsample >= nind)
    {
        randInd = new int[nind];
        for (int i = 0; i < nind; i++) randInd[i] = i;
    }
    else
    {
        int* indIndex = new int[nind];
        for (int i = 0; i < nind; i++) indIndex[i] = i;
        randInd = new int[subsample];
        gsl_ran_choose(r, randInd, subsample, indIndex, nind, sizeof(int));
        delete [] indIndex;
        nind = subsample;
    }

    LOG.logn("Individuals used for KDE: ");
    for (int ind = 0; ind < nind; ind++) {
        LOG.logn(indData->indID[randInd[ind]]);
        LOG.logn(" ");
    }
    LOG.logn("\n");

    DoubleData *rawWinData;
    //int nmiss = 0;
    //int ncols = 0;
    //int nrows = 0;
    int size = 0;
    for (unsigned int chr = 0; chr < winDataByChr->size(); chr++)
    {
        for (int ind = 0; ind < nind; ind++)
        {
            for (int locus = 0; locus < winDataByChr->at(chr)->nloci; locus+=step)
            {
                //if (winDataByChr->at(chr)->data[randInd[ind]][locus] == MISSING) nmiss++;
                if (winDataByChr->at(chr)->data[randInd[ind]][locus] != MISSING) size++;
            }
        }

        //ncols += winDataByChr->at(chr)->nloci;
        //nrows = nind;
    }
    //rawWinData = initDoubleData(ncols * nrows - nmiss);
    rawWinData = initDoubleData(size);

    int i = 0;
    for (unsigned int chr = 0; chr < winDataByChr->size(); chr++)
    {
        for (int ind = 0; ind < nind; ind++)
        {
            for (int locus = 0; locus < winDataByChr->at(chr)->nloci; locus+=step)
            {
                if (winDataByChr->at(chr)->data[randInd[ind]][locus] != MISSING)
                {
                    rawWinData->data[i] = winDataByChr->at(chr)->data[randInd[ind]][locus];
                    i++;
                }
            }
        }
    }

    gsl_rng_free(r);


    delete [] randInd;

    return rawWinData;
}

void releaseDoubleData(DoubleData *data)
{
    delete [] data->data;
    delete data;
    return;
}

void releaseDoubleData(vector < DoubleData * > *rawWinDataByPop)
{
    for (unsigned int pop = 0; pop < rawWinDataByPop->size(); pop++)
    {
        releaseDoubleData(rawWinDataByPop->at(pop));
    }
    rawWinDataByPop->clear();
    delete rawWinDataByPop;
    rawWinDataByPop = NULL;
    return;
}

void subsetData(vector< HapData * > *hapDataByChr,
                vector< GenoLikeData *> *GLDataByChr,
                IndData *indData,
                vector< HapData * > **subsetHapDataByChr,
                vector< GenoLikeData *> **subsetGLDataByChr,
                IndData **subsetIndData,
                int subsample, bool USE_GL, bool PHASED)
{
    const gsl_rng_type *T;
    gsl_rng *r;
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r, time(NULL));

    int nind = hapDataByChr->at(0)->nind;
    int *randInd;
    if (subsample >= nind)
    {
        randInd = new int[nind];
        for (int i = 0; i < nind; i++) randInd[i] = i;
    }
    else
    {
        int* indIndex = new int[nind];
        for (int i = 0; i < nind; i++) indIndex[i] = i;
        randInd = new int[subsample];
        gsl_ran_choose(r, randInd, subsample, indIndex, nind, sizeof(int));
        delete [] indIndex;
        nind = subsample;
    }

    IndData *newIndData = initIndData(nind);
    newIndData->pop = indData->pop;

    for (int ind = 0; ind < nind; ind++) {
        newIndData->indID[ind] = indData->indID[randInd[ind]];
    }
    LOG.loga("Individuals used for KDE:", newIndData->indID, nind);

    vector< HapData * > *newHapDataByChr = new vector< HapData * >;
    vector< GenoLikeData * > *newGLDataByChr;
    if (USE_GL) newGLDataByChr = new vector< GenoLikeData * >;

    int nchr = hapDataByChr->size();
    HapData *hapData;
    GenoLikeData *GLData;
    for (int chr = 0; chr < nchr; chr++)
    {
        int nloci = hapDataByChr->at(chr)->nloci;
        hapData = initHapData(nind, nloci, PHASED);
        GLData = initGLData(nind, nloci);

        for (int locus = 0; locus < nloci; locus++)
        {
            for (int ind = 0; ind < nind; ind++)
            {
                hapData->data[locus][ind] = hapDataByChr->at(chr)->data[locus][randInd[ind]];
                if (PHASED) hapData->firstCopy[locus][ind] = hapDataByChr->at(chr)->firstCopy[locus][randInd[ind]];
                if (USE_GL) GLData->data[locus][ind] = GLDataByChr->at(chr)->data[locus][randInd[ind]];
            }
        }
        newHapDataByChr->push_back(hapData);
        hapData = NULL;
        if (USE_GL) newGLDataByChr->push_back(GLData);
        GLData = NULL;
    }
    delete [] randInd;

    *(subsetHapDataByChr) = newHapDataByChr;
    *(subsetGLDataByChr) = newGLDataByChr;
    *(subsetIndData) = newIndData;
    gsl_rng_free(r);
    return;
}


