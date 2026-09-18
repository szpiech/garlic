#include "garlic-data.h"
#include <cmath>
#include <random>

static GarlicRNG *GARLIC_RNG = NULL;

void initRNG(unsigned long int seed)
{
    if (GARLIC_RNG != NULL) delete GARLIC_RNG;
    //GarlicRNG reproduces gsl_rng_mt19937 (which was gsl_rng_default) stream
    //for stream, so a seeded run gives exactly the results it did with GSL.
    GARLIC_RNG = new GarlicRNG(seed);
    return;
}

GarlicRNG *getRNG()
{
    if (GARLIC_RNG == NULL)
    {
        //Should not happen: main calls initRNG before any consumer runs.
        //Fall back to a reported random seed rather than a silent time(NULL).
        unsigned long int seed = drawRandomSeed();
        LOG.err("WARNING: RNG used before initialisation; seeding with", int(seed));
        initRNG(seed);
    }
    return GARLIC_RNG;
}

void freeRNG()
{
    if (GARLIC_RNG != NULL) delete GARLIC_RNG;
    GARLIC_RNG = NULL;
    return;
}

unsigned long int drawRandomSeed()
{
    std::random_device rd;
    //Must stay inside the range --seed can accept (param_t parses it with
    //atoi), otherwise the seed we report in the log cannot be replayed.
    //0 is the "choose for me" sentinel, so map into [1, INT_MAX].
    unsigned long int seed = (unsigned long int)(rd() % (unsigned int)(0x7FFFFFFF)) + 1UL;
    return seed;
}

//---- fast whitespace-delimited scanning -------------------------------------
//std::istream extraction (operator>>) constructs a sentry and performs locale
//lookups on every single token.  On TPED input that machinery, not the actual
//parsing, dominated runtime.  These helpers walk the line buffer directly.
static inline const char *skipSpace(const char *p, const char *e)
{
    while (p < e && (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\f' || *p == '\v')) ++p;
    return p;
}

static inline const char *tokenEnd(const char *p, const char *e)
{
    while (p < e && !(*p == ' ' || *p == '\t' || *p == '\r' || *p == '\f' || *p == '\v')) ++p;
    return p;
}

static double AUTO_OVERLAP_SLOPE = 6.375;
static double AUTO_OVERLAP_INTERCEPT = 63.888;

void setAutoOverlapCoef(double slope, double intercept) { AUTO_OVERLAP_SLOPE = slope; AUTO_OVERLAP_INTERCEPT = intercept; }

bool isSexChromosome(const string &chr)
{
    return (chr == "chrX" || chr == "chrY" || chr == "chr23" || chr == "chr24");
}

bool warnSexChromosomes(vector< MapData * > *mapDataByChr, IndData *indData)
{
    vector<string> found;
    for (unsigned int chr = 0; chr < mapDataByChr->size(); chr++)
        if (isSexChromosome(mapDataByChr->at(chr)->chr))
            found.push_back(mapDataByChr->at(chr)->chr);

    if (found.empty()) return false;

    int males = 0, females = 0, unknown = 0;
    for (int i = 0; i < indData->nind; i++)
    {
        if (indData->sex[i] == 1) males++;
        else if (indData->sex[i] == 2) females++;
        else unknown++;
    }

    string list;
    for (unsigned int i = 0; i < found.size(); i++)
    {
        if (i) list += ", ";
        list += found[i];
    }

    LOG.err("WARNING: the data contains sex chromosomes:", list);
    LOG.err("WARNING: a hemizygous male genotype is written as a homozygous call in a TPED, so");
    LOG.err("WARNING: it is indistinguishable from true autozygosity. Male X chromosomes will");
    LOG.err("WARNING: therefore be called as one run spanning the whole chromosome, and any");
    LOG.err("WARNING: FROH computed from that will be inflated. Use --autosomes-only to drop");
    LOG.err("WARNING: these chromosomes, or restrict to females using the sex column of the");
    LOG.err("WARNING: TFAM or of --pop.");

    //Reporting only the male count was misleading whenever sex was recorded
    //for PART of the cohort: 3 of 45 coded male with 42 unknown printed
    //"individuals coded male: 3", and 3 reads as the answer when up to 45
    //could be affected.  Sex is optional in a TFAM and in a --pop file, so
    //partial coverage is normal.  Report all three counts, and never state a
    //number of affected individuals that the metadata cannot support.
    if (males + females == 0)
    {
        LOG.err("WARNING: sex is not recorded for any of the", indData->nind, false);
        LOG.err(" individuals, so the number affected cannot be determined.");
    }
    else if (unknown == 0)
    {
        LOG.err("WARNING: sex recorded for all", indData->nind, false);
        LOG.err(" individuals:", males, false);
        LOG.err(" male,", females, false);
        LOG.err(" female.");
        LOG.err("WARNING: affected individuals:", males);
    }
    else
    {
        LOG.err("WARNING: sex recorded for", males + females, false);
        LOG.err(" of", indData->nind, false);
        LOG.err(" individuals:", males, false);
        LOG.err(" male,", females, false);
        LOG.err(" female,", unknown, false);
        LOG.err(" unknown.");
        LOG.err("WARNING: at least", males, false);
        LOG.err(" individuals are affected; with", unknown, false);
        LOG.err(" of unknown sex the true number");
        LOG.err("WARNING: cannot be determined and may be as high as", males + unknown, false);
        LOG.err(".");
    }

    return true;
}

int filterChromosomes(vector<string> &keep,
                      vector< MapData * > **mapDataByChr,
                      vector< HapData * > **hapDataByChr,
                      vector< FreqData * > **freqDataByChr,
                      vector< GenoLikeData * > **GLDataByChr,
                      bool USE_GL)
{
    vector<string> want;
    for (unsigned int i = 0; i < keep.size(); i++) want.push_back(checkChrName(keep[i]));

    vector< MapData * > *newMap = new vector< MapData * >;
    vector< HapData * > *newHap = new vector< HapData * >;
    vector< FreqData * > *newFreq = new vector< FreqData * >;
    vector< GenoLikeData * > *newGL = USE_GL ? new vector< GenoLikeData * > : NULL;

    vector<bool> found(want.size(), false);

    for (unsigned int chr = 0; chr < (*mapDataByChr)->size(); chr++)
    {
        string have = checkChrName((*mapDataByChr)->at(chr)->chr);
        int match = -1;
        for (unsigned int w = 0; w < want.size(); w++)
            if (want[w].compare(have) == 0) { match = int(w); break; }

        if (match >= 0)
        {
            found[match] = true;
            newMap->push_back((*mapDataByChr)->at(chr));
            newHap->push_back((*hapDataByChr)->at(chr));
            newFreq->push_back((*freqDataByChr)->at(chr));
            if (USE_GL) newGL->push_back((*GLDataByChr)->at(chr));
        }
        else
        {
            releaseMapData((*mapDataByChr)->at(chr));
            releaseHapData((*hapDataByChr)->at(chr));
            releaseFreqData((*freqDataByChr)->at(chr));
            if (USE_GL) releaseGLData((*GLDataByChr)->at(chr));
        }
    }

    for (unsigned int w = 0; w < want.size(); w++)
    {
        if (!found[w])
        {
            LOG.err("ERROR: --chr", want[w], false);
            LOG.err(" is not present in the data.");
            return -1;
        }
    }

    (*mapDataByChr)->clear();  delete *mapDataByChr;  *mapDataByChr = newMap;
    (*hapDataByChr)->clear();  delete *hapDataByChr;  *hapDataByChr = newHap;
    (*freqDataByChr)->clear(); delete *freqDataByChr; *freqDataByChr = newFreq;
    if (USE_GL) { (*GLDataByChr)->clear(); delete *GLDataByChr; *GLDataByChr = newGL; }

    return int(newMap->size());
}

double selectOverlapFrac(double variantDensity, int winsize){
    double frac = (AUTO_OVERLAP_SLOPE*log(variantDensity)+AUTO_OVERLAP_INTERCEPT)/100.0;
    if(frac > 1) frac = 1.0;
    if(frac <= 0) frac = 1.0/double(winsize);
    return frac;
}

void loadTPEDData(string tpedfile, int &numLoci, int &numInd,
                                   vector< HapData * > **hapDataByChr,
                                   vector< MapData * > **mapDataByChr,
                                   vector< FreqData * > **freqDataByChr,
                                   char TPED_MISSING, int nresample, bool PHASED, bool AUTO_FREQ)
{
    GarlicRNG *r = getRNG();

    igzstream fin;
    fin.open(tpedfile.c_str());

    if (fin.fail())
    {
        LOG.err("ERROR: Failed to open", tpedfile);
        throw 0;
    }

    string line;
    stringstream ss;
    char oneAllele = TPED_MISSING;
    int currChrLoci = 0;
    int ncols = 0;
    int nalleles = 0;
    int total = 0;
    string chr, locusName;
    string emptyChr = "_nochr";
    string prevChr = emptyChr;
    double gpos, ppos;

    vector<double> geneticPos;
    vector<pos_t> physicalPos;
    vector<string> locusNames;
    vector<char> allele;
    
    char alleleStr1, alleleStr2;

    vector< double > freq;
    vector< geno_t * > hap;
    vector< bool * > fc;
    geno_t *data;
    bool *firstCopy;

    numLoci = 0;

    bool firstLine = true;
    while (getline(fin, line)) {
        numLoci++;
        const char *p    = line.c_str();
        const char *pEnd = p + line.size();
        const char *tEnd;

        //Column count is a property of the file, not of each line.  It used to
        //be recomputed (and silently overwritten) on every line, so a ragged
        //TPED produced rows of differing length and out-of-bounds reads later.
        if (firstLine) {
            ncols = countFields(line) - 4;
            if (ncols <= 0 || (ncols % 2) != 0) {
                LOG.err("ERROR: Expected 4 leading columns and an even number of allele columns in", tpedfile);
                throw 0;
            }
            numInd = ncols / 2;
            firstLine = false;
        }

        //--- chromosome ---
        p = skipSpace(p, pEnd);
        tEnd = tokenEnd(p, pEnd);
        if (tEnd == p) {
            LOG.err("ERROR: Empty line or missing chromosome field at line", numLoci, false);
            LOG.err(" of", tpedfile);
            throw 0;
        }
        chr.assign(p, tEnd - p);
        p = tEnd;

        if (prevChr.compare(emptyChr) == 0 && numLoci == 1) prevChr = chr;

        if (chr.compare(prevChr) != 0){

            LOG.log("Chromosome",checkChrName(prevChr),false);
            LOG.log(":",currChrLoci,false);
            LOG.log(" sites.");

            (*mapDataByChr)->push_back(initMapData(geneticPos, physicalPos,locusNames, allele, currChrLoci, checkChrName(prevChr)));
            geneticPos.clear();
            physicalPos.clear();
            allele.clear();
            locusNames.clear();

            (*hapDataByChr)->push_back(initHapData(hap, fc, currChrLoci, numInd, PHASED));
            hap.clear();
            fc.clear();

            if(AUTO_FREQ){
                (*freqDataByChr)->push_back(initFreqData(freq, currChrLoci));
                freq.clear();
            }

            prevChr = chr;
            currChrLoci = 0;
        }

        currChrLoci++;

        //--- locus name ---
        p = skipSpace(p, pEnd);
        tEnd = tokenEnd(p, pEnd);
        if (tEnd == p) {
            LOG.err("ERROR: Missing locus ID at line", numLoci, false);
            LOG.err(" of", tpedfile);
            throw 0;
        }
        locusName.assign(p, tEnd - p);
        locusNames.push_back(locusName);
        p = tEnd;

        //--- genetic and physical position ---
        {
            char *q;
            p = skipSpace(p, pEnd);
            gpos = strtod(p, &q);
            if (q == p) {
                LOG.err("ERROR: Could not parse genetic position at line", numLoci, false);
                LOG.err(" of", tpedfile);
                throw 0;
            }
            p = q;
            geneticPos.push_back(gpos);

            p = skipSpace(p, pEnd);
            ppos = strtod(p, &q);
            if (q == p) {
                LOG.err("ERROR: Could not parse physical position at line", numLoci, false);
                LOG.err(" of", tpedfile);
                throw 0;
            }
            p = q;
            //strtod is kept so any decimal form still parses; the cast
            //truncates exactly as the old double-to-int assignment did, but
            //into 64 bits.
            physicalPos.push_back(pos_t(ppos));
        }

        //--- genotypes ---
        //Read 2*numInd individual non-whitespace characters, which is exactly
        //what `ss >> char` did (TPED alleles are single characters).
        nalleles = 0;
        total = 0;
        //Deliberately still raw: these rows are an ownership HANDOFF.  The
        //reader cannot size the whole block because it does not know nloci
        //yet, so it allocates a row per locus and initHapData consumes and
        //frees each one as it appends it (see appendRow in garlic-matrix.h).
        //Making them vectors would mean either copying twice or moving them
        //into a vector<vector<>>, which gives up the contiguity the LOD stage
        //depends on.
        data = new geno_t[numInd];
        if(PHASED) firstCopy = new bool[numInd];
        oneAllele = TPED_MISSING;

        for(int i = 0; i < numInd; i++){
            data[i] = 0;
            p = skipSpace(p, pEnd);
            if (p >= pEnd) {
                delete [] data;
                if (PHASED) delete [] firstCopy;
                LOG.err("ERROR: Too few genotype fields at line", numLoci, false);
                LOG.err(" of", tpedfile);
                throw 0;
            }
            alleleStr1 = *p++;
            p = skipSpace(p, pEnd);
            if (p >= pEnd) {
                delete [] data;
                if (PHASED) delete [] firstCopy;
                LOG.err("ERROR: Too few genotype fields at line", numLoci, false);
                LOG.err(" of", tpedfile);
                throw 0;
            }
            alleleStr2 = *p++;

            if(oneAllele == TPED_MISSING && alleleStr1 != TPED_MISSING) oneAllele = alleleStr1;
            if(oneAllele == TPED_MISSING && alleleStr2 != TPED_MISSING) oneAllele = alleleStr2;

            //Counted explicitly rather than accumulated into data[i] and then
            //clamped.  The old form added -9 per missing allele and clamped any
            //negative to -9, which collapsed a HALF call onto a fully missing
            //one -- losing the fact that one real allele was observed, even
            //though the frequency counters below had already used it.  The
            //frequency arithmetic here is unchanged; what changes is that the
            //matrix now records which case this was.
            int observed = 0, counted = 0;
            if (alleleStr1 != TPED_MISSING) { observed++; if (alleleStr1 == oneAllele) counted++; }
            if (alleleStr2 != TPED_MISSING) { observed++; if (alleleStr2 == oneAllele) counted++; }
            nalleles += counted;
            total    += observed;

            if (observed == 2)      data[i] = geno_t(counted);
            else if (observed == 1) data[i] = (counted == 1 ? GENO_HALF_COUNTED : GENO_HALF_OTHER);
            else                    data[i] = GENO_MISSING;

            if(PHASED) firstCopy[i] = (alleleStr1 == oneAllele);
        }

        p = skipSpace(p, pEnd);
        if (p != pEnd) {
            delete [] data;
            if (PHASED) delete [] firstCopy;
            LOG.err("ERROR: Too many columns at line", numLoci, false);
            LOG.err(" of", tpedfile, false);
            LOG.err(". Expected genotypes for", numInd, false);
            LOG.err(" individuals.");
            throw 0;
        }

        allele.push_back(oneAllele);
        hap.push_back(data);
        data = NULL;
        if(PHASED){
            fc.push_back(firstCopy);
            firstCopy = NULL;
        }

        if(AUTO_FREQ){
            freq.push_back(alleleFrequency(nalleles, total, nresample, r));
        }
    }

    LOG.log("Chromosome",checkChrName(chr),false);
    LOG.log(":",currChrLoci,false);
    LOG.log(" sites.");

    (*mapDataByChr)->push_back(initMapData(geneticPos,physicalPos,locusNames, allele, currChrLoci, checkChrName(chr)));
    geneticPos.clear();
    physicalPos.clear();
    allele.clear();
    locusNames.clear();

    (*hapDataByChr)->push_back(initHapData(hap, fc, currChrLoci, numInd, PHASED));
    hap.clear();
    fc.clear();

    if(AUTO_FREQ){
        (*freqDataByChr)->push_back(initFreqData(freq, currChrLoci));
        freq.clear();
    }
    
    return;
}

GenoLikeData *initGLData(const vector< double * > &GL, int nloci, int nind){
    GenoLikeData *glData = new GenoLikeData;
    glData->nind = nind;
    glData->nloci = nloci;
    //reserveRows rather than resize: see garlic-matrix.h, the zero pass would
    //commit every page before the copy overwrote it.
    glData->data.reserveRows(nloci, nind);
    for(int i = 0; i < nloci; i++){
        glData->data.appendRow(GL[i]);
        delete [] GL[i];
    }
    glData->data.finishRows();

    return glData;
}

double alleleFrequency(double nalleles, double total, int nresample, GarlicRNG *r)
{
    double freq = (total == 0) ? 0 : (nalleles / total);
    if (nresample > 0 && total != 0)
    {
        int count = 0;
        for (int i = 0; i < nresample; i++) if (r->uniform() <= freq) count++;
        freq = double(count) / double(nresample);
    }
    return freq;
}

vector< FreqData * > *calcFreqDataForIndices(vector< HapData * > *hapDataByChr,
                                             const vector<int> &keepInd,
                                             int nresample)
{
    GarlicRNG *r = getRNG();
    vector< FreqData * > *freqDataByChr = new vector< FreqData * >;

    for (unsigned int chr = 0; chr < hapDataByChr->size(); chr++)
    {
        HapData *hap = hapDataByChr->at(chr);
        int nloci = hap->nloci;
        vector<double> freq(nloci);
        for (int locus = 0; locus < nloci; locus++)
        {
            double nalleles = 0, total = 0;
            for (unsigned int k = 0; k < keepInd.size(); k++)
            {
                //Half calls contribute their one observed allele, matching both
                //loaders.  See the genotype encoding in garlic-data.h.
                addAlleleCounts(hap->data[locus][keepInd[k]], nalleles, total);
            }
            freq[locus] = alleleFrequency(nalleles, total, nresample, r);
        }
        freqDataByChr->push_back(initFreqData(freq, nloci));
    }
    return freqDataByChr;
}

FreqData *initFreqData(const vector<double> &freq, int nloci){
    FreqData *freqData = new FreqData;
    freqData->freq.resize(nloci);
    freqData->nloci = nloci;

    for(int i = 0; i < nloci; i++){
        freqData->freq[i] = freq[i];
    }

    return freqData;
}

HapData *initHapData(const vector< geno_t * > &hap, const vector< bool * > &fc, int nloci, int nind, bool PHASED){
    HapData *hapData = new HapData;
    hapData->nind = nind;
    hapData->nloci = nloci;
    //The reader allocates a row per locus because it does not know the locus
    //count in advance.  Gather them into one block here and release the rows,
    //so the transient overlap is one chromosome rather than the whole genome.
    hapData->data.reserveRows(nloci, nind);
    if(PHASED) hapData->firstCopy.reserveRows(nloci, nind);
    vector<unsigned char> fcRow(PHASED ? nind : 0);
    for(int i = 0; i < nloci; i++){
        hapData->data.appendRow(hap[i]);
        delete [] hap[i];
        if(PHASED){
            for(int j = 0; j < nind; j++) fcRow[j] = fc[i][j] ? 1 : 0;
            hapData->firstCopy.appendRow(&fcRow[0]);
            delete [] fc[i];
        }
    }
    hapData->data.finishRows();
    if(PHASED) hapData->firstCopy.finishRows();

    return hapData;
}

MapData *initMapData(const vector<double> &geneticPos, const vector<pos_t> &physicalPos, const vector<string> &locusNames, const vector<char> &allele, int nloci, string chr){

    MapData *mapData = new MapData;
    mapData->physicalPos.resize(nloci);
    mapData->geneticPos.resize(nloci);
    mapData->locusName.resize(nloci);
    mapData->allele.resize(nloci);
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

//Distinct labels in order of first appearance -- the same rule
//enumeratePopulations uses -- plus, for each individual, which one it is.
//Returns an empty name list when there is nothing to split on, which is the
//signal to write the single-column format.
static vector<int> groupPopulations(const vector<string> &popOfInd,
                                    vector<string> &names)
{
    names.clear();
    vector<int> idx(popOfInd.size(), 0);
    for (unsigned int i = 0; i < popOfInd.size(); i++)
    {
        unsigned int k = 0;
        for (; k < names.size(); k++) if (names[k] == popOfInd[i]) break;
        if (k == names.size()) names.push_back(popOfInd[i]);
        idx[i] = int(k);
    }
    if (names.size() <= 1) names.clear();
    return idx;
}

void freqOnlyVCF(string vcffile, string outfile, int nresample, bool PASS_ONLY,
                 const string &popfile)
{
    GarlicRNG *r = getRNG();

    string freqoutfile = outfile + ".freq.gz";
    ogzstream fout;
    fout.open(freqoutfile.c_str());
    if (fout.fail())
    {
        LOG.err("ERROR: Failed to open", freqoutfile);
        throw 0;
    }
    //Filled from the #CHROM line and --pop, below: a VCF carries no
    //population labels of its own.
    vector<string> popNames;
    vector<int> popOf;
    int npop = 1;


    igzstream fin;
    fin.open(vcffile.c_str());
    if (fin.fail())
    {
        LOG.err("ERROR: Failed to open", vcffile);
        throw 0;
    }

    const int VCF_FIXED = 9;
    string line, chr, locusName, ref, alt, filt, fmt;
    long long lineno = 0, nwritten = 0;
    int numInd = 0;
    bool haveHeader = false;

    while (getline(fin, line))
    {
        lineno++;
        if (line.size() >= 2 && line[0] == '#' && line[1] == '#') continue;

        const char *p    = line.c_str();
        const char *pEnd = p + line.size();
        const char *tEnd;

        if (!line.empty() && line[0] == '#')
        {
            int ncols = countFields(line);
            if (ncols < VCF_FIXED + 1)
            {
                LOG.err("ERROR:", vcffile, false);
                LOG.err(" has no sample columns.");
                throw 0;
            }
            numInd = ncols - VCF_FIXED;
            haveHeader = true;
            LOG.log("Samples in the VCF:", numInd);

            if (!popfile.empty() && popfile.compare("none") != 0)
            {
                //The sample IDs, so --pop can be applied to them.  Reusing
                //applyPopFile rather than re-parsing here keeps one reader
                //for the file and one set of error messages.
                vector<string> sampleIDs;
                {
                    stringstream hs(line);
                    string tok;
                    for (int c = 0; c < VCF_FIXED && hs >> tok; c++) ;
                    while (hs >> tok) sampleIDs.push_back(tok);
                }
                IndData *tmp = initIndData(int(sampleIDs.size()));
                for (unsigned int c = 0; c < sampleIDs.size(); c++)
                {
                    tmp->indID[c] = sampleIDs[c];
                    tmp->pop[c]   = "unknown";
                    tmp->sex[c]   = 0;
                }
                applyPopFile(popfile, tmp);
                popOf = groupPopulations(tmp->pop, popNames);
                releaseIndData(tmp);
                if (!popNames.empty()) npop = int(popNames.size());
            }

            fout << "CHR\tSNP\tPOS\tALLELE";
            if (popNames.empty()) fout << "\tFREQ";
            else for (int p = 0; p < npop; p++) fout << "\t" << popNames[p];
            fout << "\n";
            continue;
        }
        if (!haveHeader)
        {
            LOG.err("ERROR: data before the #CHROM header at line", lineno, false);
            LOG.err(" of", vcffile);
            throw 0;
        }
        if (line.empty()) continue;

        p = skipSpace(p, pEnd); tEnd = tokenEnd(p, pEnd); chr.assign(p, tEnd - p); p = tEnd;
        pos_t ppos;
        {
            char *q;
            p = skipSpace(p, pEnd);
            double v = strtod(p, &q);
            if (q == p)
            {
                LOG.err("ERROR: could not parse POS at line", lineno, false);
                LOG.err(" of", vcffile);
                throw 0;
            }
            p = q;
            ppos = pos_t(v);
        }
        p = skipSpace(p, pEnd); tEnd = tokenEnd(p, pEnd); locusName.assign(p, tEnd - p); p = tEnd;
        p = skipSpace(p, pEnd); tEnd = tokenEnd(p, pEnd); ref.assign(p, tEnd - p);       p = tEnd;
        p = skipSpace(p, pEnd); tEnd = tokenEnd(p, pEnd); alt.assign(p, tEnd - p);       p = tEnd;
        p = skipSpace(p, pEnd); tEnd = tokenEnd(p, pEnd);                                p = tEnd;
        p = skipSpace(p, pEnd); tEnd = tokenEnd(p, pEnd); filt.assign(p, tEnd - p);      p = tEnd;
        p = skipSpace(p, pEnd); tEnd = tokenEnd(p, pEnd);                                p = tEnd;
        p = skipSpace(p, pEnd); tEnd = tokenEnd(p, pEnd); fmt.assign(p, tEnd - p);       p = tEnd;

        //The same site rules as loadVCFData, in the same order, so the two
        //readers keep the same set of sites.
        bool pass = (filt.compare("PASS") == 0 || filt.compare(".") == 0);
        if (PASS_ONLY && !pass) continue;
        if (alt.find(',') != string::npos) continue;
        if (!isSNV(ref) || !isSNV(alt))    continue;
        int gtIndex = gtIndexOf(fmt);
        if (gtIndex < 0) continue;

        vector<double> nallelesBy(npop, 0.0), totalBy(npop, 0.0);
        for (int i = 0; i < numInd; i++)
        {
            p = skipSpace(p, pEnd);
            tEnd = tokenEnd(p, pEnd);
            if (tEnd == p)
            {
                LOG.err("ERROR: line", lineno, false);
                LOG.err(" of", vcffile, false);
                LOG.err(" has genotypes for fewer than", numInd, false);
                LOG.err(" samples.");
                throw 0;
            }
            int dosage, ploidy; bool fcopy, isPhased;
            if (!parseGT(p, tEnd, gtIndex, 1, dosage, fcopy, ploidy, isPhased))
            {
                LOG.err("ERROR: could not parse a genotype at line", lineno, false);
                LOG.err(" of", vcffile);
                throw 0;
            }
            if (ploidy != 2)
            {
                LOG.err("ERROR: a genotype at line", lineno, false);
                LOG.err(" of", vcffile, false);
                LOG.err(" has ploidy", ploidy, false);
                LOG.err("; garlic calls ROH from diploid genotypes only.");
                throw 0;
            }
            p = tEnd;
            //Two allele columns per sample in a TPED; here one genotype per
            //sample, so the sample index is the column index directly.
            {
                //One genotype per sample here, so the loop index IS the
                //sample index -- unlike a TPED, which has two columns each.
                int pp = (popNames.empty() || i >= int(popOf.size())) ? 0 : popOf[i];
                addAlleleCounts(geno_t(dosage), nallelesBy[pp], totalBy[pp]);
            }
        }

        vector<double> freqBy(npop);
        for (int pp = 0; pp < npop; pp++)
            freqBy[pp] = alleleFrequency(nallelesBy[pp], totalBy[pp], nresample, r);

        //ALT in the ALLELE column.  readFreqData compares it against
        //mapData->allele and flips to 1-f on a mismatch, so this file is
        //interchangeable with one written from a TPED.
        //An absent ID becomes CHROM:POS, matching loadVCFData, so a freq file
        //and a run over the same VCF agree on locus names -- readFreqData
        //errors on a locus-name mismatch.
        string nm = locusName;
        if (nm.compare(".") == 0)
        {
            stringstream ns;
            ns << chr << ":" << ppos;
            nm = ns.str();
        }
        fout << checkChrName(chr) << "\t" << nm
             << "\t" << ppos << "\t" << alt[0];
        for (int pp = 0; pp < npop; pp++) fout << "\t" << freqBy[pp];
        fout << "\n";
        nwritten++;
    }

    if (!haveHeader)
    {
        LOG.err("ERROR: no #CHROM header found in", vcffile, false);
        LOG.err(". Is it a VCF?");
        throw 0;
    }
    LOG.log("Sites written to the frequency file:", nwritten);

    fin.close();
    fout.close();
}

void freqOnly(string filename, string outfile, int nresample, char TPED_MISSING,
              const vector<string> &popOfInd){
    
    GarlicRNG *r = getRNG();

    string freqoutfile = outfile + ".freq.gz";

    ogzstream fout;
    fout.open(freqoutfile.c_str());
    if (fout.fail()){
        LOG.err("ERROR: Failed to open", freqoutfile);
        throw 0;
    }

    vector<string> popNames;
    vector<int> popOf = groupPopulations(popOfInd, popNames);
    const int npop = popNames.empty() ? 1 : int(popNames.size());

    //One column when there is one population, headed FREQ, exactly as before.
    fout << "CHR\tSNP\tPOS\tALLELE";
    if (popNames.empty()) fout << "\tFREQ";
    else for (int p = 0; p < npop; p++) fout << "\t" << popNames[p];
    fout << "\n";

    igzstream fin;
    fin.open(filename.c_str());

    if (fin.fail()){
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
        //Counted per population, but against ONE reference allele for the
        //locus -- the first non-missing one seen, as before -- so the columns
        //describe the same allele and the ALLELE field still means something.
        vector<double> nallelesBy(npop, 0.0), totalBy(npop, 0.0);
        for(count = 0; count < ncols-4; count++){
            ss >> junk;
            if(junk[0] != TPED_MISSING){
                //Two allele columns per individual.
                int who = count / 2;
                int p = (popNames.empty() || who >= int(popOf.size())) ? 0 : popOf[who];
                totalBy[p]++;
                if(oneAllele == TPED_MISSING) oneAllele = junk.c_str()[0]; 
                if(junk[0] == oneAllele) nallelesBy[p]++;
            }
        }

        vector<double> freqBy(npop);
        for (int p = 0; p < npop; p++)
            freqBy[p] = alleleFrequency(nallelesBy[p], totalBy[p], nresample, r);
        ss.clear();
        
        //pos_t, not int: this was the last 32-bit truncation of a physical
        //position left after 7899f42, and it silently wrapped past 2.147 Gb.
        fout << checkChrName(chr) << "\t" << locusName << "\t" << pos_t(ppos) << "\t" << oneAllele;
        for (int p = 0; p < npop; p++) fout << "\t" << freqBy[p];
        fout << "\n";
    }

    fin.close();
    fout.close();
}


double calcDensity(int numLoci, vector< MapData * > *mapDataByChr, centromere *centro){
    double density = numLoci;
    double length = 0;
    for(unsigned int i = 0; i < mapDataByChr->size(); i++){
        string str = mapDataByChr->at(i)->chr; 
        string chrstr = checkChrName(str);
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
                               int numThreads,
                               int ldSubsample)
{

    GarlicRNG *r = getRNG();

    //to hold the indicies of the randomly selected individuals
    int nind = hapDataByChr->at(0)->nind;
    vector<int> randInd;
    if (ldSubsample >= nind || ldSubsample <= 0)
    {
        ldSubsample = nind;
        randInd.resize(nind);
        for (int i = 0; i < nind; i++) randInd[i] = i;
    }
    else
    {
        vector<int> indIndex(nind);
        for (int i = 0; i < nind; i++) indIndex[i] = i;
        randInd.resize(ldSubsample);
        r->choose(randInd.data(), ldSubsample, indIndex.data(), nind);
        nind = ldSubsample;
    }


    vector< LDData * > *ldDataByChr = new vector< LDData * >;
    for(unsigned int chr = 0; chr < hapDataByChr->size(); chr++){
        cerr << mapDataByChr->at(chr)->chr << "    ";
        if(!PHASED) ldDataByChr->push_back(calcHR2LD(hapDataByChr->at(chr), genoFreqDataByChr->at(chr), winsize, numThreads, randInd.data(), ldSubsample));
        else ldDataByChr->push_back(calcR2LD(hapDataByChr->at(chr), freqDataByChr->at(chr), winsize, numThreads, randInd.data(), ldSubsample));
    }
    return ldDataByChr;
}

//Two phases: (1) evaluate every distinct pairwise LD in the band exactly once,
//(2) assemble the per-window sums from the band with a sliding update.
//Previously both were fused, so each pair was re-evaluated once per window
//containing it -- O(nloci * winsize^2 * nind) instead of O(nloci * winsize * nind).
//bandFn/bandOrders were a void *(*)(void *) pair that this function never
//used -- they existed so a caller could pass the cast worker pointer through,
//and were silenced with (void) casts at the bottom.  With std::thread the
//worker is named directly, so both parameters are gone.
static void runLDPhases(LDData *LD, int nloci, int winsize, int numThreads,
                        double *band, Bar *bar)
{
    //phase 2 only; phase 1 is type-specific and done by the caller
    int nWindows = nloci - winsize + 1;
    if (nWindows <= 0) return;

    int nt = numThreads;
    vector<unsigned int> NUM_PER_THREAD = make_thread_partition(nt, nWindows);
    vector<std::thread> peer;
    peer.reserve(nt);
    vector<BAND_work_order_t *> orders(nt);
    unsigned int previous = 0;
    for (int i = 0; i < nt; i++)
    {
        BAND_work_order_t *o = new BAND_work_order_t;
        o->winsize = winsize;
        o->nloci = nloci;
        o->band = band;
        o->LD = LD;
        o->bar = bar;
        o->start = previous;
        previous += NUM_PER_THREAD[i];
        o->stop = previous;
        orders[i] = o;
        peer.emplace_back(parallelLDFromBand, o);
    }
    for (int i = 0; i < nt; i++){
        peer[i].join();
        delete orders[i];
    }
}

LDData *calcHR2LD(HapData *hapData, GenoFreqData *genoFreqData, int winsize, int numThreads, int *indIndex, int ldSubsample){

    LDData *LD = initLDData(hapData->nloci, winsize);
    int nloci = hapData->nloci;
    int nWindows = nloci - winsize + 1;
    if (nWindows < 0) nWindows = 0;

    vector<double> bandBuf((size_t)nloci * (size_t)winsize);
    double *band = bandBuf.data();

    Bar bar;
    barInit(bar, double(nloci) + double(nWindows), 100);

    //--- phase 1: banded pairwise hr2 ---
    int nt = numThreads;
    vector<unsigned int> NUM_PER_THREAD = make_thread_partition(nt, nloci);
    vector<std::thread> peer;
    peer.reserve(nt);
    vector< HR2_work_order_t * > orders;
    unsigned int previous = 0;
    for (int i = 0; i < nt; i++)
    {
        HR2_work_order_t *order = new HR2_work_order_t;
        order->hapData = hapData;
        order->genoFreqData = genoFreqData;
        order->LD = LD;
        order->winsize = winsize;
        order->bar = &bar;
        order->start = previous;
        previous += NUM_PER_THREAD[i];
        order->stop = previous;
        order->indIndex = indIndex;
        order->ldSubsample = ldSubsample;
        order->band = band;

        peer.emplace_back(parallelHR2, order);
        orders.push_back(order);
    }
    for (int i = 0; i < nt; i++){
        peer[i].join();
        delete orders[i];
    }
    orders.clear();

    //--- phase 2: window sums from the band ---
    runLDPhases(LD, nloci, winsize, numThreads, band, &bar);

    finalize(bar);

    return LD;
}

LDData *calcR2LD(HapData *hapData, FreqData *freqData, int winsize, int numThreads, int *indIndex, int ldSubsample){

    LDData *LD = initLDData(hapData->nloci, winsize);
    int nloci = hapData->nloci;
    int nWindows = nloci - winsize + 1;
    if (nWindows < 0) nWindows = 0;

    vector<double> bandBuf((size_t)nloci * (size_t)winsize);
    double *band = bandBuf.data();

    Bar bar;
    barInit(bar, double(nloci) + double(nWindows), 100);

    //--- phase 1: banded pairwise r2 ---
    int nt = numThreads;
    vector<unsigned int> NUM_PER_THREAD = make_thread_partition(nt, nloci);
    vector<std::thread> peer;
    peer.reserve(nt);
    vector< R2_work_order_t * > orders;
    unsigned int previous = 0;
    for (int i = 0; i < nt; i++)
    {
        R2_work_order_t *order = new R2_work_order_t;
        order->hapData = hapData;
        order->freqData = freqData;
        order->LD = LD;
        order->winsize = winsize;
        order->bar = &bar;
        order->start = previous;
        previous += NUM_PER_THREAD[i];
        order->stop = previous;
        order->indIndex = indIndex;
        order->ldSubsample = ldSubsample;
        order->band = band;

        peer.emplace_back(parallelR2, order);
        orders.push_back(order);
    }
    for (int i = 0; i < nt; i++){
        peer[i].join();
        delete orders[i];
    }
    orders.clear();

    //--- phase 2: window sums from the band ---
    runLDPhases(LD, nloci, winsize, numThreads, band, &bar);

    finalize(bar);

    return LD;
}

//Was launched through (void *(*)(void *))parallelHR2 -- a cast between
//incompatible function pointer types, since this returns void.  std::thread
//calls it with its real signature, so the cast and the void * are both gone.
void parallelHR2(HR2_work_order_t *p){
    //advanceBar takes a global mutex and writes to cerr on every call; at one
    //call per locus that is ~577k lock/unlock pairs genome-wide.  Batch them.
    int barPending = 0;
    HapData *hapData = p->hapData;
    GenoFreqData *genoFreqData = p->genoFreqData;
    int winsize = p->winsize;
    int start = p->start;
    int stop = p->stop;
    Bar *bar = p->bar;
    int *indIndex = p->indIndex;
    int ldSubsample = p->ldSubsample;
    double *band = p->band;
    int nloci = hapData->nloci;

    //Fill the band rows owned by this thread: one hr2 evaluation per distinct
    //pair, rather than one per (window, pair) as before.
    for (int i = start; i < stop; i++){
        if ((++barPending) == 256) { advanceBar(*bar, 256); barPending = 0; }
        band[(size_t)i * winsize] = 1.0;
        int dMax = (winsize < nloci - i) ? winsize : (nloci - i);
        for (int d = 1; d < dMax; d++){
            band[(size_t)i * winsize + d] = hr2(hapData, genoFreqData, i, i + d, indIndex, ldSubsample);
        }
        for (int d = dMax; d < winsize; d++) band[(size_t)i * winsize + d] = 0.0;
    }
    if (barPending) advanceBar(*bar, barPending);

}

void parallelR2(R2_work_order_t *p){
    //advanceBar takes a global mutex and writes to cerr on every call; at one
    //call per locus that is ~577k lock/unlock pairs genome-wide.  Batch them.
    int barPending = 0;
    HapData *hapData = p->hapData;
    FreqData *freqData = p->freqData;
    int winsize = p->winsize;
    int start = p->start;
    int stop = p->stop;
    Bar *bar = p->bar;
    int *indIndex = p->indIndex;
    int ldSubsample = p->ldSubsample;
    double *band = p->band;
    int nloci = hapData->nloci;

    for (int i = start; i < stop; i++){
        if ((++barPending) == 256) { advanceBar(*bar, 256); barPending = 0; }
        band[(size_t)i * winsize] = 1.0;
        int dMax = (winsize < nloci - i) ? winsize : (nloci - i);
        for (int d = 1; d < dMax; d++){
            band[(size_t)i * winsize + d] = r2(hapData, freqData, i, i + d, indIndex, ldSubsample);
        }
        for (int d = dMax; d < winsize; d++) band[(size_t)i * winsize + d] = 0.0;
    }
    if (barPending) advanceBar(*bar, barPending);

}

//rho(a,b) for |a-b| < winsize, read out of the band (which stores a <= b).
static inline double bandLD(const double *band, int winsize, int a, int b)
{
    return (a <= b) ? band[(size_t)a * winsize + (b - a)]
                    : band[(size_t)b * winsize + (a - b)];
}

//LD[p][k] = sum over the window [p, p+winsize-1] of rho(p+k, j).
//Computing each row from scratch costs O(winsize^2); instead note that
//    LD[p][k] = LD[p-1][k+1] - rho(p+k, p-1) + rho(p+k, p+winsize-1)
//because both refer to the same site s = p+k and the window just slid by one.
//Only the last entry of each row needs a fresh O(winsize) sum, so a row costs
//O(winsize) rather than O(winsize^2).
void ldRowsFromBand(double *band, LDData *LD, int nloci, int winsize, int start, int stop, Bar *bar)
{
    if (start >= stop) return;
    int w = winsize;

    //First row of this thread's chunk: computed directly so that chunks are
    //independent and the recurrence never crosses a thread boundary.
    {
        int p0 = start;
        for (int k = 0; k < w; k++){
            int site = p0 + k;
            double sum = 0;
            for (int j = p0; j < p0 + w; j++) sum += bandLD(band, w, site, j);
            LD->LD[p0][k] = sum;
        }
        advanceBar(*bar,1);
    }

    for (int p = start + 1; p < stop; p++){
        advanceBar(*bar,1);
        double *prev = LD->LD[p - 1];
        double *cur  = LD->LD[p];
        int leaving  = p - 1;
        int entering = p + w - 1;
        for (int k = 0; k < w - 1; k++){
            int site = p + k;
            cur[k] = prev[k + 1]
                     - bandLD(band, w, site, leaving)
                     + bandLD(band, w, site, entering);
        }
        //last entry: site == entering, no predecessor to slide from
        double sum = 0;
        for (int j = p; j < p + w; j++) sum += bandLD(band, w, entering, j);
        cur[w - 1] = sum;
    }
}

void parallelLDFromBand(BAND_work_order_t *p){
    ldRowsFromBand(p->band, p->LD, p->nloci, p->winsize, p->start, p->stop, p->bar);
}

//Superseded by the banded implementation above; retained for reference.
void ldHR2(LDData *LD, HapData *hapData, GenoFreqData *genoFreqData, int site, int start, int end, int *indIndex, int ldSubsample) {
    for (int i = start; i <= end; i++) {
        if (i != site) LD->LD[start][site-start] += hr2(hapData, genoFreqData, i, site, indIndex, ldSubsample);
        else LD->LD[start][site-start] += 1;
    }
    return;
}

void ldR2(LDData *LD, HapData *hapData, FreqData *freqData, int site, int start, int end, int *indIndex, int ldSubsample) {
    for (int i = start; i <= end; i++) {
        if (i != site) LD->LD[start][site-start] += r2(hapData, freqData, i, site, indIndex, ldSubsample);
        else LD->LD[start][site-start] += 1;
    }
    return;
}


vector<unsigned int> make_thread_partition(int &numThreads, int ncols) {
    if (numThreads > ncols) numThreads = ncols;
    vector<unsigned int> NUM_PER_THREAD(numThreads);
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


double hr2(HapData *hapData, GenoFreqData *genoFreqData, int i, int j, int *indIndex, int ldSubsample) {
    double HA = genoFreqData->homFreq[i];
    double HB = genoFreqData->homFreq[j];
    if(HA > 0 && HA < 1 && HB > 0 && HB < 1){
        double HAB = 0;
        double total = 0;
        //for (int ind = 0; ind < hapData->nind; ind++) {
        for (int k = 0; k < ldSubsample; k++) {
            int ind = indIndex[k];
            if (genoIsCalled(hapData->data[i][ind]) && genoIsCalled(hapData->data[j][ind])) {
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

double r2(HapData *hapData, FreqData *freqData, int i, int j, int *indIndex, int ldSubsample) {
    
    double pi = freqData->freq[i];
    double pj = freqData->freq[j];

    if(pi > 0 && pi < 1 && pj > 0 && pj < 1){
        double x11 = 0;
        double total = 0;
        //for (int ind = 0; ind < hapData->nind; ind++) {
        for (int k = 0; k < ldSubsample; k++) {
            int ind = indIndex[k];
            if (genoIsCalled(hapData->data[i][ind]) && genoIsCalled(hapData->data[j][ind])) {
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
    data->LD.assign(nloci, winsize, 0);
    return data;
}
void releaseLDData(LDData *data){
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
            if (genoIsCalled(hapData->data[locus][ind]))
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
    data->homFreq.resize(nloci);
    data->nloci = nloci;
    return data;
}

void releaseGenoFreq(GenoFreqData *genoFreqData){
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

bool alignMapScaffold(vector< GenMapScaffold * > *scaffoldMapByChr, vector< MapData * > *mapDataByChr) {
    map<string, GenMapScaffold *> byName;
    for (unsigned int i = 0; i < scaffoldMapByChr->size(); i++) {
        GenMapScaffold *sc = scaffoldMapByChr->at(i);
        if (byName.count(sc->chr) > 0) {
            LOG.err("ERROR: Genetic map contains more than one block for", sc->chr);
            LOG.err("\tAll records for a chromosome must be contiguous in the map file.");
            return false;
        }
        byName[sc->chr] = sc;
    }

    vector< GenMapScaffold * > ordered;
    for (unsigned int i = 0; i < mapDataByChr->size(); i++) {
        string want = mapDataByChr->at(i)->chr;
        if (byName.count(want) == 0) {
            LOG.err("ERROR: Genetic map has no data for", want);
            return false;
        }
        ordered.push_back(byName[want]);
        byName.erase(want);
    }
    if (!byName.empty()) {
        LOG.err("ERROR: Genetic map contains chromosomes absent from the data, e.g.", byName.begin()->first);
        return false;
    }

    for (unsigned int i = 0; i < ordered.size(); i++) scaffoldMapByChr->at(i) = ordered[i];
    return true;
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

double getMapInfo(pos_t queryPos, GenMapScaffold *scaffold, int &count) {
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
        //currentIndex is a forward-only cursor, which is correct only while
        //queries arrive in ascending order.  If a query sits behind the cursor
        //(TPED positions not sorted within a chromosome) the scan used to fall
        //through and interpolate on uninitialised indices -- -Wall flags this
        //as "used uninitialized whenever 'for' loop exits because its condition
        //is false".  Fall back to a rescan from the start, and fail loudly if
        //the position really is not bracketed.
        int startIndex = -1;
        int endIndex = -1;
        for (/*scaffold->currentIndex*/; scaffold->currentIndex < scaffold->nloci - 1; scaffold->currentIndex++)
        {
            if (queryPos > scaffold->physicalPos[scaffold->currentIndex] && queryPos < scaffold->physicalPos[scaffold->currentIndex + 1])
            {
                startIndex = scaffold->currentIndex;
                endIndex = scaffold->currentIndex + 1;
                break;
            }
        }
        if (startIndex < 0)
        {
            for (scaffold->currentIndex = 0; scaffold->currentIndex < scaffold->nloci - 1; scaffold->currentIndex++)
            {
                if (queryPos > scaffold->physicalPos[scaffold->currentIndex] && queryPos < scaffold->physicalPos[scaffold->currentIndex + 1])
                {
                    startIndex = scaffold->currentIndex;
                    endIndex = scaffold->currentIndex + 1;
                    break;
                }
            }
        }
        if (startIndex < 0)
        {
            LOG.err("ERROR: Physical position", queryPos, false);
            LOG.err(" on", scaffold->chr, false);
            LOG.err(" is not bracketed by the genetic map scaffold.");
            throw 0;
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
    if (!LOG.isQuiet()) cerr << "Opening " << mapfile << "...\n";
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

    if (!LOG.isQuiet()) cerr << "Loading genetic map scaffold for " << nloci << " loci.\n";

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
    scaffoldMap->physicalPos.resize(nloci);
    scaffoldMap->geneticPos.resize(nloci);
    scaffoldMap->nloci = nloci;
    scaffoldMap->currentIndex = 0;
    return scaffoldMap;
}

void releaseGenMapScaffold(GenMapScaffold *scaffoldMap) {
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

//The site-retention predicate, in ONE place.  A NULL scaffold means "drop
//monomorphic sites only"; a non-NULL scaffold additionally drops sites outside
//the genetic map's span and inside the assembled centromere gap.
//
//keep[k] is the index in the ORIGINAL arrays of the k-th retained site, so
//keep.size() is the new locus count AND every gather below is driven by the
//same list.  That coupling is the point.  Previously this predicate was
//re-derived inline in all ten functions of this family -- 16 copies of the
//monomorphic test, 8 of the out-of-bounds clause -- and the arrangement was
//worse than plain duplication: the driver called the MapData overload first
//with newLoci = 0 so THAT copy sized every allocation, then handed the count
//to the HapData/GenoLikeData/FreqData overloads, which skipped their own
//counting pass and filled the arrays using their own copies of the predicate.
//No gather loop bounded its write index against the count.  So any divergence
//between the copy that counted and a copy that gathered was either a write
//past the end of the destination or, worse because it is silent, arrays of
//different lengths with loci misaligned between MapData and HapData.
//
//B3 was exactly this class of bug: locusName[index] was assigned
//physicalPos[i], and because the code was duplicated the same one-line error
//had to be found and fixed in two separate copies.
//
//Declared in the header rather than static so it can be unit-tested directly
//against a hand-built scaffold, which the inline version could not be.
vector<int> keepSites(const MapData *mapData, const FreqData *freqData,
                      const GenMapScaffold *scaffold)
{
    vector<int> keep;
    keep.reserve(freqData->nloci);

    //Bounded by freqData->nloci, as every one of the ten predecessors was.
    for (int i = 0; i < freqData->nloci; i++)
    {
        if (!(freqData->freq[i] > 0 && freqData->freq[i] < 1)) continue;

        if (scaffold != NULL)
        {
            if (mapData->physicalPos[i] < scaffold->physicalPos[0]) continue;
            if (mapData->physicalPos[i] > scaffold->physicalPos[scaffold->nloci - 1]) continue;
            if (mapData->physicalPos[i] > scaffold->centroStart &&
                    mapData->physicalPos[i] < scaffold->centroEnd) continue;
        }

        keep.push_back(i);
    }

    return keep;
}

//Gather the retained elements of one array.  dst must already have
//keep.size() elements, which is how every caller allocates it.
template <typename T>
static void gather(vector<T> &dst, const vector<T> &src, const vector<int> &keep)
{
    for (size_t k = 0; k < keep.size(); k++) dst[k] = src[keep[k]];
}

//Same for the matrix payloads.  Rows are loci after the contiguous layout
//change, so each retained locus is a single contiguous row copy rather than
//an element loop over individuals.
template <typename T>
static void gatherRows(Matrix<T> &dst, const Matrix<T> &src,
                       const vector<int> &keep, int ncol)
{
    for (size_t k = 0; k < keep.size(); k++)
        memcpy(dst[k], src[keep[k]], sizeof(T) * size_t(ncol));
}

//One implementation for both public entry points.  scaffoldMapByChr == NULL
//selects the monomorphic-only filter; otherwise the out-of-bounds and
//centromere clauses apply too.
static int filterSites(vector< MapData * > **mapDataByChr,
                       vector< HapData * > **hapDataByChr,
                       vector< FreqData * > **freqDataByChr,
                       vector< GenoLikeData * > **GLDataByChr,
                       vector< GenMapScaffold * > *scaffoldMapByChr,
                       bool USE_GL, bool PHASED)
{
    vector< MapData * > *mapDataByChr2 = new vector< MapData * >;
    vector< HapData * > *hapDataByChr2 = new vector< HapData * >;
    vector< FreqData * > *freqDataByChr2 = new vector< FreqData * >;
    vector< GenoLikeData * > *GLDataByChr2 = NULL;
    if (USE_GL) GLDataByChr2 = new vector< GenoLikeData * >;

    int numLoci = 0;
    for (unsigned int i = 0; i < (*mapDataByChr)->size(); i++)
    {
        MapData  *mapData  = (*mapDataByChr)->at(i);
        HapData  *hapData  = (*hapDataByChr)->at(i);
        FreqData *freqData = (*freqDataByChr)->at(i);
        const GenMapScaffold *scaffold =
            (scaffoldMapByChr != NULL) ? scaffoldMapByChr->at(i) : NULL;

        vector<int> keep = keepSites(mapData, freqData, scaffold);
        int newLoci = int(keep.size());

        MapData *mapData2 = initMapData(newLoci);
        mapData2->chr = mapData->chr;
        gather(mapData2->physicalPos, mapData->physicalPos, keep);
        gather(mapData2->geneticPos,  mapData->geneticPos,  keep);
        gather(mapData2->locusName,   mapData->locusName,   keep);
        gather(mapData2->allele,      mapData->allele,      keep);

        HapData *hapData2 = initHapData(hapData->nind, newLoci, PHASED);
        gatherRows(hapData2->data, hapData->data, keep, hapData->nind);
        if (PHASED)
            gatherRows(hapData2->firstCopy, hapData->firstCopy, keep, hapData->nind);

        GenoLikeData *GLData2 = NULL;
        if (USE_GL)
        {
            GenoLikeData *GLData = (*GLDataByChr)->at(i);
            GLData2 = initGLData(GLData->nind, newLoci);
            gatherRows(GLData2->data, GLData->data, keep, GLData->nind);
        }

        FreqData *freqData2 = initFreqData(newLoci);
        gather(freqData2->freq, freqData->freq, keep);

        mapDataByChr2->push_back(mapData2);
        hapDataByChr2->push_back(hapData2);
        if (USE_GL) GLDataByChr2->push_back(GLData2);
        freqDataByChr2->push_back(freqData2);
        numLoci += newLoci;
    }

    releaseMapData(*mapDataByChr);
    releaseHapData(*hapDataByChr);
    releaseFreqData(*freqDataByChr);
    if (USE_GL) releaseGLData(*GLDataByChr);

    *mapDataByChr  = mapDataByChr2;
    *hapDataByChr  = hapDataByChr2;
    if (USE_GL) *GLDataByChr = GLDataByChr2;
    *freqDataByChr = freqDataByChr2;

    return numLoci;
}

int filterMonomorphicSites(vector< MapData * > **mapDataByChr,
                           vector< HapData * > **hapDataByChr,
                           vector< FreqData * > **freqDataByChr,
                           vector< GenoLikeData * > **GLDataByChr,
                           bool USE_GL, bool PHASED)
{
    return filterSites(mapDataByChr, hapDataByChr, freqDataByChr, GLDataByChr,
                       NULL, USE_GL, PHASED);
}

int filterMonomorphicAndOOBSites(vector< MapData * > **mapDataByChr,
                                 vector< HapData * > **hapDataByChr,
                                 vector< FreqData * > **freqDataByChr,
                                 vector< GenoLikeData * > **GLDataByChr,
                                 vector< GenMapScaffold * > *scaffoldMapByChr,
                                 bool USE_GL, bool PHASED)
{
    return filterSites(mapDataByChr, hapDataByChr, freqDataByChr, GLDataByChr,
                       scaffoldMapByChr, USE_GL, PHASED);
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

//calcFreqData/calcFreqData2 removed: commented out since before this branch,
//their only call site in main was commented out too, and they were the last
//mention of gsl_rng in the tree.

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
    data->freq.resize(nloci);
    //data->allele.resize(nloci);

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

void writeFreqDataWide(string freqOutfile,
                       const vector< vector< FreqData * >* > &freqByPop,
                       const vector<string> &popNames,
                       vector< MapData * > *mapDataByChr)
{
    freqOutfile += ".gz";
    ogzstream fout;
    fout.open(freqOutfile.c_str());

    if (fout.fail())
    {
        LOG.err("ERROR: Failed to open", freqOutfile);
        throw 0;
    }

    fout << "CHR\tSNP\tPOS\tALLELE";
    for (unsigned int p = 0; p < popNames.size(); p++) fout << "\t" << popNames[p];
    fout << "\n";

    for (unsigned int chr = 0; chr < mapDataByChr->size(); chr++)
    {
        for (int locus = 0; locus < mapDataByChr->at(chr)->nloci; locus++)
        {
            fout << mapDataByChr->at(chr)->chr << "\t"
                 << mapDataByChr->at(chr)->locusName[locus] << "\t"
                 << pos_t(mapDataByChr->at(chr)->physicalPos[locus]) << "\t"
                 << mapDataByChr->at(chr)->allele[locus];
            for (unsigned int p = 0; p < freqByPop.size(); p++)
                fout << "\t" << freqByPop[p]->at(chr)->freq[locus];
            fout << "\n";
        }
    }

    fout.close();
    LOG.log("Wrote allele frequencies for", int(popNames.size()), false);
    LOG.log(" populations to", freqOutfile);
    return;
}

void writeFreqData(string freqOutfile,
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
                 //pos_t, not int: the last 32-bit truncation of a physical
                 //position left after 7899f42.  A position past 2.147 Gb wrote
                 //a negative number into the frequency file, and readFreqData
                 //would then fail the locus-name check on the way back in.
                 << pos_t(mapDataByChr->at(chr)->physicalPos[locus]) << "\t"
                 << mapDataByChr->at(chr)->allele[locus] << "\t"
                 << freqDataByChr->at(chr)->freq[locus] << "\n";
        }
    }
    cout << "Wrote allele frequency data to " << freqOutfile << endl;
    fout.close();
    return;
}

//Split a whitespace-separated header line into fields.
static vector<string> splitFields(const string &line)
{
    vector<string> out;
    stringstream ss(line);
    string f;
    while (ss >> f) out.push_back(f);
    return out;
}

vector< vector< FreqData * >* > *readFreqData(string freqfile,
                                              vector< MapData * > *mapDataByChr,
                                              const vector<string> &populations)
{
    //Pooled analysis asks for nothing by name and gets one set.
    bool pooled = populations.empty();
    unsigned int nwanted = (pooled ? 1 : (unsigned int)populations.size());

    int expectedRows = 1;
    for (unsigned int chr = 0; chr < mapDataByChr->size(); chr++)
        expectedRows += mapDataByChr->at(chr)->nloci;

    igzstream fin;
    fin.open(freqfile.c_str());
    if (fin.fail()){
        LOG.err("ERROR: Failed to open", freqfile);
        throw 0;
    }

    if (!LOG.isQuiet()) cerr << "Reading " << freqfile << "\n";
    string line;
    int currentCols = 0;
    int previousCols = -1;

    string header;
    getline(fin, header);
    vector<string> headerFields = splitFields(header);
    if (headerFields.size() < 5)
    {
        LOG.err("ERROR: header of", freqfile, false);
        LOG.err(" has fewer than 5 fields; expected CHR SNP POS ALLELE <population>...");
        throw 0;
    }
    //Everything after ALLELE names a population.
    vector<string> colName(headerFields.begin() + 4, headerFields.end());
    int minCols = 4 + int(colName.size());

    //Which column feeds which requested population.
    vector<int> columnFor(nwanted, 0);
    if (colName.size() == 1)
    {
        //Every frequency file garlic has ever written, and any file written
        //for a single population.  Applies to everyone, whatever it is named
        //-- the name is not checked, so a file headed FREQ, MAF or anything
        //else keeps working.
        for (unsigned int k = 0; k < nwanted; k++) columnFor[k] = 0;
    }
    else if (pooled)
    {
        LOG.err("ERROR:", freqfile, false);
        LOG.err(" has frequencies for", int(colName.size()), false);
        LOG.err(" populations, but this run analyses every sample as one.");
        LOG.err("\tThere is no correct way to choose between them.");
        LOG.err("\tSupply a file with a single frequency column, which applies to everyone.");
        throw 0;
    }
    else
    {
        for (unsigned int k = 0; k < nwanted; k++)
        {
            int found = -1;
            for (unsigned int c = 0; c < colName.size(); c++)
                if (colName[c] == populations[k]) { found = int(c); break; }
            if (found < 0)
            {
                LOG.err("ERROR:", freqfile, false);
                LOG.err(" has no frequency column for population", populations[k]);
                LOG.err("\tColumns present:", header);
                throw 0;
            }
            columnFor[k] = found;
        }
        //Extra columns are fine: one file may carry every population the
        //project has, and a run may analyse a subset of them.
    }

    //Allocated only now: every check above throws, and an allocation made
    //before them would be leaked on the way out.
    vector< vector< FreqData * >* > *out = new vector< vector< FreqData * >* >;
    for (unsigned int k = 0; k < nwanted; k++)
    {
        vector< FreqData * > *byChr = new vector< FreqData * >;
        for (unsigned int chr = 0; chr < mapDataByChr->size(); chr++)
            byChr->push_back(initFreqData(mapDataByChr->at(chr)->nloci));
        out->push_back(byChr);
    }

    stringstream ss;
    string locusID, chromosome;
    char allele;
    double position;
    int lineNum = 1;
    vector<double> value(colName.size());
    try
    {
    for (unsigned int chr = 0; chr < mapDataByChr->size(); chr++)
    {
        for (int locus = 0; locus < mapDataByChr->at(chr)->nloci; locus++)
        {
            lineNum++;
            getline(fin, line);
            if(fin.fail()){
                LOG.err("ERROR: at line", lineNum, false);
                LOG.err(" in", freqfile, false);
                LOG.err(". Perhaps too few lines?");
                throw 0;
            }
            ss.str(line);

            currentCols = countFields(line);
            if (currentCols < minCols){
                LOG.err("ERROR: Found", currentCols, false);
                LOG.err(" in", freqfile, false);
                LOG.err(" on line", lineNum, false);
                LOG.err(" but expected at least", minCols);
                throw 0;
            }
            if (currentCols != previousCols && previousCols != -1){
                LOG.err("ERROR: Differing number of columns across rows found in", freqfile);
                throw 0;
            }
            previousCols = currentCols;

            ss >> chromosome >> locusID >> position >> allele;
            for (unsigned int c = 0; c < colName.size(); c++) ss >> value[c];
            if (ss.fail())
            {
                LOG.err("ERROR: could not read", int(colName.size()), false);
                LOG.err(" frequencies on line", lineNum, false);
                LOG.err(" of", freqfile);
                throw 0;
            }
            if (mapDataByChr->at(chr)->locusName[locus].compare(locusID) != 0){
                LOG.err("ERROR: Loci appear mismatched in:", freqfile);
                LOG.err("ERROR: at line:", lineNum);
                LOG.err("ERROR: freq file locus name:", locusID);
                LOG.err("ERROR: tped file locus name:", mapDataByChr->at(chr)->locusName[locus]);
                throw 0;
            }
            //Does the internal coding of the '1' allele for this dataset match
            //the one in the freq file?  If not, take 1-freq.  ONE decision per
            //locus, applied to every population -- which is why the format is
            //wide: the counted allele belongs to the site, not the population.
            bool flip = (mapDataByChr->at(chr)->allele[locus] != allele);
            for (unsigned int k = 0; k < nwanted; k++)
            {
                double f = value[columnFor[k]];
                out->at(k)->at(chr)->freq[locus] = (flip ? 1 - f : f);
            }
            ss.clear();
        }
    }

    if (lineNum != expectedRows)
    {
        LOG.err("ERROR:", freqfile, false);
        LOG.err(" has", lineNum - 1, false);
        LOG.err(" rows but expected", expectedRows - 1);
        throw 0;
    }
    }
    catch (...)
    {
        //The rows are read into `out`; a malformed row must not leak it.
        for (unsigned int k = 0; k < out->size(); k++) releaseFreqData(out->at(k));
        out->clear();
        delete out;
        throw;
    }

    fin.close();

    return out;
}

//allocates the arrays and populates them with MISSING or "--" depending on type
vector< MapData * > *cloneMapData(const vector< MapData * > *mapDataByChr)
{
    vector< MapData * > *copy = new vector< MapData * >;
    for (unsigned int chr = 0; chr < mapDataByChr->size(); chr++)
    {
        const MapData *src = mapDataByChr->at(chr);
        MapData *dst = new MapData;
        dst->physicalPos = src->physicalPos;
        dst->geneticPos  = src->geneticPos;
        dst->locusName   = src->locusName;
        dst->allele      = src->allele;
        dst->nloci       = src->nloci;
        dst->chr         = src->chr;
        copy->push_back(dst);
    }
    return copy;
}

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
    data->locusName.resize(nloci);
    data->physicalPos.resize(nloci);
    data->geneticPos.resize(nloci);
    data->allele.resize(nloci);
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
    //The members own themselves now; there is nothing here to forget.
    delete data;
    return;
}

//Names a value's location for a diagnostic: "chr21:14834892 (sample HGDP00521)".
static string tglsWhere(vector< MapData * > *mapDataByChr, IndData *indData,
                        unsigned int chr, int locus, int ind)
{
    stringstream ss;
    ss << mapDataByChr->at(chr)->chr << ":" << mapDataByChr->at(chr)->physicalPos[locus];
    if (indData != NULL && ind < indData->nind) ss << " (sample " << indData->indID[ind] << ")";
    else ss << " (column " << (ind + 5) << ")";
    return ss.str();
}

vector< GenoLikeData * > *readTGLSData(string filename,
                                       int expectedLoci,
                                       int expectedInd,
                                       vector< MapData * > *mapDataByChr,
                                       string GL_TYPE,
                                       IndData *indData)
{
    igzstream fin;
    if (!LOG.isQuiet()) cerr << "Loading genotype likelihoods from " << filename << "\n";
    fin.open(filename.c_str());

    if (fin.fail()){
         LOG.err("ERROR: Failed to open", filename);
        throw 0;
    }
    
    string junk;
    double gl;
    string line;
    stringstream ss;
    vector< GenoLikeData * > *GLDataByChr = new vector< GenoLikeData * >;

    //A tgls file states P(genotype CORRECT) for PL and GL (see HELP_GL_TYPE),
    //which is NOT the quantity a VCF's PL/GL fields carry: VCF normalises them
    //so the called genotype is exactly 0.  Feeding those in yields an error
    //rate of 0, clamped to 1e-16 by glToError, which makes every heterozygote
    //contribute -16 to the window LOD instead of the -3 an error of 0.001
    //gives -- 113 ROH instead of 171 on the bundled chr21 data, with no other
    //symptom.  The two checks below reject that input instead of computing
    //with it.  Neither applies to GQ, where 0 legitimately means "no
    //confidence" and yields an error rate of 1.
    bool rejectZero = (GL_TYPE.compare("PL") == 0 || GL_TYPE.compare("GL") == 0);
    bool allInteger = true;     //verdict deferred: a property of the whole file
    long nvalues = 0;

    //For each chromosome
    for (unsigned int chr = 0; chr < mapDataByChr->size(); chr++){
        GenoLikeData *data = initGLData(expectedInd, mapDataByChr->at(chr)->nloci);
        GLDataByChr->push_back(data);
        //For each locus on the chromosome
        for (int locus = 0; locus < mapDataByChr->at(chr)->nloci; locus++){
            getline(fin,line);
            int num = countFields(line); 
            if(num != expectedInd+4){
                LOG.err("ERROR: Incorrect number of columns in tgls file: ", num, false);
                LOG.err(". Expected: ", expectedInd);
                throw 0;
            }
            ss.str(line);
            ss >> junk;
            ss >> junk;
            ss >> junk;
            ss >> junk;
            for (int ind = 0; ind < expectedInd; ind++){
                //An unparseable token used to leave gl at 0 and put the stream
                //in a failed state, so the rest of the line silently became
                //zeros -- error rate 1 under GQ, 1e-16 under PL/GL.
                if (!(ss >> gl)){
                    LOG.err("ERROR: could not read a " + GL_TYPE + " value at",
                            tglsWhere(mapDataByChr, indData, chr, locus, ind), false);
                    LOG.err(" in", filename);
                    throw 0;
                }
                nvalues++;

                //Values that cannot occur on their own scale.
                if (GL_TYPE.compare("GL") == 0 && gl > 0){
                    LOG.err("ERROR: GL is log10 of a probability and cannot exceed 0, but the value at",
                            tglsWhere(mapDataByChr, indData, chr, locus, ind), false);
                    LOG.err(" is", gl);
                    throw 0;
                }
                if (GL_TYPE.compare("GL") != 0 && gl < 0){
                    LOG.err("ERROR: " + GL_TYPE + " is phred-scaled and cannot be negative, but the value at",
                            tglsWhere(mapDataByChr, indData, chr, locus, ind), false);
                    LOG.err(" is", gl);
                    throw 0;
                }

                if (rejectZero && gl == 0){
                    LOG.err("ERROR: --gl-type " + GL_TYPE + ", but the value at",
                            tglsWhere(mapDataByChr, indData, chr, locus, ind), false);
                    LOG.err(" in", filename, false);
                    LOG.err(" is 0.");
                    LOG.err("\ta " + GL_TYPE + " of 0 asserts P(genotype correct) = 1 exactly, which becomes a");
                    LOG.err("\tper-genotype error of 0 -- clamped to 1e-16, making every heterozygote");
                    LOG.err("\tcontribute -16 to the window LOD instead of the -3 that an error rate of");
                    LOG.err("\t0.001 gives.  On the bundled chr21 data that is 113 ROH instead of 171,");
                    LOG.err("\twith no other symptom.");
                    LOG.err("\tIf these values came from a VCF: VCF normalises PL and GL so the CALLED");
                    LOG.err("\tgenotype is exactly 0, which is not the quantity --tgls expects.  Either");
                    LOG.err("\tuse --gl-type GQ, a phred-scaled P(call is wrong) that needs no");
                    LOG.err("\tnormalisation, or supply values that really are P(genotype correct).");
                    throw 0;
                }

                if (gl != floor(gl)) allInteger = false;

                GLDataByChr->at(chr)->data[locus][ind] = glToError(gl, GL_TYPE);

            }
            ss.clear();
        }
    }

    //Deferred verdict: no INTEGER PL expresses a realistic genotype error rate
    //under this convention -- the integers map to 0, 0.206, 0.369, 0.499, ...
    //and an error of 0.001 would need PL = 0.00435.  An all-integer PL file is
    //therefore a VCF-derived file however it was produced.
    if (GL_TYPE.compare("PL") == 0 && allInteger && nvalues > 0){
        LOG.err("ERROR: --gl-type PL, but every value in", filename, false);
        LOG.err(" is an integer.");
        LOG.err("\t--tgls expects PL as a phred-scaled P(genotype CORRECT), which requires");
        LOG.err("\tfractional values -- an error rate of 0.001 is PL = 0.00435.  Integers map");
        LOG.err("\tto error 0 (PL=0), 0.206 (1), 0.369 (2), 0.499 (3), so no integer PL");
        LOG.err("\texpresses a realistic genotype error rate.");
        LOG.err("\tInteger PLs almost certainly came from a VCF.  Use --gl-type GQ instead,");
        LOG.err("\twhich is a phred-scaled P(call is wrong) and needs no normalisation.");
        throw 0;
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

    //MISSING, not GENO_UNSET: this array holds LOD scores, and the sentinel IS
    //load-bearing here -- tail windows retain it, writeWinData prints NA for it,
    //and the KDE input skips it with x != MISSING.  -9 would be a plausible LOD
    //score and would be taken as real data.
    data->data.assign(nind, nloci, MISSING);

    return data;
}

void releaseWinData(WinData *data)
{
    if (data == NULL) return;
    delete data;
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
    //string popName = indData->pop;

    for (int chr = 0; chr < numChr; chr++)
    {
        string rawWinOutfile = outfile;
        //rawWinOutfile += ".";
        //rawWinOutfile += popName;
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
        if (!LOG.isQuiet()) cerr << "Wrote " << rawWinOutfile << "\n";
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

    data->data.assign(nloci, nind, GENO_UNSET);
    if(PHASED) data->firstCopy.assign(nloci, nind, 0);

    return data;
}

void releaseHapData(HapData *data)
{
    if (data == NULL) return;

    //The Matrix members free their own contiguous block; there is no row-0
    //bookkeeping left to get wrong.
    delete data;
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

    data->data.assign(nloci, nind, 1);

    return data;
}

void releaseGLData(GenoLikeData *data) {
    if (data == NULL) return;
    delete data;
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

string lc(string str) {
    char c[2] = {' ', '\0'};
    for (unsigned int i = 0; i < str.size(); i++) {
        c[0] = char(tolower(str[i]));
        str.replace(i, 1, c);
    }
    return str;
}

int formatIndexOf(const string &fmt, const string &key)
{
    int k = 0;
    size_t a = 0;
    while (a <= fmt.size())
    {
        size_t b = fmt.find(':', a);
        size_t len = (b == string::npos ? fmt.size() : b) - a;
        if (len == key.size() && fmt.compare(a, len, key) == 0) return k;
        if (b == string::npos) break;
        a = b + 1; k++;
    }
    return -1;
}

int gtIndexOf(const string &fmt) { return formatIndexOf(fmt, "GT"); }

bool isSNV(const string &s)
{
    if (s.size() != 1) return false;
    char c = s[0];
    if (c >= 'a' && c <= 'z') c = char(c - 'a' + 'A');
    return (c == 'A' || c == 'C' || c == 'G' || c == 'T');
}

//A VCF is self-describing where a TPED is not: it names its samples, states
//REF and ALT, states ploidy per call, and marks phase per call.  Everything
//loadTPEDData has to infer, this can check.  The outputs are deliberately the
//same structures loadTPEDData produces, so nothing downstream knows which
//reader ran.
//
//Two differences from the TPED path that are decisions, not omissions:
//
//  - the counted allele is ALT, not "the first allele observed at the locus".
//    allele[locus] is set to the ALT character and freq to the ALT frequency
//    over non-missing calls.  The LOD is symmetric in the counted allele, so
//    this changes no call; it makes the .freq.gz interchangeable with the TPED
//    path's through readFreqData.
//  - geneticPos is 0 for every site, because a VCF has no genetic-position
//    column.  Both consumers of geneticPos require --map: --weighted through
//    checkMapFile, and --cm through the mapfile check in main.  interpolate-
//    Geneticmap then overwrites it, so the zeros are never read.
//Free everything a reader has allocated but not yet handed off, on a throw.
//
//This replaces an inline
//    delete [] data; if (PHASED) delete [] firstCopy; if (USE_GL) delete [] glrow;
//that was repeated at seven throw sites.  A textual edit duplicated the last
//clause on four of them, which is a double free -- caught by ASan on the
//--gl-type PL --phased path.  One idempotent function instead: `delete [] NULL`
//is well defined, and the pointers are nulled, so calling it twice is safe.
//
//It also frees the rows ALREADY accumulated, which the inline cleanup did not.
//
//That does NOT make the error paths leak-free, and measurement says so: the
//residual is the partially built per-chromosome HapData/MapData/FreqData that
//main never releases when the reader throws.  loadTPEDData has exactly the
//same residual -- 409 allocations / 46,400 bytes on a malformed TPED against
//405 / 46,176 on a malformed VCF -- so it is a pre-existing property of the
//caller, not of this reader, and is left alone here.
//The accumulated rows alone, for throw sites reached before the per-locus
//pointers exist.
static void freeAccumulatedRows(vector< geno_t * > &hap, vector< bool * > &fc,
                                vector< double * > &gl)
{
    for (unsigned int i = 0; i < hap.size(); i++) delete [] hap[i];
    for (unsigned int i = 0; i < fc.size();  i++) delete [] fc[i];
    for (unsigned int i = 0; i < gl.size();  i++) delete [] gl[i];
    hap.clear(); fc.clear(); gl.clear();
}

static void abortRowRead(vector< geno_t * > &hap, vector< bool * > &fc,
                         vector< double * > &gl,
                         geno_t *&data, bool *&firstCopy, double *&glrow)
{
    delete [] data;      data      = NULL;
    delete [] firstCopy; firstCopy = NULL;
    delete [] glrow;     glrow     = NULL;
    freeAccumulatedRows(hap, fc, gl);
}

void loadVCFData(string vcffile, int &numLoci, int &numInd,
                 vector< HapData * > **hapDataByChr,
                 vector< MapData * > **mapDataByChr,
                 vector< FreqData * > **freqDataByChr,
                 vector< GenoLikeData * > **GLDataByChr,
                 int nresample, bool PHASED, bool AUTO_FREQ, bool PASS_ONLY,
                 string GL_TYPE, vector<string> &sampleIDs)
{
    igzstream fin;
    fin.open(vcffile.c_str());
    if (fin.fail())
    {
        LOG.err("ERROR: Failed to open", vcffile);
        throw 0;
    }
    if (!LOG.isQuiet()) cerr << "Reading " << vcffile << "\n";

    GarlicRNG *r = getRNG();

    const int VCF_FIXED = 9;   //CHROM POS ID REF ALT QUAL FILTER INFO FORMAT

    string line;
    //long long, not long: errlog has int/double/long long overloads and a
    //plain long is ambiguous between them.
    long long lineno = 0;
    bool haveHeader = false;

    string chr, locusName, ref, alt, filt, fmt;
    string emptyChr = "_nochr";
    string prevChr  = emptyChr;
    int currChrLoci = 0;

    vector<double> geneticPos;
    vector<pos_t>  physicalPos;
    vector<string> locusNames;
    vector<char>   allele;
    vector<double> freq;
    vector< geno_t * > hap;
    vector< bool * >   fc;
    vector< double * > gl;          //per-genotype error rates, when GL_TYPE is set

    const bool USE_GL = (GL_TYPE.compare("none") != 0 && !GL_TYPE.empty());

    map<string, int> chrSeen;

    long long nSkipIndel = 0, nSkipMulti = 0, nNonPass = 0, nSkipNonPass = 0, nSkipNoGT = 0;
    numLoci = 0;
    numInd  = 0;

    while (getline(fin, line))
    {
        lineno++;
        const char *p    = line.c_str();
        const char *pEnd = p + line.size();
        const char *tEnd;

        if (line.size() >= 2 && line[0] == '#' && line[1] == '#') continue;

        //--- the #CHROM header names the samples ---
        if (!line.empty() && line[0] == '#')
        {
            if (haveHeader)
            {
                LOG.err("ERROR: a second #CHROM header at line", lineno, false);
                LOG.err(" of", vcffile);
                freeAccumulatedRows(hap, fc, gl);
                throw 0;
            }
            int ncols = countFields(line);
            if (ncols < VCF_FIXED + 1)
            {
                LOG.err("ERROR:", vcffile, false);
                LOG.err(" has no sample columns; its #CHROM line has", ncols, false);
                LOG.err(" fields and at least", VCF_FIXED + 1, false);
                LOG.err(" are needed.");
                throw 0;
            }
            numInd = ncols - VCF_FIXED;
            //Sample IDs, in file order.  Every later per-sample diagnostic
            //names the sample rather than its column number.
            for (int k = 0; k < ncols; k++)
            {
                p = skipSpace(p, pEnd);
                tEnd = tokenEnd(p, pEnd);
                if (k >= VCF_FIXED) sampleIDs.push_back(string(p, tEnd - p));
                p = tEnd;
            }
            haveHeader = true;
            LOG.log("Samples in the VCF:", numInd);
            continue;
        }

        if (!haveHeader)
        {
            LOG.err("ERROR: data before the #CHROM header at line", lineno, false);
            LOG.err(" of", vcffile);
            freeAccumulatedRows(hap, fc, gl);
            throw 0;
        }
        if (line.empty()) continue;

        //--- CHROM ---
        p = skipSpace(p, pEnd);
        tEnd = tokenEnd(p, pEnd);
        if (tEnd == p)
        {
            LOG.err("ERROR: missing CHROM at line", lineno, false);
            LOG.err(" of", vcffile);
            throw 0;
        }
        chr.assign(p, tEnd - p);
        p = tEnd;

        //--- POS ---
        pos_t ppos;
        {
            char *q;
            p = skipSpace(p, pEnd);
            double v = strtod(p, &q);
            if (q == p)
            {
                LOG.err("ERROR: could not parse POS at line", lineno, false);
                LOG.err(" of", vcffile);
                freeAccumulatedRows(hap, fc, gl);
                throw 0;
            }
            p = q;
            ppos = pos_t(v);
        }

        //--- ID, REF, ALT, QUAL, FILTER, INFO, FORMAT ---
        p = skipSpace(p, pEnd); tEnd = tokenEnd(p, pEnd); locusName.assign(p, tEnd - p); p = tEnd;
        p = skipSpace(p, pEnd); tEnd = tokenEnd(p, pEnd); ref.assign(p, tEnd - p);       p = tEnd;
        p = skipSpace(p, pEnd); tEnd = tokenEnd(p, pEnd); alt.assign(p, tEnd - p);       p = tEnd;
        p = skipSpace(p, pEnd); tEnd = tokenEnd(p, pEnd);                                p = tEnd; //QUAL
        p = skipSpace(p, pEnd); tEnd = tokenEnd(p, pEnd); filt.assign(p, tEnd - p);      p = tEnd;
        p = skipSpace(p, pEnd); tEnd = tokenEnd(p, pEnd);                                p = tEnd; //INFO
        //errlog's (string, value) overloads insert a space, so composing a site
        //label from two calls rendered as "21: 13865210".  Build it once.
        string site;
        {
            stringstream ss;
            ss << chr << ":" << ppos;
            site = ss.str();
        }

        p = skipSpace(p, pEnd); tEnd = tokenEnd(p, pEnd); fmt.assign(p, tEnd - p);       p = tEnd;
        if (fmt.empty())
        {
            LOG.err("ERROR: line", lineno, false);
            LOG.err(" of", vcffile, false);
            LOG.err(" has fewer than", VCF_FIXED, false);
            LOG.err(" fixed fields.");
            freeAccumulatedRows(hap, fc, gl);
            throw 0;
        }

        //--- site filters, counted so the skips are never silent ---
        bool pass = (filt.compare("PASS") == 0 || filt.compare(".") == 0);
        if (!pass) nNonPass++;
        if (PASS_ONLY && !pass) { nSkipNonPass++; continue; }

        //Multiallelic first: it is the more specific reason, and an ALT of
        //"A,T" would otherwise be reported as an indel.
        if (alt.find(',') != string::npos) { nSkipMulti++; continue; }
        if (!isSNV(ref) || !isSNV(alt))    { nSkipIndel++; continue; }

        //--- GT's position in FORMAT.  Per-line: FORMAT may vary by site. ---
        int gtIndex = gtIndexOf(fmt);
        if (gtIndex < 0) { nSkipNoGT++; continue; }

        //The requested likelihood field.  Absent is an ERROR rather than a
        //skip: the user asked for per-genotype error rates, and quietly
        //dropping to --error for some sites would be exactly the silent
        //substitution the tgls guards exist to prevent.  The message names
        //whichever alternatives this site does carry.
        int glIndex = -1;
        if (USE_GL)
        {
            glIndex = formatIndexOf(fmt, GL_TYPE);
            if (glIndex < 0)
            {
                //Composed in one string rather than a chain of LOG.err calls:
                //the (string, value) overloads insert a space, and a trailing
                //`false` meant as the newline flag binds to err(string, bool)
                //and prints "FALSE".  Both bit this message before.
                string alt;
                const char *cand[3] = {"GQ", "PL", "GL"};
                int nalt = 0;
                for (int c = 0; c < 3; c++)
                {
                    if (GL_TYPE.compare(cand[c]) == 0) continue;
                    if (formatIndexOf(fmt, cand[c]) >= 0)
                    {
                        if (!alt.empty()) alt += " and ";
                        alt += cand[c];
                        nalt++;
                    }
                }
                stringstream ss;
                ss << "ERROR: --gl-type " << GL_TYPE << " was given, but FORMAT at "
                   << site << " is '" << fmt << "', which has no " << GL_TYPE << " field.";
                LOG.err(ss.str());
                stringstream s2;
                if (nalt == 1)
                    s2 << "\tThis site carries " << alt << "; use --gl-type " << alt
                       << " instead, or drop --gl-type to use a flat error rate from --error.";
                else if (nalt > 1)
                    s2 << "\tThis site carries " << alt << "; use --gl-type with one of those,"
                       << " or drop --gl-type to use a flat error rate from --error.";
                else
                    s2 << "\tNo genotype-quality field is present at all; drop --gl-type to use"
                       << " a flat error rate from --error.";
                LOG.err(s2.str());
                freeAccumulatedRows(hap, fc, gl);
                throw 0;
            }
        }

        //--- chromosome blocks.  Downstream code zips per-chromosome vectors
        //positionally and alignMapScaffold errors on a chromosome split into
        //non-contiguous blocks, so an interleaved VCF is rejected here with a
        //message that says what to do about it. ---
        if (prevChr.compare(emptyChr) == 0) prevChr = chr;

        if (chr.compare(prevChr) != 0)
        {
            if (chrSeen.count(chr) > 0)
            {
                LOG.err("ERROR: chromosome", chr, false);
                LOG.err(" reappears at line", lineno, false);
                LOG.err(" of", vcffile, false);
                LOG.err(" after another chromosome.");
                LOG.err("\tSites must be grouped by chromosome. Sort the VCF, e.g. with");
                LOG.err("\tbcftools sort, and try again.");
                freeAccumulatedRows(hap, fc, gl);
                throw 0;
            }
            chrSeen[prevChr] = 1;

            LOG.log("Chromosome", checkChrName(prevChr), false);
            LOG.log(":", currChrLoci, false);
            LOG.log(" sites.");

            (*mapDataByChr)->push_back(initMapData(geneticPos, physicalPos, locusNames, allele, currChrLoci, checkChrName(prevChr)));
            geneticPos.clear(); physicalPos.clear(); allele.clear(); locusNames.clear();

            (*hapDataByChr)->push_back(initHapData(hap, fc, currChrLoci, numInd, PHASED));
            hap.clear(); fc.clear();

            if (USE_GL)
            {
                (*GLDataByChr)->push_back(initGLData(gl, currChrLoci, numInd));
                gl.clear();
            }

            if (AUTO_FREQ)
            {
                (*freqDataByChr)->push_back(initFreqData(freq, currChrLoci));
                freq.clear();
            }

            prevChr = chr;
            currChrLoci = 0;
        }

        numLoci++;
        currChrLoci++;

        geneticPos.push_back(0.0);
        physicalPos.push_back(ppos);
        //A VCF's ID column is '.' at most sites.  A locus name is only used to
        //label output, and B3 was locus names being silently replaced by
        //positions, so an absent ID becomes an explicit CHROM:POS rather than
        //a file full of '.'.
        if (locusName.compare(".") == 0) locusNames.push_back(site);
        else locusNames.push_back(locusName);
        allele.push_back(alt[0]);

        //--- genotypes.  Raw rows, handed off to initHapData, exactly as in
        //loadTPEDData: the locus count is not known until the file ends. ---
        geno_t *data = new geno_t[numInd];
        bool *firstCopy = NULL;
        double *glrow = NULL;
        if (PHASED) firstCopy = new bool[numInd];
        if (USE_GL)  glrow = new double[numInd];

        double nalleles = 0, total = 0;

        for (int i = 0; i < numInd; i++)
        {
            p = skipSpace(p, pEnd);
            tEnd = tokenEnd(p, pEnd);
            if (tEnd == p)
            {
                abortRowRead(hap, fc, gl, data, firstCopy, glrow);
                LOG.err("ERROR: line", lineno, false);
                LOG.err(" of", vcffile, false);
                LOG.err(" has genotypes for fewer than", numInd, false);
                LOG.err(" samples.");
                throw 0;
            }

            int dosage, ploidy; bool fcopy, isPhased;
            if (!parseGT(p, tEnd, gtIndex, 1, dosage, fcopy, ploidy, isPhased))
            {
                string bad(p, tEnd - p);
                abortRowRead(hap, fc, gl, data, firstCopy, glrow);
                LOG.err("ERROR: could not parse the genotype of sample", sampleIDs[i], false);
                LOG.err(" at", site, false);
                LOG.err(" (field '", bad, false);
                LOG.err("', FORMAT '", fmt, false);
                LOG.err("').");
                throw 0;
            }
            if (ploidy != 2)
            {
                abortRowRead(hap, fc, gl, data, firstCopy, glrow);
                LOG.err("ERROR: sample", sampleIDs[i], false);
                LOG.err(" at", site, false);
                LOG.err(" has ploidy", ploidy, false);
                LOG.err("; garlic calls ROH from diploid genotypes only.");
                if (isSexChromosome(checkChrName(chr)))
                    LOG.err("\tThis is a sex chromosome: see --autosomes-only.");
                throw 0;
            }
            if (PHASED && !isPhased)
            {
                abortRowRead(hap, fc, gl, data, firstCopy, glrow);
                LOG.err("ERROR: --phased was given but sample", sampleIDs[i], false);
                LOG.err(" at", site, false);
                LOG.err(" is unphased ('/' rather than '|').");
                throw 0;
            }

            data[i] = geno_t(dosage);
            if (PHASED) firstCopy[i] = fcopy;

            if (USE_GL)
            {
                //The sub-field, located by index within this sample's column.
                const char *f = p;
                for (int k = 0; k < glIndex; k++)
                {
                    while (f < tEnd && *f != ':') f++;
                    if (f >= tEnd) break;
                    f++;
                }
                const char *fe = f;
                while (fe < tEnd && *fe != ':') fe++;
                string val(f, fe > f ? fe - f : 0);

                //genoIsCalled, not != GENO_MISSING: a HALF call has a genotype
                //code of its own now, is not usable as a genotype, and so needs
                //no per-genotype error rate.
                if (!genoIsCalled(geno_t(dosage)) || val.empty() || val.compare(".") == 0)
                {
                    //lod() takes its default branch for a missing genotype, so
                    //the error value is never read there.  1.0 is the honest
                    //placeholder: maximum uncertainty.
                    if (genoIsCalled(geno_t(dosage)))
                    {
                        abortRowRead(hap, fc, gl, data, firstCopy, glrow);
                        LOG.err("ERROR: sample", sampleIDs[i], false);
                        LOG.err(" at", site, false);
                        LOG.err(" has a called genotype but no", GL_TYPE, false);
                        LOG.err(" value. Drop --gl-type to use --error instead.");
                        throw 0;
                    }
                    glrow[i] = 1.0;
                }
                else if (GL_TYPE.compare("GQ") == 0)
                {
                    glrow[i] = glToError(atof(val.c_str()), "GQ");
                }
                else
                {
                    //PL or GL: the whole comma-separated array, converted to a
                    //posterior.  GL is log10 of a likelihood, so PL = -10*GL.
                    vector<double> arr;
                    size_t a = 0;
                    while (a <= val.size())
                    {
                        size_t b = val.find(',', a);
                        string tok = val.substr(a, (b == string::npos ? val.size() : b) - a);
                        double v = atof(tok.c_str());
                        arr.push_back(GL_TYPE.compare("GL") == 0 ? -10.0 * v : v);
                        if (b == string::npos) break;
                        a = b + 1;
                    }
                    //For a biallelic diploid site the VCF genotype ordering is
                    //0/0, 0/1, 1/1, so the called genotype's index IS the ALT
                    //dosage.  Only biallelic SNVs reach here.
                    if (int(arr.size()) != 3)
                    {
                        abortRowRead(hap, fc, gl, data, firstCopy, glrow);
                        LOG.err("ERROR: sample", sampleIDs[i], false);
                        LOG.err(" at", site, false);
                        LOG.err(" has", int(arr.size()), false);
                        LOG.err(" values in", GL_TYPE, false);
                        LOG.err("; a biallelic diploid site must have 3.");
                        throw 0;
                    }
                    glrow[i] = plToError(arr, dosage);
                }
            }

            p = tEnd;

            //Frequency over NON-MISSING alleles only, matching loadTPEDData,
            //which increments total only for an allele it could read.
            addAlleleCounts(geno_t(dosage), nalleles, total);
        }

        p = skipSpace(p, pEnd);
        if (p != pEnd)
        {
            abortRowRead(hap, fc, gl, data, firstCopy, glrow);
            LOG.err("ERROR: line", lineno, false);
            LOG.err(" of", vcffile, false);
            LOG.err(" has more columns than the", numInd, false);
            LOG.err(" samples named in the header.");
            throw 0;
        }

        hap.push_back(data);
        data = NULL;
        if (PHASED) { fc.push_back(firstCopy); firstCopy = NULL; }
        if (USE_GL) { gl.push_back(glrow); glrow = NULL; }

        if (AUTO_FREQ)
        {
                freq.push_back(alleleFrequency(nalleles, total, nresample, r));
        }
    }

    if (!haveHeader)
    {
        LOG.err("ERROR: no #CHROM header found in", vcffile, false);
        LOG.err(". Is it a VCF?");
        freeAccumulatedRows(hap, fc, gl);
        throw 0;
    }
    if (numLoci == 0)
    {
        LOG.err("ERROR: no usable sites in", vcffile, false);
        LOG.err(".");
        //Composed in one string: errlog's (string, bool) overload prints the
        //bool, so a trailing `false` meant as the newline flag rendered as
        //"Skipped: FALSE".
        {
            stringstream ss;
            ss << "\tgarlic uses biallelic SNVs only. Skipped: multiallelic " << nSkipMulti
               << ", non-SNV " << nSkipIndel << ", no GT " << nSkipNoGT
               << ", non-PASS " << nSkipNonPass << ".";
            LOG.err(ss.str());
        }
        freeAccumulatedRows(hap, fc, gl);
        throw 0;
    }

    LOG.log("Chromosome", checkChrName(prevChr), false);
    LOG.log(":", currChrLoci, false);
    LOG.log(" sites.");

    (*mapDataByChr)->push_back(initMapData(geneticPos, physicalPos, locusNames, allele, currChrLoci, checkChrName(prevChr)));
    geneticPos.clear(); physicalPos.clear(); allele.clear(); locusNames.clear();

    (*hapDataByChr)->push_back(initHapData(hap, fc, currChrLoci, numInd, PHASED));
    hap.clear(); fc.clear();

    if (USE_GL)
    {
        (*GLDataByChr)->push_back(initGLData(gl, currChrLoci, numInd));
        gl.clear();
    }

    if (AUTO_FREQ)
    {
        (*freqDataByChr)->push_back(initFreqData(freq, currChrLoci));
        freq.clear();
    }

    //Every skip is reported.  A VCF that is mostly indels, or whose FILTER
    //column is populated, silently losing most of its sites is exactly the
    //failure this tool should not have.
    if (nSkipMulti    > 0) LOG.log("Sites skipped, multiallelic:", nSkipMulti);
    if (nSkipIndel    > 0) LOG.log("Sites skipped, not a SNV:", nSkipIndel);
    if (nSkipNoGT     > 0) LOG.log("Sites skipped, no GT in FORMAT:", nSkipNoGT);
    if (nSkipNonPass  > 0) LOG.log("Sites skipped, FILTER not PASS:", nSkipNonPass);
    else if (nNonPass > 0) LOG.log("Sites with FILTER not PASS, KEPT (see --vcf-pass-only):", nNonPass);

    return;
}

bool parseGT(const char *sample, const char *sampleEnd, int gtIndex, int altIndex,
             int &dosage, bool &firstCopy, int &ploidy, bool &phased)
{
    //Walk to the gtIndex'th colon-separated sub-field.
    const char *f = sample;
    for (int k = 0; k < gtIndex; k++)
    {
        while (f < sampleEnd && *f != ':') f++;
        if (f >= sampleEnd) return false;       //fewer sub-fields than FORMAT promised
        f++;
    }
    const char *fEnd = f;
    while (fEnd < sampleEnd && *fEnd != ':') fEnd++;
    if (f == fEnd) return false;                //empty GT

    dosage    = 0;
    ploidy    = 0;
    phased    = true;
    firstCopy = false;

    bool anyMissing = false;
    int nObserved = 0;
    const char *q = f;
    while (q < fEnd)
    {
        //One allele: '.' or a run of digits.  The index may exceed one digit at
        //a multiallelic site.
        int idx;
        if (*q == '.') { idx = -1; q++; }
        else if (*q >= '0' && *q <= '9')
        {
            idx = 0;
            while (q < fEnd && *q >= '0' && *q <= '9') { idx = idx * 10 + (*q - '0'); q++; }
        }
        else return false;

        if (ploidy == 0 && idx >= 0) firstCopy = (idx == altIndex);
        ploidy++;

        if (idx < 0) anyMissing = true;
        else { nObserved++; if (idx == altIndex) dosage++; }

        if (q < fEnd)
        {
            if (*q == '|') q++;
            else if (*q == '/') { phased = false; q++; }
            else return false;
        }
        else break;
    }

    if (ploidy == 0) return false;
    if (anyMissing)
    {
        //A diploid HALF call keeps the allele that was observed: it is real
        //data and counts towards the allele frequency, even though the
        //genotype itself is unusable.  This is what loadTPEDData has always
        //done; the VCF path used to collapse it onto GENO_MISSING and throw
        //the observed allele away.
        if (ploidy == 2 && nObserved == 1)
            dosage = (dosage == 1 ? GENO_HALF_COUNTED : GENO_HALF_OTHER);
        else
            dosage = GENO_MISSING;
    }
    return true;
}

double plToError(const vector<double> &pl, int calledIndex)
{
    if (pl.empty()) return 1.0;
    if (calledIndex < 0 || calledIndex >= int(pl.size())) return 1.0;

    double m = pl[0];
    for (unsigned int i = 1; i < pl.size(); i++) if (pl[i] < m) m = pl[i];

    double sum = 0, best = 0;
    for (unsigned int i = 0; i < pl.size(); i++)
    {
        double p = pow(10.0, -(pl[i] - m) / 10.0);
        sum += p;
        if (int(i) == calledIndex) best = p;
    }
    if (sum <= 0) return 1.0;

    double e = 1.0 - best / sum;
    //Rounding can put e just outside [0,1].  The 1e-16 floor matches
    //glToError, so both paths hand lod() the same kind of value.
    if (e <= 0) e = 0.0000000000000001;
    if (e > 1)  e = 1;
    return e;
}

double glToError(double value, string glType)
{
    double e;
    if (glType.compare("GQ") == 0)
    {
        value /= (-10.0);
        value = (value > -10) ? value : -10;
        e = pow(10, value);
    }
    else if (glType.compare("GL") == 0)
    {
        value = (value > -10) ? value : -10;
        e = 1 - pow(10, value);
    }
    else if (glType.compare("PL") == 0)
    {
        value /= (-10.0);
        value = (value > -10) ? value : -10;
        e = 1 - pow(10, value);
    }
    else
    {
        //Unreachable: --gl-type is validated against {GQ, GL, PL} before any
        //file is read.  Return a neutral error rate rather than an
        //uninitialised value if that ever stops being true.
        LOG.err("ERROR: unknown --gl-type reached glToError:", glType);
        return 1.0;
    }

    if (e <= 0) e = 0.0000000000000001;
    if (e > 1) e = 1;
    return e;
}

string checkChrName(string chr) {
    if (chr[0] != 'c') {
        chr = "chr" + chr;
    }
    return chr;
}

void scanIndData3(string filename, int &numInd) {
    igzstream fin;
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        LOG.err("ERROR: Failed to open", filename);
        throw 0;
    }

    cout << "Reading " << filename << "\n";

    string line;
    int nind = 0;
    int min_cols = 2;
    int current_cols = 0;
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
        //The duplicate-ID check and the pooled-population warning that used to
        //be here are now in checkIndData, over the assembled IndData, so they
        //apply to every input path and not only to --tfam.
    }

    fin.close();

    numInd = nind;

    return;
}


void applyPopFile(const string &filename, IndData *indData)
{
    if (indData == NULL) return;

    igzstream fin;
    fin.open(filename.c_str());
    if (fin.fail())
    {
        LOG.err("ERROR: Failed to open", filename);
        throw 0;
    }
    if (!LOG.isQuiet()) cerr << "Reading population labels from " << filename << "\n";

    map<string, string> popOf;
    map<string, int>    sexOf;      //only for rows that supplied a third column
    vector<string>      col1, col2; //every row, in file order
    vector<string>      sexRaw;     //"" where the row had no third column
    vector<int>         lineOf;
    string line;
    int lineno = 0;

    while (getline(fin, line))
    {
        lineno++;
        //Trim and skip blanks and comments.
        size_t b = line.find_first_not_of(" \t\r\n");
        if (b == string::npos || line[b] == '#') continue;

        stringstream ss(line);
        string id, pop, sexField;
        if (!(ss >> id >> pop))
        {
            LOG.err("ERROR: line", lineno, false);
            LOG.err(" of", filename, false);
            LOG.err(" has fewer than 2 fields; --pop takes <sample_id> <population> [sex].");
            throw 0;
        }

        //pop is interpolated into a quoted UCSC track name in writeROHData, so
        //a double quote there would break the BED.
        if (pop.find('"') != string::npos)
        {
            LOG.err("ERROR: population label", pop, false);
            LOG.err(" on line", lineno, false);
            LOG.err(" of", filename, false);
            LOG.err(" contains a double quote, which would break the BED track name.");
            throw 0;
        }

        //The duplicate check is deliberately NOT here.  A TFAM-ordered file --
        //population first -- has the same population on every row, so checking
        //duplicates during the read would report "duplicate sample ID ( 36 )"
        //for the swapped-columns case and hide the diagnostic that explains it.
        //Rows are collected first; the swap check runs before the duplicate
        //check below.
        col1.push_back(id);
        col2.push_back(pop);
        sexRaw.push_back("");
        lineOf.push_back(lineno);

        if (ss >> sexField)
        {
            sexRaw[sexRaw.size() - 1] = sexField;
            if      (sexField.compare("1")  == 0) {}
            else if (sexField.compare("2")  == 0) {}
            else if (sexField.compare("0")  == 0 || sexField.compare("-9") == 0) {}
            else
            {
                LOG.err("ERROR: sex", sexField, false);
                LOG.err(" on line", lineno, false);
                LOG.err(" of", filename, false);
                LOG.err(" is not 1 (male), 2 (female), 0 or -9 (unknown).");
                throw 0;
            }
        }
    }
    fin.close();

    if (col1.empty())
    {
        LOG.err("ERROR: no usable rows in", filename);
        throw 0;
    }

    //Which of the data's samples does each column account for?  A TFAM-ordered
    //file -- population first -- is the mistake this format invites, and it
    //produces a confusing "sample not found" for every sample otherwise.
    map<string, int> haveID;
    for (int i = 0; i < indData->nind; i++) haveID[indData->indID[i]] = 1;

    int m1 = 0, m2 = 0;
    for (unsigned int i = 0; i < col1.size(); i++)
    {
        if (haveID.count(col1[i]) > 0) m1++;
        if (haveID.count(col2[i]) > 0) m2++;
    }
    if (m1 == 0 && m2 == int(col2.size()))
    {
        LOG.err("ERROR:", filename, false);
        LOG.err(" looks like its columns are swapped: no value in column 1 is a sample");
        LOG.err("\tin the data, but every value in column 2 is.");
        LOG.err("\t--pop takes <sample_id> <population>, unlike a TFAM, which is");
        LOG.err("\t<population> <sample_id>.");
        throw 0;
    }

    //Now that the swapped-columns case has been ruled out, duplicates in
    //column 1 really are duplicate sample IDs.
    for (unsigned int i = 0; i < col1.size(); i++)
    {
        if (popOf.count(col1[i]) > 0)
        {
            LOG.err("ERROR: Found duplicate sample ID ( ", col1[i], false);
            LOG.err(" ) in", filename);
            throw 0;
        }
        popOf[col1[i]] = col2[i];
        if (!sexRaw[i].empty())
            sexOf[col1[i]] = (sexRaw[i].compare("1") == 0) ? 1
                           : (sexRaw[i].compare("2") == 0) ? 2 : 0;
    }

    //Every sample must have a row.  Assigning a default would be the silent
    //pooling that the pooled-population warning exists to catch.
    int missing = 0;
    string firstMissing;
    for (int i = 0; i < indData->nind; i++)
    {
        if (popOf.count(indData->indID[i]) == 0)
        {
            if (missing == 0) firstMissing = indData->indID[i];
            missing++;
        }
    }
    if (missing > 0)
    {
        LOG.err("ERROR:", missing, false);
        LOG.err(" of", indData->nind, false);
        LOG.err(" samples have no row in", filename, false);
        LOG.err("; the first is", firstMissing, false);
        LOG.err(".");
        throw 0;
    }

    int changedPop = 0, changedSex = 0, setSex = 0;
    for (int i = 0; i < indData->nind; i++)
    {
        const string &id = indData->indID[i];
        if (indData->pop[i].compare(popOf[id]) != 0) changedPop++;
        indData->pop[i] = popOf[id];

        if (sexOf.count(id) > 0)
        {
            int s = sexOf[id];
            if (indData->sex[i] == 0 && s != 0) setSex++;
            else if (indData->sex[i] != s)      changedSex++;
            indData->sex[i] = s;
        }
    }

    LOG.log("Population labels read from", filename);
    if (changedPop > 0) LOG.log("--pop overrode the population label of", changedPop);
    if (setSex > 0)     LOG.log("--pop supplied a sex for", setSex);
    if (changedSex > 0) LOG.log("--pop overrode the sex of", changedSex);
    int extra = int(popOf.size()) - indData->nind;
    if (extra > 0)      LOG.log("rows in the file with no matching sample (ignored):", extra);

    return;
}

vector< pair<string, int> > enumeratePopulations(IndData *indData)
{
    vector< pair<string, int> > pops;
    if (indData == NULL) return pops;

    //map only to find a label again in O(log n); the ORDER is pops', which is
    //order of first appearance, not the map's (which would be sorted).
    map<string, unsigned int> indexOf;
    for (int i = 0; i < indData->nind; i++)
    {
        const string &p = indData->pop[i];
        map<string, unsigned int>::iterator it = indexOf.find(p);
        if (it == indexOf.end())
        {
            indexOf[p] = (unsigned int)(pops.size());
            pops.push_back(make_pair(p, 1));
        }
        else pops[it->second].second++;
    }
    return pops;
}

void checkIndData(IndData *indData, const string &source)
{
    if (indData == NULL) return;

    //Duplicate individual IDs.  Messages preserved verbatim from where this
    //lived in scanIndData3, including the spacing, so anything parsing them
    //sees no change.
    map<string, int> indList;
    for (int i = 0; i < indData->nind; i++)
    {
        string ind = indData->indID[i];
        if (indList.count(ind) > 0)
        {
            cerr << "ERROR: Found duplicate individual ID (" << ind << ") in " << source << endl;
            LOG.err("ERROR: Found duplicate individual ID ( ", ind, false);
            LOG.err(" ) in", source);
            throw 0;
        }
        else indList[ind] = 1;
    }

    //Nothing to warn about any more: several populations in one file are
    //analysed separately, each with its own allele frequencies, and
    //enumeratePopulations logs which ones were found.  The warning that used
    //to live here told the user to "run each population separately", which is
    //now what garlic does by itself.

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
        indData->pop[i] = pop;
        //Columns 3-5 of a TFAM are father, mother, sex.  Absent in some files,
        //so read them only if they are there and leave sex unknown otherwise.
        string pat, mat, sexField;
        if (ss >> pat >> mat >> sexField)
        {
            if (sexField == "1") indData->sex[i] = 1;
            else if (sexField == "2") indData->sex[i] = 2;
        }
        ss.clear();
    }
    fin.close();
    
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
    data->indID.assign(nind, "--");
    data->pop.resize(nind);
    data->sex.assign(nind, 0);

    return data;
}

DoubleData *initDoubleData(int n)
{
    DoubleData *data = new DoubleData;

    data->size = n;
    data->data.resize(n);

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
                double x = winDataByChr->at(chr)->data[ind][locus];
                //if (winDataByChr->at(chr)->data[ind][locus] == MISSING) nmiss++;
                //std::isnan, not isnan: isnan is a C99 MACRO in <math.h>, and
                //libstdc++'s <cmath> undefines it and provides only the
                //std:: form.  The unqualified spelling compiled here only
                //because devel's garlic-data.h included gsl/gsl_rng.h, which
                //pulled in the C <math.h>; dropping GSL in 0bdc179 removed
                //that transitive include and broke every GCC build.  libc++
                //keeps the global name, which is why macOS never noticed.
                if (x != MISSING && !(std::isnan(x)) ) size++;
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
                double x = winDataByChr->at(chr)->data[ind][locus];
                if (x != MISSING && !(std::isnan(x)) )
                {
                    rawWinData->data[i] = x;
                    i++;
                }
            }
        }
    }

    return rawWinData;
}

DoubleData *convertSubsetWinData2DoubleData(vector< WinData * > *winDataByChr, IndData *indData, int subsample, int step)
{
    GarlicRNG *r = getRNG();

    //to hold the indicies of the randomly selected individuals
    int nind = winDataByChr->at(0)->nind;
    vector<int> randInd;
    if (subsample >= nind)
    {
        randInd.resize(nind);
        for (int i = 0; i < nind; i++) randInd[i] = i;
    }
    else
    {
        vector<int> indIndex(nind);
        for (int i = 0; i < nind; i++) indIndex[i] = i;
        randInd.resize(subsample);
        r->choose(randInd.data(), subsample, indIndex.data(), nind);
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
                double x = winDataByChr->at(chr)->data[randInd[ind]][locus];
                if (x != MISSING && !(std::isnan(x))) size++;
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
                double x = winDataByChr->at(chr)->data[randInd[ind]][locus];
                if (x != MISSING && !(std::isnan(x)))
                {
                    rawWinData->data[i] = x;
                    i++;
                }
            }
        }
    }



    return rawWinData;
}

void releaseDoubleData(DoubleData *data)
{
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

void subsetDataByIndex(vector< HapData * > *hapDataByChr,
                       vector< GenoLikeData *> *GLDataByChr,
                       IndData *indData,
                       const vector<int> &keepInd,
                       vector< HapData * > **subsetHapDataByChr,
                       vector< GenoLikeData *> **subsetGLDataByChr,
                       IndData **subsetIndData,
                       bool USE_GL, bool PHASED)
{
    int nind = int(keepInd.size());

    //Checked once here rather than trusted per element: every read below is
    //data[locus][keepInd[ind]], so one bad index is an out-of-bounds read
    //repeated nloci times, and silent.
    int navail = (hapDataByChr->size() > 0 ? hapDataByChr->at(0)->nind : 0);
    for (int i = 0; i < nind; i++)
    {
        if (keepInd[i] < 0 || keepInd[i] >= navail || keepInd[i] >= indData->nind)
        {
            LOG.err("ERROR: individual index out of range in subsetDataByIndex:", keepInd[i]);
            throw 0;
        }
    }

    IndData *newIndData = initIndData(nind);

    //NB: assigning newIndData->pop = indData->pop here (as this used to do)
    //leaked the array initIndData just allocated AND aliased the caller's,
    //so the subsequent releaseIndData(subset) freed memory that main still
    //owned and freed again -> abort under --auto-winsize.  Deep copy instead.
    for (int ind = 0; ind < nind; ind++) {
        newIndData->indID[ind] = indData->indID[keepInd[ind]];
        newIndData->pop[ind]   = indData->pop[keepInd[ind]];
        newIndData->sex[ind]   = indData->sex[keepInd[ind]];
    }

    vector< HapData * > *newHapDataByChr = new vector< HapData * >;
    vector< GenoLikeData * > *newGLDataByChr = NULL;
    if (USE_GL) newGLDataByChr = new vector< GenoLikeData * >;

    int nchr = hapDataByChr->size();
    HapData *hapData;
    GenoLikeData *GLData = NULL;
    for (int chr = 0; chr < nchr; chr++)
    {
        int nloci = hapDataByChr->at(chr)->nloci;
        hapData = initHapData(nind, nloci, PHASED);
        if (USE_GL) GLData = initGLData(nind, nloci);

        for (int locus = 0; locus < nloci; locus++)
        {
            for (int ind = 0; ind < nind; ind++)
            {
                hapData->data[locus][ind] = hapDataByChr->at(chr)->data[locus][keepInd[ind]];
                if (PHASED) hapData->firstCopy[locus][ind] = hapDataByChr->at(chr)->firstCopy[locus][keepInd[ind]];
                if (USE_GL) GLData->data[locus][ind] = GLDataByChr->at(chr)->data[locus][keepInd[ind]];
            }
        }
        newHapDataByChr->push_back(hapData);
        hapData = NULL;
        if (USE_GL) newGLDataByChr->push_back(GLData);
        GLData = NULL;
    }

    *(subsetHapDataByChr) = newHapDataByChr;
    if(USE_GL) *(subsetGLDataByChr) = newGLDataByChr;
    *(subsetIndData) = newIndData;
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
    int nind = hapDataByChr->at(0)->nind;
    vector<int> randInd;
    if (subsample >= nind)
    {
        randInd.resize(nind);
        for (int i = 0; i < nind; i++) randInd[i] = i;
    }
    else
    {
        GarlicRNG *r = getRNG();
        vector<int> indIndex(nind);
        for (int i = 0; i < nind; i++) indIndex[i] = i;
        randInd.resize(subsample);
        r->choose(randInd.data(), subsample, indIndex.data(), nind);
    }

    subsetDataByIndex(hapDataByChr, GLDataByChr, indData, randInd,
                      subsetHapDataByChr, subsetGLDataByChr, subsetIndData,
                      USE_GL, PHASED);

    //Logged here, not in the gather: the gather is also how a population is
    //selected, and that has nothing to do with the KDE.
    LOG.loga("Individuals used for KDE:", (*subsetIndData)->indID.data(), int(randInd.size()));
    return;
}


