#include "garlic-data.h"

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
        cerr << "ERROR: number of loci (" << nloci << ") must be positive.\n";
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
        cerr << "ERROR: Failed to open " << freqOutfile << " for writing.\n";
        LOG.err("ERROR: Failed to open", freqOutfile);
        throw 0;
    }

    fout << "SNP\tALLELE\t" << popName << endl;

    for (unsigned int chr = 0; chr < mapDataByChr->size(); chr++)
    {
        for (int locus = 0; locus < mapDataByChr->at(chr)->nloci; locus++)
        {
            fout << mapDataByChr->at(chr)->locusName[locus] << "\t"
                 << mapDataByChr->at(chr)->allele[locus] << "\t"
                 << freqDataByChr->at(chr)->freq[locus] << "\n";
        }
    }
    cout << "Wrote allele frequency data to " << freqOutfile << endl;
    fout.close();
    return;
}

vector< FreqData * > *readFreqData(string freqfile, string popName,
                                   vector< int_pair_t > *chrCoordList,
                                   vector< MapData * > *mapDataByChr)
{
    //scan file for format integrity
    int minCols = 3;
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
        cerr << "ERROR: Failed to open " << freqfile << " for reading.\n";
        LOG.err("ERROR: Failed to open", freqfile);
        throw 0;
    }

    //cerr << "Loading frequencies...\n";

    string header, junk;
    //int popsFound = 0;
    getline(fin, header);
    int headerSize = countFields(header) - 2;
    stringstream ss;
    ss.str(header);
    ss >> junk >> junk;
    int popLocation = -1;
    string popStr;
    for (int i = 0; i < headerSize; i++)
    {
        ss >> popStr;
        if (popStr.compare(popName) == 0) popLocation = i;
    }

    if (popLocation < 0) {
        cerr << "ERROR: Could not find " << popName << " in " << freqfile << endl;
        LOG.err("ERROR: Could not find", popName, false);
        LOG.err(" in", freqfile);
        throw 0;
    }

    string locusID;//, allele;
    char allele;
    double freq;
    for (unsigned int chr = 0; chr < mapDataByChr->size(); chr++)
    {
        for (int locus = 0; locus < mapDataByChr->at(chr)->nloci; locus++)
        {
            fin >> locusID >> allele;
            if (mapDataByChr->at(chr)->locusName[locus].compare(locusID) != 0)
            {
                cerr << "ERROR: Loci appear out of order in " << freqfile << " relative to other files.\n";
                LOG.err("ERROR: Loci appear out of order in", freqfile);
                throw 0;
            }
            else
            {
                mapDataByChr->at(chr)->allele[locus] = allele;
            }

            for (int pop = 0; pop < headerSize; pop++)
            {
                fin >> freq;
                if (pop == popLocation)
                {
                    freqDataByChr->at(chr)->freq[locus] = freq;
                }
            }
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
    data->physicalPos = new int[nloci];
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

/*
vector< vector< HapData * >* > *readTPEDHapData(string filename,
        int expectedLoci,
        int expectedInd,
        vector< int_pair_t > *chrCoordList,
        vector< int_pair_t > *indCoordList)
{
    int expectedHaps = 2 * expectedInd;
    igzstream fin;
    cerr << "Checking " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

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
            throw 0;
        }
    }
    if (nloci != expectedLoci)
    {
        cerr << "ERROR: " << filename << " has " << nloci
             << " loci, but expected " << expectedLoci << ".\n";
        throw 0;
    }

    fin.close();
    fin.clear(); // clear error flags
    //fin.seekg(fileStart);
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    cerr << "Loading genotypes " << filename << "...\n";

    vector< vector< HapData * >* > *hapDataByPopByChr = new vector< vector< HapData * >* >;

    for (unsigned int pop = 0; pop < indCoordList->size(); pop++)
    {
        int totalHaps = 2 * (indCoordList->at(pop).second - indCoordList->at(pop).first + 1);
        vector< HapData * > *hapDataByChr = new vector< HapData * >;
        for (unsigned int chr = 0; chr < chrCoordList->size(); chr++)
        {
            int totalLoci = chrCoordList->at(chr).second - chrCoordList->at(chr).first + 1;
            HapData *data = initHapData(totalHaps, totalLoci);
            hapDataByChr->push_back(data);
        }
        hapDataByPopByChr->push_back(hapDataByChr);
    }

    string junk, oneAllele;
    string TPED_MISSING = "0";
    //For each chromosome
    for (unsigned int chr = 0; chr < chrCoordList->size(); chr++)
    {
        int totalLoci = chrCoordList->at(chr).second - chrCoordList->at(chr).first + 1;
        //For each locus on the chromosome
        for (int locus = 0; locus < totalLoci; locus++)
        {
            stringstream ss;
            oneAllele = TPED_MISSING;
            string genotypes;
            getline(fin, genotypes);
            genotypes += " ";
            ss.str(genotypes);
            ss >> junk;
            ss >> junk;
            ss >> junk;
            ss >> junk;
            //For each population
            for (unsigned int pop = 0; pop < indCoordList->size(); pop++)
            {
                int totalHaps = 2 * (indCoordList->at(pop).second - indCoordList->at(pop).first + 1);
                //For each haplotype in the population
                for (int hap = 0; hap < totalHaps; hap++)
                {
                    string alleleStr;
                    ss >> alleleStr;
                    short allele = -1;

                    if (alleleStr.compare(TPED_MISSING) == 0)
                    {
                        allele = -9;
                    }
                    else if (oneAllele.compare(TPED_MISSING) == 0)
                    {
                        oneAllele = alleleStr;
                        allele = 1;
                    }
                    else if (alleleStr.compare(oneAllele) == 0)
                    {
                        allele = 1;
                    }
                    else
                    {
                        allele = 0;
                    }

                    if (allele != 0 && allele != 1 && allele != -9)
                    {
                        string hapPost, popPost, chrPost, locPost;
                        hapPost = getPost(hap + 1);
                        popPost = getPost(pop + 1);
                        chrPost = getPost(chr + 1);
                        locPost = getPost(locus + 1);

                        cerr << "ERROR: The " << hap + 1 << hapPost << " haplotype in the "
                             << pop + 1 << popPost << " population at the "
                             << locus + 1 << locPost << " locus on the "
                             << chr + 1 << chrPost << " chromosome has an illegal value.\n";

                        throw 0;
                    }
                    hapDataByPopByChr->at(pop)->at(chr)->data[hap][locus] = allele;
                }
            }
        }
    }
    fin.close();

    return hapDataByPopByChr;
}
*/

vector< HapData * > *readTPEDHapData3(string filename,
                                      int expectedLoci,
                                      int expectedInd,
                                      char TPED_MISSING,
                                      vector< MapData * > *mapDataByChr)
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
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        LOG.err("ERROR: Failed to open", filename);
        throw 0;
    }

    vector< HapData * > *hapDataByChr = new vector< HapData * >;
    for (unsigned int chr = 0; chr < mapDataByChr->size(); chr++)
    {
        HapData *data = initHapData(expectedInd, mapDataByChr->at(chr)->nloci);
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
            }
        }
    }

    fin.close();
    return hapDataByChr;
}

/*
void writeTPEDDataByPop(string outfile, vector< vector< HapData * >* > *hapDataByPopByChr, vector< MapData * > *mapDataByChr, map<string, int> &pop2index)
{
    ogzstream fout;
    string tpedFile = outfile;
    //Each population
    for (map<string, int>::iterator it = pop2index.begin(); it != pop2index.end(); it++)
    {
        int pop = it->second;
        tpedFile = outfile + "." + it->first + ".tped.gz";
        fout.open(tpedFile.c_str());
        if (fout.fail())
        {
            cerr << "ERROR: Failed to open " << tpedFile << " for writing.\n";
            throw 0;
        }
        cerr << "Writing to " << tpedFile << endl;
        //each chromosome
        for (unsigned int chr = 0; chr < mapDataByChr->size(); chr++)
        {
            for (unsigned int locus = 0; locus < mapDataByChr->at(chr)->nloci; locus++)
            {
                fout << mapDataByChr->at(chr)->chr << "\t";
                fout << mapDataByChr->at(chr)->locusName[locus] << "\t";
                fout << mapDataByChr->at(chr)->geneticPos[locus] << "\t";
                fout << mapDataByChr->at(chr)->physicalPos[locus] << "\t";
                for (unsigned int ind = 0; ind < hapDataByPopByChr->at(pop)->at(chr)->nind; ind++)
                {
                    if (hapDataByPopByChr->at(pop)->at(chr)->data[locus][ind] == 2)
                    {
                        fout << mapDataByChr->at(chr)->allele[locus] << "\t" << mapDataByChr->at(chr)->allele[locus] << "\t";
                    }
                    else if (hapDataByPopByChr->at(pop)->at(chr)->data[locus][ind] == 1)
                    {
                        fout << mapDataByChr->at(chr)->allele0[locus] << "\t" << mapDataByChr->at(chr)->allele[locus] << "\t";
                    }
                    else
                    {
                        fout << mapDataByChr->at(chr)->allele0[locus] << "\t" << mapDataByChr->at(chr)->allele0[locus] << "\t";
                    }
                }
                fout << endl;
            }
        }

        fout.close();
    }
    return;
}
*/
/*
void writeTFAMDataByPop(string outfile, vector< IndData * > *indDataByPop, map<string, int> &pop2index)
{
    ogzstream fout;
    string tfamFile = outfile;
    //Each population
    for (map<string, int>::iterator it = pop2index.begin(); it != pop2index.end(); it++)
    {
        int pop = it->second;
        tfamFile = outfile + "." + it->first + ".tfam.gz";
        fout.open(tfamFile.c_str());
        if (fout.fail())
        {
            cerr << "ERROR: Failed to open " << tfamFile << " for writing.\n";
            throw 0;
        }
        cerr << "Writing to " << tfamFile << endl;
        //each chromosome
        for (unsigned int ind = 0; ind < indDataByPop->at(pop)->nind; ind++)
        {
            fout << indDataByPop->at(pop)->pop << "\t";
            fout << indDataByPop->at(pop)->indID[ind] << "\t0\t0\t0\t0";
            fout << endl;
        }

        fout.close();
    }
    return;
}
*/

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

vector< WinData * > *initWinData(vector< MapData * > *mapDataByChr, IndData *indData)
{
    int nind = indData->nind;
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

HapData *initHapData(unsigned int nind, unsigned int nloci)
{
    if (nind < 1 || nloci < 1)
    {
        cerr << "ERROR: Can not allocate HapData object.  Number of haplotypes (" << nind
             << ") and number of loci (" << nloci
             << ") must be positive.\n";
        LOG.err("ERROR: Can not allocate HapData object. Number of haplotypes (", int(nind), false);
        LOG.err(" ) and number of loci (", int(nloci), false);
        LOG.err(" ) must be positive.");
        throw 0;
    }

    HapData *data = new HapData;
    data->nind = nind;
    data->nloci = nloci;

    data->data = new short*[nloci];
    for (unsigned int i = 0; i < nloci; i++)
    {
        data->data[i] = new short[nind];
        for (unsigned int j = 0; j < nind; j++)
        {
            data->data[i][j] = MISSING;
        }
    }

    return data;
}

void releaseHapData(HapData *data)
{
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

DoubleData *convertWinData2DoubleData(vector< WinData * > *winDataByChr)
{
    int nmiss = 0;
    int ncols = 0;
    int nrows = 0;
    DoubleData *rawWinData;
    for (unsigned int chr = 0; chr < winDataByChr->size(); chr++)
    {
        for (int ind = 0; ind < winDataByChr->at(chr)->nind; ind++)
        {
            for (int locus = 0; locus < winDataByChr->at(chr)->nloci; locus++)
            {
                if (winDataByChr->at(chr)->data[ind][locus] == MISSING) nmiss++;
            }
        }

        ncols += winDataByChr->at(chr)->nloci;
        nrows = winDataByChr->at(chr)->nind;
    }
    rawWinData = initDoubleData(ncols * nrows - nmiss);

    int i = 0;
    for (unsigned int chr = 0; chr < winDataByChr->size(); chr++)
    {
        for (int ind = 0; ind < winDataByChr->at(chr)->nind; ind++)
        {
            for (int locus = 0; locus < winDataByChr->at(chr)->nloci; locus++)
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

DoubleData *convertSubsetWinData2DoubleData(vector< WinData * > *winDataByChr, IndData *indData, int subsample)
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

    cout << "Individuals used for KDE: ";
    LOG.logn("Individuals used for KDE: ");
    for (int ind = 0; ind < nind; ind++) {
        cout << indData->indID[randInd[ind]] << " ";
        LOG.logn(indData->indID[randInd[ind]]);
        LOG.logn(" ");
    }
    cout << "\n";
    LOG.logn("\n");

    DoubleData *rawWinData;
    int nmiss = 0;
    int ncols = 0;
    int nrows = 0;

    for (unsigned int chr = 0; chr < winDataByChr->size(); chr++)
    {
        for (int ind = 0; ind < nind; ind++)
        {
            for (int locus = 0; locus < winDataByChr->at(chr)->nloci; locus++)
            {
                if (winDataByChr->at(chr)->data[randInd[ind]][locus] == MISSING) nmiss++;
            }
        }

        ncols += winDataByChr->at(chr)->nloci;
        nrows = nind;
    }
    rawWinData = initDoubleData(ncols * nrows - nmiss);

    int i = 0;
    for (unsigned int chr = 0; chr < winDataByChr->size(); chr++)
    {
        for (int ind = 0; ind < nind; ind++)
        {
            for (int locus = 0; locus < winDataByChr->at(chr)->nloci; locus++)
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
                IndData *indData,
                vector< HapData * > **subsetHapDataByChr,
                IndData **subsetIndData,
                int subsample)
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

    cout << "Individuals used for KDE: ";
    for (int ind = 0; ind < nind; ind++) {
        newIndData->indID[ind] = indData->indID[randInd[ind]];
        cout << newIndData->indID[ind] << " ";
    }
    cout << "\n";
    LOG.loga("Individuals used for KDE:", newIndData->indID, nind);

    vector< HapData * > *newHapDataByChr = new vector< HapData * >;

    int nchr = hapDataByChr->size();
    HapData *hapData;
    for (int chr = 0; chr < nchr; chr++)
    {
        int nloci = hapDataByChr->at(chr)->nloci;
        hapData = initHapData(nind, nloci);

        for (int locus = 0; locus < nloci; locus++)
        {
            for (int ind = 0; ind < nind; ind++)
            {
                hapData->data[locus][ind] = hapDataByChr->at(chr)->data[locus][randInd[ind]];
            }
        }
        newHapDataByChr->push_back(hapData);
        hapData = NULL;
    }
    delete [] randInd;

    *(subsetHapDataByChr) = newHapDataByChr;
    *(subsetIndData) = newIndData;
    gsl_rng_free(r);
    return;
}


