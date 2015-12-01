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


map<string, double> readLODCutoff(string lodCutoffFile, map<string, int> &pop2size)
{
    igzstream fin;
    fin.open(lodCutoffFile.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << lodCutoffFile << "\n";
        throw 0;
    }

    map<string, double> pop2lodcutoff;
    string line, popID, cutoffStr;
    double cutoff;
    stringstream ss;
    while (getline(fin, line))
    {
        ss.str(line);
        int fields = countFields(line);
        if (fields != 2)
        {
            cerr << "ERROR: Found " << fields << " fields but expected 2 in " << lodCutoffFile << endl;
            cerr << "\tFormat is <popID> <LOD cutoff>.\n";
            throw 0;
        }

        ss >> popID >> cutoffStr;

        if (pop2lodcutoff.count(popID) > 0)
        {
            cerr << "ERROR: Duplicate population ID (" << popID << ") found in " << lodCutoffFile << endl;
            throw 0;
        }

        if (!goodDouble(cutoffStr))
        {
            cerr << "ERROR: " << cutoffStr << " is not a valid double.\n";
            throw 0;
        }

        cutoff = atof(cutoffStr.c_str());

        if (pop2size.count(popID) > 0)
        {
            pop2lodcutoff[popID] = cutoff;
        }
        ss.clear();
    }

    if (pop2lodcutoff.size() != pop2size.size())
    {
        cerr << "ERROR: " << lodCutoffFile << " must provide one LOD score cutoff per population.\n";
        cerr << "\tExpected cutoffs for\n";

        for (map<string, int>::iterator it = pop2size.begin(); it != pop2size.end(); it++)
        {
            cerr << "\t" << it->first << endl;
        }
        cerr << "\tbut found only\n";
        for (map<string, double>::iterator it = pop2lodcutoff.begin(); it != pop2lodcutoff.end(); it++)
        {
            cerr << "\t" << it->first << endl;
        }
        throw 0;
    }

    fin.close();

    return pop2lodcutoff;
}

void readBoundSizes(string boundSizeFile, map<string, double> &pop2SMbound, map<string, double> &pop2MLbound, map<string, int> &pop2size)
{
    igzstream fin;
    fin.open(boundSizeFile.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << boundSizeFile << "\n";
        throw 0;
    }

    string line, popID, SMboundStr, MLboundStr;
    double SMbound, MLbound;
    stringstream ss;
    while (getline(fin, line))
    {
        //stringstream ss;
        ss.str(line);
        int fields = countFields(line);
        if (fields != 3)
        {
            cerr << "ERROR: Found " << fields << " fields but expected 3 in " << boundSizeFile << endl;
            cerr << "\tFormat is <popID> <small/medium size bound> <medium/long size bound>.\n";
            throw 0;
        }

        ss >> popID >> SMboundStr >> MLboundStr;

        if (pop2SMbound.count(popID) > 0)
        {
            cerr << "ERROR: Duplicate population ID (" << popID << ") found in " << boundSizeFile << endl;
            throw 0;
        }

        if (!goodDouble(SMboundStr))
        {
            cerr << "ERROR: " << SMboundStr << " is not a valid double.\n";
            throw 0;
        }

        SMbound = atof(SMboundStr.c_str());

        if (!goodDouble(MLboundStr))
        {
            cerr << "ERROR: " << MLboundStr << " is not a valid double.\n";
            throw 0;
        }

        MLbound = atof(MLboundStr.c_str());

        if (pop2size.count(popID) > 0)
        {
            double tmp;
            if (SMbound <= 0 || MLbound <= 0)
            {
                cerr << "ERROR: User provided size boundaries must be positive.\n";
                cerr << "\t" << popID << " " << SMbound << " " << MLbound << endl;
                throw 0;
            }
            else if (SMbound == MLbound)
            {
                cerr << "ERROR: Size boundaries must be different.\n";
                cerr << "\t" << popID << " " << SMbound << " " << MLbound << endl;
                throw 0;
            }
            else if (SMbound > MLbound)
            {
                tmp = MLbound;
                MLbound = SMbound;
                SMbound = tmp;
            }

            pop2SMbound[popID] = SMbound;
            pop2MLbound[popID] = MLbound;
        }
        ss.clear();
    }

    if (pop2SMbound.size() != pop2size.size())
    {
        cerr << "ERROR: " << boundSizeFile << " must provide size boundaries for each population.\n";
        cerr << "\tExpected cutoffs for\n";
        for (map<string, int>::iterator it = pop2size.begin(); it != pop2size.end(); it++)
        {
            cerr << "\t" << it->first << endl;
        }
        cerr << "\tbut found only\n";
        for (map<string, double>::iterator it = pop2SMbound.begin(); it != pop2SMbound.end(); it++)
        {
            cerr << "\t" << it->first << endl;
        }
        throw 0;
    }

    fin.close();

    return;
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

/*
vector< vector< FreqData * >* > *calcFreqData(vector< vector< HapData * >* > *hapDataByPopByChr, int nresample)
{
    const gsl_rng_type *T;
    gsl_rng *r;
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r, time(NULL));

    vector< vector< FreqData * >* > *freqDataByPopByChr = new vector< vector< FreqData * >* >;

    for (int pop = 0; pop < hapDataByPopByChr->size(); pop++)
    {
        vector< FreqData * > *freqDataByChr = new vector< FreqData * >;
        for (int chr = 0; chr < hapDataByPopByChr->at(pop)->size(); chr++)
        {
            FreqData *data = calcFreqData(hapDataByPopByChr->at(pop)->at(chr), nresample, r);
            freqDataByChr->push_back(data);
        }
        freqDataByPopByChr->push_back(freqDataByChr);
    }

    gsl_rng_free(r);

    return freqDataByPopByChr;
}
*/
vector< FreqData * > *calcFreqData2(vector< HapData * > *hapDataByChr, int nresample)
{
    const gsl_rng_type *T;
    gsl_rng *r;
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r, time(NULL));

    vector< FreqData * > *freqDataByChr = new vector< FreqData * >;
    for (int chr = 0; chr < hapDataByChr->size(); chr++)
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

void releaseFreqData(vector< vector< FreqData * >* > *freqDataByPopByChr)
{
    for (int pop = 0; pop < freqDataByPopByChr->size(); pop++)
    {
        for (int chr = 0; chr < freqDataByPopByChr->at(pop)->size(); chr++)
        {
            releaseFreqData(freqDataByPopByChr->at(pop)->at(chr));
        }
        freqDataByPopByChr->at(pop)->clear();
        delete freqDataByPopByChr->at(pop);
    }
    freqDataByPopByChr->clear();
    delete freqDataByPopByChr;
    return;
}

/*
void writeFreqData(string freqOutfile,
                   vector< vector< FreqData * >* > *freqDataByPopByChr,
                   vector< MapData * > *mapDataByChr,
                   vector< IndData * > *indDataByPop)
{
    freqOutfile += ".gz";
    ogzstream fout;
    fout.open(freqOutfile.c_str());

    if (fout.fail())
    {
        cerr << "ERROR: Failed to open " << freqOutfile << " for writing.\n";
        throw 0;
    }

    fout << "SNP\tALLELE\t";

    for (int pop = 0; pop < indDataByPop->size(); pop++)
    {
        fout << indDataByPop->at(pop)->pop << "\t";
    }

    fout << "\n";

    for (int chr = 0; chr < mapDataByChr->size(); chr++)
    {
        for (int locus = 0; locus < mapDataByChr->at(chr)->nloci; locus++)
        {
            fout << mapDataByChr->at(chr)->locusName[locus] << "\t"
                 << mapDataByChr->at(chr)->allele[locus] << "\t";
            for (int pop = 0; pop < indDataByPop->size(); pop++)
            {
                fout << freqDataByPopByChr->at(pop)->at(chr)->freq[locus] << "\t";
            }
            fout << "\n";
        }
    }
    cerr << "Wrote " << freqOutfile << endl;
    fout.close();
    return;
}
*/
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
        throw 0;
    }

    fout << "SNP\tALLELE\t" << popName << endl;

    for (int chr = 0; chr < mapDataByChr->size(); chr++)
    {
        for (int locus = 0; locus < mapDataByChr->at(chr)->nloci; locus++)
        {
            fout << mapDataByChr->at(chr)->locusName[locus] << "\t"
                 << mapDataByChr->at(chr)->allele[locus] << "\t"
                 << freqDataByChr->at(chr)->freq[locus] << "\n";
        }
    }
    cerr << "Wrote " << freqOutfile << endl;
    fout.close();
    return;
}

/*
vector< vector< FreqData * >* > *readFreqData(string freqfile,
        vector< int_pair_t > *chrCoordList,
        vector< MapData * > *mapDataByChr,
        map<string, int> &pop2index)
{
    //scan file for format integrity
    int expectedCols = pop2index.size() + 2;
    int expectedRows = 1;

    for (int chr = 0; chr < mapDataByChr->size(); chr++)
    {
        expectedRows += mapDataByChr->at(chr)->nloci;
    }

    igzstream fin;
    fin.open(freqfile.c_str());
    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << freqfile << " for reading.\n";
        throw 0;
    }
    cerr << "Checking " << freqfile << " integrity...\n";
    string line;
    int currentRows = 0;
    int currentCols = 0;
    int previousCols = -1;
    while (getline(fin, line))
    {
        currentRows++;
        currentCols = countFields(line);
        if (currentCols < expectedCols)
        {
            cerr << "ERROR: Found " << currentCols << " in " << freqfile
                 << " on line " << currentRows << " but expected at least "
                 << expectedCols << ".\n";
            throw 0;
        }
        if (currentCols != previousCols && previousCols != -1)
        {
            cerr << "ERROR: " << freqfile << " has differing number of columns across rows.\n";
            throw 0;
        }
        previousCols = currentCols;
    }

    if (currentRows != expectedRows)
    {
        cerr << "ERROR: " << freqfile << " has " << currentRows << " rows but expected "
             << expectedRows << ".\n";
        throw 0;
    }

    fin.close();
    fin.clear();

    //allocate
    vector< vector< FreqData * >* > *freqDataByPopByChr = new vector< vector< FreqData * >* >;

    for (int pop = 0; pop < pop2index.size(); pop++)
    {
        vector< FreqData * > *freqDataByChr = new vector< FreqData * >;
        for (int chr = 0; chr < mapDataByChr->size(); chr++)
        {
            FreqData *data = initFreqData(mapDataByChr->at(chr)->nloci);
            freqDataByChr->push_back(data);
        }
        freqDataByPopByChr->push_back(freqDataByChr);
    }

    fin.open(freqfile.c_str());
    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << freqfile << " for reading.\n";
        throw 0;
    }

    cerr << "Loading frequencies...\n";

    string header, junk;
    int popsFound = 0;
    getline(fin, header);
    int headerSize = countFields(header) - 2;
    if (headerSize <= 0)
    {
        cerr << "ERROR: " << freqfile << " header line has too few columns.\n";
        throw 0;
    }
    stringstream ss;
    ss.str(header);
    ss >> junk >> junk;
    string *popList = new string[headerSize];
    for (int i = 0; i < headerSize; i++)
    {
        ss >> popList[i];
        if (pop2index.count(popList[i]) > 0) popsFound++;
    }

    if (popsFound != pop2index.size())
    {
        cerr << "ERROR: " << freqfile << " must provide allele frequencies for all populations.\n";
        cerr << "\tExpected frequencies for\n";

        for (map<string, int>::iterator it = pop2index.begin(); it != pop2index.end(); it++)
        {
            cerr << "\t" << it->first << endl;
        }
        cerr << "\tbut found only\n";
        for (int i = 0; i < headerSize; i++)
        {
            if (pop2index.count(popList[i]) > 0)
            {
                cerr << "\t" << popList[i] << endl;
            }
        }
        throw 0;
    }


    string locusID;//, allele;
    char allele;
    double freq;
    for (int chr = 0; chr < mapDataByChr->size(); chr++)
    {
        for (int locus = 0; locus < mapDataByChr->at(chr)->nloci; locus++)
        {
            getline(fin, line);
            ss.clear();
            ss.str(line);
            ss >> locusID >> allele;
            if (mapDataByChr->at(chr)->locusName[locus].compare(locusID) != 0)
            {
                cerr << "ERROR: Loci appear out of order in " << freqfile << " relative to other files.\n";
                throw 0;
            }
            else
            {
                mapDataByChr->at(chr)->allele[locus] = allele;
            }
            for (int pop = 0; pop < headerSize; pop++)
            {
                ss >> freq;
                if (pop2index.count(popList[pop]) > 0)
                {
                    freqDataByPopByChr->at(pop2index[popList[pop]])->at(chr)->freq[locus] = freq;
                }
            }
        }
    }

    fin.close();

    return freqDataByPopByChr;
}
*/
vector< FreqData * > *readFreqData(string freqfile, string popName,
                                   vector< int_pair_t > *chrCoordList,
                                   vector< MapData * > *mapDataByChr)
{
    //scan file for format integrity
    int minCols = 3;
    int expectedRows = 1;

    for (int chr = 0; chr < mapDataByChr->size(); chr++)
    {
        expectedRows += mapDataByChr->at(chr)->nloci;
    }

    igzstream fin;
    fin.open(freqfile.c_str());
    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << freqfile << " for reading.\n";
        throw 0;
    }
    cerr << "Checking " << freqfile << " integrity...\n";
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
            throw 0;
        }
        if (currentCols != previousCols && previousCols != -1)
        {
            cerr << "ERROR: " << freqfile << " has differing number of columns across rows.\n";
            throw 0;
        }
        previousCols = currentCols;
    }

    if (currentRows != expectedRows)
    {
        cerr << "ERROR: " << freqfile << " has " << currentRows << " rows but expected "
             << expectedRows << ".\n";
        throw 0;
    }

    fin.close();
    fin.clear();

    //allocate
    vector< FreqData * > *freqDataByChr = new vector< FreqData * >;
    for (int chr = 0; chr < mapDataByChr->size(); chr++)
    {
        FreqData *data = initFreqData(mapDataByChr->at(chr)->nloci);
        freqDataByChr->push_back(data);
    }

    fin.open(freqfile.c_str());
    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << freqfile << " for reading.\n";
        throw 0;
    }

    cerr << "Loading frequencies...\n";

    string header, junk;
    int popsFound = 0;
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
        throw 0;
    }

    string locusID;//, allele;
    char allele;
    double freq;
    for (int chr = 0; chr < mapDataByChr->size(); chr++)
    {
        for (int locus = 0; locus < mapDataByChr->at(chr)->nloci; locus++)
        {
            fin >> locusID >> allele;
            if (mapDataByChr->at(chr)->locusName[locus].compare(locusID) != 0)
            {
                cerr << "ERROR: Loci appear out of order in " << freqfile << " relative to other files.\n";
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
                    freqDataByPopByChr->at(chr)->freq[locus] = freq;
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
    for (int i = 0; i < mapDataByChr->size(); i++)
    {
        releaseMapData(mapDataByChr->at(i));
    }
    mapDataByChr->clear();
    delete mapDataByChr;
    return;
}

void releaseIndData(vector< IndData * > *indDataByPop)
{
    for (int i = 0; i < indDataByPop->size(); i++)
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

    for (int pop = 0; pop < indCoordList->size(); pop++)
    {
        int totalHaps = 2 * (indCoordList->at(pop).second - indCoordList->at(pop).first + 1);
        vector< HapData * > *hapDataByChr = new vector< HapData * >;
        for (int chr = 0; chr < chrCoordList->size(); chr++)
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
    for (int chr = 0; chr < chrCoordList->size(); chr++)
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
            for (int pop = 0; pop < indCoordList->size(); pop++)
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
/*
vector< vector< HapData * >* > *readTPEDHapData2(string filename,
        int expectedLoci,
        int expectedInd,
        vector< int_pair_t > *chrCoordList,
        string *indList,
        map<string, string> &ind2pop,
        map<string, int> &pop2size,
        map<string, int> &pop2index,
        char TPED_MISSING,
        vector< MapData * > *mapDataByChr)
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
    //track index for individuals in hapDataByPopByChr->at(pop)->at(chr)->data
    map<string, int> pop2currind;

    //For each population
    for (map<string, int>::iterator it = pop2size.begin(); it != pop2size.end(); it++)
    {
        pop2currind[it->first] = 0; //init
        //number of haplotypes in current population
        int totalInd = it->second;
        vector< HapData * > *hapDataByChr = new vector< HapData * >;
        for (int chr = 0; chr < chrCoordList->size(); chr++)
        {
            int totalLoci = chrCoordList->at(chr).second - chrCoordList->at(chr).first + 1;
            HapData *data = initHapData(totalInd, totalLoci);
            hapDataByChr->push_back(data);
        }
        hapDataByPopByChr->push_back(hapDataByChr);
    }

    //string alleleStr1, alleleStr2;
    char alleleStr1, alleleStr2;
    string junk;//, oneAllele;
    char oneAllele, zeroAllele;
    //stringstream ss;
    //For each chromosome
    for (int chr = 0; chr < chrCoordList->size(); chr++)
    {
        int totalLoci = chrCoordList->at(chr).second - chrCoordList->at(chr).first + 1;
        //For each locus on the chromosome
        for (int locus = 0; locus < totalLoci; locus++)
        {
            for (map<string, int>::iterator it = pop2size.begin(); it != pop2size.end(); it++)
            {
                pop2currind[it->first] = 0;
            }
            //stringstream ss;
            oneAllele = mapDataByChr->at(chr)->allele[locus];
            zeroAllele = mapDataByChr->at(chr)->allele0[locus];
            //string genotypes;
            //getline(fin, genotypes);
            //ss.str(genotypes);
            /*
            ss >> junk;
            ss >> junk;
            ss >> junk;
            ss >> junk;
            */
fin >> junk;
fin >> junk;
fin >> junk;
fin >> junk;
for (int ind = 0; ind < expectedInd; ind++)
{
    string popStr = ind2pop[indList[ind]];
    int pop = pop2index[popStr];
    int index = pop2currind[popStr];
    //ss >> alleleStr1 >> alleleStr2;
    fin >> alleleStr1 >> alleleStr2;
    short allele1, allele2;
    allele1 = -1;
    allele2 = -1;

    //if (alleleStr1.compare(TPED_MISSING) == 0)
    if (alleleStr1 == TPED_MISSING)
    {
        allele1 = -9;
    }
    //else if (oneAllele.compare(TPED_MISSING) == 0)
    else if (oneAllele == TPED_MISSING)
    {
        oneAllele = alleleStr1;
        allele1 = 1;
    }
    //else if (alleleStr1.compare(oneAllele) == 0)
    else if (alleleStr1 == oneAllele)
    {
        allele1 = 1;
    }
    else
    {
        zeroAllele = alleleStr1;
        allele1 = 0;
    }

    //if (alleleStr2.compare(TPED_MISSING) == 0)
    if (alleleStr2 == TPED_MISSING)
    {
        allele2 = -9;
    }
    //else if (oneAllele.compare(TPED_MISSING) == 0)
    else if (oneAllele == TPED_MISSING)
    {
        oneAllele = alleleStr2;
        allele2 = 1;
    }
    //else if (alleleStr2.compare(oneAllele) == 0)
    else if (alleleStr2 == oneAllele)
    {
        allele2 = 1;
    }
    else
    {
        zeroAllele = alleleStr2;
        allele2 = 0;
    }

    //load into data stru
    if (allele1 + allele2 < 0)
    {
        hapDataByPopByChr->at(pop)->at(chr)->data[locus][index] = -9;
    }
    else
    {
        hapDataByPopByChr->at(pop)->at(chr)->data[locus][index] = allele1 + allele2;
    }
    //hapDataByPopByChr->at(pop)->at(chr)->data[2 * index + 1][locus] = allele2;
    pop2currind[ind2pop[indList[ind]]]++;
    //ss.clear();
}

mapDataByChr->at(chr)->allele[locus] = oneAllele;
mapDataByChr->at(chr)->allele0[locus] = zeroAllele;
}
}

fin.close();
return hapDataByPopByChr;
}
* /
vector< HapData * > *readTPEDHapData3(string filename,
                                      int expectedLoci,
                                      int expectedInd,
                                      char TPED_MISSING,
                                      vector< MapData * > *mapDataByChr)
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
    fin.clear();

    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    vector< HapData * > *hapDataByChr = new vector< HapData * >;
    for (int chr = 0; chr < mapDataByChr->size(); chr++)
    {
        HapData *data = initHapData(expectedInd, mapDataByChr->at(chr)->nloci);
        hapDataByChr->push_back(data);
    }

    cerr << "Loading genotypes " << filename << "...\n";

    //string alleleStr1, alleleStr2;
    char alleleStr1, alleleStr2;
    string junk;//, oneAllele;
    char oneAllele, zeroAllele;
    short allele;
    //For each chromosome
    for (int chr = 0; chr < mapDataByChr->size(); chr++)
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
        for (int chr = 0; chr < mapDataByChr->size(); chr++)
        {
            for (int locus = 0; locus < mapDataByChr->at(chr)->nloci; locus++)
            {
                fout << mapDataByChr->at(chr)->chr << "\t";
                fout << mapDataByChr->at(chr)->locusName[locus] << "\t";
                fout << mapDataByChr->at(chr)->geneticPos[locus] << "\t";
                fout << mapDataByChr->at(chr)->physicalPos[locus] << "\t";
                for (int ind = 0; ind < hapDataByPopByChr->at(pop)->at(chr)->nind; ind++)
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
        for (int ind = 0; ind < indDataByPop->at(pop)->nind; ind++)
        {
            fout << indDataByPop->at(pop)->pop << "\t";
            fout << indDataByPop->at(pop)->indID[ind] << "\t0\t0\t0\t0";
            fout << endl;
        }

        fout.close();
    }
    return;
}

/*
vector< vector< HapData * >* > *readHapData(string filename,
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
    int nhaps = 0;
    int nloci = -1;
    while (getline(fin, line))
    {
        nhaps++;
        nloci = countFields(line);
        //cout << "nloci: " << current_nloci << endl;
        if (nloci != expectedLoci)
        {
            cerr << "ERROR: line " << nhaps << " of " << filename << " has " << nloci
                 << ", but expected " << expectedLoci << ".\n";
            throw 0;
        }
    }
    if (nhaps != expectedHaps)
    {
        cerr << "ERROR: " << filename << " has " << nhaps
             << " haplotypes, but expected " << expectedHaps << ".\n";
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

    for (int pop = 0; pop < indCoordList->size(); pop++)
    {
        int totalHaps = 2 * (indCoordList->at(pop).second - indCoordList->at(pop).first + 1);
        vector< HapData * > *hapDataByChr = new vector< HapData * >;
        for (int chr = 0; chr < chrCoordList->size(); chr++)
        {
            int totalLoci = chrCoordList->at(chr).second - chrCoordList->at(chr).first + 1;
            HapData *data = initHapData(totalHaps, totalLoci);
            hapDataByChr->push_back(data);
        }
        hapDataByPopByChr->push_back(hapDataByChr);
    }

    //For each population
    for (int pop = 0; pop < indCoordList->size(); pop++)
    {
        int totalHaps = 2 * (indCoordList->at(pop).second - indCoordList->at(pop).first + 1);
        //For each haplotype in the population
        for (int hap = 0; hap < totalHaps; hap++)
        {
            //For each chromosome
            for (int chr = 0; chr < chrCoordList->size(); chr++)
            {
                int totalLoci = chrCoordList->at(chr).second - chrCoordList->at(chr).first + 1;
                //For each locus on the chromosome
                for (int locus = 0; locus < totalLoci; locus++)
                {
                    short allele;
                    fin >> allele;
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
                             << chr + 1 << chrPost << " chromosome has an illegal value.  Must be 0/1/-9.\n";

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
/*
vector< vector< HapData * >* > *readHapData2(string filename,
        int expectedLoci,
        int expectedInd,
        vector< int_pair_t > *chrCoordList,
        string *indList,
        map<string, string> &ind2pop,
        map<string, int> &pop2size,
        map<string, int> &pop2index)
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
    int nhaps = 0;
    int nloci = -1;
    while (getline(fin, line))
    {
        nhaps++;
        nloci = countFields(line);
        //cout << "nloci: " << current_nloci << endl;
        if (nloci != expectedLoci)
        {
            cerr << "ERROR: line " << nhaps << " of " << filename << " has " << nloci
                 << ", but expected " << expectedLoci << ".\n";
            throw 0;
        }
    }
    if (nhaps != expectedHaps)
    {
        cerr << "ERROR: " << filename << " has " << nhaps
             << " haplotypes, but expected " << expectedHaps << ".\n";
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

    //track index for individuals in hapDataByPopByChr->at(pop)->at(chr)->data
    map<string, int> pop2currind;
    vector< vector< HapData * >* > *hapDataByPopByChr = new vector< vector< HapData * >* >;

    //For each population
    for (map<string, int>::iterator it = pop2size.begin(); it != pop2size.end(); it++)
    {
        pop2currind[it->first] = 0; //init
        //number of haplotypes in current population
        int totalHaps = 2 * (it->second);
        vector< HapData * > *hapDataByChr = new vector< HapData * >;
        for (int chr = 0; chr < chrCoordList->size(); chr++)
        {
            int totalLoci = chrCoordList->at(chr).second - chrCoordList->at(chr).first + 1;
            HapData *data = initHapData(totalHaps, totalLoci);
            hapDataByChr->push_back(data);
        }
        hapDataByPopByChr->push_back(hapDataByChr);
    }

    //for each individual
    string hap1, hap2;
    short allele1, allele2;
    for (int ind = 0; ind < expectedInd; ind++)
    {
        stringstream ss1, ss2;
        getline(fin, hap1);
        getline(fin, hap2);
        ss1.str(hap1);
        ss2.str(hap2);

        int pop = pop2index[ind2pop[indList[ind]]];
        int index = pop2currind[ind2pop[indList[ind]]];
        //For each chromosome
        for (int chr = 0; chr < chrCoordList->size(); chr++)
        {
            //int totalLoci = chrCoordList->at(chr).second - chrCoordList->at(chr).first + 1;
            //For each locus on the chromosome
            for (int locus = 0; locus < hapDataByPopByChr->at(pop)->at(chr)->nloci; locus++)
            {
                ss1 >> allele1;
                ss2 >> allele2;
                //error checking
                if (allele1 != 0 && allele1 != 1 && allele1 != -9)
                {
                    cerr << "ERROR: Individual " << indList[ind] << " has an illegal value ("
                         << allele1 << ").  Must be 0/1/-9.\n";

                    throw 0;
                }
                if (allele2 != 0 && allele2 != 1 && allele2 != -9)
                {
                    cerr << "ERROR: Individual " << indList[ind] << " has an illegal value ("
                         << allele2 << ").  Must be 0/1/-9.\n";

                    throw 0;
                }
                //load into data stru
                hapDataByPopByChr->at(pop)->at(chr)->data[2 * index][locus] = allele1;
                hapDataByPopByChr->at(pop)->at(chr)->data[2 * index + 1][locus] = allele2;
            }
        }
        pop2currind[ind2pop[indList[ind]]]++;
    }

    fin.close();

    return hapDataByPopByChr;
}
*/

void releaseHapData(vector< vector< HapData * >* > *hapDataByPopByChr)
{
    for (int pop = 0; pop < hapDataByPopByChr->size(); pop++)
    {
        for (int chr = 0; chr < hapDataByPopByChr->at(pop)->size(); chr++)
        {
            releaseHapData(hapDataByPopByChr->at(pop)->at(chr));
        }
        hapDataByPopByChr->at(pop)->clear();
        delete hapDataByPopByChr->at(pop);
    }
    hapDataByPopByChr->clear();
    delete hapDataByPopByChr;
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
        cerr << "ERROR: Can't allocate WinData object.  Number of individuals (" << nind << ") and number of loci (" << nloci << ") must be positive.\n";
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
    for (int chr = 0; chr < winDataByChr->size(); chr++)
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

    for (int pop = 0; pop < indDataByPop->size(); pop++)
    {
        int nind = indDataByPop->at(pop)->nind;
        vector< WinData * > *winDataByChr = new vector< WinData * >;
        for (int chr = 0; chr < mapDataByChr->size(); chr++)
        {
            int nloci = mapDataByChr->at(chr)->nloci;
            WinData *data = initWinData(nind, nloci);
            winDataByChr->push_back(data);
        }
        winDataByPopByChr->push_back(winDataByChr);
    }

    return winDataByPopByChr;
}

vector< vector< WinData * >* > *initWinData(vector< MapData * > *mapDataByChr,
        vector< IndData * > *indDataByPop, int pop)
{
    vector< vector< WinData * >* > *winDataByPopByChr = new vector< vector< WinData * >* >;

    int nind = indDataByPop->at(pop)->nind;
    vector< WinData * > *winDataByChr = new vector< WinData * >;
    for (int chr = 0; chr < mapDataByChr->size(); chr++)
    {
        int nloci = mapDataByChr->at(chr)->nloci;
        WinData *data = initWinData(nind, nloci);
        winDataByChr->push_back(data);
    }
    winDataByPopByChr->push_back(winDataByChr);

    return winDataByPopByChr;
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
        cerr << "ERROR: Can't allocate HapData object.  Number of haplotypes (" << nind << ") and number of loci (" << nloci << ") must be positive.\n";
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
    cerr << "Scanning " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

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
    for (int i = 0; i < str.size(); i++) {
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
    cerr << "Loading map from " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    string junk;
    map<char, int> alleles;
    map<char, int>::iterator it;
    int count = 4;
    //For each chromosome
    for (int i = 0; i < chrCoordList->size(); i++)
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
            alleles.clear();
            count = 4;
            while (alleles.size() < 2 && count < numCols) {
                fin >> junk;
                if (junk.compare(TPED_MISSING) != 0) alleles[junk.c_str()[0]] = 0;
                count++;
            }

            if (count >= numCols) {
                cerr << "ERROR: locus " << data->locusName[locus] << " appears monomorphic.\n";
                throw 0;
            }
            it = alleles.begin();
            data->allele[locus] = it->first;
            /*
            it++;
            data->allele0[locus] = it->first;
            */
            getline(fin, junk);
        }
        data->chr = lc(data->chr);
        data->chr = checkChrName(data->chr);
        cerr << size << " loci on chromosome " << data->chr << endl;
        mapDataByChr->push_back(data);
    }

    fin.close();

    return mapDataByChr;
}


void scanIndData3(string filename, int &numInd, string &popName) {
    igzstream fin;
    cerr << "Scanning " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

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
            throw 0;
        }
        //stringstream ss;
        ss.str(line);
        ss >> pop >> ind;
        if (indList.count(ind) > 0)
        {
            cerr << "ERROR: Found duplicate individual ID (" << ind << ") in " << filename << endl;
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
            throw 0;
        }
        ss.clear();
    }

    fin.close();

    numInd = nind;

    return;
}


IndData *readIndData3(string filename, int numInd, string *indList)
{
    igzstream fin;
    cerr << "Loading IDs...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    IndData *data = initIndData(numInd);

    string line, pop, ind;
    stringstream ss;
    for (int i = 0; i < numInd; i++)
    {
        getline(fin, line);
        ss.str(line);
        ss >> pop >> ind;
        indList[i] = ind;
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
    for (int chr = 0; chr < winDataByChr->size(); chr++)
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
    for (int chr = 0; chr < winDataByChr->size(); chr++)
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

DoubleData *convertSubsetWinData2DoubleData(vector< WinData * > *winDataByChr, int subsample)
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

    DoubleData *rawWinData;
    int nmiss = 0;
    int ncols = 0;
    int nrows = 0;

    for (int chr = 0; chr < winDataByChr->size(); chr++)
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
    for (int chr = 0; chr < winDataByChr->size(); chr++)
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
    for (int pop = 0; pop < rawWinDataByPop->size(); pop++)
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

    for (int ind = 0; ind < nind; ind++) {
        newIndData->indID[ind] = indData->indID[randInd[ind]];
    }

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


