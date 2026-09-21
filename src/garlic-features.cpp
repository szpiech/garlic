#include "garlic-features.h"
#include "garlic-errlog.h"
#include "gzstream.h"
#include <algorithm>
#include <sstream>
#include <fstream>
#include <set>
#include <limits>
#include <cstdlib>
#include <cerrno>

using namespace std;

//---- the classification table ----------------------------------------------

static string upperAscii(const string &s)
{
    string t;
    t.reserve(s.size());
    for (unsigned int i = 0; i < s.size(); i++)
    {
        char c = s[i];
        //By hand over 'a'..'z' for the same reason canonChrKey does it:
        //toupper() is locale-dependent and std::toupper(char) is undefined for
        //the negative char values every non-ASCII UTF-8 byte produces.
        if (c >= 'a' && c <= 'z') c = char(c - 'a' + 'A');
        t.push_back(c);
    }
    return t;
}

int FeatureTable::read(const string &filename)
{
    igzstream fin;
    fin.open(filename.c_str());
    if (fin.fail())
    {
        LOG.err("ERROR: Failed to open", filename);
        return -1;
    }

    string line;
    long long lineno = 0;
    //Class indices are handed out in first-appearance order here and remapped
    //to sorted order at the end, so the output columns are in a stable order
    //that does not depend on which class happens to appear first in the file.
    while (getline(fin, line))
    {
        lineno++;
        //A '#' anywhere before content is a comment, and that includes the
        //header line a user is likely to write.
        unsigned int i = 0;
        while (i < line.size() && (line[i] == ' ' || line[i] == '\t' ||
                                   line[i] == '\r')) i++;
        if (i >= line.size() || line[i] == '#') continue;

        int ncols = countFields(line);
        if (ncols != 4)
        {
            LOG.err("ERROR: line", lineno, false);
            LOG.err(" of", filename, false);
            LOG.err(" has", ncols, false);
            LOG.err(" columns; a feature file has 4 (chr, pos, allele, class).");
            //The scripts this replaces took "chr:pos ref alt class", so the
            //first field of an unconverted file carries a colon.  Say so
            //rather than letting the position parse fail three lines down.
            if (line.find(':') != string::npos)
                LOG.err("\tThis looks like the old count_features_in_roh.pl format. Convert it with:\n"
                        "\tawk '{split($1,a,\":\"); print a[1], a[2], $3, $4}' old.txt > new.txt");
            return -1;
        }

        string chr, posStr, allele, cls;
        {
            stringstream ss(line);
            ss >> chr >> posStr >> allele >> cls;
        }

        if (line.find(':') != string::npos && chr.find(':') != string::npos)
        {
            LOG.err("ERROR: the chromosome at line", lineno, false);
            LOG.err(" of", filename, false);
            LOG.err(" is", chr, false);
            LOG.err(", which contains a colon.");
            LOG.err("\tThis looks like the old count_features_in_roh.pl format. Convert it with:\n"
                    "\tawk '{split($1,a,\":\"); print a[1], a[2], $3, $4}' old.txt > new.txt");
            return -1;
        }

        pos_t pos;
        {
            errno = 0;
            char *end = NULL;
            long long v = strtoll(posStr.c_str(), &end, 10);
            if (end == posStr.c_str() || *end != '\0' || errno == ERANGE || v < 1)
            {
                LOG.err("ERROR: could not read a position from", posStr, false);
                LOG.err(" at line", lineno, false);
                LOG.err(" of", filename, false);
                LOG.err(". Positions are 1-based and positive.");
                return -1;
            }
            pos = pos_t(v);
        }

        int ci;
        {
            map<string, int>::iterator it = classIndex.find(cls);
            if (it != classIndex.end()) ci = it->second;
            else
            {
                ci = int(classNames.size());
                classNames.push_back(cls);
                classIndex[cls] = ci;
            }
        }

        vector<FeatureAllele> &v = byChr[canonChrKey(chr)][pos];
        string up = upperAscii(allele);
        for (unsigned int k = 0; k < v.size(); k++)
        {
            if (v[k].cls == ci && v[k].allele.compare(up) == 0)
            {
                LOG.err("ERROR: line", lineno, false);
                LOG.err(" of", filename, false);
                LOG.err(" repeats", chr, false);
                LOG.err(" ", posStr, false);
                LOG.err(" ", allele, false);
                LOG.err(" ", cls, false);
                LOG.err(".");
                LOG.err("\tA site may carry several classes, but not the same one twice.");
                return -1;
            }
        }
        FeatureAllele fa;
        fa.allele = up;
        fa.cls    = ci;
        v.push_back(fa);
        rows++;
    }
    fin.close();

    if (rows == 0)
    {
        LOG.err("ERROR:", filename, false);
        LOG.err(" gave no classified sites.");
        return -1;
    }

    //Sorted class order, so two cohorts annotated with the same scheme produce
    //the same columns whatever order the rows happened to be written in.
    vector<string> sorted = classNames;
    sort(sorted.begin(), sorted.end());
    vector<int> remap(classNames.size());
    for (unsigned int k = 0; k < classNames.size(); k++)
    {
        for (unsigned int j = 0; j < sorted.size(); j++)
            if (sorted[j].compare(classNames[k]) == 0) { remap[k] = int(j); break; }
    }
    map<string, map<pos_t, vector<FeatureAllele> > >::iterator c;
    for (c = byChr.begin(); c != byChr.end(); ++c)
    {
        map<pos_t, vector<FeatureAllele> >::iterator s;
        for (s = c->second.begin(); s != c->second.end(); ++s)
            for (unsigned int k = 0; k < s->second.size(); k++)
                s->second[k].cls = remap[s->second[k].cls];
    }
    classNames = sorted;
    classIndex.clear();
    for (unsigned int k = 0; k < classNames.size(); k++) classIndex[classNames[k]] = int(k);

    return 0;
}

const map<pos_t, vector<FeatureAllele> > *FeatureTable::chromosome(const string &chrKey) const
{
    map<string, map<pos_t, vector<FeatureAllele> > >::const_iterator it = byChr.find(chrKey);
    if (it == byChr.end()) return NULL;
    return &(it->second);
}

long long FeatureTable::nsites() const
{
    long long n = 0;
    map<string, map<pos_t, vector<FeatureAllele> > >::const_iterator it;
    for (it = byChr.begin(); it != byChr.end(); ++it) n += (long long)it->second.size();
    return n;
}

//---- the calls, indexed for lookup -----------------------------------------

int ROHIndex::addPopulation(const string &popLabel,
                            const vector< ROHData * > *rohDataByInd,
                            const vector< MapData * > *mapDataByChr,
                            const vector<double> &bounds,
                            const IndData *indData,
                            centromere *centro,
                            const ExcludedRegions *excluded,
                            const vector<ChrRole> *role)
{
    const int p = int(pops.size());
    pops.push_back(PopEntry());
    pops[p].label   = popLabel;
    pops[p].bounds  = bounds;
    pops[p].nclassN = int(bounds.size()) + 1;

    //What was analysed, per chromosome.  A chromosome carrying fewer than two
    //loci is not analysed -- writeFROH leaves it out of the denominator for
    //the same reason -- so it is left out here and its sites come back
    //UNASSESSED rather than NONE.
    for (unsigned int c = 0; c < mapDataByChr->size(); c++)
    {
        const MapData *md = mapDataByChr->at(c);
        if (md->nloci < 2) continue;
        ChrInfo ci;
        ci.lo   = md->physicalPos[0];
        ci.hi   = md->physicalPos[md->nloci - 1];
        ci.role = (role != NULL && c < role->size()) ? role->at(c) : CHR_AUTOSOME;
        if (centro != NULL)
        {
            ci.gapLo = centro->centromereStart(md->chr);
            ci.gapHi = centro->centromereEnd(md->chr);
        }
        if (excluded != NULL)
        {
            const vector<Interval> *iv = excluded->get(md->chr);
            if (iv != NULL) ci.excluded = *iv;
        }
        pops[p].chrInfo[canonChrKey(md->chr)] = ci;
    }

    for (unsigned int i = 0; i < rohDataByInd->size(); i++)
    {
        const ROHData *rd = rohDataByInd->at(i);
        const string &id = rd->indID;

        if (indOf.count(id) != 0)
        {
            //Duplicate IDs within one file are already refused by
            //checkIndData; this catches the same ID in two populations, which
            //nothing else looks for and which would make the counts of one
            //overwrite the other.
            LOG.err("ERROR: individual", id, false);
            LOG.err(" appears in more than one population.");
            LOG.err("\tCounting classified genotypes needs individual IDs to be unique across the run.");
            return -1;
        }

        const int ii = int(inds.size());
        inds.push_back(IndEntry());
        inds[ii].id   = id;
        inds[ii].pop  = p;
        inds[ii].zygo = (indData != NULL && i < (unsigned int)indData->nind)
                        ? indData->zygo[i] : ZYG_UNKNOWN;
        inds[ii].popDisplay = (indData != NULL && i < (unsigned int)indData->nind)
                              ? indData->pop[i] : popLabel;
        indOf[id] = ii;
        pops[p].members.push_back(ii);

        //Gathered per chromosome and sorted, rather than trusted to arrive in
        //order: bucketAt binary-searches, so the ordering is load-bearing.
        map<string, vector< pair<pos_t, pair<pos_t, int> > > > gather;
        for (unsigned int k = 0; k < rd->chr.size(); k++)
        {
            const MapData *md = mapDataByChr->at(rd->chr[k]);
            const int cls = rohSizeClassIndex(rd->length[k], bounds);
            gather[canonChrKey(md->chr)].push_back(
                make_pair(rd->start[k], make_pair(rd->stop[k], cls)));
        }
        map<string, vector< pair<pos_t, pair<pos_t, int> > > >::iterator g;
        for (g = gather.begin(); g != gather.end(); ++g)
        {
            sort(g->second.begin(), g->second.end());
            ChrTracts &t = inds[ii].byChr[g->first];
            t.start.reserve(g->second.size());
            t.stop.reserve(g->second.size());
            t.cls.reserve(g->second.size());
            for (unsigned int k = 0; k < g->second.size(); k++)
            {
                t.start.push_back(g->second[k].first);
                t.stop.push_back(g->second[k].second.first);
                t.cls.push_back(g->second[k].second.second);
            }
        }
    }

    return 0;
}

int ROHIndex::indexOf(const string &id) const
{
    map<string, int>::const_iterator it = indOf.find(id);
    return (it == indOf.end()) ? -1 : it->second;
}

void ROHIndex::resolveChr(const string &chrKey,
                          vector<const ChrTracts *> &tracts,
                          vector<const ChrInfo *> &info) const
{
    tracts.assign(inds.size(), (const ChrTracts *)NULL);
    info.assign(inds.size(), (const ChrInfo *)NULL);

    //One lookup per population rather than one per individual.
    vector<const ChrInfo *> byPop(pops.size(), (const ChrInfo *)NULL);
    for (unsigned int p = 0; p < pops.size(); p++)
    {
        map<string, ChrInfo>::const_iterator it = pops[p].chrInfo.find(chrKey);
        if (it != pops[p].chrInfo.end()) byPop[p] = &(it->second);
    }

    for (unsigned int i = 0; i < inds.size(); i++)
    {
        info[i] = byPop[inds[i].pop];
        map<string, ChrTracts>::const_iterator it = inds[i].byChr.find(chrKey);
        if (it != inds[i].byChr.end()) tracts[i] = &(it->second);
    }
}

//---- counting against calls that already exist -----------------------------

int ROHIndex::sizeClassIndexFromLabel(const string &label)
{
    if (label.empty()) return -1;
    long long idx = 0;
    for (unsigned int i = 0; i < label.size(); i++)
    {
        const char c = label[i];
        if (c < 'A' || c > 'Z') return -1;
        idx = idx * 26 + (c - 'A' + 1);
        if (idx > 1000000) return -1;       //not a label garlic would write
    }
    return int(idx - 1);
}

//Pulls the individual and population out of a track line.  The Perl's regex
//was /^track .+Ind: (.+) Pop:(.+) ROH.+/, which is greedy on both captures
//and misparses an ID that itself contains " Pop:"; this anchors on the LAST
//occurrence instead.  Returns false when the line is not a track line garlic
//wrote.
static bool parseTrackLine(const string &line, string &id, string &pop)
{
    const string key = "name=\"Ind: ";
    size_t a = line.find(key);
    if (a == string::npos) return false;
    a += key.size();
    size_t b = line.find('"', a);
    if (b == string::npos) return false;

    string inner = line.substr(a, b - a);           //<id> Pop:<pop> ROH
    const string tail = " ROH";
    if (inner.size() >= tail.size() &&
        inner.compare(inner.size() - tail.size(), tail.size(), tail) == 0)
        inner.erase(inner.size() - tail.size());

    const string sep = " Pop:";
    size_t c = inner.rfind(sep);
    if (c == string::npos) return false;
    id  = inner.substr(0, c);
    pop = inner.substr(c + sep.size());
    return !id.empty();
}

int ROHIndex::addFromBed(const string &bedfile, string &coordNote)
{
    igzstream fin;
    fin.open(bedfile.c_str());
    if (fin.fail())
    {
        LOG.err("ERROR: Failed to open", bedfile);
        return -1;
    }

    //One row as read, before the coordinate convention is known.
    struct Row { int ind; string chrKey; pos_t s, e; int cls; };
    vector<Row> rows;
    vector<string> idOrder, popOfInd;
    map<string, int> idIndex;

    //Which convention the file uses, counted rather than trusted.  A
    //conforming BED has chromEnd - chromStart == the length in column 5;
    //files written before d5946c4 wrote the 1-based start verbatim, so the
    //difference is one less.  The version in the track line cannot be used:
    //v1.1.6a was released with both conventions.
    long long nConforming = 0, nLegacy = 0, nUndecidable = 0;

    string line, id, pop;
    long long lineno = 0;
    int curInd = -1;
    int maxCls = 0;

    while (getline(fin, line))
    {
        lineno++;
        if (line.empty()) continue;
        if (line.compare(0, 6, "track ") == 0 || line.compare(0, 6, "track\t") == 0)
        {
            if (!parseTrackLine(line, id, pop))
            {
                LOG.err("ERROR: could not read the individual from the track line at line",
                        lineno, false);
                LOG.err(" of", bedfile);
                LOG.err("\tExpected: track name=\"Ind: <id> Pop:<pop> ROH\" ...");
                fin.close();
                return -1;
            }
            map<string, int>::iterator it = idIndex.find(id);
            if (it != idIndex.end())
            {
                //Two track blocks for one individual would make the second
                //silently replace or duplicate the first.
                LOG.err("ERROR:", bedfile, false);
                LOG.err(" has more than one track for individual", id, false);
                LOG.err(".");
                fin.close();
                return -1;
            }
            curInd = int(idOrder.size());
            idIndex[id] = curInd;
            idOrder.push_back(id);
            popOfInd.push_back(pop);
            continue;
        }
        if (line[0] == '#' || line.compare(0, 7, "browser") == 0) continue;
        if (curInd < 0)
        {
            LOG.err("ERROR: a call at line", lineno, false);
            LOG.err(" of", bedfile, false);
            LOG.err(" comes before any track line, so it belongs to no individual.");
            fin.close();
            return -1;
        }

        string chr, cls;
        double startv = 0, endv = 0, lenv = 0;
        {
            stringstream ss(line);
            if (!(ss >> chr >> startv >> endv >> cls >> lenv))
            {
                LOG.err("ERROR: line", lineno, false);
                LOG.err(" of", bedfile, false);
                LOG.err(" does not have the five columns a garlic .roh.bed has");
                LOG.err("\t(chrom, chromStart, chromEnd, class, length).");
                fin.close();
                return -1;
            }
        }

        const int k = sizeClassIndexFromLabel(cls);
        if (k < 0)
        {
            LOG.err("ERROR:", cls, false);
            LOG.err(" at line", lineno, false);
            LOG.err(" of", bedfile, false);
            LOG.err(" is not a size class label garlic writes (A, B, ... Z, AA, ...).");
            fin.close();
            return -1;
        }
        if (k > maxCls) maxCls = k;

        Row r;
        r.ind    = curInd;
        r.chrKey = canonChrKey(chr);
        r.s      = pos_t(startv);
        r.e      = pos_t(endv);
        r.cls    = k;
        rows.push_back(r);

        //Only a physical length can settle the convention.  Under --cm
        //column 5 is a genetic distance and says nothing about the
        //coordinates, which is why the answer is counted over rows rather
        //than taken from the first one.
        const pos_t span = r.e - r.s;
        if (lenv >= 1.0 && lenv == double(pos_t(lenv)))
        {
            const pos_t len = pos_t(lenv);
            if (span == len)          nConforming++;
            else if (span == len - 1) nLegacy++;
            else                      nUndecidable++;
        }
        else nUndecidable++;
    }
    fin.close();

    if (rows.empty())
    {
        LOG.err("ERROR:", bedfile, false);
        LOG.err(" contains no calls.");
        return -1;
    }

    //A 0-based start needs +1 to become the 1-based position a genotype file
    //reports; a legacy file already carries that position.
    int shift = 1;
    if (nConforming > 0 && nLegacy == 0)
    {
        shift = 1;
        coordNote = "0-based chromStart (garlic 2.0.0 and the later 1.1.6a builds)";
    }
    else if (nLegacy > 0 && nConforming == 0)
    {
        shift = 0;
        coordNote = "1-based chromStart (garlic 1.1.6a and earlier); read as written";
        LOG.log("The ROH file uses the pre-2.0.0 coordinate convention (chromEnd - chromStart is length - 1).");
    }
    else if (nConforming > 0 && nLegacy > 0)
    {
        LOG.err("ERROR:", bedfile, false);
        LOG.err(" is inconsistent about its coordinates:", nConforming, false);
        LOG.err(" call(s) have chromEnd - chromStart equal to the length column and", nLegacy, false);
        LOG.err(" have one less.");
        LOG.err("\tIt cannot be read without knowing which is meant; it looks like two files joined.");
        return -1;
    }
    else
    {
        //Every row was --cm, or the length column disagrees with both.
        shift = 1;
        coordNote = "assumed 0-based chromStart; the length column could not confirm it";
        LOG.log("Calls whose length column settles no convention:", nUndecidable);
        LOG.err("WARNING: the coordinate convention of " + bedfile + " could not be determined from");
        LOG.err("\tits length column, which is what a --cm .roh.bed looks like.  Assuming the");
        LOG.err("\tconforming 0-based chromStart that garlic 2.0.0 writes.  If the file came from");
        LOG.err("\ta build older than d5946c4, every interval is one base out.");
    }

    //One population per distinct label, in order of first appearance.  A
    //single label is unlabelled for naming, exactly as a single-population
    //run is.
    vector<string> popNames;
    vector<int> popOfIdx(idOrder.size(), 0);
    for (unsigned int i = 0; i < popOfInd.size(); i++)
    {
        int found = -1;
        for (unsigned int k = 0; k < popNames.size(); k++)
            if (popNames[k].compare(popOfInd[i]) == 0) { found = int(k); break; }
        if (found < 0) { found = int(popNames.size()); popNames.push_back(popOfInd[i]); }
        popOfIdx[i] = found;
    }

    for (unsigned int k = 0; k < popNames.size(); k++)
    {
        pops.push_back(PopEntry());
        pops.back().label   = (popNames.size() == 1) ? string("") : popNames[k];
        pops.back().nclassN = maxCls + 1;
        //No boundaries: a .roh.bed records the class letters, not what they
        //mean.  The header says so rather than inventing numbers.
    }

    for (unsigned int i = 0; i < idOrder.size(); i++)
    {
        inds.push_back(IndEntry());
        inds.back().id         = idOrder[i];
        inds.back().popDisplay = popOfInd[i];
        inds.back().pop        = popOfIdx[i];
        inds.back().zygo       = ZYG_UNKNOWN;
        indOf[idOrder[i]]      = int(i);
        pops[popOfIdx[i]].members.push_back(int(i));
    }

    //Every chromosome named in the file is "analysed" as far as anything here
    //can tell, over a span that covers everything: a .roh.bed does not record
    //where the markers reached, so no site can be shown to be UNASSESSED and
    //everything outside a run is NONE.  The header states this.
    //Gathered per individual so each one's tracts can be sorted.
    vector< map<string, vector< pair<pos_t, pair<pos_t, int> > > > > byInd(idOrder.size());
    set<string> chrSeen;
    for (unsigned int r = 0; r < rows.size(); r++)
    {
        const Row &row = rows[r];
        byInd[row.ind][row.chrKey].push_back(
            make_pair(row.s + shift, make_pair(row.e, row.cls)));
        chrSeen.insert(row.chrKey);
    }

    for (unsigned int p = 0; p < pops.size(); p++)
    {
        for (set<string>::const_iterator c = chrSeen.begin(); c != chrSeen.end(); ++c)
        {
            ChrInfo ci;
            ci.lo = 0;
            ci.hi = numeric_limits<pos_t>::max();
            ci.role = CHR_AUTOSOME;
            pops[p].chrInfo[*c] = ci;
        }
    }

    for (unsigned int i = 0; i < idOrder.size(); i++)
    {
        map<string, vector< pair<pos_t, pair<pos_t, int> > > >::iterator g;
        for (g = byInd[i].begin(); g != byInd[i].end(); ++g)
        {
            sort(g->second.begin(), g->second.end());
            ChrTracts &t = inds[i].byChr[g->first];
            for (unsigned int k = 0; k < g->second.size(); k++)
            {
                t.start.push_back(g->second[k].first);
                t.stop.push_back(g->second[k].second.first);
                t.cls.push_back(g->second[k].second.second);
            }
        }
    }

    builtFromBed = true;
    LOG.log("Individuals in the ROH file:", int(idOrder.size()));
    LOG.log("Calls read:", (long long)rows.size());
    LOG.log("ROH coordinates:", coordNote);
    return 0;
}

//---- the counts ------------------------------------------------------------

void FeatureCounts::init(int nind, int nclasses, int nsizeclass)
{
    nclass  = nsizeclass;
    nbucket = nsizeclass + 2;          //size classes, NONE, UNASSESSED
    hom.assign(nind, vector< vector<long long> >(nclasses, vector<long long>(nbucket, 0)));
    het = hom;
    n   = hom;
    missing.assign(nind, 0);
    hemi.assign(nind, 0);
    seen.assign(nind, false);
    rowsSeen = 0;
    rowsUnusable = 0;
    sitesNotFound = 0;
}

//Reports the individuals each file has and the other does not.  The Perl
//reported neither: an individual absent from the calls autovivified an empty
//interval list and every one of its homozygotes was counted as outside ROH.
static void reportUnmatched(const vector<string> &names, const string &what)
{
    if (names.empty()) return;
    ostringstream ss;
    ss << "WARNING: " << names.size() << " individual(s) " << what << ":";
    for (unsigned int i = 0; i < names.size() && i < 5; i++) ss << " " << names[i];
    if (names.size() > 5) ss << " ...";
    LOG.err(ss.str());
}

//A token of a whitespace-separated line.  garlic-data.cpp has the same two
//helpers as static inlines; they are three lines and not worth a header.
static inline const char *skipBlank(const char *p, const char *e)
{
    while (p < e && (*p == ' ' || *p == '\t' || *p == '\r')) p++;
    return p;
}

static inline const char *tokEnd(const char *p, const char *e)
{
    while (p < e && *p != ' ' && *p != '\t' && *p != '\r') p++;
    return p;
}

static inline char upChar(char c)
{
    return (c >= 'a' && c <= 'z') ? char(c - 'a' + 'A') : c;
}

int countFeaturesTPED(const string &tpedfile,
                      const string &tfamfile,
                      char TPED_MISSING,
                      const FeatureTable &features,
                      const ROHIndex &index,
                      FeatureCounts &counts)
{
    //The TFAM is read with garlic's own reader, so a malformed file and a
    //duplicate ID are refused with the messages they are refused with
    //everywhere else.
    IndData *tfam = NULL;
    try
    {
        int numInd = 0;
        scanIndData3(tfamfile, numInd);
        tfam = readIndData3(tfamfile, numInd);
    }
    catch (...)
    {
        logCurrentException("reading " + tfamfile);
        if (tfam != NULL) releaseIndData(tfam);
        return -1;
    }

    //Column -> individual in the call set.  This is the join the Perl did by
    //hash lookup with no miss handling.
    vector<int> colToInd(tfam->nind, -1);
    vector<string> notCalled;
    int nmatched = 0;
    for (int i = 0; i < tfam->nind; i++)
    {
        const int ri = index.indexOf(tfam->indID[i]);
        colToInd[i] = ri;
        if (ri < 0) notCalled.push_back(tfam->indID[i]);
        else nmatched++;
    }
    releaseIndData(tfam);

    if (nmatched == 0)
    {
        LOG.err("ERROR: no individual in", tfamfile, false);
        LOG.err(" has a run-of-homozygosity call in this run.");
        LOG.err("\tThe two files name no individual in common, so every classified genotype");
        LOG.err("\twould be counted as outside a run.  Check that the sample IDs match.");
        return -1;
    }
    reportUnmatched(notCalled, "in " + tfamfile + " have no ROH calls and are not counted");

    //Bucket columns are sized by the largest population's class count; a
    //population with fewer simply leaves the last columns at zero.
    int nsize = 0;
    for (int p = 0; p < index.npop(); p++)
        if (index.nclass(p) > nsize) nsize = index.nclass(p);
    counts.init(index.nind(), features.nclass(), nsize);
    //"The counting file carried a column for this individual", which is what
    //decides whether a row is written -- not "a classified site was reached",
    //which would drop everyone when the two files share no site.
    for (unsigned int i = 0; i < colToInd.size(); i++)
        if (colToInd[i] >= 0) counts.seen[colToInd[i]] = true;

    igzstream fin;
    fin.open(tpedfile.c_str());
    if (fin.fail())
    {
        LOG.err("ERROR: Failed to open", tpedfile);
        return -1;
    }

    string line, chr, lastChr;
    const map<pos_t, vector<FeatureAllele> > *featChr = NULL;
    vector<const ROHIndex::ChrTracts *> tracts;
    vector<const ROHIndex::ChrInfo *> info;
    long long lineno = 0, sitesFound = 0;
    const int nind = index.nind();

    while (getline(fin, line))
    {
        lineno++;
        const char *p    = line.c_str();
        const char *pEnd = p + line.size();
        const char *tEnd;

        p = skipBlank(p, pEnd);
        if (p == pEnd || *p == '#') continue;
        tEnd = tokEnd(p, pEnd);
        chr.assign(p, tEnd - p);
        p = tEnd;

        //Resolved once per chromosome rather than once per line.
        if (chr.compare(lastChr) != 0)
        {
            lastChr = chr;
            const string key = canonChrKey(chr);
            featChr = features.chromosome(key);
            index.resolveChr(key, tracts, info);
        }
        //Nothing classified on this chromosome: the genotype columns are
        //never touched, which is what keeps a whole-genome file cheap.
        if (featChr == NULL) continue;

        p = skipBlank(p, pEnd); tEnd = tokEnd(p, pEnd); p = tEnd;   //locus name
        p = skipBlank(p, pEnd); tEnd = tokEnd(p, pEnd); p = tEnd;   //genetic position

        pos_t ppos;
        {
            p = skipBlank(p, pEnd);
            char *q = NULL;
            const double v = strtod(p, &q);
            if (q == p)
            {
                LOG.err("ERROR: could not read a physical position at line", lineno, false);
                LOG.err(" of", tpedfile);
                fin.close();
                return -1;
            }
            p = q;
            ppos = pos_t(v);
        }

        map<pos_t, vector<FeatureAllele> >::const_iterator site = featChr->find(ppos);
        if (site == featChr->end()) continue;
        sitesFound++;

        //A TPED allele is a single character, as everywhere else in garlic
        //(loadTPEDData reads junk[0]).  A longer allele in the feature file --
        //an indel from a VCF-derived annotation -- cannot be matched against
        //one, and is reported rather than counted as absent.
        const vector<FeatureAllele> &rows = site->second;
        int usable = 0;
        for (unsigned int k = 0; k < rows.size(); k++)
            if (rows[k].allele.size() == 1) usable++;
        counts.rowsUnusable += (long long)(rows.size() - usable);
        counts.rowsSeen     += usable;
        if (usable == 0) continue;

        for (unsigned int i = 0; i < colToInd.size(); i++)
        {
            //Both allele columns are consumed whatever becomes of them, so
            //the columns stay in step with the samples.
            p = skipBlank(p, pEnd);
            if (p == pEnd)
            {
                LOG.err("ERROR: line", lineno, false);
                LOG.err(" of", tpedfile, false);
                LOG.err(" has genotypes for fewer individuals than", tfamfile, false);
                LOG.err(" names.");
                fin.close();
                return -1;
            }
            tEnd = tokEnd(p, pEnd);
            const char a1 = *p;
            p = tEnd;
            p = skipBlank(p, pEnd);
            if (p == pEnd)
            {
                LOG.err("ERROR: line", lineno, false);
                LOG.err(" of", tpedfile, false);
                LOG.err(" has an odd number of allele columns.");
                fin.close();
                return -1;
            }
            tEnd = tokEnd(p, pEnd);
            const char a2 = *p;
            p = tEnd;

            const int ri = colToInd[i];
            if (ri < 0) continue;

            const int bucket = ROHIndex::bucketAt(tracts[ri], info[ri],
                                                  index.zygoOf(ri), ppos);
            //One copy of the chromosome: the pair is one allele written
            //twice, so reading it as a homozygote is exactly the error this
            //branch exists to avoid.
            if (bucket == BUCKET_INELIGIBLE) { counts.hemi[ri] += usable; continue; }
            if (a1 == TPED_MISSING || a2 == TPED_MISSING)
            { counts.missing[ri] += usable; continue; }

            const int col = (bucket >= 0) ? bucket
                          : (bucket == BUCKET_NONE ? counts.colNone()
                                                   : counts.colUnassessed());
            const char u1 = upChar(a1), u2 = upChar(a2);
            for (unsigned int k = 0; k < rows.size(); k++)
            {
                if (rows[k].allele.size() != 1) continue;
                const char f = rows[k].allele[0];
                const int copies = (u1 == f ? 1 : 0) + (u2 == f ? 1 : 0);
                counts.n[ri][rows[k].cls][col]++;
                if (copies == 2)      counts.hom[ri][rows[k].cls][col]++;
                else if (copies == 1) counts.het[ri][rows[k].cls][col]++;
            }
        }

        p = skipBlank(p, pEnd);
        if (p != pEnd)
        {
            LOG.err("ERROR: line", lineno, false);
            LOG.err(" of", tpedfile, false);
            LOG.err(" has genotypes for more individuals than", tfamfile, false);
            LOG.err(" names.");
            fin.close();
            return -1;
        }
    }
    fin.close();

    counts.sitesNotFound = features.nsites() - sitesFound;

    vector<string> notCounted;
    for (int i = 0; i < nind; i++)
        if (!counts.seen[i]) notCounted.push_back(index.indID(i));
    reportUnmatched(notCounted, "have ROH calls but no genotypes in " + tpedfile);

    LOG.log("Classified sites found in the counting genotypes:", (long long)sitesFound);
    LOG.log("Classified sites not found:", counts.sitesNotFound);
    if (counts.rowsUnusable > 0)
        LOG.log("Classified rows skipped (allele longer than a TPED allele):", counts.rowsUnusable);

    return 0;
}

int countFeaturesVCF(const string &vcffile,
                     bool PASS_ONLY,
                     const FeatureTable &features,
                     const ROHIndex &index,
                     FeatureCounts &counts)
{
    igzstream fin;
    fin.open(vcffile.c_str());
    if (fin.fail())
    {
        LOG.err("ERROR: Failed to open", vcffile);
        return -1;
    }

    const int VCF_FIXED = 9;
    string line, chr, lastChr, ref, alt, filt, fmt;
    const map<pos_t, vector<FeatureAllele> > *featChr = NULL;
    vector<const ROHIndex::ChrTracts *> tracts;
    vector<const ROHIndex::ChrInfo *> info;
    vector<int> colToInd;
    vector<string> alleles;      //REF then each ALT, upper-cased
    vector<int> rowAllele;       //allele index per feature row at this site
    vector<int> needIdx;         //the distinct ones, usually a single entry
    long long lineno = 0, sitesFound = 0;
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
            vector<string> sampleIDs;
            {
                stringstream hs(line);
                string tok;
                for (int c = 0; c < VCF_FIXED && hs >> tok; c++) ;
                while (hs >> tok) sampleIDs.push_back(tok);
            }
            if (sampleIDs.empty())
            {
                LOG.err("ERROR:", vcffile, false);
                LOG.err(" has no sample columns.");
                fin.close();
                return -1;
            }

            colToInd.assign(sampleIDs.size(), -1);
            vector<string> notCalled;
            int nmatched = 0;
            for (unsigned int i = 0; i < sampleIDs.size(); i++)
            {
                const int ri = index.indexOf(sampleIDs[i]);
                colToInd[i] = ri;
                if (ri < 0) notCalled.push_back(sampleIDs[i]);
                else nmatched++;
            }
            if (nmatched == 0)
            {
                LOG.err("ERROR: no sample in", vcffile, false);
                LOG.err(" has a run-of-homozygosity call in this run.");
                LOG.err("\tThe two files name no individual in common, so every classified genotype");
                LOG.err("\twould be counted as outside a run.  Check that the sample IDs match.");
                fin.close();
                return -1;
            }
            reportUnmatched(notCalled, "in " + vcffile + " have no ROH calls and are not counted");

            int nsize = 0;
            for (int pp = 0; pp < index.npop(); pp++)
                if (index.nclass(pp) > nsize) nsize = index.nclass(pp);
            counts.init(index.nind(), features.nclass(), nsize);
            for (unsigned int i = 0; i < colToInd.size(); i++)
                if (colToInd[i] >= 0) counts.seen[colToInd[i]] = true;

            haveHeader = true;
            continue;
        }
        if (line.empty()) continue;
        if (!haveHeader)
        {
            LOG.err("ERROR: data before the #CHROM header at line", lineno, false);
            LOG.err(" of", vcffile);
            fin.close();
            return -1;
        }

        p = skipBlank(p, pEnd); tEnd = tokEnd(p, pEnd); chr.assign(p, tEnd - p); p = tEnd;
        if (chr.compare(lastChr) != 0)
        {
            lastChr = chr;
            const string key = canonChrKey(chr);
            featChr = features.chromosome(key);
            index.resolveChr(key, tracts, info);
        }
        if (featChr == NULL) continue;

        pos_t ppos;
        {
            p = skipBlank(p, pEnd);
            char *q = NULL;
            const double v = strtod(p, &q);
            if (q == p)
            {
                LOG.err("ERROR: could not parse POS at line", lineno, false);
                LOG.err(" of", vcffile);
                fin.close();
                return -1;
            }
            p = q;
            ppos = pos_t(v);
        }

        map<pos_t, vector<FeatureAllele> >::const_iterator site = featChr->find(ppos);
        if (site == featChr->end()) continue;

        p = skipBlank(p, pEnd); tEnd = tokEnd(p, pEnd);                          p = tEnd;  //ID
        p = skipBlank(p, pEnd); tEnd = tokEnd(p, pEnd); ref.assign(p, tEnd - p); p = tEnd;
        p = skipBlank(p, pEnd); tEnd = tokEnd(p, pEnd); alt.assign(p, tEnd - p); p = tEnd;
        p = skipBlank(p, pEnd); tEnd = tokEnd(p, pEnd);                          p = tEnd;  //QUAL
        p = skipBlank(p, pEnd); tEnd = tokEnd(p, pEnd); filt.assign(p, tEnd - p);p = tEnd;
        p = skipBlank(p, pEnd); tEnd = tokEnd(p, pEnd);                          p = tEnd;  //INFO
        p = skipBlank(p, pEnd); tEnd = tokEnd(p, pEnd); fmt.assign(p, tEnd - p); p = tEnd;

        //Honoured here as well as when calling, so one flag means one thing
        //about how this run treats a VCF record.  A site dropped by it is
        //reported as a classified site that was not found.
        if (PASS_ONLY && filt.compare("PASS") != 0 && filt.compare(".") != 0) continue;

        sitesFound++;

        //REF is allele 0 and each ALT follows, which is what a GT names.
        alleles.clear();
        alleles.push_back(upperAscii(ref));
        {
            string cur;
            for (unsigned int k = 0; k <= alt.size(); k++)
            {
                if (k == alt.size() || alt[k] == ',') { alleles.push_back(upperAscii(cur)); cur.clear(); }
                else cur.push_back(alt[k]);
            }
        }

        const vector<FeatureAllele> &rows = site->second;
        rowAllele.assign(rows.size(), -1);
        needIdx.clear();
        int usable = 0;
        for (unsigned int k = 0; k < rows.size(); k++)
        {
            for (unsigned int a = 0; a < alleles.size(); a++)
                if (alleles[a].compare(rows[k].allele) == 0) { rowAllele[k] = int(a); break; }
            if (rowAllele[k] < 0) continue;
            usable++;
            bool have = false;
            for (unsigned int j = 0; j < needIdx.size(); j++)
                if (needIdx[j] == rowAllele[k]) { have = true; break; }
            if (!have) needIdx.push_back(rowAllele[k]);
        }
        counts.rowsUnusable += (long long)(rows.size() - usable);
        counts.rowsSeen     += usable;

        const int gtIndex = gtIndexOf(fmt);
        if (gtIndex < 0)
        {
            //No genotypes to read, so the rows just counted cannot be
            //evaluated for anybody.  Undo rather than leave the per-individual
            //accounting unable to add up.
            counts.rowsSeen -= usable;
            counts.rowsUnusable += usable;
            continue;
        }
        if (usable == 0) continue;

        for (unsigned int i = 0; i < colToInd.size(); i++)
        {
            p = skipBlank(p, pEnd);
            if (p == pEnd)
            {
                LOG.err("ERROR: line", lineno, false);
                LOG.err(" of", vcffile, false);
                LOG.err(" has genotypes for fewer than the header's samples.");
                fin.close();
                return -1;
            }
            tEnd = tokEnd(p, pEnd);
            const char *sBeg = p, *sEnd = tEnd;
            p = tEnd;

            const int ri = colToInd[i];
            if (ri < 0) continue;

            const int bucket = ROHIndex::bucketAt(tracts[ri], info[ri],
                                                  index.zygoOf(ri), ppos);
            if (bucket == BUCKET_INELIGIBLE) { counts.hemi[ri] += usable; continue; }

            int dosage, ploidy; bool fcopy, phased;
            if (!parseGT(sBeg, sEnd, gtIndex, needIdx[0], dosage, fcopy, ploidy, phased))
            {
                LOG.err("ERROR: could not parse a genotype at line", lineno, false);
                LOG.err(" of", vcffile);
                fin.close();
                return -1;
            }
            //Not diploid here: one allele is not a genotype, and reading it
            //as one is how a hemizygous call becomes a spurious homozygote.
            if (ploidy != 2)     { counts.hemi[ri]    += usable; continue; }
            //Any allele of the call missing, including a half call: parseGT
            //returns one of the negative codes and there is no genotype.
            if (dosage < 0)      { counts.missing[ri] += usable; continue; }

            const int col = (bucket >= 0) ? bucket
                          : (bucket == BUCKET_NONE ? counts.colNone()
                                                   : counts.colUnassessed());
            for (unsigned int k = 0; k < rows.size(); k++)
            {
                if (rowAllele[k] < 0) continue;
                int copies = dosage;
                if (rowAllele[k] != needIdx[0])
                {
                    int d2, pl2; bool fc2, ph2;
                    if (!parseGT(sBeg, sEnd, gtIndex, rowAllele[k], d2, fc2, pl2, ph2))
                    {
                        LOG.err("ERROR: could not parse a genotype at line", lineno, false);
                        LOG.err(" of", vcffile);
                        fin.close();
                        return -1;
                    }
                    copies = d2;
                }
                counts.n[ri][rows[k].cls][col]++;
                if (copies == 2)      counts.hom[ri][rows[k].cls][col]++;
                else if (copies == 1) counts.het[ri][rows[k].cls][col]++;
            }
        }
    }
    fin.close();

    if (!haveHeader)
    {
        LOG.err("ERROR: no #CHROM header found in", vcffile, false);
        LOG.err(". Is it a VCF?");
        return -1;
    }

    counts.sitesNotFound = features.nsites() - sitesFound;

    vector<string> notCounted;
    for (int i = 0; i < index.nind(); i++)
        if (!counts.seen[i]) notCounted.push_back(index.indID(i));
    reportUnmatched(notCounted, "have ROH calls but no genotypes in " + vcffile);

    LOG.log("Classified sites found in the counting genotypes:", (long long)sitesFound);
    LOG.log("Classified sites not found:", counts.sitesNotFound);
    if (counts.rowsUnusable > 0)
        LOG.log("Classified rows skipped (allele at neither REF nor ALT, or no GT):",
                counts.rowsUnusable);

    return 0;
}

//---- the table -------------------------------------------------------------

static string bucketLabel(int col, int nclass)
{
    if (col < nclass) return sizeClassLabel(col);
    if (col == nclass) return "NONE";
    return "UNASSESSED";
}

int writeFeatureCounts(const string &outfile,
                       const FeatureTable &features,
                       const ROHIndex &index,
                       int pop,
                       const FeatureCounts &counts,
                       const FeatureRunInfo &runInfo)
{
    ofstream out;
    //ios::binary for the reason writeROHData gives: without it the Windows
    //CRT turns every \n into \r\n and the same run produces different bytes
    //on different platforms.
    out.open(outfile.c_str(), ios::binary);
    if (out.fail())
    {
        LOG.err("ERROR: Failed to open", outfile);
        return -1;
    }

    const int nclass = index.nclass(pop);
    const vector<string> &cls = features.classes();
    const vector<double> &bounds = index.bounds(pop);

    out << "## garlic feature counts\n";
    if (runInfo.pooled) out << "## populations_pooled\ttrue\n";
    out << "## feature_file\t" << runInfo.featureFile << "\n";
    out << "## feature_classes";
    for (unsigned int k = 0; k < cls.size(); k++) out << (k ? "," : "\t") << cls[k];
    out << "\n";
    out << "## feature_sites\t" << features.nsites() << "\tsites, "
        << features.nrows() << " classified rows\n";
    out << "## feature_rows_seen\t" << counts.rowsSeen
        << "\tclassified rows found in the genotypes; "
        << counts.sitesNotFound << " sites were not found\n";
    //Said plainly because it is the difference between this table and one
    //made from garlic's internal matrix, which would have dropped exactly the
    //rare sites this counts.
    out << "## genotype_source\t" << runInfo.genotypeSource
        << "\traw calls; no garlic site filter applied\n";
    if (!runInfo.rohSource.empty())
    {
        out << "## roh_source\t" << runInfo.rohSource
            << "\tcalls read from this file; not made by this run\n";
        out << "## roh_coordinates\t" << runInfo.coordNote << "\n";
    }
    if (index.fromBed())
        //Not left blank: a reader comparing two tables has to know that a
        //class letter here is only a letter, with no length behind it.
        out << "## size_class_boundaries\tunknown\ta .roh.bed records the class "
            << "letters, not the boundaries they came from\n";
    else
    {
        out << "## size_class_boundaries";
        for (unsigned int k = 0; k < bounds.size(); k++) out << "\t" << bounds[k];
        out << "\n";
    }
    out << "## buckets";
    for (int c = 0; c < nclass + 2; c++) out << (c ? "," : "\t") << bucketLabel(c, nclass);
    out << ",ALL\n";
    out << "## metrics\thom = homozygous for the classified allele; "
        << "het = one copy; n = genotyped sites\n";
    if (index.fromBed())
        //The honest statement of what is lost by counting against a file.
        //In a run that makes its own calls these are different things.
        out << "## unassessed\tnot determinable from a .roh.bed, which records where the runs "
            << "are and not which chromosomes were analysed; every site outside a run is NONE\n";
    else
        out << "## unassessed\ta site where no run could have been called: a chromosome "
            << "not analysed, an assembly gap, an excluded region, or outside the analysed span\n";
    out << "## in_roh\tALL minus NONE minus UNASSESSED\n";
    //An individual with one copy of a chromosome has no diploid genotype
    //there, so those sites are in no bucket.  Counted rather than dropped
    //silently, because in a TPED a hemizygous call looks like a homozygote.
    out << "## n_hemi\tclassified rows where the individual is not diploid; counted in no bucket\n";
    out << "## n_missing\tclassified rows the individual has no call for\n";

    out << "ind\tpop\tzygosity\tn_missing\tn_hemi";
    for (unsigned int k = 0; k < cls.size(); k++)
        for (int c = 0; c < nclass + 3; c++)
        {
            const string b = (c < nclass + 2) ? bucketLabel(c, nclass) : string("ALL");
            out << "\t" << cls[k] << "_" << b << "_hom"
                << "\t" << cls[k] << "_" << b << "_het"
                << "\t" << cls[k] << "_" << b << "_n";
        }
    out << "\n";

    const vector<int> &members = index.indsOf(pop);
    for (unsigned int m = 0; m < members.size(); m++)
    {
        const int i = members[m];
        //No genotypes for this individual: a row of zeros would say it has no
        //classified homozygotes, which is not what was measured.
        if (!counts.seen[i]) continue;

        const int z = index.zygoOf(i);
        out << index.indID(i) << "\t" << index.popDisplay(i) << "\t"
            << (z == ZYG_HOMOGAMETIC ? "homogametic"
               : (z == ZYG_HETEROGAMETIC ? "heterogametic" : "unknown"))
            << "\t" << counts.missing[i] << "\t" << counts.hemi[i];

        for (unsigned int k = 0; k < cls.size(); k++)
        {
            long long aHom = 0, aHet = 0, aN = 0;
            for (int c = 0; c < nclass + 2; c++)
            {
                out << "\t" << counts.hom[i][k][c]
                    << "\t" << counts.het[i][k][c]
                    << "\t" << counts.n[i][k][c];
                aHom += counts.hom[i][k][c];
                aHet += counts.het[i][k][c];
                aN   += counts.n[i][k][c];
            }
            out << "\t" << aHom << "\t" << aHet << "\t" << aN;
        }
        out << "\n";
    }

    out.close();
    LOG.log("Classified genotype counts:", outfile);
    return 0;
}

int ROHIndex::bucketAt(const ChrTracts *t, const ChrInfo *ci, int zygo, pos_t pos)
{
    //Not analysed at all: a chromosome removed by --chr or --autosomes-only,
    //one carrying fewer than two loci, or one the calling data never had.
    if (ci == NULL) return BUCKET_UNASSESSED;

    //A run of homozygosity cannot be called for an individual who has one
    //copy of the chromosome, so nothing here is evidence about them either
    //way.  Kept distinct from UNASSESSED because the caller must also refuse
    //to read the genotype: a hemizygous call written as a doubled allele
    //looks exactly like a homozygote.
    if (!eligibleForCalling(ci->role, zygo)) return BUCKET_INELIGIBLE;

    //Outside the span the markers cover, inside the assembly gap, or inside a
    //region dropped before calling: no run could have been called here.
    if (pos < ci->lo || pos > ci->hi) return BUCKET_UNASSESSED;
    if (ci->gapHi > ci->gapLo && pos >= ci->gapLo && pos <= ci->gapHi)
        return BUCKET_UNASSESSED;
    for (unsigned int k = 0; k < ci->excluded.size(); k++)
        if (pos >= ci->excluded[k].start && pos <= ci->excluded[k].end)
            return BUCKET_UNASSESSED;

    if (t == NULL || t->start.empty()) return BUCKET_NONE;

    //The tracts are sorted by start and disjoint, so the only candidate is the
    //last one that begins at or before pos.
    vector<pos_t>::const_iterator it =
        upper_bound(t->start.begin(), t->start.end(), pos);
    if (it == t->start.begin()) return BUCKET_NONE;
    const size_t k = size_t(it - t->start.begin()) - 1;
    if (pos <= t->stop[k]) return t->cls[k];
    return BUCKET_NONE;
}
