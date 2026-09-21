#include "garlic-features.h"
#include "garlic-errlog.h"
#include "gzstream.h"
#include <algorithm>
#include <sstream>
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
    pops[p].label  = popLabel;
    pops[p].bounds = bounds;

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
