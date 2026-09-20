#include "garlic-errlog.h"
#include "countFeatures-cli.h"
#include <iostream>
#include <fstream>
#include "garlic-data.h"
//#include "garlic-roh.h"
#include "param_t.h"

using namespace std;

int int main(int argc, char const *argv[])
{
	param_t *params = getCLI(argc, argv);
	if (params == NULL) return 0;

	bool argerr = false;

	string outfile = params->getStringFlag(ARG_OUTFILE);

	string tpedfile = params->getStringFlag(ARG_TPED);
	string tfamfile = params->getStringFlag(ARG_TFAM);
	string rohfile = params->getStringFlag(ARG_ROHFILE);
	string featurefile = params->getStringFlag(ARG_FEATURES);

	argerr = checkRequiredFiles(tpedfile, tfamfile, rohfile, featurefile);
	LOG.log("TPED file:", tpedfile);
	char TPED_MISSING = params->getCharFlag(ARG_TPED_MISSING);
    LOG.log("TPED missing data code:", TPED_MISSING);
    LOG.log("TFAM file:", tfamfile);
    LOG.log("ROH file:", rohfile);
    LOG.log("Feature file:", featurefile);

    vector< int_pair_t > *chrCoordList;
    vector< MapData * > *mapDataByChr;
    IndData *indData;
    vector< HapData * > *hapDataByChr;

    try
    {
    	/*
        chrCoordList = scanTPEDMapData(tpedfile, numLoci, numCols);
        mapDataByChr = readTPEDMapData(tpedfile, numCols, chrCoordList, TPED_MISSING);

        LOG.log("Total loci:", numLoci);

        scanIndData3(tfamfile, numInd, popName);
        indData = readIndData3(tfamfile, numInd);

        LOG.log("Population:", popName);
        LOG.log("Total diploid individuals:", numInd);

        hapDataByChr = readTPEDHapData3(tpedfile, numLoci, numInd, TPED_MISSING, mapDataByChr);
		*/
        //Read in feature file

        //Read in ROH file


    }
    catch (...) { return 1; }

	return 0;
}