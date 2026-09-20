#include "countFeatures-cli.h"
#include <iostream>
#include <string>
#include "param_t.h"
#include "garlic-errlog.h"

const string VERSION = "1.0.0";

const string ARG_OUTFILE = "--out";
const string DEFAULT_OUTFILE = "outfile";
const string HELP_OUTFILE = "The base name for all output files.";

const string ARG_TPED = "--tped";
const string DEFAULT_TPED = "none";
const string HELP_TPED = "A tped formatted file containing map and genotype information.";

const string ARG_TFAM = "--tfam";
const string DEFAULT_TFAM = "none";
const string HELP_TFAM = "A tfam formatted file containing population and individual IDs.";

const string ARG_TPED_MISSING = "--tped-missing";
const char DEFAULT_TPED_MISSING = '0';
const string HELP_TPED_MISSING = "Single character missing data code for TPED files.";

const string ARG_FEATURES = "--features";
const string DEFAULT_FEATURES = "none";
const string HELP_FEATURES = "A feature file formatted <chr> <pos> <allele> <class>.";

const string ARG_ROHFILE = "--roh";
const string DEFAULT_ROHFILE = "none";
const string HELP_ROHFILE = "An roh bed file output from GARLIC.";

param_t *getCLI(int argc, char *argv[])
{
	param_t *params = new param_t;
	params->addFlag(ARG_OUTFILE, DEFAULT_OUTFILE, "", HELP_OUTFILE);
	params->addFlag(ARG_TPED, DEFAULT_TPED, "", HELP_TPED);
	params->addFlag(ARG_TFAM, DEFAULT_TFAM, "", HELP_TFAM);
	params->addFlag(ARG_TPED_MISSING, DEFAULT_TPED_MISSING, "", HELP_TPED_MISSING);
	params->addFlag(ARG_FEATURES, DEFAULT_FEATURES, "", HELP_FEATURES);
	params->addFlag(ARG_ROHFILE, DEFAULT_ROHFILE, "", HELP_ROHFILE);

	params->setPreamble("countFeatures v" + VERSION);

	if (!params->parseCommandLine(argc, argv))
	{
		delete params; 
		return NULL;
	}
	return params;
}

bool checkRequiredFiles(string tpedfile, string tfamfile, string featurefile, string rohfile)
{
	if (tpedfile.compare(DEFAULT_TPED) == 0 || 
		tfamfile.compare(DEFAULT_TFAM) == 0 ||
		featurefile.compare(DEFAULT_FEATURE) == 0 ||
		rohfile.compare(DEFAULT_ROHFILE) == 0)
	{
		//cerr << "ERROR: Must provide both a tped, a tfam file.\n";
		LOG.err("ERROR: Must provide both a tped, a tfam, an roh, and a feature file.");
		return true;
	}
	return false;
}
