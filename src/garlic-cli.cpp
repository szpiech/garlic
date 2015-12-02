#include "garlic-cli.h"

param_t *getCLI(int argc, char *argv[]) {
	param_t *params = new param_t;
	params->addFlag(ARG_OUTFILE, DEFAULT_OUTFILE, "", HELP_OUTFILE);
	params->addFlag(ARG_THREADS, DEFAULT_THREADS, "", HELP_THREADS);
	params->addFlag(ARG_ERROR, DEFAULT_ERROR, "", HELP_ERROR);
	params->addFlag(ARG_WINSIZE, DEFAULT_WINSIZE, "", HELP_WINSIZE);
	//params->addFlag(ARG_POINTS, DEFAULT_POINTS, "", HELP_POINTS);
	//params->addFlag(ARG_BW, DEFAULT_BW, "", HELP_BW);
	params->addFlag(ARG_MAX_GAP, DEFAULT_MAX_GAP, "", HELP_MAX_GAP);
	params->addFlag(ARG_RESAMPLE, DEFAULT_RESAMPLE, "", HELP_RESAMPLE);
	params->addFlag(ARG_TPED, DEFAULT_TPED, "", HELP_TPED);
	params->addFlag(ARG_TFAM, DEFAULT_TFAM, "", HELP_TFAM);
	params->addFlag(ARG_RAW_LOD, DEFAULT_RAW_LOD, "", HELP_RAW_LOD);
	params->addListFlag(ARG_BOUND_SIZE, DEFAULT_BOUND_SIZE, "", HELP_BOUND_SIZE);
	params->addFlag(ARG_LOD_CUTOFF, DEFAULT_LOD_CUTOFF, "", HELP_LOD_CUTOFF);
	//params->addFlag(ARG_LOD_CUTOFF_FILE, DEFAULT_LOD_CUTOFF_FILE, "", HELP_LOD_CUTOFF_FILE);
	//params->addFlag(ARG_BOUND_SIZE_FILE, DEFAULT_BOUND_SIZE_FILE, "", HELP_BOUND_SIZE_FILE);
	params->addFlag(ARG_TPED_MISSING, DEFAULT_TPED_MISSING, "", HELP_TPED_MISSING);
	params->addFlag(ARG_FREQ_FILE, DEFAULT_FREQ_FILE, "", HELP_FREQ_FILE);
	params->addFlag(ARG_FREQ_ONLY, DEFAULT_FREQ_ONLY, "", HELP_FREQ_ONLY);
	params->addListFlag(ARG_WINSIZE_MULTI, DEFAULT_WINSIZE_MULTI, "", HELP_WINSIZE_MULTI);
	params->addFlag(ARG_POP_SPLIT, DEFAULT_POP_SPLIT , "", HELP_POP_SPLIT);
	params->addFlag(ARG_KDE_SUBSAMPLE, DEFAULT_KDE_SUBSAMPLE , "", HELP_KDE_SUBSAMPLE);
	params->addFlag(ARG_AUTO_WINSIZE, DEFAULT_AUTO_WINSIZE, "", HELP_AUTO_WINSIZE);
	params->addFlag(ARG_BUILD, DEFAULT_BUILD, "", HELP_BUILD);
	params->addFlag(ARG_CENTROMERE_FILE, DEFAULT_CENTROMERE_FILE, "", HELP_CENTROMERE_FILE);

	try { params->parseCommandLine(argc, argv); }
	catch (...) { throw 0; }

	return params;
}

bool checkBuild(string BUILD)
{
	if (BUILD.compare("hg18") != 0 &&
	        BUILD.compare("hg19") != 0 &&
	        BUILD.compare("hg38") != 0 &&
	        BUILD.compare(DEFAULT_BUILD) != 0) {
		cerr << "ERROR: Must choose hg18/hg19/hg38/none for build version.\n";
		return 1;
	}
	return 0;
}

bool checkMultiWinsizes(vector<int> &multiWinsizes, bool &WINSIZE_EXPLORE)
{
	if (multiWinsizes[0] != DEFAULT_WINSIZE_MULTI)
	{
		for (int i = 0; i < multiWinsizes.size(); i++)
		{
			if (multiWinsizes[i] <= 0)
			{
				cerr << "ERROR: SNP window sizes must be > 1.\n";
				return 1;
			}
		}
		WINSIZE_EXPLORE = true;
	}
	return 0;
}

bool checkAutoFreq(string freqfile, bool FREQ_ONLY, bool &AUTO_FREQ)
{
	if (freqfile.compare(DEFAULT_FREQ_FILE) != 0)
	{
		AUTO_FREQ = false;
		if (FREQ_ONLY)
		{
			cerr << "ERROR: Specifying both " << ARG_FREQ_ONLY << " and " << ARG_FREQ_FILE << " accomplishes nothing.\n";
			return 1;
		}
	}
	return 0;
}

bool checkAutoWinsize(bool WINSIZE_EXPLORE, bool AUTO_WINSIZE)
{
	//Check if both AUTO_WINSIZE and WINSIZE_EXPLORE are set
	//If so, exit with error.
	if (WINSIZE_EXPLORE && AUTO_WINSIZE)
	{
		cerr << "ERROR: Must set only one of " << ARG_WINSIZE_MULTI << " and " << ARG_AUTO_WINSIZE << ".\n";
		return 1;
	}
	return 0;
}

bool checkAutoCutoff(bool LOD_CUTOFF, bool &AUTO_CUTOFF)
{
	if (LOD_CUTOFF != DEFAULT_LOD_CUTOFF) {
		AUTO_CUTOFF = false;
	}
	return 0;
}

bool checkBoundSizes(vector<double> &boundSizes, bool &AUTO_BOUNDS)
{
	if (boundSizes[0] != DEFAULT_BOUND_SIZE && boundSizes.size() != 2) {
		cerr << "ERROR: Must provide two bounds.\n";
		return 1;
	}
	else if (boundSizes.size() == 2)
	{
		double tmp;
		AUTO_BOUNDS = false;
		if (boundSizes[0] <= 0 || boundSizes[1] <= 0)
		{
			cerr << "ERROR: User provided size boundaries must be positive.\n";
			return 1;
		}
		else if (boundSizes[0] > boundSizes[1])
		{
			tmp = boundSizes[0];
			boundSizes[0] = boundSizes[1];
			boundSizes[1] = tmp;
		}
		else if (boundSizes[0] == boundSizes[1])
		{
			cerr << "ERROR: Size boundaries must be different.\n";
			return 1;
		}
	}
	return 0;
}

bool checkRequiredFiles(string tpedfile, string tfamfile)
{
	if (tpedfile.compare(DEFAULT_TPED) == 0 || tfamfile.compare(DEFAULT_TFAM) == 0)
	{
		cerr << "ERROR: Must provide both a tped and tfam file.\n";
		return 1;
	}
	return 0;
}

bool checkThreads(int numThreads)
{
	if (numThreads <= 0)
	{
		cerr << "ERROR: Number of threads must be > 0.\n";
		return 1;
	}
	return 0;
}

bool checkError(double error)
{
	if (error <= 0 || error >= 1)
    {
        cerr << "ERROR: Genotype error rate must be > 0 and < 1.\n";
        return 1;
    }
    return 0;
}

bool checkWinsize(int winsize)
{
	if (winsize <= 1)
    {
        cerr << "ERROR: SNP window size must be > 1.\n";
        return 1;
    }
    return 0;
}

bool checkMaxGap(int MAX_GAP)
{
	if(MAX_GAP < 0)
	{
		cerr << "ERROR: max gap must be > 0.\n";
		return 1;
	}
	else if(MAX_GAP < 10000)
	{
		cerr << "WARNING: max gap set very low: " << MAX_GAP << endl;
	}
	return 0;
}