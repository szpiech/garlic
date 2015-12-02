#include "garlic-cli.h"

param_t *getCLI(int argc, char *argv[]){
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