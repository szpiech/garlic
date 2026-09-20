#ifndef __GARLIC_CLI_H__
#define __GARLIC_CLI_H__

#include <iostream>
#include <string>
#include "param_t.h"
#include "garlic-errlog.h"

using namespace std;

extern const string VERSION;

extern const string ARG_OUTFILE;
extern const string DEFAULT_OUTFILE;
extern const string HELP_OUTFILE;

extern const string ARG_TPED;
extern const string DEFAULT_TPED;
extern const string HELP_TPED;

extern const string ARG_TFAM;
extern const string DEFAULT_TFAM;
extern const string HELP_TFAM;

extern const string ARG_FEATURES;
extern const string DEFAULT_FEATURES;
extern const string HELP_FEATURES;

extern const string ARG_ROHFILE;
extern const string DEFAULT_ROHFILE;
extern const string HELP_ROHFILE;

extern const string ARG_TPED_MISSING;
extern const char DEFAULT_TPED_MISSING;
extern const string HELP_TPED_MISSING;

param_t *getCLI(int argc, char *argv[]);
bool checkRequiredFiles(string tpedfile, string tfamfile, string featurefile, string rohfile);

#endif
