#ifndef __GARLIC_CENTROMERES_H__
#define __GARLIC_CENTROMERES_H__
#include <map>
#include "garlic-pos.h"
#include <string>
#include "gzstream.h"
#include <iostream>
#include "garlic-errlog.h"

using namespace std;

//Defined in garlic-data.cpp.  Declared again here rather than included from
//garlic-data.h, which includes THIS header: chromosome-name canonicalisation
//has to be available on both sides of that dependency, and a repeated
//declaration is the smaller price.
string canonChrKey(const string &name);

class centromere {

public:

	centromere(string arg, string file, string defaultFileName);
	centromere();

	pos_t centromereStart(string chr);
	pos_t centromereEnd(string chr);

	//--no-centromere: there is no assembled gap to skip.  Every lookup
	//returns 0 and the missing-chromosome warning is suppressed, because the
	//absence is the user's explicit choice rather than a name mismatch.
	void suppressMissingWarnings();

	void readCustomCentromeres(string filename);
	void makeHG18();
	void makeHG19();
	void makeHG38();
	void makeT2TCHM13();
	void makeWarning();

private:

	map <string, pos_t> gapStart;
	map <string, pos_t> gapEnd;
	bool quietMissing;

	int countFields(const string &str);
	map <string, int> chrWarning;
	string checkChrName(string chr);

	//Rekeys gapStart/gapEnd by canonChrKey once the table is built, so a
	//lookup matches whatever spelling the data file uses.
	void canonicaliseKeys();
};

#endif