#ifndef __GARLIC_CENTROMERES_H__
#define __GARLIC_CENTROMERES_H__
#include <map>
#include <string>
#include "gzstream.h"
#include <iostream>
#include "garlic-errlog.h"

using namespace std;

class centromere {

public:

	centromere(string arg, string file, string defaultFileName);
	centromere();

	int centromereStart(string chr);
	int centromereEnd(string chr);

	void readCustomCentromeres(string filename);
	void makeHG18();
	void makeHG19();
	void makeHG38();

private:

	map <string, int> gapStart;
	map <string, int> gapEnd;

	int countFields(const string &str);
};

#endif