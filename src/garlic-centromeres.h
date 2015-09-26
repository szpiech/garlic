#ifndef __GARLIC_DATA_H__
#define __GARLIC_DATA_H__
#include <map>
#include <string>
#include "garlic-data.h"
#include "gzstream.h"
#include <iostream>

class centromere{

public:

	centromere(string arg);

	int centromereStart(string chr);
	int centromereEnd(string chr);

private:
	
	map <string, int_pair_t> gap;

	void readCustomCentromeres(string filename);
	void makeHG18();
	void makeHG19();
//	void makeHG38();

	int countFields(const string &str);
};

#endif