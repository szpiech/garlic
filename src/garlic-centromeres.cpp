#include "garlic-centromeres.h"

using namespace std;

centromere::centromere(string arg) {
	if (arg.compare("hg18") == 0) {
		makeHG18();
	}
	else if (arg.compare("hg19") == 0) {
		makeHG19()
	}
	else if (arg.compare("hg38") == 0) {
//		makeHG38();
	}
	else{
		readCustomCentromeres(arg);
	}
}
int centromere::centromereStart(string chr) {
	if (gap.count(chr) == 0) return 0;
	return gap[chr].first;
}
int centromere::centromereEnd(string chr) {
	if (gap.count(chr) == 0) return 0;
	return gap[chr].second;
}

/*
	<chr name> <start> <end>
*/
void centromere::readCustomCentromeres(string filename){
	igzstream fin;
	fin.open(filename.c_str());
	if(fin.fail()){
		cerr << "ERROR: Could not open " << filename << " for reading.\n";
		throw -1;
	}

	string line;
	int numChr = 0;
	int curCols = 0;
	while(getline(fin,line)){
		numChr++;
		curCols = countFields(line);
		if(curCols != 3){
			cerr << "ERORR: Custom centromere file requires three columns.  Found " << curCols << ".\n";
		}
	}
	fin.close();

	string chrname;
	fin.open(filename.c_str());
	for(int row = 0; row < numChr; row++){
		fin >> chrname;
		fin >> gap[chrname].first >> gap[chrname].second;
	}
	cerr << "Loaded custom centromere limits for " << numChr << " chromosomes.\n";
	return;
}
void centromere::makeHG18() {
	//UCSC build hg18 centromere boundaries
	gap["chr1"].first = 121236957;
	gap["chr1"].second = 123476957;

	gap["chr2"].first = 91689898;
	gap["chr2"].second = 94689898;

	gap["chr3"].first = 90587544;
	gap["chr3"].second = 93487544;

	gap["chr4"].first = 49354874;
	gap["chr4"].second = 52354874;

	gap["chr5"].first = 46441398;
	gap["chr5"].second = 49441398;

	gap["chr6"].first = 58938125;
	gap["chr6"].second = 61938125;

	gap["chr7"].first = 58058273;
	gap["chr7"].second = 61058273;

	gap["chr8"].first = 43958052;
	gap["chr8"].second = 46958052;

	gap["chr9"].first = 47107499;
	gap["chr9"].second = 50107499;

	gap["chr10"].first = 39244941;
	gap["chr10"].second = 41624941;

	gap["chr11"].first = 51450781;
	gap["chr11"].second = 54450781;

	gap["chr12"].first = 34747961;
	gap["chr12"].second = 36142961;

	gap["chr13"].first = 16000000;
	gap["chr13"].second = 17868000;

	gap["chr14"].first = 15070000;
	gap["chr14"].second = 18070000;

	gap["chr15"].first = 15260000;
	gap["chr15"].second = 18260000;

	gap["chr16"].first = 35143302;
	gap["chr16"].second = 36943302;

	gap["chr17"].first = 22187133;
	gap["chr17"].second = 22287133;

	gap["chr18"].first = 15400898;
	gap["chr18"].second = 16764896;

	gap["chr19"].first = 26923622;
	gap["chr19"].second = 29923622;

	gap["chr20"].first = 26267569;
	gap["chr20"].second = 28033230;

	gap["chr21"].first = 10260000;
	gap["chr21"].second = 13260000;

	gap["chr22"].first = 11330000;
	gap["chr22"].second = 14330000;

	gap["chrX"].first = 58598737;
	gap["chrX"].second = 61598737;

}

void centromere::makeHG19() {
	//UCSC build hg19 centromere boundaries
	gap["chr1"].first = 121535434;
	gap["chr1"].second = 124535434;

	gap["chr2"].first = 92326171;
	gap["chr2"].second = 95326171;

	gap["chr3"].first = 90504854;
	gap["chr3"].second = 93504854;

	gap["chr4"].first = 49660117;
	gap["chr4"].second = 52660117;

	gap["chr5"].first = 46405641;
	gap["chr5"].second = 49405641;

	gap["chr6"].first = 58830166;
	gap["chr6"].second = 61830166;

	gap["chr7"].first = 58054331;
	gap["chr7"].second = 61054331;

	gap["chr8"].first = 43838887;
	gap["chr8"].second = 46838887;

	gap["chr9"].first = 47367679;
	gap["chr9"].second = 50367679;

	gap["chr10"].first = 39254935;
	gap["chr10"].second = 42254935;

	gap["chr11"].first = 51644205;
	gap["chr11"].second = 54644205;

	gap["chr12"].first = 34856694;
	gap["chr12"].second = 37856694;

	gap["chr13"].first = 16000000;
	gap["chr13"].second = 19000000;

	gap["chr14"].first = 16000000;
	gap["chr14"].second = 19000000;

	gap["chr15"].first = 17000000;
	gap["chr15"].second = 20000000;

	gap["chr16"].first = 35335801;
	gap["chr16"].second = 38335801;

	gap["chr17"].first = 22263006;
	gap["chr17"].second = 25263006;

	gap["chr18"].first = 15460898;
	gap["chr18"].second = 18460898;

	gap["chr19"].first = 24681782;
	gap["chr19"].second = 27681782;

	gap["chr20"].first = 26369569;
	gap["chr20"].second = 29369569;

	gap["chr21"].first = 11288129;
	gap["chr21"].second = 14288129;

	gap["chr22"].first = 13000000;
	gap["chr22"].second = 16000000;

	gap["chrX"].first = 58632012;
	gap["chrX"].second = 61632012;
}


int centromere::countFields(const string &str)
{
    string::const_iterator it;
    int result;
    int numFields = 0;
    int seenChar = 0;
    for (it = str.begin() ; it < str.end(); it++)
    {
        result = isspace(*it);
        if (result == 0 && seenChar == 0)
        {
            numFields++;
            seenChar = 1;
        }
        else if (result != 0)
        {
            seenChar = 0;
        }
    }
    return numFields;
}