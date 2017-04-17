#include "garlic-centromeres.h"

centromere::centromere(string arg, string file, string defaultFileName) {
	if (arg.compare("hg18") == 0) {
		makeHG18();
		makeWarning();
	}
	else if (arg.compare("hg19") == 0) {
		makeHG19();
		makeWarning();
	}
	else if (arg.compare("hg38") == 0) {
		makeHG38();
		makeWarning();
	}
	else if (file.compare(defaultFileName) != 0) {
		readCustomCentromeres(file);
	}
	else if (arg.compare("none") == 0) {
		gapStart.clear();
		gapEnd.clear();
		makeWarning();
	}
	return;
}

centromere::centromere() {
	gapStart.clear();
	gapEnd.clear();
	return;
}

int centromere::centromereStart(string chr) {
	if (gapStart.count(chr) == 0) {
		if (chrWarning[chr] == 0) {
			LOG.err("WARNING: No centromere start information for chr:", chr);
			LOG.err("WARNING: If you provided custom centromeres check that chromosome names match between data files.");
			//LOG.log("WARNING: No centromere start information for chr:", chr);
			//LOG.log("WARNING: If you provided custom centromeres check that chromosome names match between data files.");
			chrWarning[chr]++;
		}
		return 0;
	}
	return gapStart[chr];
}

int centromere::centromereEnd(string chr) {
	if (gapEnd.count(chr) == 0) {
		if (chrWarning[chr] == 0) {
			LOG.err("WARNING: No centromere end information for chr:", chr);
			LOG.err("WARNING: If you provided custom centromeres check that chromosome names match between data files.");
			//LOG.log("WARNING: No centromere end information for chr:", chr);
			//LOG.log("WARNING: If you provided custom centromeres check that chromosome names match between data files.");
			chrWarning[chr]++;
		}
		return 0;
	}
	return gapEnd[chr];
}

/*
	<chr name> <start> <end>
*/
void centromere::readCustomCentromeres(string filename) {
	igzstream fin;
	fin.open(filename.c_str());
	if (fin.fail()) {
		//cerr << "ERROR: Could not open " << filename << " for reading.\n";
		LOG.err("ERROR: Could not open", filename);
		throw - 1;
	}

	string line;
	int numChr = 0;
	int curCols = 0;
	while (getline(fin, line)) {
		numChr++;
		curCols = countFields(line);
		if (curCols != 3) {
			//cerr << "ERROR: Custom centromere file requires three columns.  Found " << curCols << ".\n";
			LOG.err("ERROR: Custom centromere file requires three columns.  Found", curCols);
		}
	}
	fin.close();
	fin.clear();
	string chrname;
	int start;
	int end;
	fin.open(filename.c_str());
	for (int row = 0; row < numChr; row++) {
		fin >> chrname;
		fin >> start >> end;
		gapStart[chrname] = start;
		gapEnd[chrname] = end;
		chrWarning[chrname] = 0;
	}
	cerr << "Loaded custom centromere limits for " << numChr << " chromosomes.\n";
	fin.close();
	return;
}

void centromere::makeWarning() {

	chrWarning.clear();

	chrWarning["chr1"] = 0;
	chrWarning["1"] = 0;

	chrWarning["chr2"] = 0;
	chrWarning["2"] = 0;

	chrWarning["chr3"] = 0;
	chrWarning["3"] = 0;

	chrWarning["chr4"] = 0;
	chrWarning["4"] = 0;

	chrWarning["chr5"] = 0;
	chrWarning["5"] = 0;

	chrWarning["chr6"] = 0;
	chrWarning["6"] = 0;

	chrWarning["chr7"] = 0;
	chrWarning["7"] = 0;

	chrWarning["chr8"] = 0;
	chrWarning["8"] = 0;

	chrWarning["chr9"] = 0;
	chrWarning["9"] = 0;

	chrWarning["chr10"] = 0;
	chrWarning["10"] = 0;

	chrWarning["chr11"] = 0;
	chrWarning["11"] = 0;

	chrWarning["chr12"] = 0;
	chrWarning["12"] = 0;

	chrWarning["chr13"] = 0;
	chrWarning["13"] = 0;

	chrWarning["chr14"] = 0;
	chrWarning["14"] = 0;

	chrWarning["chr15"] = 0;
	chrWarning["15"] = 0;

	chrWarning["chr16"] = 0;
	chrWarning["16"] = 0;

	chrWarning["chr17"] = 0;
	chrWarning["17"] = 0;

	chrWarning["chr18"] = 0;
	chrWarning["18"] = 0;

	chrWarning["chr19"] = 0;
	chrWarning["19"] = 0;

	chrWarning["chr20"] = 0;
	chrWarning["20"] = 0;

	chrWarning["chr21"] = 0;
	chrWarning["21"] = 0;

	chrWarning["chr22"] = 0;
	chrWarning["22"] = 0;

	chrWarning["chrX"] = 0;
	chrWarning["X"] = 0;
}


void centromere::makeHG18() {
	gapStart.clear();
	gapEnd.clear();

	//UCSC build hg18 centromere boundaries
	gapStart["chr1"] = 121236957;
	gapEnd["chr1"] = 123476957;

	gapStart["chr2"] = 91689898;
	gapEnd["chr2"] = 94689898;

	gapStart["chr3"] = 90587544;
	gapEnd["chr3"] = 93487544;

	gapStart["chr4"] = 49354874;
	gapEnd["chr4"] = 52354874;

	gapStart["chr5"] = 46441398;
	gapEnd["chr5"] = 49441398;

	gapStart["chr6"] = 58938125;
	gapEnd["chr6"] = 61938125;

	gapStart["chr7"] = 58058273;
	gapEnd["chr7"] = 61058273;

	gapStart["chr8"] = 43958052;
	gapEnd["chr8"] = 46958052;

	gapStart["chr9"] = 47107499;
	gapEnd["chr9"] = 50107499;

	gapStart["chr10"] = 39244941;
	gapEnd["chr10"] = 41624941;

	gapStart["chr11"] = 51450781;
	gapEnd["chr11"] = 54450781;

	gapStart["chr12"] = 34747961;
	gapEnd["chr12"] = 36142961;

	gapStart["chr13"] = 16000000;
	gapEnd["chr13"] = 17868000;

	gapStart["chr14"] = 15070000;
	gapEnd["chr14"] = 18070000;

	gapStart["chr15"] = 15260000;
	gapEnd["chr15"] = 18260000;

	gapStart["chr16"] = 35143302;
	gapEnd["chr16"] = 36943302;

	gapStart["chr17"] = 22187133;
	gapEnd["chr17"] = 22287133;

	gapStart["chr18"] = 15400898;
	gapEnd["chr18"] = 16764896;

	gapStart["chr19"] = 26923622;
	gapEnd["chr19"] = 29923622;

	gapStart["chr20"] = 26267569;
	gapEnd["chr20"] = 28033230;

	gapStart["chr21"] = 10260000;
	gapEnd["chr21"] = 13260000;

	gapStart["chr22"] = 11330000;
	gapEnd["chr22"] = 14330000;

	gapStart["chrX"] = 58598737;
	gapEnd["chrX"] = 61598737;

	gapStart["1"] = 121236957;
	gapEnd["1"] = 123476957;

	gapStart["2"] = 91689898;
	gapEnd["2"] = 94689898;

	gapStart["3"] = 90587544;
	gapEnd["3"] = 93487544;

	gapStart["4"] = 49354874;
	gapEnd["4"] = 52354874;

	gapStart["5"] = 46441398;
	gapEnd["5"] = 49441398;

	gapStart["6"] = 58938125;
	gapEnd["6"] = 61938125;

	gapStart["7"] = 58058273;
	gapEnd["7"] = 61058273;

	gapStart["8"] = 43958052;
	gapEnd["8"] = 46958052;

	gapStart["9"] = 47107499;
	gapEnd["9"] = 50107499;

	gapStart["10"] = 39244941;
	gapEnd["10"] = 41624941;

	gapStart["11"] = 51450781;
	gapEnd["11"] = 54450781;

	gapStart["12"] = 34747961;
	gapEnd["12"] = 36142961;

	gapStart["13"] = 16000000;
	gapEnd["13"] = 17868000;

	gapStart["14"] = 15070000;
	gapEnd["14"] = 18070000;

	gapStart["15"] = 15260000;
	gapEnd["15"] = 18260000;

	gapStart["16"] = 35143302;
	gapEnd["16"] = 36943302;

	gapStart["17"] = 22187133;
	gapEnd["17"] = 22287133;

	gapStart["18"] = 15400898;
	gapEnd["18"] = 16764896;

	gapStart["19"] = 26923622;
	gapEnd["19"] = 29923622;

	gapStart["20"] = 26267569;
	gapEnd["20"] = 28033230;

	gapStart["21"] = 10260000;
	gapEnd["21"] = 13260000;

	gapStart["22"] = 11330000;
	gapEnd["22"] = 14330000;

	gapStart["X"] = 58598737;
	gapEnd["X"] = 61598737;
}

void centromere::makeHG19() {
	gapStart.clear();
	gapEnd.clear();

	//UCSC build hg19 centromere boundaries
	gapStart["chr1"] = 121535434;
	gapEnd["chr1"] = 124535434;

	gapStart["chr2"] = 92326171;
	gapEnd["chr2"] = 95326171;

	gapStart["chr3"] = 90504854;
	gapEnd["chr3"] = 93504854;

	gapStart["chr4"] = 49660117;
	gapEnd["chr4"] = 52660117;

	gapStart["chr5"] = 46405641;
	gapEnd["chr5"] = 49405641;

	gapStart["chr6"] = 58830166;
	gapEnd["chr6"] = 61830166;

	gapStart["chr7"] = 58054331;
	gapEnd["chr7"] = 61054331;

	gapStart["chr8"] = 43838887;
	gapEnd["chr8"] = 46838887;

	gapStart["chr9"] = 47367679;
	gapEnd["chr9"] = 50367679;

	gapStart["chr10"] = 39254935;
	gapEnd["chr10"] = 42254935;

	gapStart["chr11"] = 51644205;
	gapEnd["chr11"] = 54644205;

	gapStart["chr12"] = 34856694;
	gapEnd["chr12"] = 37856694;

	gapStart["chr13"] = 16000000;
	gapEnd["chr13"] = 19000000;

	gapStart["chr14"] = 16000000;
	gapEnd["chr14"] = 19000000;

	gapStart["chr15"] = 17000000;
	gapEnd["chr15"] = 20000000;

	gapStart["chr16"] = 35335801;
	gapEnd["chr16"] = 38335801;

	gapStart["chr17"] = 22263006;
	gapEnd["chr17"] = 25263006;

	gapStart["chr18"] = 15460898;
	gapEnd["chr18"] = 18460898;

	gapStart["chr19"] = 24681782;
	gapEnd["chr19"] = 27681782;

	gapStart["chr20"] = 26369569;
	gapEnd["chr20"] = 29369569;

	gapStart["chr21"] = 11288129;
	gapEnd["chr21"] = 14288129;

	gapStart["chr22"] = 13000000;
	gapEnd["chr22"] = 16000000;

	gapStart["chrX"] = 58632012;
	gapEnd["chrX"] = 61632012;


	gapStart["1"] = 121535434;
	gapEnd["1"] = 124535434;

	gapStart["2"] = 92326171;
	gapEnd["2"] = 95326171;

	gapStart["3"] = 90504854;
	gapEnd["3"] = 93504854;

	gapStart["4"] = 49660117;
	gapEnd["4"] = 52660117;

	gapStart["5"] = 46405641;
	gapEnd["5"] = 49405641;

	gapStart["6"] = 58830166;
	gapEnd["6"] = 61830166;

	gapStart["7"] = 58054331;
	gapEnd["7"] = 61054331;

	gapStart["8"] = 43838887;
	gapEnd["8"] = 46838887;

	gapStart["9"] = 47367679;
	gapEnd["9"] = 50367679;

	gapStart["10"] = 39254935;
	gapEnd["10"] = 42254935;

	gapStart["11"] = 51644205;
	gapEnd["11"] = 54644205;

	gapStart["12"] = 34856694;
	gapEnd["12"] = 37856694;

	gapStart["13"] = 16000000;
	gapEnd["13"] = 19000000;

	gapStart["14"] = 16000000;
	gapEnd["14"] = 19000000;

	gapStart["15"] = 17000000;
	gapEnd["15"] = 20000000;

	gapStart["16"] = 35335801;
	gapEnd["16"] = 38335801;

	gapStart["17"] = 22263006;
	gapEnd["17"] = 25263006;

	gapStart["18"] = 15460898;
	gapEnd["18"] = 18460898;

	gapStart["19"] = 24681782;
	gapEnd["19"] = 27681782;

	gapStart["20"] = 26369569;
	gapEnd["20"] = 29369569;

	gapStart["21"] = 11288129;
	gapEnd["21"] = 14288129;

	gapStart["22"] = 13000000;
	gapEnd["22"] = 16000000;

	gapStart["X"] = 58632012;
	gapEnd["X"] = 61632012;
}

void centromere::makeHG38() {
	gapStart.clear();
	gapEnd.clear();

	//UCSC build hg38 centromere boundaries

	gapStart["chr1"] = 122026459;
	gapEnd["chr1"] = 124932724;

	gapStart["chr2"] = 92188145;
	gapEnd["chr2"] = 94090557;

	gapStart["chr3"] = 90772458;
	gapEnd["chr3"] = 93655574;

	gapStart["chr4"] = 49712061;
	gapEnd["chr4"] = 51743951;

	gapStart["chr5"] = 46485900;
	gapEnd["chr5"] = 50059807;

	gapStart["chr6"] = 58553888;
	gapEnd["chr6"] = 59829934;

	gapStart["chr7"] = 58169653;
	gapEnd["chr7"] = 61528020;

	gapStart["chr8"] = 44033744;
	gapEnd["chr8"] = 45877265;

	gapStart["chr9"] = 43389635;
	gapEnd["chr9"] = 45518558;

	gapStart["chr10"] = 39686682;
	gapEnd["chr10"] = 41593521;

	gapStart["chr11"] = 51078348;
	gapEnd["chr11"] = 54425074;

	gapStart["chr12"] = 34769407;
	gapEnd["chr12"] = 37185252;

	gapStart["chr13"] = 16000000;
	gapEnd["chr13"] = 18051248;

	gapStart["chr14"] = 16000000;
	gapEnd["chr14"] = 18173523;

	gapStart["chr15"] = 17083673;
	gapEnd["chr15"] = 19725254;

	gapStart["chr16"] = 36311158;
	gapEnd["chr16"] = 38265669;

	gapStart["chr17"] = 22813679;
	gapEnd["chr17"] = 26616164;

	gapStart["chr18"] = 15460899;
	gapEnd["chr18"] = 20861206;

	gapStart["chr19"] = 24498980;
	gapEnd["chr19"] = 27190874;

	gapStart["chr20"] = 26436232;
	gapEnd["chr20"] = 30038348;

	gapStart["chr21"] = 10864560;
	gapEnd["chr21"] = 12915808;

	gapStart["chr22"] = 12954788;
	gapEnd["chr22"] = 15054318;

	gapStart["chrX"] = 58605579;
	gapEnd["chrX"] = 62412542;


	gapStart["1"] = 122026459;
	gapEnd["1"] = 124932724;

	gapStart["2"] = 92188145;
	gapEnd["2"] = 94090557;

	gapStart["3"] = 90772458;
	gapEnd["3"] = 93655574;

	gapStart["4"] = 49712061;
	gapEnd["4"] = 51743951;

	gapStart["5"] = 46485900;
	gapEnd["5"] = 50059807;

	gapStart["6"] = 58553888;
	gapEnd["6"] = 59829934;

	gapStart["7"] = 58169653;
	gapEnd["7"] = 61528020;

	gapStart["8"] = 44033744;
	gapEnd["8"] = 45877265;

	gapStart["9"] = 43389635;
	gapEnd["9"] = 45518558;

	gapStart["10"] = 39686682;
	gapEnd["10"] = 41593521;

	gapStart["11"] = 51078348;
	gapEnd["11"] = 54425074;

	gapStart["12"] = 34769407;
	gapEnd["12"] = 37185252;

	gapStart["13"] = 16000000;
	gapEnd["13"] = 18051248;

	gapStart["14"] = 16000000;
	gapEnd["14"] = 18173523;

	gapStart["15"] = 17083673;
	gapEnd["15"] = 19725254;

	gapStart["16"] = 36311158;
	gapEnd["16"] = 38265669;

	gapStart["17"] = 22813679;
	gapEnd["17"] = 26616164;

	gapStart["18"] = 15460899;
	gapEnd["18"] = 20861206;

	gapStart["19"] = 24498980;
	gapEnd["19"] = 27190874;

	gapStart["20"] = 26436232;
	gapEnd["20"] = 30038348;

	gapStart["21"] = 10864560;
	gapEnd["21"] = 12915808;

	gapStart["22"] = 12954788;
	gapEnd["22"] = 15054318;

	gapStart["X"] = 58605579;
	gapEnd["X"] = 62412542;

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