#include "rohscan-data.h"

FreqData* calcFreqData(HapData* hapData)
{
  FreqData* freqData = initFreqData(hapData->nloci);
  double total, freq;

  for(int locus = 0; locus < hapData->nloci; locus++)
    {
      total = 0;
      freq = 0;
      for(int hap = 0; hap < hapData->nhaps; hap++)
	{
	  if(hapData->data[hap][locus] != -9)
	    {
	      freq += hapData->data[hap][locus];
	      total++;
	    }
	}
      freqData->freq[locus] = freq/total;
    }
  return freqData;
}

vector< vector< FreqData* >* >* calcFreqData(vector< vector< HapData* >* >* hapDataByPopByChr)
{
  vector< vector< FreqData* >* >* freqDataByPopByChr = new vector< vector< FreqData* >* >;
  
  for(int pop = 0; pop < hapDataByPopByChr->size(); pop++)
    {
      vector< FreqData* >* freqDataByChr = new vector< FreqData* >;
      for (int chr = 0; chr < hapDataByPopByChr->at(pop)->size(); chr++)
	{
	  FreqData* data = calcFreqData(hapDataByPopByChr->at(pop)->at(chr));
	  freqDataByChr->push_back(data);
	}
      freqDataByPopByChr->push_back(freqDataByChr);
    }
  return freqDataByPopByChr;
}


//allocates the arrays and populates them with MISSING
FreqData* initFreqData(int nloci)
{
  if(nloci < 1)
    {
      cerr << "ERROR: number of loci ("<< nloci <<") must be positive.\n";
      throw 0;
    }

  FreqData *data = new FreqData;
  data->nloci = nloci;
  data->freq = new double[nloci];

  for (int locus = 0; locus < nloci; locus++)
    {
      data->freq[locus] = MISSING;
    }

  return data;
}

void releaseFreqData(FreqData *data)
{
  if(data == NULL) return;
  data->nloci = -9;
  delete [] data->freq;
  delete data;
  data = NULL;
  return;
}

void releaseFreqData(vector< vector< FreqData* >* >* freqDataByPopByChr)
{
  for(int pop = 0; pop < freqDataByPopByChr->size(); pop++)
    {
      for(int chr = 0; chr < freqDataByPopByChr->at(pop)->size(); chr++)
	{
	  releaseFreqData(freqDataByPopByChr->at(pop)->at(chr));
	}
      freqDataByPopByChr->at(pop)->clear();
      delete freqDataByPopByChr->at(pop);
    }
  freqDataByPopByChr->clear();
  delete freqDataByPopByChr;
  return;
}

//allocates the arrays and populates them with MISSING or "--" depending on type
MapData* initMapData(int nloci)
{
  if(nloci < 1)
    {
      cerr << "ERROR: number of loci ("<< nloci <<") must be positive.\n";
      throw 0;
    }

  MapData *data = new MapData;
  data->nloci = nloci;
  data->locusName = new string[nloci];
  data->physicalPos = new int[nloci];
  data->geneticPos = new double[nloci];

  for (int locus = 0; locus < nloci; locus++)
    {
      data->locusName[locus] = "--";
      data->physicalPos[locus] = MISSING;
      data->geneticPos[locus] = MISSING;
    }

  return data;
}

void releaseMapData(MapData *data)
{
  if(data == NULL) return;
  data->nloci = -9;
  delete [] data->locusName;
  delete [] data->physicalPos;
  delete [] data->geneticPos;
  delete data;
  data = NULL;
  return;
}

void releaseMapData(vector< MapData* >* mapDataByChr)
{
  for (int i = 0; i < mapDataByChr->size(); i++)
    {
      releaseMapData(mapDataByChr->at(i));
    }
  mapDataByChr->clear();
  delete mapDataByChr;
  return;
}

void releaseIndData(vector< IndData* >* indDataByPop)
{
  for (int i = 0; i < indDataByPop->size(); i++)
    {
      releaseIndData(indDataByPop->at(i));
    }
  indDataByPop->clear();
  delete indDataByPop;
  return;
}

void releaseIndData(IndData *data)
{
  delete [] data->indID;
  delete data;
  return;
}

vector< vector< HapData* >* >* readHapData(string filename, 
						       int expectedLoci, 
						       int expectedInd,
						       vector< int_pair_t >* chrCoordList,
						       vector< int_pair_t >* indCoordList)
{
  int expectedHaps = 2*expectedInd;
  ifstream fin;
  cerr << "Checking " << filename << "...\n";
  fin.open(filename.c_str());
  
  if(fin.fail())
    {
      cerr << "ERROR: Failed to open " << filename << " for reading.\n";
      throw 0;
    }
  
  int fileStart = fin.tellg();
  string line;
  int nhaps = 0;
  int nloci = -1;
  while(getline(fin,line))
    {
      nhaps++;
      nloci = countFields(line);
      //cout << "nloci: " << current_nloci << endl;
      if(nloci != expectedLoci)
	{
	  cerr << "ERROR: line " << nhaps << " of " << filename << " has " << nloci 
	       << ", but expected " << expectedLoci << ".\n";
	  throw 0;
	}
    }
  if(nhaps != expectedHaps)
    {
      cerr << "ERROR: " << filename << " has " << nhaps 
	   << " haplotypes, but expected " << expectedHaps << ".\n";
      throw 0;
    }
  
  fin.clear(); // clear error flags
  fin.seekg(fileStart);
  
  cerr << "Loading " << filename << "...\n";

  vector< vector< HapData* >* >* hapDataByPopByChr = new vector< vector< HapData* >* >;

  for(int pop = 0; pop < indCoordList->size(); pop++)
    {
      int totalHaps = 2*(indCoordList->at(pop).second - indCoordList->at(pop).first + 1);
      vector< HapData* >* hapDataByChr = new vector< HapData* >;
      for (int chr = 0; chr < chrCoordList->size(); chr++)
	{
	  int totalLoci = chrCoordList->at(chr).second - chrCoordList->at(chr).first + 1;
	  HapData* data = initHapData(totalHaps,totalLoci);
	  hapDataByChr->push_back(data);
	}
      hapDataByPopByChr->push_back(hapDataByChr);
    }

  //For each population
  for(int pop = 0; pop < indCoordList->size(); pop++)
    {
      int totalHaps = 2*(indCoordList->at(pop).second - indCoordList->at(pop).first + 1);
      //For each haplotype in the population
      for(int hap = 0; hap < totalHaps; hap++)
	{
	  //For each chromosome
	  for (int chr = 0; chr < chrCoordList->size(); chr++)
	    {
	      int totalLoci = chrCoordList->at(chr).second - chrCoordList->at(chr).first + 1;
	      //For each locus on the chromosome
	      for (int locus = 0; locus < totalLoci; locus++)
		{
		  short allele;
		  fin >> allele;
		  if(allele != 0 && allele != 1 && allele != -9)
		    {
		      string hapPost, popPost, chrPost, locPost;
		      hapPost = getPost(hap+1);
		      popPost = getPost(pop+1);
		      chrPost = getPost(chr+1);
		      locPost = getPost(locus+1);

		      cerr << "ERROR: The " << hap+1 << hapPost << " haplotype in the " 
			   << pop+1 << popPost << " population at the " 
			   << locus+1 << locPost << " locus on the "
			   << chr+1 << chrPost << " chromosome has an illegal value.  Must be 0/1/-9.\n";
		      
		      throw 0;
		    }
		  hapDataByPopByChr->at(pop)->at(chr)->data[hap][locus] = allele;
		}
	    }
	}
      
    }
  fin.close();

  return hapDataByPopByChr;
}

void releaseHapData(vector< vector< HapData* >* >* hapDataByPopByChr)
{
  for(int pop = 0; pop < hapDataByPopByChr->size(); pop++)
    {
      for (int chr = 0; chr < hapDataByPopByChr->at(pop)->size(); chr++)
	{
	  releaseHapData(hapDataByPopByChr->at(pop)->at(chr));
	}
      hapDataByPopByChr->at(pop)->clear();
      delete hapDataByPopByChr->at(pop);
    }
  hapDataByPopByChr->clear();
  delete hapDataByPopByChr;
}

string getPost(int num)
{
  string post;
  if(num == 1) post = "st";
  else if(num == 2) post = "nd";
  else if(num == 3) post = "rd";
  else post = "th";
  return post;
}

WinData* initWinData(unsigned int nind, unsigned int nloci)
{
  if(nind < 1 || nloci < 1)
    {
      cerr << "ERROR: Can't allocate WinData object.  Number of individuals ("<< nind <<") and number of loci ("<< nloci <<") must be positive.\n";
      throw 0;
    }
  
  WinData* data = new WinData;
  data->nind = nind;
  data->nloci = nloci;
  //data->nmiss = 0;
  
  data->data = new double*[nind];
  for (unsigned int i = 0; i < nind; i++)
    {
      data->data[i] = new double[nloci];
      for (unsigned int j = 0; j < nloci; j++)
	{
	  data->data[i][j] = MISSING;
	}
    }
  
  return data;
}

void releaseWinData(WinData *data)
{
  if(data == NULL) return;
  for (int i = 0; i < data->nind; i++)
    {
      delete [] data->data[i];
    }
  
  delete [] data->data;
  
  data->data = NULL;
  data->nind = -9;
  data->nloci = -9;
  //data->nmiss = -9;
  delete data;
  data = NULL;
  return;
}

void releaseWinData(vector< vector< WinData* >* >* winDataByPopByChr)
{
  for(int pop = 0; pop < winDataByPopByChr->size(); pop++)
    {
      for (int chr = 0; chr < winDataByPopByChr->at(pop)->size(); chr++)
	{
	  releaseWinData(winDataByPopByChr->at(pop)->at(chr));
	}
      winDataByPopByChr->at(pop)->clear();
      delete winDataByPopByChr->at(pop);
    }
  winDataByPopByChr->clear();
  delete winDataByPopByChr;
}

vector< vector< WinData* >* >* initWinData(vector< MapData* >* mapDataByChr,
					   vector< IndData* >* indDataByPop)
{
  vector< vector< WinData* >* >* winDataByPopByChr = new vector< vector< WinData* >* >;
  
  for(int pop = 0; pop < indDataByPop->size(); pop++)
    {
      int nind = indDataByPop->at(pop)->nind;
      vector< WinData* >* winDataByChr = new vector< WinData* >;
      for (int chr = 0; chr < mapDataByChr->size(); chr++)
	{
	  int nloci = mapDataByChr->at(chr)->nloci;
	  WinData* data = initWinData(nind,nloci);
	  winDataByChr->push_back(data);
	}
      winDataByPopByChr->push_back(winDataByChr);
    }
 
  return winDataByPopByChr;
}

void writeWinData(vector< vector< WinData* >* >* winDataByPopByChr,
		  vector< IndData* >* indDataByPop,
		  vector< MapData* >* mapDataByChr,
		  string outfile)
{
  ofstream fout;
  int numPop = indDataByPop->size();
  int numChr = mapDataByChr->size();
  for (int pop = 0; pop < numPop; pop++)
    {
      string popName = indDataByPop->at(pop)->pop;

      for (int chr = 0; chr < numChr; chr++)
	{
	  char chrnum[5];
	  sprintf(chrnum,"%d",mapDataByChr->at(chr)->chr);
	  string rawWinOutfile = outfile;
	  rawWinOutfile += ".";
	  rawWinOutfile += popName;
	  rawWinOutfile += ".chr";
	  rawWinOutfile += chrnum;
	  rawWinOutfile += ".raw.lod.windows";

	  fout.open(rawWinOutfile.c_str());
	  if(fout.fail())
	    {
	      cerr << "ERROR: Failed to open " << rawWinOutfile << " for writing.\n"; 
	      throw -1;
	    }
	  
	  WinData* winData = winDataByPopByChr->at(pop)->at(chr);

	  for (int ind = 0; ind < winData->nind; ind++)
	    {
	      for (int locus = 0; locus < winData->nloci; locus++)
		{
		  if(winData->data[ind][locus] == MISSING) fout << "NA";
		    else fout << winData->data[ind][locus];
		  if(locus < winData->nloci-1) fout << " ";
		}
	      fout << endl;
	    }
	  cerr << "Wrote " << rawWinOutfile << "\n";
	  fout.close();
	}
      
    }

  return;
}




HapData* initHapData(unsigned int nhaps, unsigned int nloci)
{
  if(nhaps < 1 || nloci < 1)
    {
      cerr << "ERROR: Can't allocate WinData object.  Number of haplotypes ("<< nhaps <<") and number of loci ("<< nloci <<") must be positive.\n";
      throw 0;
    }

  HapData* data = new HapData;
  data->nhaps = nhaps;
  data->nloci = nloci;

  data->data = new short*[nhaps];
  for (unsigned int i = 0; i < nhaps; i++)
    {
      data->data[i] = new short[nloci];
      for (unsigned int j = 0; j < nloci; j++)
	{
	  data->data[i][j] = MISSING;
	}
    }

  return data;
}

void releaseHapData(HapData *data)
{
  if(data == NULL) return;
  for (int i = 0; i < data->nhaps; i++)
    {
      delete [] data->data[i];
    }
  
  delete [] data->data;
  
  data->data = NULL;
  data->nhaps = -9;
  data->nloci = -9;
  delete data;
  data = NULL;
  return;
}


int countFields(const string &str)
{
  string::const_iterator it;
  int result;
  int numFields = 0;
  int seenChar = 0;
  for (it = str.begin() ; it < str.end(); it++)
    {
      result = isspace(*it);
      if(result == 0 && seenChar == 0)
	{
	  numFields++;
	  seenChar = 1;
	}
      else if(result != 0)
	{
	  seenChar = 0;
	}
    }
  return numFields;
}



vector< int_pair_t >* scanMapData(string filename, int &numLoci)
{
  ifstream fin;
  cerr << "Scanning " << filename << "...\n";
  fin.open(filename.c_str());
  
  if(fin.fail())
    {
      cerr << "ERROR: Failed to open " << filename << " for reading.\n";
      throw 0;
    }
  
  vector< int_pair_t >* chrStartStop = new vector< int_pair_t >;
  stringstream ss;
  string line;
  int nloci = 0;
  int index;
  int num_cols = 4;
  int current_cols = 0;
  int prevChr = -1;
  int currChr = prevChr;
  int_pair_t currChrCoordinates;
  while(getline(fin,line))
    {
      nloci++;
      index = nloci-1;
      current_cols = countFields(line);
      if(current_cols != num_cols)
	{
	  cerr << "ERROR: line " << nloci << " of " << filename << " has " << current_cols 
	       << ", but expected " << num_cols << ".\n";
	  throw 0;
	}
      ss.str(line);
      ss >> currChr;
      if(prevChr == -1 && index == 0)
	{
	  prevChr = currChr;
	  currChrCoordinates.first = index;
	} 

      if(currChr != prevChr)
	{
	  currChrCoordinates.second = index-1;
	  chrStartStop->push_back(currChrCoordinates);
	  currChrCoordinates.first = index;
	  prevChr = currChr;
	}
    }

  fin.close();

  numLoci = nloci;

  currChrCoordinates.second = index;
  chrStartStop->push_back(currChrCoordinates);


  return chrStartStop;
}

vector< MapData* >* readMapData(string filename, vector< int_pair_t >* chrCoordList)
{
  vector< MapData* >* mapDataByChr = new vector< MapData* >;  

  ifstream fin;
  cerr << "Loading " << filename << "...\n";
  fin.open(filename.c_str());
  
  if(fin.fail())
    {
      cerr << "ERROR: Failed to open " << filename << " for reading.\n";
      throw 0;
    }

  //For each chromosome
  for (int i = 0; i < chrCoordList->size(); i++)
    {
      int size = chrCoordList->at(i).second - chrCoordList->at(i).first + 1;
      MapData* data = initMapData(size);
      for(int locus = 0; locus < size; locus++)
	{
	  fin >> data->chr;
	  fin >> data->locusName[locus];
	  fin >> data->geneticPos[locus];
	  fin >> data->physicalPos[locus];
	}
      cerr << size << " loci on chromosome " << data->chr << endl;
      mapDataByChr->push_back(data);
    }
  
  fin.close();

  return mapDataByChr;
}

vector< int_pair_t >* scanIndData(string filename, int &numInd)
{
  ifstream fin;
  cerr << "Scanning " << filename << "...\n";
  fin.open(filename.c_str());
  
  if(fin.fail())
    {
      cerr << "ERROR: Failed to open " << filename << " for reading.\n";
      throw 0;
    }
  
  vector< int_pair_t >* indStartStop = new vector< int_pair_t >;
  stringstream ss;
  string line;
  int nind = 0;
  int index;
  int num_cols = 2;
  int current_cols = 0;
  string prevPop = "__BLANK__";
  string currPop = prevPop;
  int_pair_t currPopCoordinates;
  while(getline(fin,line))
    {
      nind++;
      index = nind-1;
      current_cols = countFields(line);
      if(current_cols != num_cols)
	{
	  cerr << "ERROR: line " << nind << " of " << filename << " has " << current_cols 
	       << ", but expected " << num_cols << ".\n";
	  throw 0;
	}
      ss.str(line);
      ss >> currPop;
      if(prevPop.compare("__BLANK__") == 0 && index == 0)
	{
	  prevPop = currPop;
	  currPopCoordinates.first = index;
	} 

      if(currPop.compare(prevPop) != 0)
	{
	  currPopCoordinates.second = index-1;
	  indStartStop->push_back(currPopCoordinates);
	  currPopCoordinates.first = index;
	  prevPop = currPop;
	}
    }

  fin.close();

  numInd = nind;

  currPopCoordinates.second = index;
  indStartStop->push_back(currPopCoordinates);

  return indStartStop;
}

vector< IndData* >* readIndData(string filename, vector< int_pair_t >* indCoordList)
{
  vector< IndData* >* indDataByPop = new vector< IndData* >;  

  ifstream fin;
  cerr << "Loading " << filename << "...\n";
  fin.open(filename.c_str());
  
  if(fin.fail())
    {
      cerr << "ERROR: Failed to open " << filename << " for reading.\n";
      throw 0;
    }

  //For each chromosome
  for (int i = 0; i < indCoordList->size(); i++)
    {
      int size = indCoordList->at(i).second - indCoordList->at(i).first + 1;
      IndData* data = initIndData(size);
      for(int ind = 0; ind < size; ind++)
	{
	  fin >> data->pop;
	  fin >> data->indID[ind];
	}
      cerr << size << " individuals in population " << data->pop << endl;
      indDataByPop->push_back(data);
    }
  
  fin.close();

  return indDataByPop;
}

IndData* initIndData(int nind)
{
  if(nind < 1)
    {
      cerr << "ERROR: number of individuals ("<< nind <<") must be positive.\n";
      throw 0;
    }
  
  IndData *data = new IndData;
  data->nind = nind;
  data->indID = new string[nind];

  for (int ind = 0; ind < nind; ind++)
    {
      data->indID[ind] = "--";
    }

  return data;
}

DoubleData* initDoubleData(int n)
{
  DoubleData* data = new DoubleData;
  
  data->size = n;
  data->data = new double[n];
  
  return data;
}

vector < DoubleData* >* convertWinData2DoubleData(vector< vector< WinData* >* >* winDataByPopByChr)
{
  vector < DoubleData* >* rawWinDataByPop = new vector < DoubleData* >;
  double val;
  for (int pop = 0; pop < winDataByPopByChr->size(); pop++)
    {
      int nmiss = 0;
      int ncols = 0;
      int nrows = 0;
      DoubleData* data;
      for(int chr = 0; chr < winDataByPopByChr->at(pop)->size();chr++)
	{
	  for(int ind = 0; ind < winDataByPopByChr->at(pop)->at(chr)->nind; ind++)
	    {
	      for(int locus = 0; locus < winDataByPopByChr->at(pop)->at(chr)->nloci; locus++)
		{
		  val = winDataByPopByChr->at(pop)->at(chr)->data[ind][locus];
		  if(val == MISSING) nmiss++;
		}
	    }

	  ncols += winDataByPopByChr->at(pop)->at(chr)->nloci;
	  nrows = winDataByPopByChr->at(pop)->at(chr)->nind;
	}
      data = initDoubleData(ncols*nrows-nmiss);
      rawWinDataByPop->push_back(data);

      //cerr << "missing: " << nmiss << endl;
    }

  int i;
  for (int pop = 0; pop < winDataByPopByChr->size(); pop++)
    {
      i = 0;
      for(int chr = 0; chr < winDataByPopByChr->at(pop)->size();chr++)
	{
	  for(int ind = 0; ind < winDataByPopByChr->at(pop)->at(chr)->nind; ind++)
	    {
	      for(int locus = 0; locus < winDataByPopByChr->at(pop)->at(chr)->nloci; locus++)
		{
		  val = winDataByPopByChr->at(pop)->at(chr)->data[ind][locus];
		  if(val != MISSING)
		    {
		      rawWinDataByPop->at(pop)->data[i] = val;
		      i++;
		    }
		}
	    }
	}
    }

  return rawWinDataByPop;
}

void releaseDoubleData(DoubleData* data)
{
  delete [] data->data;
  delete data;
  return;
}

void releaseDoubleData(vector < DoubleData* >* rawWinDataByPop)
{
  for(int pop = 0; pop < rawWinDataByPop->size(); pop++)
    {
      releaseDoubleData(rawWinDataByPop->at(pop));
    }
  rawWinDataByPop->clear();
  delete rawWinDataByPop;
  rawWinDataByPop = NULL;
  return;
}
