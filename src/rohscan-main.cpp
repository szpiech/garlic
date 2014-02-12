#include <iostream>
#include <cstdio>
#include <fstream>
#include <cmath>
#include "rohscan-data.h"
#include "rohscan-roh.h"
#include "rohscan-kde.h";
#include "param_t.h"
#include "gmm.h"
#include "gsl/gsl_statistics.h"
#include "gsl/gsl_sort.h"
#include "BoundFinder.h"

using namespace std;

const string ARG_MAPFILE = "--map";
const string DEFAULT_MAPFILE = "__mapfile";
const string HELP_MAPFILE = "A mapfile with one row per variant site.  Formatted <chr#> <locusID> <genetic pos> <physical pos>";

const string ARG_HAPFILE = "--hap";
const string DEFAULT_HAPFILE = "__hapfile";
const string HELP_HAPFILE = "A hapfile with one row per individual, and one column per variant.  Variants should be coded 0/1/-9.";

const string ARG_INDFILE = "--ind";
const string DEFAULT_INDFILE = "__indfile";
const string HELP_INDFILE = "An indfile containing population and individual IDs.  One row per individual, formatted <popID> <indID>";

const string ARG_OUTFILE = "--out";
const string DEFAULT_OUTFILE = "outfile";
const string HELP_OUTFILE = "The base name for all output files.";

const string ARG_THREADS = "--threads";
const int DEFAULT_THREADS = 1;
const string HELP_THREADS = "The number of threads to spawn during calculations.";

const string ARG_ERROR = "--error";
const double DEFAULT_ERROR = 0.001;
const string HELP_ERROR = "The assumed genotyping error rate.";

const string ARG_WINSIZE = "--winsize";
const int DEFAULT_WINSIZE = 60;
const string HELP_WINSIZE = "The window size in # of SNPs in which to calculate LOD scores.";

const string ARG_POINTS = "--kde-points";
const int DEFAULT_POINTS = 512;
const string HELP_POINTS = "The number of equally spaced points at which to do KDE.";

const string ARG_BW = "--kde-bw";
const double DEFAULT_BW = -1;
const string HELP_BW = "Manually set the bandwidth for the KDE of lod scores.  By default, the nrd0 rule of thumb is used.";

const string ARG_MAX_GAP = "--max-gap";
const int DEFAULT_MAX_GAP = 200000;
const string HELP_MAX_GAP = "A LOD score window is not calculated if the gap in bps between two loci is greater than this parameter.";

const string ARG_RESAMPLE = "--resample";
const int DEFAULT_RESAMPLE = 0;
const string HELP_RESAMPLE = "Number of resamples for estimating allele frequencies.  When set to 0 (default), rohscan will use allele frequencies as calculated from the data.";

int main(int argc, char *argv[])
{

  param_t params;

  params.addFlag(ARG_MAPFILE,DEFAULT_MAPFILE,"mapfileLabel",HELP_MAPFILE);
  params.addFlag(ARG_HAPFILE,DEFAULT_HAPFILE,"hapfileLabel",HELP_HAPFILE);
  params.addFlag(ARG_INDFILE,DEFAULT_INDFILE,"indfileLabel",HELP_INDFILE);
  params.addFlag(ARG_OUTFILE,DEFAULT_OUTFILE,"outfileLabel",HELP_OUTFILE);
  params.addFlag(ARG_THREADS,DEFAULT_THREADS,"threadsLabel",HELP_THREADS);
  params.addFlag(ARG_ERROR,DEFAULT_ERROR,"errorLabel",HELP_ERROR);
  params.addFlag(ARG_WINSIZE,DEFAULT_WINSIZE,"winsizeLabel",HELP_WINSIZE);
  params.addFlag(ARG_POINTS,DEFAULT_POINTS,"pointsLabel",HELP_POINTS);
  params.addFlag(ARG_BW,DEFAULT_BW,"bwLabel",HELP_BW);
  params.addFlag(ARG_MAX_GAP,DEFAULT_MAX_GAP,"maxGapLabel",HELP_MAX_GAP);
  params.addFlag(ARG_RESAMPLE,DEFAULT_RESAMPLE,"resampleLabel",HELP_RESAMPLE);

  try
    {
      params.parseCommandLine(argc,argv);
    }
  catch (...)
    {
      return -1;
    }


  int winsize = params.getIntFlag(ARG_WINSIZE);
  string mapfile = params.getStringFlag(ARG_MAPFILE);
  string hapfile = params.getStringFlag(ARG_HAPFILE);
  string indfile = params.getStringFlag(ARG_INDFILE);
  string outfile = params.getStringFlag(ARG_OUTFILE);
  int numThreads = params.getIntFlag(ARG_THREADS);
  double error = params.getDoubleFlag(ARG_ERROR);
  int MAX_GAP = params.getIntFlag(ARG_MAX_GAP);
  int nresample = params.getIntFlag(ARG_RESAMPLE);
  
  

  if(mapfile.compare(DEFAULT_MAPFILE) == 0)
    {
      cerr << "ERROR: Must specify a mapfile.\n";
      return 1;
    }

  if(hapfile.compare(DEFAULT_HAPFILE) == 0)
    {
      cerr << "ERROR: Must specify a hapfile.\n";
      return 1;
    }

  if(indfile.compare(DEFAULT_INDFILE) == 0)
    {
      cerr << "ERROR: Must specify an indfile.\n";
      return 1;
    }

  if(numThreads <= 0)
    {
      cerr << "ERROR: Number of threads must be > 0.\n";
      return 1;
    }

  if(error <= 0 || error >= 1)
    {
      cerr << "ERROR: Genotype error rate must be > 0 and < 1.\n";
      return 1;
    }

  if(winsize <= 1)
    {
      cerr << "ERROR: SNP window size must be > 1.\n";
      return 1;
    }

  int numLoci, numInd;
  vector< int_pair_t >* chrCoordList;
  vector< MapData* >* mapDataByChr;
    
  vector< int_pair_t >* indCoordList;
  vector< IndData* >* indDataByPop;
    
  vector< vector< HapData* >* >* hapDataByPopByChr;
  vector< vector< FreqData* >* >* freqDataByPopByChr;

  vector< vector< WinData* >* >* winDataByPopByChr;

  try
    {
      chrCoordList = scanMapData(mapfile,numLoci);
      mapDataByChr = readMapData(mapfile,chrCoordList);
      
      indCoordList = scanIndData(indfile,numInd);
      indDataByPop = readIndData(indfile,indCoordList);
      
      hapDataByPopByChr = readHapData(hapfile,numLoci,numInd,chrCoordList,indCoordList);
      
      freqDataByPopByChr = calcFreqData(hapDataByPopByChr,nresample);

      winDataByPopByChr = initWinData(mapDataByChr,indDataByPop);
    }
  catch(...)
    {
      return 1;
    }

  /*
  if(winsize > hapData->nloci)
    {
      cerr << "ERROR: SNP window size (" << winsize << ") must be <= numloci (" << hapData->nloci << ").\n";
      return 1;
    }
  
  if(hapData->nhaps % 2 != 0)
    {
      cerr << "ERROR: There are an odd number of haplotypes.\n";
      return 1;
    }
  */
  int numChr = chrCoordList->size();
  int numPop = indCoordList->size();
  chrCoordList->clear();
  indCoordList->clear();
  delete chrCoordList;
  delete indCoordList;
  
  //Create a vector of pop/chr pairs
  //These will be distributed across threads
  vector<int_pair_t> *popChrPairs = new vector<int_pair_t>;
  int_pair_t pair;
  for(int pop = 0; pop < numPop; pop++)
    {
      for(int chr = 0; chr < numChr; chr++)
	{
	  pair.first = pop;
	  pair.second = chr;
	  popChrPairs->push_back(pair);
	}
    }

  cerr << "There are " << popChrPairs->size() << " pop/chr combinations to compute.\n";

  if(popChrPairs->size() < numThreads)
    {
      numThreads = popChrPairs->size();
      cerr << "WARNING: there are fewer pop/chr pairs than threads requested.  Running with "
	   << numThreads << " threads instead.\n";
    }

  //Partition pop/chr pairs amongst the specified threads
  unsigned long int *NUM_PER_THREAD = new unsigned long int[numThreads];
  unsigned long int div = popChrPairs->size()/numThreads;
  for(int i = 0; i < numThreads; i++) NUM_PER_THREAD[i] = div;
  for(int i = 0; i < (popChrPairs->size())%numThreads;i++) NUM_PER_THREAD[i]++;


  work_order_t *order;
  pthread_t *peer = new pthread_t[numThreads];
  int prev_index = 0;
  for(int i = 0; i < numThreads; i++)
    {
      order = new work_order_t;
      order->first_index = prev_index;
      order->last_index = prev_index+NUM_PER_THREAD[i];
      prev_index += NUM_PER_THREAD[i];

      order->winsize = winsize;
      order->error = error;
      order->MAX_GAP = MAX_GAP;

      order->popChrPairs = popChrPairs;
      order->indDataByPop = indDataByPop;
      order->mapDataByChr = mapDataByChr;
      order->hapDataByPopByChr = hapDataByPopByChr;
      order->freqDataByPopByChr = freqDataByPopByChr;
      order->winDataByPopByChr = winDataByPopByChr;

      order->id = i;
      pthread_create(&(peer[i]),
		     NULL,
		     (void *(*)(void*))scan,
		     (void *)order);
      
    }
 
  for(int i = 0; i < numThreads; i++)
    {
      pthread_join(peer[i],NULL);
    }
  
  releaseHapData(hapDataByPopByChr);
  releaseFreqData(freqDataByPopByChr);

  //Output raw windows
  try
    {
      writeWinData(winDataByPopByChr,indDataByPop,mapDataByChr,outfile);
    }
  catch (...)
    {
      return -1;
    }
  
  ofstream fout;

  //Format the LOD window data into a single array per pop with no missing data
  //Prepped for KDE
  vector < DoubleData* >* rawWinDataByPop = convertWinData2DoubleData(winDataByPopByChr);
  
  //Compute KDE of LOD score distribution
  cerr << "Estimating distribution of raw LOD score windows:\n";
  vector < KDEResult* >* kdeResultByPop = computeKDE(rawWinDataByPop,indDataByPop);
  releaseDoubleData(rawWinDataByPop);

  string boundaryOutfile = outfile;
  boundaryOutfile += ".lod.cutoff";
  
  fout.open(boundaryOutfile.c_str());
  if(fout.fail())
    {
      cerr << "ERROR: Could not open " << boundaryOutfile << " for writing.\n";
      return -1;
    }

  //For each population, find the LOD score cutoff
  double *lodScoreCutoffByPop = new double[numPop];
  for (int pop = 0; pop < numPop; pop++)
    {
      try
	{
	  lodScoreCutoffByPop[pop] = get_min_btw_modes(kdeResultByPop->at(pop)->x,kdeResultByPop->at(pop)->y,512);
	}
      catch(...)
	{
	  cerr << "ERROR: Failed to find the minimum between modes in the LOD score density.\n";
	  return -1;
	}
      
      string popName = indDataByPop->at(pop)->pop;
      
      fout << popName << " " << lodScoreCutoffByPop[pop] << "\n";   
      cerr << popName << " LOD score cutoff: " << lodScoreCutoffByPop[pop] << "\n";     
    }
  cerr << "Wrote " << boundaryOutfile << "\n";
  fout.close();
  
  //Output kde points
  try
    {
      writeKDEResult(kdeResultByPop,indDataByPop,outfile);
    }
  catch (...)
    {
      return -1;
    }
  releaseKDEResult(kdeResultByPop);
  
  cerr << "Begin ROH window assembly.\n";
  //Assemble ROH for each individual in each pop
  vector< ROHLength* >* rohLengthByPop = new vector< ROHLength* >;
  vector< vector< ROHData* >* >* rohDataByPopByInd = assembleROHWindows(winDataByPopByChr,
									mapDataByChr,
									indDataByPop,
									lodScoreCutoffByPop,
									&rohLengthByPop,
									winsize);

  cerr << "Complete.\n";

  
  /*
  ofstream out;
  out.open("test.points.mix");

  int pop = 0;
  //for(int pop = 0; pop < rohLengthByPop->size();pop++)
    {
      for(int i = 0; i < rohLengthByPop->at(pop)->size; i++)
	{
	  out << double(rohLengthByPop->at(pop)->length[i])/1000000.0 << endl;
	}
    }
    
    out.close();
  */



    //return 0;

  //GMM set up for size classifications
  //There might be some benefit in doing this across a range of ngaussians and choosing the classification 
  //That has highest BIC
  int ngaussians = 3;
  size_t maxIter = 1000;
  double tolerance = 1e-8;
  double *W = new double[ngaussians];
  double *Mu = new double[ngaussians];
  double *Sigma = new double[ngaussians];
  double *shortMedBound = new double[rohLengthByPop->size()];
  double *medLongBound = new double[rohLengthByPop->size()];
  size_t *sortIndex = new size_t[ngaussians];

  for(int pop = 0; pop < rohLengthByPop->size();pop++)
    {
      //calculate mean and var for the population size distribution to use for initial guess
      double var = gsl_stats_variance(rohLengthByPop->at(pop)->length,1,rohLengthByPop->at(pop)->size);
      double mu = gsl_stats_mean(rohLengthByPop->at(pop)->length,1,rohLengthByPop->at(pop)->size);
      for(int n = 0; n < ngaussians; n++)
	{
	  W[n] = 1.0/double(ngaussians);
	  Mu[n] = mu*double(n+1)/double(ngaussians+1);
	  //gsl_stats_quantile_from_sorted_data(rohLengthByPop->at(pop)->length, 1, rohLengthByPop->at(pop)->size, double(n+1)/double(ngaussians+1));
	  Sigma[n] = var*(n+1)/double(ngaussians);
	}

      GMM gmm(ngaussians,W,Mu,Sigma,maxIter,tolerance,false);

      gmm.estimate(rohLengthByPop->at(pop)->length,rohLengthByPop->at(pop)->size);
      
      for(int n = 0; n < ngaussians; n++)
	{
	  W[n] = gmm.getMixCoefficient(n);
	  Mu[n] = gmm.getMean(n);
	  Sigma[n] = gmm.getVar(n);
	  sortIndex[n] = n;
	}
      
      gsl_sort_index(sortIndex,Mu,1,ngaussians);

      cerr << rohLengthByPop->at(pop)->pop << " [" 
	   << W[sortIndex[0]] << " ," << W[sortIndex[1]] << " ," << W[sortIndex[2]] << "] [" 
	   << Mu[sortIndex[0]] << " ," << Mu[sortIndex[1]] << " ," << Mu[sortIndex[2]] << "] [" 
	   << Sigma[sortIndex[0]] << " ," << Sigma[sortIndex[1]] << " ," << Sigma[sortIndex[2]] << "]\n";
      
      //Find boundaries, there are ngaussians-1 of them, but for the moment this is defined to be 2
      //This finds the 'first' root of the difference between two gaussians
      BoundFinder SM(Mu[sortIndex[0]],Sigma[sortIndex[0]],W[sortIndex[0]],Mu[sortIndex[1]],Sigma[sortIndex[1]],W[sortIndex[1]],1000,1e-4,false);
      shortMedBound[pop] = SM.findBoundary();
      BoundFinder ML(Mu[sortIndex[1]],Sigma[sortIndex[1]],W[sortIndex[1]],Mu[sortIndex[2]],Sigma[sortIndex[2]],W[sortIndex[2]],1000,1e-4,false);
      medLongBound[pop] = ML.findBoundary();

      cerr << "A/B: " << shortMedBound[pop] << " B/C: " << medLongBound[pop] << endl;

      /*
      cerr << rohLengthByPop->at(pop)->pop << " [" 
	   << W[0] << " ," << W[1] << " ," << W[2] << "] [" 
	   << Mu[0] << " ," << Mu[1] << " ," << Mu[2] << "] [" 
	   << Sigma[0] << " ," << Sigma[1] << " ," << Sigma[2] << "]\n";
      */
    }

  //Output ROH calls to file, one for each individual
  //includes A/B/C size classifications
  //Could be modified to allow for arbitrary number of size classifications
  for (int pop = 0; pop < rohDataByPopByInd->size(); pop++)
    {
      vector< ROHData* >* rohDataByInd = rohDataByPopByInd->at(pop);
      for(int ind = 0; ind < rohDataByInd->size(); ind++)
	{
	  ROHData* rohData = rohDataByInd->at(ind);
	  cerr << "Writing ROH tracts for " << rohData->indID << endl;
	  string rohOutfile = outfile;
	  rohOutfile += ".";
	  rohOutfile += rohData->indID;
	  rohOutfile += ".roh";
	  ofstream out;
	  out.open(rohOutfile.c_str());
	  if(out.fail())
	    {
	      cerr << "ERROR: Failed to open " << rohOutfile << " for writing.\n";
	      return -1;
	    }
	  
	  for (int roh = 0; roh < rohData->chr.size(); roh++)
	    {
	      int size = (rohData->stop[roh] - rohData->start[roh]);
	      char sizeClass = 'C';
	      if(size < shortMedBound[pop]) sizeClass = 'A';
	      else if(size < medLongBound[pop]) sizeClass = 'B';
	      out << rohData->chr[roh] << " " << rohData->start[roh] << " " << rohData->stop[roh] << " " << size << " " << sizeClass << endl;
	    }
	  
	  out.close();
	}
    }


  return 0;
}

