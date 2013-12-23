#include "rohscan-kde.h"

KDEResult* computeKDE(double* data,int size)
{
  
 //Used, as in the R function density, to extend the range of the fixed width points
  //used to compute the KDE
  double CUT = 3;

 // The dimensionality of each sample vector.
  int d = 1;

  // The number of targets (vectors at which gauss transform is evaluated).
  int M = 512;

  // The number of sources which will be used for the gauss transform. 
  int n = size;

  // Desired maximum absolute error after normalizing output by sum of weights.
  // If the weights, q_i (see below), add up to 1, then this is will be the 
  // maximum absolute error.
  // The smaller epsilon is, the more accurate the results will be, at the
  // expense of increased computational complexity.
  double epsilon = 1e-2;
  
  // Number of weights.  For each set of weights a different Gauss Transform is computed, 
  // but by giving multiple sets of weights at once some overhead can be shared.
  int W = 1;

  double h = nrd0(data,size); // bandwitdh
  double min, max;
  gsl_stats_minmax(&min,&max,data,1,n);
  max += CUT*h;
  min -= CUT*h;

  cerr << "\n\tMin LOD: " << min << "\n\tMax LOD: " << max
       << "\n\th: " << h << "\n";

  //Results
  double *kde_points = new double[M];
  //Create target array
  double *targets = new double [M];
  
  memset(kde_points, 0, sizeof(double)*M);
  memset(targets, 0, sizeof(double)*M);

  //Initialize the equally spaced target points
  for(int i = 0; i < M; i++)
    {
      double obs = (double(i+1)/double(M)) * ( max - min ) + min;
      targets[i] = obs;
      //cout << obs << " " << kde(data,obs,n,h) << endl;
    }
  
  double targetPointSpacing = targets[1] - targets[0];
  
  //Weights, man what a waste of memory...
  double *q = new double[n];
  for(int i = 0; i < n; i++)
    {
      q[i] = 1.0/double(n);
    }
  
  figtree( d, n, M, W, data, h, q, targets, epsilon, kde_points );

  delete [] q;

  double sum = 0;
  for(int i = 0; i < M; i++)
    {
      sum += kde_points[i];
    }
  
  for(int i = 0; i < M; i++)
    {
      kde_points[i]/=(sum*targetPointSpacing);
    }
  
  KDEResult* kdeResult = new KDEResult;
  kdeResult->x = targets;
  kdeResult->y = kde_points;
  kdeResult->size = M;

  return kdeResult;
}

void releaseKDEResult(KDEResult *data)
{
  delete [] data->x;
  delete [] data->y;
  delete data;
  data = NULL;
  return;
}

vector < KDEResult* >* computeKDE(vector < DoubleData* >* rawWinDataByPop, vector< IndData* >* indDataByPop)
{
  vector < KDEResult* >* kdeResultByPop = new  vector < KDEResult* >;

  for(int pop = 0; pop < rawWinDataByPop->size(); pop++)
    {
      cerr << "\t" << indDataByPop->at(pop)->pop << "...";
      KDEResult* result = computeKDE(rawWinDataByPop->at(pop)->data,rawWinDataByPop->at(pop)->size);
      cerr << "\tdone.\n";
      kdeResultByPop->push_back(result);
    }
  return kdeResultByPop;  
}


void releaseKDEResult(vector < KDEResult* >* kdeResultByPop)
{
  for(int pop = 0; pop < kdeResultByPop->size(); pop++)
    {
      releaseKDEResult(kdeResultByPop->at(pop));
    }
  delete kdeResultByPop;
  kdeResultByPop = NULL;
  return;
}

double nrd0(double x[], const int N)
{
	gsl_sort(x, 1, N);
	double hi = gsl_stats_sd(x, 1, N);
	double iqr =
		gsl_stats_quantile_from_sorted_data (x,1, N,0.75) - 
        gsl_stats_quantile_from_sorted_data (x,1, N,0.25);
	double lo = GSL_MIN(hi, iqr/1.34);
	double bw = 0.9 * lo * pow(N,-0.2);
	return(bw);
}
/*
double get_min_btw_modes(double *x, double *y, int size)
{
  vector<int> local_max_index;
  
  //The goal here is to first 'smooth' the distribution
  //And then find the two local maxima of the (presumably) bimodal distribution
  //Using these as bounds, search for the minimum between the two maxima
  //This should be the cutoff we desire

  //Stride is the 'smoothing' parameter
  //Too small and we find numerous maxima because of the wobbly KDE
  //Too large and we will only find the global max
  //So here we cycle over a range of options,
  //breaking when we find only two local maxima
  //***I don't think we can guarantee that this is the one we want
  //***will have to figure out a way to determine that
  //***however it ought to be the case the majority of cases
  for (int stride = 1; stride <= size/2; stride++)
    {
      vector<int> max_index_list;
      int current_value = -1;
      
      for (int i = 0; i < size-stride; i++)
	{
	  int arg_max = get_arg_max(&(y[i]),stride);
	  
	  if(arg_max+i != current_value)
	    {
	      current_value = arg_max+i;
	      max_index_list.push_back(current_value);
	    }
	}
     
      int previous_sign = (slope(x[max_index_list[0]],y[max_index_list[0]],x[max_index_list[1]],y[max_index_list[1]]) > 0) ? 1 : 0;
      int current_sign = 0;
      
      for(int i = 1; i < max_index_list.size()-1; i++)
	{
	  int index = max_index_list[i];
	  int index1 = max_index_list[i+1];
	  current_sign = (slope(x[index],y[index],x[index1],y[index1]) > 0) ? 1 : 0;
	  
	  if(previous_sign-current_sign == 1)//then it is a max
	    { 
	      //cerr << "local max at: " << index << " " << x[index] << " " << y[index] << endl;
	      local_max_index.push_back(index);
	    }
	  
	  previous_sign = current_sign;
	}

      if(local_max_index.size() == 2)
	{
	  break;
	}
      else
	{
	  local_max_index.clear();
	}
    }

  if(local_max_index.size() != 2)
    {
      throw -1;
    }

  int arg_min = get_arg_min(&(y[local_max_index[0]]),(local_max_index[1]-local_max_index[0])+1) + local_max_index[0];

  return x[arg_min];
}
*/

double get_min_btw_modes(double *x, double *y, int size)
{
  //double initialGuess = 0;
  double *x_abs = new double[size];

  for(int i = 0; i < size; i++)
    {
      x_abs[i] = fabs(x[i]);
    }

  int minGuessIndex = get_arg_min(x_abs,size);

  delete [] x_abs;

  int rightMaxIndex = get_arg_max(&(y[minGuessIndex]),size-minGuessIndex) + minGuessIndex;
  int leftMaxIndex = get_arg_max(&(y[0]),minGuessIndex+1);
  int minIndex = get_arg_min(&(y[leftMaxIndex]),rightMaxIndex-leftMaxIndex+1) + leftMaxIndex;
  
  return x[minIndex];

}

double slope(double x0, double y0, double x1, double y1)
{
  return (y1-y0)/(x1-x0);  
}

int get_arg_max(double* nums, int size)
{
  double max = numeric_limits<double>::min();
  int arg_max = -1;

  for (int i = 0; i < size; i++)
    {
      if(max < nums[i])
	{
	  max = nums[i];
	  arg_max = i;
	}
    }

  return arg_max;
}

int get_arg_min(double* nums, int size)
{
  double min = numeric_limits<double>::max();
  int arg_min = -1;
  for (int i = 0; i < size; i++)
    {
      if(min > nums[i])
	{
	  min = nums[i];
	  arg_min = i;
	}
    }

  return arg_min;
}

void writeKDEResult(vector < KDEResult* >* kdeResultByPop, vector< IndData* >* indDataByPop, string outfile)
{
  ofstream fout;
  int numPop = indDataByPop->size();
  for (int pop = 0; pop < numPop; pop++)
    { 
      string popName = indDataByPop->at(pop)->pop;
      string kdeOutfile = outfile;
      kdeOutfile += ".";
      kdeOutfile += popName;
      kdeOutfile += ".kde";
      
      fout.open(kdeOutfile.c_str());
      if(fout.fail())
	{
	  cerr << "ERROR: Failed to open " << kdeOutfile << " for writing.\n";
	  throw -1;
	}

      for(int i = 0; i < kdeResultByPop->at(pop)->size; i++)
	{
	  fout << kdeResultByPop->at(pop)->x[i] << " " << kdeResultByPop->at(pop)->y[i] << endl;
	}
      cerr << "Wrote " << kdeOutfile << endl;
      fout.close();
    }

  return;
}
