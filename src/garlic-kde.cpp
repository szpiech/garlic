#include "garlic-kde.h"

double calculateWiggle(KDEResult *kdeResult, int winsize) {
    double tot = 0;
    for (int i = 0; i < kdeResult->size; i++) kdeResult->y[i] = kdeResult->y[i] * 100;
    for (int i = 0; i < kdeResult->size - winsize; i++) {
        double c0, c1, cov00, cov01, cov11, sumsq;
        gsl_fit_linear (&(kdeResult->x[i]), 1, &(kdeResult->y[i]), 1, winsize, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
        tot += sumsq / double(winsize);
    }
    return tot;
}

KDEResult *computeKDE(double *data, int size)
{
    LOG.log("KDE with", size, false);
    LOG.log(" points.");

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

    double h = nrd0(data, size); // bandwitdh
    double min, max;
    gsl_stats_minmax(&min, &max, data, 1, n);
    max += CUT * h;
    min -= CUT * h;


    //cout << "\n\tMin LOD: " << min << "\n\tMax LOD: " << max
    //<< "\n\th: " << h << "\n\tsize: " << size << endl;

    //Results
    double *kde_points = new double[M];
    //Create target array
    double *targets = new double [M];

    //memset(kde_points, 0, sizeof(double)*M);
    //memset(targets, 0, sizeof(double)*M);

    //Initialize the equally spaced target points
    for (int i = 0; i < M; i++)
    {
        kde_points[i] = 0;
        double obs = (double(i + 1) / double(M)) * ( max - min ) + min;
        targets[i] = obs;
        //cout << obs << " " << kde(data,obs,n,h) << endl;
    }

    double targetPointSpacing = targets[1] - targets[0];

    //Weights, man what a waste of memory...
    double *q = new double[n];
    for (int i = 0; i < n; i++)
    {
        q[i] = 1.0 / double(n);
    }

    //figtree not thread safe due to dependent libraries using global vars somewhere
    //pthread_mutex_lock(&kde_mutex);
    figtree( d, n, M, W, data, h, q, targets, epsilon, kde_points );
    //pthread_mutex_unlock(&kde_mutex);

    delete [] q;

    double sum = 0;
    for (int i = 0; i < M; i++)
    {
        sum += kde_points[i];
    }

    for (int i = 0; i < M; i++)
    {
        kde_points[i] /= (sum * targetPointSpacing);
    }

    KDEResult *kdeResult = new KDEResult;
    kdeResult->x = targets;
    kdeResult->y = kde_points;
    kdeResult->size = M;

    return kdeResult;
}


KDEResult *cloneKDEResult(KDEResult *data)
{
    KDEResult *kdeResult = new KDEResult;
    kdeResult->size = data->size;
    kdeResult->x = new double[kdeResult->size];
    kdeResult->y = new double[kdeResult->size];

    for (int i = 0; i < kdeResult->size; i++)
    {
        kdeResult->x[i] = data->x[i];
        kdeResult->y[i] = data->y[i];
    }
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

double nrd0(double x[], const int N)
{
    gsl_sort(x, 1, N);
    double hi = gsl_stats_sd(x, 1, N);
    double iqr =
        gsl_stats_quantile_from_sorted_data (x, 1, N, 0.75) -
        gsl_stats_quantile_from_sorted_data (x, 1, N, 0.25);
    double lo = GSL_MIN(hi, iqr / 1.34);
    double bw = 0.9 * lo * pow(N, -0.2);
    return (bw);
}

double get_min_btw_modes(double *x, double *y, int size, int wsize)
{
    //double initialGuess = 0;
    
    int winsize = 20;
    double maxes;
    double *uniq_maxes = new double[size-winsize];
    double *uniq_counts = new double[size-winsize];
    for(int i = 0; i < size-winsize; i++){
        uniq_counts[i] = 0;
        uniq_maxes[i] = 0;
    }

    int index = 0;
    for(int i = 0; i < size-winsize; i++){
        maxes = y[get_arg_max(&(y[i]),winsize)+i];
        if(i == 1){
            uniq_maxes[i] = maxes;
            uniq_counts[i]++;
        }
        else if(uniq_maxes[index] == maxes){
            uniq_counts[index]++;    
        }
        else if(uniq_maxes[index] != maxes){
            index++;
            uniq_maxes[index] = maxes;
            uniq_counts[index]++;    
        }
    }

    int maxCount = uniq_counts[0];
    int secondMaxCount = 0;
    for(int i = 1; i < size-winsize; i++){
        if(maxCount <= uniq_counts[i]){
            secondMaxCount = maxCount;
            maxCount = uniq_counts[i];
        }
        else if (secondMaxCount <= uniq_counts[i]){
            secondMaxCount = uniq_counts[i];
        }
    }

    vector<double> values;
    for(int i = 0; i < size-winsize; i++){
        if(maxCount == uniq_counts[i] || secondMaxCount == uniq_counts[i]){
            values.push_back(uniq_maxes[i]);
        }
    }

    int maxIndex = -1;
    int secondMaxIndex = -1;
    double firstMax = -1;
    double secondMax = -1;
    for(unsigned int i = 0; i < values.size(); i++){
        if(firstMax <= values[i]){
            secondMax = firstMax;
            firstMax = values[i];
        }
        else if (secondMax <= values[i]){
            secondMax = values[i];
        }
    }

    int leftMaxIndex = -1;
    int rightMaxIndex = -1;

    for(int i = 0; i < size; i++){
        if(y[i] == firstMax){
            leftMaxIndex = i;
        }
        if(y[i] == secondMax){
            rightMaxIndex = i;
        }
    }

    int tmp;
    if(rightMaxIndex < leftMaxIndex){
        tmp = rightMaxIndex;
        rightMaxIndex = leftMaxIndex;
        leftMaxIndex = tmp;
    }

    
    delete [] uniq_maxes;
    delete [] uniq_counts;

    int minIndex = get_arg_min(&(y[leftMaxIndex]), rightMaxIndex - leftMaxIndex + 1) + leftMaxIndex;


    if(abs(x[minIndex]/wsize) < 1) return x[minIndex];
    else return 0;

}

double slope(double x0, double y0, double x1, double y1)
{
    return (y1 - y0) / (x1 - x0);
}

int get_arg_max(double *nums, int size)
{
    double max = numeric_limits<double>::min();
    int arg_max = -1;

    for (int i = 0; i < size; i++)
    {
        if (max < nums[i])
        {
            max = nums[i];
            arg_max = i;
        }
    }

    return arg_max;
}

int get_arg_min(double *nums, int size)
{
    double min = numeric_limits<double>::max();
    int arg_min = -1;
    for (int i = 0; i < size; i++)
    {
        if (min > nums[i])
        {
            min = nums[i];
            arg_min = i;
        }
    }

    return arg_min;
}

void writeKDEResult(KDEResult *kdeResult, string outfile)
{
    ofstream fout;
    fout.open(outfile.c_str());
    if (fout.fail())
    {
        LOG.err("ERROR: Failed to open", outfile);
        throw - 1;
    }

    for (int i = 0; i < kdeResult->size; i++)
    {
        fout << kdeResult->x[i] << " " << kdeResult->y[i] << endl;
    }
    LOG.log("Wrote KDE results to", outfile);
    fout.close();

    return;
}

string makeKDEFilename(string basename, int winsize)
{
    char winStr[10];
    sprintf(winStr, "%d", winsize);
    basename += ".";
    basename += winStr;
    basename += "SNPs.kde";
    return basename;
}