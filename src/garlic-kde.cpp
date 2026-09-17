#include <algorithm>
#include <thread>
#include "garlic-kde.h"

static int KDE_POINTS = 512;
static double KDE_CUT = 3;
static int KDE_MODE_SPAN = 20;

void setKDEGrid(int points, double cut) { KDE_POINTS = points; KDE_CUT = cut; }
void setModeSpan(int span) { KDE_MODE_SPAN = span; }


double calculateWiggle(KDEResult *kdeResult, int winsize) {
    double tot = 0;
    for (int i = 0; i < kdeResult->size; i++) kdeResult->y[i] = kdeResult->y[i] * 100;
    for (int i = 0; i < kdeResult->size - winsize; i++) {
        double sumsq;
        sumsq = garlicFitSumsq(&(kdeResult->x[i]), &(kdeResult->y[i]), winsize);
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
    double CUT = KDE_CUT;

    // The number of targets (vectors at which gauss transform is evaluated).
    int M = KDE_POINTS;

    // The number of sources which will be used for the gauss transform.
    int n = size;

    double h = nrd0(data, size); // bandwitdh
    double min, max;
    garlicMinMax(&min, &max, data, n);
    max += CUT * h;
    min -= CUT * h;


    //cout << "\n\tMin LOD: " << min << "\n\tMax LOD: " << max
    //<< "\n\th: " << h << "\n\tsize: " << size << endl;

    //Built in place in the result rather than as raw arrays handed over to it,
    //so nothing owns these but the KDEResult.
    KDEResult *kdeResult = new KDEResult;
    kdeResult->size = M;
    kdeResult->y.assign(M, 0.0);          //kde_points
    kdeResult->x.resize(M);               //targets
    double *kde_points = kdeResult->y.data();
    double *targets = kdeResult->x.data();

    //Initialize the equally spaced target points
    for (int i = 0; i < M; i++)
    {
        double obs = (double(i + 1) / double(M)) * ( max - min ) + min;
        targets[i] = obs;
    }

    double targetPointSpacing = targets[1] - targets[0];

    kdeGaussian(data, n, h, targets, M, kde_points);

    double sum = 0;
    for (int i = 0; i < M; i++)
    {
        sum += kde_points[i];
    }

    for (int i = 0; i < M; i++)
    {
        kde_points[i] /= (sum * targetPointSpacing);
    }

    return kdeResult;
}


KDEResult *cloneKDEResult(KDEResult *data)
{
    KDEResult *kdeResult = new KDEResult;
    kdeResult->size = data->size;
    kdeResult->x.resize(kdeResult->size);
    kdeResult->y.resize(kdeResult->size);

    for (int i = 0; i < kdeResult->size; i++)
    {
        kdeResult->x[i] = data->x[i];
        kdeResult->y[i] = data->y[i];
    }
    return kdeResult;
}

void releaseKDEResult(KDEResult *data)
{
    delete data;
    data = NULL;
    return;
}

double nrd0(double x[], const int N)
{
    garlicSort(x, N);
    double hi = garlicSD(x, N);
    double iqr =
        garlicQuantileFromSorted(x, N, 0.75) -
        garlicQuantileFromSorted(x, N, 0.25);
    double lo = (hi < iqr / 1.34) ? hi : iqr / 1.34;
    double bw = 0.9 * lo * pow(N, -0.2);
    return (bw);
}

//============================ exact Gaussian KDE =============================
//Replaces the figtree call that used to live in computeKDE.
//
//figtree's default (FIGTREE_EVAL_AUTO) was non-deterministic -- its k-center
//clustering reseeds itself with srand(time(NULL)) -- and FIGTREE_EVAL_DIRECT,
//which is deterministic and exact, is a flat O(n*M) double loop.  At the
//defaults that is 0.67 s of a 2.87 s run on the bundled example, and ~31 s
//under --no-kde-thinning (24.5 M points); it grows linearly with cohort size.
//
//Two observations make the exact sum much cheaper:
//
//  1. exp(-z^2/h^2) UNDERFLOWS TO EXACTLY 0.0 in IEEE double once
//     z^2/h^2 > 745.2.  So every source further than h*sqrt(746) from a target
//     contributes exactly 0.0, and skipping it is not an approximation --
//     adding 0.0 cannot change a sum.  This is exactness by construction, not
//     a tolerance.
//  2. nrd0() sorts `data` in place (garlicSort) immediately before we need it, so
//     the contributing range for each target is one binary search.
//
//Targets are independent, so the 512 grid points parallelise with no sharing.
static int KDE_NUM_THREADS = 1;

void setKDEThreads(int n) { KDE_NUM_THREADS = (n > 0 ? n : 1); }

struct KDE_work_order_t
{
    const double *data;
    int n;
    double h;
    const double *targets;
    double *out;
    int start;
    int stop;
    int stride;
    double R;
};

//Takes its order type directly.  Under pthread_create every worker had to be
//void *(void *) and cast the argument back, and pthread_join's return value
//was discarded at all seven call sites -- so the void * round trip carried no
//information and only cost the compiler its type checking.
static void parallelKDE(KDE_work_order_t *p)
{
    const double *data = p->data;
    const double *targets = p->targets;
    const int n = p->n;
    const double invh2 = 1.0 / (p->h * p->h);
    const double R = p->R;
    const double invn = 1.0 / double(n);

    //Strided rather than blocked: the contributing window is widest in the
    //dense centre of the grid, so contiguous blocks hand one thread most of
    //the work.
    for (int i = p->start; i < p->stop; i += p->stride)
    {
        const double c = targets[i];
        //only sources in [c-R, c+R] can contribute a nonzero term
        const double *lo = std::lower_bound(data, data + n, c - R);
        const double *hi = std::upper_bound(data, data + n, c + R);
        double sum = 0.0;
        for (const double *q = lo; q < hi; q++)
        {
            const double z = c - *q;
            sum += exp(-z * z * invh2);
        }
        p->out[i] = sum * invn;
    }
}

//out[i] = (1/n) * sum_j exp(-(targets[i]-data[j])^2 / h^2), the same quantity
//figtree was computing with unit weights q_j = 1/n.
void kdeGaussian(double *data, int n, double h, const double *targets, int M, double *out)
{
    if (n < 1 || M < 1) return;

    //nrd0() leaves data sorted; verify rather than assume, and sort if some
    //future caller changes that.
    bool sorted = true;
    for (int i = 1; i < n; i++) { if (data[i] < data[i - 1]) { sorted = false; break; } }
    if (!sorted) garlicSort(data, n);

    //exp() underflows to exactly 0.0 beyond this separation, so terms outside
    //contribute nothing at all -- not merely something small.
    const double R = h * sqrt(746.0);

    int nt = KDE_NUM_THREADS;
    if (nt > M) nt = M;
    if (nt < 1) nt = 1;

    KDE_work_order_t proto;
    proto.data = data; proto.n = n; proto.h = h; proto.targets = targets;
    proto.out = out; proto.R = R; proto.start = 0; proto.stop = M; proto.stride = 1;

    if (nt == 1)
    {
        parallelKDE(&proto);
        return;
    }

    //orders is filled completely BEFORE any thread is launched: it is a vector,
    //so emplacing into it while a thread holds &orders[i] would dangle on
    //reallocation.
    vector<KDE_work_order_t> orders(nt);
    for (int i = 0; i < nt; i++)
    {
        orders[i] = proto;
        orders[i].start = i;
        orders[i].stride = nt;
    }
    vector<std::thread> peer;
    peer.reserve(nt);
    for (int i = 0; i < nt; i++) peer.emplace_back(parallelKDE, &(orders[i]));
    for (int i = 0; i < nt; i++) peer[i].join();
}


double get_min_btw_modes(double *x, double *y, int size, int wsize)
{
    //double initialGuess = 0;
    
    int winsize = KDE_MODE_SPAN;
    double maxes;
    vector<double> uniq_maxes(size-winsize);
    vector<double> uniq_counts(size-winsize);
    for(int i = 0; i < size-winsize; i++){
        uniq_counts[i] = 0;
        uniq_maxes[i] = 0;
    }

    //Record WHICH grid point was the windowed maximum alongside its value, so
    //the mode's location never has to be recovered by scanning y for a float
    //that compares exactly equal (see below).
    vector<int> uniq_argmax(size-winsize);
    for(int i = 0; i < size-winsize; i++) uniq_argmax[i] = -1;

    int index = 0;
    for(int i = 0; i < size-winsize; i++){
        int argmax_i = get_arg_max(&(y[i]),winsize)+i;
        maxes = y[argmax_i];
        if(i == 1){
            uniq_maxes[i] = maxes;
            uniq_argmax[i] = argmax_i;
            uniq_counts[i]++;
        }
        else if(uniq_maxes[index] == maxes){
            uniq_counts[index]++;    
        }
        else if(uniq_maxes[index] != maxes){
            index++;
            uniq_maxes[index] = maxes;
            uniq_argmax[index] = argmax_i;
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
    vector<int> valueIdx;
    for(int i = 0; i < size-winsize; i++){
        if((maxCount == uniq_counts[i] || secondMaxCount == uniq_counts[i]) && uniq_argmax[i] >= 0){
            values.push_back(uniq_maxes[i]);
            valueIdx.push_back(uniq_argmax[i]);
        }
    }

    //Two distinct candidate modes are required; with fewer there is no
    //"between the modes" to minimise over.  This used to leave
    //leftMaxIndex/rightMaxIndex at -1 and index y[-1].
    if(values.size() < 2){
        LOG.err("ERROR: Found fewer than two modes in the LOD score density (", int(values.size()), false);
        LOG.err(" candidate(s)).");
        throw 0;
    }

    double firstMax = -1;
    double secondMax = -1;
    int leftMaxIndex = -1;
    int rightMaxIndex = -1;
    for(unsigned int i = 0; i < values.size(); i++){
        if(firstMax <= values[i]){
            secondMax = firstMax;
            rightMaxIndex = leftMaxIndex;
            firstMax = values[i];
            leftMaxIndex = valueIdx[i];
        }
        else if (secondMax <= values[i]){
            secondMax = values[i];
            rightMaxIndex = valueIdx[i];
        }
    }


    if(leftMaxIndex < 0 || rightMaxIndex < 0){
        LOG.err("ERROR: Could not locate two modes in the LOD score density.");
        throw 0;
    }

    if(rightMaxIndex < leftMaxIndex){
        int tmp = rightMaxIndex;
        rightMaxIndex = leftMaxIndex;
        leftMaxIndex = tmp;
    }

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
    //::min() is the smallest POSITIVE normal double (+2.2e-308), not -inf, so
    //an all-nonpositive input returned arg_max = -1 and the caller indexed
    //y[-1].  ::lowest() is the intended value.  (get_arg_min below correctly
    //uses ::max(), which is what makes the asymmetry easy to miss.)
    double max = numeric_limits<double>::lowest();
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