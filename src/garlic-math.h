#ifndef __GARLIC_MATH_H__
#define __GARLIC_MATH_H__

#include <cstddef>

//Replacements for the handful of GSL routines garlic used, so that building it
//needs no external numerical library.  Only ~10 GSL functions were ever called,
//against four per-platform static archives (libgsl.a, libgslcblas.a) that were
//the largest source of build friction for a tool whose users compile from
//source on clusters.
//
//These are written to reproduce GSL's results BIT FOR BIT, not merely to
//compute the same quantity, because garlic's regression suite compares output
//checksums.  That matters most for the mean and variance: GSL does not sum and
//divide, it uses incremental recurrences
//
//    mean_{i+1}     = mean_i + (x_i - mean_i) / (i + 1)
//    variance_{i+1} = variance_i + (delta*delta - variance_i) / (i + 1)
//
//which give different last digits from the naive formulas.  The comments below
//name the GSL routine each function stands in for.

//gsl_stats_mean
double garlicMean(const double *data, size_t n);
//gsl_stats_variance (unbiased; the recurrence above, scaled by n/(n-1))
double garlicVariance(const double *data, size_t n);
//gsl_stats_sd
double garlicSD(const double *data, size_t n);
//gsl_stats_minmax
void garlicMinMax(double *min, double *max, const double *data, size_t n);
//gsl_stats_quantile_from_sorted_data.  GSL's definition: index = f*(n-1), then
//linear interpolation between the two neighbouring order statistics.
double garlicQuantileFromSorted(const double *sorted, size_t n, double f);
//gsl_sort
void garlicSort(double *data, size_t n);
//gsl_sort_index (ascending permutation)
void garlicSortIndex(size_t *p, const double *data, size_t n);
//gsl_ran_gaussian_pdf(x, sigma)
double garlicGaussianPDF(double x, double sigma);
//The only output of gsl_fit_linear garlic used was sumsq, the residual sum of
//squares of the least-squares line.  Computed GSL's way: recurrence means, then
//sum (dy - slope*dx)^2.
double garlicFitSumsq(const double *x, const double *y, size_t n);

//gsl_sf_log's domain check, which garlic depended on without meaning to.  GSL
//raised a domain error for x <= 0, and garlic installed a handler that turned
//that into a clean exit-2 with a diagnosis.  Plain log() returns NaN instead,
//which used to let a collapsed GMM component run all 1000 EM iterations and
//surface much later as "root finder failed to converge" -- a symptom, not the
//cause.  Throws int 1, which is what the GMM's callers already catch.
double garlicLogChecked(double x, const char *what);

//M_PI is NOT standard C++.  It is an X/Open extension: glibc and Apple's libc
//define it from <cmath> unconditionally, but MinGW follows MSVC and defines it
//only under _USE_MATH_DEFINES, so a MinGW build fails with "'M_PI' was not
//declared in this scope".  devel got it from gsl/gsl_math.h, which defines it;
//dropping GSL in 0bdc179 removed that, the same way it removed the transitive
//<math.h> behind the unqualified isnan (9724178).
//
//Defining it here rather than setting _USE_MATH_DEFINES: that macro only works
//if it precedes the FIRST <cmath> anywhere in the translation unit, which is a
//property of include order that nothing enforces.  garlic-math.h is the
//replacement for gsl_math.h, so the constant belongs with it.
const double GARLIC_PI = 3.14159265358979323846;

//----------------------------------------------------------------- root finding
//Brent's method on a bracketing interval, reproducing gsl_root_fsolver_brent
//driven by gsl_root_test_interval(lo, hi, 0, epsrel).  Reproduced rather than
//written afresh so the GMM size-class boundaries are unchanged: a different
//root finder converges to a different point inside the tolerance, which would
//move every published boundary slightly.
//
//Returns true on convergence and writes the root to *root; false if the
//endpoints do not straddle zero (which gsl_root_fsolver_set treated as a GSL
//error) or if maxIter is reached.  verbose mirrors the old per-iteration trace.
bool garlicRootBrent(double (*f)(double, void *), void *params,
                     double xLower, double xUpper,
                     int maxIter, double epsrel, bool verbose,
                     double *root);

//----------------------------------------------------------------------- random
//Stands in for gsl_rng with gsl_rng_default, which is gsl_rng_mt19937.  std's
//mersenne_twister_engine is the same algorithm with the same seeding
//recurrence, and GSL's uniform for mt19937 is next32 / 2^32, so the stream is
//identical -- seeded subsampling gives exactly the results it did before.
//GSL mapped seed 0 to 4357; that is reproduced too.
class GarlicRNG
{
public:
    explicit GarlicRNG(unsigned long int seed);
    //gsl_rng_uniform
    double uniform();
    //gsl_ran_choose: selection sampling of k of n ints, order preserved
    bool choose(int *dest, size_t k, const int *src, size_t n);
private:
    //MT19937 state, written out rather than using std::mt19937 so the mapping
    //from seed to stream cannot change with a standard-library version.
    unsigned long mt[624];
    int mti;
    unsigned long next32();
};

#endif
