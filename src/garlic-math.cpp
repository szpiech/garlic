#include "garlic-math.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <float.h>
#include <string>
#include "garlic-errlog.h"

using std::string;

double garlicMean(const double *data, size_t n)
{
    //gsl_stats_mean: recurrence, not sum/n.
    double mean = 0.0;
    for (size_t i = 0; i < n; i++) mean += (data[i] - mean) / double(i + 1);
    return mean;
}

double garlicVariance(const double *data, size_t n)
{
    //gsl_stats_variance -> gsl_stats_variance_m with the computed mean, which
    //accumulates by the same style of recurrence and then applies the n/(n-1)
    //correction for the unbiased estimate.
    const double mean = garlicMean(data, n);
    double variance = 0.0;
    for (size_t i = 0; i < n; i++)
    {
        const double delta = data[i] - mean;
        variance += (delta * delta - variance) / double(i + 1);
    }
    return variance * (double(n) / double(n - 1));
}

double garlicSD(const double *data, size_t n)
{
    return sqrt(garlicVariance(data, n));
}

void garlicMinMax(double *min, double *max, const double *data, size_t n)
{
    double mn = data[0];
    double mx = data[0];
    for (size_t i = 0; i < n; i++)
    {
        if (data[i] < mn) mn = data[i];
        if (data[i] > mx) mx = data[i];
    }
    *min = mn;
    *max = mx;
}

double garlicQuantileFromSorted(const double *sorted, size_t n, double f)
{
    const double index = f * double(n - 1);
    const size_t lhs = size_t(index);
    const double delta = index - double(lhs);
    if (n == 0) return 0.0;
    if (lhs == n - 1) return sorted[lhs];
    return (1.0 - delta) * sorted[lhs] + delta * sorted[lhs + 1];
}

void garlicSort(double *data, size_t n)
{
    std::sort(data, data + n);
}

void garlicSortIndex(size_t *p, const double *data, size_t n)
{
    //stable_sort rather than sort so that equal values keep input order and the
    //permutation is reproducible; gsl_sort_index leaves ties unspecified.
    for (size_t i = 0; i < n; i++) p[i] = i;
    std::stable_sort(p, p + n,
                     [data](size_t a, size_t b) { return data[a] < data[b]; });
}

double garlicGaussianPDF(double x, double sigma)
{
    const double u = x / fabs(sigma);
    return (1.0 / (sqrt(2.0 * M_PI) * fabs(sigma))) * exp(-u * u / 2.0);
}

double garlicFitSumsq(const double *x, const double *y, size_t n)
{
    //gsl_fit_linear, keeping only sumsq.  Means by recurrence, as GSL does.
    double m_x = 0.0, m_y = 0.0, m_dx2 = 0.0, m_dxdy = 0.0;
    for (size_t i = 0; i < n; i++)
    {
        m_x += (x[i] - m_x) / double(i + 1);
        m_y += (y[i] - m_y) / double(i + 1);
    }
    for (size_t i = 0; i < n; i++)
    {
        const double dx = x[i] - m_x;
        const double dy = y[i] - m_y;
        m_dx2  += (dx * dx - m_dx2)  / double(i + 1);
        m_dxdy += (dx * dy - m_dxdy) / double(i + 1);
    }
    const double b = m_dxdy / m_dx2;
    double s2 = 0.0;
    for (size_t i = 0; i < n; i++)
    {
        const double dx = x[i] - m_x;
        const double dy = y[i] - m_y;
        const double d = dy - b * dx;
        s2 += d * d;
    }
    return s2;
}

double garlicLogChecked(double x, const char *what)
{
    if (!(x > 0.0))
    {
        LOG.err("ERROR: numerical failure: log of a non-positive value in", string(what), false);
        LOG.err(", value", x);
        LOG.err("\tThis usually means degenerate input to the size-class GMM (for example");
        LOG.err("\ta LOD cutoff so low that every window is called, giving near-identical");
        LOG.err("\tROH lengths). Pass --size-bounds to set the boundaries explicitly.");
        throw 1;
    }
    return log(x);
}

//----------------------------------------------------------------- root finding

//A faithful transcription of GSL's brent_init / brent_iterate / and
//gsl_root_test_interval.  Deliberately kept in GSL's variable names and
//operation order: any algebraic "tidying" changes the last digits of the root.
bool garlicRootBrent(double (*f)(double, void *), void *params,
                     double xLower, double xUpper,
                     int maxIter, double epsrel, bool verbose,
                     double *root)
{
    double a = xLower, b = xUpper, c = xUpper;
    double d = xUpper - xLower, e = xUpper - xLower;
    double fa = f(xLower, params), fb = f(xUpper, params), fc = fb;

    if ((fa < 0.0 && fb < 0.0) || (fa > 0.0 && fb > 0.0)) return false;

    *root = 0.5 * (xLower + xUpper);
    double xLo = xLower, xHi = xUpper;

    for (int iter = 1; iter <= maxIter; iter++)
    {
        //---- brent_iterate ----
        int ac_equal = 0;
        if ((fb < 0 && fc < 0) || (fb > 0 && fc > 0))
        {
            ac_equal = 1;
            c = a; fc = fa; d = b - a; e = b - a;
        }
        if (fabs(fc) < fabs(fb))
        {
            ac_equal = 1;
            a = b; b = c; c = a;
            fa = fb; fb = fc; fc = fa;
        }

        const double tol = 0.5 * DBL_EPSILON * fabs(b);
        const double m = 0.5 * (c - b);

        if (fb == 0)
        {
            *root = b; xLo = b; xHi = b;
            return true;
        }
        if (fabs(m) <= tol)
        {
            *root = b;
            if (b < c) { xLo = b; xHi = c; } else { xLo = c; xHi = b; }
            return true;
        }

        if (fabs(e) < tol || fabs(fa) <= fabs(fb))
        {
            d = m; e = m;                       //bisection
        }
        else
        {
            double p, q, r;
            const double s = fb / fa;
            if (ac_equal) { p = 2 * m * s; q = 1 - s; }
            else
            {
                q = fa / fc;
                r = fb / fc;
                p = s * (2 * m * q * (q - r) - (b - a) * (r - 1));
                q = (q - 1) * (r - 1) * (s - 1);
            }
            if (p > 0) q = -q; else p = -p;

            const double lim = (3 * m * q - fabs(tol * q) < fabs(e * q))
                               ? 3 * m * q - fabs(tol * q) : fabs(e * q);
            if (2 * p < lim) { e = d; d = p / q; }
            else             { d = m; e = m; }  //interpolation failed
        }

        a = b; fa = fb;
        if (fabs(d) > tol) b += d;
        else               b += (m > 0 ? +tol : -tol);
        fb = f(b, params);

        *root = b;
        double cc = c;
        if ((fb < 0 && fc < 0) || (fb > 0 && fc > 0)) cc = a;
        if (b < cc) { xLo = b; xHi = cc; } else { xLo = cc; xHi = b; }

        if (verbose)
            std::cerr << iter << " " << xLo << " " << xHi << " " << *root
                      << " " << (xHi - xLo) << std::endl;

        //---- gsl_root_test_interval(xLo, xHi, 0, epsrel) ----
        double min_abs;
        if ((xLo > 0.0 && xHi > 0.0) || (xLo < 0.0 && xHi < 0.0))
            min_abs = (fabs(xLo) < fabs(xHi)) ? fabs(xLo) : fabs(xHi);
        else
            min_abs = 0.0;
        if (fabs(xHi - xLo) < epsrel * min_abs) return true;
    }
    return false;
}

//----------------------------------------------------------------------- random

GarlicRNG::GarlicRNG(unsigned long int seed)
{
    //gsl_rng_mt19937's mt_set, including its remapping of seed 0.
    if (seed == 0) seed = 4357;
    mt[0] = seed & 0xffffffffUL;
    for (int i = 1; i < 624; i++)
        mt[i] = (1812433253UL * (mt[i - 1] ^ (mt[i - 1] >> 30)) + i) & 0xffffffffUL;
    mti = 624;
}

unsigned long GarlicRNG::next32()
{
    const unsigned long UPPER_MASK = 0x80000000UL;
    const unsigned long LOWER_MASK = 0x7fffffffUL;
    static const unsigned long mag01[2] = {0x0UL, 0x9908b0dfUL};
    unsigned long y;

    if (mti >= 624)
    {
        int kk;
        for (kk = 0; kk < 624 - 397; kk++)
        {
            y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + 397] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (; kk < 624 - 1; kk++)
        {
            y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + (397 - 624)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[624 - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
        mt[624 - 1] = mt[397 - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];
        mti = 0;
    }

    y = mt[mti++];
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);
    return y & 0xffffffffUL;
}

double GarlicRNG::uniform()
{
    //GSL's mt19937 get_double
    return next32() / 4294967296.0;
}

bool GarlicRNG::choose(int *dest, size_t k, const int *src, size_t n)
{
    //gsl_ran_choose: selection sampling, one uniform per candidate.
    if (k > n) return false;
    size_t j = 0;
    for (size_t i = 0; i < n && j < k; i++)
    {
        if ((double(n) - double(i)) * uniform() < (double(k) - double(j)))
        {
            dest[j] = src[i];
            j++;
        }
    }
    return true;
}
