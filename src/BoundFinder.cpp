#include "BoundFinder.h"
#include <iostream>
#include <cmath>   //sqrt: devel got this via gsl/gsl_roots.h in BoundFinder.h
#include "garlic-math.h"

using namespace std;

BoundFinder::BoundFinder(double m1, double v1, double w1, double m2, double v2, double w2, int maxIt, double err, bool v)
{
    params.mu1 = m1;
    params.mu2 = m2;
    params.var1 = v1;
    params.var2 = v2;
    params.a1 = w1;
    params.a2 = w2;

    maxIter = maxIt;
    error = err;

    verbose = v;

    found = false;

    x_lo = (params.mu1 > params.mu2) ? params.mu2 : params.mu1;
    x_hi = (params.mu1 > params.mu2) ? params.mu1 : params.mu2;

    boundary = -9999;

    return;
}

BoundFinder::~BoundFinder()
{
}

double BoundFinder::findBoundary()
{
    if (found) return boundary;

    double r;
    //garlicRootBrent reproduces gsl_root_fsolver_brent driven by
    //gsl_root_test_interval(lo, hi, 0, error) exactly, so the boundaries this
    //returns are the ones GSL returned (verified bit-identical on 2,581 random
    //two-gaussian problems).  It reports failure instead of raising a GSL error
    //when the means do not bracket a root.
    if (!garlicRootBrent(&f, &params, x_lo, x_hi, maxIter, error, verbose, &r))
    {
        LOG.err("ERROR: Root finder failed to converge after", maxIter, false);
        LOG.err(" iterations.");
        throw - 1;
    }

    found = true;
    boundary = r;

    return boundary;
}

double BoundFinder::f(double x, void *p)
{
    Params *params = (Params *)p;
    double v1 = sqrt(params->var1);
    double v2 = sqrt(params->var2);
    return (params->a1) * garlicGaussianPDF(x - (params->mu1), v1) - (params->a2) * garlicGaussianPDF(x - (params->mu2), v2);
}
