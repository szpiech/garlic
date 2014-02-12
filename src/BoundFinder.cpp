#include "BoundFinder.h"
#include <iostream>
#include "gsl/gsl_randist.h"

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

  F.function = &f;
  F.params = &params;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);

  boundary = -9999;

  return;
}

BoundFinder::~BoundFinder()
{
  gsl_root_fsolver_free (s);
}

double BoundFinder::findBoundary()
{
  if(found) return boundary;

  int status;
  int iter = 0;
  double r;
  
  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      r = gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi, 0, error);
      
      if(verbose)
	{
	  cerr << iter << " " << x_lo << " " 
	       << x_hi << " " << r << " " << (x_hi - x_lo) << endl;
	}

    } while (status == GSL_CONTINUE && iter < maxIter);
  
  if(status != GSL_SUCCESS)
    {
      cerr << "Root finder failed to converge after " << maxIter << " iterations.\n";
      throw -1;
    }
  
  found = true;
  boundary = r;
  
  return boundary;

}

double BoundFinder::f(double x, void *p)
{
  Params *params = (Params*)p;
  double v1 = sqrt(params->var1);
  double v2 = sqrt(params->var2);
  return (params->a1)*gsl_ran_gaussian_pdf(x-(params->mu1),v1) - (params->a2)*gsl_ran_gaussian_pdf(x-(params->mu2),v2);
}
