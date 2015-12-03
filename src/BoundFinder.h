/*
  A class to handle finding the 'boundary' between two gaussians



  Copyright (C) 2013  Zachary A Szpiech (szpiech@gmail.com)

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef __BOUNDFINDER_H__
#define __BOUNDFINDER_H__

#include "gsl/gsl_roots.h"
#include "garlic-errlog.h"

class BoundFinder
{
public:

    double findBoundary();

    BoundFinder(double mu1, double var1, double a1, double mu2, double var2, double a2, int maxIt = 1000, double err = 1e-5, bool v = true);
    ~BoundFinder();

private:

    struct Params
    {
        double mu1;
        double mu2;
        double var1;
        double var2;
        //mixture coefficients
        double a1;
        double a2;
    };

    Params params;

    double x_hi;
    double x_lo;

    double boundary;

    int maxIter;
    double error;
    bool verbose;

    bool found;

    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;

    gsl_function F;
    static double f(double x, void *p);

};

#endif
