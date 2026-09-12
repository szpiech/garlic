#ifndef __GARLIC_KDE_H__
#define __GARLIC_KDE_H__

#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <cstring>
#include <limits>
//#include <pthread.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_fit.h>
#include "garlic-data.h"
#include "garlic-errlog.h"

struct KDEResult
{
    int size;
    double *x;
    double *y;
};

double nrd0(double *data, const int n);

KDEResult *computeKDE(double *data, int size);
//Exact Gaussian KDE on a fixed target grid; replaces the figtree dependency.
void kdeGaussian(double *data, int n, double h, const double *targets, int M, double *out);
void setKDEThreads(int n);
//Grid resolution, range extension (in bandwidths) and the span used to smooth
//the density when locating its modes.  All three were hardcoded (512, 3, 20)
//and all three affect the selected LOD cutoff.
void setKDEGrid(int points, double cut);
void setModeSpan(int span);
KDEResult *cloneKDEResult(KDEResult *data);
void releaseKDEResult(KDEResult *data);
void writeKDEResult(KDEResult *kdeResult, string outfile);
string makeKDEFilename(string basename, int winsize);

int get_arg_max(double *nums, int size);
int get_arg_min(double *nums, int size);
double get_min_btw_modes(double *x, double *y, int size, int wsize);
double slope(double x0, double y0, double x1, double y1);

double calculateWiggle(KDEResult *kdeResult, int size = 20);

/*
KDEWinsizeReport *initKDEWinsizeReport();
void releaseKDEWinsizeReport(KDEWinsizeReport *winsizeReport);
*/
//extern pthread_mutex_t io_mutex;

#endif
