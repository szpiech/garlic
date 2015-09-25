#ifndef __GARLIC_KDE_H__
#define __GARLIC_KDE_H__

#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <cstring>
#include <limits>
#include <pthread.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_fit.h>
#include "figtree.h"
#include "garlic-data.h"

struct KDEResult
{
    int size;
    double *x;
    double *y;
};

struct KDEWork
{
	vector < DoubleData * > *rawWinDataByPop;
	int id;
	int numThreads;
	vector< IndData * > *indDataByPop;
	vector < KDEResult * > *kdeResultByPop;
};

double nrd0(double *data, const int n);

KDEResult *computeKDE(double *data, int size);
vector < KDEResult * > *computeKDE(vector < DoubleData * > *rawWinDataByPop, vector< IndData * > *indDataByPop, int numThreads);
void releaseKDEResult(KDEResult *data);
void releaseKDEResult(vector < KDEResult * > *kdeResultByPop);
void writeKDEResult(vector < KDEResult * > *kdeResultByPop, vector< IndData * > *indDataByPop, string outfile, int *winSizeByPop);

int get_arg_max(double *nums, int size);
int get_arg_min(double *nums, int size);
double get_min_btw_modes(double *x, double *y, int size);
double slope(double x0, double y0, double x1, double y1);

void doKDE(void *kdework);

double calculateWiggle(KDEResult *kdeResultByPop, int size = 20);

#endif
