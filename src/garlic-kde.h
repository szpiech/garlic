#ifndef __GARLIC_KDE_H__
#define __GARLIC_KDE_H__

#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <cstring>
#include <limits>
#include "garlic-math.h"
#include "garlic-data.h"
#include "garlic-errlog.h"

struct KDEResult
{
    int size;
    vector<double> x;
    vector<double> y;
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

//Which modes a LOD density has, by the sign anchor.  E[LOD] per site is
//-KL(N||A) < 0 for a non-autozygous window and +KL(A||N) > 0 for an
//autozygous one, whatever the depth or error rate, so the two modes lie on
//opposite sides of zero: the non-autozygous mode is the dominant peak at
//LOD < 0, the autozygous mode the dominant peak at LOD > 0, and the cutoff
//is the density minimum between them.  A density with no peak on one side
//has no cutoff, and says so instead of inventing one.
enum LODModeVerdict
{
    LOD_MODES_BIMODAL = 0,           //both modes found; cutoff is the valley between them
    LOD_MODES_UNIMODAL_NULL,         //a non-autozygous mode only: no resolvable autozygous excess
    LOD_MODES_UNIMODAL_AUTOZYGOUS,   //an autozygous mode only: the non-autozygous state is not represented
    LOD_MODES_NONE                   //no mode passed the filters (degenerate density)
};
struct LODModes
{
    LODModeVerdict verdict;
    double modeN;        //location of the non-autozygous mode (NaN if absent)
    double modeA;        //location of the autozygous mode (NaN if absent)
    double cutoff;       //valley between them (NaN unless bimodal)
    double valleyHeight; //density at the cutoff, as a fraction of the shorter mode's density (NaN unless bimodal)
    int    nCandidates;  //peaks that passed the height and persistence filters
    LODModes();
};
LODModes identifyLODModes(const double *x, const double *y, int size);
const char *lodModeVerdictName(LODModeVerdict v);
//The relative height below which a local maximum is density noise or a tail
//bump rather than a mode, and the number of consecutive windowed-argmax
//positions a maximum must hold to count as a mode (half the smoothing span).
void setModeFilters(double minFrac, int minPersist);
//Cutoff or throw: the valley of a bimodal density, an exception otherwise.
//Kept for callers that only need a number and treat its absence as
//"no second mode" (the sex-chromosome report, the strata).
double get_min_btw_modes(double *x, double *y, int size, int wsize);
double slope(double x0, double y0, double x1, double y1);

double calculateWiggle(KDEResult *kdeResult, int size = 20);

/*
KDEWinsizeReport *initKDEWinsizeReport();
void releaseKDEWinsizeReport(KDEWinsizeReport *winsizeReport);
*/

#endif
