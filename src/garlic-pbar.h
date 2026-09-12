#ifndef __PBAR_H__
#define __PBAR_H__
#include <iostream>
#include <iomanip>

using namespace std;

struct Bar
{
    double total;
    double current;
    int totalTicks;
    int currentTick;
};

//The bar emits backspaces, so it is only meaningful on a terminal.  Redirected
//stderr used to collect literal "\b\b\b 0%\b\b\b100%" for every chromosome.
void setProgressEnabled(bool on);
bool progressEnabled();

void advanceBar(Bar &bar, double inc);
void barInit(Bar &bar, double total, int totalTicks);
void finalize(Bar &bar);

#endif
