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

void advanceBar(Bar &bar, double inc);
void barInit(Bar &bar, double total, int totalTicks);
void finalize(Bar &bar);

#endif
