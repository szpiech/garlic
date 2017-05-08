#include "garlic-pbar.h"
#include <pthread.h>

pthread_mutex_t mutex_progress = PTHREAD_MUTEX_INITIALIZER;

void advanceBar(Bar &bar, double inc)
{
    pthread_mutex_lock(&mutex_progress);
    bar.current += inc;
    if (bar.current / bar.total >= double(bar.currentTick) / double(bar.totalTicks))
    {
        bar.currentTick++;
        for (int i = 0; i < 3; i++) cerr << '\b';
        if(int((bar.current / bar.total) * 100) < 10) cerr << " ";
        cerr /*<< setprecision(4)*/ << int((bar.current / bar.total) * 100) << '%';
        cerr.flush();
    }
    pthread_mutex_unlock(&mutex_progress);
    return;
}

void barInit(Bar &bar, double total, int totalTicks)
{
    bar.total = total;
    bar.current = 0;
    bar.totalTicks = totalTicks;
    bar.currentTick = 0;
    return;
}

void finalize(Bar &bar){
    for (int i = 0; i < 3; i++) cerr << '\b';
    cerr << "100%" << endl;
}
