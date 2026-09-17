#include "garlic-pbar.h"
#include <mutex>

//std::mutex rather than a pthread_mutex_t with a static initialiser: the
//initialiser macro is the only thing this file needed pthreads for, and a
//std::mutex is default-constructed with no macro and no platform variant.
static std::mutex mutex_progress;

static bool BAR_ON = true;

void setProgressEnabled(bool on) { BAR_ON = on; }
bool progressEnabled() { return BAR_ON; }

void advanceBar(Bar &bar, double inc)
{
    if (!BAR_ON) { std::lock_guard<std::mutex> lock(mutex_progress); bar.current += inc; return; }
    std::unique_lock<std::mutex> lock(mutex_progress);
    bar.current += inc;
    if (bar.current / bar.total >= double(bar.currentTick) / double(bar.totalTicks))
    {
        bar.currentTick++;
        for (int i = 0; i < 3; i++) cerr << '\b';
        if(int((bar.current / bar.total) * 100) < 10) cerr << " ";
        cerr /*<< setprecision(4)*/ << int((bar.current / bar.total) * 100) << '%';
        cerr.flush();
    }
    lock.unlock();
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
    if (!BAR_ON) return;
    for (int i = 0; i < 3; i++) cerr << '\b';
    cerr << "100%" << endl;
}
