#ifndef __STOPWATCH_H
#define __STOPWATCH_H

#include <stdbool.h>

#ifdef _MSC_VER
#include <Windows.h>
#include <stdint.h> // portable: uint64_t   MSVC: __int64 
#else
#include <sys/time.h>
#endif

/* simple stopwatch class */
struct stop_watch {
    struct timeval tv;
    bool running;

};

struct stop_watch* new_stop_watch();
void init_stop_watch(struct stop_watch *sw);
void start_stop_watch(struct stop_watch *sw); /* start the stopwatch */
long stop_stop_watch(struct stop_watch *sw); /* stop the stopwatch and get the value in usecs */

#endif /* __STOPWATCH_H */
