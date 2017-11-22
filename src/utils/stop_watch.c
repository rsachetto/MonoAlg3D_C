/*
 * StopWatch.cpp
 *
 *  Created on: 29/07/2011
 *      Author: sachetto
 */

#include "stop_watch.h"

#include <stdlib.h>
#include <stdio.h>

#ifdef _MSC_VER
int gettimeofday(struct timeval * tp, struct timezone * tzp)
{
	// Note: some broken versions only have 8 trailing zero's, the correct epoch has 9 trailing zero's
	// This magic number is the number of 100 nanosecond intervals since January 1, 1601 (UTC)
	// until 00:00:00 January 1, 1970 
	static const uint64_t EPOCH = ((uint64_t)116444736000000000ULL);

	SYSTEMTIME  system_time;
	FILETIME    file_time;
	uint64_t    time;

	GetSystemTime(&system_time);
	SystemTimeToFileTime(&system_time, &file_time);
	time = ((uint64_t)file_time.dwLowDateTime);
	time += ((uint64_t)file_time.dwHighDateTime) << 32;

	tp->tv_sec = (long)((time - EPOCH) / 10000000L);
	tp->tv_usec = (long)(system_time.wMilliseconds * 1000);
	return 0;
}
#endif

struct stop_watch* new_stop_watch() {
	return (struct stop_watch*) malloc(sizeof(struct stop_watch));
}

void init_stop_watch(struct stop_watch *sw) {
	sw->running = false;
}

void start_stop_watch(struct stop_watch *sw) {
	if (gettimeofday(&(sw->tv), NULL) < 0) {
		perror("gettimeofday");
		return;
	}
	sw->running = true;
}

long stop_stop_watch(struct stop_watch *sw) {

    struct timeval tv_stop;

	if (!sw->running) {
		fprintf(stderr,"Stopwatch not running.");
		return -1;
	}

	if (gettimeofday(&tv_stop, NULL) < 0) {
		perror("gettimeofday");
		return -1;
	}

	sw->running = false;
	return ((tv_stop.tv_sec - sw->tv.tv_sec) * 1000000
				+ tv_stop.tv_usec - sw->tv.tv_usec);
}
