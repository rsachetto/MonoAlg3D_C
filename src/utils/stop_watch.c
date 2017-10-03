/*
 * StopWatch.cpp
 *
 *  Created on: 29/07/2011
 *      Author: sachetto
 */

#include "stop_watch.h"

#include <stdlib.h>
#include <stdio.h>

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
