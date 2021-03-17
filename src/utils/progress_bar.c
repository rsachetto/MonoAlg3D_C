//
// Created by sachetto on 11/03/2021.
//

#include <stdio.h>

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void print_progress(double percentage) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\rProgress - %3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}