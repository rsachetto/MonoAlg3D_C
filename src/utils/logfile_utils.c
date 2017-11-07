//
// Created by sachetto on 18/10/17.
//

#include <stdarg.h>

#include "logfile_utils.h"
#include <stdio.h>

static FILE *logfile = NULL;

void print_to_stdout_and_file(char const *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vprintf(fmt, ap);
    va_end(ap);
    va_start(ap, fmt);
    if(logfile)
        vfprintf(logfile, fmt, ap);
    va_end(ap);
}

void open_logfile(const char *path) {
    logfile = fopen(path, "w");

    if(logfile == NULL) {
        fprintf(stderr, "Error opening %s, printing output only in the sdtout (Terminal)\n", path);
    }
    else {
        printf("Log will be saved in %s\n", path);
    }
}

void close_logfile() {
    if(logfile) fclose(logfile);
}
