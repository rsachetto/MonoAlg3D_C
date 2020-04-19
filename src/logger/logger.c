//
// Created by sachetto on 07/11/2019.
//

#include "logger.h"
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

static struct logt L;

void set_no_stdout(bool val) {
    L.quiet = val;
}

void log_to_stdout_and_file(char const *fmt, ...) {

    va_list ap;

    if (!L.quiet) {
        va_start(ap, fmt);
        vprintf(fmt, ap);
        fflush(stdout);
        va_end(ap);
    }

    if (L.fp) {
        va_start(ap, fmt);
        vfprintf(L.fp, fmt, ap);
        fflush(L.fp);
        va_end(ap);
    }
}

void log_to_stderr_and_file(char const *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vprintf(fmt, ap);
    fflush(stderr);
    va_end(ap);
    va_start(ap, fmt);
    if (L.fp) {
        vfprintf(L.fp, fmt, ap);
        fflush(L.fp);
    }
    va_end(ap);
}

void log_to_stderr_and_file_and_exit(char const *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    log_to_stderr_and_file(fmt, ap);
    va_end(ap);
    exit(EXIT_FAILURE);
}

void open_logfile(const char *path) {
    L.fp = fopen(path, "w");

    if (L.fp == NULL) {
        fprintf(stderr, "Error opening %s, printing output only in the sdtout (Terminal)\n", path);
    } else {
        printf("Log will be saved in %s\n", path);
    }
}

void close_logfile() {
    if (L.fp) fclose(L.fp);
}
