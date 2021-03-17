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

void log_to_stdout_and_file(const char *prefix, char const *fmt, ...) {

    char *fmt_plus_prefix;

    if(prefix) {
        fmt_plus_prefix = malloc(strlen(fmt) + strlen(prefix) + 2);
        sprintf(fmt_plus_prefix, "%s %s", prefix, fmt);
    }
    else {
        fmt_plus_prefix = (char*) fmt;
    }

    va_list ap;

    if (!L.quiet) {
        va_start(ap, fmt);
        vprintf(fmt_plus_prefix, ap);
        fflush(stdout);
        va_end(ap);
    }

    if (L.fp) {
        va_start(ap, fmt);
        vfprintf(L.fp, fmt_plus_prefix, ap);
        fflush(L.fp);
        va_end(ap);
    }

    if(prefix) {
        free(fmt_plus_prefix);
    }
}

void log_to_stderr_and_file(bool exit_program, const char *prefix, char const *fmt, ...) {

    char *fmt_plus_prefix;

    if(prefix) {
        fmt_plus_prefix = malloc(strlen(fmt) + strlen(prefix) + 2);
        sprintf(fmt_plus_prefix, "%s %s", prefix, fmt);
    }
    else {
        fmt_plus_prefix = (char*) fmt;
    }

    va_list ap;
    va_start(ap, fmt);
    vprintf(fmt_plus_prefix, ap);
    fflush(stderr);
    va_end(ap);

    va_start(ap, fmt);
    if (L.fp) {
        vfprintf(L.fp, fmt_plus_prefix, ap);
        fflush(L.fp);
    }
    va_end(ap);

    if(prefix) {
        free(fmt_plus_prefix);
    }

    if(exit_program)
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
