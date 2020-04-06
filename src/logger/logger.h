//
// Created by sachetto on 07/11/2019.
//

#ifndef MONOALG3D_C_LOGGER_H
#define MONOALG3D_C_LOGGER_H

#include <stdbool.h>
#include <stdio.h>

struct logt {
    FILE *fp;
    bool quiet;
};

void set_no_stdout(bool val);
void log_to_stdout_and_file(char const *fmt, ...);
void log_to_stderr_and_file(char const *fmt, ...);
void log_to_stderr_and_file_and_exit(char const *fmt, ...);
void open_logfile(const char *path);
void close_logfile();

#endif //MONOALG3D_C_LOGGER_H
