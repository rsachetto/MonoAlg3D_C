//
// Created by sachetto on 07/11/2019.
//

#ifndef MONOALG3D_C_LOGGER_H
#define MONOALG3D_C_LOGGER_H

#include <stdbool.h>
#include <stdio.h>

#define LOG_LINE_SEPARATOR "======================================================================\n"

struct logt {
    FILE *fp;
    bool quiet;
};

#define INFO_PREFIX "[INFO]"
#define WARN_PREFIX "[WARN]"
#define ERROR_PREFIX "[ERR]"
#define NO_PREFIX NULL

#define log_msg(fmt, ...) log_to_stdout_and_file(NO_PREFIX, (fmt), ##__VA_ARGS__)
#define log_info(fmt, ...) log_to_stdout_and_file(NO_PREFIX, (fmt), ##__VA_ARGS__)
#define log_warn(fmt, ...) log_to_stdout_and_file(WARN_PREFIX, (fmt), ##__VA_ARGS__)
#define log_error(fmt, ...) log_to_stderr_and_file(false, ERROR_PREFIX, (fmt), ##__VA_ARGS__)
#define log_error_and_exit(fmt, ...) log_to_stderr_and_file(true, ERROR_PREFIX, (fmt), ##__VA_ARGS__)

#ifdef __cplusplus
extern "C" {
#endif

void set_no_stdout(bool val);
void open_logfile(const char *path);
void close_logfile();
void log_to_stdout_and_file(const char *prefix, char const *fmt, ...);
void log_to_stderr_and_file(bool exit, const char *prefix, char const *fmt, ...);

#ifdef __cplusplus
}
#endif

#endif // MONOALG3D_C_LOGGER_H
