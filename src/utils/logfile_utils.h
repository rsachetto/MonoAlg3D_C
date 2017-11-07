//
// Created by sachetto on 18/10/17.
//

#ifndef MONOALG3D_LOGFILE_UTILS_H
#define MONOALG3D_LOGFILE_UTILS_H

#define LOG_LINE_SEPARATOR "======================================================================\n"


void print_to_stdout_and_file(char const *fmt, ...);
void open_logfile(const char *path);
void close_logfile();

#endif //MONOALG3D_LOGFILE_UTILS_H
