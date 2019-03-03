//
// Created by sachetto on 18/10/17.
//

#ifndef MONOALG3D_FILE_UTILS_H
#define MONOALG3D_FILE_UTILS_H

#define LOG_LINE_SEPARATOR "======================================================================\n"


#include <stdbool.h>

bool no_stdout;

void print_to_stdout_and_file(char const *fmt, ...);
void print_to_stderr_and_file_and_exit(char const *fmt, ...);
void open_logfile(const char *path);
void close_logfile();
int cp_file(const char *to, const char *from);
char * read_entire_file(char *filename, long *size);
char ** list_files_from_dir(const char *dir, const char *prefix);
char **read_lines(const char *filename);
bool dir_exists(const char *path);
void create_dir(const char *out_dir);
int remove_directory(const char *path);

#endif //MONOALG3D_LOGFILE_UTILS_H
