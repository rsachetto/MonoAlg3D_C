//
// Created by sachetto on 18/10/17.
//

#ifndef MONOALG3D_FILE_UTILS_H
#define MONOALG3D_FILE_UTILS_H

#define LOG_LINE_SEPARATOR "======================================================================\n"


#include <stdbool.h>
#include "../common_types/common_types.h"

bool no_stdout;

void print_to_stdout_and_file(char const *fmt, ...);
void print_to_stderr_and_file_and_exit(char const *fmt, ...);
void open_logfile(const char *path);
void close_logfile();
int cp_file(const char *to, const char *from);
char * read_entire_file(char *filename, long *size);
string_array list_files_from_dir(const char *dir, const char *prefix);
string_array read_lines(const char *filename);
bool dir_exists(const char *path);
void fixpath(char *path);
void create_dir(char *out_dir);
int remove_directory(const char *path);

#endif //MONOALG3D_FILE_UTILS_H
