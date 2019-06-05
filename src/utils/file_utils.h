//
// Created by sachetto on 18/10/17.
//

#ifndef MONOALG3D_FILE_UTILS_H
#define MONOALG3D_FILE_UTILS_H

#define LOG_LINE_SEPARATOR "======================================================================\n"


#include <stdbool.h>
#include <stddef.h>
#include <sys/mman.h>

#include "../common_types/common_types.h"

bool no_stdout;

void print_to_stdout_and_file(char const *fmt, ...);
void print_to_stderr_and_file_and_exit(char const *fmt, ...);
void open_logfile(const char *path);
void close_logfile();
int cp_file(const char *to, const char *from);
char * read_entire_file(const char *filename, long *size);
char *read_entire_file_with_mmap(const char *filename, size_t *size);
string_array list_files_from_dir(const char *dir, const char *prefix);
string_array list_files_from_dir_sorted(const char *dir, const char *prefix);
string_array read_lines(const char *filename);
bool dir_exists(const char *path);
void fixpath(char *path);
void create_dir(char *out_dir);
int remove_directory(const char *path);
size_t base64_decode(unsigned char* out, const char *src, size_t len, size_t *bytes_read);
char * get_dir_from_path(const char * path);

#endif //MONOALG3D_FILE_UTILS_H
