//
// Created by sachetto on 18/10/17.
//

#ifndef MONOALG3D_FILE_UTILS_H
#define MONOALG3D_FILE_UTILS_H

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <sys/mman.h>

#include "../common_types/common_types.h"

#define FILE_EXTENSION_EQUALS(file_ext__, ext__) (strcmp((file_ext__), (ext__)) == 0)
#define FILE_HAS_EXTENSION(file__, ext__) FILE_EXTENSION_EQUALS((file__).file_extension, ext__)

#define FILE_EXTENSION_EQUALS_PREFIX(file_ext__, ext__) (strncmp((file_ext__), (ext__), strlen((ext__))) == 0)
#define FILE_HAS_EXTENSION_PREFIX(file__, ext__) FILE_EXTENSION_EQUALS_PREFIX((file__).file_extension, ext__)
#define ENDS_WITH_SLASH(path__) (path__[strlen(path__) - 1] == '/')

struct path_information {
    bool exists, is_file, is_dir;
    char *file_extension;
    char *filename_without_extension;
    char *dir_name;
};

#ifdef __cplusplus
extern "C" {
#endif

float *read_timesteps_from_case_file(sds case_file_path);
char *get_current_directory();
const char *get_filename_ext(const char *filename);
int cp_file(const char *to, const char *from);
char *read_entire_file(const char *filename, size_t *size);
char *read_entire_file_with_mmap(const char *filename, size_t *size);
string_array list_files_from_dir(const char *dir, const char *prefix, const char *extension, string_array ignore_extensions, bool sort);
string_array read_lines(const char *filename);
bool dir_exists(const char *path);
bool file_exists(const char *path);

int get_step_from_filename(char *filename);
void get_path_information(const char *path, struct path_information *input_info);
void free_path_information(struct path_information *input_info);

char *get_filename_without_ext(const char *filename);
void create_dir(char *out_dir);
int remove_directory(const char *path);
size_t base64_decode(unsigned char *out, const char *src, size_t len, size_t *bytes_read);
char *get_dir_from_path(const char *path);
char *get_file_from_path(const char *path);
bool check_simulation_completed(char *simulation_dir);

real_cpu **read_octave_mat_file_to_array(FILE *matrix_file, uint64_t *num_lines, uint64_t *nnz);
real_cpu *read_octave_vector_file_to_array(FILE *vec_file, uint64_t *num_lines);
FILE *open_file_or_exit(char *filename, char *mode);
char *get_timestamped_dir_name(char *dir_name);

#ifdef __cplusplus
}
#endif

#endif // MONOALG3D_FILE_UTILS_H
