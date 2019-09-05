//
// Created by sachetto on 14/10/17.
//

#ifndef MONOALG3D_CONFIG_COMMON_H
#define MONOALG3D_CONFIG_COMMON_H

#include <stdbool.h>
#include "../common_types/common_types.h"

struct config {
    void *handle;

    char *main_function_name;
    bool main_function_name_was_set;

    char *init_function_name;
    bool init_function_name_was_set;

    char *end_function_name;
    bool end_function_name_was_set;

    char *library_file_path;
    bool library_file_path_was_set;

    struct string_hash_entry *config_data;

    void *persistent_data;

    void *main_function;
    void *init_function;
    void *end_function;
};

struct config *alloc_and_init_config_data();
void free_config_data(struct config *cm);
void init_config_functions(struct config *config, char *default_lib, char *config_type);

#endif //MONOALG3D_CONFIG_COMMON_H
