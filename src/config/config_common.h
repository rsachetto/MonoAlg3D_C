//
// Created by sachetto on 14/10/17.
//

#ifndef MONOALG3D_CONFIG_COMMON_H
#define MONOALG3D_CONFIG_COMMON_H

#include <stdbool.h>
#include "../common_types/common_types.h"

struct config_common {
    void *handle;
    char *function_name;
    bool function_name_was_set;
    char *library_file_path;
    bool library_file_path_was_set;
    struct string_hash_entry *config;

};

struct generic_config {
    //This is used to cast any config that has only a config_common and a
    //function pointer
    struct config_common config_data;
};

void init_config_common_data(struct config_common *cm);
void free_config_common_data(struct config_common *cm);


#endif //MONOALG3D_CONFIG_COMMON_H
