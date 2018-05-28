//
// Created by sachetto on 14/10/17.
//

#ifndef MONOALG3D_CONFIG_COMMON_H
#define MONOALG3D_CONFIG_COMMON_H

#include <stdbool.h>
#include "../../hash/string_hash.h"

struct config_common {
    void *handle;
    char *function_name;
    bool function_name_was_set;
    char *library_file_path;
    bool library_file_path_was_set;
    struct string_hash *config;
};

void init_config_common_data(struct config_common *cm);


#endif //MONOALG3D_CONFIG_COMMON_H
