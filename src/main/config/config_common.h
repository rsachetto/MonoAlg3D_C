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
    char *library_file_path;
    bool configured;
    struct string_hash *config;
};

void init_config_common_data(struct config_common *cm);


#endif //MONOALG3D_CONFIG_COMMON_H
