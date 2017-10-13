//
// Created by sachetto on 03/10/17.
//

#ifndef MONOALG3D_OUTPUT_INFO_H
#define MONOALG3D_OUTPUT_INFO_H

#include <stdbool.h>
#include "string/sds.h"
#include "config_parser.h"

struct output_utils {
    int print_rate;
    sds output_dir_name;
};

struct output_utils *new_output_utils();
void free_output_utils(struct output_utils* info);
bool dir_exists(const char *path);
void configure_output_from_options(struct output_utils *output_utils,
                                   struct user_options *options);

#endif //MONOALG3D_OUTPUT_INFO_H
