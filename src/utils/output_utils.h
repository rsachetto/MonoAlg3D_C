//
// Created by sachetto on 03/10/17.
//

#ifndef MONOALG3D_OUTPUT_INFO_H
#define MONOALG3D_OUTPUT_INFO_H

#include <stdbool.h>

struct output_utils {
    int print_rate;

    //TODO: we should be very carefull when allocating this strings
    char *output_dir_name;
};

bool dir_exists(const char *path);

#endif //MONOALG3D_OUTPUT_INFO_H
