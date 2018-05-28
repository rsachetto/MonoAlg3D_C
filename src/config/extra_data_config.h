//
// Created by sachetto on 13/10/17.
//

#ifndef MONOALG3D_EXTRA_DATA_CONFIG_H
#define MONOALG3D_EXTRA_DATA_CONFIG_H

#include "../../alg/grid/grid.h"
#include "../../hash/string_hash.h"
#include "config_common.h"
#include "../constants.h"

#define SET_EXTRA_DATA(name) EXPORT_FN void *name (struct grid * the_grid, struct string_hash *config, size_t *extra_data_size)
typedef SET_EXTRA_DATA(set_extra_data_fn);

struct extra_data_config {
    struct config_common config_data;
    set_extra_data_fn *set_extra_data;
};

void init_extra_data_functions(struct extra_data_config *config);
struct extra_data_config* new_extra_data_config();
void print_extra_data_config_values(struct extra_data_config* s);
void free_extra_data_config(struct extra_data_config* s);

#endif //MONOALG3D_EXTRA_DATA_CONFIG_H
