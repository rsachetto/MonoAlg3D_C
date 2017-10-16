//
// Created by sachetto on 13/10/17.
//

#ifndef MONOALG3D_EXTRA_DATA_CONFIG_H
#define MONOALG3D_EXTRA_DATA_CONFIG_H

#include "../../alg/grid/grid.h"
#include "../../hash/string_hash.h"
#include "config_common.h"

typedef void * (*set_extra_data_fn_pt)(struct grid *, struct string_hash *, size_t *);

struct extra_data_config {
    struct config_common config_data;
    set_extra_data_fn_pt set_extra_data_fn;
};

void init_extra_data_functions(struct extra_data_config *config);
struct extra_data_config* new_extra_data_config();
void print_extra_data_config_values(struct extra_data_config* s);
void free_extra_data_config(struct extra_data_config* s);

#endif //MONOALG3D_EXTRA_DATA_CONFIG_H
