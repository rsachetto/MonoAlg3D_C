//
// Created by sachetto on 13/10/17.
//

#ifndef MONOALG3D_EXTRA_DATA_CONFIG_H
#define MONOALG3D_EXTRA_DATA_CONFIG_H

#include "../alg/grid/grid.h"
#include "config_common.h"
#include "../monodomain/constants.h"
#include "../common_types/common_types.h"

#define SET_EXTRA_DATA(name) void *name (struct time_info *time_info, struct config *config, struct grid * the_grid, size_t *extra_data_size)
typedef SET_EXTRA_DATA(set_extra_data_fn);

void print_extra_data_config_values(struct config* s);

#endif //MONOALG3D_EXTRA_DATA_CONFIG_H
