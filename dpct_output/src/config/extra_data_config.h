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

#define FREE_EXTRA_DATA(name) void *name (void *extra_data)
typedef FREE_EXTRA_DATA(free_extra_data_fn);

#define print_extra_data_config_values(s) LOG_COMMON_CONFIG("[extra_data]", s)
#define print_purkinje_extra_data_config_values(s) LOG_COMMON_CONFIG("[purkinje_extra_data]", s)

#define CALL_FREE_EXTRA_DATA(config, extra_data)                                                                      \
    do {                                                                                                              \
        if(config->end_function) {                                                                                    \
            ((free_extra_data_fn *)config->end_function)(extra_data);                                                 \
        }                                                                                                             \
        else {                                                                                                        \
            free(extra_data);                                                                                         \
        }                                                                                                             \
    } while(0)

#endif //MONOALG3D_EXTRA_DATA_CONFIG_H
