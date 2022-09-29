//
// Created by sachetto on 13/10/17.
//

#ifndef MONOALG3D_STIM_CONFIG_H
#define MONOALG3D_STIM_CONFIG_H

#include "../alg/grid/grid.h"
#include "../monodomain/constants.h"
#include "config_common.h"

#define SET_STIM_VALUE(i, stim_value) ((real *)(config->persistent_data))[i] = stim_value

#define ALLOCATE_STIMS()                                                                                                                                       \
    do {                                                                                                                                                       \
        if(the_grid->adaptive) {                                                                                                                               \
            if(config->persistent_data) {                                                                                                                      \
                free(config->persistent_data);                                                                                                                 \
            }                                                                                                                                                  \
            config->persistent_data = (real *)malloc(n_active * sizeof(real));                                                                                 \
        } else {                                                                                                                                               \
            if(!config->persistent_data) {                                                                                                                     \
                config->persistent_data = (real *)malloc(n_active * sizeof(real));                                                                             \
            }                                                                                                                                                  \
        }                                                                                                                                                      \
    } while(0)
#define SET_SPATIAL_STIM(name) void name(struct time_info *time_info, struct config *config, struct grid *the_grid, bool is_purkinje)
typedef SET_SPATIAL_STIM(set_spatial_stim_fn);

#define print_stim_config_values(s) LOG_COMMON_CONFIG("[stim]", (s))

#endif // MONOALG3D_STIM_CONFIG_H
