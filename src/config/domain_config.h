//
// Created by sachetto on 13/10/17.
//

#ifndef MONOALG3D_DOMAIN_CONFIG_H
#define MONOALG3D_DOMAIN_CONFIG_H

#include "../alg/grid/grid.h"
#include "config_common.h"
#include "../monodomain/constants.h"
#include "../common_types/common_types.h"

#define SET_SPATIAL_DOMAIN(name) int name(struct config *config, struct grid *the_grid)
typedef SET_SPATIAL_DOMAIN(set_spatial_domain_fn);

#define SET_CUSTOM_DATA_FOR_MESH(name) void name(struct cell_node *cell, real_cpu *custom_data)
typedef SET_CUSTOM_DATA_FOR_MESH(set_custom_data_for_mesh_fn);

#define print_domain_config_values(s) LOG_COMMON_CONFIG("[domain]", s)

#endif //MONOALG3D_STIM_CONFIG_H
