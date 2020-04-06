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

void print_domain_config_values(struct config* s);

#endif //MONOALG3D_STIM_CONFIG_H
