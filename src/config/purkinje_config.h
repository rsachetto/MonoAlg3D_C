//
// Created by bergolho on 19/07/18.
//

#ifndef MONOALG3D_PURKINJE_CONFIG_H
#define MONOALG3D_PURKINJE_CONFIG_H

#include "../alg/grid/grid.h"
#include "config_common.h"
#include "../monodomain/constants.h"

#define SET_SPATIAL_PURKINJE(name) EXPORT_FN int name(struct time_info *time_info, struct config *config, struct grid *the_grid)
typedef SET_SPATIAL_PURKINJE(set_spatial_purkinje_fn);

void print_purkinje_config_values(struct config* s);

#endif //MONOALG3D_STIM_CONFIG_H
