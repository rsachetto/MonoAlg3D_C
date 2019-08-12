//
// Created by sachetto on 13/10/17.
//

#ifndef MONOALG3D_SAVE_MESH_CONFIG_H
#define MONOALG3D_SAVE_MESH_CONFIG_H

#include "config_common.h"
#include "../monodomain/constants.h"
#include "../alg/grid/grid.h"

#define SAVE_MESH(name) EXPORT_FN void name(struct config *config, struct grid *the_grid, int iteration_count, real_cpu current_t, real_cpu last_t, real dt)
typedef SAVE_MESH(save_mesh_fn);

void print_save_mesh_config_values(struct config* s);

#endif //MONOALG3D_SAVE_MESH_CONFIG_H
