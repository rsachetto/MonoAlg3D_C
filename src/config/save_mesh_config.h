//
// Created by sachetto on 13/10/17.
//

#ifndef MONOALG3D_SAVE_MESH_CONFIG_H
#define MONOALG3D_SAVE_MESH_CONFIG_H

#include "../alg/grid/grid.h"
#include "../ode_solver/ode_solver.h"
#include "../monodomain/constants.h"
#include "config_common.h"

#define SAVE_MESH(name) void name(struct time_info *time_info, struct config *config, struct grid *the_grid, struct ode_solver *ode_solver)
typedef SAVE_MESH(save_mesh_fn);

#define INIT_SAVE_MESH(name)  void name(struct config *config)
typedef INIT_SAVE_MESH(init_save_mesh_fn);

#define END_SAVE_MESH(name)  void name(struct config *config, struct grid *the_grid)
typedef END_SAVE_MESH(end_save_mesh_fn);

#define CALL_INIT_SAVE_MESH(config)                                                                                    \
    do {                                                                                                               \
        if(config && config->init_function) {                                                                          \
            ((init_save_mesh_fn *)config->init_function)(config);                                                      \
        }                                                                                                              \
    } while(0)

#define CALL_END_SAVE_MESH(config,grid)                                                                                     \
    do {                                                                                                               \
        if(config && config->end_function && grid) {                                                                           \
            ((end_save_mesh_fn *)config->end_function)(config,grid);                                                        \
        }                                                                                                              \
    } while(0)

void print_save_mesh_config_values(struct config *s);

#endif // MONOALG3D_SAVE_MESH_CONFIG_H
