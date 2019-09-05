//
// Created by sachetto on 13/10/17.
//

#ifndef MONOALG3D_SAVE_MESH_CONFIG_H
#define MONOALG3D_SAVE_MESH_CONFIG_H

#include "config_common.h"
#include "../monodomain/constants.h"
#include "../alg/grid/grid.h"

//Forward declaration
struct save_mesh_config;

#define SAVE_MESH(name) EXPORT_FN void name(struct config *config, struct grid *the_grid, int iteration_count, real_cpu current_t, real_cpu last_t, real dt, const char scalar_name)
typedef SAVE_MESH(save_mesh_fn);

#define INIT_SAVE_MESH(name) EXPORT_FN void name(struct config *config)
typedef INIT_SAVE_MESH(init_save_mesh_fn);

#define END_SAVE_MESH(name) EXPORT_FN void name(struct config *config)
typedef END_SAVE_MESH(end_save_mesh_fn);

#define CALL_INIT_SAVE_MESH(config)                                                                          \
    do {                                                                                                               \
        if(config->init_function) {                                                                                    \
            ((init_save_mesh_fn*) config->init_function)(config);                                     \
        }                                                                                                              \
    } while(0)

#define CALL_END_SAVE_MESH(config)                                                                                 \
    do {                                                                                                               \
        if(config->end_function) {                                                                                     \
            ((end_save_mesh_fn*)config->end_function)(config);                                              \
        }                                                                                                              \
    } while(0)

void print_save_mesh_config_values(struct config* s);

#endif //MONOALG3D_SAVE_MESH_CONFIG_H
