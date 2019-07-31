//
// Created by sachetto on 13/10/17.
//

#ifndef MONOALG3D_LINEAR_SOLVER_CONFIG_H
#define MONOALG3D_LINEAR_SOLVER_CONFIG_H

#include "../alg/grid/grid.h"
#include "../monodomain/constants.h"
#include "config_common.h"

struct linear_system_solver_config; // Forward declaration

#define SOLVE_LINEAR_SYSTEM(name)                                                                                      \
    EXPORT_FN void name(struct linear_system_solver_config *config, struct grid *the_grid,                             \
                        uint32_t *number_of_iterations, real_cpu *error)
typedef SOLVE_LINEAR_SYSTEM(linear_system_solver_fn);

#define INIT_LINEAR_SYSTEM(name) EXPORT_FN void name(struct linear_system_solver_config *config, struct grid *the_grid)
typedef INIT_LINEAR_SYSTEM(init_linear_system_solver_fn);

#define END_LINEAR_SYSTEM(name) EXPORT_FN void name(struct linear_system_solver_config *config)
typedef END_LINEAR_SYSTEM(end_linear_system_solver_fn);

struct linear_system_solver_config {
    struct config_common config_data;
    linear_system_solver_fn *solve_linear_system;
    init_linear_system_solver_fn *init_linear_system;
    end_linear_system_solver_fn *end_linear_system;
};

#define CALL_INIT_LINEAR_SYSTEM(config, grid)                                                                          \
    do {                                                                                                               \
        if(config->init_linear_system) {                                                                               \
            config->init_linear_system(config, grid);                                                                  \
        }                                                                                                              \
    } while(0)

#define CALL_END_LINEAR_SYSTEM(config)                                                                                 \
    do {                                                                                                               \
        if(config->end_linear_system) {                                                                                \
            config->end_linear_system(config);                                                                               \
        }                                                                                                              \
    } while(0)

void init_linear_system_solver_functions(struct linear_system_solver_config *config);
struct linear_system_solver_config *new_linear_system_solver_config();
void print_linear_system_solver_config_values(struct linear_system_solver_config *s);
void free_linear_system_solver_config(struct linear_system_solver_config *s);

#endif // MONOALG3D_ASSEMBLY_CONFIG_H
