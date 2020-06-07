//
// Created by sachetto on 13/10/17.
//

#ifndef MONOALG3D_LINEAR_SOLVER_CONFIG_H
#define MONOALG3D_LINEAR_SOLVER_CONFIG_H

#include "../alg/grid/grid.h"
#include "../monodomain/constants.h"
#include "config_common.h"

#define SOLVE_LINEAR_SYSTEM(name)                                                                                      \
              void name(struct time_info *time_info, struct config *config, struct grid *the_grid,                     \
                        uint32_t num_active_cells, struct cell_node **active_cells,                                    \
                        uint32_t *number_of_iterations, real_cpu *error)
typedef SOLVE_LINEAR_SYSTEM(linear_system_solver_fn);

#define INIT_LINEAR_SYSTEM(name)  void name(struct config *config, struct grid *the_grid)
typedef INIT_LINEAR_SYSTEM(init_linear_system_solver_fn);

#define END_LINEAR_SYSTEM(name)  void name(struct config *config)
typedef END_LINEAR_SYSTEM(end_linear_system_solver_fn);

#define CALL_INIT_LINEAR_SYSTEM(config, grid)                                                                          \
    do {                                                                                                               \
        if(config && config->init_function) {                                                                          \
            ((init_linear_system_solver_fn *)config->init_function)(config, grid);                                     \
        }                                                                                                              \
    } while(0)

#define CALL_END_LINEAR_SYSTEM(config)                                                                                 \
    do {                                                                                                               \
        if(config && config->end_function) {                                                                           \
            ((end_linear_system_solver_fn *)config->end_function)(config);                                             \
        }                                                                                                              \
    } while(0)

void print_linear_system_solver_config_values(struct config *s);

#endif // MONOALG3D_ASSEMBLY_CONFIG_H
