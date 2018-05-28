//
// Created by sachetto on 13/10/17.
//

#ifndef MONOALG3D_LINEAR_SOLVER_CONFIG_H
#define MONOALG3D_LINEAR_SOLVER_CONFIG_H

#include "config_common.h"
#include "../monodomain/constants.h"
#include "../alg/grid/grid.h"

struct linear_system_solver_config; //Forward declaration

#define SOLVE_LINEAR_SYSTEM(name) EXPORT_FN void name(struct linear_system_solver_config *config, struct grid *the_grid, uint32_t *number_of_iterations, double *error)
typedef SOLVE_LINEAR_SYSTEM(linear_system_solver_fn);

struct linear_system_solver_config {
    struct config_common config_data;
    linear_system_solver_fn *solve_linear_system;
};

void init_linear_system_solver_functions(struct linear_system_solver_config *config);
struct linear_system_solver_config* new_linear_system_solver_config();
void print_linear_system_solver_config_values(struct linear_system_solver_config* s);
void free_linear_system_solver_config(struct linear_system_solver_config* s);

#endif //MONOALG3D_ASSEMBLY_CONFIG_H
