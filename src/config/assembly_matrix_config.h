//
// Created by sachetto on 13/10/17.
//

#ifndef MONOALG3D_ASSEMBLY_CONFIG_H
#define MONOALG3D_ASSEMBLY_CONFIG_H

#include "config_common.h"
#include "../monodomain/constants.h"
#include "../monodomain/monodomain_solver.h"

struct monodomain_solver;

//TODO: split into two different configs
#define ASSEMBLY_MATRIX(name)  void name(struct config *config, struct monodomain_solver *the_solver, struct grid *the_grid)
typedef ASSEMBLY_MATRIX(assembly_matrix_fn);

#define INIT_ASSEMBLY_MATRIX(name)  void name(struct config *config, struct monodomain_solver *the_solver, \
                                                       struct grid *the_grid,  real_cpu initial_v, real_cpu purkinje_initial_v)

typedef INIT_ASSEMBLY_MATRIX(set_pde_initial_condition_fn);

void print_assembly_matrix_config_values(struct config* s);

#endif //MONOALG3D_ASSEMBLY_CONFIG_H
