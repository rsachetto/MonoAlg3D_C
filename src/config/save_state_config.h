//
// Created by sachetto on 13/10/17.
//

#ifndef MONOALG3D_SAVE_STATE_CONFIG_H
#define MONOALG3D_SAVE_STATE_CONFIG_H

#include "config_common.h"
#include "../monodomain/constants.h"
#include "../monodomain/ode_solver.h"
#include "../monodomain/monodomain_solver.h"
#include "../alg/grid/grid.h"

//Forward declaration
struct ode_solver;
struct monodomain_solver;

#define SAVE_STATE(name) EXPORT_FN void name(char *output_dir, \
                                             struct config *config, \
                                             struct grid *the_grid, \
                                             struct monodomain_solver *the_monodomain_solver,  \
                                             struct ode_solver *the_ode_solver)

typedef SAVE_STATE(save_state_fn);


void print_save_state_config_values(struct config* s);


#endif //MONOALG3D_SAVE_CONFIG_H
