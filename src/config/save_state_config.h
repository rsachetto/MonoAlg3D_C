//
// Created by sachetto on 13/10/17.
//

#ifndef MONOALG3D_SAVE_STATE_CONFIG_H
#define MONOALG3D_SAVE_STATE_CONFIG_H

#include "../alg/grid/grid.h"
#include "../monodomain/constants.h"
#include "../monodomain/monodomain_solver.h"
#include "../ode_solver/ode_solver.h"
#include "config_common.h"

// Forward declaration
struct ode_solver;
struct monodomain_solver;

#define SAVE_STATE(name)                                                                                               \
     void name(struct time_info *time_info, struct config *config, struct grid *the_grid,                      \
                        struct monodomain_solver *the_monodomain_solver, struct ode_solver *the_ode_solver,            \
                        char *output_dir)

typedef SAVE_STATE(save_state_fn);

void print_save_state_config_values(struct config *s);

#endif // MONOALG3D_SAVE_CONFIG_H
