//
// Created by sachetto on 13/10/17.
//

#ifndef MONOALG3D_RESTORE_STATE_CONFIG_H
#define MONOALG3D_RESTORE_STATE_CONFIG_H

#include "config_common.h"
#include "../monodomain/constants.h"
#include "../monodomain/monodomain_solver.h"
#include "../ode_solver/ode_solver.h"
#include "../alg/grid/grid.h"

//Forward declaration
struct monodomain_solver;
struct ode_solver;

#define RESTORE_STATE(name)  bool name(struct time_info *time_info,                      \
                                                struct config *config,                            \
                                                struct grid *the_grid,                            \
                                                struct monodomain_solver *the_monodomain_solver,  \
                                                struct ode_solver *the_ode_solver, char* input_dir)
typedef RESTORE_STATE(restore_state_fn);


void print_restore_state_config_values(struct config* s);


#endif //MONOALG3D_SAVE_CONFIG_H
