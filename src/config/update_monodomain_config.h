//
// Created by sachetto on 13/10/17.
//

#ifndef MONOALG3D_UPDATE_MONODOMAIN_CONFIG_H
#define MONOALG3D_UPDATE_MONODOMAIN_CONFIG_H

#include "config_common.h"
#include "../monodomain/constants.h"
#include "../monodomain/monodomain_solver.h"

struct ode_solver;
struct monodomain_solver;

#define UPDATE_MONODOMAIN(name) EXPORT_FN void name(struct config *config, uint32_t initial_number_of_cells, struct monodomain_solver *the_solver, const uint32_t num_active_cells, struct cell_node **active_cells, struct ode_solver *the_ode_solver)
typedef UPDATE_MONODOMAIN(update_monodomain_fn);

void print_update_monodomain_config_values(struct config* s);


#endif //MONOALG3D_UPDATE_MONODOMAIN_CONFIG_H
