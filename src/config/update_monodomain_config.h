//
// Created by sachetto on 13/10/17.
//

#ifndef MONOALG3D_UPDATE_MONODOMAIN_CONFIG_H
#define MONOALG3D_UPDATE_MONODOMAIN_CONFIG_H

#include "config_common.h"
#include "../monodomain/constants.h"
#include "../monodomain/monodomain_solver.h"

struct update_monodomain_config;
struct ode_solver;
struct monodomain_solver;

#define UPDATE_MONODOMAIN(name) EXPORT_FN void name(struct update_monodomain_config *config, uint32_t initial_number_of_cells, struct monodomain_solver *the_solver, struct grid *the_grid, struct ode_solver *the_ode_solver)
typedef UPDATE_MONODOMAIN(update_monodomain_fn);

struct update_monodomain_config {
    struct config_common config_data;
    update_monodomain_fn *update_monodomain;
};

void init_update_monodomain_functions(struct update_monodomain_config *config);
struct update_monodomain_config* new_update_monodomain_config();
void free_update_monodomain_config(struct update_monodomain_config* s);
void print_update_monodomain_config_values(struct update_monodomain_config* s);


#endif //MONOALG3D_UPDATE_MONODOMAIN_CONFIG_H
