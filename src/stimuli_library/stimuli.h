//
// Created by sachetto on 13/10/17.
//

#ifndef MONOALG3D_STIMS_H
#define MONOALG3D_STIMS_H

#include "../main/ode_solver.h"
#include "../alg/grid/grid.h"

void set_benchmark_spatial_stim (struct grid *the_grid, Real stim_current, Real *spatial_currents, struct string_hash *config);
void stim_if_x_less_than (struct grid *the_grid, Real stim_current, Real *spatial_currents, struct string_hash *config);

#endif //MONOALG3D_STIMS_H
