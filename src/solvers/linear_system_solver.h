//
// Created by sachetto on 04/10/17.
//

#ifndef MONOALG3D_LINEAR_SYSTEM_SOLVER_C_H
#define MONOALG3D_LINEAR_SYSTEM_SOLVER_C_H

#include "../grid/grid.h"
#include <stdbool.h>

uint64_t conjugate_gradient(struct grid* the_grid, int max_its, double tol, bool use_jacobi, double *error);

#endif //MONOALG3D_LINEAR_SYSTEM_SOLVER_C_H
