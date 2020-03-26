//
// Created by sachetto on 27/09/18.
//

#ifndef MONOALG3D_COMMON_DATA_STRUCTURES_H
#define MONOALG3D_COMMON_DATA_STRUCTURES_H

#include <stdbool.h>
#include "../monodomain/constants.h"

#define Pragma(x) _Pragma(#x)
#define OMP(directive) Pragma(omp directive)

// This is used when we are dealing with fibrotic meshes
struct fibrotic_mesh_info {
    bool fibrotic;
    bool border_zone;
    char scar_type;
};

#define FIBROTIC_INFO(grid_cell) (struct fibrotic_mesh_info *)grid_cell->mesh_extra_info
#define FIBROTIC(grid_cell) (FIBROTIC_INFO(grid_cell))->fibrotic
#define BORDER_ZONE(grid_cell) (FIBROTIC_INFO(grid_cell))->border_zone
#define SCAR_TYPE(grid_cell) (FIBROTIC_INFO(grid_cell))->scar_type

#define INITIALIZE_FIBROTIC_INFO(grid_cell)                                                                            \
    do {                                                                                                               \
        size_t __size__ = sizeof (struct fibrotic_mesh_info);                                                          \
        (grid_cell)->mesh_extra_info = malloc (__size__);                                                              \
        (grid_cell)->mesh_extra_info_size = __size__;                                                                  \
        FIBROTIC ((grid_cell)) = false;                                                                                \
        BORDER_ZONE (grid_cell) = false;                                                                               \
        SCAR_TYPE ((grid_cell)) = 'n';                                                                                 \
} while (0)

struct conjugate_gradient_info {
    real_cpu r;  /* Element of the int_vector r = b - Ax associated to this cell. */
    real_cpu p;  /* Element of the search direction int_vector in the conjugate gradient algorithm. */
    real_cpu p1; /* p's upgrade in the conjugate gradient algorithm. */
    real_cpu z;  // Jacobi preconditioner
};

struct biconjugate_gradient_info {
    real_cpu r;  /* Element of the int_vector r = b - Ax associated to this cell. */
    real_cpu p;  /* Element of the search direction int_vector in the conjugate gradient algorithm. */
    real_cpu p1; /* p's upgrade in the conjugate gradient algorithm. */
    real_cpu z;  // Jacobi preconditioner

    real_cpu x_aux;
    real_cpu r_aux;
    real_cpu z_aux;
    real_cpu p_aux;
    real_cpu p1_aux;
    real_cpu xA;
};

struct jacobi_info {
    real_cpu x_aux; // jacobi variable
};

#define CG_INFO(grid_cell) (struct conjugate_gradient_info *)grid_cell->linear_system_solver_extra_info
#define CG_R(grid_cell) (CG_INFO(grid_cell))->r
#define CG_P(grid_cell) (CG_INFO(grid_cell))->p
#define CG_P1(grid_cell) (CG_INFO(grid_cell))->p1
#define CG_Z(grid_cell) (CG_INFO(grid_cell))->z

#define BCG_INFO(grid_cell) (struct biconjugate_gradient_info *)grid_cell->linear_system_solver_extra_info
#define BCG_R(grid_cell) (BCG_INFO(grid_cell))->r
#define BCG_P(grid_cell) (BCG_INFO(grid_cell))->p
#define BCG_P1(grid_cell) (BCG_INFO(grid_cell))->p1
#define BCG_Z(grid_cell) (BCG_INFO(grid_cell))->z
#define BCG_X_AUX(grid_cell) (BCG_INFO(grid_cell))->x_aux
#define BCG_R_AUX(grid_cell) (BCG_INFO(grid_cell))->r_aux
#define BCG_Z_AUX(grid_cell) (BCG_INFO(grid_cell))->z_aux
#define BCG_P_AUX(grid_cell) (BCG_INFO(grid_cell))->p_aux
#define BCG_P1_AUX(grid_cell) (BCG_INFO(grid_cell))->p1_aux
#define BCG_XA(grid_cell) (BCG_INFO(grid_cell))->xA

#define JACOBI_INFO(grid_cell) (struct jacobi_info *)grid_cell->linear_system_solver_extra_info
#define JACOBI_X_AUX(grid_cell) (JACOBI_INFO(grid_cell))->x_aux

#define INITIALIZE_LINEAR_SYSTEM_SOLVER_INFO(grid_cell, info_struct)                                                   \
    do {                                                                                                               \
        size_t __size__ = sizeof (struct info_struct);                                                                 \
        (grid_cell)->linear_system_solver_extra_info = malloc (__size__);                                              \
        (grid_cell)->linear_system_solver_extra_info_size = __size__;                                                  \
    } while (0)

#define INITIALIZE_CONJUGATE_GRADIENT_INFO(grid_cell) INITIALIZE_LINEAR_SYSTEM_SOLVER_INFO (grid_cell, conjugate_gradient_info)
#define INITIALIZE_BICONJUGATE_GRADIENT_INFO(grid_cell) INITIALIZE_LINEAR_SYSTEM_SOLVER_INFO (grid_cell, biconjugate_gradient_info)
#define INITIALIZE_JACOBI_INFO(grid_cell) INITIALIZE_LINEAR_SYSTEM_SOLVER_INFO (grid_cell, jacobi_info)

#endif // MONOALG3D_COMMON_DATA_STRUCTURES_H
