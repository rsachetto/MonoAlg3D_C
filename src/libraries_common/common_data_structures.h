//
// Created by sachetto on 27/09/18.
//

#ifndef MONOALG3D_COMMON_DATA_STRUCTURES_H
#define MONOALG3D_COMMON_DATA_STRUCTURES_H

#include <stdbool.h>

// This is used when we are dealing with fibrotic meshes
struct fibrotic_mesh_info {
    bool fibrotic;
    bool border_zone;
    char scar_type;
};

#define FIBROTIC(grid_cell) ((struct fibrotic_mesh_info *)(grid_cell)->mesh_extra_info)->fibrotic
#define BORDER_ZONE(grid_cell) ((struct fibrotic_mesh_info *)(grid_cell)->mesh_extra_info)->border_zone
#define SCAR_TYPE(grid_cell) ((struct fibrotic_mesh_info *)(grid_cell)->mesh_extra_info)->scar_type
#define FIBROTIC_INFO(grid_cell) (struct fibrotic_mesh_info *)grid_cell->mesh_extra_info

#define INITIALIZE_FIBROTIC_INFO(grid_cell)                                                                            \
    do {                                                                                                               \
        size_t size = sizeof (struct fibrotic_mesh_info);                                                              \
        (grid_cell)->mesh_extra_info = malloc (size);                                                                  \
        (grid_cell)->mesh_extra_info_size = size;                                                                      \
        FIBROTIC ((grid_cell)) = false;                                                                                \
        BORDER_ZONE (grid_cell) = false;                                                                               \
        SCAR_TYPE ((grid_cell)) = 'n';                                                                                 \
    } while (0)

struct conjugate_gradient_info {
    double r;  /* Element of the int_vector r = b - Ax associated to this cell. */
    double p;  /* Element of the search direction int_vector in the conjugate gradient algorithm. */
    double p1; /* p's upgrade in the conjugate gradient algorithm. */
    double z;  // Jacobi preconditioner
};

struct biconjugate_gradient_info {
    double r;  /* Element of the int_vector r = b - Ax associated to this cell. */
    double p;  /* Element of the search direction int_vector in the conjugate gradient algorithm. */
    double p1; /* p's upgrade in the conjugate gradient algorithm. */
    double z;  // Jacobi preconditioner

    double x_aux;
    double r_aux;
    double z_aux;
    double p_aux;
    double p1_aux;
    double xA;
};

struct jacobi_info {
    double x_aux; // jacobi variable
};

#define CG_INFO(grid_cell) (struct conjugate_gradient_info *)grid_cell->linear_system_solver_extra_info
#define CG_R(grid_cell) ((struct conjugate_gradient_info *)(grid_cell)->linear_system_solver_extra_info)->r
#define CG_P(grid_cell) ((struct conjugate_gradient_info *)(grid_cell)->linear_system_solver_extra_info)->p
#define CG_P1(grid_cell) ((struct conjugate_gradient_info *)(grid_cell)->linear_system_solver_extra_info)->p1
#define CG_Z(grid_cell) ((struct conjugate_gradient_info *)(grid_cell)->linear_system_solver_extra_info)->z

#define BCG_INFO(grid_cell) (struct biconjugate_gradient_info *)grid_cell->linear_system_solver_extra_info
#define BCG_R(grid_cell) ((struct biconjugate_gradient_info *)(grid_cell)->linear_system_solver_extra_info)->r
#define BCG_P(grid_cell) ((struct biconjugate_gradient_info *)(grid_cell)->linear_system_solver_extra_info)->p
#define BCG_P1(grid_cell) ((struct biconjugate_gradient_info *)(grid_cell)->linear_system_solver_extra_info)->p1
#define BCG_Z(grid_cell) ((struct biconjugate_gradient_info *)(grid_cell)->linear_system_solver_extra_info)->z
#define BCG_X_AUX(grid_cell) ((struct biconjugate_gradient_info *)(grid_cell)->linear_system_solver_extra_info)->x_aux
#define BCG_R_AUX(grid_cell) ((struct biconjugate_gradient_info *)(grid_cell)->linear_system_solver_extra_info)->r_aux
#define BCG_Z_AUX(grid_cell) ((struct biconjugate_gradient_info *)(grid_cell)->linear_system_solver_extra_info)->z_aux
#define BCG_P_AUX(grid_cell) ((struct biconjugate_gradient_info *)(grid_cell)->linear_system_solver_extra_info)->p_aux
#define BCG_P1_AUX(grid_cell) ((struct biconjugate_gradient_info *)(grid_cell)->linear_system_solver_extra_info)->p1_aux
#define BCG_XA(grid_cell) ((struct biconjugate_gradient_info *)(grid_cell)->linear_system_solver_extra_info)->xA

#define JACOBI_INFO(grid_cell) (struct jacobi_info *)grid_cell->linear_system_solver_extra_info
#define JACOBI_X_AUX(grid_cell) ((struct jacobi_info *)(grid_cell)->linear_system_solver_extra_info)->x_aux

#define INITIALIZE_LINEAR_SYSTEM_SOLVER_INFO(grid_cell, info_struct)                                                                 \
    do {                                                                                                               \
        size_t size = sizeof (struct info_struct);                                                                     \
        (grid_cell)->linear_system_solver_extra_info = malloc (size);                                                  \
        (grid_cell)->linear_system_solver_extra_info_size = size;                                                      \
    } while (0)

#define INITIALIZE_CONJUGATE_GRADIENT_INFO(grid_cell) INITIALIZE_LINEAR_SYSTEM_SOLVER_INFO (grid_cell, conjugate_gradient_info)
#define INITIALIZE_BICONJUGATE_GRADIENT_INFO(grid_cell) INITIALIZE_LINEAR_SYSTEM_SOLVER_INFO (grid_cell, biconjugate_gradient_info)
#define INITIALIZE_JACOBI_INFO(grid_cell) INITIALIZE_LINEAR_SYSTEM_SOLVER_INFO (grid_cell, jacobi_info)

#endif // MONOALG3D_COMMON_DATA_STRUCTURES_H
