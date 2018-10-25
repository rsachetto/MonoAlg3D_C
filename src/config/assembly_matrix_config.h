//
// Created by sachetto on 13/10/17.
//

#ifndef MONOALG3D_ASSEMBLY_CONFIG_H
#define MONOALG3D_ASSEMBLY_CONFIG_H

#include "config_common.h"
#include "../monodomain/constants.h"
#include "../monodomain/monodomain_solver.h"

struct assembly_matrix_config;
struct monodomain_solver;

#define ASSEMBLY_MATRIX(name) EXPORT_FN void name(struct assembly_matrix_config *config, struct monodomain_solver *the_solver, struct grid *the_grid)
typedef ASSEMBLY_MATRIX(assembly_matrix_fn);

struct assembly_matrix_config {
    struct config_common config_data;
    assembly_matrix_fn *assembly_matrix;
};

void init_assembly_matrix_functions(struct assembly_matrix_config *config);
struct assembly_matrix_config* new_assembly_matrix_config();
void free_assembly_matrix_config(struct assembly_matrix_config* s);
void print_assembly_matrix_config_values(struct assembly_matrix_config* s);


#endif //MONOALG3D_ASSEMBLY_CONFIG_H
