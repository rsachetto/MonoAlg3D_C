#ifndef MONOALG3D_GRID_DOMAIN_H
#define MONOALG3D_GRID_DOMAIN_H

#include "../main/config/domain_config.h"
#include "../alg/grid/grid.h"


void set_plain_sphere_fibrosis(struct grid* the_grid, double phi,  double plain_center, double sphere_radius, double bz_size,
                               double bz_radius);
void set_plain_fibrosis(struct grid* the_grid, double phi);
int get_num_refinement_steps_to_discretization (double side_len, double h);

void initialize_grid_with_mouse_mesh (struct grid *the_grid, struct domain_config *domain_config);
void initialize_grid_with_rabbit_mesh (struct grid *the_grid, struct domain_config *domain_config);
void initialize_grid_with_benchmark_mesh (struct grid *the_grid, struct domain_config *domain_config);

void initialize_grid_with_plain_mesh (struct grid *the_grid, struct domain_config *domain_config);
void initialize_grid_with_plain_fibrotic_mesh(struct grid *the_grid, double side_length, double start_h, int num_layers,
                                              double phi);
void initialize_grid_with_plain_and_sphere_fibrotic_mesh(struct grid *the_grid, double side_length,
                                                         double start_h, int num_layers, double phi,
                                                         double plain_center, double sphere_radius, double bz_size,
                                                         double bz_radius);

#endif //MONOALG3D_GRID_DOMAIN_H