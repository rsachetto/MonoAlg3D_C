#ifndef MONOALG3D_GRID_DOMAIN_H
#define MONOALG3D_GRID_DOMAIN_H

#include "../main/config/domain_config.h"
#include "../alg/grid/grid.h"

void initialize_grid_with_mouse_mesh (struct grid *the_grid, struct domain_config *domain_config);
void initialize_grid_with_rabbit_mesh (struct grid *the_grid, struct domain_config *domain_config);
void initialize_grid_with_benchmark_mesh (struct grid *the_grid, struct domain_config *domain_config);

void initialize_grid_with_plain_mesh (struct grid *the_grid, struct domain_config *domain_config);
void initialize_grid_with_plain_fibrotic_mesh(struct grid *the_grid, struct domain_config* config);
void initialize_grid_with_plain_and_sphere_fibrotic_mesh(struct grid *the_grid, struct domain_config* config);

#endif //MONOALG3D_GRID_DOMAIN_H