//
// Created by sachetto on 19/10/17.
//

#ifndef MONOALG3D_DOMAIN_HELPERS_H
#define MONOALG3D_DOMAIN_HELPERS_H

#include "../alg/grid/grid.h"

void set_benchmark_domain (struct grid *the_grid);
void set_cuboid_domain(struct grid *the_grid, real_cpu sizeX, real_cpu sizeY, real_cpu sizeZ);

void set_custom_mesh(struct grid *the_grid, const char *file_name, size_t size,char *read_format);

void set_custom_mesh_with_bounds (struct grid *the_grid, const char *file_name, size_t size,
                                  real_cpu minx, real_cpu maxx, real_cpu miny, real_cpu maxy, real_cpu minz,
                                  real_cpu maxz,  bool read_fibrosis);

void set_cell_not_changeable (struct cell_node *c, real_cpu initial_discretization);

void set_plain_fibrosis (struct grid *the_grid, real_cpu phi, unsigned fib_seed);

void set_plain_fibrosis_and_write_positions_to_file (struct grid *the_grid, real_cpu phi, unsigned fib_seed);
void set_plain_fibrosis_using_file (struct grid *the_grid, const char filename[]);

void set_plain_fibrosis_inside_region (struct grid *the_grid, real_cpu phi, unsigned fib_seed,\
                        const double min_x, const double max_x,\
                        const double min_y, const double max_y,\
                        const double min_z, const double max_z);

void set_plain_source_sink_fibrosis (struct grid *the_grid, real_cpu channel_width, real_cpu channel_length);

void set_plain_sphere_fibrosis (struct grid *the_grid, real_cpu phi, real_cpu plain_center, real_cpu sphere_radius,
                                real_cpu bz_size, real_cpu bz_radius, unsigned fib_seed);

void set_plain_sphere_fibrosis_without_inactivating(struct grid *the_grid, real_cpu plain_center, real_cpu sphere_radius, real_cpu bz_radius);


void set_human_mesh_fibrosis(struct grid *grid, real_cpu phi, unsigned seed, real_cpu big_scar_center_x,
                     real_cpu big_scar_center_y, real_cpu big_scar_center_z, real_cpu small_scar_center_x,
                     real_cpu small_scar_center_y, real_cpu small_scar_center_z);

void set_human_mesh_fibrosis_from_file(struct grid *grid, char type, const char *filename, int size);

int calculate_cuboid_side_lengths(real_cpu start_dx, real_cpu start_dy, real_cpu start_dz, real_cpu side_length_x,
                                   real_cpu side_length_y, real_cpu side_length_z, real_cpu *real_side_length_x,
                                   real_cpu *real_side_length_y, real_cpu *real_side_length_z);

void set_fibrosis_from_file(struct grid *grid, const char *filename, int size);

void set_plain_fibrosis_and_write_positions_to_file (struct grid *the_grid, real_cpu phi, unsigned fib_seed);
void set_plain_fibrosis_using_file (struct grid *the_grid, const char filename[]);
// ---------------------------------------------------------------------------------------------------------------

void set_plain_fibrosis_inside_region (struct grid *the_grid, real_cpu phi, unsigned fib_seed,\
                        const double min_x, const double max_x,\
                        const double min_y, const double max_y,\
                        const double min_z, const double max_z);


#endif // MONOALG3D_DOMAIN_HELPERS_H
