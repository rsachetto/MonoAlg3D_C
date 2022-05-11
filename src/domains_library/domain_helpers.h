//
// Created by sachetto on 19/10/17.
//

#ifndef MONOALG3D_DOMAIN_HELPERS_H
#define MONOALG3D_DOMAIN_HELPERS_H

#include "../alg/grid/grid.h"
#include "../config/domain_config.h"

struct custom_mesh_basic_data_hash_entry {
    struct point_3d key;
    int32_t value;
};

int set_cuboid_domain_mesh(struct grid *the_grid, real_cpu start_dx, real_cpu start_dy, real_cpu start_dz, real_cpu side_length_x, real_cpu side_length_y,
                           real_cpu side_length_z);
int set_square_mesh(struct config *config, struct grid *the_grid);

void set_benchmark_domain(struct grid *the_grid);
void set_cuboid_domain(struct grid *the_grid, real_cpu sizeX, real_cpu sizeY, real_cpu sizeZ);

void set_custom_mesh(struct grid *the_grid, const char *file_name, size_t size, char *read_format);

void set_custom_mesh_with_bounds(struct grid *the_grid, const char *file_name, size_t size, real_cpu minx, real_cpu maxx, real_cpu miny, real_cpu maxy,
                                 real_cpu minz, real_cpu maxz, bool read_fibrosis);

void set_cell_not_changeable(struct cell_node *c, real_cpu initial_discretization);

void set_plain_fibrosis(struct grid *the_grid, real_cpu phi, unsigned fib_seed);

void set_plain_source_sink_fibrosis(struct grid *the_grid, real_cpu channel_width, real_cpu channel_length);

void set_plain_sphere_fibrosis(struct grid *the_grid, real_cpu phi, real_cpu plain_center, real_cpu sphere_radius, real_cpu bz_size, real_cpu bz_radius,
                               unsigned fib_seed);

void set_plain_sphere_fibrosis_without_inactivating(struct grid *the_grid, real_cpu plain_center, real_cpu sphere_radius, real_cpu bz_radius);

int calculate_cuboid_side_lengths(real_cpu start_dx, real_cpu start_dy, real_cpu start_dz, real_cpu side_length_x, real_cpu side_length_y,
                                  real_cpu side_length_z, real_cpu *real_side_length_x, real_cpu *real_side_length_y, real_cpu *real_side_length_z);

void set_fibrosis_from_file(struct grid *grid, const char *filename, int size);

void set_plain_fibrosis_using_file(struct grid *the_grid, const char filename[]);

void set_plain_fibrosis_inside_region(struct grid *the_grid, real_cpu phi, unsigned fib_seed, double min_x, double max_x, double min_y, double max_y,
                                      double min_z, double max_z);

uint32_t set_custom_mesh_from_file(struct grid *the_grid, const char *mesh_file, uint32_t num_volumes, double start_h, uint8_t num_extra_fields,
                                   set_custom_data_for_mesh_fn set_custom_data_for_mesh);


void set_cube_sphere_fibrosis(struct grid *the_grid, real_cpu phi, real_cpu sphere_center[3], real_cpu sphere_radius, unsigned fib_seed);

int calc_num_refs(real_cpu start_h, real_cpu desired_h);

#endif // MONOALG3D_DOMAIN_HELPERS_H
