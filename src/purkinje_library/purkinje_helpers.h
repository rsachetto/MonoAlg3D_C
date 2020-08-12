//
// Created by bergolho on 19/07/18.
//

#ifndef MONOALG3D_PURKINJE_HELPERS_H
#define MONOALG3D_PURKINJE_HELPERS_H

#include "../alg/grid/grid.h"

#define UM_TO_CM 0.0001

struct point
{
    real_cpu x, y, z;
};

struct branch
{
    uint32_t source;
    uint32_t destination;
};

void set_custom_purkinje_network (struct grid_purkinje *the_purkinje, const char *file_name, const real_cpu side_length);
void set_purkinje_network_from_file (struct graph *the_purkinje_network, const char *file_name, const real_cpu side_length);
void set_purkinje_coupling_parameters(struct graph *the_purkinje_network, const real_cpu rpmj, const real_cpu pmj_scale, const real_cpu assym_ratio,\
                                    const uint32_t nmin_pmj, const uint32_t nmax_pmj, const bool retro_propagation, const char pmj_filename[]);

void build_skeleton_purkinje (const char *filename, struct graph *skeleton_network);
void build_mesh_purkinje (struct graph *the_purkinje_network, struct graph *skeleton_network, const real_cpu side_length);

void depth_first_search (struct graph *the_purkinje_network, struct node *u, int level, uint32_t *map_skeleton_to_mesh);
void grow_segment (struct graph *the_purkinje_network, struct node *u, struct edge *v, uint32_t *map_skeleton_to_mesh);

void calc_unitary_vector (real_cpu d_ori[], struct node *u, struct node *v);
void calculate_number_of_terminals (struct graph *the_purkinje_network);

void read_purkinje_network_from_file (const char *filename, struct point **points, struct branch **branches, int *N, int *E);
void write_purkinje_network_to_vtk (struct graph *the_purkinje_network);

bool is_terminal (const struct node *n);

int check_purkinje_mesh_for_errors (struct graph *the_purkinje_network);

// TO DO: Other types of network will be implemented here ...

#endif // MONOALG3D_DOMAIN_HELPERS_H
