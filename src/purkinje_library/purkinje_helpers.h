//
// Created by bergolho on 19/07/18.
//

#ifndef MONOALG3D_PURKINJE_HELPERS_H
#define MONOALG3D_PURKINJE_HELPERS_H

#include "../alg/grid/grid.h"

#define UM_TO_CM 0.0001

struct point
{
    double x, y, z;
};

struct branch
{
    uint32_t source;
    uint32_t destination;
};

void set_custom_purkinje_network (struct grid *the_grid, const char *file_name, const double side_length);
void set_purkinje_network_from_file (struct graph *the_purkinje_network, const char *file_name, const double side_length);
void build_skeleton_purkinje (const char *filename, struct graph *skeleton_network);
void build_mesh_purkinje (struct graph *the_purkinje_network, struct graph *skeleton_network, const double side_length);
void depth_first_search (struct graph *the_purkinje_network, struct node *u, int level, uint32_t *map_skeleton_to_mesh);
void grow_segment (struct graph *the_purkinje_network, struct node *u, struct edge *v, uint32_t *map_skeleton_to_mesh);
void calc_unitary_vector (double d_ori[], struct node *u, struct node *v);
void read_purkinje_network_from_file (const char *filename, struct point **points, struct branch **branches, int *N, int *E);
void write_purkinje_network_to_vtk (struct graph *the_purkinje_network);

int check_purkinje_input (const double side_length);

// TO DO: Other types of network will be implemented here ...

#endif // MONOALG3D_DOMAIN_HELPERS_H
