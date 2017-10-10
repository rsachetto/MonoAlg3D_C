//
// Created by sachetto on 29/09/17.
//

#ifndef MONOALG3D_GRID_H
#define MONOALG3D_GRID_H

#include "../cell/cell.h"
#include "../../utils/vector/int_vector.h"
#include "../../utils/vector/uint32_vector.h"

#include <stdlib.h>
#include <stdio.h>

struct grid {

    struct cell_node *first_cell;     // First cell of grid.
    double side_length;        // Length of cube grid. Default = 1.0.
    uint64_t number_of_cells;  // Number of cells of grid.

    uint64_t num_active_cells;
    uint8_t num_cell_neighbours;

    uint32_vector *free_sv_positions;
    uint32_vector *refined_this_step;

    struct cell_node* *active_cells;

    bool init_ode;
    bool adaptive;

};



void initialize_grid(struct grid *the_grid, double side_length, uint8_t num_cell_neighbours);
void clean_and_free_grid(struct grid* the_grid);
void construct_grid(struct grid *the_grid);
void initialize_and_construct_grid(struct grid *the_grid, double side_length, uint8_t num_cell_neighbours);

void print_grid(struct grid* the_grid, FILE *output_file);
bool print_grid_and_check_for_activity (struct grid *the_grid, FILE *output_file, int count);

void clean_grid(struct grid *the_grid);
void order_grid_cells (struct grid *the_grid);

void set_grid_flux(struct grid *the_grid);

bool refine_grid_with_bound(struct grid* the_grid, double refinement_bound, double min_h);
void refine_grid(struct grid* the_grid, int num_steps);
void refine_grid_cell_at(struct grid* the_grid, uint64_t cell_number );

bool derefine_grid_with_bound(struct grid *the_grid, double derefinement_bound, double max_h);
void derefine_all_grid (struct grid* the_grid);
void derefine_grid_inactive_cells (struct grid* the_grid);

void print_grid_matrix(struct grid *the_grid, FILE* output_file);
void print_grid_vector(struct grid* the_grid, FILE *output_file, char name);
double * grid_vector_to_array(struct grid* the_grid, char name, uint64_t *num_lines);

void initialize_grid_with_mouse_mesh (struct grid *the_grid, const char *mesh_file);
void initialize_grid_with_rabbit_mesh (struct grid *the_grid, const char *mesh_file);
void initialize_grid_with_benchmark_mesh (struct grid *the_grid, double start_h);

void initialize_grid_with_plain_mesh (struct grid *the_grid, double desired_side_lenght, double start_h, int num_layers);
void initialize_grid_with_plain_fibrotic_mesh (struct grid *the_grid, double side_length, double start_h, double num_layers, double phi);
void initialize_grid_with_plain_and_sphere_fibrotic_mesh (struct grid *the_grid, double side_length,
                                                          double start_h, double num_layers, double phi,
                                                          double plain_center, double sphere_radius, double bz_size,
                                                          double bz_radius);
void set_plain_sphere_fibrosis(struct grid* the_grid, double phi,  double plain_center, double sphere_radius, double bz_size,
                               double bz_radius);
void set_plain_fibrosis(struct grid* the_grid, double phi);

#endif //MONOALG3D_GRID_H
