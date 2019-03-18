//
// Created by sachetto on 29/09/17.
//

#ifndef MONOALG3D_GRID_H
#define MONOALG3D_GRID_H

#include "../cell/cell.h"
#include "../../graph/graph.h"
#include "../../common_types/common_types.h"

#include <stdlib.h>
#include <stdio.h>

struct grid {

    struct cell_node *first_cell;     // First cell of grid.

    real_cpu side_length_x;
    real_cpu side_length_y;
    real_cpu side_length_z;

    uint32_t number_of_cells;  // Number of cells of grid.

    uint32_t num_active_cells;

    //dynamic arrays, need to point to NULL
    ui32_array free_sv_positions;
    ui32_array refined_this_step;

    struct cell_node **active_cells;
    bool adaptive;

    struct graph *the_purkinje_network;

};


struct grid* new_grid();
void initialize_grid(struct grid *the_grid, real_cpu side_length_x, real_cpu side_length_y, real_cpu side_length_z);
void clean_and_free_grid(struct grid* the_grid);
void construct_grid(struct grid *the_grid);
void initialize_and_construct_grid(struct grid *the_grid, real_cpu side_length_x, real_cpu side_length_y, real_cpu side_length_z);

void print_grid(struct grid* the_grid, FILE *output_file);
void print_grid_with_scar_info(struct grid *the_grid, FILE *output_file, bool binary);

void clean_grid(struct grid *the_grid);
void order_grid_cells (struct grid *the_grid);

void set_grid_flux(struct grid *the_grid);

bool refine_grid_with_bound(struct grid* the_grid, real_cpu refinement_bound,  real_cpu min_dx, real_cpu min_dy, real_cpu min_dz);
void refine_grid(struct grid* the_grid, int num_steps);
void refine_grid_cell(struct grid *the_grid, struct cell_node* grid_cell);
void refine_fibrotic_cells(struct grid *the_grid);
void refine_border_zone_cells(struct grid *the_grid);

bool derefine_grid_with_bound (struct grid *the_grid, real_cpu derefinement_bound, real_cpu max_dx, real_cpu max_dy, real_cpu max_dz);
void derefine_all_grid (struct grid* the_grid);
void derefine_grid_inactive_cells (struct grid* the_grid);
bool derefine_grid_cell (struct grid *the_grid, struct cell_node *grid_cell);


void print_grid_matrix(struct grid *the_grid, FILE* output_file);
void print_grid_vector(struct grid* the_grid, FILE *output_file, char name);
real_cpu * grid_vector_to_array(struct grid *the_grid, char name, uint32_t *num_lines);
void print_grid_matrix_as_octave_matrix(struct grid *the_grid, FILE *output_file);

int get_num_refinement_steps_to_discretization (real_cpu side_len, real_cpu h);
void save_grid_domain (struct grid * the_grid, const char *file_name);

void initialize_and_construct_grid_purkinje (struct grid *the_grid);
void initialize_grid_purkinje (struct grid *the_grid);
void construct_grid_purkinje (struct grid *the_grid);

#endif //MONOALG3D_GRID_H
