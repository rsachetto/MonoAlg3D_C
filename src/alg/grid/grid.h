//
// Created by sachetto on 29/09/17.
//

#ifndef MONOALG3D_GRID_H
#define MONOALG3D_GRID_H

#include "../cell/cell.h"

#include <stdlib.h>
#include <stdio.h>

struct grid {

    struct cell_node *first_cell;     // First cell of grid.
    double side_length;        // Length of cube grid. Default = 1.0.
    uint32_t number_of_cells;  // Number of cells of grid.

    uint32_t num_active_cells;

    //uint32_vector *free_sv_positions;
    //uint32_vector *refined_this_step;

    //dynamic arrays, need to point to NULL
    uint32_t *free_sv_positions; 
    uint32_t *refined_this_step;

    struct cell_node* *active_cells;
    bool adaptive;

};


struct grid* new_grid();
void initialize_grid(struct grid *the_grid, double side_length);
void clean_and_free_grid(struct grid* the_grid);
void construct_grid(struct grid *the_grid);
void initialize_and_construct_grid(struct grid *the_grid, double side_length);

void print_grid(struct grid* the_grid, FILE *output_file);
void print_grid_with_scar_info(struct grid *the_grid, FILE *output_file, bool binary);

bool print_grid_and_check_for_activity(const struct grid *the_grid, FILE *output_file, const int count, const bool binary);

void clean_grid(struct grid *the_grid);
void order_grid_cells (struct grid *the_grid);

void set_grid_flux(struct grid *the_grid);

bool refine_grid_with_bound(struct grid* the_grid, double refinement_bound, double min_h);
void refine_grid(struct grid* the_grid, int num_steps);
void refine_grid_cell_at(struct grid *the_grid, uint32_t cell_number);
void refine_fibrotic_cells(struct grid *the_grid);
void refine_border_zone_cells(struct grid *the_grid);

bool derefine_grid_with_bound(struct grid *the_grid, double derefinement_bound, double max_h);
void derefine_all_grid (struct grid* the_grid);
void derefine_grid_inactive_cells (struct grid* the_grid);

void print_grid_matrix(struct grid *the_grid, FILE* output_file);
void print_grid_vector(struct grid* the_grid, FILE *output_file, char name);
double * grid_vector_to_array(struct grid *the_grid, char name, uint32_t *num_lines);

void save_grid_domain (struct grid * the_grid, const char *file_name);

void lock_grid(struct grid *the_grid);

void unlock_grid(struct grid *the_grid);



#endif //MONOALG3D_GRID_H
