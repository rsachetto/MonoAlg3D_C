//
// Created by sachetto on 29/09/17.
//

#ifndef MONOALG3D_GRID_H
#define MONOALG3D_GRID_H

#include "../../common_types/common_types.h"
#include "../cell/cell.h"
#include "../grid_purkinje/grid_purkinje.h"

#include <stdio.h>
#include <stdlib.h>

#define FOR_EACH_CELL(grid) for(struct cell_node *cell = grid->first_cell; cell != NULL; cell = cell->next)

#define FOR_EACH_PURKINJE_CELL(grid) for(struct cell_node *cell = grid->purkinje->first_cell; cell != NULL; cell = cell->next)

struct grid {
    struct cell_node *first_cell; // First cell of grid.
    struct point_3d cube_side_length;
    struct point_3d mesh_side_length;

    uint32_t number_of_cells; // Number of cells of grid.
    uint32_t num_active_cells;

    // dynamic arrays, need to point to NULL
    ui32_array free_sv_positions;
    ui32_array refined_this_step;

    struct cell_node **active_cells;
    bool adaptive;

    // Purkinje section
    struct grid_purkinje *purkinje;

    struct point_3d start_discretization;
    struct point_3d max_discretization;

    void *extra_info;
};

#ifdef __cplusplus
extern "C" {
#endif

struct grid *new_grid();
void initialize_grid(struct grid *the_grid, struct point_3d side_length);
void clean_and_free_grid(struct grid *the_grid);
void construct_grid(struct grid *the_grid);
void initialize_and_construct_grid(struct grid *the_grid, struct point_3d side_length);

void print_grid(struct grid *the_grid, FILE *output_file);

void clean_grid(struct grid *the_grid);
void order_grid_cells(struct grid *the_grid);

void set_grid_flux(struct grid *the_grid);

void grid_to_csr(struct grid *the_grid, float **A, int **IA, int **JA, bool is_purkinje);
void grid_to_csr_for_ecg(struct grid *the_grid, float **A, int **IA, int **JA, bool is_purkinje, bool for_ecg);

bool refine_grid_with_bound(struct grid *the_grid, real_cpu refinement_bound, real_cpu min_dx, real_cpu min_dy, real_cpu min_dz);
void refine_grid(struct grid *the_grid, int num_steps);
void refine_grid_with_bounds(struct grid *the_grid, int num_steps, struct point_3d min_bounds, struct point_3d max_bounds);
void refine_grid_cell(struct grid *the_grid, struct cell_node *grid_cell);

bool derefine_grid_with_bound(struct grid *the_grid, real_cpu derefinement_bound, real_cpu max_dx, real_cpu max_dy, real_cpu max_dz);
void derefine_all_grid(struct grid *the_grid);
void derefine_grid_inactive_cells(struct grid *the_grid);

void print_grid_matrix(struct grid *the_grid, FILE *output_file);
void print_grid_vector(struct grid *the_grid, FILE *output_file, char name);
real_cpu *grid_vector_to_array(struct grid *the_grid, char name, uint32_t *num_lines);
void print_grid_matrix_as_octave_matrix(struct grid *the_grid, FILE *output_file);

int get_num_refinement_steps_to_discretization(real_cpu side_len, real_cpu h);
void save_grid_domain(struct grid *the_grid, const char *file_name);

void initialize_and_construct_grid_purkinje(struct grid *the_grid);
void initialize_grid_purkinje(struct grid *the_grid);
void construct_grid_purkinje(struct grid *the_grid);

void construct_grid_from_file(struct grid *grid, FILE *matrix_a, FILE *vector_b);

struct terminal *link_purkinje_to_tissue(struct grid *the_grid);
struct terminal *link_purkinje_to_tissue_default(struct grid *the_grid);
struct terminal *link_purkinje_to_tissue_using_pmj_locations(struct grid *the_grid);
void update_link_purkinje_to_endocardium(struct grid *the_grid, struct terminal *the_terminals);
void set_active_terminals(struct terminal *the_terminals, const uint32_t number_of_terminals, const char filename[]);
void print_terminals(struct terminal *the_terminals, const uint32_t number_of_terminals);
void free_terminals(struct terminal *the_terminals, const uint32_t number_of_terminals);

#ifdef __cplusplus
}
#endif

#endif // MONOALG3D_GRID_H
