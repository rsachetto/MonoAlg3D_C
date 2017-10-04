//
// Created by sachetto on 29/09/17.
//

#ifndef MONOALG3D_SOLVER_H
#define MONOALG3D_SOLVER_H

#include "../grid/grid.h"
#include "../utils/output_utils.h"
#include "ode_solver.h"
#include <stdbool.h>
#include <unitypes.h>

// TODO: should we separate edp solver from monodomain?
struct monodomain_solver {
    int numer_of_iterations;

    float tolerance;
    int num_threads;
    int max_iterations;
    double final_time;

    // TODO: create a file for stimulus definition!!
    double stim_start;
    double stim_duration;
    double stim_current;
    uint64_t *cells_to_solve;

    // TODO: this should be in the monodomain solver maybe???
    double beta, cm; // micrometers
    double initial_V;
    double sigma_x;
    double sigma_y;
    double sigma_z;

    double start_h, max_h, min_h;
    bool adaptive;
    int refine_each;
    int derefine_each;
    float refinement_bound;
    float derefinement_bound;

    bool abort_on_no_activity;

    // Time used for solving wave equation.
    double dt;
    bool use_jacobi;
};

void solve_monodomain (struct grid *the_grid, struct monodomain_solver *the_monodomain_solver,
                       struct ode_solver *the_edo_solver, struct output_utils *output_info);
void init_solver (struct monodomain_solver *the_solver);

void save_old_cell_positions (struct grid *the_grid);
void update_cells_to_solve (struct grid *the_grid, struct monodomain_solver *solver);
void set_initial_conditions (struct monodomain_solver *the_solver, struct grid *the_grid);

void initialize_diagonal_elements (struct monodomain_solver *the_solver, struct grid *the_grid);

void fill_discretization_matrix_elements (struct monodomain_solver *the_solver,
                                          struct cell_node *grid_cell, void *neighbor_grid_cell,
                                          char direction);

void set_discretization_matrix (struct monodomain_solver *the_solver, struct grid *the_grid);

void print_solver_info (struct monodomain_solver *the_monodomain_solver,
                        struct ode_solver *the_ode_solver, struct grid *the_grid,
                        struct output_utils *output_info);

void print_grid_matrix(struct grid *the_grid, FILE* output_file);

#endif // MONOALG3D_SOLVER_H
