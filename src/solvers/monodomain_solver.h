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

struct monodomain_solver {

    //TODO: @Check: maybe this solver need a alg grid?

    Real tolerance;
    int num_threads;
    int max_iterations;
    Real final_time;

    Real beta, cm; // micrometers
    Real sigma_x;
    Real sigma_y;
    Real sigma_z;

    Real start_h, max_h, min_h;
    bool adaptive;
    int refine_each;
    int derefine_each;
    Real refinement_bound;
    Real derefinement_bound;

    bool abort_on_no_activity;

    // Time used for solving wave equation.
    Real dt;
    bool use_jacobi;


};

void solve_monodomain (struct grid *the_grid, struct monodomain_solver *the_monodomain_solver,
                       struct ode_solver *the_edo_solver, struct output_utils *output_info);
void init_solver (struct monodomain_solver *the_solver);

void save_old_cell_positions (struct grid *the_grid);
void update_cells_to_solve (struct grid *the_grid, struct ode_solver *solver);
void set_initial_conditions (struct monodomain_solver *the_solver, struct grid *the_grid, Real initial_v);

void initialize_diagonal_elements (struct monodomain_solver *the_solver, struct grid *the_grid);

void fill_discretization_matrix_elements (struct monodomain_solver *the_solver,
                                          struct cell_node *grid_cell, void *neighbor_grid_cell,
                                          char direction);

void set_discretization_matrix (struct monodomain_solver *the_solver, struct grid *the_grid);

void print_solver_info (struct monodomain_solver *the_monodomain_solver,
                        struct ode_solver *the_ode_solver, struct grid *the_grid,
                        struct output_utils *output_info);

void print_grid_matrix(struct grid *the_grid, FILE* output_file);

void update_ode_state_vector(struct ode_solver *the_ode_solver, struct grid *the_grid, bool adaptive);

void set_ode_extra_data(struct grid* the_grid, struct ode_solver *the_ode_solver);
void set_spatial_stim(struct grid* the_grid, struct ode_solver *the_ode_solver);

#endif // MONOALG3D_SOLVER_H
