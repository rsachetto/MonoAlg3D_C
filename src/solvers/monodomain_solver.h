//
// Created by sachetto on 29/09/17.
//

#ifndef MONOALG3D_SOLVER_H
#define MONOALG3D_SOLVER_H

#include "../grid/grid/grid.h"
#include "../utils/output_utils.h"
#include "ode_solver.h"
#include <stdbool.h>
#include <unitypes.h>

struct monodomain_solver {

    //TODO: @Check: maybe this solver need a alg grid?

    double tolerance;
    int num_threads;
    int max_iterations;
    double final_time;

    double beta, cm; // micrometers
    double sigma_x;
    double sigma_y;
    double sigma_z;

    double start_h, max_h, min_h;
    int refine_each;
    int derefine_each;
    double refinement_bound;
    double derefinement_bound;

    bool abort_on_no_activity;

    // Time used for solving wave equation.
    double dt;
    bool use_jacobi;


};

void solve_monodomain (struct grid *the_grid, struct monodomain_solver *the_monodomain_solver,
                       struct ode_solver *the_edo_solver, struct output_utils *output_info);
void init_solver (struct monodomain_solver *the_solver);

void save_old_cell_positions (struct grid *the_grid);
void update_cells_to_solve (struct grid *the_grid, struct ode_solver *solver);
void set_initial_conditions (struct monodomain_solver *the_solver, struct grid *the_grid, double initial_v);

void initialize_diagonal_elements (struct monodomain_solver *the_solver, struct grid *the_grid);

void fill_discretization_matrix_elements(struct monodomain_solver *the_solver, struct cell_node *grid_cell,
                                         void *neighbor_grid_cell, char direction, uint8_t max_elements);

void set_discretization_matrix (struct monodomain_solver *the_solver, struct grid *the_grid);

void print_solver_info (struct monodomain_solver *the_monodomain_solver,
                        struct ode_solver *the_ode_solver, struct grid *the_grid,
                        struct output_utils *output_info);

void update_ode_state_vector(struct ode_solver *the_ode_solver, struct grid *the_grid, uint32_t max_number_of_cells);

void set_ode_extra_data(struct grid* the_grid, struct ode_solver *the_ode_solver);
void set_spatial_stim(struct grid* the_grid, struct ode_solver *the_ode_solver);

void update_monodomain(uint32_t initial_number_of_cells, uint32_t num_active_cells, struct cell_node **active_cells,
                       double beta,
                       double cm, double dt_edp, Real *sv, int n_equations_cell_model, bool use_gpu);


#endif // MONOALG3D_SOLVER_H
