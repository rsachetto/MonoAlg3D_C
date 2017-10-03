//
// Created by sachetto on 29/09/17.
//

#ifndef MONOALG3D_SOLVER_H
#define MONOALG3D_SOLVER_H

#include <unitypes.h>
#include <stdbool.h>
#include "../utils/output_info.h"
#include "edo_solver.h"
#include "../grid/grid.h"

//TODO: should we separate edp solver from monodomain?
struct monodomain_solver {
    int numer_of_iterations;
    int refine_each;
    int derefine_each;
    float tolerance;
    int num_threads;
    int max_iterations;
    double final_time;

    //TODO: create a file for stimulus definition!!
    double stim_start;
    double stim_dur;
    uint64_t *cells_to_solve;

    //TODO: this should be in the monodomain solver maybe???
    double beta, cm; //micrometers
    double initial_V;
    double sigma_x;
    double sigma_y;
    double sigma_z;



    // Time used for solving wave equation.
    double dt;

};

void solve_monodomain( struct grid *the_grid,  struct monodomain_solver *the_monodomain_solver, struct edo_solver *the_edo_solver, struct output_info *output_info );
void init_solver(struct monodomain_solver* the_solver);

void save_old_cell_positions(struct grid* the_grid);
void update_cells_to_solve(struct grid* the_grid, struct monodomain_solver *solver);
void set_initial_conditions(struct monodomain_solver *the_solver, struct grid *the_grid);

void initialize_diagonal_elements(struct monodomain_solver *the_solver, struct grid *the_grid);
void fill_discretization_matrix_elements(struct monodomain_solver *the_solver, struct cell_node *grid_cell, void* neighbor_grid_cell, char direction);
#endif //MONOALG3D_SOLVER_H
