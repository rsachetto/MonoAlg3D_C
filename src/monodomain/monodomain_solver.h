#ifndef MONOALG3D_SOLVER_H
#define MONOALG3D_SOLVER_H

#include "../alg/grid/grid.h"
#include "ode_solver.h"
#include "../config/config_parser.h"

#include <stdbool.h>
#include <stdint.h>

#define RESTART_SIMULATION 42
#define END_SIMULATION 43
#define SIMULATION_FINISHED 44

//Forward declarations
struct user_options;
struct ode_solver;

struct monodomain_solver {

    int num_threads;
    real_cpu final_time;

    real_cpu beta, cm; // micrometers

    //TODO: maybe use an extra data variable as we did on the alg cell
    real_cpu kappa_x, kappa_y, kappa_z;
    
    int refine_each;
    int derefine_each;
    real_cpu refinement_bound;
    real_cpu derefinement_bound;

    real_cpu start_adapting_at;

    bool abort_on_no_activity;

    // NEW VARIABLE!
    bool calc_activation_time;
    bool print_conductivity;
    bool print_min_vm;
    bool print_max_vm;
    bool print_apd;

    // Time used for solving wave equation.
    real_cpu dt;
    real_cpu current_time;
    int current_count;

};

struct monodomain_solver *new_monodomain_solver ();

int solve_monodomain(struct monodomain_solver *the_monodomain_solver, struct ode_solver *the_ode_solver,
                      struct grid *the_grid, struct user_options *configs);

void save_old_cell_positions (struct grid *the_grid);
void update_cells_to_solve (struct grid *the_grid, struct ode_solver *solver);
void set_initial_conditions (struct monodomain_solver *the_solver, struct grid *the_grid, real_cpu initial_v);

void print_solver_info(struct monodomain_solver *the_monodomain_solver,
                        struct ode_solver *the_ode_solver, struct ode_solver *the_purkinje_ode_solver,
                        struct grid *the_grid, struct user_options *options);

bool update_ode_state_vector_and_check_for_activity(real_cpu vm_thresold, struct ode_solver *the_ode_solver, struct ode_solver *the_purkinje_ode_solver, struct grid *the_grid);

void set_ode_extra_data(struct config *config, struct grid *the_grid, struct ode_solver *the_ode_solver);
void set_spatial_stim(struct string_voidp_hash_entry *stim_configs, struct grid *the_grid);
                    
void configure_monodomain_solver_from_options(struct monodomain_solver *the_monodomain_solver, struct user_options *options);

bool print_result(const struct grid *the_grid, const struct user_options *configs, int count);

void debug_print_and_leave ();

// NEW FUNCTIONS !
void calculate_activation_time (const real_cpu cur_time, const real_cpu dt, const uint32_t n_active, struct cell_node **ac,\
                        struct ode_solver *the_ode_solver);
void print_activation_time (struct grid *the_grid);
void print_conductivity (struct grid *the_grid);


void linear_system_solver_purkinje (struct config *config, struct grid *the_grid, uint32_t *number_of_iterations, real_cpu *error);
void compute_pmj_current_purkinje_to_tissue (struct ode_solver *the_ode_solver, struct grid *the_grid, struct terminal *the_terminals);
void compute_pmj_current_tissue_to_purkinje (struct ode_solver *the_purkinje_ode_solver, struct grid *the_grid, struct terminal *the_terminals);

void solve_explicit_purkinje_diffusion (struct grid *the_grid, struct ode_solver *the_purkinje_ode_solver, struct monodomain_solver *the_monodomain_solver);

double** allocate_matrix_LU (const uint32_t num_rows, const uint32_t num_cols);
void lu_decomposition (double **lu, uint32_t pivot[], struct cell_node **ac, const uint32_t n_active);
void choose_pivot (double **lu, const uint32_t num_rows, uint32_t *pivot_line, double *Amax, const uint32_t i);
void switch_lines (double **lu, const uint32_t num_cols, uint32_t pivot[], const int pivot_line, const int i);

#endif // MONOALG3D_SOLVER_H
