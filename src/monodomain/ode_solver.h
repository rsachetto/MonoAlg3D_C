//
// Created by sachetto on 02/10/17.
//

#ifndef MONOALG3D_EDO_SOLVER_H
#define MONOALG3D_EDO_SOLVER_H

#include <stdbool.h>
#include <stdint.h>
#include "config/stim_config_hash.h"
#include "config/config_parser.h"

#include "../models_library/model_common.h"

struct ode_solver {

    void *handle;

    real max_dt;
    real min_dt;
    real rel_tol;
    real abs_tol;

    //used for the adaptive time step solver
    real previous_dt;
    real time_new;

    uint32_t *cells_to_solve;

    bool gpu;
    int gpu_id;

    real *sv;
    void *edo_extra_data;
    size_t extra_data_size;
    struct cell_model_data model_data;

    size_t pitch;

    //User provided functions
    get_cell_model_data_fn *get_cell_model_data;
    set_ode_initial_conditions_cpu_fn *set_ode_initial_conditions_cpu;
    set_ode_initial_conditions_gpu_fn *set_ode_initial_conditions_gpu;
    solve_model_ode_cpu_fn *solve_model_ode_cpu;
    solve_model_ode_gpu_fn *solve_model_ode_gpu;
    //update_gpu_fn_pt update_gpu_fn;


};

void set_ode_initial_conditions_for_all_volumes(struct ode_solver *solver, uint32_t num_volumes);

void update_state_vectors_after_refinement(struct ode_solver *ode_solver, const uint32_t *refined_this_step);
struct ode_solver* new_ode_solver();
void free_ode_solver(struct ode_solver *solver);
void init_ode_solver_with_cell_model(struct ode_solver* solver);
void solve_all_volumes_odes(struct ode_solver *the_ode_solver, uint32_t n_active, double cur_time, int num_steps,
                            struct stim_config_hash *stim_configs);
void configure_ode_solver_from_options(struct ode_solver *solver, struct user_options *options);


#endif //MONOALG3D_EDO_SOLVER_H
