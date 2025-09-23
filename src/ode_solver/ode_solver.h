//
// Created by sachetto on 02/10/17.
//

#ifndef MONOALG3D_EDO_SOLVER_H
#define MONOALG3D_EDO_SOLVER_H

#include "../common_types/common_types.h"
#include "../config/config_parser.h"
#include <stdbool.h>
#include <stdint.h>

#ifdef __CUDACC__
extern "C" {
#include "../config_helpers/config_helpers.h"
#include "../logger/logger.h"
};
#else
#include "../config_helpers/config_helpers.h"
#include "../logger/logger.h"
#endif

// Forward declaration
struct user_options;
struct ode_solver;

struct cell_model_data {
    int number_of_ode_equations;
    real initial_v;
    char *model_library_path;
};

#define GET_CELL_MODEL_DATA(name) void name(struct cell_model_data *cell_model, bool get_initial_v, bool get_neq)
typedef GET_CELL_MODEL_DATA(get_cell_model_data_fn);

// COMMON FUNCTION
#define SOLVE_MODEL_ODES(name) void name(struct ode_solver *ode_solver, struct string_hash_entry *ode_extra_config, real current_t, real *stim_currents)
typedef SOLVE_MODEL_ODES(solve_model_ode_gpu_fn);
typedef SOLVE_MODEL_ODES(solve_model_ode_cpu_fn);

// CPU FUNCTIONS
#define SET_ODE_INITIAL_CONDITIONS_CPU(name) void name(struct ode_solver *solver, struct string_hash_entry *ode_extra_config)
typedef SET_ODE_INITIAL_CONDITIONS_CPU(set_ode_initial_conditions_cpu_fn);

// GPU FUNCTIONS
#define SET_ODE_INITIAL_CONDITIONS_GPU(name) size_t name(struct ode_solver *solver, struct string_hash_entry *ode_extra_config)
typedef SET_ODE_INITIAL_CONDITIONS_GPU(set_ode_initial_conditions_gpu_fn);

struct ode_solver {

    void *handle;

    real max_dt;
    real min_dt;
    bool auto_dt;

    bool adaptive;
    real rel_tol;
    real abs_tol;

    size_t num_cells_to_solve;
    uint32_t *cells_to_solve;
    uint32_t num_steps;

    bool gpu;
    int gpu_id;

    uint32_t original_num_cells;
    real *sv;
    void *ode_extra_data;
    size_t extra_data_size;
    struct cell_model_data model_data;

    size_t pitch;

    real *ode_dt, *ode_previous_dt, *ode_time_new;

    // User provided functions
    get_cell_model_data_fn *get_cell_model_data;
    set_ode_initial_conditions_cpu_fn *set_ode_initial_conditions_cpu;
    set_ode_initial_conditions_gpu_fn *set_ode_initial_conditions_gpu;
    solve_model_ode_cpu_fn *solve_model_ode_cpu;
    solve_model_ode_gpu_fn *solve_model_ode_gpu;
    // update_gpu_fn_pt update_gpu_fn;
};

#ifdef __cplusplus
extern "C" {
#endif

void set_ode_initial_conditions_for_all_volumes(struct ode_solver *solver, struct string_hash_entry *ode_extra_config);

void update_state_vectors_after_refinement(struct ode_solver *ode_solver, const uint32_t *refined_this_step);
struct ode_solver *new_ode_solver();
void free_ode_solver(struct ode_solver *solver);
void init_ode_solver_with_cell_model(struct ode_solver *solver);

void solve_all_volumes_odes(struct ode_solver *the_ode_solver, real_cpu cur_time, struct string_voidp_hash_entry *stim_configs,
                            struct string_hash_entry *ode_extra_config);

void configure_ode_solver_from_options(struct ode_solver *solver, struct user_options *options);

void configure_purkinje_ode_solver_from_options(struct ode_solver *purkinje_solver, struct user_options *options);
void configure_purkinje_ode_solver_from_ode_solver(struct ode_solver *purkinje_solver, struct ode_solver *solver);

#ifdef __cplusplus
}
#endif

#endif // MONOALG3D_EDO_SOLVER_H
