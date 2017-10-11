//
// Created by sachetto on 02/10/17.
//

#ifndef MONOALG3D_EDO_SOLVER_H
#define MONOALG3D_EDO_SOLVER_H

#include <stdbool.h>
#include <unitypes.h>
#include "../utils/vector/uint32_vector.h"
#include "../models/model_common.h"

#define EULER_METHOD 0
#define EULER_METHOD_ADPT 1

typedef void (*get_cell_model_data_fn_pt)(struct cell_model_data*, bool, bool);

//CPU FUNCTIONS
typedef void (*set_ode_initial_conditions_cpu_fn_pt)(Real *);
typedef void (*solve_model_ode_cpu_fn_pt)(Real, Real *, Real , Real , Real , Real , int , void *);

//GPU FUNCTIONS
typedef void (*set_ode_initial_conditions_gpu_fn_pt)(Real **, int , uint32_t );
typedef void (*solve_model_ode_gpu_fn_pt)(Real, Real *, Real *, uint32_t *, size_t , Real, Real, Real, int, int, void *);


struct ode_solver {

    void *handle;

    Real max_dt;
    Real min_dt;
    Real rel_tol;
    Real abs_tol;

    uint8_t method;

    //used for the adaptive time step solver
    Real previous_dt;
    Real time_new;


    // TODO: create a file for stimulus definition!!
    Real stim_start;
    Real stim_duration;
    Real stim_current;
    uint32_t *cells_to_solve;

    bool gpu;
    int gpu_id;

    Real *sv;
    Real *stim_currents;
    void *edo_extra_data;
    struct cell_model_data model_data;

    //User provided functions
    get_cell_model_data_fn_pt get_cell_model_data_fn;
    set_ode_initial_conditions_cpu_fn_pt set_ode_initial_conditions_cpu_fn;
    set_ode_initial_conditions_gpu_fn_pt set_ode_initial_conditions_gpu_fn;
    solve_model_ode_cpu_fn_pt solve_model_ode_cpu_fn;
    solve_model_ode_gpu_fn_pt solve_model_ode_gpu_fn;


};

void set_ode_initial_conditions_for_all_volumes(struct ode_solver *solver, uint32_t num_volumes);
const char* get_ode_method_name(int met);

void update_state_vectors_after_refinement(Real *sv, uint32_vector *refined_this_step, int neq);
struct ode_solver* new_ode_solver();
void free_ode_solver(struct ode_solver *solver);
void init_ode_solver_with_cell_model(struct ode_solver* solver);
void solve_odes_cpu(struct ode_solver *the_ode_solver, uint64_t  n_active, Real cur_time, int num_steps);
int parse_ode_ini_file(void* user, const char* section, const char* name, const char* value);



#endif //MONOALG3D_EDO_SOLVER_H
