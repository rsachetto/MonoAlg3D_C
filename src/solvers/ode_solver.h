//
// Created by sachetto on 02/10/17.
//

#ifndef MONOALG3D_EDO_SOLVER_H
#define MONOALG3D_EDO_SOLVER_H

#include <stdbool.h>
#include <unitypes.h>
#include "../models/model_common.h"

#define EULER_METHOD 0
#define EULER_METHOD_ADPT 1

struct ode_solver {

    Real max_dt;
    Real min_dt;
    Real rel_tol;
    Real abs_tol;
    char *model_library_path;

    uint8_t method;

    //used for the adaptive time step solver
    Real previous_dt;
    Real time_new;


    // TODO: create a file for stimulus definition!!
    Real stim_start;
    Real stim_duration;
    Real stim_current;
    uint64_t *cells_to_solve;

    bool gpu;
    int gpu_id;

    Real *sv;
    Real *stim_currents;
    void *edo_extra_data;
    struct cell_model_data model_data;

    //User provided functions
    void (*get_cell_model_data_fn)(struct cell_model_data*);
    void (*set_ode_initial_conditions_fn)(Real *);

    //Use dinamic libraries to load from a .so model file
    //https://www.dwheeler.com/program-library/Program-Library-HOWTO/x172.html
    //https://www.cprogramming.com/tutorial/shared-libraries-linux-gcc.html


};

void set_ode_initial_conditions_for_all_volumes(struct ode_solver *solver, uint64_t num_volumes);
const char* get_ode_method_name(int met);

struct ode_solver* new_ode_solver(const char *model_library_path);
void init_ode_solver_with_cell_model(struct ode_solver* solver);
void solve_odes_cpu(struct ode_solver *the_ode_solver, uint64_t  n_active, Real cur_time, int num_steps);



#endif //MONOALG3D_EDO_SOLVER_H
