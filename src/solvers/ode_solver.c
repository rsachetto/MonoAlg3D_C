//
// Created by sachetto on 02/10/17.
//

#include "ode_solver.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dlfcn.h>

struct ode_solver* new_ode_solver(const char *model_library_path) {
    struct ode_solver* result = (struct ode_solver *) malloc(sizeof(struct ode_solver));
    result->sv = NULL;
    result->stim_currents = NULL;
    result->edo_extra_data = NULL;
    result->cells_to_solve = NULL;
    result->handle = NULL;

    result->get_cell_model_data_fn = NULL;
    result->set_ode_initial_conditions_fn = NULL;
    result->solve_model_ode_cpu_fn = NULL;

    result->model_library_path = strdup(model_library_path);
    init_ode_solver_with_cell_model(result);
    return result;
}

void free_ode_solver(struct ode_solver *solver) {
    if(solver->sv) {
        free(solver->sv);
    }

    if(solver->stim_currents) {
        free(solver->stim_currents);
    }

    if(solver->edo_extra_data) {
        free(solver->edo_extra_data);
    }

    if(solver->cells_to_solve) {
        free(solver->cells_to_solve);
    }

    if(solver->model_library_path) {
        free(solver->model_library_path);
    }

    if(solver->handle) {
        dlclose(solver->handle);
    }

    free(solver);

}


void init_ode_solver_with_cell_model(struct ode_solver* solver) {

    char *error;
    void *handle = solver->handle;

    printf("Opening %s as model lib\n", solver->model_library_path);

    handle = dlopen (solver->model_library_path, RTLD_LAZY);
    if (!handle) {
        fputs (dlerror(), stderr);
        exit(1);
    }

    solver->get_cell_model_data_fn = dlsym(handle, "init_cell_model_data");
    if ((error = dlerror()) != NULL)  {
        fputs(error, stderr);
        fprintf(stderr, "init_cell_model_data function not found in the provided model library");
        exit(1);
    }

    solver->set_ode_initial_conditions_fn = dlsym(handle, "set_model_initial_conditions");
    if ((error = dlerror()) != NULL)  {
        fputs(error, stderr);
        fprintf(stderr, "set_model_initial_conditions function not found in the provided model library");
        exit(1);
    }

    solver->solve_model_ode_cpu_fn = dlsym(handle, "solve_model_ode_cpu");
    if ((error = dlerror()) != NULL)  {
        fputs(error, stderr);
        fprintf(stderr, "solve_model_ode_cpu function not found in the provided model library");
        exit(1);
    }


    //TODO: we need to think about this handler stuff
    //dlclose(handle);

}

void set_ode_initial_conditions_for_all_volumes(struct ode_solver *solver, uint64_t num_cells) {

    (*(solver->get_cell_model_data_fn))(&(solver->model_data));

    set_ode_initial_conditions_fn_pt soi_fn_pt = solver->set_ode_initial_conditions_fn;

    if (solver->gpu) {
        if (solver->method == EULER_METHOD) {
            //TODO: @Incomplete
            //pitch = setIC_ode_gpu(&sv, originalNumCells);
        } else {
        //TODO: @Incomplete
        //pitch = setIC_ode_gpu_adapt(&sv, originalNumCells, dt_edo);
        }
    } else {
        int n_odes = solver->model_data.number_of_ode_equations;

        if(solver->sv != NULL) {
            free(solver->sv);
        }

        solver->sv = (Real*)malloc(n_odes*num_cells*sizeof(Real));
        for(u_int64_t i = 0; i < num_cells; i++) {
            soi_fn_pt(solver->sv + (i*n_odes));
        }
    }
}

void solve_odes_cpu(struct ode_solver *the_ode_solver, uint64_t  n_active, Real cur_time, int num_steps) {

    solve_model_ode_cpu_fn_pt solve_odes_pt = the_ode_solver->solve_model_ode_cpu_fn;
    uint64_t sv_id;

    Real dt = the_ode_solver->min_dt;
    int n_odes = the_ode_solver->model_data.number_of_ode_equations;
    Real *stims = the_ode_solver->stim_currents;
    Real *sv = the_ode_solver->sv;

    Real stim_start = the_ode_solver->stim_start;
    Real stim_dur = the_ode_solver->stim_duration;
    void *extra_data = the_ode_solver->edo_extra_data;

    for(int i = 0; i < n_active; i++) {
        sv_id = the_ode_solver->cells_to_solve[i];

        for (int j = 0; j < num_steps; ++j) {
            solve_odes_pt(dt, sv + (sv_id*n_odes), stims[i], stim_start, stim_dur, cur_time, n_odes, extra_data);
            cur_time += dt;
        }

    }

}


const char* get_ode_method_name(int met) {

    switch(met) {
        case 0: return "Euler Method";
        case 1: return "Euler Method with adaptive time step (Formula)";
        default: printf("Invalid Method!!\n");  exit(0);
    }
}