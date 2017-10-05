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

    result->get_cell_model_data_fn = NULL;
    result->model_library_path = strdup(model_library_path);
    init_ode_solver_with_cell_model(result);
    return result;
}

void init_ode_solver_with_cell_model(struct ode_solver* solver) {

    void *handle;
    char *error;

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

    //TODO: we need to think about this handler stuff
    //dlclose(handle);

}

void set_ode_initial_conditions_for_all_volumes(struct ode_solver *solver, uint64_t num_cells) {

    (*(solver->get_cell_model_data_fn))(&(solver->model_data));


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

        solver->sv = (Real*)malloc(n_odes*num_cells);
        for(u_int64_t i = 0; i < num_cells; i++) {
            (*(solver->set_ode_initial_conditions_fn))(solver->sv + (i*n_odes));
        }
    }
}

void solve_odes_cpu(struct ode_solver *the_ode_solver, uint64_t  n_active, Real cur_time, int num_steps) {

}


const char* get_ode_method_name(int met) {

    switch(met) {
        case 0: return "Euler Method";
        case 1: return "Euler Method with adaptive time step (Formula)";
        default: printf("Invalid Method!!\n");  exit(0);
    }
}