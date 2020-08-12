//
// Created by sachetto on 02/10/17.
//

#include "ode_solver.h"

#include <string.h>
#include <dlfcn.h>
#include <assert.h>
#include "../utils/file_utils.h"
#include "../config/stim_config.h"

#ifdef COMPILE_CUDA
#include "../gpu_utils/gpu_utils.h"
#endif

#include "../3dparty/stb_ds.h"
#include "../config_helpers/config_helpers.h"
#include "../3dparty/sds/sds.h"


struct ode_solver* new_ode_solver() {
    struct ode_solver* result = (struct ode_solver *) malloc(sizeof(struct ode_solver));
    result->sv = NULL;
    result->cells_to_solve = NULL;
    result->handle = NULL;

    result->get_cell_model_data = NULL;
    result->set_ode_initial_conditions_cpu = NULL;
    result->solve_model_ode_cpu = NULL;

    result->set_ode_initial_conditions_gpu = NULL;
    result->solve_model_ode_gpu = NULL;

    result->model_data.initial_v = INFINITY;
    result->model_data.number_of_ode_equations = -1;
    result->model_data.model_library_path = NULL;

    result->ode_extra_data = NULL;
    result->ode_extra_data = 0;
    
    return result;

}

void free_ode_solver(struct ode_solver *solver) {
    if(solver->sv) {
        if(solver->gpu) {
#ifdef COMPILE_CUDA
            cudaFree(solver->sv);
#endif
        }
        else {
            free(solver->sv);
        }
    }

    if(solver->ode_extra_data) {
        free(solver->ode_extra_data);
    }

    if(solver->cells_to_solve) {
        free(solver->cells_to_solve);
    }

    if(solver->model_data.model_library_path) {
        free(solver->model_data.model_library_path);
    }

    if(solver->handle) {
        dlclose(solver->handle);
    }

    free(solver);

}


void init_ode_solver_with_cell_model(struct ode_solver* solver) {

    char *error;

    if(!solver->model_data.model_library_path) {
        fprintf(stderr, "model_library_path not provided. Exiting!\n");
        exit(1);
    }

    solver->handle = dlopen (solver->model_data.model_library_path, RTLD_LAZY);
    if (!solver->handle) {
        fprintf(stderr, "%s\n", dlerror());
        exit(1);
    }

    solver->get_cell_model_data = dlsym(solver->handle, "init_cell_model_data");
    if ((error = dlerror()) != NULL)  {
        fprintf(stderr, "%s\n", error);
        fprintf(stderr, "init_cell_model_data function not found in the provided model library\n");
        if(!isfinite(solver->model_data.initial_v)) {
            fprintf(stderr, "intial_v not provided in the [cell_model] of the config file! Exiting\n");
            exit(1);
        }

    }

    solver->set_ode_initial_conditions_cpu = dlsym(solver->handle, "set_model_initial_conditions_cpu");
    if ((error = dlerror()) != NULL)  {
        fprintf(stderr, "%s\n", error);
        fprintf(stderr, "set_model_initial_conditions function not found in the provided model library\n");
        exit(1);
    }

    solver->solve_model_ode_cpu = dlsym(solver->handle, "solve_model_odes_cpu");
    if ((error = dlerror()) != NULL)  {
        fprintf(stderr, "%s\n", error);
        fprintf(stderr, "solve_model_odes_cpu function not found in the provided model library\n");
        exit(1);
    }

#ifdef COMPILE_CUDA
    solver->set_ode_initial_conditions_gpu = dlsym(solver->handle, "set_model_initial_conditions_gpu");
    if ((error = dlerror()) != NULL)  {
        fputs(error, stderr);
        fprintf(stderr, "set_model_initial_conditions_gpu function not found in the provided model library\n");
        exit(1);
    }

    solver->solve_model_ode_gpu = dlsym(solver->handle, "solve_model_odes_gpu");
    if ((error = dlerror()) != NULL)  {
        fputs(error, stderr);
        fprintf(stderr, "\nsolve_model_odes_gpu function not found in the provided model library\n");
        exit(1);
    }
#endif

}

void set_ode_initial_conditions_for_all_volumes(struct ode_solver *solver, struct string_hash_entry *ode_extra_config) 
{

    bool get_initial_v = !isfinite(solver->model_data.initial_v);
    bool get_neq = solver->model_data.number_of_ode_equations == -1;

    (*(solver->get_cell_model_data))(&(solver->model_data), get_initial_v, get_neq);

    if (solver->gpu) 
    {
#ifdef COMPILE_CUDA

        set_ode_initial_conditions_gpu_fn *soicg_fn_pt = solver->set_ode_initial_conditions_gpu;

        if(!soicg_fn_pt) 
        {
            fprintf(stderr, "The ode solver was set to use the GPU, \n "
                    "but no function called set_model_initial_conditions_gpu "
                    "was provided in the %s shared library file\n", solver->model_data.model_library_path);
            exit(11);
        }

        if(solver->sv != NULL) 
        {
            check_cuda_errors(cudaFree(solver->sv));
        }

        solver->pitch = soicg_fn_pt(solver, ode_extra_config);
#endif
    } 
    else 
    {

        set_ode_initial_conditions_cpu_fn *soicc_fn_pt = solver->set_ode_initial_conditions_cpu;

        if(!soicc_fn_pt) 
        {
            fprintf(stderr, "The ode solver was set to use the CPU, \n "
                    "but no function called set_model_initial_conditions_cpu "
                    "was provided in the %s shared library file\n", solver->model_data.model_library_path);
            exit(11);
        }

        if(solver->sv != NULL) 
        {
            free(solver->sv);
        }

        //We do not malloc here sv anymore. This have to be done in the model solver
        soicc_fn_pt(solver, ode_extra_config);
    }

    assert(solver->sv);
}

void solve_all_volumes_odes(struct ode_solver *the_ode_solver, real_cpu cur_time, struct string_voidp_hash_entry *stim_configs,
                            struct string_hash_entry *ode_extra_config) {

    assert(the_ode_solver->sv);

    size_t n_active = the_ode_solver->num_cells_to_solve;

    real dt = the_ode_solver->min_dt;
    real_cpu time = cur_time;

    real *merged_stims = (real*)calloc(sizeof(real), n_active);

    struct config *tmp = NULL;

	uint32_t i;
    ptrdiff_t n = hmlen(stim_configs);

    uint32_t num_steps = the_ode_solver->num_steps;

    if(stim_configs) {

        real stim_start = 0.0;
        real  stim_dur = 0.0;
        real  stim_period = 0.0;

        for (long k = 0; k < n; k++) {
            tmp = (struct config*) stim_configs[k].value;
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_start, tmp->config_data, "start");
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_dur, tmp->config_data, "duration");
            GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, stim_period, tmp->config_data, "period");

            for (int j = 0; j < num_steps; ++j) {
                if ((time >= stim_start) && (time <= stim_start + stim_dur)) {
                    OMP(parallel for)
                    for (i = 0; i < n_active; i++) {
                        merged_stims[i] += ((real*)(tmp->persistent_data))[i];
                    }
                }
                time += dt;
            }

            if(stim_period > 0.0) {
                if (time >= stim_start + stim_period) {
                    stim_start = stim_start + stim_period;
                    sds stim_start_char = sdscatprintf(sdsempty(), "%lf", stim_start);
                    shput_dup_value(tmp->config_data, "start", stim_start_char);
                    sdsfree(stim_start_char);
                }
            }

            time = cur_time;
        }
    }

    if(the_ode_solver->gpu) {
        #ifdef COMPILE_CUDA
        solve_model_ode_gpu_fn *solve_odes_pt = the_ode_solver->solve_model_ode_gpu;
       // solve_odes_pt(ode_extra_config, cur_time, dt, sv, merged_stims, the_ode_solver->cells_to_solve, n_active, num_steps, extra_data,
       //               extra_data_size);

      solve_odes_pt(the_ode_solver, ode_extra_config, cur_time, merged_stims) ;

        #endif
    }
    else {
        solve_model_ode_cpu_fn *solve_odes_pt = the_ode_solver->solve_model_ode_cpu;
        //solve_odes_pt(ode_extra_config, cur_time, dt, sv, merged_stims, the_ode_solver->cells_to_solve, n_active, num_steps, extra_data);
        solve_odes_pt(the_ode_solver, ode_extra_config, cur_time, merged_stims);
    }

    free(merged_stims);
}

void update_state_vectors_after_refinement(struct ode_solver *ode_solver, const uint32_t *refined_this_step) {

    assert(ode_solver);
    assert(ode_solver->sv);

    size_t num_refined_cells = (size_t) arrlen(refined_this_step)/8;

    real *sv = ode_solver->sv;
    int neq = ode_solver->model_data.number_of_ode_equations;
    real *sv_src;
    real *sv_dst;
    int  i;


    if(ode_solver->gpu) {
        #ifdef COMPILE_CUDA
        size_t pitch_h = ode_solver->pitch;

        OMP(parallel for private(sv_src, sv_dst))
        for (i = 0; i < num_refined_cells; i++) {

            size_t index_id = i * (size_t )8;

            uint32_t index = refined_this_step[index_id];
            sv_src = &sv[index];

            for (int j = 1; j < 8; j++) {
                index = refined_this_step[index_id + j];
                sv_dst = &sv[index];
                check_cuda_errors(cudaMemcpy2D(sv_dst, pitch_h, sv_src, pitch_h, sizeof(real), (size_t )neq, cudaMemcpyDeviceToDevice));
            }
        }
        #endif
    }
    else {

        OMP(parallel for private(sv_src, sv_dst))
        for (i = 0; i < num_refined_cells; i++) {

            size_t index_id = i * (size_t )8;

            uint32_t index = refined_this_step[index_id];
            sv_src = &sv[index * neq];

            for (int j = 1; j < 8; j++) {
                index = refined_this_step[index_id + j];
                sv_dst = &sv[index * neq];
                memcpy(sv_dst, sv_src, neq * sizeof(real));
            }


        }
    }

}

void configure_ode_solver_from_options(struct ode_solver *solver, struct user_options *options) {
    solver->gpu_id = options->gpu_id;
    solver->adaptive = options->ode_adaptive;

    if(solver->adaptive) {
        solver->max_dt = (real)options->dt_pde;
    }

    solver->abs_tol = options->ode_abstol;
    solver->rel_tol = options->ode_reltol;

    solver->min_dt = (real)options->dt_ode;

    solver->gpu = options->gpu;

    if(options->model_file_path) {
        free(solver->model_data.model_library_path);
        solver->model_data.model_library_path = strdup(options->model_file_path);
    }

}

void configure_purkinje_ode_solver_from_options (struct ode_solver *purkinje_solver, struct user_options *options) {

    purkinje_solver->gpu_id = options->purkinje_gpu_id;
    purkinje_solver->adaptive = options->purkinje_ode_adaptive;

    if(purkinje_solver->adaptive) {
        purkinje_solver->max_dt = (real)options->dt_pde;
    }

    purkinje_solver->abs_tol = options->ode_abstol;
    purkinje_solver->rel_tol = options->ode_reltol;

    purkinje_solver->min_dt = (real)options->purkinje_dt_ode;
    
    purkinje_solver->gpu = options->purkinje_gpu;

    if(options->purkinje_model_file_path) 
    {
        free(purkinje_solver->model_data.model_library_path);
        purkinje_solver->model_data.model_library_path = strdup(options->purkinje_model_file_path);
    }
}

void configure_purkinje_ode_solver_from_ode_solver (struct ode_solver *purkinje_solver, struct ode_solver *solver) {

    purkinje_solver->gpu_id = solver->gpu_id;
    purkinje_solver->min_dt = (real)solver->min_dt;
    purkinje_solver->max_dt = (real)solver->max_dt;
    purkinje_solver->gpu = solver->gpu;
    purkinje_solver->adaptive = solver->adaptive;
    purkinje_solver->abs_tol = solver->abs_tol;
    purkinje_solver->rel_tol = solver->rel_tol;

    if(solver->model_data.model_library_path) 
    {
        purkinje_solver->model_data.model_library_path = strdup(solver->model_data.model_library_path);
    }
}

void solve_purkinje_volumes_odes(struct ode_solver *the_ode_solver, real_cpu cur_time,
                                 struct string_voidp_hash_entry *stim_configs,
                                 struct string_hash_entry *ode_extra_config) {


    size_t n_active = the_ode_solver->num_cells_to_solve;
    uint32_t num_steps = the_ode_solver->num_steps;

    assert(the_ode_solver->sv);

    real dt = the_ode_solver->min_dt;
    real_cpu time = cur_time;

    real *merged_stims = (real*)calloc(sizeof(real), n_active);

    struct config *tmp = NULL;

	uint32_t i;
    ptrdiff_t n = hmlen(stim_configs);

    if(stim_configs) 
    {

        real stim_start = 0.0;
        real  stim_dur = 0.0;
        real  stim_period = 0.0;

        for (long k = 0; k < n; k++) 
        {
            tmp = (struct config*) stim_configs[k].value;
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_start, tmp->config_data, "start");
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_dur, tmp->config_data, "duration");
            GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, stim_period, tmp->config_data, "period");

            for (int j = 0; j < num_steps; ++j) 
            {
                if ((time >= stim_start) && (time <= stim_start + stim_dur)) 
                {
                    OMP(parallel for)
                    for (i = 0; i < n_active; i++) 
                    {
                        // This variable should be an accumulator to allow multiple stimulus
                        merged_stims[i] += ((real*)(tmp->persistent_data))[i];
                    }
                }
                time += dt;
            }

            if(stim_period > 0.0) 
            {
                if (time >= stim_start + stim_period) 
                {
                    stim_start = stim_start + stim_period;
                    sds stim_start_char = sdscatprintf(sdsempty(), "%lf", stim_start);
                    shput_dup_value(tmp->config_data, "start", stim_start_char);
                    sdsfree(stim_start_char);
                }
            }

            time = cur_time;
        }
    }

    if(the_ode_solver->gpu) 
    {
        #ifdef COMPILE_CUDA
        solve_model_ode_gpu_fn *solve_odes_pt = the_ode_solver->solve_model_ode_gpu;
//        solve_odes_pt(ode_extra_config, cur_time, dt, sv, merged_stims, the_ode_solver->cells_to_solve, n_active, num_steps, extra_data,
//                      extra_data_size);

        solve_odes_pt(the_ode_solver, ode_extra_config, cur_time, merged_stims);

#endif
    }
    else 
    {
        solve_model_ode_cpu_fn *solve_odes_pt = the_ode_solver->solve_model_ode_cpu;
        //solve_odes_pt(ode_extra_config, cur_time, dt, sv, merged_stims, the_ode_solver->cells_to_solve, n_active, num_steps, extra_data);
        solve_odes_pt(the_ode_solver, ode_extra_config, cur_time, merged_stims);
    }

    free(merged_stims);
}

void init_purkinje_ode_solver_with_cell_model (struct ode_solver* solver, const char *model_library_path) {

    char *error;

    if(model_library_path)
    {
        free(solver->model_data.model_library_path);
        solver->model_data.model_library_path = strdup(model_library_path);
    }

    if(!solver->model_data.model_library_path) 
    {
        fprintf(stderr, "model_library_path not provided. Exiting!\n");
        exit(1);
    }

    solver->handle = dlopen (solver->model_data.model_library_path, RTLD_LAZY);
    if (!solver->handle) 
    {
        fprintf(stderr, "%s\n", dlerror());
        exit(1);
    }

    solver->get_cell_model_data = dlsym(solver->handle, "init_cell_model_data");
    if ((error = dlerror()) != NULL)  
    {
        fprintf(stderr, "%s\n", error);
        fprintf(stderr, "init_cell_model_data function not found in the provided model library\n");
        if(!isfinite(solver->model_data.initial_v)) 
        {
            fprintf(stderr, "intial_v not provided in the [cell_model] of the config file! Exiting\n");
            exit(1);
        }

    }

    solver->set_ode_initial_conditions_cpu = dlsym(solver->handle, "set_model_initial_conditions_cpu");
    if ((error = dlerror()) != NULL)  
    {
        fprintf(stderr, "%s\n", error);
        fprintf(stderr, "set_model_initial_conditions function not found in the provided model library\n");
        exit(1);
    }

    solver->solve_model_ode_cpu = dlsym(solver->handle, "solve_model_odes_cpu");
    if ((error = dlerror()) != NULL)  
    {
        fprintf(stderr, "%s\n", error);
        fprintf(stderr, "solve_model_odes_cpu function not found in the provided model library\n");
        exit(1);
    }

#ifdef COMPILE_CUDA
    solver->set_ode_initial_conditions_gpu = dlsym(solver->handle, "set_model_initial_conditions_gpu");
    if ((error = dlerror()) != NULL)  
    {
        fputs(error, stderr);
        fprintf(stderr, "set_model_initial_conditions_gpu function not found in the provided model library\n");
        exit(1);
    }

    solver->solve_model_ode_gpu = dlsym(solver->handle, "solve_model_odes_gpu");
    if ((error = dlerror()) != NULL)  
    {
        fputs(error, stderr);
        fprintf(stderr, "\nsolve_model_odes_gpu function not found in the provided model library\n");
        exit(1);
    }
#endif

}
