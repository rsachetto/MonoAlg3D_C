//
// Created by sachetto on 02/10/17.
//

#include "ode_solver.h"

#include <assert.h>
#include <dlfcn.h>
#include <string.h>

#if defined(COMPILE_CUDA) || defined(COMPILE_SYCL)
#define COMPILE_GPU
#endif

#ifdef COMPILE_GPU
#include "../gpu_utils/accel_utils.h"
#endif

#include "../3dparty/sds/sds.h"
#include "../3dparty/stb_ds.h"
#include "../config_helpers/config_helpers.h"

struct ode_solver *new_ode_solver() {
    struct ode_solver *result = (struct ode_solver *)malloc(sizeof(struct ode_solver));
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
    result->extra_data_size = 0;

    result->auto_dt = false;

    return result;
}

void free_ode_solver(struct ode_solver *solver) {
    if(solver->sv) {
        if(solver->gpu) {
#ifdef COMPILE_GPU
            free_device(solver->sv);
#endif
        } else {
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

void init_ode_solver_with_cell_model(struct ode_solver *solver) {

    char *error;

    if(!solver->model_data.model_library_path) {
        fprintf(stderr, "model_library_path not provided. Exiting!\n");
        exit(1);
    }

    solver->handle = dlopen(solver->model_data.model_library_path, RTLD_LAZY);
    if(!solver->handle) {
        fprintf(stderr, "%s\n", dlerror());
        exit(1);
    }

    solver->get_cell_model_data = (get_cell_model_data_fn *)dlsym(solver->handle, "init_cell_model_data");
    if((error = dlerror()) != NULL) {
        fprintf(stderr, "%s\n", error);
        fprintf(stderr, "init_cell_model_data function not found in the provided model library\n");
        if(!isfinite(solver->model_data.initial_v)) {
            fprintf(stderr, "intial_v not provided in the [cell_model] of the config file! Exiting\n");
            exit(1);
        }
    }

    solver->set_ode_initial_conditions_cpu = (set_ode_initial_conditions_cpu_fn *)dlsym(solver->handle, "set_model_initial_conditions_cpu");
    if((error = dlerror()) != NULL) {
        fprintf(stderr, "%s\n", error);
        fprintf(stderr, "set_model_initial_conditions function not found in the provided model library\n");
        exit(1);
    }

    solver->solve_model_ode_cpu = (solve_model_ode_cpu_fn *)dlsym(solver->handle, "solve_model_odes_cpu");
    if((error = dlerror()) != NULL) {
        fprintf(stderr, "%s\n", error);
        fprintf(stderr, "solve_model_odes_cpu function not found in the provided model library\n");
        exit(1);
    }

#ifdef COMPILE_GPU
    solver->set_ode_initial_conditions_gpu = (set_ode_initial_conditions_gpu_fn *)dlsym(solver->handle, "set_model_initial_conditions_gpu");
    if((error = dlerror()) != NULL) {
        fputs(error, stderr);
        fprintf(stderr, "set_model_initial_conditions_gpu function not found in the provided model library\n");
        exit(1);
    }

    solver->solve_model_ode_gpu = (solve_model_ode_gpu_fn *)dlsym(solver->handle, "solve_model_odes_gpu");
    if((error = dlerror()) != NULL) {
        fputs(error, stderr);
        fprintf(stderr, "\nsolve_model_odes_gpu function not found in the provided model library\n");
        exit(1);
    }
#endif
}

void set_ode_initial_conditions_for_all_volumes(struct ode_solver *solver, struct string_hash_entry *ode_extra_config) {

    bool get_initial_v = !isfinite(solver->model_data.initial_v);
    bool get_neq = solver->model_data.number_of_ode_equations == -1;

    (*(solver->get_cell_model_data))(&(solver->model_data), get_initial_v, get_neq);

    if(solver->gpu) {
#ifdef COMPILE_GPU

        set_ode_initial_conditions_gpu_fn *soicg_fn_pt = solver->set_ode_initial_conditions_gpu;

        if(!soicg_fn_pt) {
            fprintf(stderr,
                    "The ode solver was set to use the GPU, \n "
                    "but no function called set_model_initial_conditions_gpu "
                    "was provided in the %s shared library file\n",
                    solver->model_data.model_library_path);
            exit(11);
        }

        if(solver->sv != NULL) {
            free_device(solver->sv);
            solver->sv = NULL;
        }

        solver->pitch = soicg_fn_pt(solver, ode_extra_config);
#endif
    } else {

        set_ode_initial_conditions_cpu_fn *soicc_fn_pt = solver->set_ode_initial_conditions_cpu;

        if(!soicc_fn_pt) {
            fprintf(stderr,
                    "The ode solver was set to use the CPU, \n "
                    "but no function called set_model_initial_conditions_cpu "
                    "was provided in the %s shared library file\n",
                    solver->model_data.model_library_path);
            exit(11);
        }

        if(solver->sv != NULL) {
            free(solver->sv);
        }

        // We do not malloc here sv anymore. This have to be done in the model solver
        soicc_fn_pt(solver, ode_extra_config);
    }

    if(solver->sv == NULL) {
        log_error_and_exit("Error allocating memory for the ODE's state vector. Exiting!\n");
    }
}

void solve_all_volumes_odes(struct ode_solver *the_ode_solver, real_cpu cur_time, struct string_voidp_hash_entry *stim_configs,
                            struct string_hash_entry *ode_extra_config) {

    assert(the_ode_solver->sv);

    size_t n_active = the_ode_solver->num_cells_to_solve;

    real dt = the_ode_solver->min_dt;
    real_cpu time = cur_time;

    real *merged_stims = (real *)calloc(sizeof(real), n_active);

    struct config *tmp = NULL;

    uint32_t i;
    size_t n = hmlen(stim_configs);
    uint32_t num_steps = the_ode_solver->num_steps;

    if(stim_configs) {

        for(size_t k = 0; k < n; k++) {

            real stim_start = 0.0;
            real stim_dur = 0.0;
            real stim_period = 0.0;

            tmp = (struct config *)stim_configs[k].value;
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_start, tmp, "start");
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_dur, tmp, "duration");
            GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, stim_period, tmp, "period");

            for(int j = 0; j < num_steps; ++j) {
                if((time >= stim_start) && (time <= stim_start + stim_dur)) {
                    OMP(parallel for)
                    for(i = 0; i < n_active; i++) {
                        merged_stims[i] += ((real *)(tmp->persistent_data))[i];
                    }
                }
                time += dt;
            }

            if(stim_period > 0.0) {
                if(time >= stim_start + stim_period) {
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
#ifdef COMPILE_GPU
        solve_model_ode_gpu_fn *solve_odes_fn = the_ode_solver->solve_model_ode_gpu;
        solve_odes_fn(the_ode_solver, ode_extra_config, cur_time, merged_stims);
#endif
    } else {
        solve_model_ode_cpu_fn *solve_odes_fn = the_ode_solver->solve_model_ode_cpu;
        solve_odes_fn(the_ode_solver, ode_extra_config, cur_time, merged_stims);
    }

    free(merged_stims);
}

void update_state_vectors_after_refinement(struct ode_solver *ode_solver, const uint32_t *refined_this_step) {

    assert(ode_solver);
    assert(ode_solver->sv);

    size_t num_refined_cells = (size_t)arrlen(refined_this_step) / 8;

    real *sv = ode_solver->sv;
    int neq = ode_solver->model_data.number_of_ode_equations;
    real *sv_src;
    real *sv_dst;
    int i;

    if(ode_solver->adaptive) {
        neq += 3;
    }

    const size_t max_index = 8;

    if(ode_solver->gpu) {
#ifdef COMPILE_GPU
        size_t pitch_h = ode_solver->pitch;

        for(i = 0; i < num_refined_cells; i++) {

            size_t index_id = i * max_index;

            uint32_t index = refined_this_step[index_id];
            sv_src = &sv[index];

            for(int j = 1; j < max_index; j++) {
                index = refined_this_step[index_id + j];
                sv_dst = &sv[index];
                memcpy2d_device(sv_dst, pitch_h, sv_src, pitch_h, sizeof(real), (size_t)neq, DEVICE_TO_DEVICE);
            }
        }
#endif
    } else {

        OMP(parallel for private(sv_src, sv_dst))
        for(i = 0; i < num_refined_cells; i++) {

            size_t index_id = i * max_index;

            uint32_t index = refined_this_step[index_id];
            sv_src = &sv[index * neq];

            for(int j = 1; j < max_index; j++) {
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
    solver->auto_dt = options->auto_dt_ode;

    if(solver->adaptive) {
        solver->max_dt = (real)options->dt_pde;

        if(solver->auto_dt || (options->dt_ode == 0.0)) {
            real min_dt = 1e-10;

            // This is highly unlikely
            if(min_dt > solver->max_dt) {
                min_dt = min_dt / 1.1;
            }

            solver->min_dt = min_dt;

        } else {
            solver->min_dt = (real)options->dt_ode;
        }
        solver->abs_tol = options->ode_abstol;
        solver->rel_tol = options->ode_reltol;

    } else {
        if(options->dt_ode == 0.0) {
            solver->min_dt = 0.01;
        } else {
            solver->min_dt = (real)options->dt_ode;
        }
    }

    solver->gpu = options->gpu;

    if(options->model_file_path) {
        free(solver->model_data.model_library_path);
        solver->model_data.model_library_path = strdup(options->model_file_path);
    }
}

void configure_purkinje_ode_solver_from_options(struct ode_solver *purkinje_solver, struct user_options *options) {

    purkinje_solver->gpu_id = options->purkinje_gpu_id;
    purkinje_solver->adaptive = options->purkinje_ode_adaptive;

    if(purkinje_solver->adaptive) {
        purkinje_solver->max_dt = (real)options->dt_pde;
    }

    purkinje_solver->abs_tol = options->ode_abstol;
    purkinje_solver->rel_tol = options->ode_reltol;
    purkinje_solver->min_dt = (real)options->purkinje_dt_ode;
    purkinje_solver->gpu = options->purkinje_gpu;

    if(options->purkinje_model_file_path) {
        free(purkinje_solver->model_data.model_library_path);
        purkinje_solver->model_data.model_library_path = strdup(options->purkinje_model_file_path);
    }
}

void configure_purkinje_ode_solver_from_ode_solver(struct ode_solver *purkinje_solver, struct ode_solver *solver) {

    purkinje_solver->gpu_id = solver->gpu_id;
    purkinje_solver->min_dt = (real)solver->min_dt;
    purkinje_solver->max_dt = (real)solver->max_dt;
    purkinje_solver->gpu = solver->gpu;
    purkinje_solver->adaptive = solver->adaptive;
    purkinje_solver->abs_tol = solver->abs_tol;
    purkinje_solver->rel_tol = solver->rel_tol;

    if(solver->model_data.model_library_path) {
        purkinje_solver->model_data.model_library_path = strdup(solver->model_data.model_library_path);
    }
}
