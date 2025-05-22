//
// Created by sachetto on 03/10/17.
//

#ifdef COMPILE_GUI
#include "../gui/gui.h"
#endif

#if defined(COMPILE_CUDA) || defined(COMPILE_SYCL)
#define COMPILE_GPU
#endif

#ifdef COMPILE_GPU
#include "../gpu_utils/accel_utils.h"
#endif

#include "../3dparty/stb_ds.h"
#include "../config/config_common.h"
#include "../config/modify_current_domain_config.h"
#include "../config/stim_config.h"
#include "../libraries_common/common_data_structures.h"
#include "../save_mesh_library/save_mesh_helper.h"
#include "../utils/file_utils.h"
#include "../utils/stop_watch.h"
#include "monodomain_solver.h"
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <unistd.h>

#define MAX_TMP 1024

struct monodomain_solver *new_monodomain_solver() {

    struct monodomain_solver *result = (struct monodomain_solver *)malloc(sizeof(struct monodomain_solver));

    result->beta = 0.14;
    result->cm = 1.0;
    result->only_abort_after_dt = 0;

    return result;
}

int solve_monodomain(struct monodomain_solver *the_monodomain_solver, struct ode_solver *the_ode_solver, struct grid *the_grid, struct user_options *configs,
                     struct gui_shared_info *gui_config) {

    assert(configs);
    assert(the_grid);
    assert(the_monodomain_solver);
    assert(the_ode_solver);

    log_msg(LOG_LINE_SEPARATOR);

    uint64_t ode_total_time = 0, cg_total_time = 0, total_write_time = 0, total_mat_time = 0, total_ref_time = 0, total_deref_time = 0, cg_partial,
             total_ecg_time = 0, total_order_time = 0;
    uint64_t total_update_sv_time = 0;
    uint64_t total_update_cells_time = 0;

    uint64_t purkinje_ode_total_time = 0, purkinje_cg_total_time = 0, purkinje_cg_partial = 0;

    uint32_t total_cg_it = 0, purkinje_total_cg_it = 0;

    struct stop_watch stop_watch, solver_time;

    init_stop_watch(&stop_watch);
    init_stop_watch(&solver_time);

    ///////MAIN CONFIGURATION BEGIN//////////////////
    init_ode_solver_with_cell_model(the_ode_solver);

    struct string_voidp_hash_entry *stimuli_configs = configs->stim_configs;
    struct string_voidp_hash_entry *modify_domain_configs = configs->modify_domain_configs;
    struct string_voidp_hash_entry *purkinje_stimuli_configs = configs->purkinje_stim_configs;
    struct config *extra_data_config = configs->extra_data_config;
    struct config *purkinje_extra_data_config = configs->purkinje_extra_data_config;
    struct config *domain_config = configs->domain_config;
    struct config *purkinje_config = configs->purkinje_config;
    struct config *assembly_matrix_config = configs->assembly_matrix_config;
    struct config *linear_system_solver_config = configs->linear_system_solver_config;
    struct config *purkinje_linear_system_solver_config = configs->purkinje_linear_system_solver_config;
    struct config *save_mesh_config = configs->save_mesh_config;
    struct config *save_state_config = configs->save_state_config;
    struct config *restore_state_config = configs->restore_state_config;
    struct config *update_monodomain_config = configs->update_monodomain_config;
    struct config *calc_ecg_config = configs->calc_ecg_config;

    bool has_extra_data = (extra_data_config != NULL);
    bool has_purkinje_extra_data = (purkinje_extra_data_config != NULL);

    real_cpu last_stimulus_time = -1.0;
    bool has_any_periodic_stim = false;

    real_cpu dt_pde = the_monodomain_solver->dt;
    real_cpu finalT = the_monodomain_solver->final_time;
    real_cpu dt_ode = the_ode_solver->min_dt;

    size_t num_stims = shlen(stimuli_configs);
    if(num_stims) {
        // Init all stimuli
        STIM_CONFIG_HASH_FOR_INIT_FUNCTIONS(stimuli_configs);

        // Find last stimuli
        real_cpu s_end;
        real_cpu stim_start = 0.0;
        real_cpu stim_duration = 0.0;
        real_cpu stim_period = 0;

        for(size_t i = 0; i < num_stims; i++) {

            struct config *sconfig = (struct config *)stimuli_configs[i].value;

            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, stim_start, sconfig, "start");
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, stim_duration, sconfig, "duration");
            GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, stim_period, sconfig, "period");

            s_end = stim_start + stim_duration;

            has_any_periodic_stim |= (bool)(stim_period > 0.0);

            if(s_end > last_stimulus_time) {
                last_stimulus_time = s_end;
            }
        }
    }

    size_t num_purkinje_stims = shlen(purkinje_stimuli_configs);
    if(num_purkinje_stims) {
        // Init all stimuli
        STIM_CONFIG_HASH_FOR_INIT_FUNCTIONS(purkinje_stimuli_configs);

        // Find last stimuli
        real_cpu s_end;
        real_cpu stim_start = 0.0;
        real_cpu stim_duration = 0.0;
        real_cpu stim_period = 0;

        for(size_t i = 0; i < num_purkinje_stims; i++) {

            struct config *sconfig = (struct config *)purkinje_stimuli_configs[i].value;

            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, stim_start, sconfig, "start");
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, stim_duration, sconfig, "duration");
            GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, stim_period, sconfig, "period");

            s_end = stim_start + stim_duration;

            has_any_periodic_stim |= (bool)(stim_period > 0.0);

            if(s_end > last_stimulus_time) {
                last_stimulus_time = s_end;
            }
        }
    }

    real_cpu modify_at_end = 0.0;

    size_t num_modify_domains = shlen(modify_domain_configs);

    if(num_modify_domains) {
        MODIFY_DOMAIN_CONFIG_HASH_FOR_INIT_FUNCTIONS(modify_domain_configs);

        real_cpu modify_at = 0.0;

        for(size_t i = 0; i < num_modify_domains; i++) {

            struct config *dconfig = (struct config *)modify_domain_configs[i].value;

            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, modify_at, dconfig, "modify_after_dt");

            if(modify_at > modify_at_end)
                modify_at_end = modify_at;
        }
    }

    if(!purkinje_config && !domain_config) {
        log_error_and_exit("Error configuring the domain! No Purkinje or tissue configuration was provided!\n");
    }

    // Configure the functions and set the Purkinje mesh domain
    if(purkinje_config) {
        init_config_functions(purkinje_config, "shared_libs/libdefault_purkinje.so", "purkinje");
    }

    // Configure the functions and set the mesh domain
    if(domain_config) {
        init_config_functions(domain_config, "./shared_libs/libdefault_domains.so", "domain");
    }

    if(assembly_matrix_config != NULL) {
        init_config_functions(assembly_matrix_config, "./shared_libs/libdefault_matrix_assembly.so", "assembly_matrix");
    } else {
        log_error_and_exit("No assembly matrix configuration provided! Exiting!\n");
    }

    if(linear_system_solver_config) {
        init_config_functions(linear_system_solver_config, "./shared_libs/libdefault_linear_system_solver.so", "linear_system_solver");
    } else {
        log_error_and_exit("No linear solver configuration provided! Exiting!\n");
    }

    if(purkinje_linear_system_solver_config) {
        init_config_functions(purkinje_linear_system_solver_config, "./shared_libs/libdefault_linear_system_solver.so", "purkinje_linear_system_solver");
    }

    bool calc_ecg = (calc_ecg_config != NULL);
    int calc_ecg_rate = 1;

    if(calc_ecg) {
        init_config_functions(calc_ecg_config, "./shared_libs/libdefault_calc_ecg.so", "calc_ecg");
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, calc_ecg_rate, calc_ecg_config, "calc_rate");
        calc_ecg &= (calc_ecg_rate > 0);
    }

    int print_rate = 0;
    int output_print_rate = 1;
    char *out_dir_name = strdup("./");

    bool save_to_file = (save_mesh_config != NULL);
    real_cpu start_saving_after_dt = 0.0;

    if(save_to_file) {
        init_config_functions(save_mesh_config, "./shared_libs/libdefault_save_mesh.so", "save_result");
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, print_rate, save_mesh_config, "print_rate");
        GET_PARAMETER_STRING_VALUE_OR_USE_DEFAULT(out_dir_name, save_mesh_config, "output_dir");
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, start_saving_after_dt, save_mesh_config, "start_saving_after_dt");
        save_to_file &= (print_rate > 0) && (out_dir_name);

        if(print_rate > 0)
            output_print_rate = print_rate;

    } else {
        log_info("No configuration provided to save the results! The results will not be saved!\n");
        free(out_dir_name);
    }

    bool save_checkpoint = (save_state_config != NULL);

    if(save_checkpoint && the_grid->adaptive) {
        log_info("Saving checkpoint is not implemented for adaptive grids yet!\n");
        save_checkpoint = false;
    }

    if(save_checkpoint) {
        init_config_functions(save_state_config, "./shared_libs/libdefault_save_state.so", "save_state");
    } else {
        log_info("No configuration provided to make simulation checkpoints! Chekpoints will not be created!\n");
    }

    bool restore_checkpoint = (restore_state_config != NULL);
    char *restore_in_dir_name = NULL;

    if(restore_checkpoint && the_grid->adaptive) {
        log_warn("Restoring checkpoint is not implemented for adaptive grids yet!\n");
        restore_checkpoint = false;
    }

    if(restore_state_config) {
        init_config_functions(restore_state_config, "./shared_libs/libdefault_restore_state.so", "restore_state");
        GET_PARAMETER_STRING_VALUE_OR_USE_DEFAULT(restore_in_dir_name, restore_state_config, "input_dir");
        if(restore_in_dir_name == NULL) {
            restore_in_dir_name = out_dir_name;
        }
    }

    if(has_extra_data) {
        init_config_functions(extra_data_config, "./shared_libs/libdefault_extra_data.so", "extra_data");
    }

    if(has_purkinje_extra_data) {
        init_config_functions(purkinje_extra_data_config, "./shared_libs/libdefault_extra_data.so", "extra_data");
    }

    log_msg(LOG_LINE_SEPARATOR);

    bool restore_success = false;

    struct time_info time_info = {0.0, finalT, dt_pde, 0};

    if(restore_checkpoint) {
        // Here we only restore the monodomain_solver_state...
        restore_success = ((restore_state_fn *)restore_state_config->main_function)(&time_info, restore_state_config, save_mesh_config, NULL,
                                                                                    the_monodomain_solver, NULL, NULL, restore_in_dir_name);
    }

    // HACK: we have to restore the last_t time info as the restore state function changes it to a wrong value
    time_info.final_t = finalT;

    if(update_monodomain_config) {
        init_config_functions(update_monodomain_config, "./shared_libs/libdefault_update_monodomain.so", "update_monodomain");
    } else {
        log_error_and_exit("No update monodomain configuration provided! Exiting!\n");
    }
    ///////MAIN CONFIGURATION END//////////////////

    int refine_each = the_monodomain_solver->refine_each;
    int derefine_each = the_monodomain_solver->derefine_each;

    bool redo_matrix = false;

    bool activity;

    real_cpu refinement_bound = the_monodomain_solver->refinement_bound;
    real_cpu derefinement_bound = the_monodomain_solver->derefinement_bound;

    bool adaptive = the_grid->adaptive;
    real_cpu start_adpt_at = the_monodomain_solver->start_adapting_at;

    bool calc_retropropagation = true;

#ifdef COMPILE_GUI
    bool show_gui = configs->show_gui;
    if(show_gui) {
        gui_config->grid_info.alg_grid = the_grid;
        gui_config->simulating = true;
        gui_config->paused = !configs->start_visualization_unpaused;
    }
#endif

#ifdef COMPILE_GPU
    bool linear_solver_on_gpu = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(linear_solver_on_gpu, linear_system_solver_config, "use_gpu");

    bool init_gpu = linear_solver_on_gpu || the_ode_solver->gpu;

    if(init_gpu) {
        char *device_info = get_device_info(the_ode_solver, linear_solver_on_gpu);
        log_info("%s\n", device_info);
        free(device_info);
    }
#endif

    int success;
    struct ode_solver *the_purkinje_ode_solver = NULL;

    if(purkinje_config) {

        // Allocate a new 'ode_solver' for the Purkinje
        the_purkinje_ode_solver = new_ode_solver();

        // Here we configure the Purkinje ode_solver using the [purkinje_ode_solver] parameters
        // If there is no [purkinje_ode_solver] section we configure the Purkinje ODE solver using the input from the [ode_solver] section
        if(!domain_config) {
            // Only Purkinje simulation
            configure_purkinje_ode_solver_from_ode_solver(the_purkinje_ode_solver, the_ode_solver);
        } else {
            // Purkinje + Tissue simulation
            configure_purkinje_ode_solver_from_options(the_purkinje_ode_solver, configs);
        }

        init_ode_solver_with_cell_model(the_purkinje_ode_solver);

        success = ((set_spatial_purkinje_fn *)purkinje_config->main_function)(purkinje_config, the_grid, the_purkinje_ode_solver);
        if(!success) {
            log_error_and_exit("Error configuring the Purkinje domain!\n");
        }

        calc_retropropagation = the_grid->purkinje->network->calc_retropropagation;
    }

    if(domain_config) {
        success = ((set_spatial_domain_fn *)domain_config->main_function)(domain_config, the_grid);

        if(!success) {
            log_error_and_exit("Error configuring the tissue domain!\n");
        }
    }

    if(restore_checkpoint) {
        restore_success &= ((restore_state_fn *)restore_state_config->main_function)(&time_info, restore_state_config, save_mesh_config, the_grid, NULL, NULL,
                                                                                     NULL, restore_in_dir_name);
    }

    real_cpu start_dx, start_dy, start_dz;
    real_cpu max_dx, max_dy, max_dz;

    start_dx = start_dy = start_dz = 100.0;
    max_dx = max_dy = max_dz = 100.0;

    uint32_t original_num_cells = 0;
    uint32_t original_num_purkinje_cells = 0;

    if(domain_config) {

        start_dx = the_grid->start_discretization.x;
        start_dy = the_grid->start_discretization.y;
        start_dz = the_grid->start_discretization.z;

        max_dx = the_grid->max_discretization.x;
        max_dy = the_grid->max_discretization.y;
        max_dz = the_grid->max_discretization.z;

        // This is used only to print information about the domain:
        char tmp[MAX_TMP];

        // TODO: change how this is handled
        if(!shget(domain_config->config_data, "start_dx")) {
            snprintf(tmp, MAX_TMP, "%lf", start_dx);
            shput_dup_value(domain_config->config_data, "start_dx", tmp);
        }

        if(!shget(domain_config->config_data, "start_dy")) {
            snprintf(tmp, MAX_TMP, "%lf", start_dy);
            shput_dup_value(domain_config->config_data, "start_dy", tmp);
        }

        if(!shget(domain_config->config_data, "start_dz")) {
            snprintf(tmp, MAX_TMP, "%lf", start_dz);
            shput_dup_value(domain_config->config_data, "start_dz", tmp);
        }

        order_grid_cells(the_grid);

        original_num_cells = the_grid->num_active_cells;
        the_ode_solver->original_num_cells = original_num_cells;
        the_ode_solver->num_cells_to_solve = original_num_cells;
    }

    // Purkinje section
    if(purkinje_config) {
        original_num_purkinje_cells = the_grid->purkinje->number_of_purkinje_cells;
        the_purkinje_ode_solver->original_num_cells = original_num_purkinje_cells;
        the_purkinje_ode_solver->num_cells_to_solve = original_num_purkinje_cells;
    }

    save_old_cell_positions(the_grid);

    if(adaptive) {
        update_cells_to_solve(the_grid, the_ode_solver);
    }

    struct terminal *the_terminals = NULL;
    if(domain_config && purkinje_config) {

        log_info("Start - link_purkinje_to_tissue\n");
        the_terminals = link_purkinje_to_tissue(the_grid);
        log_info("End - link_purkinje_to_tissue\n");
    }

    if(has_extra_data) {
        the_ode_solver->ode_extra_data =
            ((set_extra_data_fn *)extra_data_config->main_function)(&time_info, extra_data_config, the_grid, &(the_ode_solver->extra_data_size));

        if(the_ode_solver->ode_extra_data == NULL) {
            log_warn("set_extra_data function was called but returned NULL!\n");
        } else {
            if(the_ode_solver->extra_data_size == 0) {
                log_warn("set_extra_data function was called but extra_data_size is 0!\n");
                log_warn("Maybe you forgot to call the SET_EXTRA_DATA_SIZE(size) macro in your extra_data_function!\n");
            }
        }
    }

    if(has_purkinje_extra_data) {
        the_purkinje_ode_solver->ode_extra_data = ((set_extra_data_fn *)purkinje_extra_data_config->main_function)(
            &time_info, purkinje_extra_data_config, the_grid, &(the_purkinje_ode_solver->extra_data_size));

        if(the_purkinje_ode_solver->ode_extra_data == NULL) {
            log_warn("set_extra_data function was called but returned NULL!\n");
        } else {
            if(the_purkinje_ode_solver->extra_data_size == 0) {
                log_warn("set_extra_data function was called but extra_data_size is 0!\n");
                log_warn("Maybe you forgot to call the SET_EXTRA_DATA_SIZE(size) macro in your extra_data_function!\n");
            }
        }
    }

    log_info("Setting ODE's initial conditions\n");

    if(domain_config) {
        set_ode_initial_conditions_for_all_volumes(the_ode_solver, configs->ode_extra_config);
    }
    if(purkinje_config) {
        set_ode_initial_conditions_for_all_volumes(the_purkinje_ode_solver, configs->purkinje_ode_extra_config);
    }

    // We need to call this function after because of the pitch.... maybe we have to change the way
    // we pass this parameters to the cell model....
    if(restore_checkpoint) {
        restore_success &= ((restore_state_fn *)restore_state_config->main_function)(&time_info, restore_state_config, save_mesh_config, NULL, NULL,
                                                                                     the_ode_solver, the_purkinje_ode_solver, restore_in_dir_name);
    }

    real_cpu initial_v, purkinje_initial_v = 0;
    initial_v = the_ode_solver->model_data.initial_v;

    if(purkinje_config) {
        purkinje_initial_v = the_purkinje_ode_solver->model_data.initial_v;
    }

    print_solver_info(the_monodomain_solver, the_ode_solver, the_purkinje_ode_solver, the_grid, configs);

    the_ode_solver->num_steps = 1;

    if(!the_ode_solver->adaptive) {
        if(dt_pde >= dt_ode) {
            the_ode_solver->num_steps = (int)(dt_pde / dt_ode);
            log_info("Solving EDO %d times before solving PDE\n", the_ode_solver->num_steps);
        } else {
            log_info("WARNING: EDO time step is greater than PDE time step. Adjusting EDP time step %lf to EDO time step %lf\n", dt_pde, dt_ode);
            dt_pde = dt_ode;
        }
    }

    if(purkinje_config) {
        the_purkinje_ode_solver->num_steps = 1;

        if(!the_purkinje_ode_solver->adaptive) {
            if(dt_pde >= dt_ode) {
                the_purkinje_ode_solver->num_steps = (int)(dt_pde / dt_ode);
                log_info("Solving Purkinje EDO %d times before solving PDE\n", the_purkinje_ode_solver->num_steps);
            } else {
                log_info("WARNING: EDO time step is greater than PDE time step. Adjusting EDP time step %lf to EDO time step %lf\n", dt_pde, dt_ode);
                dt_pde = dt_ode;
            }
        }
    }

    fflush(stdout);

    start_stop_watch(&stop_watch);

    if(!restore_checkpoint || !restore_success) {
        if(assembly_matrix_config->init_function) {
            ((set_pde_initial_condition_fn *)assembly_matrix_config->init_function)(assembly_matrix_config, the_monodomain_solver, the_ode_solver, the_grid,
                                                                                    initial_v, purkinje_initial_v);
        } else {
            log_error_and_exit("Function for the Monodomain initial conditions not provided (init_function on [assembly_matrix] section)!\n");
        }
    }

    ((assembly_matrix_fn *)assembly_matrix_config->main_function)(assembly_matrix_config, the_monodomain_solver, the_grid);

    total_mat_time = stop_stop_watch(&stop_watch);

    start_stop_watch(&solver_time);

    int save_state_rate = 0;
    char *save_checkpoint_out_dir = NULL;

    if(save_checkpoint) {
        GET_PARAMETER_STRING_VALUE_OR_USE_DEFAULT(save_checkpoint_out_dir, save_state_config, "output_dir");
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, save_state_rate, save_state_config, "save_rate");

        if(save_checkpoint_out_dir == NULL) {
            save_checkpoint_out_dir = out_dir_name;
        }
    }

    real_cpu vm_threshold = configs->vm_threshold;

    bool abort_on_no_activity = the_monodomain_solver->abort_on_no_activity;

    real_cpu solver_error = 0, purkinje_solver_error = 0;
    uint32_t solver_iterations = 0, purkinje_solver_iterations = 0;

    if(num_stims > 0) {
        set_spatial_stim(&time_info, stimuli_configs, the_grid, false);
    }

    if(num_purkinje_stims > 0) {
        set_spatial_stim(&time_info, purkinje_stimuli_configs, the_grid, true);
    }

    real_cpu cur_time = time_info.current_t;

    struct stop_watch iteration_time_watch;
    uint64_t iteration_time;

    CALL_INIT_LINEAR_SYSTEM(linear_system_solver_config, the_grid, false || !domain_config);
    CALL_INIT_SAVE_MESH(save_mesh_config);
    CALL_INIT_CALC_ECG(calc_ecg_config, the_ode_solver, the_grid, out_dir_name);

    if(purkinje_linear_system_solver_config) {
        CALL_INIT_LINEAR_SYSTEM(purkinje_linear_system_solver_config, the_grid, true);
    }

#ifdef COMPILE_GUI
    if(show_gui) {
        gui_config->grid_info.loaded = true;
    }
#endif

    real_cpu only_abort_after_dt = the_monodomain_solver->only_abort_after_dt;

    if(out_dir_name) {
        sds buffer_ini = sdscatfmt(sdsempty(), "%s/original_configuration.ini", out_dir_name);
        options_to_ini_file(configs, buffer_ini);
        sdsfree(buffer_ini);
    }

    int count = time_info.iteration;

    init_stop_watch(&iteration_time_watch);

    log_info("Starting simulation\n");

    // Main simulation loop start
    while(cur_time - finalT <= dt_pde) {

        start_stop_watch(&iteration_time_watch);

        time_info.current_t = cur_time;
        time_info.iteration = count;

#ifdef COMPILE_GUI
        if(show_gui) {
            omp_set_nest_lock(&gui_config->sleep_lock);
            if(gui_config->restart) {

                gui_config->time = 0.0f;

                CALL_END_LINEAR_SYSTEM(linear_system_solver_config);
                CALL_END_SAVE_MESH(save_mesh_config, the_grid);
                CALL_END_LINEAR_SYSTEM(purkinje_linear_system_solver_config);
                CALL_END_CALC_ECG(calc_ecg_config);

                return RESTART_SIMULATION;
            }
            if(gui_config->exit) {
                CALL_END_LINEAR_SYSTEM(linear_system_solver_config);
                CALL_END_SAVE_MESH(save_mesh_config, the_grid);
                CALL_END_LINEAR_SYSTEM(purkinje_linear_system_solver_config);
                CALL_END_CALC_ECG(calc_ecg_config);

                return END_SIMULATION;
            }
        }
#endif

        if(save_to_file && (count % print_rate == 0) && (cur_time >= start_saving_after_dt)) {
            start_stop_watch(&stop_watch);
            ((save_mesh_fn *)save_mesh_config->main_function)(&time_info, save_mesh_config, the_grid, the_ode_solver, the_purkinje_ode_solver);
            total_write_time += stop_stop_watch(&stop_watch);
        }

        if(calc_ecg && (count % calc_ecg_rate == 0)) {
            start_stop_watch(&stop_watch);
            ((calc_ecg_fn *)calc_ecg_config->main_function)(&time_info, calc_ecg_config, the_grid, out_dir_name);
            total_ecg_time += stop_stop_watch(&stop_watch);
        }

        if(cur_time > 0.0) {
            activity = update_ode_state_vector_and_check_for_activity(vm_threshold, the_ode_solver, the_purkinje_ode_solver, the_grid);

            if(abort_on_no_activity && cur_time > last_stimulus_time && cur_time > only_abort_after_dt) {
                if(!activity) {
                    log_info("No activity, aborting simulation\n");
                    break;
                }
            }
        }

        if(purkinje_config) {
            start_stop_watch(&stop_watch);

            // REACTION: Purkinje
            solve_all_volumes_odes(the_purkinje_ode_solver, cur_time, purkinje_stimuli_configs, configs->purkinje_ode_extra_config);

            purkinje_ode_total_time += stop_stop_watch(&stop_watch);

            start_stop_watch(&stop_watch);

            // UPDATE: Purkinje
            ((update_monodomain_fn *)update_monodomain_config->main_function)(&time_info, update_monodomain_config, the_grid, the_monodomain_solver,
                                                                              the_grid->purkinje->num_active_purkinje_cells, the_grid->purkinje->purkinje_cells,
                                                                              the_purkinje_ode_solver, original_num_purkinje_cells);

            purkinje_ode_total_time += stop_stop_watch(&stop_watch);

            start_stop_watch(&stop_watch);

            // TODO: show the purkinje fibers in the visualization tool
            //            #ifdef COMPILE_GUI
            //            if (show_gui) {
            //                omp_set_lock(&gui_config->draw_lock);
            //            }
            //            #endif

            // COUPLING: Calculate the PMJ current from the Tissue to the Purkinje
            if(domain_config && calc_retropropagation)
                compute_pmj_current_tissue_to_purkinje(the_purkinje_ode_solver, the_grid, the_terminals);

            // DIFUSION: Purkinje
            if(purkinje_linear_system_solver_config) // Purkinje-coupled
                ((linear_system_solver_fn *)purkinje_linear_system_solver_config->main_function)(
                    &time_info, purkinje_linear_system_solver_config, the_grid, the_grid->purkinje->num_active_purkinje_cells,
                    the_grid->purkinje->purkinje_cells, &purkinje_solver_iterations, &purkinje_solver_error);
            else // Only-Purkinje
                ((linear_system_solver_fn *)linear_system_solver_config->main_function)(
                    &time_info, linear_system_solver_config, the_grid, the_grid->purkinje->num_active_purkinje_cells, the_grid->purkinje->purkinje_cells,
                    &purkinje_solver_iterations, &purkinje_solver_error);

            purkinje_cg_partial = stop_stop_watch(&stop_watch);

            purkinje_cg_total_time += purkinje_cg_partial;

            purkinje_total_cg_it += purkinje_solver_iterations;
        }

        if(domain_config) {

            start_stop_watch(&stop_watch);

            // REACTION
            solve_all_volumes_odes(the_ode_solver, cur_time, stimuli_configs, configs->ode_extra_config);
            ((update_monodomain_fn *)update_monodomain_config->main_function)(&time_info, update_monodomain_config, the_grid, the_monodomain_solver,
                                                                              the_grid->num_active_cells, the_grid->active_cells, the_ode_solver,
                                                                              original_num_cells);

            ode_total_time += stop_stop_watch(&stop_watch);

            start_stop_watch(&stop_watch);

#ifdef COMPILE_GUI
            if(show_gui) {
                omp_set_nest_lock(&gui_config->draw_lock);
            }
#endif

            // COUPLING: Calculate the PMJ current from the Purkinje to the Tissue
            if(purkinje_config) {
                compute_pmj_current_purkinje_to_tissue(the_ode_solver, the_grid, the_terminals);
            }

            // DIFUSION: Tissue
            ((linear_system_solver_fn *)linear_system_solver_config->main_function)(
                &time_info, linear_system_solver_config, the_grid, the_grid->num_active_cells, the_grid->active_cells, &solver_iterations, &solver_error);
            if(isnan(solver_error)) {
                log_error("\nSimulation stopped due to NaN on time %lf. This is probably a problem with the cellular model solver.\n.", cur_time);
#ifdef COMPILE_GUI
                if(show_gui) {
                    omp_unset_nest_lock(&gui_config->draw_lock);
                }
#endif
                return SIMULATION_FINISHED;
            }

            cg_partial = stop_stop_watch(&stop_watch);

            cg_total_time += cg_partial;

            total_cg_it += solver_iterations;
        }

        if(count % output_print_rate == 0) {
            if(purkinje_config && domain_config) {
                log_info("t = %.5lf, Iterations = "
                         "%" PRIu32 ", Error Norm = %e, Number of Tissue Cells:"
                         "%" PRIu32 ", Tissue CG Iterations time: %ld us\n"
                         "            , Iterations = "
                         "%" PRIu32 ", Error Norm = %e, Number of Purkinje Cells:"
                         "%" PRIu32 ", Purkinje CG Iterations time: %ld us",
                         cur_time, solver_iterations, solver_error, the_grid->num_active_cells, cg_partial, purkinje_solver_iterations, purkinje_solver_error,
                         the_grid->purkinje->num_active_purkinje_cells, purkinje_cg_partial);
            } else if(domain_config) {
                log_info("t = %lf, Iterations = "
                         "%" PRIu32 ", Error Norm = %e, Number of Cells:"
                         "%" PRIu32 ", CG Iterations time: %ld us",
                         cur_time, solver_iterations, solver_error, the_grid->num_active_cells, cg_partial);
            } else {
                log_info("t = %lf, Iterations = "
                         "%" PRIu32 ", Error Norm = %e, Number of Purkinje Cells:"
                         "%" PRIu32 ", Purkinje CG Iterations time: %ld us",
                         cur_time, purkinje_solver_iterations, purkinje_solver_error, the_grid->purkinje->num_active_purkinje_cells, purkinje_cg_partial);
            }
        }

        if(adaptive) {
            redo_matrix = false;
            if(cur_time >= start_adpt_at) {
                if(count % refine_each == 0) {
                    start_stop_watch(&stop_watch);
                    redo_matrix = refine_grid_with_bound(the_grid, refinement_bound, start_dx, start_dy, start_dz);
                    total_ref_time += stop_stop_watch(&stop_watch);
                }

                if(count % derefine_each == 0) {
                    start_stop_watch(&stop_watch);
                    redo_matrix |= derefine_grid_with_bound(the_grid, derefinement_bound, max_dx, max_dy, max_dz);
                    total_deref_time += stop_stop_watch(&stop_watch);
                }

                if(redo_matrix) {

                    start_stop_watch(&stop_watch);
                    order_grid_cells(the_grid);
                    total_order_time += stop_stop_watch(&stop_watch);

                    if(stimuli_configs) {
                        if(cur_time <= last_stimulus_time || has_any_periodic_stim) {
                            set_spatial_stim(&time_info, stimuli_configs, the_grid, false);
                        }
                    }
                    if(has_extra_data) {
                        free(the_ode_solver->ode_extra_data);
                        the_ode_solver->ode_extra_data = ((set_extra_data_fn *)extra_data_config->main_function)(&time_info, extra_data_config, the_grid,
                                                                                                                 &(the_ode_solver->extra_data_size));
                    }

                    start_stop_watch(&stop_watch);
                    update_cells_to_solve(the_grid, the_ode_solver);
                    total_update_cells_time += stop_stop_watch(&stop_watch);

                    start_stop_watch(&stop_watch);
                    if(arrlen(the_grid->refined_this_step) > 0) {
                        update_state_vectors_after_refinement(the_ode_solver, the_grid->refined_this_step);
                    }
                    total_update_sv_time += stop_stop_watch(&stop_watch);

                    start_stop_watch(&stop_watch);
                    ((assembly_matrix_fn *)assembly_matrix_config->main_function)(assembly_matrix_config, the_monodomain_solver, the_grid);
                    total_mat_time += stop_stop_watch(&stop_watch);

                    // MAPPING: Update the mapping between the Purkinje mesh and the refined/derefined grid
                    if(purkinje_config && domain_config) {
                        update_link_purkinje_to_endocardium(the_grid, the_terminals);
                    }

                    // Freeing resources and reconfiguring the linear system solver
                    CALL_END_LINEAR_SYSTEM(linear_system_solver_config);
                    CALL_INIT_LINEAR_SYSTEM(linear_system_solver_config, the_grid, false);
                }
            }
        }

        if(num_modify_domains) {

            for(size_t i = 0; i < num_modify_domains; i++) {

                bool modification_applied = false;

                struct config *dconfig = (struct config *)modify_domain_configs[i].value;

                // TODO: if only one modification can be applied this does not make sense.
                // We will need a boolean for each modification
                GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(modification_applied, dconfig, "modification_applied");

                if(!modification_applied) {

                    real_cpu modify_at = 0.0;
                    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, modify_at, dconfig, "modify_after_dt");

                    if(cur_time >= modify_at) {

                        ((modify_current_domain_fn *)dconfig->main_function)(&time_info, dconfig, the_grid);

                        order_grid_cells(the_grid);

                        if(stimuli_configs) {
                            if(cur_time <= last_stimulus_time || has_any_periodic_stim) {
                                set_spatial_stim(&time_info, stimuli_configs, the_grid, false);
                            }
                        }

                        if(has_extra_data) {
                            CALL_FREE_EXTRA_DATA(extra_data_config, the_ode_solver->ode_extra_data);
                            the_ode_solver->ode_extra_data = ((set_extra_data_fn *)extra_data_config->main_function)(&time_info, extra_data_config, the_grid,
                                                                                                                     &(the_ode_solver->extra_data_size));
                        }

                        if(arrlen(the_grid->refined_this_step) > 0) {
                            update_state_vectors_after_refinement(the_ode_solver, the_grid->refined_this_step);
                        }

                        update_cells_to_solve(the_grid, the_ode_solver);

                        start_stop_watch(&stop_watch);
                        ((assembly_matrix_fn *)assembly_matrix_config->main_function)(assembly_matrix_config, the_monodomain_solver, the_grid);
                        total_mat_time += stop_stop_watch(&stop_watch);

                        CALL_END_LINEAR_SYSTEM(linear_system_solver_config);
                        CALL_INIT_LINEAR_SYSTEM(linear_system_solver_config, the_grid, false);

                        CALL_END_SAVE_MESH(save_mesh_config, the_grid);
                        CALL_INIT_SAVE_MESH(save_mesh_config);

                        shput_dup_value(dconfig->config_data, "modification_applied", "true");
                    }
                }
            }
        }

#ifdef COMPILE_GUI
        if(configs->show_gui) {
            omp_unset_nest_lock(&gui_config->draw_lock);
            gui_config->time = (float)cur_time;
        }
#endif
        count++;
        cur_time += dt_pde;

        if(save_checkpoint) {
            if(save_state_rate && count != 0 && (count % save_state_rate == 0)) {
                time_info.iteration = count;
                time_info.current_t = cur_time;
                printf("Saving state with time = %lf, and count = %d\n", time_info.current_t, time_info.iteration);
                ((save_state_fn *)save_state_config->main_function)(&time_info, save_state_config, save_mesh_config, the_grid, the_monodomain_solver,
                                                                    the_ode_solver, the_purkinje_ode_solver, save_checkpoint_out_dir);
            }
        }

        iteration_time = stop_stop_watch(&iteration_time_watch);

        if((count - 1) % output_print_rate == 0) {
            log_msg(", Total Iteration time: %ld us\n", iteration_time);
        }
    }

    // if no save_rate is passed we only save at the end of the simulation;
    if(save_checkpoint && save_state_rate == 0) {
        time_info.iteration = count;
        time_info.current_t = cur_time;
        printf("Saving state with time = %lf, and count = %d\n", time_info.current_t, time_info.iteration);
        ((save_state_fn *)save_state_config->main_function)(&time_info, save_state_config, save_mesh_config, the_grid, the_monodomain_solver, the_ode_solver,
                                                            the_purkinje_ode_solver, save_checkpoint_out_dir);
    }

    uint64_t res_time = stop_stop_watch(&solver_time);

    double conv_rate = 1000.0 * 1000.0 * 60.0;
    log_info("Resolution Time: %ld μs (%lf min)\n", res_time, res_time / conv_rate);

    if(calc_ecg) {
        log_info("ECG calculation Time: %ld μs (%lf min)\n", total_ecg_time, total_ecg_time / conv_rate);
    }

    if(domain_config) {
        log_info("Total Write Time: %ld μs (%lf min)\n", total_write_time, total_write_time / conv_rate);
        log_info("ODE Total Time: %ld μs (%lf min)\n", ode_total_time, ode_total_time / conv_rate);
        log_info("CG Total Time: %ld μs (%lf min)\n", cg_total_time, cg_total_time / conv_rate);

        if(adaptive) {
            log_info("Assemble matrix time: %ld μs (%lf min)\n", total_mat_time, total_mat_time / conv_rate);
            log_info("Refine time: %ld μs (%lf min)\n", total_ref_time, total_ref_time / conv_rate);
            log_info("Derefine time: %ld μs (%lf min)\n", total_deref_time, total_deref_time / conv_rate);
            log_info("Order grid time: %ld μs (%lf min)\n", total_order_time, total_order_time / conv_rate);
            log_info("Update cells time: %ld μs (%lf min)\n", total_update_cells_time, total_update_cells_time / conv_rate);
            log_info("Update SV time: %ld μs (%lf min)\n", total_update_sv_time, total_update_sv_time / conv_rate);
        }

        log_info("CG Total Iterations: %u\n", total_cg_it);

        uint64_t u_time = res_time - (total_ecg_time + total_write_time + ode_total_time + cg_total_time + total_mat_time + total_ref_time + total_deref_time +
                                      total_order_time + total_update_sv_time + total_update_cells_time);

        log_info("Unmeasured time: %ld μs (%lf min)\n", u_time, u_time / conv_rate);
    }

    if(purkinje_config) {
        log_info("Purkinje ODE Total Time: %ld μs (%lf min)\n", purkinje_ode_total_time, purkinje_ode_total_time / conv_rate);
        log_info("Purkinje CG Total Time: %ld μs (%lf min)\n", purkinje_cg_total_time, purkinje_cg_total_time / conv_rate);
        log_info("Purkinje CG Total Iterations: %u\n", purkinje_total_cg_it);
    }

    if(purkinje_config && domain_config) {
        write_pmj_delay(the_grid, save_mesh_config, the_terminals);
        write_terminals_info(the_grid, save_mesh_config, the_terminals);
        free_terminals(the_terminals, the_grid->purkinje->network->number_of_terminals);
    }

#ifdef COMPILE_GUI
    if(show_gui) {
        gui_config->solver_time = res_time;
        gui_config->ode_total_time = ode_total_time;
        gui_config->cg_total_time = cg_total_time;
        gui_config->total_mat_time = total_mat_time;
        gui_config->total_ref_time = total_ref_time;
        gui_config->total_deref_time = total_deref_time;
        gui_config->total_write_time = total_write_time;
        gui_config->total_cg_it = total_cg_it;
        gui_config->simulating = false;
    }
#endif

    CALL_END_LINEAR_SYSTEM(linear_system_solver_config);
    CALL_END_SAVE_MESH(save_mesh_config, the_grid);
    CALL_END_CALC_ECG(calc_ecg_config);
    if(purkinje_linear_system_solver_config)
        CALL_END_LINEAR_SYSTEM(purkinje_linear_system_solver_config);

    return SIMULATION_FINISHED;
}

void set_spatial_stim(struct time_info *time_info, struct string_voidp_hash_entry *stim_configs, struct grid *the_grid, bool purkinje) {

    struct config *tmp = NULL;
    size_t n = shlen(stim_configs);
    for(size_t i = 0; i < n; i++) {
        tmp = (struct config *)stim_configs[i].value;
        ((set_spatial_stim_fn *)tmp->main_function)(time_info, tmp, the_grid, purkinje);
    }
}

bool update_ode_state_vector_and_check_for_activity(real_cpu vm_threshold, struct ode_solver *the_ode_solver, struct ode_solver *the_purkinje_ode_solver,
                                                    struct grid *the_grid) {
    bool act = false;

    // Tissue section
    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    if(the_ode_solver) {
        int n_odes = the_ode_solver->model_data.number_of_ode_equations;

        real *sv = the_ode_solver->sv;

        if(the_ode_solver->gpu) {
#ifdef COMPILE_GPU
            uint32_t max_number_of_cells = the_ode_solver->original_num_cells;
            real *vms;
            size_t mem_size = max_number_of_cells * sizeof(real);

            vms = (real *)malloc(mem_size);

            if(the_grid->adaptive) {
                // check_cuda_error(cudaMemcpy(vms, sv, mem_size, cudaMemcpyDeviceToHost));
                memcpy_device(vms, sv, mem_size, DEVICE_TO_HOST);
            }

            OMP(parallel for)
            for(uint32_t i = 0; i < n_active; i++) {
                vms[ac[i]->sv_position] = (real)ac[i]->v;

                if(ac[i]->v > vm_threshold) {
                    act = true;
                }
            }

            // check_cuda_error(cudaMemcpy(sv, vms, mem_size, cudaMemcpyHostToDevice));
            memcpy_device(sv, vms, mem_size, HOST_TO_DEVICE);
            free(vms);
#endif
        } else {
            OMP(parallel for)
            for(uint32_t i = 0; i < n_active; i++) {
                sv[ac[i]->sv_position * n_odes] = (real)ac[i]->v;

                if(ac[i]->v > vm_threshold) {
                    act = true;
                }
            }
        }
    }

    if(the_purkinje_ode_solver) {
        // Purkinje section
        uint32_t n_active_purkinje = the_grid->purkinje->number_of_purkinje_cells;
        struct cell_node **ac_purkinje = the_grid->purkinje->purkinje_cells;

        int n_odes_purkinje = the_purkinje_ode_solver->model_data.number_of_ode_equations;

        real *sv_purkinje = the_purkinje_ode_solver->sv;

        if(the_purkinje_ode_solver->gpu) {
#ifdef COMPILE_GPU
            uint32_t max_number_of_purkinje_cells = the_purkinje_ode_solver->original_num_cells;
            real *vms_purkinje;
            size_t mem_size_purkinje = max_number_of_purkinje_cells * sizeof(real);

            vms_purkinje = (real *)malloc(mem_size_purkinje);

            if(the_grid->adaptive) {
                memcpy_device(vms_purkinje, sv_purkinje, mem_size_purkinje, DEVICE_TO_HOST);
            }

            OMP(parallel for)
            for(uint32_t i = 0; i < n_active_purkinje; i++) {
                vms_purkinje[ac_purkinje[i]->sv_position] = (real)ac_purkinje[i]->v;

                if(ac_purkinje[i]->v > vm_threshold) {
                    act = true;
                }
            }

            memcpy_device(sv_purkinje, vms_purkinje, mem_size_purkinje, HOST_TO_DEVICE);
            free(vms_purkinje);
#endif
        } else {
            OMP(parallel for)
            for(uint32_t i = 0; i < n_active_purkinje; i++) {
                sv_purkinje[ac_purkinje[i]->sv_position * n_odes_purkinje] = (real)ac_purkinje[i]->v;

                if(ac_purkinje[i]->v > vm_threshold) {
                    act = true;
                }
            }
        }
    }

    return act;
}

void save_old_cell_positions(struct grid *the_grid) {

    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    uint32_t i;

    OMP(parallel for)
    for(i = 0; i < n_active; i++) {
        ac[i]->sv_position = ac[i]->grid_position;
    }

    // Purkinje section
    struct grid_purkinje *the_purkinje = the_grid->purkinje;

    if(the_purkinje && the_purkinje->first_cell) {
        uint32_t n_purkinje_active = the_purkinje->num_active_purkinje_cells;
        struct cell_node **ac_purkinje = the_purkinje->purkinje_cells;
        OMP(parallel for)
        for(i = 0; i < n_purkinje_active; i++) {
            ac_purkinje[i]->sv_position = ac_purkinje[i]->grid_position;
        }
    }
}

void update_cells_to_solve(struct grid *the_grid, struct ode_solver *solver) {

    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    if(solver->cells_to_solve) {
        free(solver->cells_to_solve);
    }

    solver->num_cells_to_solve = n_active;
    solver->cells_to_solve = (uint32_t *)malloc(n_active * sizeof(uint32_t));
    uint32_t *cts = solver->cells_to_solve;
    int i;

    OMP(parallel for)
    for(i = 0; i < n_active; i++) {
        cts[i] = ac[i]->sv_position;
    }
}

void print_solver_info(struct monodomain_solver *the_monodomain_solver, struct ode_solver *the_ode_solver, struct ode_solver *the_purkinje_ode_solver,
                       struct grid *the_grid, struct user_options *options) {

    log_msg(LOG_LINE_SEPARATOR);

    log_info("System parameters: \n");

    log_msg(LOG_LINE_SEPARATOR);

#if defined(_OPENMP)
    log_info("[main] Using OpenMP with %d threads\n", omp_get_max_threads());
#endif

    log_info("[monodomain_solver] Beta = %.10lf, Cm = %.10lf\n", the_monodomain_solver->beta, the_monodomain_solver->cm);
    log_info("[monodomain_solver] PDE time step = %.10lf\n", the_monodomain_solver->dt);
    log_info("[monodomain_solver] ODE min time step = %.10lf\n", the_ode_solver->min_dt);
    log_info("[monodomain_solver] Simulation Final Time = %lf\n", the_monodomain_solver->final_time);

    log_msg(LOG_LINE_SEPARATOR);

    if(the_ode_solver->gpu) {
        log_info("[ode_solver] Using GPU to solve ODEs\n");
    }

    log_info("[ode_solver] Using %s as model lib\n", the_ode_solver->model_data.model_library_path);
    log_info("[ode_solver] Initial V: %lf\n", the_ode_solver->model_data.initial_v);
    log_info("[ode_solver] Number of ODEs in cell model: %d\n", the_ode_solver->model_data.number_of_ode_equations);

    size_t len = shlen(options->ode_extra_config);

    if(len) {
        log_msg(LOG_LINE_SEPARATOR);
        STRING_HASH_PRINT_KEY_VALUE_LOG("[ode_solver]", options->ode_extra_config);
    }

    if(the_purkinje_ode_solver) {
        log_msg(LOG_LINE_SEPARATOR);

        if(the_purkinje_ode_solver->gpu) {
            log_info("[purkinje_ode_solver] Using GPU to solve ODEs\n");
        }

        log_info("[purkinje_ode_solver] Using %s as model lib\n", the_purkinje_ode_solver->model_data.model_library_path);
        log_info("[purkinje_ode_solver] Initial V: %lf\n", the_purkinje_ode_solver->model_data.initial_v);
        log_info("[purkinje_ode_solver] Number of ODEs in cell model: %d\n", the_purkinje_ode_solver->model_data.number_of_ode_equations);
    }

    len = shlen(options->purkinje_ode_extra_config);

    if(len) {
        log_msg(LOG_LINE_SEPARATOR);
        STRING_HASH_PRINT_KEY_VALUE_LOG("[purkinje_ode_solver]", options->purkinje_ode_extra_config);
    }

    log_msg(LOG_LINE_SEPARATOR);

    log_info("[grid] Initial N. of Elements = "
             "%" PRIu32 "\n",
             the_grid->num_active_cells);

    if(the_grid->adaptive) {
        log_info("Using adaptativity\n");
        log_info("[monodomain_solver] Refinement Bound = %lf\n", the_monodomain_solver->refinement_bound);
        log_info("[monodomain_solver] Derefinement Bound = %lf\n", the_monodomain_solver->derefinement_bound);
        log_info("[monodomain_solver] Refining each %d time steps\n", the_monodomain_solver->refine_each);
        log_info("[monodomain_solver] Derefining each %d time steps\n", the_monodomain_solver->derefine_each);

        log_info("[domain] Domain maximum Space Discretization: dx %lf um, dy %lf um, dz %lf um\n", the_grid->max_discretization.x,
                 the_grid->max_discretization.y, the_grid->max_discretization.z);

        log_info("[monodomain_solver] The adaptivity will start in time: %lf ms\n", the_monodomain_solver->start_adapting_at);
    }

    log_msg(LOG_LINE_SEPARATOR);

    if(options->save_state_config) {
        print_save_state_config_values(options->save_state_config);
        log_msg(LOG_LINE_SEPARATOR);
    }

    if(options->linear_system_solver_config) {
        print_linear_system_solver_config_values(options->linear_system_solver_config);
        log_msg(LOG_LINE_SEPARATOR);
    }

    if(options->purkinje_linear_system_solver_config) {
        print_purkinje_linear_system_solver_config_values(options->purkinje_linear_system_solver_config);
        log_msg(LOG_LINE_SEPARATOR);
    }

    if(options->save_mesh_config) {
        print_save_mesh_config_values(options->save_mesh_config);
        log_msg(LOG_LINE_SEPARATOR);
    }

    if(options->stim_configs) {

        size_t num_stims = shlen(options->stim_configs);

        for(int i = 0; i < num_stims; i++) {

            struct string_voidp_hash_entry e = options->stim_configs[i];
            log_info("Stimulus name: %s\n", e.key);
            struct config *tmp = (struct config *)e.value;
            print_stim_config_values(tmp);
            log_msg(LOG_LINE_SEPARATOR);
        }
    }

    if(options->purkinje_stim_configs) {

        size_t num_stims = shlen(options->purkinje_stim_configs);

        for(int i = 0; i < num_stims; i++) {

            struct string_voidp_hash_entry e = options->purkinje_stim_configs[i];
            log_info("Stimulus name: %s\n", e.key);
            struct config *tmp = (struct config *)e.value;
            print_purkinje_config_values(tmp);
            log_msg(LOG_LINE_SEPARATOR);
        }
    }

    if(options->domain_config) {
        print_domain_config_values(options->domain_config);
        log_msg(LOG_LINE_SEPARATOR);
    }

    if(options->purkinje_config) {
        print_purkinje_config_values(options->purkinje_config);
        log_msg(LOG_LINE_SEPARATOR);
    }

    if(options->extra_data_config) {
        print_extra_data_config_values(options->extra_data_config);
        log_msg(LOG_LINE_SEPARATOR);
    }

    if(options->purkinje_extra_data_config) {
        print_extra_data_config_values(options->purkinje_extra_data_config);
        log_msg(LOG_LINE_SEPARATOR);
    }

    if(options->update_monodomain_config) {
        print_update_monodomain_config_values(options->update_monodomain_config);
        log_msg(LOG_LINE_SEPARATOR);
    }

    if(options->assembly_matrix_config) {
        print_assembly_matrix_config_values(options->assembly_matrix_config);
        log_msg(LOG_LINE_SEPARATOR);
    }

    if(options->calc_ecg_config) {
        print_calc_ecg_config_values(options->calc_ecg_config);
        log_msg(LOG_LINE_SEPARATOR);
    }
}

void configure_monodomain_solver_from_options(struct monodomain_solver *the_monodomain_solver, struct user_options *options) {

    assert(the_monodomain_solver);
    assert(options);

    the_monodomain_solver->num_threads = options->num_threads;
    the_monodomain_solver->final_time = options->final_time;

    the_monodomain_solver->refine_each = options->refine_each;
    the_monodomain_solver->derefine_each = options->derefine_each;
    the_monodomain_solver->refinement_bound = options->ref_bound;
    the_monodomain_solver->derefinement_bound = options->deref_bound;

    the_monodomain_solver->abort_on_no_activity = options->abort_no_activity;

    the_monodomain_solver->dt = options->dt_pde;

    the_monodomain_solver->beta = options->beta;
    the_monodomain_solver->cm = options->cm;
    the_monodomain_solver->start_adapting_at = options->start_adapting_at;
}

void compute_pmj_current_purkinje_to_tissue(struct ode_solver *the_ode_solver, struct grid *the_grid, struct terminal *the_terminals) {
    assert(the_ode_solver);
    assert(the_grid);
    assert(the_terminals);

    // Tissue solution
    struct cell_node **ac = the_grid->active_cells;
    uint32_t n_active = the_grid->num_active_cells;
    real *sv = the_ode_solver->sv;
    uint32_t nodes = the_ode_solver->model_data.number_of_ode_equations;

    // Purkinje solution
    struct cell_node **ac_purkinje = the_grid->purkinje->purkinje_cells;

    // Purkinje coupling parameters
    real rpmj = the_grid->purkinje->network->rpmj;
    real pmj_scale = the_grid->purkinje->network->pmj_scale;

    real Gpmj = 1.0 / rpmj;

    if(the_ode_solver->gpu) {
#ifdef COMPILE_GPU
        real *vms;
        uint32_t max_number_of_cells = the_ode_solver->original_num_cells;
        size_t mem_size = max_number_of_cells * sizeof(real);

        vms = (real *)malloc(mem_size);

        if(the_grid->adaptive) {
            // check_cuda_error(cudaMemcpy(vms, sv, mem_size, cudaMemcpyDeviceToHost));
            memcpy_device(vms, sv, mem_size, DEVICE_TO_HOST);
        }

        OMP(parallel for)
        for(uint32_t i = 0; i < n_active; i++) {
            vms[ac[i]->sv_position] = (real)ac[i]->v;
        }

        uint32_t num_of_purkinje_terminals = the_grid->purkinje->network->number_of_terminals;
        for(uint32_t i = 0; i < num_of_purkinje_terminals; i++) {

            // Compute the PMJ current
            real Ipmj = 0.0;
            uint32_t num_tissue_cells = arrlen(the_terminals[i].tissue_cells);
            uint32_t purkinje_index = the_terminals[i].purkinje_cell->id;
            rpmj = the_terminals[i].purkinje_cell->rpmj;
            Gpmj = 1.0 / rpmj;
            for(uint32_t j = 0; j < num_tissue_cells; j++) {
                uint32_t tissue_index = the_terminals[i].tissue_cells[j]->sv_position;
                Ipmj += (vms[tissue_index] - ac_purkinje[purkinje_index]->v);
            }
            Ipmj *= (Gpmj / pmj_scale);

            // Add this current to the RHS from each tissue cell
            for(uint32_t j = 0; j < num_tissue_cells; j++) {
                uint32_t tissue_index = the_terminals[i].tissue_cells[j]->sv_position;

                ac[tissue_index]->b -= Ipmj;
            }
        }

        // check_cuda_error(cudaMemcpy(sv, vms, mem_size, cudaMemcpyHostToDevice));
        memcpy_device(sv, vms, mem_size, HOST_TO_DEVICE);
        free(vms);
#endif
    } else {
        uint32_t num_of_purkinje_terminals = the_grid->purkinje->network->number_of_terminals;
        for(uint32_t i = 0; i < num_of_purkinje_terminals; i++) {

            // Compute the PMJ current
            real Ipmj = 0.0;
            uint32_t num_tissue_cells = arrlen(the_terminals[i].tissue_cells);
            uint32_t purkinje_index = the_terminals[i].purkinje_cell->id;
            rpmj = the_terminals[i].purkinje_cell->rpmj;
            Gpmj = 1.0 / rpmj;
            for(uint32_t j = 0; j < num_tissue_cells; j++) {
                uint32_t tissue_index = the_terminals[i].tissue_cells[j]->sv_position;
                Ipmj += (sv[tissue_index * nodes] - ac_purkinje[purkinje_index]->v);
            }
            Ipmj *= (Gpmj / pmj_scale);

            // Add this current to the RHS from each tissue cell
            for(uint32_t j = 0; j < num_tissue_cells; j++) {
                uint32_t tissue_index = the_terminals[i].tissue_cells[j]->sv_position;

                ac[tissue_index]->b -= Ipmj;
            }
        }
    }
}

void compute_pmj_current_tissue_to_purkinje(struct ode_solver *the_purkinje_ode_solver, struct grid *the_grid, struct terminal *the_terminals) {
    assert(the_purkinje_ode_solver);
    assert(the_grid);
    assert(the_terminals);

    // Tissue solution
    struct cell_node **ac = the_grid->active_cells;
    // Purkinje solution
    struct cell_node **ac_purkinje = the_grid->purkinje->purkinje_cells;

    real *sv = the_purkinje_ode_solver->sv;
    uint32_t nodes = the_purkinje_ode_solver->model_data.number_of_ode_equations;

    // Purkinje coupling parameters
    real rpmj = the_grid->purkinje->network->rpmj;
    real pmj_scale = the_grid->purkinje->network->pmj_scale;
    real asymm_ratio = the_grid->purkinje->network->asymm_ratio;

    real Gpmj = 1.0 / rpmj;

    if(the_purkinje_ode_solver->gpu) {
#ifdef COMPILE_GPU

        real *vms;
        uint32_t max_number_of_cells = the_purkinje_ode_solver->original_num_cells;
        size_t mem_size = max_number_of_cells * sizeof(real);

        vms = (real *)malloc(mem_size);

        // check_cuda_error(cudaMemcpy(vms, sv, mem_size, cudaMemcpyDeviceToHost));
        memcpy_device(vms, sv, mem_size, DEVICE_TO_HOST);

        uint32_t num_of_purkinje_terminals = the_grid->purkinje->network->number_of_terminals;
        for(uint32_t i = 0; i < num_of_purkinje_terminals; i++) {

            // Compute the PMJ current
            real Ipmj = 0.0;
            uint32_t num_tissue_cells = arrlen(the_terminals[i].tissue_cells);
            uint32_t purkinje_index = the_terminals[i].purkinje_cell->id;
            rpmj = the_terminals[i].purkinje_cell->rpmj;
            Gpmj = 1.0 / rpmj;
            for(uint32_t j = 0; j < num_tissue_cells; j++) {
                uint32_t tissue_index = the_terminals[i].tissue_cells[j]->sv_position;
                Ipmj += (vms[purkinje_index] - ac[tissue_index]->v);
            }
            // Asymmetry of conduction across the PMJ
            Ipmj *= (Gpmj / (pmj_scale * asymm_ratio));

            // Add this current to the RHS of the Purkinje cell
            ac_purkinje[purkinje_index]->b -= Ipmj;
        }

        free(vms);
#endif
    } else {
        uint32_t num_of_purkinje_terminals = the_grid->purkinje->network->number_of_terminals;
        for(uint32_t i = 0; i < num_of_purkinje_terminals; i++) {

            // Compute the PMJ current
            real Ipmj = 0.0;
            uint32_t num_tissue_cells = arrlen(the_terminals[i].tissue_cells);
            uint32_t purkinje_index = the_terminals[i].purkinje_cell->id;
            rpmj = the_terminals[i].purkinje_cell->rpmj;
            Gpmj = 1.0 / rpmj;
            for(uint32_t j = 0; j < num_tissue_cells; j++) {
                uint32_t tissue_index = the_terminals[i].tissue_cells[j]->sv_position;
                Ipmj += (sv[purkinje_index * nodes] - ac[tissue_index]->v);
            }
            Ipmj *= (Gpmj / (pmj_scale * asymm_ratio));

            // Add this current to the RHS of the Purkinje cell
            ac_purkinje[purkinje_index]->b -= Ipmj;
        }
    }
}

// TODO: Maybe write a library to the PMJ coupling ...
void write_pmj_delay(struct grid *the_grid, struct config *config, struct terminal *the_terminals) {
    assert(the_grid);
    assert(config);
    assert(the_terminals);

    struct save_coupling_with_activation_times_persistent_data *persistent_data =
        (struct save_coupling_with_activation_times_persistent_data *)config->persistent_data;
    char *main_function_name = config->main_function_name;

    if(strcmp(main_function_name, "save_purkinje_coupling_with_activation_times") == 0) {

        char *output_dir = NULL;
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config, "output_dir");

        sds output_dir_with_file = sdsnew(output_dir);
        output_dir_with_file = sdscat(output_dir_with_file, "/pmj_delay.csv");
        log_warn("PMJ delay information will be saved at:> '%s'\n", output_dir_with_file);

        FILE *output_file = NULL;
        output_file = fopen(output_dir_with_file, "w");
        fprintf(output_file, "curPulse,curTerm,pkLAT,meanTissLAT,pmjDelay,isActive,hasBlock\n");

        uint32_t num_terminals = the_grid->purkinje->network->number_of_terminals;

        bool has_block;
        uint32_t purkinje_index;
        struct node *purkinje_cell;
        struct point_3d cell_coordinates;
        real_cpu center_x, center_y, center_z;
        struct cell_node **purkinje_cells = the_grid->purkinje->purkinje_cells;

        // Get the total number of pulses
        purkinje_cell = the_terminals[0].purkinje_cell;
        purkinje_index = purkinje_cell->id;

        center_x = purkinje_cells[purkinje_index]->center.x;
        center_y = purkinje_cells[purkinje_index]->center.y;
        center_z = purkinje_cells[purkinje_index]->center.z;

        cell_coordinates.x = center_x;
        cell_coordinates.y = center_y;
        cell_coordinates.z = center_z;

        int n_pulses = 0;
        n_pulses = (int)hmget(persistent_data->purkinje_num_activations, cell_coordinates);

        // For each pulses calculate its PMJ delay
        for(int k = 0; k < n_pulses; k++) {

            // fprintf(output_file,"====================== PULSE %u ======================\n", k+1);
            for(uint32_t i = 0; i < num_terminals; i++) {

                has_block = false;
                uint32_t term_id = i;
                bool is_terminal_active = the_terminals[i].active;

                // [PURKINJE] Get the informaion from the Purkinje cell
                purkinje_cell = the_terminals[i].purkinje_cell;
                purkinje_index = purkinje_cell->id;

                center_x = purkinje_cells[purkinje_index]->center.x;
                center_y = purkinje_cells[purkinje_index]->center.y;
                center_z = purkinje_cells[purkinje_index]->center.z;

                cell_coordinates.x = center_x;
                cell_coordinates.y = center_y;
                cell_coordinates.z = center_z;

                int n_activations_purkinje = 0;
                float *activation_times_array_purkinje = NULL;

                n_activations_purkinje = (int)hmget(persistent_data->purkinje_num_activations, cell_coordinates);
                activation_times_array_purkinje = (float *)hmget(persistent_data->purkinje_activation_times, cell_coordinates);

                real_cpu purkinje_lat = activation_times_array_purkinje[k];
                // fprintf(output_file,"Terminal %u --> Purkinje cell %u --> LAT = %g\n", i, purkinje_index, purkinje_lat);

                // [TISSUE] Get the information from the Tissue cells
                struct cell_node **tissue_cells = the_terminals[i].tissue_cells;
                uint32_t number_tissue_cells = arrlen(tissue_cells);

                // Calculate the mean LAT of the tissue cells surrounding the Purkinje cell
                real_cpu mean_tissue_lat = 0.0;
                real_cpu min_tissue_lat = __DBL_MAX__;
                uint32_t cur_pulse = k;
                for(uint32_t j = 0; j < number_tissue_cells; j++) {

                    cell_coordinates.x = tissue_cells[j]->center.x;
                    cell_coordinates.y = tissue_cells[j]->center.y;
                    cell_coordinates.z = tissue_cells[j]->center.z;

                    int n_activations_tissue = 0;
                    float *activation_times_array_tissue = NULL;

                    n_activations_tissue = (int)hmget(persistent_data->tissue_num_activations, cell_coordinates);
                    activation_times_array_tissue = (float *)hmget(persistent_data->tissue_activation_times, cell_coordinates);

                    // Check if the number of activations from the tissue and Purkinje cell are equal
                    if(n_activations_purkinje > n_activations_tissue) {
                        // log_error("[purkinje_coupling] ERROR! The number of activations of the tissue and Purkinje cells are different!\n");
                        // log_error("[purkinje_coupling] Probably there was a block on the anterograde direction!\n");
                        // log_error("[purkinje_coupling] Consider only the result from the second pulse! (retrograde direction)!\n");
                        // fprintf(output_file,"ERROR! Probably there was a block on the anterograde direction!\n");
                        has_block = true;
                        cur_pulse = 0;
                        return;
                    }
                    mean_tissue_lat += activation_times_array_tissue[cur_pulse];
                    if(activation_times_array_tissue[cur_pulse] < min_tissue_lat) {
                        min_tissue_lat = activation_times_array_tissue[cur_pulse];
                    }
                }

                if(is_terminal_active) {
                    mean_tissue_lat /= (real_cpu)number_tissue_cells;

                    // PMJ delay is calculated using the mean LAT of the coupled tissue cells minus the LAT of the terminal Purkinje cell
                    real_cpu pmj_delay = (mean_tissue_lat - purkinje_lat);

                    // pulse_id, terminal_id, purkinje_lat, mean_tissue_lat, pmj_delay, is_active, has_block
                    fprintf(output_file, "%d,%u,%g,%g,%g,%d,%d\n", k, term_id, purkinje_lat, mean_tissue_lat, pmj_delay, (int)is_terminal_active,
                            (int)has_block);
                }
            }
        }
        fclose(output_file);
    } else {
        log_error("[purkinje_coupling] ERROR! No 'persistant_data' was found!\n");
        log_error("[purkinje_coupling] You must use the 'save_purkinje_coupling_with_activation_times' function to print the PMJ delay!\n");
    }
}

// TODO: Find a better place to this function ...
void write_terminals_info(struct grid *the_grid, struct config *config, struct terminal *the_terminals) {
    assert(the_grid);
    assert(config);
    assert(the_terminals);

    char *main_function_name = config->main_function_name;

    if(strcmp(main_function_name, "save_purkinje_coupling_with_activation_times") == 0) {

        char *output_dir = NULL;
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config, "output_dir");

        // Purkinje terminal output file
        uint32_t num_coupled_tissue_cells = 0;
        uint32_t num_terminals = the_grid->purkinje->network->number_of_terminals;
        uint32_t purkinje_index;
        struct node *purkinje_cell;
        real_cpu center_x, center_y, center_z;
        struct cell_node **purkinje_cells = the_grid->purkinje->purkinje_cells;

        sds purkinje_term_file = sdsnew(output_dir);
        purkinje_term_file = sdscat(purkinje_term_file, "/purkinje_terminals.vtk");
        log_warn("Purkinje terminal information will be saved at:> '%s'\n", purkinje_term_file);

        FILE *output_pk_term_file = NULL;
        output_pk_term_file = fopen(purkinje_term_file, "w");

        fprintf(output_pk_term_file, "# vtk DataFile Version 4.2\n");
        fprintf(output_pk_term_file, "vtk output\n");
        fprintf(output_pk_term_file, "ASCII\n");
        fprintf(output_pk_term_file, "DATASET POLYDATA\n");
        fprintf(output_pk_term_file, "POINTS %u float\n", num_terminals);

        for(uint32_t i = 0; i < num_terminals; i++) {

            struct cell_node **tissue_cells = the_terminals[i].tissue_cells;
            num_coupled_tissue_cells += arrlen(tissue_cells);

            purkinje_cell = the_terminals[i].purkinje_cell;
            purkinje_index = purkinje_cell->id;

            center_x = purkinje_cells[purkinje_index]->center.x;
            center_y = purkinje_cells[purkinje_index]->center.y;
            center_z = purkinje_cells[purkinje_index]->center.z;

            fprintf(output_pk_term_file, "%g %g %g\n", center_x, center_y, center_z);
        }
        fprintf(output_pk_term_file, "VERTICES %u %u\n", num_terminals, num_terminals * 2);
        for(uint32_t i = 0; i < num_terminals; i++) {
            fprintf(output_pk_term_file, "1 %u\n", i);
        }
        fprintf(output_pk_term_file, "POINT_DATA %u\n", num_terminals);
        fprintf(output_pk_term_file, "FIELD FieldData 1\n");
        fprintf(output_pk_term_file, "isActive 1 %u float\n", num_terminals);
        for(uint32_t i = 0; i < num_terminals; i++) {
            bool is_terminal_active = the_terminals[i].active;
            fprintf(output_pk_term_file, "%d\n", (int)is_terminal_active);
        }
        fclose(output_pk_term_file);

        // Purkinje-tissue coupled file
        sds coupled_tissue_cells_file = sdsnew(output_dir);
        coupled_tissue_cells_file = sdscat(coupled_tissue_cells_file, "/coupled_tissue_cells.vtk");
        log_warn("Coupled tissue information will be saved at:> '%s'\n", coupled_tissue_cells_file);

        FILE *output_coupled_tiss_file = NULL;
        output_coupled_tiss_file = fopen(coupled_tissue_cells_file, "w");

        fprintf(output_coupled_tiss_file, "# vtk DataFile Version 4.2\n");
        fprintf(output_coupled_tiss_file, "vtk output\n");
        fprintf(output_coupled_tiss_file, "ASCII\n");
        fprintf(output_coupled_tiss_file, "DATASET POLYDATA\n");
        fprintf(output_coupled_tiss_file, "POINTS %u float\n", num_coupled_tissue_cells);

        for(uint32_t i = 0; i < num_terminals; i++) {
            struct cell_node **tissue_cells = the_terminals[i].tissue_cells;
            for(uint32_t j = 0; j < arrlen(tissue_cells); j++) {
                center_x = tissue_cells[j]->center.x;
                center_y = tissue_cells[j]->center.y;
                center_z = tissue_cells[j]->center.z;
                fprintf(output_coupled_tiss_file, "%g %g %g\n", center_x, center_y, center_z);
            }
        }
        fprintf(output_coupled_tiss_file, "VERTICES %u %u\n", num_coupled_tissue_cells, num_coupled_tissue_cells * 2);
        for(uint32_t i = 0; i < num_coupled_tissue_cells; i++) {
            fprintf(output_coupled_tiss_file, "1 %u\n", i);
        }
        fprintf(output_coupled_tiss_file, "POINT_DATA %u\n", num_coupled_tissue_cells);
        fprintf(output_coupled_tiss_file, "FIELD FieldData 2\n");
        fprintf(output_coupled_tiss_file, "isActive 1 %u float\n", num_coupled_tissue_cells);
        for(uint32_t i = 0; i < num_terminals; i++) {
            bool is_terminal_active = the_terminals[i].active;
            struct cell_node **tissue_cells = the_terminals[i].tissue_cells;
            for(uint32_t j = 0; j < arrlen(tissue_cells); j++) {
                fprintf(output_coupled_tiss_file, "%d\n", (int)is_terminal_active);
            }
        }
        fprintf(output_coupled_tiss_file, "METADATA\n");
        fprintf(output_coupled_tiss_file, "INFORMATION 0\n\n");
        fprintf(output_coupled_tiss_file, "associatedTerminal 1 %u float\n", num_coupled_tissue_cells);
        for(uint32_t i = 0; i < num_terminals; i++) {
            struct cell_node **tissue_cells = the_terminals[i].tissue_cells;
            for(uint32_t j = 0; j < arrlen(tissue_cells); j++) {
                fprintf(output_coupled_tiss_file, "%u\n", i);
            }
        }
        fprintf(output_coupled_tiss_file, "METADATA\n");
        fprintf(output_coupled_tiss_file, "INFORMATION 0\n\n");
        fclose(output_coupled_tiss_file);
    } else {
        log_error("[purkinje_coupling] ERROR! No 'persistant_data' was found!\n");
        log_error("[purkinje_coupling] You must use the 'save_purkinje_coupling_with_activation_times' function to print the PMJ delay!\n");
    }
}
