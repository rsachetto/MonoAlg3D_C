//
// Created by sachetto on 03/10/17.
//

#include "monodomain_solver.h"
#include "../utils/file_utils.h"
#include "../utils/stop_watch.h"
#include "../libraries_common/common_data_structures.h"

#ifdef COMPILE_CUDA
#include "../gpu_utils/gpu_utils.h"
#endif

#ifdef COMPILE_OPENGL
#include "../draw/draw.h"
#endif

#include "../string/sds.h"
#include <assert.h>
#include <inttypes.h>

#include "../config/assembly_matrix_config.h"
#include "../config/domain_config.h"
#include "../config/purkinje_config.h"
#include "../config/stim_config.h"
#include "../config/linear_system_solver_config.h"

#include "../single_file_libraries/stb_ds.h"
#include "../config_helpers/config_helpers.h"

#include <unistd.h>

#include <stdio.h>
#include <float.h>

struct monodomain_solver *new_monodomain_solver() {

    struct monodomain_solver *result = (struct monodomain_solver *)malloc(sizeof(struct monodomain_solver));

    result->beta = 0.14;
    result->cm = 1.0;
    result->current_time = 0.0;
    result->current_count = 0;

    result->kappa_x = 0.0;
    result->kappa_y = 0.0;
    result->kappa_z = 0.0;

    result->calc_activation_time = false;
    result->print_conductivity = false;
    result->print_min_vm = false;
    result->print_max_vm = false;
    result->print_apd = false;

    return result;
}

int solve_monodomain(struct monodomain_solver *the_monodomain_solver, struct ode_solver *the_ode_solver,
                      struct grid *the_grid, struct user_options *configs)
{

    assert(configs);

    assert(the_grid);
    assert(the_monodomain_solver);
    assert(the_ode_solver);

    print_to_stdout_and_file(LOG_LINE_SEPARATOR);

    long ode_total_time = 0, cg_total_time = 0, total_write_time = 0, total_mat_time = 0, total_ref_time = 0,
         total_deref_time = 0, cg_partial, total_config_time = 0;
    long purkinje_ode_total_time = 0, purkinje_cg_total_time = 0, purkinje_total_write_time = 0, purkinje_total_mat_time = 0,
         purkinje_cg_partial;

    uint32_t total_cg_it = 0, purkinje_total_cg_it = 0;

    struct stop_watch solver_time, ode_time, cg_time, part_solver, part_mat, write_time, ref_time, deref_time,
        config_time;
    struct stop_watch purkinje_solver_time, purkinje_ode_time, purkinje_cg_time, purkinje_part_solver, purkinje_part_mat;

    init_stop_watch(&config_time);

    start_stop_watch(&config_time);

    ///////MAIN CONFIGURATION BEGIN//////////////////
    init_ode_solver_with_cell_model(the_ode_solver);
    struct string_voidp_hash_entry *stimuli_configs = configs->stim_configs;
    struct string_voidp_hash_entry *purkinje_stimuli_configs = configs->purkinje_stim_configs;
    struct config *extra_data_config = configs->extra_data_config;
    struct config *domain_config = configs->domain_config;
    struct config *purkinje_config = configs->purkinje_config;
    struct config *assembly_matrix_config = configs->assembly_matrix_config;
    struct config *linear_system_solver_config = configs->linear_system_solver_config;
    struct config *save_mesh_config = configs->save_mesh_config;
    struct config *save_state_config = configs->save_state_config;
    struct config *restore_state_config = configs->restore_state_config;
    struct config *update_monodomain_config = configs->update_monodomain_config;

    bool has_extra_data = (extra_data_config != NULL);

    real_cpu last_stimulus_time = -1.0;
    bool has_any_periodic_stim = false;

    // Tissue stimuli
    if(stimuli_configs)
    {
        // Init all stimuli
        STIM_CONFIG_HASH_FOR_INIT_FUNCTIONS(stimuli_configs);

        // Find last stimuli
        size_t s_size = shlen(stimuli_configs);
        real_cpu s_end;
        real_cpu stim_start = 0.0;
        real_cpu stim_duration = 0.0;
        real_cpu stim_period = 0;
        bool unnused;

        for(unsigned long i = 0; i < s_size; i++)
        {

            struct config *sconfig = (struct config*) stimuli_configs[i].value;

            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, stim_start, sconfig->config_data, "start");
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, stim_duration, sconfig->config_data, "duration");
            GET_PARAMETER_NUMERIC_VALUE(real_cpu, stim_period, sconfig->config_data, "period", unnused);

            s_end = stim_start + stim_duration;

            has_any_periodic_stim |= (bool)(stim_period > 0.0);

            if(s_end > last_stimulus_time)
            {
                last_stimulus_time = s_end;
            }

        }
    }

    // Purkinje stimuli
    if (purkinje_stimuli_configs)
    {
        // Init all stimuli
        STIM_CONFIG_HASH_FOR_INIT_FUNCTIONS(purkinje_stimuli_configs);

        // Find last stimuli
        size_t s_size = shlen(purkinje_stimuli_configs);
        real_cpu s_end;
        real_cpu stim_start = 0.0;
        real_cpu stim_duration = 0.0;
        real_cpu stim_period = 0;
        bool unnused;

        for(unsigned long i = 0; i < s_size; i++)
        {

            struct config *sconfig = (struct config*) purkinje_stimuli_configs[i].value;

            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, stim_start, sconfig->config_data, "start");
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, stim_duration, sconfig->config_data, "duration");
            GET_PARAMETER_NUMERIC_VALUE(real_cpu, stim_period, sconfig->config_data, "period", unnused);

            s_end = stim_start + stim_duration;

            has_any_periodic_stim |= (bool)(stim_period > 0.0);

            if(s_end > last_stimulus_time)
            {
                last_stimulus_time = s_end;
            }

        }
    }

    // Configure the functions and set the Purkinje mesh domain
    if (purkinje_config)
    {
        init_config_functions(purkinje_config, "shared_libs/libdefault_purkinje.so", "purkinje");
    }

    // Configure the functions and set the mesh domain
    if(domain_config)
    {
        init_config_functions(domain_config, "./shared_libs/libdefault_domains.so", "domain");
    }

    if( !purkinje_config && !domain_config )
    {
        print_to_stderr_and_file_and_exit("Error configuring the domain! No Purkinje or tissue configuration was provided!\n");
    }


    if(assembly_matrix_config)
    {
        init_config_functions(assembly_matrix_config, "./shared_libs/libdefault_matrix_assembly.so", "assembly_matrix");
    }
    else
    {
        print_to_stderr_and_file_and_exit("No assembly matrix configuration provided! Exiting!\n");
    }

    if(linear_system_solver_config)
    {
        init_config_functions(linear_system_solver_config, "./shared_libs/libdefault_linear_system_solver.so", "linear_system_solver");
    }
    else
    {
        print_to_stderr_and_file_and_exit("No linear solver configuration provided! Exiting!\n");
    }

    int print_rate = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, print_rate, save_mesh_config->config_data, "print_rate");

    char *out_dir_name = NULL;
    GET_PARAMETER_VALUE_CHAR_OR_USE_DEFAULT(out_dir_name, save_mesh_config->config_data, "output_dir");

    bool save_to_file = (save_mesh_config != NULL) && (print_rate > 0) && (out_dir_name);

    if(save_to_file)
    {
        init_config_functions(save_mesh_config, "./shared_libs/libdefault_save_mesh.so", "save_result");
    }
    else
    {
        print_to_stdout_and_file("No configuration provided to save the results! The results will not be saved!\n");
    }

    bool save_checkpoint = (save_state_config != NULL);

    if(save_checkpoint && the_grid->adaptive) {
        print_to_stdout_and_file("Saving checkpoint is not implemented for adaptive grids yet!\n");
        save_checkpoint = false;
    }

    if(save_checkpoint)
    {
        init_config_functions(save_state_config, "./shared_libs/libdefault_save_state.so", "save_state");
    }
    else
    {
        print_to_stdout_and_file(
            "No configuration provided to make simulation checkpoints! Chekpoints will not be created!\n");
    }

    bool restore_checkpoint = (restore_state_config != NULL);

    if(restore_checkpoint && the_grid->adaptive) {
        print_to_stdout_and_file("Restoring checkpoint is not implemented for adaptive grids yet!\n");
        restore_checkpoint = false;
    }

    if(restore_state_config)
    {
        init_config_functions(restore_state_config, "./shared_libs/libdefault_restore_state.so", "restore_state");
    }

    if(has_extra_data)
    {
        init_config_functions(extra_data_config, "./shared_libs/libdefault_extra_data.so", "extra_data");
    }

    print_to_stdout_and_file(LOG_LINE_SEPARATOR);

    bool restore_success = false;

    if(restore_checkpoint)
    {
        // Here we only restore the monodomain_solver_state...
        restore_success = ((restore_state_fn *)restore_state_config->main_function)(out_dir_name, restore_state_config, NULL,
                                            the_monodomain_solver, NULL);
    }

    if(update_monodomain_config)
    {
        init_config_functions(update_monodomain_config, "./shared_libs/libdefault_update_monodomain.so", "update_monodomain");
    }
    else
    {
        print_to_stderr_and_file_and_exit("No update monodomain configuration provided! Exiting!\n");
    }

    ///////MAIN CONFIGURATION END//////////////////
    int refine_each = the_monodomain_solver->refine_each;
    int derefine_each = the_monodomain_solver->derefine_each;

    bool redo_matrix;

    bool activity;

    #ifdef COMPILE_CUDA
    bool gpu = the_ode_solver->gpu;
    #endif

    int count = the_monodomain_solver->current_count;

    real_cpu refinement_bound = the_monodomain_solver->refinement_bound;
    real_cpu derefinement_bound = the_monodomain_solver->derefinement_bound;

    bool adaptive = the_grid->adaptive;
    real_cpu start_adpt_at = the_monodomain_solver->start_adapting_at;
    real_cpu dt_pde = the_monodomain_solver->dt;
    real_cpu finalT = the_monodomain_solver->final_time;
    real_cpu dt_ode = the_ode_solver->min_dt;

#ifdef COMPILE_OPENGL
    bool draw = configs->draw;
    if (draw) {
        draw_config.grid_info.grid_to_draw = the_grid;
        draw_config.simulating = true;
        draw_config.paused = !configs->start_visualization_unpaused;
    } else {
        draw_config.paused = false;
    }
#endif

#ifdef COMPILE_CUDA
    if(gpu) {
        int device_count;
        int device = the_ode_solver->gpu_id;
        check_cuda_errors(cudaGetDeviceCount(&device_count));
        struct cudaDeviceProp prop;
        check_cuda_errors(cudaGetDeviceProperties(&prop, the_ode_solver->gpu_id));
        print_to_stdout_and_file("%d devices available, running on Device %d: %s\n", device_count, device, prop.name);
        check_cuda_errors(cudaSetDevice(device));
    }
#endif

    int success;

    struct ode_solver *the_purkinje_ode_solver = NULL;
    if (purkinje_config)
    {
        // Allocate a new 'ode_solver' for the Purkinje
        the_purkinje_ode_solver = new_ode_solver();

        // Here we configure the Purkinje ode_solver using the [purkinje_ode_solver] parameters
        // If there is no [purkinje_ode_solver] section we configure the Purkinje ODE solver using the input from the [ode_solver] section
        if (!domain_config)     // ONLY Purkinje simulation
            configure_purkinje_ode_solver_from_ode_solver(the_purkinje_ode_solver,the_ode_solver);
        // Otherwise, there is a [purkinje_ode_solver] section and we are doing a coupled simulation
        else                    // Purkinje + Tissue simulation
            configure_purkinje_ode_solver_from_options(the_purkinje_ode_solver,configs);

        success = ((set_spatial_purkinje_fn*) purkinje_config->main_function)(purkinje_config,the_grid,the_purkinje_ode_solver);
        if(!success)
        {
            print_to_stderr_and_file_and_exit("Error configuring the Purkinje domain!\n");
        }

    }

    if (domain_config)
    {
        success = ((set_spatial_domain_fn *) domain_config->main_function)(domain_config, the_grid);

        if (!success)
        {
            print_to_stderr_and_file_and_exit("Error configuring the tissue domain!\n");
        }
    }

    if (!purkinje_config && !domain_config)
    {
        print_to_stderr_and_file_and_exit("Error configuring the domain! No Purkinje or tissue configuration was provided!\n");
    }

    if(restore_checkpoint)
    {
        // TODO: Create a Purkinje restore function in the 'restore_library' and put here ...
        restore_success &= ((restore_state_fn*)restore_state_config->main_function)(out_dir_name, restore_state_config, the_grid, NULL, NULL);
    }

    real_cpu start_dx, start_dy, start_dz;
    real_cpu max_dx, max_dy, max_dz;

    start_dx = start_dy = start_dz = 100.0;
    max_dx = max_dy = max_dz = 100.0;

/*
    // TODO: Eliminate this ...
    if (purkinje_config)
    {
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, start_dx, purkinje_config->config_data, "start_discretization");
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, start_dy, purkinje_config->config_data, "start_discretization");
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, start_dz, purkinje_config->config_data, "start_discretization");

        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, max_dx, purkinje_config->config_data, "start_discretization");
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, max_dy, purkinje_config->config_data, "start_discretization");
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, max_dz, purkinje_config->config_data, "start_discretization");
    }
*/

    if (domain_config)
    {
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, start_dx, domain_config->config_data, "start_dx");
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, start_dy, domain_config->config_data, "start_dy");
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, start_dz, domain_config->config_data, "start_dz");

        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, max_dx, domain_config->config_data, "maximum_dx");
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, max_dy, domain_config->config_data, "maximum_dy");
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, max_dz, domain_config->config_data, "maximum_dz");
    }

    uint32_t original_num_cells, original_num_purkinje_cells;
    if (domain_config)
    {
        order_grid_cells(the_grid);
        original_num_cells = the_grid->num_active_cells;
        the_ode_solver->original_num_cells = original_num_cells;
        the_ode_solver->num_cells_to_solve = original_num_cells;
    }

    // Purkinje section
    if (purkinje_config)
    {
        original_num_purkinje_cells = the_grid->the_purkinje->number_of_purkinje_cells;
        the_purkinje_ode_solver->original_num_cells = original_num_purkinje_cells;
        the_purkinje_ode_solver->num_cells_to_solve = original_num_purkinje_cells;
    }

    save_old_cell_positions(the_grid);

    if(adaptive)
    {
        update_cells_to_solve(the_grid, the_ode_solver);
    }

    // NEW FUNCTION !!!
    // Map the indexes from the closest endocardium cells that are next to the Purkinje terminals
    // =============================================
    struct terminal *the_terminals = NULL;
    if (domain_config && purkinje_config)
        the_terminals = link_purkinje_to_endocardium(the_grid);
    // =============================================

    print_to_stdout_and_file("Setting ODE's initial conditions\n");

    // TODO: Include the Purkinje extra data here ...
    if (has_extra_data)
    {
      if (purkinje_config)
        set_ode_extra_data(extra_data_config,the_grid,the_purkinje_ode_solver);
      if (domain_config)
        set_ode_extra_data(extra_data_config, the_grid, the_ode_solver);
    }

    if (domain_config)
        set_ode_initial_conditions_for_all_volumes(the_ode_solver, configs->ode_extra_config);
    if (purkinje_config)
        set_ode_initial_conditions_for_all_volumes(the_purkinje_ode_solver, configs->ode_extra_config);

    // We need to call this function after because of the pitch.... maybe we have to change the way
    // we pass this parameters to the cell model....
    if(restore_checkpoint)
    {
        restore_success &= ((restore_state_fn*)restore_state_config->main_function)(out_dir_name, restore_state_config, NULL, NULL, the_ode_solver);
    }

    real_cpu initial_v, purkinje_initial_v;
    initial_v = the_ode_solver->model_data.initial_v;
    if (purkinje_config)
        purkinje_initial_v = the_purkinje_ode_solver->model_data.initial_v;

    total_config_time = stop_stop_watch(&config_time);

    print_solver_info(the_monodomain_solver, the_ode_solver, the_purkinje_ode_solver, the_grid, configs);

    int ode_step = 1;

    if(dt_pde >= dt_ode)
    {
        ode_step = (int)(dt_pde / dt_ode);
        print_to_stdout_and_file("Solving EDO %d times before solving PDE\n", ode_step);
    }
    else
    {
        print_to_stdout_and_file("WARNING: EDO time step is greater than PDE time step. Adjusting to EDO time "
                                 "step: %lf\n",
                                 dt_ode);
        dt_pde = dt_ode;
    }

    fflush(stdout);

    init_stop_watch(&solver_time);
    init_stop_watch(&ode_time);
    init_stop_watch(&cg_time);
    init_stop_watch(&part_solver);
    init_stop_watch(&part_mat);
    init_stop_watch(&write_time);
    init_stop_watch(&ref_time);
    init_stop_watch(&deref_time);
    if (purkinje_config)
    {
        init_stop_watch(&purkinje_ode_time);
        init_stop_watch(&purkinje_cg_time);
        init_stop_watch(&purkinje_part_solver);
    }

    start_stop_watch(&part_mat);

    if(!restore_checkpoint || !restore_success)
    {
        ((set_pde_initial_condition_fn*)assembly_matrix_config->init_function)(assembly_matrix_config, the_monodomain_solver, the_grid, initial_v,purkinje_initial_v);
    }

    ((assembly_matrix_fn*) assembly_matrix_config->main_function)(assembly_matrix_config, the_monodomain_solver, the_grid);

    // TESTING LU DECOMPOSITION
/*
    real_cpu **lu = NULL;
    uint32_t *pivot = NULL;
    if (purkinje_config)
    {
        uint32_t num_rows = the_grid->the_purkinje->num_active_purkinje_cells;
        uint32_t num_cols = the_grid->the_purkinje->num_active_purkinje_cells;

        pivot = (uint32_t*)calloc(num_rows,sizeof(uint32_t));

        lu = allocate_matrix_LU(num_rows,num_cols);
        lu_decomposition(lu,pivot,the_grid->the_purkinje->purkinje_cells,the_grid->the_purkinje->num_active_purkinje_cells);
    }
 */
    // TESTING LU DECOMPOSITION



    total_mat_time = stop_stop_watch(&part_mat);
    start_stop_watch(&solver_time);

    int save_state_rate = 0;

    if(save_checkpoint)
    {
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, save_state_rate, save_state_config->config_data, "save_rate");
    }

    real_cpu vm_threshold = configs->vm_threshold;

    bool abort_on_no_activity = the_monodomain_solver->abort_on_no_activity;
    bool calc_activation_time = the_monodomain_solver->calc_activation_time;
    bool print_conductivity = the_monodomain_solver->print_conductivity;
    bool print_min_vm = the_monodomain_solver->print_min_vm;
    bool print_max_vm = the_monodomain_solver->print_max_vm;
    bool print_apd = the_monodomain_solver->print_apd;
    bool calc_retropropagation = the_grid->the_purkinje->the_network->calc_retropropagation;

    real_cpu solver_error, purkinje_solver_error;
    uint32_t solver_iterations = 0, purkinje_solver_iterations = 0;

    real *spatial_stim_currents = NULL;
    real *purkinje_spatial_stim_currents = NULL;

    if(stimuli_configs)
    {
        spatial_stim_currents = (real*)malloc(sizeof(real)*original_num_cells);
        set_spatial_stim(stimuli_configs, the_grid);
    }

    if (purkinje_stimuli_configs)
    {
        purkinje_spatial_stim_currents = (real*)malloc(sizeof(real)*original_num_purkinje_cells);
        set_spatial_purkinje_stim(purkinje_stimuli_configs, the_grid);
    }

    real_cpu cur_time = the_monodomain_solver->current_time;

    if(save_mesh_config != NULL)
    {
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, print_rate, save_mesh_config->config_data, "print_rate");
    }

    print_to_stdout_and_file("Starting simulation\n");

    struct stop_watch iteration_time_watch;
    long iteration_time;

    init_stop_watch(&iteration_time_watch);

    CALL_INIT_LINEAR_SYSTEM(linear_system_solver_config, the_grid);
    CALL_INIT_SAVE_MESH(save_mesh_config);


#ifdef COMPILE_OPENGL
    if(configs->draw) {
        translate_mesh_to_origin(the_grid);
    }
    draw_config.grid_info.loaded = true;
#endif

    // Main simulation loop start
    while(cur_time <= finalT)
    {
        start_stop_watch(&iteration_time_watch);

        #ifdef COMPILE_OPENGL
        if(draw) {
            omp_set_lock(&draw_config.sleep_lock);
            if (draw_config.restart) {
                draw_config.time = 0.0;

                CALL_END_LINEAR_SYSTEM(linear_system_solver_config);
                CALL_END_SAVE_MESH(save_mesh_config);
                return RESTART_SIMULATION;
            }
            if (draw_config.exit)  {
                CALL_END_LINEAR_SYSTEM(linear_system_solver_config);
                CALL_END_SAVE_MESH(save_mesh_config);
                return END_SIMULATION;
            }
        }
        #endif

        if (save_to_file && (count % print_rate == 0))
        {
            start_stop_watch(&write_time);
            ((save_mesh_fn *)save_mesh_config->main_function)(save_mesh_config, the_grid, count, cur_time, finalT, dt_pde,'v');
            total_write_time += stop_stop_watch(&write_time);
        }

        if (cur_time > 0.0)
        {
            // UPDATE: Tissue and Purkinje (if it exists)
            activity = update_ode_state_vector_and_check_for_activity(vm_threshold, the_ode_solver, the_purkinje_ode_solver, the_grid);

            // MAPPING: Tissue -> Purkinje
            //if (domain_config && purkinje_config)
            //    map_tissue_solution_to_purkinje(the_purkinje_ode_solver,the_grid,the_terminals);


            if (abort_on_no_activity && cur_time > last_stimulus_time)
            {
                if (!activity)
                {
                    print_to_stdout_and_file("No activity, aborting simulation\n");
                    break;
                }

            }
        }

        if (purkinje_config)
        {
            start_stop_watch(&purkinje_ode_time);

            // REACTION: Purkinje
            solve_purkinje_volumes_odes(the_purkinje_ode_solver, the_grid->the_purkinje->number_of_purkinje_cells, cur_time, ode_step, purkinje_stimuli_configs, configs->ode_extra_config);

            purkinje_ode_total_time += stop_stop_watch(&purkinje_ode_time);

            // NEW FEATURE! Calculate the Purkinje activation time here !
            if (cur_time > 0.0 && calc_activation_time)
            {
                calculate_activation_time(cur_time,dt_pde,the_grid->the_purkinje->num_active_purkinje_cells,the_grid->the_purkinje->purkinje_cells,the_purkinje_ode_solver);
            }

            start_stop_watch(&purkinje_ode_time);

    // Implicito
            // UPDATE: Purkinje
            ((update_monodomain_fn*)update_monodomain_config->main_function)(update_monodomain_config, original_num_purkinje_cells, the_monodomain_solver, the_grid->the_purkinje->num_active_purkinje_cells, the_grid->the_purkinje->purkinje_cells, the_purkinje_ode_solver);

            purkinje_ode_total_time += stop_stop_watch(&purkinje_ode_time);

            start_stop_watch(&purkinje_cg_time);

            #ifdef COMPILE_OPENGL
            if (draw) {
                omp_set_lock(&draw_config.draw_lock);
            }
            #endif

            // COUPLING: Calculate the PMJ current from the Tissue to the Purkinje
            if (domain_config)
                compute_pmj_current_tissue_to_purkinje(the_purkinje_ode_solver,the_grid,the_terminals);

            // DIFUSION: Purkinje
            linear_system_solver_purkinje(linear_system_solver_config, the_grid, &purkinje_solver_iterations, &purkinje_solver_error);
            //linear_system_solver_purkinje_lu(lu,pivot,the_grid);

            purkinje_cg_partial = stop_stop_watch(&purkinje_cg_time);

            purkinje_cg_total_time += purkinje_cg_partial;

            purkinje_total_cg_it += purkinje_solver_iterations;

            // MAPPING: Purkinje -> Tissue
            //if (domain_config)
            //    map_purkinje_solution_to_tissue(the_ode_solver,the_grid,the_terminals);

        }

        if (domain_config)
        {

            start_stop_watch(&ode_time);

            // REACTION: Tissue
            solve_all_volumes_odes(the_ode_solver, the_grid->num_active_cells, cur_time, ode_step, stimuli_configs, configs->ode_extra_config);

            ode_total_time += stop_stop_watch(&ode_time);

            // NEW FEATURE! Calculate the tissue activation time here !
            if (cur_time > 0.0 && calc_activation_time)
            {
                calculate_activation_time(cur_time,dt_pde,the_grid->num_active_cells,the_grid->active_cells,the_ode_solver);
            }

            start_stop_watch(&ode_time);

            // UPDATE: Tissue
            ((update_monodomain_fn*)update_monodomain_config->main_function)(update_monodomain_config, original_num_cells, the_monodomain_solver, the_grid->num_active_cells, the_grid->active_cells, the_ode_solver);

            ode_total_time += stop_stop_watch(&ode_time);

            start_stop_watch(&cg_time);

            #ifdef COMPILE_OPENGL
            if (draw) {
                omp_set_lock(&draw_config.draw_lock);
            }
            #endif

            // COUPLING: Calculate the PMJ current from the Purkinje to the Tissue
            if (purkinje_config)
                compute_pmj_current_purkinje_to_tissue(the_ode_solver,the_grid,the_terminals);

            // DIFUSION: Tissue
            ((linear_system_solver_fn *)linear_system_solver_config->main_function)(linear_system_solver_config, the_grid, &solver_iterations, &solver_error);

            cg_partial = stop_stop_watch(&cg_time);

            cg_total_time += cg_partial;

            total_cg_it += solver_iterations;

        }


        if (count % print_rate == 0)
        {
            if (purkinje_config && domain_config)
                print_to_stdout_and_file("t = %.5lf, Iterations = "
                                     "%" PRIu32 ", Error Norm = %e, Number of Tissue Cells:"
                                     "%" PRIu32 ", Tissue CG Iterations time: %ld us\n"
                                         "            , Iterations = "
                                     "%" PRIu32 ", Error Norm = %e, Number of Purkinje Cells:"
                                     "%" PRIu32 ", Purkinje CG Iterations time: %ld us",
                                     cur_time, solver_iterations, solver_error, the_grid->num_active_cells, cg_partial,
                                     purkinje_solver_iterations, purkinje_solver_error, the_grid->the_purkinje->num_active_purkinje_cells,purkinje_cg_partial);
            else if (domain_config)
                print_to_stdout_and_file("t = %lf, Iterations = "
                                     "%" PRIu32 ", Error Norm = %e, Number of Cells:"
                                     "%" PRIu32 ", CG Iterations time: %ld us",
                                     cur_time, solver_iterations, solver_error, the_grid->num_active_cells,
                                     cg_partial);
            else
                print_to_stdout_and_file("t = %lf, Iterations = "
                                     "%" PRIu32 ", Error Norm = %e, Number of Purkinje Cells:"
                                     "%" PRIu32 ", Purkinje CG Iterations time: %ld us",
                                     cur_time, purkinje_solver_iterations, purkinje_solver_error, the_grid->the_purkinje->num_active_purkinje_cells,
                                     purkinje_cg_partial);
        }

        if (adaptive)
        {
            redo_matrix = false;
            if (cur_time >= start_adpt_at)
            {
                if (count % refine_each == 0)
                {

                    start_stop_watch(&ref_time);
                    redo_matrix = refine_grid_with_bound(the_grid, refinement_bound, start_dx, start_dy, start_dz);
                    total_ref_time += stop_stop_watch(&ref_time);
                }

                if (count % derefine_each == 0)
                {
                    start_stop_watch(&deref_time);
                    redo_matrix |= derefine_grid_with_bound(the_grid, derefinement_bound, max_dx, max_dy, max_dz);
                    total_deref_time += stop_stop_watch(&deref_time);
                }

            }
            if (redo_matrix)
            {
                order_grid_cells(the_grid);

                if (stimuli_configs) {
                    if (cur_time <= last_stimulus_time || has_any_periodic_stim) {
                        free(spatial_stim_currents);
                        spatial_stim_currents = (real*)malloc(sizeof(real)*the_grid->num_active_cells);
                        set_spatial_stim(stimuli_configs, the_grid);
                    }
                }
                if (has_extra_data)
                {
                    set_ode_extra_data(extra_data_config, the_grid, the_ode_solver);
                }

                update_cells_to_solve(the_grid, the_ode_solver);

                if (arrlen(the_grid->refined_this_step) > 0)
                {
                    update_state_vectors_after_refinement(the_ode_solver, the_grid->refined_this_step);
                }

                start_stop_watch(&part_mat);
                ((assembly_matrix_fn *)assembly_matrix_config->main_function)(assembly_matrix_config, the_monodomain_solver, the_grid);

                // MAPPING: Update the mapping between the Purkinje mesh and the refined/derefined grid
                if (purkinje_config && domain_config)
                    update_link_purkinje_to_endocardium(the_grid,the_terminals);

                total_mat_time += stop_stop_watch(&part_mat);

            }

        }

        #ifdef COMPILE_OPENGL
        if (configs->draw) {
            omp_unset_lock(&draw_config.draw_lock);
            draw_config.time = cur_time;
        }
        #endif
        count++;
        cur_time += dt_pde;

        if (save_checkpoint)
        {
            if (count != 0 && (count % save_state_rate == 0))
            {
                the_monodomain_solver->current_count = count;
                the_monodomain_solver->current_time = cur_time;
                printf("Saving state with time = %lf, and count = %d\n", the_monodomain_solver->current_time,
                       the_monodomain_solver->current_count);
                ((save_state_fn *)save_state_config->main_function)(out_dir_name, save_state_config, the_grid, the_monodomain_solver, the_ode_solver);
            }
        }


        iteration_time = stop_stop_watch(&iteration_time_watch);

        if ( (count - 1) % print_rate == 0)
        {
            print_to_stdout_and_file(", Total Iteration time: %ld us\n", iteration_time);
        }

    }

    // ------------------------------------------------------------
    // NEW FEATURE ! Save the activation map in a VTK-format file
    if (calc_activation_time)
    {
        print_to_stdout_and_file("Saving activation map!\n");
        ((save_mesh_fn *)save_mesh_config->main_function)(save_mesh_config, the_grid, count, cur_time, finalT, dt_pde,'a');
    }
    // NEW FEATURE ! Save the conductivity map in a VTK-format file
    if (print_conductivity)
    {
        print_to_stdout_and_file("Saving conductivity map!\n");
        ((save_mesh_fn *)save_mesh_config->main_function)(save_mesh_config, the_grid, count, cur_time, finalT, dt_pde,'c');
    }
    if (print_min_vm)
    {
        print_to_stdout_and_file("Saving Minimum Vm map!\n");
        ((save_mesh_fn *)save_mesh_config->main_function)(save_mesh_config, the_grid, count, cur_time, finalT, dt_pde,'m');
    }
    if (print_max_vm)
    {
        print_to_stdout_and_file("Saving Maximum Vm map!\n");
        ((save_mesh_fn *)save_mesh_config->main_function)(save_mesh_config, the_grid, count, cur_time, finalT, dt_pde,'M');
    }
    if (print_apd)
    {
        print_to_stdout_and_file("Saving APD map!\n");
        ((save_mesh_fn *)save_mesh_config->main_function)(save_mesh_config, the_grid, count, cur_time, finalT, dt_pde,'d');
    }
    // ------------------------------------------------------------

    long res_time = stop_stop_watch(&solver_time);
    print_to_stdout_and_file("Resolution Time: %ld μs\n", res_time);
    print_to_stdout_and_file("Assembly matrix time: %ld μs\n", total_mat_time);
    print_to_stdout_and_file("Write time: %ld μs\n", total_write_time);
    print_to_stdout_and_file("Initial configuration time: %ld μs\n", total_config_time);
    if (domain_config)
    {
        print_to_stdout_and_file("ODE Total Time: %ld μs\n", ode_total_time);
        print_to_stdout_and_file("CG Total Time: %ld μs\n", cg_total_time);
        print_to_stdout_and_file("CG Total Iterations: %u\n", total_cg_it);
        print_to_stdout_and_file("Refine time: %ld μs\n", total_ref_time);
        print_to_stdout_and_file("Derefine time: %ld μs\n", total_deref_time);
    }
    if (purkinje_config)
    {
        print_to_stdout_and_file("Purkinje ODE Total Time: %ld μs\n", purkinje_ode_total_time);
        print_to_stdout_and_file("Purkinje CG Total Time: %ld μs\n", purkinje_cg_total_time);
        print_to_stdout_and_file("Purkinje CG Total Iterations: %u\n", purkinje_total_cg_it);
    }

#ifdef COMPILE_OPENGL
   draw_config.solver_time = res_time;
   draw_config.ode_total_time = ode_total_time;
   draw_config.cg_total_time = cg_total_time;
   draw_config.total_mat_time = total_mat_time;
   draw_config.total_ref_time = total_ref_time;
   draw_config.total_deref_time = total_deref_time;
   draw_config.total_write_time = total_write_time;
   draw_config.total_config_time = total_config_time;
   draw_config.total_cg_it  = total_cg_it;
   draw_config.simulating = false;
#endif

    CALL_END_LINEAR_SYSTEM(linear_system_solver_config);
    CALL_END_SAVE_MESH(save_mesh_config);

    return SIMULATION_FINISHED;

}

void set_spatial_stim(struct string_voidp_hash_entry *stim_configs, struct grid *the_grid)
{

    struct config *tmp = NULL;
    size_t n = shlen(stim_configs);

    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;
    bool adaptive = the_grid->adaptive;

    for(size_t i = 0; i < n; i++)
    {
        tmp = (struct config *)stim_configs[i].value;
        ((set_spatial_stim_fn*)tmp->main_function)(tmp,n_active,ac,adaptive);
    }
}

void set_spatial_purkinje_stim(struct string_voidp_hash_entry *stim_configs, struct grid *the_grid)
{
    assert(the_grid->the_purkinje);

    struct config *tmp = NULL;
    size_t n = shlen(stim_configs);

    uint32_t n_active = the_grid->the_purkinje->num_active_purkinje_cells;
    struct cell_node **ac = the_grid->the_purkinje->purkinje_cells;
    bool adaptive = false;

    for(size_t i = 0; i < n; i++)
    {
        tmp = (struct config *)stim_configs[i].value;
        ((set_spatial_stim_fn*)tmp->main_function)(tmp,n_active,ac,adaptive);
    }
}

void set_ode_extra_data(struct config *config, struct grid *the_grid, struct ode_solver *the_ode_solver) {

    free(the_ode_solver->ode_extra_data);
    the_ode_solver->ode_extra_data =
            ((set_extra_data_fn*)config->main_function)(the_grid, config, &(the_ode_solver->extra_data_size));
}

bool update_ode_state_vector_and_check_for_activity(real_cpu vm_threshold, struct ode_solver *the_ode_solver, struct ode_solver *the_purkinje_ode_solver, struct grid *the_grid)
{
    bool act = false;

    // Tissue section
    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    if (the_ode_solver)
    {
        int n_odes = the_ode_solver->model_data.number_of_ode_equations;

        real *sv = the_ode_solver->sv;

        if(the_ode_solver->gpu)
        {
        #ifdef COMPILE_CUDA
            uint32_t max_number_of_cells = the_ode_solver->original_num_cells;
            real *vms;
            size_t mem_size = max_number_of_cells * sizeof(real);

            vms = (real *)malloc(mem_size);

            if(the_grid->adaptive)
                check_cuda_errors(cudaMemcpy(vms, sv, mem_size, cudaMemcpyDeviceToHost));

            #pragma omp parallel for
            for(uint32_t i = 0; i < n_active; i++)
            {
                vms[ac[i]->sv_position] = (real)ac[i]->v;

                if(ac[i]->v > vm_threshold)
                {
                    act = true;
                }
            }

            check_cuda_errors(cudaMemcpy(sv, vms, mem_size, cudaMemcpyHostToDevice));
            free(vms);
        #endif
        }
        else
        {
            #pragma omp parallel for
            for(uint32_t i = 0; i < n_active; i++)
            {
                sv[ac[i]->sv_position * n_odes] = (real)ac[i]->v;

                if(ac[i]->v > vm_threshold)
                {
                    act = true;
                }
            }
        }
    }

    if (the_purkinje_ode_solver)
    {
        // Purkinje section
        uint32_t n_active_purkinje = the_grid->the_purkinje->number_of_purkinje_cells;
        struct cell_node **ac_purkinje = the_grid->the_purkinje->purkinje_cells;

        int n_odes_purkinje = the_purkinje_ode_solver->model_data.number_of_ode_equations;

        real *sv_purkinje = the_purkinje_ode_solver->sv;

        if(the_purkinje_ode_solver->gpu)
        {
        #ifdef COMPILE_CUDA
            uint32_t max_number_of_purkinje_cells = the_purkinje_ode_solver->original_num_cells;
            real *vms_purkinje;
            size_t mem_size_purkinje = max_number_of_purkinje_cells * sizeof(real);

            vms_purkinje = (real *)malloc(mem_size_purkinje);

            if(the_grid->adaptive)
                check_cuda_errors(cudaMemcpy(vms_purkinje, sv_purkinje, mem_size_purkinje, cudaMemcpyDeviceToHost));

            #pragma omp parallel for
            for(uint32_t i = 0; i < n_active_purkinje; i++)
            {
                vms_purkinje[ac_purkinje[i]->sv_position] = (real)ac_purkinje[i]->v;

                if(ac_purkinje[i]->v > vm_threshold)
                {
                    act = true;
                }
            }

            check_cuda_errors(cudaMemcpy(sv_purkinje, vms_purkinje, mem_size_purkinje, cudaMemcpyHostToDevice));
            free(vms_purkinje);
        #endif
        }
        else
        {
            #pragma omp parallel for
            for(uint32_t i = 0; i < n_active_purkinje; i++)
            {
                sv_purkinje[ac_purkinje[i]->sv_position * n_odes_purkinje] = (real)ac_purkinje[i]->v;

                if(ac_purkinje[i]->v > vm_threshold)
                {
                    act = true;
                }
            }
        }
    }

    return act;
}

void save_old_cell_positions(struct grid *the_grid)
{

    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    int i;

    #pragma omp parallel for
    for(i = 0; i < n_active; i++)
    {
        ac[i]->sv_position = ac[i]->grid_position;
    }

    // Purkinje section
    struct grid_purkinje *the_purkinje = the_grid->the_purkinje;

    if (the_purkinje->first_cell)
    {
        uint32_t n_purkinje_active = the_purkinje->num_active_purkinje_cells;
        struct cell_node **ac_purkinje = the_purkinje->purkinje_cells;
        #pragma omp parallel for
        for(i = 0; i < n_purkinje_active; i++)
        {
            //print_to_stdout_and_file("Cell %u -- grid_position = %u\n",i,ac_purkinje[i]->grid_position);
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

    #pragma omp parallel for
    for(i = 0; i < n_active; i++) {
        cts[i] = ac[i]->sv_position;
    }
}


void print_solver_info(struct monodomain_solver *the_monodomain_solver,
                        struct ode_solver *the_ode_solver, struct ode_solver *the_purkinje_ode_solver,
                        struct grid *the_grid, struct user_options *options)
{

    print_to_stdout_and_file(LOG_LINE_SEPARATOR);

    print_to_stdout_and_file("System parameters: \n");
    #if defined(_OPENMP)
    print_to_stdout_and_file("[main] Using OpenMP with %d threads\n", omp_get_max_threads());
    #endif

    print_to_stdout_and_file("[monodomain_solver] Beta = %.10lf, Cm = %.10lf\n", the_monodomain_solver->beta, the_monodomain_solver->cm);
    print_to_stdout_and_file("[monodomain_solver] PDE time step = %lf\n", the_monodomain_solver->dt);
    print_to_stdout_and_file("[monodomain_solver] ODE min time step = %lf\n", the_ode_solver->min_dt);
    print_to_stdout_and_file("[monodomain_solver] Simulation Final Time = %lf\n", the_monodomain_solver->final_time);

    print_to_stdout_and_file(LOG_LINE_SEPARATOR);

    if(the_ode_solver->gpu)
        print_to_stdout_and_file("[ode_solver] Using GPU to solve ODEs\n");
    print_to_stdout_and_file("[ode_solver] Using %s as model lib\n", the_ode_solver->model_data.model_library_path);
    print_to_stdout_and_file("[ode_solver] Initial V: %lf\n", the_ode_solver->model_data.initial_v);
    print_to_stdout_and_file("[ode_solver] Number of ODEs in cell model: %d\n", the_ode_solver->model_data.number_of_ode_equations);

    print_to_stdout_and_file(LOG_LINE_SEPARATOR);

    if(options->ode_extra_config)
    {

        if (shlen(options->ode_extra_config) == 1)
        {
            print_to_stdout_and_file("[ode_solver] Extra ODE Solver parameter:\n");
        }
        else if (shlen(options->ode_extra_config) > 1)
        {
            print_to_stdout_and_file("[ode_solver] Extra ODE Solver parameters:\n");
        }

        STRING_HASH_PRINT_KEY_VALUE_LOG(options->ode_extra_config);
    }

    print_to_stdout_and_file(LOG_LINE_SEPARATOR);

    if (the_purkinje_ode_solver)
    {
        if(the_purkinje_ode_solver->gpu)
            print_to_stdout_and_file("[purkinje_ode_solver] Using GPU to solve ODEs\n");
        print_to_stdout_and_file("[purkinje_ode_solver] Using %s as model lib\n", the_purkinje_ode_solver->model_data.model_library_path);
        print_to_stdout_and_file("[purkinje_ode_solver] Initial V: %lf\n", the_purkinje_ode_solver->model_data.initial_v);
        print_to_stdout_and_file("[purkinje_ode_solver] Number of ODEs in cell model: %d\n", the_purkinje_ode_solver->model_data.number_of_ode_equations);
    }


    print_to_stdout_and_file(LOG_LINE_SEPARATOR);

    if(options->purkinje_ode_extra_config)
    {

        if (shlen(options->purkinje_ode_extra_config) == 1)
        {
            print_to_stdout_and_file("[purkinje_ode_solver] Extra ODE Solver parameter:\n");
        }
        else if (shlen(options->purkinje_ode_extra_config) > 1)
        {
            print_to_stdout_and_file("[purkinje_ode_solver] Extra ODE Solver parameters:\n");
        }

        STRING_HASH_PRINT_KEY_VALUE_LOG(options->purkinje_ode_extra_config);
    }

    print_to_stdout_and_file(LOG_LINE_SEPARATOR);

    print_to_stdout_and_file("[grid] Initial N. of Elements = "
                             "%" PRIu32 "\n",
                             the_grid->num_active_cells);

    if(the_grid->adaptive)
    {
        print_to_stdout_and_file("Using adaptativity\n");
        print_to_stdout_and_file("[monodomain_solver] Refinement Bound = %lf\n", the_monodomain_solver->refinement_bound);
        print_to_stdout_and_file("[monodomain_solver] Derefinement Bound = %lf\n", the_monodomain_solver->derefinement_bound);
        print_to_stdout_and_file("[monodomain_solver] Refining each %d time steps\n", the_monodomain_solver->refine_each);
        print_to_stdout_and_file("[monodomain_solver] Derefining each %d time steps\n", the_monodomain_solver->derefine_each);

        char *max_dx, *max_dy, *max_dz;

        max_dx = shget(options->domain_config->config_data, "maximum_dx");
        max_dy = shget(options->domain_config->config_data, "maximum_dy");
        max_dz = shget(options->domain_config->config_data, "maximum_dz");

        print_to_stdout_and_file("[domain] Domain maximum Space Discretization: dx %s um, dy %s um, dz %s um\n", max_dx, max_dy, max_dz);


        print_to_stdout_and_file("[monodomain_solver] The adaptivity will start in time: %lf ms\n",
                                 the_monodomain_solver->start_adapting_at);
    }

    if(options->linear_system_solver_config)
    {
        print_linear_system_solver_config_values(options->linear_system_solver_config);
        print_to_stdout_and_file(LOG_LINE_SEPARATOR);
    }

    if(options->save_mesh_config)
    {
        print_save_mesh_config_values(options->save_mesh_config);
        print_to_stdout_and_file(LOG_LINE_SEPARATOR);
    }

    if(options->stim_configs)
    {

        size_t num_stims = shlen(options->stim_configs);

        if(num_stims == 1)
            print_to_stdout_and_file("[stim] Stimulus configuration:\n");
        else
            print_to_stdout_and_file("[stim] Stimuli configuration:\n");

        for(int i = 0; i < num_stims; i++) {

            struct string_voidp_hash_entry e = options->stim_configs[i];
            print_to_stdout_and_file("Stimulus name: %s\n", e.key);
            print_stim_config_values((struct config*) e.value);
            print_to_stdout_and_file(LOG_LINE_SEPARATOR);

        }
    }

    if(options->purkinje_stim_configs)
    {

        size_t num_stims = shlen(options->purkinje_stim_configs);

        if(num_stims == 1)
            print_to_stdout_and_file("[stim_purkinje] Stimulus configuration:\n");
        else
            print_to_stdout_and_file("[stim_purkinje] Stimuli configuration:\n");

        for(int i = 0; i < num_stims; i++) {

            struct string_voidp_hash_entry e = options->purkinje_stim_configs[i];
            print_to_stdout_and_file("Stimulus name: %s\n", e.key);
            print_stim_config_values((struct config*) e.value);
            print_to_stdout_and_file(LOG_LINE_SEPARATOR);

        }
    }

    if (options->domain_config)
    {
        print_domain_config_values(options->domain_config);
        print_to_stdout_and_file(LOG_LINE_SEPARATOR);
    }

    if (options->purkinje_config)
    {
        print_purkinje_config_values(options->purkinje_config);
        print_to_stdout_and_file (LOG_LINE_SEPARATOR);
    }

    if(options->extra_data_config)
    {
        print_extra_data_config_values(options->extra_data_config);
        print_to_stdout_and_file(LOG_LINE_SEPARATOR);
    }

    if(options->update_monodomain_config)
    {
        print_update_monodomain_config_values(options->update_monodomain_config);
        print_to_stdout_and_file(LOG_LINE_SEPARATOR);
    }

    if(options->assembly_matrix_config)
    {
        print_assembly_matrix_config_values(options->assembly_matrix_config);
        print_to_stdout_and_file(LOG_LINE_SEPARATOR);
    }
}

void configure_monodomain_solver_from_options(struct monodomain_solver *the_monodomain_solver,
                                              struct user_options *options) {

    assert(the_monodomain_solver);
    assert(options);

    the_monodomain_solver->num_threads = options->num_threads;
    the_monodomain_solver->final_time = options->final_time;

    the_monodomain_solver->refine_each = options->refine_each;
    the_monodomain_solver->derefine_each = options->derefine_each;
    the_monodomain_solver->refinement_bound = options->ref_bound;
    the_monodomain_solver->derefinement_bound = options->deref_bound;

    the_monodomain_solver->abort_on_no_activity = options->abort_no_activity;

    the_monodomain_solver->calc_activation_time = options->calc_activation_time;
    the_monodomain_solver->print_conductivity = options->print_conductivity_map;
    the_monodomain_solver->print_min_vm = options->print_min_vm_map;
    the_monodomain_solver->print_max_vm = options->print_max_vm_map;
    the_monodomain_solver->print_apd = options->print_apd_map;

    the_monodomain_solver->dt = options->dt_pde;

    the_monodomain_solver->beta = options->beta;
    the_monodomain_solver->cm = options->cm;
    the_monodomain_solver->start_adapting_at = options->start_adapting_at;
}

// ----------------------------------------------------------------------------------------------------------------------------------------
// NEW FUNCTIONS !
void calculate_activation_time (const real_cpu cur_time, const real_cpu dt, const uint32_t n_active, struct cell_node **ac,\
                        struct ode_solver *the_ode_solver)
{
    // V^n+1
    //uint32_t n_active;
    //struct cell_node **ac;

    int n_odes = the_ode_solver->model_data.number_of_ode_equations;
    const double apd_percentage = 0.9;

    // V^n+1/2
    real *sv = the_ode_solver->sv;

    int i;

    if(the_ode_solver->gpu)
    {
#ifdef COMPILE_CUDA
        uint32_t max_number_of_cells = the_ode_solver->original_num_cells;
        real *vms;
        size_t mem_size = max_number_of_cells * sizeof(real);

        vms = (real *)malloc(mem_size);
        check_cuda_errors(cudaMemcpy(vms, sv, mem_size, cudaMemcpyDeviceToHost));

#pragma omp parallel for
        for(i = 0; i < n_active; i++)
        {
            real v_new = vms[ac[i]->sv_position];
            real v_old = (real)ac[i]->v;

            real dvdt = (v_new - v_old) / dt;

            // Activation time
            if(dvdt > ac[i]->max_dvdt)
            {
                ac[i]->max_dvdt = dvdt;
                ac[i]->activation_time = cur_time;
            }

            // APD
            if (v_old < ac[i]->min_v)
                ac[i]->min_v = v_old;
            
            if (v_old > ac[i]->max_v)
            {
             	ac[i]->max_v = v_old;

                ac[i]->v_threashold = ac[i]->min_v + (ac[i]->max_v - ac[i]->min_v)*(1.0-apd_percentage);

                ac[i]->apd = ac[i]->threashold_time - ac[i]->activation_time;

                ac[i]->after_peak = true;

            }
            if (v_old < ac[i]->v_threashold && ac[i]->after_peak)
            {
             	ac[i]->threashold_time = cur_time;
                ac[i]->apd = ac[i]->threashold_time - ac[i]->activation_time;
                ac[i]->after_peak = false;

            }

        }

        free(vms);
#endif
    }
    else
    {
#pragma omp parallel for
        for(i = 0; i < n_active; i++)
        {
            real v_new = sv[ac[i]->sv_position * n_odes];
            real v_old = (real)ac[i]->v;

            real dvdt = (v_new - v_old) / dt;

            // Activation time
            if ( (dvdt > ac[i]->max_dvdt) )
            {

                ac[i]->max_dvdt = dvdt;
                ac[i]->activation_time = cur_time;
            }

            // APD
            if (v_old < ac[i]->min_v)
                ac[i]->min_v = v_old;
            
            if (v_old > ac[i]->max_v)
            {
             	ac[i]->max_v = v_old;

                ac[i]->v_threashold = ac[i]->min_v + (ac[i]->max_v - ac[i]->min_v)*(1.0-apd_percentage);

                ac[i]->apd = ac[i]->threashold_time - ac[i]->activation_time;

                ac[i]->after_peak = true;

            }
            if (v_old < ac[i]->v_threashold && ac[i]->after_peak)
            {
             	ac[i]->threashold_time = cur_time;
                ac[i]->apd = ac[i]->threashold_time - ac[i]->activation_time;
                ac[i]->after_peak = false;

            }

        }
    }

}

void print_activation_time (struct grid *the_grid)
{
    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    int i;

    for(i = 0; i < n_active; i++)
    {
        print_to_stdout_and_file("Cell %i -- Activation time = %g ms\n",i,ac[i]->activation_time);
    }
}

// TODO: Change the LINEAR_SYSTEM macro to (struct config *config, const uint32_t num_active_cells, struct cell_node **ac, uint32_t *number_of_iterations, real_cpu *error)
void linear_system_solver_purkinje (struct config *config, struct grid *the_grid, uint32_t *number_of_iterations, real_cpu *error)
{
    bool jacobi_initialized = false;
    bool bcg_initialized = false;
    bool use_preconditioner = false;
    int max_its = 50;
    real_cpu tol = 1e-16;

    real_cpu  rTr,
            r1Tr1,
            pTAp,
            alpha,
            beta,
            precision = tol,
            rTz,
            r1Tz1;

    // Use the Purkinje linked-list
    uint32_t num_active_cells = the_grid->the_purkinje->num_active_purkinje_cells;
    struct cell_node** ac = the_grid->the_purkinje->purkinje_cells;

    *error = 1.0;
    *number_of_iterations = 1;

    //__________________________________________________________________________
    //Computes int_vector A*x, residue r = b - Ax, scalar rTr = r^T * r and
    //sets initial search direction p.

    rTr = 0.0;
    rTz = 0.0;

    struct element element;
    uint32_t i;

    #pragma omp parallel for private (element) reduction(+:rTr,rTz)
    for (i = 0; i < num_active_cells; i++) {

        if(CG_INFO(ac[i]) == NULL) {
            INITIALIZE_CONJUGATE_GRADIENT_INFO(ac[i]);
        }

        struct element *cell_elements = ac[i]->elements;
        ac[i]->Ax = 0.0;

        size_t max_el = arrlen(cell_elements);

        for(int el = 0; el < max_el; el++) {
            element = cell_elements[el];
            ac[i]->Ax += element.value * element.cell->v;
        }

        CG_R(ac[i]) = ac[i]->b - ac[i]->Ax;
        if(use_preconditioner) {
            real_cpu value = cell_elements[0].value;
            if(value == 0.0) value = 1.0;
            CG_Z(ac[i]) = (1.0/value) * CG_R(ac[i]); // preconditioner
            rTz += CG_R(ac[i]) * CG_Z(ac[i]);
            CG_P(ac[i]) = CG_Z(ac[i]);
        }
        else {
            CG_P(ac[i]) = CG_R(ac[i]);
        }

        real_cpu r = CG_R(ac[i]);

        rTr += r * r;

    }

    *error = rTr;


    //__________________________________________________________________________
    //Conjugate gradient iterations.
    if( *error >= precision ) {
        while( *number_of_iterations < max_its ) {
            //__________________________________________________________________
            // Computes Ap and pTAp. Uses Ax to store Ap.
            pTAp = 0.0;

            #pragma omp parallel for private(element) reduction(+ : pTAp)
            for (i = 0; i < num_active_cells; i++) {

                ac[i]->Ax = 0.0;
                struct element *cell_elements = ac[i]->elements;

                size_t max_el = arrlen(cell_elements);
                for(int el = 0; el < max_el; el++) {
                    element = cell_elements[el];
                    ac[i]->Ax += element.value * CG_P(element.cell);
                }

                pTAp += CG_P(ac[i]) * ac[i]->Ax;
            }

            //__________________________________________________________________
            // Computes alpha.
            if(use_preconditioner) {
                alpha = rTz/pTAp;
            }
            else {
                alpha = rTr/pTAp;
            }
            //__________________________________________________________________


            r1Tr1 = 0.0;
            r1Tz1 = 0.0;

            // Computes new value of solution: u = u + alpha*p.
            #pragma omp parallel for reduction (+:r1Tr1,r1Tz1)
            for (i = 0; i < num_active_cells; i++) {
                ac[i]->v += alpha * CG_P(ac[i]);

                CG_R(ac[i]) -= alpha * ac[i]->Ax;

                real_cpu r = CG_R(ac[i]);

                if(use_preconditioner) {
                    real_cpu value = ac[i]->elements[0].value;
                    if(value == 0.0) value = 1.0;
                    CG_Z(ac[i]) = (1.0/value) * r;
                    r1Tz1 += CG_Z(ac[i]) * r;
                }
                r1Tr1 += r * r;
            }
            //__________________________________________________________________
            //Computes beta.
            if(use_preconditioner) {
                beta = r1Tz1/rTz;
            }
            else {
                beta = r1Tr1/rTr;
            }

            *error = r1Tr1;

            *number_of_iterations = *number_of_iterations + 1;
            if( *error <= precision ) {
                break;
            }
            //__________________________________________________________________
            //Computes int_vector p1 = r1 + beta*p and uses it to upgrade p.
            #pragma omp parallel for
            for (i = 0; i < num_active_cells; i++) {
                if(use_preconditioner) {
                    CG_P1(ac[i]) = CG_Z(ac[i]) + beta * CG_P(ac[i]);
                }
                else {
                    CG_P1(ac[i]) = CG_R(ac[i]) + beta * CG_P(ac[i]);
                }
                CG_P(ac[i]) = CG_P1(ac[i]);
            }

            rTz = r1Tz1;
            rTr = r1Tr1;

        }

    }//end of conjugate gradient iterations.

}//end conjugateGradient() function.

void compute_pmj_current_purkinje_to_tissue (struct ode_solver *the_ode_solver, struct grid *the_grid, struct terminal *the_terminals)
{
    assert(the_ode_solver);
    assert(the_grid);
    assert(the_terminals);

    // Tissue solution
    struct cell_node **ac = the_grid->active_cells;
    uint32_t n_active = the_grid->num_active_cells;
    real *sv = the_ode_solver->sv;
    real initial_v = the_ode_solver->model_data.initial_v;
    uint32_t nodes = the_ode_solver->model_data.number_of_ode_equations;

    // Purkinje solution
    struct cell_node **ac_purkinje = the_grid->the_purkinje->purkinje_cells;

    // Purkinje coupling parameters
    real rpmj = the_grid->the_purkinje->the_network->rpmj;
    real pmj_scale = the_grid->the_purkinje->the_network->pmj_scale;

    // TODO: Switch this to a region around the terminal
    rpmj *= pmj_scale;
    real Gpmj = 1.0 / rpmj;

    if(the_ode_solver->gpu)
    {
    #ifdef COMPILE_CUDA

        uint32_t max_number_of_cells = the_ode_solver->original_num_cells;
        real *vms;
        size_t mem_size = max_number_of_cells * sizeof(real);

        vms = (real *)malloc(mem_size);

        if(the_grid->adaptive)
            check_cuda_errors(cudaMemcpy(vms, sv, mem_size, cudaMemcpyDeviceToHost));

        #pragma omp parallel for
        for(uint32_t i = 0; i < n_active; i++)
        {
            vms[ac[i]->sv_position] = (real)ac[i]->v;
        }

        uint32_t num_of_purkinje_terminals = the_grid->the_purkinje->the_network->number_of_terminals;
        for (uint32_t i = 0; i < num_of_purkinje_terminals; i++)
        {
            //printf("[map_purkinje_solution_to_tissue] Terminal %u -- Tissue index = %u -- Purkinje index = %u\n",i,the_terminals[i].endocardium_index,the_terminals[i].purkinje_index);

            uint32_t tissue_index = the_terminals[i].endocardium_index;
            uint32_t purkinje_index = the_terminals[i].purkinje_index;

            // Compute the PMJ current
            real Ipmj = Gpmj * (vms[tissue_index] - ac_purkinje[purkinje_index]->v);

            // Add this current to the RHS of this cell
            ac[tissue_index]->b -= Ipmj;
        }

        check_cuda_errors(cudaMemcpy(sv, vms, mem_size, cudaMemcpyHostToDevice));
        free(vms);
    #endif
    }
    else
    {
        uint32_t num_of_purkinje_terminals = the_grid->the_purkinje->the_network->number_of_terminals;
        for (uint32_t i = 0; i < num_of_purkinje_terminals; i++)
        {
            //printf("Terminal %u -- Tissue index = %u -- Purkinje index = %u\n",i,the_terminals[i].endocardium_index,the_terminals[i].purkinje_index);

            uint32_t tissue_index = the_terminals[i].endocardium_cell->sv_position;
            uint32_t purkinje_index = the_terminals[i].purkinje_index;

            // Compute the PMJ current
            real Ipmj = Gpmj * (sv[tissue_index*nodes] - ac_purkinje[purkinje_index]->v);
            //printf("ac[tissue_index]->b = %g || Ipmj = %g || ac[tissue_index]->b = %g\n",ac[tissue_index]->b,Ipmj,ac[tissue_index]->b+Ipmj);

            // Add this current to the RHS of this cell
            ac[tissue_index]->b -= Ipmj;
        }
    }
}

void compute_pmj_current_tissue_to_purkinje (struct ode_solver *the_purkinje_ode_solver, struct grid *the_grid, struct terminal *the_terminals)
{
    assert(the_purkinje_ode_solver);
    assert(the_grid);
    assert(the_terminals);

    // Tissue solution
    struct cell_node **ac = the_grid->active_cells;
    uint32_t n_active = the_grid->num_active_cells;

    // Purkinje solution
    struct cell_node **ac_purkinje = the_grid->the_purkinje->purkinje_cells;
    uint32_t n_active_purkinje = the_grid->the_purkinje->num_active_purkinje_cells;
    real *sv = the_purkinje_ode_solver->sv;
    real initial_v = the_purkinje_ode_solver->model_data.initial_v;
    uint32_t nodes = the_purkinje_ode_solver->model_data.number_of_ode_equations;

    real rpmj = the_grid->the_purkinje->the_network->rpmj;
    real Gpmj = 1.0 / rpmj;

    if(the_purkinje_ode_solver->gpu)
    {
    #ifdef COMPILE_CUDA

        uint32_t max_number_of_cells = the_purkinje_ode_solver->original_num_cells;
        real *vms;
        size_t mem_size = max_number_of_cells * sizeof(real);

        vms = (real *)malloc(mem_size);

        check_cuda_errors(cudaMemcpy(vms, sv, mem_size, cudaMemcpyDeviceToHost));

        //#pragma omp parallel for
        //for(uint32_t i = 0; i < n_active_purkinje; i++)
        //{
        //    vms[ac_purkinje[i]->sv_position] = (real)ac_purkinje[i]->v;
        //}

        uint32_t num_of_purkinje_terminals = the_grid->the_purkinje->the_network->number_of_terminals;
        for (uint32_t i = 0; i < num_of_purkinje_terminals; i++)
        {
            //printf("[map_purkinje_solution_to_tissue] Terminal %u -- Tissue index = %u -- Purkinje index = %u\n",i,the_terminals[i].endocardium_index,the_terminals[i].purkinje_index);

            uint32_t tissue_index = the_terminals[i].endocardium_index;
            uint32_t purkinje_index = the_terminals[i].purkinje_index;

            // Compute the PMJ current
            real Ipmj = Gpmj * (vms[purkinje_index] - ac[tissue_index]->v);

            // Add this current to the RHS of this cell
            ac_purkinje[purkinje_index]->b -= Ipmj;
        }

        //check_cuda_errors(cudaMemcpy(sv, vms, mem_size, cudaMemcpyHostToDevice));
        free(vms);
    #endif
    }
    else
    {
        uint32_t num_of_purkinje_terminals = the_grid->the_purkinje->the_network->number_of_terminals;
        for (uint32_t i = 0; i < num_of_purkinje_terminals; i++)
        {
            //printf("Terminal %u -- Tissue index = %u -- Purkinje index = %u\n",i,the_terminals[i].endocardium_index,the_terminals[i].purkinje_index);

            uint32_t tissue_index = the_terminals[i].endocardium_cell->sv_position;
            uint32_t purkinje_index = the_terminals[i].purkinje_index;

            // Compute the PMJ current
            real Ipmj = Gpmj * (sv[purkinje_index*nodes] - ac[tissue_index]->v);
            //printf("ac[tissue_index]->b = %g || Ipmj = %g || ac[tissue_index]->b = %g\n",ac[tissue_index]->b,Ipmj,ac[tissue_index]->b+Ipmj);

            // Add this current to the RHS of this cell
            ac_purkinje[purkinje_index]->b -= Ipmj;
        }
    }
}

double** allocate_matrix_LU (const uint32_t num_rows, const uint32_t num_cols)
{
    double **a = (double**)malloc(sizeof(double*)*num_rows);
    a[0] = (double*)malloc(sizeof(double)*num_rows*num_cols);
    // Hacking the addresses of the lines
    for (int i = 1; i < num_rows; i++)
        a[i] = a[0] + num_cols*i;
    return a;
}

void lu_decomposition (double **lu, uint32_t pivot[], struct cell_node **ac, const uint32_t n_active)
{

// Init all elements to zero
    uint32_t num_rows = n_active;
    uint32_t num_cols = n_active;
    for (uint32_t i = 0; i < num_rows; i++)
        for (uint32_t j = 0; j < num_cols; j++)
            lu[i][j] = 0.0;

// Fill the non-zero elements
    for (uint32_t i = 0; i < n_active; i++)
    {
        struct element *cell_elements = ac[i]->elements;
        size_t max_elements = arrlen(cell_elements);

        uint32_t row = i;

        for(size_t j = 0; j < max_elements; j++)
        {
            uint32_t col = cell_elements[j].column;
            double value = cell_elements[j].value;

            lu[row][col] = value;
        }
    }

// Apply the LU decomposition algorithm

    uint32_t pivot_line;
    double Amax;

    for (uint32_t i = 0; i < num_rows; i++)
        pivot[i] = i;

    // For each line
    for (int i = 0; i < num_rows-1; i++)
    {
        choose_pivot(lu,num_rows,&pivot_line,&Amax,i);
        // The pivot element changed ? If so we need to switch lines
        if (pivot_line != i)
        {
            switch_lines(lu,num_cols,pivot,pivot_line,i);
        }

    // Check if the matrix is not singular
    // If not, apply a Gaussian Elimination on the current line
    if (fabs(lu[i][i]) != 0.0)
    {
        double r = 1 / lu[i][i];
        for (int j = i+1; j < num_rows; j++)
        {
            double m = lu[j][i] * r;

            // Write the multiplier on the inferior triangular matrix L
            lu[j][i] = m;

            // Write the results on the superior triangular matrix U
            for (int k = i+1; k < num_cols; k++)
            {
                lu[j][k] -= (m * lu[i][k]);
            }
        }
    }
}

}

void choose_pivot (double **lu, const uint32_t num_rows, uint32_t *pivot_line, double *Amax, const uint32_t i)
{
    *pivot_line = i;
    *Amax = fabs(lu[i][i]);

    // For all the lines below 'i' search for the maximum value in module
    for (uint32_t j = i+1; j < num_rows; j++)
    {
        double value = fabs(lu[j][i]);
        if (value > *Amax)
        {
            *Amax = value;
            *pivot_line = j;
        }
    }
}

void switch_lines (double **lu, const uint32_t num_cols, uint32_t pivot[], const int pivot_line, const int i)
{
    // Copy the role line
    // TO DO: maybe a pointer swap will be faster ...
    for (int j = 0; j < num_cols; j++)
    {
        double tmp = lu[i][j];
        lu[i][j] = lu[pivot_line][j];
        lu[pivot_line][j] = tmp;
    }
    int m = pivot[i];
    pivot[i] = pivot[pivot_line];
    pivot[pivot_line] = m;
}

void linear_system_solver_purkinje_lu(double **lu, uint32_t pivot[], struct grid *the_grid)
{

    uint32_t n_active = the_grid->the_purkinje->num_active_purkinje_cells;
    struct cell_node **ac = the_grid->the_purkinje->purkinje_cells;

    uint32_t num_rows = n_active;
    uint32_t num_cols = n_active;

    double y[num_rows];

/*
    printf("Before\n");
    for (uint32_t i = 0; i < n_active; i++)
        printf("Cell %u -- V = %g\n",i,ac[i]->v);
*/
    // Forward substitution using the pivot array
    uint32_t k = pivot[0];
    y[0] = ac[k]->b;
    for (uint32_t i = 1; i < num_rows; i++)
    {
        double sum = 0.0;
        //for (uint32_t j = 0; j <= i; j++)
        for (uint32_t j = 0; j < i; j++)
        {
            sum += lu[i][j] * y[j];
        }
        k = pivot[i];
        y[i] = (ac[k]->b - sum);
    }

    // Backward substitution
    ac[num_rows-1]->v = y[num_rows-1] / lu[num_rows-1][num_rows-1];
    for (int i = num_rows-2; i >= 0; i--)
    {
        double sum = 0.0;
        for (int j = i+1; j < num_cols; j++)
            sum += lu[i][j] * ac[j]->v;
        ac[i]->v = (y[i] - sum) / lu[i][i];
    }

/*
    printf("After\n");
    for (uint32_t i = 0; i < n_active; i++)
        printf("Cell %u -- V = %g\n",i,ac[i]->v);

    exit(EXIT_SUCCESS);
*/
}

/*
void solve_explicit_purkinje_diffusion (struct grid *the_grid, struct ode_solver *the_purkinje_ode_solver, struct monodomain_solver *the_monodomain_solver)
{
    assert(the_grid);
    assert(the_purkinje_ode_solver);
    assert(the_monodomain_solver);

    struct graph *the_network = the_grid->the_purkinje->the_network;
    real_cpu beta = the_monodomain_solver->beta;
    real_cpu cm = the_monodomain_solver->cm;
    real_cpu dt_pde = the_monodomain_solver->dt;

    uint32_t n_active = the_grid->the_purkinje->num_active_purkinje_cells;
    struct cell_node **ac = the_grid->the_purkinje->purkinje_cells;

    int n_equations_cell_model = the_purkinje_ode_solver->model_data.number_of_ode_equations;
    real *sv = the_purkinje_ode_solver->sv;

    real_cpu alpha, multiplier;
    struct node *n;
    struct edge *e;

    n = the_network->list_nodes;
    while (n != NULL)
    {
        uint32_t src_index = n->id;
        real_cpu sigma_x = ac[src_index]->sigma.x;
        real_cpu v_old = sv[ ac[src_index]->sv_position * n_equations_cell_model ];

        real update = -( n->num_edges * sv[ac[src_index]->sv_position * n_equations_cell_model] );

        alpha = ALPHA(beta,cm,dt_pde,ac[src_index]->discretization.x,ac[src_index]->discretization.y,ac[src_index]->discretization.z);
        multiplier = (sigma_x * ac[src_index]->discretization.y * ac[src_index]->discretization.z)/(alpha * ac[src_index]->discretization.x);

        e = n->list_edges;
        while (e != NULL)
        {
            uint32_t dest_index = e->id;
            //real_cpu sigma_x2 = ac[dest_index]->sigma.x;
            //real sigma_x = (2.0f * sigma_x1 * sigma_x2) / (sigma_x1 + sigma_x2);

            update += sv[ ac[dest_index]->sv_position * n_equations_cell_model ];

            e = e->next;
        }

        update *= multiplier;

        ac[src_index]->v = v_old + update;

        n = n->next;
    }

}
*/
