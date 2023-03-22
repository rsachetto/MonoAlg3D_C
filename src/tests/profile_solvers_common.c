//
// Created by sachetto on 10/03/23.
//

void profile_solver(bool preconditioner, char *method_name, char *init_name, char *end_name, struct grid *grid, int nt, struct elapsed_times *times) {

#if defined(_OPENMP)
    omp_set_num_threads(nt);
#endif

    struct config *linear_system_solver_config;

    linear_system_solver_config = alloc_and_init_config_data();
    linear_system_solver_config->main_function_name = strdup(method_name);

    if(init_name)  linear_system_solver_config->init_function_name = strdup(init_name);
    if(end_name)  linear_system_solver_config->end_function_name = strdup(end_name);

    shput_dup_value(linear_system_solver_config->config_data, "tolerance", "1e-16");

    if (preconditioner)
        shput_dup_value(linear_system_solver_config->config_data, "use_preconditioner", "yes");
    else
        shput_dup_value(linear_system_solver_config->config_data, "use_preconditioner", "no");

    shput_dup_value(linear_system_solver_config->config_data, "max_iterations", "200");

    uint32_t n_iter;

    init_config_functions(linear_system_solver_config, "./shared_libs/libdefault_linear_system_solver.so", "linear_system_solver");

    struct time_info ti = ZERO_TIME_INFO;

    struct stop_watch etime;
    init_stop_watch(&etime);

    start_stop_watch(&etime);
    CALL_INIT_LINEAR_SYSTEM(linear_system_solver_config, grid, false);
    times->init_time = (double)stop_stop_watch(&etime);

    double error;

    start_stop_watch(&etime);
    ((linear_system_solver_fn*)linear_system_solver_config->main_function)(&ti, linear_system_solver_config, grid, grid->num_active_cells, grid->active_cells, &n_iter, &error);
    times->run_time = (double)stop_stop_watch(&etime);

    start_stop_watch(&etime);
    CALL_END_LINEAR_SYSTEM(linear_system_solver_config);
    times->end_time = (double)stop_stop_watch(&etime);

    free_config_data(linear_system_solver_config);

}