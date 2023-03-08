#include "common.h"
#include "../3dparty/ini_parser/ini.h"
#include "../utils/file_utils.h"

struct user_options *load_options_from_file(char *config_file) {
    // Here we parse the config file

    struct user_options *options = new_user_options();

    if(config_file) {
        options->config_file = strdup(config_file);

        if(ini_parse(options->config_file, parse_config_file, options) < 0) {
            fprintf(stderr, "Error: Can't load the config file %s\n", options->config_file);
            return NULL;
        }

        return options;
    }

    return NULL;

}

int run_simulation_with_config(struct user_options *options, char *out_dir) {

    struct grid *the_grid;
    the_grid = new_grid();

    struct monodomain_solver *monodomain_solver;
    monodomain_solver = new_monodomain_solver();

    struct ode_solver *ode_solver;
    ode_solver = new_ode_solver();

    shput_dup_value(options->save_mesh_config->config_data, "output_dir", out_dir);

    char *out_dir_name = NULL;
    GET_PARAMETER_STRING_VALUE_OR_USE_DEFAULT(out_dir_name, options->save_mesh_config, "output_dir");

    // Create the output dir and the logfile
    if(out_dir_name) {
        remove_directory(out_dir_name);
        create_dir(out_dir_name);

        sds buffer_log = sdsempty();
        buffer_log = sdscatfmt(buffer_log, "%s/outputlog.txt", out_dir_name);
        open_logfile(buffer_log);

        sdsfree(buffer_log);

    }
    else {
        return 0;
    }

    configure_ode_solver_from_options(ode_solver, options);
    configure_monodomain_solver_from_options(monodomain_solver, options);
    configure_grid_from_options(the_grid, options);

#ifndef COMPILE_CUDA
    if(ode_solver->gpu) {
        log_warn("Cuda runtime not found in this system. Fallbacking to CPU solver!!\n");
        ode_solver->gpu = false;
    }
#endif

    int np = monodomain_solver->num_threads;

    if(np == 0)
        np = 1;

#if defined(_OPENMP)
    omp_set_num_threads(np);
#endif

    set_no_stdout(true);

    solve_monodomain(monodomain_solver, ode_solver, the_grid, options, NULL);

    clean_and_free_grid(the_grid);
    free_ode_solver(ode_solver);
    free(monodomain_solver);

    close_logfile();

    return 1;
}



