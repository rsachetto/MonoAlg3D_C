#include <stdio.h>
#include "alg/grid/grid.h"
#include "main/ode_solver.h"
#include "main/monodomain_solver.h"
#include "ini_parser/ini.h"
#include "main/config_parser.h"

int main(int argc, char **argv) {

    struct user_options* options;

    //We load the default options first
    options = new_user_options();
    
    struct grid *the_grid;
    struct monodomain_solver *edp_solver;
    struct ode_solver *ode_solver;
    struct output_utils *output_info;

    the_grid   = new_grid();
    edp_solver = new_monodomain_solver();
    ode_solver = new_ode_solver();
    output_info = new_output_utils();

    //First we have to get the config file path
    get_config_file(argc, argv, options);

    if(options->config_file) {
        //Here we parse the config file
        if (ini_parse(options->config_file, parse_config_file, options) < 0) {
            fprintf(stderr, "Error: Can't load the config file %s\n", options->config_file);
            return EXIT_FAILURE;
        }
    }
    else {
        fprintf(stderr, "Error: no config file provided! Some configurations are only available via config file!\n");
        return EXIT_FAILURE;
    }

    //The command line options always overwrite the config file
    parse_options(argc, argv, options);

    configure_ode_solver_from_options(ode_solver, options);
    configure_monodomain_solver_from_options(edp_solver, options);
    configure_output_from_options(output_info, options);
    configure_grid_from_options(the_grid, options);

    init_ode_solver_with_cell_model(ode_solver);

    STIM_CONFIG_HASH_FOR_EACH_KEY_APPLY_FN_IN_VALUE_KEY(options->stim_configs, init_stim_functions);
    init_domain_functions(options->domain_config);

#ifndef COMPILE_CUDA
    if(ode_solver->gpu) {
        printf("Cuda runtime not found in this system. Fallbacking to CPU solver!!\n");
        ode_solver->gpu = false;
    }
#endif

    options->domain_config->set_spatial_domain_fn(the_grid, options->domain_config->config);

    solve_monodomain(the_grid, edp_solver, ode_solver, output_info, options->stim_configs);

    free_output_utils(output_info);
    clean_and_free_grid(the_grid);
    free_ode_solver(ode_solver);

    free(edp_solver);
    free(options);

    STIM_CONFIG_HASH_FOR_EACH_KEY_APPLY_FN_IN_VALUE(options->stim_configs, free_stim_config);
    free(options);

    return EXIT_SUCCESS;
}

