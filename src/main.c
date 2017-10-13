#include <stdio.h>
#include <unitypes.h>
#include "grid/grid/grid.h"
#include "solvers/ode_solver.h"
#include "solvers/monodomain_solver.h"
#include "utils/config_parser.h"
#include "utils/ini_parser/ini.h"

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
            return 1    ;
        }
    }

    //The command line options always overwrite the config file
    parse_options(argc, argv, options);

    configure_ode_solver_from_options(ode_solver, options);
    configure_monodomain_solver_from_options(edp_solver, options);
    configure_output_from_options(output_info, options);
    configure_grid_from_options(the_grid, options);

    init_ode_solver_with_cell_model(ode_solver);

    //TODO: we have to define how to handle the stimuli
    ode_solver->stim_duration = 2.0;
    ode_solver->stim_start = 0.0;
    ode_solver->stim_current = -50.0f;

#ifndef COMPILE_CUDA
    if(ode_solver->gpu) {
        printf("Cuda runtime not found in this system. Fallbacking to CPU solver!!\n");
        ode_solver->gpu = false;
    }
#endif

    //TODO: this should be an user provided function!
    initialize_grid_with_benchmark_mesh(the_grid, edp_solver->start_h);

    solve_monodomain(the_grid, edp_solver, ode_solver, output_info);

    free_output_utils(output_info);
    clean_and_free_grid(the_grid);
    free_ode_solver(ode_solver);

    free(edp_solver);
    free(options);

    return 0;
}

