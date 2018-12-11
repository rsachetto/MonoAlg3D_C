#include "alg/grid/grid.h"
#include "ini_parser/ini.h"
#include "monodomain/monodomain_solver.h"
#include "monodomain/ode_solver.h"
#include "monodomain/output_utils.h"
#include "string/sds.h"
#include "utils/file_utils.h"

#ifdef COMPILE_OPENGL
#include "draw/draw.h"
#endif

int main(int argc, char **argv) {

    struct user_options *options;
    options = new_user_options();

    struct grid *the_grid;
    the_grid = new_grid();

    struct monodomain_solver *monodomain_solver;
    monodomain_solver = new_monodomain_solver();

    struct ode_solver *ode_solver;
    ode_solver = new_ode_solver();

    // First we have to get the config file path
    get_config_file(argc, argv, options);

    if(options->config_file) {
        // Here we parse the config file
        if(ini_parse(options->config_file, parse_config_file, options) < 0) {
            fprintf(stderr, "Error: Can't load the config file %s\n", options->config_file);
            return EXIT_FAILURE;
        }
    }

    // The command line options always overwrite the config file
    parse_options(argc, argv, options);

    // Create the output dir and the logfile
    if(options->save_mesh_config && options->save_mesh_config->out_dir_name) {
        sds buffer_log = sdsnew("");
        sds buffer_ini = sdsnew("");

        create_dir(options->save_mesh_config->out_dir_name);
        buffer_log = sdscatfmt(buffer_log, "%s/outputlog.txt", options->save_mesh_config->out_dir_name);
        open_logfile(buffer_log);

        print_to_stdout_and_file("Command line to reproduce this simulation:\n");
        for(int i = 0; i < argc; i++) {
            print_to_stdout_and_file("%s ", argv[i]);
        }

        print_to_stdout_and_file("\n");

        buffer_ini = sdscatfmt(buffer_ini, "%s/original_configuration.ini", options->save_mesh_config->out_dir_name);

        print_to_stdout_and_file("For reproducibility purposes the configuration file was copied to file: %s\n",
                                 buffer_ini);

        cp_file(buffer_ini, options->config_file);

        sdsfree(buffer_log);
        sdsfree(buffer_ini);
    }

    configure_ode_solver_from_options(ode_solver, options);
    configure_monodomain_solver_from_options(monodomain_solver, options);
    configure_grid_from_options(the_grid, options);

#ifndef COMPILE_CUDA
    if(ode_solver->gpu) {
        print_to_stdout_and_file("Cuda runtime not found in this system. Fallbacking to CPU solver!!\n");
        ode_solver->gpu = false;
    }
#endif

#ifndef COMPILE_OPENGL
    if(options->draw) {
        print_to_stdout_and_file("OpenGL not found. The output will not be draw!!\n");
        options->draw = false;
    }
#endif

    int np = monodomain_solver->num_threads;

    if(np == 0)
        np = 1;

#if defined(_OPENMP)
    omp_set_num_threads(np);
#endif

    if(options->draw) {
#ifdef COMPILE_OPENGL

#if defined(_OPENMP)
        omp_set_nested(true);
#endif

#pragma omp parallel sections num_threads(2)
        {
#pragma omp section
            {
                grid_to_draw = NULL;
                init_opengl(argc, argv);
            }

#pragma omp section
            { solve_monodomain(monodomain_solver, ode_solver, the_grid, options); }
        }
#endif
    } else {
        solve_monodomain(monodomain_solver, ode_solver, the_grid, options);
    }

    clean_and_free_grid(the_grid);
    free_ode_solver(ode_solver);

    free(monodomain_solver);

    free_user_options(options);
    close_logfile();

    return EXIT_SUCCESS;
}
