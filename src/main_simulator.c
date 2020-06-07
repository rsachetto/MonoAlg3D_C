#include <string.h>

#include "alg/grid/grid.h"
#include "3dparty/ini_parser/ini.h"
#include "monodomain/monodomain_solver.h"
#include "ode_solver/ode_solver.h"
#include "3dparty/sds/sds.h"
#include "utils/file_utils.h"
#include "config_helpers/config_helpers.h"
#include "logger/logger.h"

#ifdef COMPILE_GUI
    #include "gui/gui.h"
#endif


void configure_simulation(int argc, char **argv, struct user_options **options, struct monodomain_solver **monodomain_solver,  struct ode_solver **ode_solver, struct grid **the_grid ) {

    *options = new_user_options();
    *the_grid = new_grid();
    *monodomain_solver = new_monodomain_solver();
    *ode_solver = new_ode_solver();

    // First we have to get the config file path
    get_config_file(argc, argv, *options);

    if((*(options))->config_file) {
        // Here we parse the config file
        if(ini_parse((*(options))->config_file, parse_config_file, *options) < 0) {
            fprintf(stderr, "Error: Can't load the config file %s\n", (*(options))->config_file);
            exit(EXIT_FAILURE);
        }
    }
    else {
        fprintf(stderr, "\nError: The config file is mandatory.\n\n");
        display_usage(argv);
        exit(EXIT_FAILURE);
    }

    // The command line options always overwrite the config file
    parse_options(argc, argv, *options);
    //This variable is from file_utils.h
    set_no_stdout( (*(options))->quiet);

    // Create the output dir and the logfile
    if((*(options))->save_mesh_config) {

        char *out_dir_name = NULL;
        GET_PARAMETER_VALUE_CHAR_OR_USE_DEFAULT(out_dir_name, (*(options))->save_mesh_config->config_data, "output_dir");

        if(out_dir_name) {
            sds buffer_log = sdsnew("");
            sds buffer_ini = sdsnew("");

            bool remove_older_simulation_dir = false;
            GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(remove_older_simulation_dir, (*(options))->save_mesh_config->config_data, "remove_older_simulation");

            if (remove_older_simulation_dir) {
                remove_directory(out_dir_name);
            }

            create_dir(out_dir_name);
            buffer_log = sdscatfmt(buffer_log, "%s/outputlog.txt", out_dir_name);
            open_logfile(buffer_log);

            log_to_stdout_and_file("Command line to reproduce this simulation:\n");
            for (int i = 0; i < argc; i++) {
                log_to_stdout_and_file("%s ", argv[i]);
            }

            log_to_stdout_and_file("\n");

            buffer_ini = sdscatfmt(buffer_ini, "%s/original_configuration.ini", out_dir_name);

            log_to_stdout_and_file("For reproducibility purposes the configuration file was copied to file: %s\n",
                                     buffer_ini);

            //moved to monodomain solver
            //cp_file(buffer_ini, (*(options))->config_file);

            sdsfree(buffer_log);
            sdsfree(buffer_ini);
        }
    }

    configure_ode_solver_from_options(*ode_solver, *options);
    configure_monodomain_solver_from_options(*monodomain_solver, *options);
    configure_grid_from_options(*the_grid, *options);

};

void free_current_simulation_resources(struct user_options *options, struct monodomain_solver *monodomain_solver,  struct ode_solver *ode_solver, struct grid *the_grid) {
    clean_and_free_grid(the_grid);
    free_ode_solver(ode_solver);
    free(monodomain_solver);
    free_user_options(options);
    close_logfile();
}

#ifdef COMPILE_GUI
void init_gui_config(struct gui_config *gui_config, struct user_options *options) {

    gui_config->config_name = strdup(options->config_file);
    gui_config->grid_info.alg_grid = NULL;
    gui_config->max_v = options->max_v;
    gui_config->min_v = options->min_v;

    if(gui_config->min_v == 0) gui_config->min_v = 0.1f;

    gui_config->simulating = false;
    gui_config->time = 0.0;

    gui_config->adaptive = options->adaptive;
    gui_config->final_time = options->final_time;
    gui_config->dt = options->dt_pde;

    gui_config->exit = false;
    gui_config->restart = false;

    gui_config->draw_type = DRAW_SIMULATION;
    gui_config->error_message = NULL;
    gui_config->grid_info.loaded = false;
    gui_config->int_scale = false;
}
#endif

int main(int argc, char **argv) {

    struct user_options *options = NULL;
    struct grid *the_grid;
    struct monodomain_solver *monodomain_solver = NULL;
    struct ode_solver *ode_solver = NULL;

    configure_simulation(argc, argv, &options, &monodomain_solver, &ode_solver, &the_grid);

#ifndef COMPILE_CUDA
    if(ode_solver->gpu) {
        log_to_stdout_and_file("Cuda runtime not found in this system. Fallbacking to CPU solver!!\n");
        ode_solver->gpu = false;
    }
#endif

#ifndef COMPILE_GUI
    if(options->show_gui) {
        log_to_stdout_and_file("OpenGL not found. The output will not be draw!!\n");
        options->show_gui = false;
    }
#endif

    int np = monodomain_solver->num_threads;

    if(np == 0)
        np = 1;

#if defined(_OPENMP)
    omp_set_num_threads(np);
#endif

    //If COMPILE_GUI is not set this is always false. See above.
    if(options->show_gui) {

        #ifdef COMPILE_GUI //If this is defined so OMP is also defined

        omp_set_nested(true);

        omp_init_lock(&gui_config.draw_lock);
        omp_init_lock(&gui_config.sleep_lock);
        init_gui_config(&gui_config, options);

        OMP(parallel sections num_threads(2))
        {
            OMP(section)
            {
                init_and_open_gui_window();
            }

            OMP(section)
            {
                int result = solve_monodomain(monodomain_solver, ode_solver, the_grid, options);

                while (result == RESTART_SIMULATION || result == SIMULATION_FINISHED) {
                    if(result == RESTART_SIMULATION) {
                        free_current_simulation_resources(options, monodomain_solver, ode_solver, the_grid);
                        configure_simulation(argc, argv, &options, &monodomain_solver, &ode_solver, &the_grid);
                        init_gui_config(&gui_config, options);
                        result = solve_monodomain(monodomain_solver, ode_solver, the_grid, options);
                    }

                    if(gui_config.restart) result = RESTART_SIMULATION;

                    if(gui_config.exit)  {
                        free_current_simulation_resources(options, monodomain_solver, ode_solver, the_grid);
                        break;
                    }
                }
            }
        }

        #endif //COMPILE_GUI
    } else {
        solve_monodomain(monodomain_solver, ode_solver, the_grid, options);
        free_current_simulation_resources(options, monodomain_solver, ode_solver, the_grid);
    }

    return EXIT_SUCCESS;
}
