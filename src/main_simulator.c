#include "3dparty/ini_parser/ini.h"
#include "3dparty/sds/sds.h"
#include "3dparty/stb_ds.h"
#include "alg/grid/grid.h"
#include "config_helpers/config_helpers.h"
#include "logger/logger.h"
#include "monodomain/monodomain_solver.h"
#include "ode_solver/ode_solver.h"
#include "utils/file_utils.h"
#include <stdlib.h>
#include <string.h>

#ifdef COMPILE_GUI
#include "gui/gui.h"
#endif

#if defined(COMPILE_CUDA) || defined(COMPILE_SYCL)
#define COMPILE_GPU
#endif

void configure_simulation(int argc, char **argv, struct user_options **options, struct monodomain_solver **monodomain_solver, struct ode_solver **ode_solver,
                          struct grid **the_grid) {

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
    } else {
        fprintf(stderr, "\nError: The config file is mandatory.\n\n");
        display_usage(argv);
        exit(EXIT_FAILURE);
    }

    // The command line options always overwrite the config file
    parse_options(argc, argv, *options);

    set_no_stdout((*(options))->quiet);

    // Create the output dir and the logfile
    if((*(options))->save_mesh_config) {

        char *out_dir_name = NULL;
        bool add_timestamp = false;
        GET_PARAMETER_STRING_VALUE_OR_USE_DEFAULT(out_dir_name, (*(options))->save_mesh_config, "output_dir");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(add_timestamp, (*(options))->save_mesh_config, "add_timestamp");

        if(out_dir_name) {

            if(add_timestamp) {
                char *tmp = get_timestamped_dir_name(out_dir_name);
                if(tmp) {
                    free(out_dir_name);
                    out_dir_name = strdup(tmp);
                    free(tmp);
                    ADD_STRING_PARAMETER_TO_CONFIG("output_dir", out_dir_name, (*(options))->save_mesh_config); // LEAK??
                }
            }

            sds buffer_log = sdsnew("");
            sds buffer_ini = sdsnew("");

            bool remove_older_simulation_dir = false;
            GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(remove_older_simulation_dir, (*(options))->save_mesh_config, "remove_older_simulation");

            if(remove_older_simulation_dir) {
                remove_directory(out_dir_name);
            }

            create_dir(out_dir_name);
            buffer_log = sdscatfmt(buffer_log, "%s/outputlog.txt", out_dir_name);
            open_logfile(buffer_log);

            log_info("Command line to reproduce this simulation:\n");

            for(int i = 0; i < argc; i++) {
                log_msg("%s ", argv[i]);
            }

            log_msg("\n");

            buffer_ini = sdscatfmt(buffer_ini, "%s/original_configuration.ini", out_dir_name);

            log_info("For reproducibility purposes the configuration file was copied to file: %s\n", buffer_ini);

            sdsfree(buffer_log);
            sdsfree(buffer_ini);
        }
    }

    if((*(options))->save_state_config) {
        char *save_state_out_dir_name = NULL;
        GET_PARAMETER_STRING_VALUE_OR_USE_DEFAULT(save_state_out_dir_name, (*(options))->save_state_config, "output_dir");
        if(save_state_out_dir_name != NULL) {
            create_dir(save_state_out_dir_name);
        }
    }

    configure_ode_solver_from_options(*ode_solver, *options);
    configure_monodomain_solver_from_options(*monodomain_solver, *options);
    configure_grid_from_options(*the_grid, *options);
}

void free_current_simulation_resources(struct user_options *options, struct monodomain_solver *monodomain_solver, struct ode_solver *ode_solver,
                                       struct grid *the_grid) {
    clean_and_free_grid(the_grid);
    free_ode_solver(ode_solver);
    free(monodomain_solver);
    free_user_options(options);
    close_logfile();
}

#ifdef COMPILE_GUI
static void init_gui_config_for_simulation(const struct user_options *options, struct gui_shared_info *gui_config, bool only_restart) {

    if(!only_restart) {
        omp_init_nest_lock(&gui_config->draw_lock);
        omp_init_nest_lock(&gui_config->sleep_lock);
    }

    gui_config->config_name = strdup(options->config_file);
    gui_config->grid_info.alg_grid = NULL;
    gui_config->max_v = options->max_v;
    gui_config->min_v = options->min_v;

    if(gui_config->min_v == 0.0f)
        gui_config->min_v = 0.1f;

    gui_config->simulating = false;
    gui_config->time = 0.0f;

    gui_config->adaptive = options->adaptive;
    gui_config->final_time = options->final_time;
    gui_config->dt = options->dt_pde;

    gui_config->exit = false;
    gui_config->restart = false;

    gui_config->draw_type = DRAW_SIMULATION;
    gui_config->message = NULL;
    gui_config->grid_info.loaded = false;
    gui_config->int_scale = false;
    gui_config->ui_scale = 0.0;
}
#endif

int main(int argc, char **argv) {

    struct user_options *options = NULL;

    struct grid *the_grid = NULL;
    struct monodomain_solver *monodomain_solver = NULL;
    struct ode_solver *ode_solver = NULL;

    configure_simulation(argc, argv, &options, &monodomain_solver, &ode_solver, &the_grid);

#ifndef COMPILE_GPU
    if(ode_solver->gpu) {
        log_warn("MonoAlg3D was not compiled with CUDA support. Falling back to CPU solver!!\n");
        ode_solver->gpu = false;
    }
#endif

#ifndef COMPILE_GUI
    if(options->show_gui) {
        log_warn("MonoAlg3D was not compiled with GUI support. Using CLI!!\n");
        options->show_gui = false;
    }
#endif

#ifdef DEBUG_INFO
    log_warn("Running the debug version. If you do not need debug information, build with './build -r' or 'make release' to improve the performace!\n");
#endif

#if defined(_OPENMP)
    int np = monodomain_solver->num_threads;

    if(np == 0)
        np = 1;

    omp_set_num_threads(np);
#endif

    // If COMPILE_GUI is not set this is always false. See above.
    if(options->show_gui) {

#ifdef COMPILE_GUI // If this is defined so OMP is also defined

        struct gui_shared_info *gui_config = MALLOC_ONE_TYPE(struct gui_shared_info);

        omp_set_nested(true);

        init_gui_config_for_simulation(options, gui_config, false);

        OMP(parallel sections num_threads(2)) {
            OMP(section) {
                init_and_open_gui_window(gui_config);
            }

            OMP(section) {
                int result = solve_monodomain(monodomain_solver, ode_solver, the_grid, options, gui_config);

                while(result == RESTART_SIMULATION || result == SIMULATION_FINISHED) {
                    if(result == RESTART_SIMULATION) {
                        free_current_simulation_resources(options, monodomain_solver, ode_solver, the_grid);
                        configure_simulation(argc, argv, &options, &monodomain_solver, &ode_solver, &the_grid);
                        init_gui_config_for_simulation(options, gui_config, true);
                        result = solve_monodomain(monodomain_solver, ode_solver, the_grid, options, gui_config);
                    }

                    if(gui_config->restart)
                        result = RESTART_SIMULATION;

                    if(gui_config->exit) {
                        free_current_simulation_resources(options, monodomain_solver, ode_solver, the_grid);
                        break;
                    }
                }
            }
        }

#endif // COMPILE_GUI
    } else {
        solve_monodomain(monodomain_solver, ode_solver, the_grid, options, NULL);
        free_current_simulation_resources(options, monodomain_solver, ode_solver, the_grid);
    }

    return EXIT_SUCCESS;
}
