#include "alg/grid/grid.h"
#include "ini_parser/ini.h"
#include "monodomain/monodomain_solver.h"
#include "monodomain/ode_solver.h"
#include "string/sds.h"
#include "utils/file_utils.h"
#include <mpi.h>
#include <string.h>

#ifdef COMPILE_OPENGL
#include "draw/draw.h"
#include "single_file_libraries/stb_ds.h"

#endif

struct changed_parameters {
    char *section;
    char *name;
    real_cpu *values;
};

void configure_new_parameters(sds new_out_dir_name, struct changed_parameters *changed, struct user_options *options, int n, int p) {
    new_out_dir_name = sdscatprintf(new_out_dir_name, "_%s_%lf", changed[n].name, changed[n].values[p]);

    sds char_value = sdscatprintf(sdsempty(), "%lf", changed[n].values[p]);

    if(strcmp("extra_data", changed[n].section) == 0) {
        if(options->extra_data_config) {
            shput(options->extra_data_config->config_data.config, changed[n].name,
                  strdup(char_value));
        }
    }
    else if(strcmp("domain", changed[n].section) == 0) {

        if(options->domain_config) {

            if (strcmp("start_dx", changed[n].name) == 0) {
                options->domain_config->start_dx = changed[n].values[p];
            } else if (strcmp("start_dy", changed[n].name) == 0) {
                options->domain_config->start_dy = changed[n].values[p];
            } else if (strcmp("start_dz", changed[n].name) == 0) {
                options->domain_config->start_dz = changed[n].values[p];
            } else if (strcmp("maximum_dx", changed[n].name) == 0) {
                options->domain_config->max_dx = changed[n].values[p];
            } else if (strcmp("maximum_dy", changed[n].name) == 0) {
                options->domain_config->max_dy = changed[n].values[p];
            } else if (strcmp("maximum_dz", changed[n].name) == 0) {
                options->domain_config->max_dz = changed[n].values[p];
            } else {
                shput(options->domain_config->config_data.config, changed[n].name,
                      strdup(char_value));
            }
        }
    }
    else if(strcmp("assembly_matrix", changed[n].section) == 0) {
        if(options->assembly_matrix_config) {
            shput(options->assembly_matrix_config->config_data.config,
                                            changed[n].name, strdup(char_value));
        }
    }
    else if(strcmp("linear_system_solver", changed[n].section) == 0) {
        if(options->linear_system_solver_config) {
            shput(options->linear_system_solver_config->config_data.config,
                                            changed[n].name, strdup(char_value));
        }
    }


    sdsfree(char_value);

}

int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);

    int rank, num_proccess, num_max_proc;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proccess);

    int simulation_number_start;
    int num_simulations;
    int num_par_change;

    char *output_folder;

    struct string_hash_entry *config_to_change = NULL;
    char *entire_config_file = NULL;
    long config_file_size;

    if(rank == 0) {

        struct batch_options *batch_options;
        batch_options = new_batch_options();

        // TODO: maybe we want to get the config file first...
        parse_batch_options(argc, argv, batch_options);

        if(ini_parse(batch_options->batch_config_file, parse_batch_config_file, batch_options) < 0) {
            fprintf(stderr, "Error: Can't load the config file %s\n", batch_options->batch_config_file);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        entire_config_file = read_entire_file(batch_options->initial_config, &config_file_size);

        if(entire_config_file == NULL) {
            fprintf(stderr, "Error: Can't load the config file %s\n", batch_options->initial_config);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        config_to_change = batch_options->config_to_change;

        if(num_proccess > batch_options->num_simulations) {
            num_max_proc = batch_options->num_simulations;
        } else {
            num_max_proc = num_proccess;
        }

        num_simulations = batch_options->num_simulations / num_max_proc;
        int last_rank_extra = batch_options->num_simulations % num_max_proc;

        simulation_number_start = 0;
        int size_out_folder = (int)strlen(batch_options->output_folder) + 1;

        create_dir(batch_options->output_folder);

        num_par_change = batch_options->num_par_change;

        // Send the data to all other process
        for(int i = 1; i < num_proccess; i++) {

            simulation_number_start += num_simulations;
            int num_sims = num_simulations;


            if(i == num_max_proc - 1) {
                num_sims += last_rank_extra;
            }

            if(i >= batch_options->num_simulations) {
                num_sims = 0;
            }

            int n = (int)shlen(config_to_change);

            MPI_Send(&num_sims, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&num_par_change, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&simulation_number_start, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&size_out_folder, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(batch_options->output_folder, size_out_folder, MPI_CHAR, i, 0, MPI_COMM_WORLD);
            MPI_Send(&n, 1, MPI_INT, i, 0, MPI_COMM_WORLD);

            for(int k = 0; k < n; k++) {

                struct string_hash_entry e = config_to_change[k];

                int key_value_sizes[2];
                key_value_sizes[0] = (int)strlen(e.key) + 1;
                key_value_sizes[1] = (int)strlen(e.value) + 1;
                MPI_Send(key_value_sizes, 2, MPI_INT, i, 0, MPI_COMM_WORLD);

                MPI_Send(e.key, key_value_sizes[0], MPI_CHAR, i, 0, MPI_COMM_WORLD);
                MPI_Send(e.value, key_value_sizes[1], MPI_CHAR, i, 0, MPI_COMM_WORLD);

            }

            MPI_Send(&config_file_size, 1, MPI_LONG, i, 0, MPI_COMM_WORLD);
            MPI_Send(entire_config_file, (int)config_file_size, MPI_CHAR, i, 0, MPI_COMM_WORLD);
        }

        simulation_number_start = 0;
        output_folder = batch_options->output_folder;
    }

    else {
        int size_out_folder;
        int num_options;

        MPI_Recv(&num_simulations, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&num_par_change, 1, MPI_INT,  0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&simulation_number_start, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Recv(&size_out_folder, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        output_folder = malloc((size_t)size_out_folder);

        MPI_Recv(output_folder, size_out_folder, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        create_dir(output_folder);

        MPI_Recv(&num_options, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for(int o = 0; o < num_options; o++) {
            int key_value_sizes[2];

            MPI_Recv(key_value_sizes, 2, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            char *key = malloc(key_value_sizes[0]);
            char *value = malloc(key_value_sizes[1]);

            MPI_Recv(key, key_value_sizes[0], MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(value, key_value_sizes[1], MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            shput(config_to_change, key, strdup(value));

            free(key);
            free(value);
        }

        MPI_Recv(&config_file_size, 1, MPI_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        entire_config_file = malloc((size_t)config_file_size);
        MPI_Recv(entire_config_file, (int)config_file_size, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
    }

    if(num_simulations == 0) {
        MPI_Finalize();
        return EXIT_SUCCESS;
    }

    struct user_options *options;
    options = new_user_options();

    struct grid *the_grid;


    struct monodomain_solver *monodomain_solver;
    monodomain_solver = new_monodomain_solver();

    struct ode_solver *ode_solver;


    if(entire_config_file) {
        // Here we parse the config file
        if(ini_parse_string(entire_config_file, parse_config_file, options) < 0) {
            fprintf(stderr, "Error parsing config file\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
    }

    if(!options->save_mesh_config) {
        options->save_mesh_config->out_dir_name = strdup("batch_run");
    }

    struct changed_parameters *changed = NULL;

    size_t config_n = shlen(config_to_change);

    for(int k = 0; k < config_n; k++) {

            struct string_hash_entry e = config_to_change[k];

            struct changed_parameters c;

            int count;
            sds *section_name = sdssplit(e.key, "|", &count);

            c.section = strdup(section_name[0]);
            c.name = strdup(section_name[1]);
            sdsfreesplitres(section_name, count);


            sds *value_range = sdssplit(e.value, "|", &count);

            real_cpu value = strtod(value_range[0], NULL);
            real_cpu inc  = strtod(value_range[1], NULL);
            sdsfreesplitres(value_range, count);

            c.values = malloc(sizeof(real_cpu)*num_par_change);

            for(int n = 0; n < num_par_change; n++) {
                c.values[n] = value + n*inc;
            }

            arrput(changed, c);

    }

    char *initial_out_dir_name  = strdup(options->save_mesh_config->out_dir_name);

    //Parse the modification directives
    for(int s = simulation_number_start; s < simulation_number_start+num_simulations; s++) {

        for (int p = 0; p < num_par_change; p++) {
            the_grid = new_grid();
            ode_solver = new_ode_solver();

            sds new_out_dir_name =
                    sdscatprintf(sdsempty(), "%s/%s_run_%d", output_folder, initial_out_dir_name, s);

            for(int n = 0; n < config_n; n++) {
                configure_new_parameters(new_out_dir_name, changed, options, n, p);
            }

            free(options->save_mesh_config->out_dir_name);
            options->save_mesh_config->out_dir_name = strdup(new_out_dir_name);
            sdsfree(new_out_dir_name);

            printf("Rank %d, performing simulation %d and saving in %s\n", rank, simulation_number_start,
                   options->save_mesh_config->out_dir_name);

            // Create the output dir and the logfile
            if (options->save_mesh_config && options->save_mesh_config->out_dir_name) {
                sds buffer_log = sdsnew("");
                sds buffer_ini = sdsnew("");

                create_dir(options->save_mesh_config->out_dir_name);
                buffer_log = sdscatfmt(buffer_log, "%s/outputlog.txt", options->save_mesh_config->out_dir_name);
                open_logfile(buffer_log);

                print_to_stdout_and_file("Command line to reproduce this simulation:\n");
                for (int i = 0; i < argc; i++) {
                    print_to_stdout_and_file("%s ", argv[i]);
                }

                print_to_stdout_and_file("\n");

                buffer_ini =
                        sdscatfmt(buffer_ini, "%s/original_configuration.ini", options->save_mesh_config->out_dir_name);

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

            int nt = monodomain_solver->num_threads;
            options->draw = false;

            if (nt == 0)
                nt = 1;

#if defined(_OPENMP)
            omp_set_num_threads(nt);
#endif

            solve_monodomain(monodomain_solver, ode_solver, the_grid, options);

            clean_and_free_grid(the_grid);
            free_ode_solver(ode_solver);
            close_logfile();
        }
    }

    free(monodomain_solver);
    free_user_options(options);
    free(initial_out_dir_name);

    MPI_Finalize();

    return EXIT_SUCCESS;
}
