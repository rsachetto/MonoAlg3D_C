#include "alg/grid/grid.h"
#include "ini_parser/ini.h"
#include "monodomain/monodomain_solver.h"
#include "monodomain/ode_solver.h"
#include "monodomain/output_utils.h"
#include "string/sds.h"
#include "utils/file_utils.h"
#include <mpi.h>
#include <string.h>

#ifdef COMPILE_OPENGL
#include "draw/draw.h"
#endif

int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);

    int rank, num_proccess, num_max_proc;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proccess);

    int simulation_number_start;
    int num_simulations;
    int num_par_change;

    char *output_folder;

    struct string_hash *config_to_change = string_hash_create();
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

            MPI_Send(&num_sims, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&num_par_change, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&simulation_number_start, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&size_out_folder, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(batch_options->output_folder, size_out_folder, MPI_CHAR, i, 0, MPI_COMM_WORLD);
            MPI_Send(&config_to_change->n, 1, MPI_INT, i, 0, MPI_COMM_WORLD);

            for(int k = 0; k < config_to_change->size; k++) {
                for(struct string_elt *e = config_to_change->table[k % config_to_change->size]; e != 0; e = e->next) {

                    int key_value_sizes[2];
                    key_value_sizes[0] = (int)strlen(e->key) + 1;
                    key_value_sizes[1] = (int)strlen(e->value) + 1;
                    MPI_Send(key_value_sizes, 2, MPI_INT, i, 0, MPI_COMM_WORLD);

                    MPI_Send(e->key, key_value_sizes[0], MPI_CHAR, i, 0, MPI_COMM_WORLD);
                    MPI_Send(e->value, key_value_sizes[1], MPI_CHAR, i, 0, MPI_COMM_WORLD);
                }
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

            string_hash_insert(config_to_change, key, value);

            

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


    char *changed_parameters = NULL;



    //Parse the modification directives
    for(int s = simulation_number_start; s < simulation_number_start+num_simulations; s++) {

        for (int p = 0; p < num_par_change; p++) {
            the_grid = new_grid();
            ode_solver = new_ode_solver();

            sds new_out_dir_name =
                    sdscatprintf(sdsempty(), "%s/%s_%d", output_folder, options->save_mesh_config->out_dir_name, s);

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

    MPI_Finalize();

    return EXIT_SUCCESS;
}
