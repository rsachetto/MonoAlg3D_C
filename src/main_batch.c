#include "alg/grid/grid.h"
#include "ini_parser/ini.h"
#include "monodomain/monodomain_solver.h"
#include "monodomain/ode_solver.h"
#include "monodomain/output_utils.h"
#include "string/sds.h"
#include "utils/logfile_utils.h"
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
    char *output_folder;

    if(rank == 0) {

        struct batch_options *batch_options;
        batch_options = new_batch_options();

        //TODO: maybe we want to get the config file first...
        parse_batch_options(argc, argv, batch_options);

        if(ini_parse(batch_options->batch_config_file, parse_batch_config_file, batch_options) < 0) {
            fprintf(stderr, "Error: Can't load the config file %s\n", batch_options->batch_config_file);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }


        if(num_proccess > batch_options->num_simulations) {
            num_max_proc = batch_options->num_simulations;
        }
        else {
            num_max_proc = num_proccess;
        }

        //This folder needs to be in a shared filesystem
        //TODO: send this to every other proccess....
        create_dir(batch_options->output_folder);

        num_simulations = batch_options->num_simulations/num_max_proc;
        int last_rank_extra = batch_options->num_simulations%num_max_proc;

        simulation_number_start = 0;
        int size_out_folder = (int)strlen(batch_options->output_folder);

        for(int i = 1; i < num_proccess; i++) {

            simulation_number_start += num_simulations;
            int num_sims = num_simulations;

            if(i == num_max_proc - 1) {
                num_sims += last_rank_extra;
            }

            if (i >= batch_options->num_simulations) {
                num_sims = 0;
            }

            MPI_Send(&num_sims,                1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&simulation_number_start, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
            MPI_Send(&size_out_folder, 1, MPI_INT, i, 3, MPI_COMM_WORLD);
            MPI_Send(batch_options->output_folder, size_out_folder, MPI_CHAR, i, 4, MPI_COMM_WORLD);
        }

        simulation_number_start = 0;
        output_folder = batch_options->output_folder;
    }

    else {
        int size_out_folder;

        MPI_Recv(&num_simulations, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&simulation_number_start, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&size_out_folder, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        output_folder = malloc(size_out_folder);

        MPI_Recv(output_folder, size_out_folder, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    if(num_simulations == 0) {
        MPI_Finalize();
        return EXIT_SUCCESS;
    }

    printf("Rank %d, performing %d simulations, starting at index %d and saving in %s\n", rank, num_simulations, simulation_number_start, output_folder);

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

    int nt = monodomain_solver->num_threads;

    if(nt == 0)
        nt = 1;

#if defined(_OPENMP)
    omp_set_num_threads(nt);
#endif

    solve_monodomain(monodomain_solver, ode_solver, the_grid, options);


    clean_and_free_grid(the_grid);
    free_ode_solver(ode_solver);

    free(monodomain_solver);

    free_user_options(options);
    close_logfile();

    MPI_Finalize();

    return EXIT_SUCCESS;
}
