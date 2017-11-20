#include <pthread.h>
#include "alg/grid/grid.h"
#include "monodomain/ode_solver.h"
#include "monodomain/monodomain_solver.h"
#include "ini_parser/ini.h"
#include "utils/logfile_utils.h"

#ifdef COMPILE_OPENGL
#include "draw/draw.h"
#endif

struct mono_args {
    struct user_options* options;
    struct grid *the_grid;
    struct monodomain_solver *monodomain_solver;
    struct ode_solver *ode_solver;
};

void * start_monodomain_in_thread(void *args) {
    struct mono_args* args_local = (struct mono_args*) args;
    solve_monodomain(args_local->monodomain_solver, args_local->ode_solver, args_local->the_grid, args_local->options);
    return NULL;
}

int main(int argc, char **argv) {

    struct user_options* options;

    options = new_user_options();

    struct grid *the_grid;
    struct monodomain_solver *monodomain_solver;
    struct ode_solver *ode_solver;

    the_grid   = new_grid();
    monodomain_solver = new_monodomain_solver();
    ode_solver = new_ode_solver();

    //First we have to get the config file path
    get_config_file(argc, argv, options);

    if(options->config_file) {
        //Here we parse the config file
        if (ini_parse(options->config_file, parse_config_file, options) < 0) {
            fprintf(stderr, "Error: Can't load the config file %s\n", options->config_file);
            return EXIT_FAILURE;
        }
    }

    //The command line options always overwrite the config file
    parse_options(argc, argv, options);

    if(options->out_dir_name) {
        sds buffer = sdsnew("");
        create_dir_if_no_exists(options->out_dir_name);
        buffer = sdscatfmt(buffer,"%s/outputlog.txt", options->out_dir_name);
        open_logfile(buffer);
        sdsfree(buffer);
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

    if(options->draw) {
#ifdef COMPILE_OPENGL
        grid_to_draw = NULL;
        struct mono_args *thread_args = (struct mono_args *) malloc(sizeof(struct mono_args));
        thread_args->options = options;
        thread_args->ode_solver = ode_solver;
        thread_args->the_grid = the_grid;
        thread_args->monodomain_solver = monodomain_solver;

        pthread_t solver_thread;


        /* create a second thread which executes the monodomain solver as the main thread needs to update the Window */
        if(pthread_create(&solver_thread, NULL, start_monodomain_in_thread, thread_args)) {
            fprintf(stderr, "Error creating thread\n");
            return 1;
        }

        init_opengl(argc, argv);
#endif
    }
    else {
        solve_monodomain(monodomain_solver, ode_solver, the_grid, options);
    }

    clean_and_free_grid(the_grid);
    free_ode_solver(ode_solver);

    free(monodomain_solver);

    free_user_options(options);
    close_logfile();

    return EXIT_SUCCESS;
}

