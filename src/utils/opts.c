#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>


#include "opts.h"

static const char *optString = "l:p:an:gp:s:x:m:r:d:z:e:f:l:R:D:G:h?";

/* Display program usage, and exit.
 */
void display_usage( char** argv ) {

    printf("Usage: %s [options]\n\n", argv[0]);
    printf("Options:\n");
    printf("--final-time, -f [simulation final time]. Simulation final time. Default UNDEFINED.\n");
    printf("--side-lenght, -l [length]. Tissue side lenght (micro m). Default: 12800 micro m \n");
    printf("--out-dir, -o [output-dir]. Simulation output directory. Default: ./ \n");
    printf("--adaptive, -a. No argument required. Use adaptivity. Default No use.\n");
    printf("--print-rate, -p [output-print-rate]. Output print rate (in number of iterations). Default: 1 \n");
    printf("--start-h, -s [start-h]. Initial space discretization (micro m). Default: 12.5 micro m \n");
    printf("--max-h, -x [max-h].Maximum space discretization (micro m). Default: 200.0 micro m \n");
    printf("--max-cg-its, -m [max-its]. Maximum number of CG iterations. Default: number of volumes \n");
    printf("--refinement-bound, -r [ref-bound]. ALG refinement bound (Vm variantion between time steps). Default: 5.0 \n");
    printf("--derefinement-bound, -d [deref-bound]. ALG derefinement bound (Vm variantion between time steps). Default: 1.0 \n");
    printf("--dt-edp, -z [dt]. Simulation time discretization (PDE). Default: 0.01 \n");
    printf("--dt-edo, -e [dt]. Minimum ODE time discretization (using time adaptivity. Default: 0.01 \n");
    printf("--num_threads, -n [num-threads]. Solve using OpenMP. Default: 1 \n");
    printf("--gpu, -g. Solve ODEs using GPU. Default: No \n");
    printf("--use-jacobi, -j Use Jacobi Preconditioner. Default: No \n");
    printf("--refine-each, -R [ts], Refine each ts timesteps. Default: 1 \n");
    printf("--derefine-each, -D [ts], Derefine each ts timesteps. Default: 1 \n");
    printf("--gpu-id, -G [id], ID of the GPU to be used. Default: 0 \n");
    printf("--help, -h. Shows this help and exit \n");
    exit( EXIT_FAILURE );
}

struct user_options * new_user_options() {

    struct user_options* user_args = (struct user_options*)malloc(sizeof(struct user_options));

    /* Initialize user_options before we get to work. */
    user_args->out_dir_name = "./";
    user_args->side_lenght = 12800.0;
    user_args->adaptive = false;
    user_args->num_threads = 1;
    user_args->gpu = false;
    user_args->final_time = 10.0;
    user_args->print_rate = 1;
    user_args->start_h = 100.0;
    user_args->max_h = 200.0;
    user_args->max_its = 50;
    user_args->ref_bound = 0.11;
    user_args->deref_bound = 0.10;
    user_args->dt_edp = 0.02;
    user_args->dt_edo = 0.01;
    user_args->use_jacobi = false;
    user_args->refine_each = 1;
    user_args->derefine_each = 1;
    user_args->gpu_id = 0;

    return user_args;
}

void parse_options(int argc, char**argv, struct user_options *user_args) {
    int opt = 0;

    int longIndex;

    opt = getopt_long( argc, argv, optString, longOpts, &longIndex );

    while( opt != -1 ) {
        switch( opt ) {

            case 'R':
                user_args->refine_each = atoi(optarg);
                break;
            case 'D':
                user_args->derefine_each = atoi(optarg);
                break;
            case 'e':
                user_args->dt_edo = atof(optarg);
                break;
            case 'm':
                user_args->max_its = atoi(optarg);
                break;
            case 'l':
                user_args->side_lenght = atof(optarg);
                break;
            case 'p':
                user_args->print_rate = atoi(optarg);
                break;
            case 'o':
                user_args->out_dir_name = optarg;
                break;
            case 'f':
                user_args->final_time = atof(optarg);
                break;
            case 's':
                user_args->start_h = atof(optarg);
                break;
            case 'x':
                user_args->max_h = atof(optarg);
                break;
            case 'r':
                user_args->ref_bound = atof(optarg);
                break;
            case 'd':
                user_args->deref_bound = atof(optarg);
                break;
            case 'a':
                user_args->adaptive = true;
                break;
            case 'n':
                if(atoi(optarg) > 0)
                    user_args->num_threads = atoi(optarg);
                break;
            case 'g':
                user_args->gpu = true;
                break;
            case 'z':
                user_args->dt_edp = atof(optarg);
                break;
            case 'G':
                user_args->gpu_id = atoi(optarg);
                break;
            case 'j':
                user_args->use_jacobi = true;
                break;
            case 'h':	/* fall-through is intentional */
            case '?':
                display_usage(argv);
                break;

            default:
                /* You won't actually get here. */
                break;
        }

        opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
    }


}

