/*
 * opts.h
 *
 *  Created on: 27/05/2011
 *      Author: sachetto
 */

#ifndef OPTS_H_
#define OPTS_H_

#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#define JACOBI 1000
#define STIM 2000

struct user_args {
    bool saveToFile;				/*-s option */
    double stim_dur;				/*-d option */
    double final_time;				/*-f option */
    double sideLenght;			    /*-l option */
    const char *out_dir_name;	    /*-o option */
    bool adaptive;	                /*-a option */
    int print_rate;	            	/*-p option */
    double start_h;					/*-i option*/
    double min_h;					/*-m option*/
    double max_h;					/*-x option*/
    int max_its;					/*-t option*/
    double ref_bound;				/*-r option*/
    double deref_bound;				/*-b option*/
    double stim_cur;				/*-c option*/
    double dt;						/*-z option*/
    double min_dt_edo;				/*-e option*/
    double max_dt_edo;				/*-v option*/
    double abs_tol; 				/*-g option*/
    double rel_tol;					/*-n option*/
    bool use_rabbit;				/*-w option*/
    bool use_mouse; 				/*-k option*/
    bool use_human; 				/*-H option*/
    bool use_plain; 				/*-P option*/
    bool use_plain_with_sphere;     /*-S option*/
    int method;                     /*-u option*/
    int parallel;                  /*-q option*/
    bool gpu;                       /*-j option*/
    bool benchmark;					/*-y option*/
    bool use_jacobi;
    int refine_each;
    int derefine_each;
    double phi;
    int gpu_id;                     /*-G option*/
    bool stim_file;                 /*-S option*/
    bool humanSub;                  /*-U option*/
    bool scarfile;                  /*-C option*/

};

static const struct option longOpts[] = {
        { "save-to-file", no_argument, NULL, 's' },
        { "stim-duration", required_argument, NULL, 'd' },
        { "side-lenght", required_argument, NULL, 'l' },
        { "out-dir", required_argument , NULL, 'o' },
        { "adaptive", no_argument, NULL, 'a' },
        { "parallel", required_argument, NULL, 'q' },
        { "gpu", no_argument, NULL, 'j' },
        { "print-rate",required_argument , NULL, 'p' },
        { "start-h", required_argument, NULL, 'i' },
        { "min-h", required_argument, NULL, 'm' },
        { "max-h", required_argument, NULL, 'x' },
        { "max-cg-its", required_argument, NULL, 't' },
        { "refinement-bound", required_argument, NULL, 'r' },
        { "derefinement-bound", required_argument, NULL, 'b' },
        { "stim-current", required_argument, NULL, 'c' },
        { "dt-edp", required_argument, NULL, 'z' },
        { "min-dt-edo", required_argument, NULL, 'e' },
        { "max-dt-edo", required_argument, NULL, 'v' },
        { "abs-tol", required_argument, NULL, 'g' },
        { "rel-tol", required_argument, NULL, 'n' },
        { "method", required_argument, NULL, 'u' },
        { "use-rabbit-heart", no_argument, NULL, 'w' },
        { "use-mouse-heart", no_argument, NULL, 'k' },
        { "use-human-heart", no_argument, NULL, 'H' },
        { "use-human-submesh", no_argument, NULL, 'U' },
        { "benchmark", no_argument, NULL, 'y' },
        { "plain", no_argument, NULL, 'P' },
        { "plain-sphere", no_argument, NULL, 'S' },
        { "final-time", required_argument, NULL, 'f' },
        { "use-jacobi", no_argument, NULL, JACOBI },
        { "refine-each", required_argument, NULL, 'R'},
        { "derefine-each", required_argument, NULL, 'D'},
        { "phi", required_argument, NULL, 'I'},
        { "gpu-id", required_argument, NULL, 'G'},
        { "stim-file", no_argument, NULL, STIM},
        { "scar-file", no_argument, NULL, 'C'},
        { "help", no_argument, NULL, 'h' },
        { NULL, no_argument, NULL, 0 }
};

static const char *optString = "e:v:g:n:u:z:c:r:b:t:i:m:x:sd:f:p:l:o:aykq:jyhHPICU:G:S?";

/* Display program usage, and exit.
 */
void display_usage( char** argv ) {

    printf("Usage: %s [options]\n\n", argv[0]);
    printf("Options:\n");
    printf("--final-time, -f [simulation final time]. Simulation final time. Default UNDEFINED.\n");
    printf("--save-to-file, -s. No argument required. Save the output. Default No save.\n");
    printf("--stim-duration, -d [stim_duration]. Stimulus duration (ms). Default: 2 ms \n");
    printf("--side-lenght, -l [length]. Tissue side lenght (micro m). Default: 12800 micro m \n");
    printf("--out-dir, -o [output-dir]. Simulation output directory. Default: ./ \n");
    printf("--adaptive, -a. No argument required. Use adaptivity. Default No use.\n");
    printf("--print-rate, -p [output-print-rate]. Output print rate (in number of iterations). Default: 1 \n");
    printf("--start-h, -i [start-h]. Initial space discretization (micro m). Default: 12.5 micro m \n");
    printf("--min-h, -m [min-h].Minimum space discretizati'on (micro m). Default: 12.5 micro m \n");
    printf("--max-h, -x [max-h].Maximum space discretization (micro m). Default: 200.0 micro m \n");
    printf("--max-cg-its, -t [max-its]. Maximum number of CG iterations. Default: number of volumes \n");
    printf("--refinement-bound, -r [ref-bound]. ALG refinement bound (Vm variantion between time steps). Default: 5.0 \n");
    printf("--derefinement-bound, -b [deref-bound]. ALG derefinement bound (Vm variantion between time steps). Default: 1.0 \n");
    printf("--stim-current, -c [stim-cur]. Stimulus current. Default: Cell model Vm \n");
    printf("--dt-edp, -z [dt]. Simulation time discretization (PDE). Default: 0.01 \n");
    printf("--min-dt-edo, -e [dt]. Minimum ODE time discretization (using time adaptivity. Default: 0.01 \n");
    printf("--max-dt-edo, -v [dt]. Maximum ODE time discretization (using time adaptivity. Default: 0.05 \n");
    printf("--abs-tol, -g [absolute_tolerance]. Absolute tolerance to the euler with adaptive time step methods. Default: 1e-6. \n");
    printf("--rel-tol, -n [relative_tolerance]. Relative tolerance to the euler with adaptive time step methods. Default: 1e-6. \n");
    printf("--method, -u [ode_method].EDO monodomain_solver method. Default: 0\n\tOptions:\n\t 0: Euler Method\n\t 1: Euler Method with adaptive time step (Formula)\n");
    printf("--use-rabbit-heart, -w. Use rabbit mesh file (rabheart.alg). Default false\n");
    printf("--use-mouse-heart, -k. Use mouse mesh file (mouse.alg). Default false\n");
    printf("--use-human-heart, -H. Use human mesh file (human.alg). Default false\n");
    printf("--plain, -P. Plain Fibrosis test domain. Default false\n");
    printf("--plain-sphere, -S. Plain  Fibrosis test domain with a sphere in the center. Default false\n");
    printf("--benchmark, -y. Run N-Version benchmark. Default false\n");
    printf("--parallel, -q [num-threads]. Solve using OpenMP. Default: No \n");
    printf("--gpu, -j. Solve ODEs using GPU. Default: No \n");
    printf("--use-jacobi, Use Jacobi Preconditioner. Default: No \n");
    printf("--refine-each, R [ts], Refine each ts timesteps. Default: 1 \n");
    printf("--derefine-each, D [ts], Derefine each ts timesteps. Default: 1 \n");
    printf("--phi, percentage of cells to cut in a fibrotic region. Default: 0 \n");
    printf("--gpu-id, ID of the GPU to be used. Default: 0 \n");
    printf("--stim-file, use a file named stim.pts containing the points to stimulate. \n");
    printf("-C, --scar-file, use a file named scar.pts containing the points where the scar is. \n");
    printf("-U, --use-human-submesh --stim-file, crate a sub mesh containing only the big scar. \n");
    printf("--help, -h. Shows this help and exit \n");
    exit( EXIT_FAILURE );
}

void parse_options(int argc, char**argv, struct user_args *user_args) {
    int opt = 0;

    /* Initialize user_args before we get to work. */
    user_args->saveToFile = false;
    user_args->out_dir_name = "./";
    user_args->sideLenght = 12800.0;
    user_args->adaptive = false;
    user_args->parallel = 0;
    user_args->use_rabbit = false;
    user_args->use_mouse = false;
    user_args->use_human = false;
    user_args->use_plain = false;
    user_args->use_plain_with_sphere = false;
    user_args->humanSub = false;
    user_args->scarfile = false;

    user_args->gpu = false;
    user_args->final_time = 10.0;
    user_args->stim_dur = 2.0;
    user_args->print_rate = 1;
    user_args->start_h = 12.5;
    user_args->min_h = 12.5;
    user_args->max_h = 200.0;
    user_args->max_its = 0;
    user_args->ref_bound = 4.0;
    user_args->deref_bound = 1.0;
    user_args->stim_cur = INFINITY;
    user_args->dt = 0.01;
    user_args->min_dt_edo = 0.01;
    user_args->max_dt_edo = 0.05;
    user_args->abs_tol = 1e-6;
    user_args->rel_tol = 1e-6;
    user_args->method = 0;
    user_args->benchmark = false;
    user_args->use_jacobi = false;
    user_args->refine_each = 1;
    user_args->derefine_each = 1;
    user_args->phi=0;
    user_args->gpu_id = 0;
    user_args->stim_file = false;

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
                user_args->min_dt_edo = atof(optarg);
                break;
            case 'v':
                user_args->max_dt_edo = atof(optarg);
                break;
            case 'g':
                user_args->abs_tol = atof(optarg);
                break;
            case 'n':
                user_args->rel_tol = atof(optarg);
                break;
            case 'u':
                user_args->method = atoi(optarg);
                break;
            case 't':
                user_args->max_its = atoi(optarg);
                break;
            case 's':
                user_args->saveToFile = true;
                break;
            case 'l':
                user_args->sideLenght = atof(optarg);
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
            case 'i':
                user_args->start_h = atof(optarg);
                break;
            case 'm':
                user_args->min_h = atof(optarg);
                break;
            case 'x':
                user_args->max_h = atof(optarg);
                break;
            case 'r':
                user_args->ref_bound = atof(optarg);
                break;
            case 'b':
                user_args->deref_bound = atof(optarg);
                break;
            case 'd':
                user_args->stim_dur= atof(optarg);
                break;
            case 'a':
                user_args->adaptive = true;
                break;
            case 'q':
                if(atoi(optarg) > 0)
                    user_args->parallel = atoi(optarg);
                break;
            case 'j':
                user_args->gpu = true;
                break;
            case 'c':
                user_args->stim_cur= atof(optarg);
                break;
            case 'z':
                user_args->dt = atof(optarg);
                break;
            case 'w':
                user_args->use_rabbit = true;
                break;
            case 'y':
                user_args->benchmark = true;
                break;
            case 'k':
                user_args->use_mouse = true;
                break;
            case 'H':
                user_args->use_human = true;
                break;
            case 'U':
                user_args->humanSub = true;
                break;
            case 'C':
                user_args->scarfile = true;
                break;
            case 'P':
                user_args->use_plain = true;
                break;
            case 'S':
                user_args->use_plain_with_sphere = true;
                break;
            case 'I':
                user_args->phi = atof(optarg);
                break;
            case 'G':
                user_args->gpu_id = atoi(optarg);
                break;
            case JACOBI:
                user_args->use_jacobi = true;
                break;
            case STIM:
                user_args->stim_file = true;
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


    if(user_args->humanSub && !user_args->use_human) {
        user_args->use_human = true;
    }

}

#endif /* OPTS_H_ */
