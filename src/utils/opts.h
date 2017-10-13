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

struct user_options {
    double final_time;				/*-f option */
    double side_lenght;			    /*-l option */
    const char *out_dir_name;	    /*-o option */
    bool adaptive;	                /*-a option */
    int print_rate;	            	/*-p option */
    double start_h;					/*-i option*/
    double max_h;					/*-x option*/
    int max_its;					/*-t option*/
    double ref_bound;				/*-r option*/
    double deref_bound;				/*-b option*/
    double dt_edp;						/*-z option*/
    double dt_edo;				/*-e option*/
    int num_threads;                  /*-q option*/
    bool gpu;                       /*-j option*/
    bool use_jacobi;
    int refine_each;
    int derefine_each;
    int gpu_id;                     /*-G option*/

};

static const struct option longOpts[] = {
        { "side-lenght", required_argument, NULL, 'l' },
        { "out-dir", required_argument , NULL, 'o' },
        { "adaptive", no_argument, NULL, 'a' },
        { "num_threads", required_argument, NULL, 'n' },
        { "gpu", no_argument, NULL, 'g' },
        { "print-rate",required_argument , NULL, 'p' },
        { "start-h", required_argument, NULL, 's' },
        { "max-h", required_argument, NULL, 'x' },
        { "max-cg-its", required_argument, NULL, 'm' },
        { "refinement-bound", required_argument, NULL, 'r' },
        { "derefinement-bound", required_argument, NULL, 'd' },
        { "dt-edp", required_argument, NULL, 'z' },
        { "dt-edo", required_argument, NULL, 'e' },
        { "final-time", required_argument, NULL, 'f' },
        { "use-jacobi", no_argument, NULL, 'j' },
        { "refine-each", required_argument, NULL, 'R'},
        { "derefine-each", required_argument, NULL, 'D'},
        { "gpu-id", required_argument, NULL, 'G'},
        { "help", no_argument, NULL, 'h' },
        { NULL, no_argument, NULL, 0 }
};

/* Display program usage, and exit.
 */
void display_usage( char** argv );
struct user_options * new_user_options();
void parse_options(int argc, char**argv, struct user_options *user_args);

#endif /* OPTS_H_ */
