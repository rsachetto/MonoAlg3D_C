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
    bool final_time_was_set;
    double side_lenght;			    /*-l option */
    bool side_lenght_was_set;
    char *out_dir_name;	            /*-o option */
    bool out_dir_name_was_set;
    bool adaptive;	                /*-a option */
    bool adaptive_was_set;
    int print_rate;	            	/*-p option */
    bool print_rate_was_set;
    double start_h;					/*-s option*/
    bool start_h_was_set;
    double max_h;					/*-x option*/
    bool max_h_was_set;
    int max_its;					/*-m option*/
    bool max_its_was_set;
    double ref_bound;				/*-r option*/
    bool ref_bound_was_set;
    double deref_bound;				/*-d option*/
    bool deref_bound_was_set;
    double dt_edp;					/*-z option*/
    bool dt_edp_was_set;
    double dt_edo;				    /*-e option*/
    bool dt_edo_was_set;
    int num_threads;                /*-n option*/
    bool num_threads_was_set;
    bool gpu;                       /*-g option*/
    bool gpu_was_set;
    bool use_jacobi;                /*-j option*/
    bool use_jacobi_was_set;
    int refine_each;                /*-R option*/
    bool refine_each_was_set;
    int derefine_each;              /*-D option*/
    bool derefine_each_was_set;
    int gpu_id;                     /*-G option*/
    bool gpu_id_was_set;
    bool abort_no_activity;         /*-b option*/
    bool abort_no_activity_was_set;
    double cg_tol;                  /*-t option*/
    bool cg_tol_was_set;
    char *model_file_path;          /*-k option*/
    bool model_file_path_was_set;
    char *config_file;              /*-c option*/


};

static const struct option longOpts[] = {
        { "config_file", required_argument, NULL, 'c' },
        { "side_lenght", required_argument, NULL, 'l' },
        { "output_dir", required_argument , NULL, 'o' },
        { "use_adaptivity", no_argument, NULL, 'a' },
        { "abort_on_no_activity", no_argument, NULL, 'b' },
        { "num_threads", required_argument, NULL, 'n' },
        { "use_gpu", no_argument, NULL, 'g' },
        { "print_rate",required_argument , NULL, 'p' },
        { "start_discretization", required_argument, NULL, 's' },
        { "maximum_discretization", required_argument, NULL, 'x' },
        { "max_cg_its", required_argument, NULL, 'm' },
        { "cg_tolerance", required_argument, NULL, 't' },
        { "refinement_bound", required_argument, NULL, 'r' },
        { "derefinement_bound", required_argument, NULL, 'd' },
        { "dt_edp", required_argument, NULL, 'z' },
        { "dt_edo", required_argument, NULL, 'e' },
        { "simulation_time", required_argument, NULL, 'f' },
        { "use_preconditioner", no_argument, NULL, 'j' },
        { "refine_each", required_argument, NULL, 'R'},
        { "derefine_each", required_argument, NULL, 'D'},
        { "gpu_id", required_argument, NULL, 'G'},
        { "model_file_path", required_argument, NULL, 'k'},
        { "help", no_argument, NULL, 'h' },
        { NULL, no_argument, NULL, 0 }
};

/* Display program usage, and exit.
 */
void display_usage( char** argv );
struct user_options * new_user_options();
void parse_options(int argc, char**argv, struct user_options *user_args);
void get_config_file(int argc, char**argv, struct user_options *user_args);
int parse_config_file(void* user, const char* section, const char* name, const char* value);

#endif /* OPTS_H_ */
