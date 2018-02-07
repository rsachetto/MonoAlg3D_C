/*
 * opts.h
 *
 *  Created on: 27/05/2011
 *      Author: sachetto
 */

#ifndef MONOALG3D_CONFIG_PARSER_H
#define MONOALG3D_CONFIG_PARSER_H

#ifdef _MSC_VER
#include "../../getopt/getopt.h"
#else
#include <getopt.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "../../alg/grid/grid.h"
#include "../output_utils.h"
#include "stim_config_hash.h"
#include "domain_config.h"
#include "extra_data_config.h"


#define SIGMA_X 1400
#define SIGMA_Y 1500
#define SIGMA_Z 1600
#define START_REFINING 1700
#define DOMAIN_OPT 1800
#define EXTRA_DATA_OPT 1900
#define STIM_OPT 2000
#define DRAW_OPT 3000
#define BETA 4000
#define CM 5000

struct user_options {
    double final_time;				/*-f option */
    bool final_time_was_set;
    char *out_dir_name;	            /*-o option */
    bool out_dir_name_was_set;
    bool adaptive;	                /*-a option */
    bool adaptive_was_set;
    int print_rate;	            	/*-p option */
    bool print_rate_was_set;
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
    bool binary;                    /*-y option*/
    bool binary_was_set;

    bool draw;


    double sigma_x;
    bool sigma_x_was_set;
    double sigma_y;
    bool sigma_y_was_set;
    double sigma_z;
    bool sigma_z_was_set;
    double beta;
    bool beta_was_set;
    double cm;
    bool cm_was_set;

    double start_adapting_at;
    bool start_adapting_at_was_set;
    char *config_file;              /*-c option*/

    struct stim_config_hash *stim_configs;
    struct domain_config *domain_config;
    struct extra_data_config *extra_data_config;


};


/* Display program usage, and exit.
 */
void display_usage( char** argv );
struct user_options * new_user_options();
void parse_options(int argc, char**argv, struct user_options *user_args);
void get_config_file(int argc, char**argv, struct user_options *user_args);
int parse_config_file(void* user, const char* section, const char* name, const char* value);

void configure_grid_from_options(struct grid* grid, struct user_options *options);


void free_user_options(struct user_options *s);

void issue_overwrite_warning (const char *var, const char *old_value, const char *new_value, const char *config_file);

#endif /* MONOALG3D_CONFIG_PARSER_H */
