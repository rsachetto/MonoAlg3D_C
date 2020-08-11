/*
 * opts.h
 *
 *  Created on: 27/05/2011
 *      Author: sachetto
 */

#ifndef MONOALG3D_CONFIG_PARSER_H
#define MONOALG3D_CONFIG_PARSER_H

#include <getopt.h>

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "../alg/grid/grid.h"

#include "domain_config.h"
#include "purkinje_config.h"
#include "extra_data_config.h"
#include "assembly_matrix_config.h"
#include "save_mesh_config.h"
#include "linear_system_solver_config.h"
#include "../common_types/common_types.h"
#include "save_state_config.h"
#include "restore_state_config.h"
#include "update_monodomain_config.h"

#define START_REFINING 1700
#define DOMAIN_OPT 1800
#define EXTRA_DATA_OPT 1900
#define STIM_OPT 2000
#define ASSEMBLY_MATRIX_OPT 2100
#define LINEAR_SYSTEM_SOLVER_OPT 2200
#define UPDATE_MONODOMAIN_SOLVER_OPT 2300
#define MODIFY_DOMAIN_OPT 2400
#define ODE_SOLVER_OPT 2500
#define DRAW_OPT 3000
#define SAVE_OPT 3100
#define SAVE_STATE_OPT 3200
#define RESTORE_STATE_OPT 3300
#define MAX_V_OPT 3400
#define MIN_V_OPT 3500
#define BETA 4000
#define CM 5000
#define VISUALIZATION_PAUSED_OPT 5100



struct user_options {
    real_cpu final_time;				/*-f option */
    bool final_time_was_set;
    bool adaptive;	                /*-a option */
    bool adaptive_was_set;

    real_cpu ref_bound;				/*-r option*/
    bool ref_bound_was_set;
    real_cpu deref_bound;				/*-d option*/
    bool deref_bound_was_set;
    real_cpu dt_pde;					/*-z option*/
    bool dt_pde_was_set;

    real_cpu dt_ode;
    bool dt_ode_was_set;

    bool ode_adaptive;
    bool ode_adaptive_was_set;

    real ode_reltol;
    bool ode_reltol_was_set;

    real ode_abstol;
    bool ode_abstol_was_set;

    int num_threads;                /*-n option*/
    bool num_threads_was_set;
    bool gpu;                       /*-g option*/
    bool gpu_was_set;
    int refine_each;                /*-R option*/
    bool refine_each_was_set;
    int derefine_each;              /*-D option*/
    bool derefine_each_was_set;
    int gpu_id;                     /*-G option*/
    bool gpu_id_was_set;
    bool abort_no_activity;         /*-b option*/
    bool abort_no_activity_was_set;

    real_cpu only_abort_after_dt;
    bool only_abort_after_dt_was_set;

    real_cpu vm_threshold;            /*-v option*/
    bool vm_threshold_was_set;

    bool save_activation_time_map;
    bool save_activation_time_map_was_set;

    bool save_apd_map;
    bool save_apd_map_was_set;

    real_cpu purkinje_dt_ode;
    bool purkinje_dt_ode_was_set;

    bool purkinje_gpu;
    bool purkinje_gpu_was_set;

    int purkinje_gpu_id;         
    bool purkinje_gpu_id_was_set;

    bool purkinje_ode_adaptive;
    bool purkinje_ode_adaptive_was_set;
    
    char *purkinje_model_file_path;         
    bool purkinje_model_file_path_was_set;

    char *model_file_path;          /*-k option*/
    bool model_file_path_was_set;

    bool show_gui;

    real_cpu beta;
    bool beta_was_set;

    real_cpu cm;
    bool cm_was_set;

    real_cpu start_adapting_at;
    bool start_adapting_at_was_set;
    char *config_file;              /*-c option*/

    bool quiet; /*-q option*/
    bool quiet_was_set;

    bool start_visualization_unpaused;

    struct string_voidp_hash_entry *stim_configs;
    struct string_voidp_hash_entry *purkinje_stim_configs;
    struct string_voidp_hash_entry *modify_domain_configs;
    
    struct config *domain_config;
    struct config *purkinje_config;
    struct config *extra_data_config;
    struct config *assembly_matrix_config;
    struct config *linear_system_solver_config;
    struct config *save_mesh_config;
    struct config *save_state_config;
    struct config *restore_state_config;
    struct config *update_monodomain_config;

    struct string_hash_entry *ode_extra_config;
    struct string_hash_entry *purkinje_ode_extra_config;

    real_cpu max_v, min_v;

};

struct batch_options {
    char *batch_config_file;
    char *output_folder;
    char *initial_config;
    int num_simulations;
    bool skip_existing_run;
    struct string_hash_entry *config_to_change;
};

struct visualization_options {
    char *input;
    char *files_prefix;
    bool save_activation_only;
    int start_file;
    int step;
    real_cpu max_v, min_v, dt;
};

struct conversion_options {
    char *input;
    char *output;
};

void display_usage( char** argv );
void display_batch_usage(char **argv);
void display_conversion_usage(char **argv);

struct user_options * new_user_options();
struct batch_options * new_batch_options();
struct visualization_options * new_visualization_options();
struct conversion_options * new_conversion_options();

void parse_options(int argc, char**argv, struct user_options *user_args);
void parse_batch_options(int argc, char**argv, struct batch_options *user_args);
void parse_visualization_options(int argc, char**argv, struct visualization_options *user_args);
void parse_conversion_options(int argc, char **argv, struct conversion_options *user_args);

void get_config_file(int argc, char**argv, struct user_options *user_args);

int parse_config_file(void* user, const char* section, const char* name, const char* value);
int parse_batch_config_file(void *user, const char *section, const char *name, const char *value);

void options_to_ini_file(struct user_options *config, char *ini_file_path);

void configure_grid_from_options(struct grid* grid, struct user_options *options);

void free_user_options(struct user_options *s);
void free_batch_options(struct batch_options * options);
void free_visualization_options(struct visualization_options * options);
void free_conversion_options(struct conversion_options *options);

void maybe_issue_overwrite_warning(const char *var, const char *section, const char *old_value, const char *new_value, const char *config_file);
void set_or_overwrite_common_data(struct config* config, const char *key, const char *value, const char *section, const char *config_file);


#endif /* MONOALG3D_CONFIG_PARSER_H */
