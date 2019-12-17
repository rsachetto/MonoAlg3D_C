#include <assert.h>
#include <string.h>

#include "../ini_parser/ini_file_sections.h"
#include "../string/sds.h"
#include "../utils/file_utils.h"
#include "config_parser.h"
#include "stim_config.h"

#include "../single_file_libraries/stb_ds.h"
#include "../config_helpers/config_helpers.h"

static const char *batch_opt_string = "c:h?";
static const struct option long_batch_options[] = {{"config_file", required_argument, NULL, 'c'}};


static const char *visualization_opt_string = "x:m:d:p:v:a:c:s:h?";
static const struct option long_visualization_options[] = {
        {"visualization_max_v", required_argument, NULL, 'x'},
        {"visualization_min_v", required_argument, NULL, 'm'},
        {"dt", required_argument, NULL, 'd'},
        {"show_activation_map", required_argument, NULL, 'a'},
        {"convert_activation_map", required_argument, NULL, 'c'},
        {"prefix", required_argument, NULL, 'p'},
        {"start_at", required_argument, NULL, 's'},
        {"pvd", required_argument, NULL, 'v'},
        {"help", no_argument, NULL, 'h'},
        {NULL, no_argument, NULL, 0}
};


static const char *opt_string = "c:abn:g:m:t:r:d:z:e:f:jR:D:G:k:v:qh?";
static const struct option long_options[] = {
        {"config_file", required_argument, NULL, 'c'},
        {"use_adaptivity", no_argument, NULL, 'a'},
        {"abort_on_no_activity", no_argument, NULL, 'b'},
        {"vm_threshold", required_argument, NULL, 'v'},
        {"num_threads", required_argument, NULL, 'n'},
        {"use_gpu", required_argument, NULL, 'g'},
        {"print_rate", required_argument, NULL, 'p'},
        {"max_cg_its", required_argument, NULL, 'm'},
        {"cg_tolerance", required_argument, NULL, 't'},
        {"refinement_bound", required_argument, NULL, 'r'},
        {"derefinement_bound", required_argument, NULL, 'd'},
        {"dt_pde", required_argument, NULL, 'z'},
        {"dt_ode", required_argument, NULL, 'e'},
        {"simulation_time", required_argument, NULL, 'f'},
        {"use_preconditioner", no_argument, NULL, 'j'},
        {"refine_each", required_argument, NULL, 'R'},
        {"derefine_each", required_argument, NULL, 'D'},
        {"gpu_id", required_argument, NULL, 'G'},
        {"model_file_path", required_argument, NULL, 'k'},
        {"beta", required_argument, NULL, BETA},
        {"cm", required_argument, NULL, CM},
        {"start_adapting_at", required_argument, NULL, START_REFINING},
        {"domain", required_argument, NULL, DOMAIN_OPT}, //Complex option in the format --domain "name='domain', start_dx=250.0" ...
        {"save_result", required_argument, NULL, SAVE_OPT}, //Complex option
        {"assembly_matrix", required_argument, NULL, ASSEMBLY_MATRIX_OPT}, //Complex option
        {"extra_data", required_argument, NULL, EXTRA_DATA_OPT}, //Complex option
        {"stimulus", required_argument, NULL, STIM_OPT}, //Complex option
        {"save_state", required_argument, NULL, SAVE_STATE_OPT}, //Complex option
        {"restore_state", required_argument, NULL, RESTORE_STATE_OPT}, //Complex option
        {"linear_system_solver", required_argument, NULL, LINEAR_SYSTEM_SOLVER_OPT}, //Complex option
        {"update_monodomain", required_argument, NULL, UPDATE_MONODOMAIN_SOLVER_OPT}, //Complex option
        {"visualize", no_argument, NULL, DRAW_OPT},
        {"visualization_max_v", required_argument, NULL, MAX_V_OPT},
        {"visualization_min_v", required_argument, NULL, MIN_V_OPT},
        {"start_simulation_unpaused", no_argument, NULL, VISUALIZATION_PAUSED_OPT},
        {"quiet", no_argument, NULL, 'q'},
        {"help", no_argument, NULL, 'h'},
        {NULL, no_argument, NULL, 0}};


void display_usage(char **argv) {

    printf("Usage: %s [options]\n\n", argv[0]);
    printf("Options:\n");
    printf("--config_file | -c [configuration_file_path]. Simulation configuration file. Mandatory!\n");
    printf("--simulation_time | -f [simulation final time]. Simulation final time. Default 10.\n");
    printf("--use_adaptivity | -a. No argument required. Use adaptivity. Default No use.\n");
    printf("--abort_on_no_activity | -b. No argument required. The simulation will be aborted if no activity is "
           "verified after print_rate time steps. Default false.\n");
    printf("--vm_treshold | -v. V to abort. If abort_on_no_activity is set, this will be used as abort threshold. "
           "Default -86.0.\n");
    printf("--print_rate | -p [output-print-rate]. Output print rate (in number of iterations). Default: 1 \n");
    printf("--max_cg_its | -m [max-its]. Maximum number of CG iterations. Default: number of volumes \n");
    printf("--cg_tolerance | -t [tolerance]. Conjugate Gradiente tolerance. Default: 1e-16 \n");
    printf("--refinement_bound | -r [ref-bound]. ALG refinement bound (Vm variantion between time steps). Default: "
           "5.0 \n");
    printf("--derefinement_bound | -d [deref-bound]. ALG derefinement bound (Vm variantion between time steps). "
           "Default: 1.0 \n");
    printf("--dt_pde | -z [dt]. Simulation time discretization (PDE). Default: 0.01 \n");
    printf("--dt_ode | -e [dt]. Minimum ODE time discretization (using time adaptivity. Default: 0.01 \n");
    printf("--beta. Value of beta for simulation (Default: 0.14 \n");
    printf("--cm. Value cm (Default: 1.0 \n");
    printf("--num_threads | -n [num-threads]. Solve using OpenMP. Default: 1 \n");
    printf("--use_gpu | -g [0|1|yes|no|true|false]. Solve ODEs using GPU. Default: No \n");
    printf("--use_preconditioner | -j Use Jacobi Preconditioner. Default: No \n");
    printf("--refine_each | -R [ts], Refine each ts timesteps. Default: 1 \n");
    printf("--derefine_each | -D [ts], Derefine each ts timesteps. Default: 1 \n");
    printf("--gpu_id | -G [id], ID of the GPU to be used. Default: 0 \n");
    printf("--model_file_path | -k [.so file path], Path of the .so representing the cell model and the ode solver. "
           "Default: NULL \n");
    printf("--domain, Change the domain configuration without changing the .ini file."
           " Example: --domain \"name=domain, start_dx=250.0\"\n");
    printf("--visualize, Draw a iterative 3D output of the simulation. Not recommended for big meshes. Default: "
           "not draw\n");
    printf("--visualization_max_v, maximum value for V. This is only used to configure the color map in the visualization window. Default: -86.0\n");
    printf("--visualization_min_v, minimum value for V. This is only used to configure the color map in the visualization window. Default: 40.0\n");
    printf("--start_simulation_unpaused. When visualizing, starts the simulation unpaused. Default: true (starts paused)\n");

    printf("--help | -h. Shows this help and exit \n");
    exit(EXIT_FAILURE);
}

void display_batch_usage(char **argv) {

    printf("Usage: %s [options]\n\n", argv[0]);
    printf("Options:\n");
    printf("--config_file | -c [configuration_file_path]. Batch simulations configuration file. Default NULL.\n");
    printf("--help | -h. Shows this help and exit \n");
    exit(EXIT_FAILURE);
}

void display_visualization_usage(char **argv) {

    printf("Usage: %s [options] input_folder \n\n", argv[0]);
    printf("Options:\n");
    printf("--visualization_max_v | -x, maximum value for V. Default: -86.0\n");
    printf("--visualization_min_v | -m, minimum value for V. Default: 40.0\n");
    printf("--dt | -d, dt for the simulation. Default: 0\n");
    printf("--show_activation_map | -a activation_map_file, visualize only the activation map file. Default: NULL\n");
    printf("--convert_activation_map | -c activation_map_file, only convert the activation map file to VTU without opening it. Default: NULL\n");
    printf("--prefix | -p, simulation output files prefix . Default: V_it\n");
    printf("--pvd | -v, pvd file. Default: NULL\n");
    printf("--start_at | -s, Visualize starting at file number [n]. Default: 0\n");
    printf("--help | -h. Shows this help and exit \n");
    exit(EXIT_FAILURE);
}

void issue_overwrite_warning(const char *var, const char *section, const char *old_value, const char *new_value, const char *config_file) {
    fprintf(stderr,
            "WARNING: option %s in %s was set in the file %s to %s and is being overwritten "
            "by the command line flag to %s!\n",
            var, section, config_file, old_value, new_value);
}

struct batch_options *new_batch_options() {
    struct batch_options *user_args = (struct batch_options *)calloc(1, sizeof(struct batch_options));
            sh_new_arena(user_args->config_to_change);
            shdefault(user_args->config_to_change, NULL);

    return user_args;
}

void free_batch_options(struct batch_options * options) {
            shfree(options->config_to_change);
    free(options);
}


struct visualization_options *new_visualization_options() {
    struct visualization_options *options = (struct visualization_options *)malloc(sizeof(struct visualization_options));
    options->input_folder = NULL;
    options->max_v = 40.0f;
    options->min_v = -86.0f;
    options->dt = 0.0;
    options->files_prefix = strdup("V_it");
    options->pvd_file = NULL;
    options->activation_map = NULL;
    options->save_activation_only = false;

    return options;
}

void free_visualization_options(struct visualization_options * options) {
    free(options->input_folder);
    free(options->files_prefix);
    free(options->pvd_file);
    free(options);
}

struct user_options *new_user_options() {

    struct user_options *user_args = (struct user_options *)malloc(sizeof(struct user_options));

    user_args->adaptive = false;
    user_args->adaptive_was_set = false;

    user_args->num_threads = 1;
    user_args->num_threads_was_set = false;

    user_args->gpu = false;
    user_args->gpu_was_set = false;

    user_args->final_time = 10.0;
    user_args->final_time_was_set = false;

    user_args->ref_bound = 0.11;
    user_args->ref_bound_was_set = false;

    user_args->deref_bound = 0.10;
    user_args->deref_bound_was_set = false;

    user_args->dt_pde = 0.02;
    user_args->dt_pde_was_set = false;

    user_args->dt_ode = 0.01;
    user_args->dt_ode_was_set = false;

    user_args->refine_each = 1;
    user_args->refine_each_was_set = false;

    user_args->derefine_each = 1;
    user_args->derefine_each_was_set = false;

    user_args->gpu_id = 0;
    user_args->gpu_id_was_set = false;

    user_args->abort_no_activity = false;
    user_args->abort_no_activity_was_set = false;

    user_args->model_file_path = NULL;
    user_args->model_file_path_was_set = false;

    user_args->config_file = NULL;

    user_args->beta = 0.14;
    user_args->beta_was_set = false;

    user_args->cm = 1.0;
    user_args->cm_was_set = false;

    user_args->start_adapting_at = 1.0;
    user_args->start_adapting_at_was_set = false;

    user_args->vm_threshold = -86.0f;
    user_args->vm_threshold_was_set = false;

    user_args->quiet = false;
    user_args->quiet_was_set = false;

    user_args->calc_activation_time = false;
    user_args->calc_activation_time_was_set = false;

    user_args->print_min_vm_map = false;
    user_args->print_min_vm_map_was_set = false;

    user_args->print_max_vm_map = false;
    user_args->print_max_vm_map_was_set = false;

    user_args->print_apd_map = false;
    user_args->print_apd_map_was_set = false;

    user_args->print_conductivity_map = false;
    user_args->print_conductivity_map_was_set = false;

    user_args->stim_configs = NULL;
    user_args->ode_extra_config = NULL;

    sh_new_arena(user_args->stim_configs);
    shdefault(user_args->stim_configs, NULL);

    sh_new_arena(user_args->purkinje_stim_configs);
    shdefault(user_args->purkinje_stim_configs, NULL);

    sh_new_arena(user_args->ode_extra_config);
    shdefault(user_args->ode_extra_config, NULL);

    sh_new_arena(user_args->purkinje_ode_extra_config);
    shdefault(user_args->purkinje_ode_extra_config, NULL);

    user_args->domain_config = NULL;
    user_args->purkinje_config = NULL;
    user_args->extra_data_config = NULL;
    user_args->assembly_matrix_config = NULL;
    user_args->linear_system_solver_config = NULL;
    user_args->save_mesh_config = NULL;
    user_args->save_state_config = NULL;
    user_args->restore_state_config = NULL;
    user_args->update_monodomain_config = NULL;

    user_args->draw = false;
    user_args->max_v = 40.0f;
    user_args->min_v = -86.0f;

    user_args->main_found = false;

    user_args->start_visualization_unpaused = false;

    return user_args;
}


void set_or_overwrite_common_data(struct config* config, const char *key, const char *value, const char *section, const char *config_file) {
    if(strcmp(key, "main_function") == 0) {

        if(section && config_file) {
            if (config->main_function_name_was_set) {
                issue_overwrite_warning("main_function", section, config->main_function_name, value, config_file);
            }
        }
        free(config->main_function_name);
        config->main_function_name = strdup(value);
    }
    else if(strcmp(key, "init_function") == 0) {
        if(section && config_file) {
            if (config->init_function_name_was_set) {
                issue_overwrite_warning("init_function", section, config->init_function_name, value, config_file);
            }
        }
        free(config->init_function_name);
        config->init_function_name = strdup(value);
    }

    else if(strcmp(key, "end_function") == 0) {
        if(section && config_file) {
            if (config->end_function_name_was_set) {
                issue_overwrite_warning("end_function", section, config->end_function_name, value, config_file);
            }
        }
        free(config->end_function_name);
        config->end_function_name = strdup(value);
    }
    
    else if(strcmp(key, "library_file") == 0) {
        if(section && config_file) {
            if (config->library_file_path_was_set) {
                issue_overwrite_warning("library_file", section, config->library_file_path, value, config_file);
            }
        }
        free(config->library_file_path);
        config->library_file_path = strdup(value);
    } else {
        char *opt_value = shget(config->config_data, key);
        if(section && config_file) {
            if (opt_value) {
                issue_overwrite_warning(key, section, opt_value, value, config_file);
            }
        }
        shput_dup_value(config->config_data, key, value);
    }
}

void set_common_data(struct config* config, const char *key, const char *value) {

    if(strcmp(key, "main_function") == 0) {
        config->main_function_name = strdup(value);
        config->main_function_name_was_set = true;
    }
    if(strcmp(key, "init_function") == 0) {
        config->init_function_name = strdup(value);
        config->init_function_name_was_set = true;
    }
    if(strcmp(key, "end_function") == 0) {
        config->end_function_name = strdup(value);
        config->end_function_name_was_set = true;
    }
    else if(strcmp(key, "library_file") == 0) {
        config->library_file_path = strdup(value);
    }
    else {
        shput_dup_value(config->config_data, key, value);
    }

}

void set_stim_config(const char *args, struct string_voidp_hash_entry *stim_configs, const char *config_file) {

    sds extra_config;
    sds *extra_config_tokens;
    int tokens_count;
    extra_config = sdsnew(args);
    extra_config_tokens = sdssplit(extra_config, ",", &tokens_count);
    char *stim_name = NULL;
    char old_value[32];
    char *key, *value;

    assert(stim_configs);

    for(int i = 0; i < tokens_count; i++) {
        extra_config_tokens[i] = sdstrim(extra_config_tokens[i], " ");

        int values_count;
        sds *key_value = sdssplit(extra_config_tokens[i], "=", &values_count);

        if(values_count != 2) {
            fprintf(stderr, "Invalid format for optios %s. Exiting!\n", args);
            exit(EXIT_FAILURE);
        }

        if(strcmp(key_value[0], "name") == 0) {
            stim_name = strdup(key_value[1]);
            sdsfreesplitres(key_value, values_count);
            break;
        }

        sdsfreesplitres(key_value, values_count);
    }

    if(stim_name == NULL) {
        fprintf(stderr, "The stimulus name must be passed in the stimulus option! Exiting!\n");
        exit(EXIT_FAILURE);
    }

    struct config *sc = (struct config*) shget(stim_configs, stim_name);

    if(sc == NULL) {
        sc = alloc_and_init_config_data();
        print_to_stdout_and_file("Creating new stimulus name %s from command line options!\n", stim_name);
        shput(stim_configs, stim_name, sc);
    }

    for(int i = 0; i < tokens_count; i++) {

        int values_count;
        sds *key_value = sdssplit(extra_config_tokens[i], "=", &values_count);

        key_value[0] = sdstrim(key_value[0], " ");
        key_value[1] = sdstrim(key_value[1], " ");

        key = key_value[0];
        value = key_value[1];

        if(strcmp(key, "start") == 0) {

            bool stim_start_was_set;
            real stim_start = 0.0;
            GET_PARAMETER_NUMERIC_VALUE(real, stim_start, sc->config_data, key, stim_start_was_set);

            if(stim_start_was_set) {
                sprintf(old_value, "%lf", stim_start);
                issue_overwrite_warning(key, stim_name, old_value, value, config_file);
            }
            shput(sc->config_data, key, strdup(value));
        } else if(strcmp(key, "duration") == 0) {

            bool stim_duration_was_set;
            real stim_duration = 0.0;
            GET_PARAMETER_NUMERIC_VALUE(real, stim_duration, sc->config_data, key, stim_duration_was_set);

            if(stim_duration_was_set) {
                sprintf(old_value, "%lf", stim_duration);
                issue_overwrite_warning(key, stim_name, old_value, value, config_file);
            }
            shput_dup_value(sc->config_data, key, value);
        } else if(strcmp(key, "current") == 0) {

            bool stim_current_was_set;
            real stim_current = 0.0;
            GET_PARAMETER_NUMERIC_VALUE(real, stim_current, sc->config_data, key, stim_current_was_set);

            if(stim_current_was_set) {
                sprintf(old_value, "%lf", stim_current);
                issue_overwrite_warning(key, stim_name, old_value, value, config_file);
            }
            shput_dup_value(sc->config_data, key, value);

        } else if(strcmp(key, "period") == 0) {

            bool stim_period_was_set;
            real stim_period = 0.0;
            GET_PARAMETER_NUMERIC_VALUE(real, stim_period, sc->config_data, key, stim_period_was_set);

            if(stim_period_was_set) {
                sprintf(old_value, "%lf", stim_period);
                issue_overwrite_warning(key, stim_name, old_value, value, config_file);
            }
            shput_dup_value(sc->config_data, key, value);
        }
        else {
            set_or_overwrite_common_data(sc, key, value, stim_name, config_file);
        }
        sdsfreesplitres(key_value, values_count);
    }

    sdsfreesplitres(extra_config_tokens, tokens_count);
    free(stim_name);
}

void set_domain_config(const char *args, struct config *dc, const char *config_file) {

    sds extra_config;
    sds *extra_config_tokens;
    int tokens_count;
    extra_config = sdsnew(args);
    extra_config_tokens = sdssplit(extra_config, ",", &tokens_count);
    char old_value[32];
    char *key, *value;

    assert(dc);

    for(int i = 0; i < tokens_count; i++) {
        extra_config_tokens[i] = sdstrim(extra_config_tokens[i], " ");

        int values_count;
        sds *key_value = sdssplit(extra_config_tokens[i], "=", &values_count);

        if(values_count != 2) {
            fprintf(stderr, "Invalid format for options %s. Exiting!\n", args);
            exit(EXIT_FAILURE);
        }

        key_value[0] = sdstrim(key_value[0], " ");
        key_value[1] = sdstrim(key_value[1], " ");

        key = key_value[0];
        value = key_value[1];

        if(strcmp(key, "name") == 0) {
            char *domain_name = NULL;
            GET_PARAMETER_VALUE_CHAR_OR_USE_DEFAULT(domain_name, dc->config_data, "name");
            if(domain_name) {
                issue_overwrite_warning("name", "domain", domain_name, value, config_file);
            }
            shput(dc->config_data, key, strdup(value));
        } else if(strcmp(key, "start_dx") == 0) {

            bool start_dx_was_set;
            real_cpu start_dx = 0;
            GET_PARAMETER_NUMERIC_VALUE(real_cpu, start_dx, dc->config_data, key, start_dx_was_set);

            if(start_dx_was_set) {
                sprintf(old_value, "%lf", start_dx);
                issue_overwrite_warning("start_dx", "domain", old_value, value, config_file);

            }
            shput_dup_value(dc->config_data, key, value);

        } else if(strcmp(key, "maximum_dx") == 0) {

            bool max_dx_was_set;
            real_cpu max_dx = 0;
            GET_PARAMETER_NUMERIC_VALUE(real_cpu, max_dx, dc->config_data, key, max_dx_was_set);

            if(max_dx_was_set) {

                sprintf(old_value, "%lf", max_dx);
                issue_overwrite_warning("maximum_dx", "domain", old_value, value, config_file);
            }
            shput_dup_value(dc->config_data, key, value);
        }

        else if(strcmp(key, "start_dy") == 0) {

            bool start_dy_was_set;
            real_cpu start_dy = 0;
            GET_PARAMETER_NUMERIC_VALUE(real_cpu, start_dy, dc->config_data, key, start_dy_was_set);

            if(start_dy_was_set) {
                sprintf(old_value, "%lf", start_dy);
                issue_overwrite_warning("start_dy", "domain", old_value, value, config_file);
            }
            shput_dup_value(dc->config_data, key, value);
        } else if(strcmp(key, "maximum_dy") == 0) {

            bool max_dy_was_set;
            real_cpu max_dy = 0;
            GET_PARAMETER_NUMERIC_VALUE(real_cpu, max_dy, dc->config_data, key, max_dy_was_set);

            if(max_dy_was_set) {
                sprintf(old_value, "%lf", max_dy);
                issue_overwrite_warning("maximum_dy", "domain", old_value, value, config_file);
            }
            shput_dup_value(dc->config_data, key, value);
        }

        else if(strcmp(key, "start_dz") == 0) {

            bool start_dz_was_set;
            real_cpu start_dz = 0;
            GET_PARAMETER_NUMERIC_VALUE(real_cpu, start_dz, dc->config_data, key, start_dz_was_set);

            if(start_dz_was_set) {
                sprintf(old_value, "%lf", start_dz);
                issue_overwrite_warning("start_dz", "domain", old_value, value, config_file);
            }
            shput_dup_value(dc->config_data, key, value);

        } else if(strcmp(key, "maximum_dz") == 0) {

            bool max_dz_was_set;
            real_cpu max_dz = 0;
            GET_PARAMETER_NUMERIC_VALUE(real_cpu, max_dz, dc->config_data, "maximum_dz", max_dz_was_set);

            if(max_dz_was_set) {
                sprintf(old_value, "%lf", max_dz);
                issue_overwrite_warning("maximum_dz", "domain", old_value, value, config_file);
            }
            shput_dup_value(dc->config_data, key, value);
        }
        else {
            set_or_overwrite_common_data(dc, key, value, "domain", config_file);
        }
        sdsfreesplitres(key_value, values_count);
    }

    sdsfreesplitres(extra_config_tokens, tokens_count);
}

void set_save_mesh_config(const char *args, struct config *sm, const char *config_file) {

    sds extra_config;
    sds *extra_config_tokens;
    int tokens_count;
    extra_config = sdsnew(args);
    extra_config_tokens = sdssplit(extra_config, ",", &tokens_count);
    char old_value[32];
    char *key, *value;

    assert(sm);

    for(int i = 0; i < tokens_count; i++) {
        extra_config_tokens[i] = sdstrim(extra_config_tokens[i], " ");

        int values_count;
        sds *key_value = sdssplit(extra_config_tokens[i], "=", &values_count);

        if(values_count != 2) {
            fprintf(stderr, "Invalid format for options %s. Exiting!\n", args);
            exit(EXIT_FAILURE);
        }

        key_value[0] = sdstrim(key_value[0], " ");
        key_value[1] = sdstrim(key_value[1], " ");

        key = key_value[0];
        value = key_value[1];

        if(strcmp(key, "print_rate") == 0) {

            bool print_rate_was_set;
            int print_rate = 0;
            GET_PARAMETER_NUMERIC_VALUE(int, print_rate, sm->config_data, key, print_rate_was_set);

            if(print_rate_was_set) {
                sprintf(old_value, "%d", print_rate);
                issue_overwrite_warning(key, "save_mesh", old_value, value, config_file);
            }
            shput_dup_value(sm->config_data, key, value);

        } else if(strcmp(key, "output_dir") == 0) {

            char *out_dir_name = NULL;
            GET_PARAMETER_VALUE_CHAR_OR_USE_DEFAULT(out_dir_name, sm->config_data, key);

            if(out_dir_name) {
                issue_overwrite_warning(key, "save_mesh", out_dir_name, value, config_file);
            }

            shput_dup_value(sm->config_data, key, value);
        }
        else if(strcmp(key, "remove_older_simulation") == 0) {

            bool remove_older_simulation_was_set;
            bool remove_older_simulation = false;
            GET_PARAMETER_NUMERIC_VALUE(bool, remove_older_simulation, sm->config_data, key, remove_older_simulation_was_set);

            if(remove_older_simulation_was_set) {

                if(remove_older_simulation) {
                    sprintf(old_value, "yes");
                } else {
                    sprintf(old_value, "no");
                }

                issue_overwrite_warning(key, "save_mesh", old_value, optarg, config_file);
            }
            if(IS_TRUE(optarg) || IS_FALSE(optarg)) {
                shput_dup_value(sm->config_data, key, optarg);
            } else {
                fprintf(stderr,
                        "Warning: Invalid value for remove_older_simulation option: %s! Valid options are: true, yes, false, no, 0 or 1 "
                        "Setting the value to false\n",
                        optarg);
                shput_dup_value(sm->config_data, key, "false");
            }
        }
        else {
            set_or_overwrite_common_data(sm, key, value, "save_mesh", config_file);
        }
        sdsfreesplitres(key_value, values_count);
    }

    sdsfreesplitres(extra_config_tokens, tokens_count);
}

void set_config(const char *args, struct config *config, const char *config_file, char *config_type) {
    sds extra_config;
    sds *extra_config_tokens;
    int tokens_count;
    extra_config = sdsnew(args);
    extra_config_tokens = sdssplit(extra_config, ",", &tokens_count);
    char *opt_value;
    char *key, *value;

    struct string_hash_entry *sh = config->config_data;

    for(int i = 0; i < tokens_count; i++) {
        extra_config_tokens[i] = sdstrim(extra_config_tokens[i], " ");

        int values_count;
        sds *key_value = sdssplit(extra_config_tokens[i], "=", &values_count);

        if(values_count != 2) {
            fprintf(stderr, "Invalid format for options %s. Exiting!\n", args);
            exit(EXIT_FAILURE);
        }

        key_value[0] = sdstrim(key_value[0], " ");
        key_value[1] = sdstrim(key_value[1], " ");

        key = key_value[0];
        value = key_value[1];

        if(strcmp(key, "main_function") == 0) {
            if(config->main_function_name_was_set) {
                issue_overwrite_warning("main_function", config_type, config->main_function_name, value, config_file);
            }
            free(config->main_function_name);
            config->main_function_name = strdup(value);
        }
        if(strcmp(key, "init_function") == 0) {
            if(config->init_function_name_was_set) {
                issue_overwrite_warning("init_function", config_type, config->init_function_name, value, config_file);
            }
            free(config->init_function_name);
            config->init_function_name = strdup(value);
        }
        if(strcmp(key, "end_function") == 0) {
            if(config->end_function_name_was_set) {
                issue_overwrite_warning("end_function", config_type, config->end_function_name, value, config_file);
            }
            free(config->end_function_name);
            config->end_function_name = strdup(value);
        }
        else if(strcmp(key, "library_file") == 0) {
            if(config->library_file_path_was_set) {
                issue_overwrite_warning("library_file", config_type, config->library_file_path, value, config_file);
            }
            free(config->library_file_path);
            config->library_file_path = strdup(value);
        }
        else {
            opt_value = shget(sh, key);
            if(opt_value) {
                issue_overwrite_warning(key, config_type, opt_value, value, config_file);
            }

            shput(sh, key, strdup(value));
        }
        sdsfreesplitres(key_value, values_count);

        sdsfreesplitres(extra_config_tokens, tokens_count);
    }
}

void parse_batch_options(int argc, char **argv, struct batch_options *user_args) {

    int opt = 0;
    int option_index;

    opt = getopt_long_only(argc, argv, batch_opt_string, long_batch_options, &option_index);

    while(opt != -1) {
        switch(opt) {
            case 'c':
                user_args->batch_config_file = strdup(optarg);
                break;
            case 'h': /* fall-through is intentional */
            case '?':
                display_batch_usage(argv);
                break;
            default:
                /* You won't actually get here. */
                break;
        }

        opt = getopt_long(argc, argv, batch_opt_string, long_batch_options, &option_index);
    }
}

void parse_visualization_options(int argc, char **argv, struct visualization_options *user_args) {

    int opt = 0;
    int option_index;

    opt = getopt_long_only(argc, argv, visualization_opt_string, long_visualization_options, &option_index);

    while(opt != -1) {
        switch(opt) {
            case 'x':
                user_args->max_v = strtod(optarg, NULL);
                break;
            case 'm':
                user_args->min_v = strtod(optarg, NULL);
                break;
            case 'd':
                user_args->dt = strtod(optarg, NULL);
                break;
            case 'p':
                free(user_args->files_prefix);
                user_args->files_prefix = strdup(optarg);
                break;
            case 'v':
                user_args->pvd_file = strdup(optarg);
                break;
            case 'a':
                user_args->activation_map = strdup(optarg);
                break;
            case 'c':
                user_args->activation_map = strdup(optarg);
                user_args->save_activation_only = true;
                break;
            case 's':
                user_args->start_file = (int)strtod(optarg, NULL);
                break;
            case 'h': /* fall-through is intentional */
            case '?':
                display_visualization_usage(argv);
                break;
            default:
                break;
        }

        opt = getopt_long(argc, argv, visualization_opt_string, long_visualization_options, &option_index);
    }

    if(user_args->activation_map == NULL) {
        for (int index = optind; index < argc; index++)
            user_args->input_folder = strdup(argv[index]);
    }

}

void get_config_file(int argc, char **argv, struct user_options *user_args) {

    optind = 0;

    int opt = 0;

    int option_index;

    opt = getopt_long(argc, argv, opt_string, long_options, &option_index);

    while(opt != -1) {
        if(opt == 'c') {
            user_args->config_file = optarg;
            return;
        }
        opt = getopt_long(argc, argv, opt_string, long_options, &option_index);
    }

    // We reset the index after parsing the config_file
    optind = 1;
}

void parse_options(int argc, char **argv, struct user_options *user_args) {

    int opt = 0;
    int option_index;

    opt = getopt_long_only(argc, argv, opt_string, long_options, &option_index);
    char old_value[32];

    while(opt != -1) {
        switch(opt) {
            case 'R':
                if(user_args->refine_each_was_set) {
                    sprintf(old_value, "%d", user_args->refine_each);
                    issue_overwrite_warning("refine_each", "alg", old_value, optarg, user_args->config_file);
                }
                user_args->refine_each = (int)strtol(optarg, NULL, 10);

                break;
            case 'D':
                if(user_args->derefine_each_was_set) {
                    sprintf(old_value, "%d", user_args->derefine_each);
                    issue_overwrite_warning("derefine_each", "alg", old_value, optarg, user_args->config_file);
                }
                user_args->derefine_each = (int)strtol(optarg, NULL, 10);

                break;
            case 'e':
                if(user_args->dt_ode_was_set) {
                    sprintf(old_value, "%lf", user_args->dt_ode);
                    issue_overwrite_warning("dt_ode", "ode_solver", old_value, optarg, user_args->config_file);
                }
                user_args->dt_ode = strtof(optarg, NULL);
                break;
            case 'k':
                if(user_args->model_file_path_was_set) {
                    if(user_args->model_file_path) {
                        issue_overwrite_warning("model_file_path", "ode_solver", user_args->model_file_path, optarg,
                                                user_args->config_file);
                    } else {
                        issue_overwrite_warning("model_file_path", "ode_solver", "No Save", optarg, user_args->config_file);
                    }
                }
                free(user_args->model_file_path);
                user_args->model_file_path = strdup(optarg);

                break;
            case 'f':
                if(user_args->final_time_was_set) {
                    sprintf(old_value, "%lf", user_args->final_time);
                    issue_overwrite_warning("simulation_time", "main", old_value, optarg, user_args->config_file);
                }
                user_args->final_time = strtof(optarg, NULL);

                break;
            case 'r':
                if(user_args->ref_bound_was_set) {
                    sprintf(old_value, "%lf", user_args->ref_bound);
                    issue_overwrite_warning("refinement_bound", "alg", old_value, optarg, user_args->config_file);
                }
                user_args->ref_bound = strtof(optarg, NULL);
                break;
            case 'd':
                if(user_args->deref_bound_was_set) {
                    sprintf(old_value, "%lf", user_args->deref_bound);
                    issue_overwrite_warning("derefinement_bound", "alg", old_value, optarg, user_args->config_file);
                }
                user_args->deref_bound = strtof(optarg, NULL);
                break;
            case 'a':
                if(user_args->adaptive_was_set) {
                    sprintf(old_value, "%d", user_args->adaptive);
                    issue_overwrite_warning("use_adaptivity", "main", old_value, optarg, user_args->config_file);
                }
                user_args->adaptive = true;
                break;
            case 'n':
                if(((int)strtol(optarg, NULL, 10)) > 0) {
                    if(user_args->num_threads_was_set) {
                        sprintf(old_value, "%d", user_args->num_threads);
                        issue_overwrite_warning("n"
                                                "main",
                                                "num_threads",
                                                old_value, optarg, user_args->config_file);
                    }
                    user_args->num_threads = (int)strtol(optarg, NULL, 10);
                }
                break;
            case 'g':
                if(user_args->gpu_was_set) {

                    if(user_args->gpu) {
                        sprintf(old_value, "yes");
                    } else {
                        sprintf(old_value, "no");
                    }

                    issue_overwrite_warning("use_gpu", "ode_solver", old_value, optarg, user_args->config_file);
                }
                if(IS_TRUE(optarg)) {
                    user_args->gpu = true;
                } else if(IS_FALSE(optarg)) {
                    user_args->gpu = false;
                } else {
                    fprintf(stderr,
                            "Warning: Invalid value for use_gpu option: %s! Valid options are: true, yes, false, no, 0 or 1. "
                            "Setting the value to false\n",
                            optarg);
                    user_args->gpu = false;
                }
                break;
            case 'z':
                if(user_args->dt_pde_was_set) {
                    sprintf(old_value, "%lf", user_args->dt_pde);
                    issue_overwrite_warning("dt_pde", "main", old_value, optarg, user_args->config_file);
                }
                user_args->dt_pde = strtof(optarg, NULL);
                break;
            case BETA:
                if(user_args->beta_was_set) {
                    sprintf(old_value, "%lf", user_args->beta);
                    issue_overwrite_warning("beta", "main", old_value, optarg, user_args->config_file);
                }
                user_args->beta = strtof(optarg, NULL);
                break;
            case CM:
                if(user_args->cm) {
                    sprintf(old_value, "%lf", user_args->cm);
                    issue_overwrite_warning("cm", "main", old_value, optarg, user_args->config_file);
                }
                user_args->cm = strtof(optarg, NULL);
                break;
            case START_REFINING:
                if(user_args->start_adapting_at_was_set) {
                    sprintf(old_value, "%lf", user_args->start_adapting_at);
                    issue_overwrite_warning("start_adapting_at", "main", old_value, optarg, user_args->config_file);
                }
                user_args->start_adapting_at = strtof(optarg, NULL);
                break;
            case 'G':
                if(user_args->gpu_id_was_set) {
                    sprintf(old_value, "%d", user_args->gpu_id);
                    issue_overwrite_warning("gpu_id", "ode_solver", old_value, optarg, user_args->config_file);
                }
                user_args->gpu_id = (int)strtol(optarg, NULL, 10);
                break;
            case 'b':
                if(user_args->abort_no_activity_was_set) {
                    sprintf(old_value, "%d", user_args->abort_no_activity);
                    issue_overwrite_warning("abort_on_no_activity", "main", old_value, optarg, user_args->config_file);
                }
                user_args->abort_no_activity = true;
                break;
            case 'v':
                if(user_args->vm_threshold_was_set) {
                    sprintf(old_value, "%lf", user_args->vm_threshold);
                    issue_overwrite_warning("vm_threshold", "main", old_value, optarg, user_args->config_file);
                }
                user_args->vm_threshold = strtof(optarg, NULL);
                break;

            case DOMAIN_OPT:
                if(user_args->domain_config == NULL) {
                    print_to_stdout_and_file("Creating new domain config from command line!\n");
                    user_args->domain_config = alloc_and_init_config_data();
                }
                set_domain_config(optarg, user_args->domain_config, user_args->config_file);
                break;
            case SAVE_OPT:
                if(user_args->save_mesh_config == NULL) {
                    print_to_stdout_and_file("Creating new save config from command line!\n");
                    user_args->save_mesh_config = alloc_and_init_config_data();
                }
                set_save_mesh_config(optarg, user_args->save_mesh_config, user_args->config_file);
                break;
            case ASSEMBLY_MATRIX_OPT:
                if(user_args->assembly_matrix_config == NULL) {
                    print_to_stdout_and_file("Creating new assembly_matrix config from command line!\n");
                    user_args->assembly_matrix_config = alloc_and_init_config_data();
                }
                set_config(optarg, user_args->assembly_matrix_config, user_args->config_file, "assembly_matrix");
                break;
            case LINEAR_SYSTEM_SOLVER_OPT:
                if(user_args->linear_system_solver_config == NULL) {
                    print_to_stdout_and_file("Creating new linear_system_solver config from command line!\n");
                    user_args->linear_system_solver_config = alloc_and_init_config_data();
                }
                set_config(optarg, user_args->linear_system_solver_config, user_args->config_file, "linear_system_solver");
                break;
            case UPDATE_MONODOMAIN_SOLVER_OPT:
                if(user_args->update_monodomain_config == NULL) {
                    print_to_stdout_and_file("Creating new update_monodomain config from command line!\n");
                    user_args->update_monodomain_config = alloc_and_init_config_data();
                }
                set_config(optarg, user_args->update_monodomain_config, user_args->config_file, "update_monodomain");
                break;
            case EXTRA_DATA_OPT:
                if(user_args->extra_data_config == NULL) {
                    print_to_stdout_and_file("Creating new extra data config from command line!\n");
                    user_args->extra_data_config = alloc_and_init_config_data();
                }
                set_config(optarg, user_args->extra_data_config, user_args->config_file, "extra_data");
                break;
            case STIM_OPT:
                if(user_args->stim_configs == NULL) {
                    print_to_stdout_and_file("Creating new stim config from command line!\n");
                }
                set_stim_config(optarg, user_args->stim_configs, user_args->config_file);
                break;
            case SAVE_STATE_OPT:
                if(user_args->save_state_config == NULL) {
                    print_to_stdout_and_file("Creating new save state config from command line!\n");
                    user_args->save_state_config = alloc_and_init_config_data();
                }
                set_config(optarg, user_args->save_state_config, user_args->config_file, "save_state");
                break;
            case RESTORE_STATE_OPT:
                if(user_args->restore_state_config == NULL) {
                    print_to_stdout_and_file("Creating new restore state config from command line!\n");
                    user_args->restore_state_config = alloc_and_init_config_data();
                }
                set_config(optarg, user_args->restore_state_config, user_args->config_file, "restore_state");
                break;
            case DRAW_OPT:
                user_args->draw = true;
                break;
            case MAX_V_OPT:
                user_args->max_v = strtof(optarg, NULL);
                break;
            case MIN_V_OPT:
                user_args->min_v = strtof(optarg, NULL);
                break;
            case VISUALIZATION_PAUSED_OPT:
                user_args->start_visualization_unpaused = true;
                break;
            case 'q':
                if(user_args->quiet_was_set) {
                    sprintf(old_value, "%d", user_args->quiet);
                    issue_overwrite_warning("quiet", "main", old_value, optarg, user_args->config_file);
                }
                user_args->quiet = true;
                break;
            case 'h': /* fall-through is intentional */
            case '?':
                display_usage(argv);
                break;
            default:
                /* You won't actually get here. */
                break;
        }

        opt = getopt_long(argc, argv, opt_string, long_options, &option_index);
    }
}

int parse_batch_config_file(void *user, const char *section, const char *name, const char *value) {

    struct batch_options *pconfig = (struct batch_options *)user;

    if(MATCH_SECTION(BATCH_SECTION)) {

        if(MATCH_NAME("output_folder")) {
            pconfig->output_folder = strdup(value);

        } else if(MATCH_NAME("num_simulations_per_parameter_change")) {
            pconfig->num_simulations = (int)strtol(value, NULL, 10);
        } else if(MATCH_NAME("initial_config")) {
            pconfig->initial_config = strdup(value);
        } else if(MATCH_NAME("skip_existing_run")) {

            if(IS_TRUE(value)) {
                pconfig->skip_existing_run = true;
            } else if(IS_FALSE(value)) {
                pconfig->skip_existing_run = false;
            } else {
                fprintf(stderr,
                        "Warning: Invalid value for skip_existing_run option: %s! Valid options are: true, yes, false, no, 0 or 1. "
                        "Setting the value to false\n",
                        value);
                pconfig->skip_existing_run = false;
            }
        }
    } else if(MATCH_SECTION(MODIFICATION_SECTION)) {
        shput(pconfig->config_to_change, name, strdup(value));
    } else {
        fprintf(stderr, "\033[33;5;7mInvalid name %s in section %s on the batch config file!\033[0m\n", name, section);
        return 0;
    }

    return 1;
}

int parse_config_file(void *user, const char *section, const char *name, const char *value) {
    struct user_options *pconfig = (struct user_options *)user;

    if(MATCH_SECTION(MAIN_SECTION)) {
        pconfig->main_found = true;
    }

    if(MATCH_SECTION_AND_NAME(MAIN_SECTION, "num_threads")) {
        pconfig->num_threads = (int)strtol(value, NULL, 10);
        pconfig->num_threads_was_set = true;
    } else if(MATCH_SECTION_AND_NAME(MAIN_SECTION, "dt_pde")) {
        pconfig->dt_pde = strtof(value, NULL);
        pconfig->dt_pde_was_set = true;
    } else if(MATCH_SECTION_AND_NAME(MAIN_SECTION, "simulation_time")) {
        pconfig->final_time = strtof(value, NULL);
        pconfig->final_time_was_set = true;
    } else if(MATCH_SECTION_AND_NAME(MAIN_SECTION, "beta")) {
        pconfig->beta = strtof(value, NULL);
        pconfig->beta_was_set = true;
    } else if(MATCH_SECTION_AND_NAME(MAIN_SECTION, "cm")) {
        pconfig->cm = strtof(value, NULL);
        pconfig->cm_was_set = true;
    } else if(MATCH_SECTION_AND_NAME(MAIN_SECTION, "start_adapting_at")) {
        pconfig->start_adapting_at = strtof(value, NULL);
        pconfig->start_adapting_at_was_set = true;
    } else if(MATCH_SECTION_AND_NAME(MAIN_SECTION, "abort_on_no_activity")) {
        if(IS_TRUE(value)) {
            pconfig->abort_no_activity = true;
        } else {
            pconfig->abort_no_activity = false;
        }
        pconfig->abort_no_activity_was_set = true;
    } else if(MATCH_SECTION_AND_NAME(MAIN_SECTION, "calc_activation_time")) {
        if(strcmp(value, "true") == 0 || strcmp(value, "yes") == 0) {
            pconfig->calc_activation_time = true;
        } else {
            pconfig->calc_activation_time = false;
        }
        pconfig->calc_activation_time_was_set = true;
    } 
    else if(MATCH_SECTION_AND_NAME(MAIN_SECTION, "print_conductivity_map")) {
        if(strcmp(value, "true") == 0 || strcmp(value, "yes") == 0) {
            pconfig->print_conductivity_map = true;
        } else {
            pconfig->print_conductivity_map = false;
        }
        pconfig->print_conductivity_map_was_set = true;
    } 
    else if(MATCH_SECTION_AND_NAME(MAIN_SECTION, "print_min_vm_map")) {
        if(strcmp(value, "true") == 0 || strcmp(value, "yes") == 0) {
            pconfig->print_min_vm_map = true;
        } else {
            pconfig->print_min_vm_map = false;
        }
        pconfig->print_min_vm_map_was_set = true;
    }
    else if(MATCH_SECTION_AND_NAME(MAIN_SECTION, "print_max_vm_map")) {
        if(strcmp(value, "true") == 0 || strcmp(value, "yes") == 0) {
            pconfig->print_max_vm_map = true;
        } else {
            pconfig->print_max_vm_map = false;
        }
        pconfig->print_max_vm_map_was_set = true;
    }
    else if(MATCH_SECTION_AND_NAME(MAIN_SECTION, "print_apd_map")) {
        if(strcmp(value, "true") == 0 || strcmp(value, "yes") == 0) {
            pconfig->print_apd_map = true;
        } else {
            pconfig->print_apd_map = false;
        }
        pconfig->print_apd_map_was_set = true;
    }
    else if(MATCH_SECTION_AND_NAME(MAIN_SECTION, "vm_threshold")) {
        pconfig->vm_threshold = strtof(value, NULL);
        pconfig->vm_threshold_was_set = true;
    } else if(MATCH_SECTION_AND_NAME(MAIN_SECTION, "use_adaptivity")) {
        if( IS_TRUE(value) ) {
            pconfig->adaptive = true;
        } else {
            pconfig->adaptive = false;
        }
        pconfig->adaptive_was_set = true;
    } else if(MATCH_SECTION_AND_NAME(MAIN_SECTION, "quiet")) {
        if(IS_TRUE(value)) {
            pconfig->quiet = true;
        } else {
            pconfig->quiet = false;
        }
        pconfig->quiet = true;
    } else if(MATCH_SECTION_AND_NAME(ALG_SECTION, "refinement_bound")) {
        pconfig->ref_bound = strtof(value, NULL);
        pconfig->ref_bound_was_set = true;
    } else if(MATCH_SECTION_AND_NAME(ALG_SECTION, "derefinement_bound")) {
        pconfig->deref_bound = strtof(value, NULL);
        pconfig->deref_bound_was_set = true;
    } else if(MATCH_SECTION_AND_NAME(ALG_SECTION, "refine_each")) {
        pconfig->refine_each = (int)strtol(value, NULL, 10);
        pconfig->refine_each_was_set = true;
    } else if(MATCH_SECTION_AND_NAME(ALG_SECTION, "derefine_each")) {
        pconfig->derefine_each = (int)strtol(value, NULL, 10);
        pconfig->derefine_each_was_set = true;
    } else if(MATCH_SECTION(ODE_SECTION)) {

        if(MATCH_NAME("dt_ode")) {
            pconfig->dt_ode = strtof(value, NULL);
            pconfig->dt_ode_was_set = true;
        }
        else if(MATCH_NAME("use_gpu")) {
            if(IS_TRUE(value)) {
                pconfig->gpu = true;
            } else {
                pconfig->gpu = false;
            }
            pconfig->gpu_was_set = true;
        } else if(MATCH_NAME("gpu_id")) {
            pconfig->gpu_id = (int)strtol(value, NULL, 10);
            pconfig->gpu_id_was_set = true;
        } else if(MATCH_NAME("library_file")) {
            pconfig->model_file_path = strdup(value);
            pconfig->model_file_path_was_set = true;
        }
        else {
            shput(pconfig->ode_extra_config, name, strdup(value));
        }
    } else if(MATCH_SECTION(ODE_PURKINJE_SECTION)) {

        if(MATCH_NAME("dt_ode")) {
            pconfig->purkinje_dt_ode = strtof(value, NULL);
            pconfig->purkinje_dt_ode_was_set = true;
        }
        else if(MATCH_NAME("use_gpu")) {
            if(IS_TRUE(value)) {
                pconfig->purkinje_gpu = true;
            } else {
                pconfig->purkinje_gpu = false;
            }
            pconfig->purkinje_gpu_was_set = true;
        } else if(MATCH_NAME("gpu_id")) {
            pconfig->purkinje_gpu_id = (int)strtol(value, NULL, 10);
            pconfig->purkinje_gpu_id_was_set = true;
        } else if(MATCH_NAME("library_file")) {
            pconfig->purkinje_model_file_path = strdup(value);
            pconfig->purkinje_model_file_path_was_set = true;
        }
        else {
            shput(pconfig->purkinje_ode_extra_config, name, strdup(value));
        }
    } else if(SECTION_STARTS_WITH(STIM_PURKINJE_SECTION)) {

        struct config *tmp = (struct config *) shget(pconfig->purkinje_stim_configs, section);

        if (tmp == NULL) {
            tmp = alloc_and_init_config_data();
            shput(pconfig->purkinje_stim_configs, section, tmp);
        }

        if (MATCH_NAME("name")) {
            fprintf(stderr,
                    "name is a reserved word and should not be used inside a stimulus config section. Found in %s. "
                    "Exiting!\n",
                    section);
            exit(EXIT_FAILURE);
        }
        else {
            set_common_data(tmp, name, value);
        }

    } else if(SECTION_STARTS_WITH(STIM_SECTION)) {

        struct config *tmp = (struct config *) shget(pconfig->stim_configs, section);

        if (tmp == NULL) {
            tmp = alloc_and_init_config_data();
            shput(pconfig->stim_configs, section, tmp);
        }

        if (MATCH_NAME("name")) {
            fprintf(stderr,
                    "name is a reserved word and should not be used inside a stimulus config section. Found in %s. "
                    "Exiting!\n",
                    section);
            exit(EXIT_FAILURE);
        }
        else {
            set_common_data(tmp, name, value);
        }

    } else if(MATCH_SECTION(DOMAIN_SECTION)) {

        if(pconfig->domain_config == NULL) {
            pconfig->domain_config = alloc_and_init_config_data();
        }

        set_common_data(pconfig->domain_config, name, value);
    }
    else if(MATCH_SECTION(PURKINJE_SECTION)) {

        if(pconfig->purkinje_config == NULL) {
            pconfig->purkinje_config = alloc_and_init_config_data();
        }

        set_common_data(pconfig->purkinje_config, name, value);

    } else if(MATCH_SECTION(MATRIX_ASSEMBLY_SECTION)) {

        if(pconfig->assembly_matrix_config == NULL) {
            pconfig->assembly_matrix_config = alloc_and_init_config_data();
        }

        set_common_data(pconfig->assembly_matrix_config, name, value);

    }
    else if(MATCH_SECTION(UPDATE_MONODOMAIN_SECTION)) {

        if (pconfig->update_monodomain_config == NULL) {
            pconfig->update_monodomain_config = alloc_and_init_config_data();
        }

        set_common_data(pconfig->update_monodomain_config, name, value);
    }
    else if(MATCH_SECTION(LINEAR_SYSTEM_SOLVER_SECTION)) {

        if(pconfig->linear_system_solver_config == NULL) {
            pconfig->linear_system_solver_config = alloc_and_init_config_data();
        }

        set_common_data(pconfig->linear_system_solver_config, name, value);

    } else if(MATCH_SECTION(EXTRA_DATA_SECTION)) {

        if(pconfig->extra_data_config == NULL) {
            pconfig->extra_data_config = alloc_and_init_config_data();
        }

        set_common_data(pconfig->extra_data_config, name, value);

    } else if(MATCH_SECTION(SAVE_RESULT_SECTION)) {

        if(pconfig->save_mesh_config == NULL) {
            pconfig->save_mesh_config = alloc_and_init_config_data();
        }

        set_common_data(pconfig->save_mesh_config, name, value);

    } else if(MATCH_SECTION(SAVE_STATE_SECTION)) {

        if(pconfig->save_state_config == NULL) {
            pconfig->save_state_config = alloc_and_init_config_data();
        }

        set_common_data(pconfig->save_state_config, name, value);

    } else if(MATCH_SECTION(RESTORE_STATE_SECTION)) {

        if(pconfig->restore_state_config == NULL) {
            pconfig->restore_state_config = alloc_and_init_config_data();
        }

        set_common_data(pconfig->restore_state_config, name, value);

    } else {

        fprintf(stderr, "\033[33;5;7mInvalid name %s in section %s on the config file!\033[0m\n", name, section);
        return 0;
    }

    return 1;
}

#define WRITE_INI_SECTION(SECTION) fprintf(ini_file, "[%s]\n", SECTION)

#define WRITE_NAME_VALUE_WITHOUT_CHECK(NAME, VALUE, SPECIFIER) fprintf(ini_file, "%s = %" SPECIFIER "\n", NAME, VALUE)

#define WRITE_NAME_VALUE(NAME, VALUE, SPECIFIER) \
    do { \
        if(VALUE ## _was_set) \
            WRITE_NAME_VALUE_WITHOUT_CHECK(NAME, VALUE, SPECIFIER); \
    }while(0)

#define WRITE_EXTRA_CONFIG(hash) \
    do {\
        for (long j = 0; j < hmlen(hash); j++) { \
            char *name = hash[j].key; \
            if (strcmp(name, "main_function") != 0 && strcmp(name, "library_file") != 0) { \
                char *value = hash[j].value;\
                WRITE_NAME_VALUE_WITHOUT_CHECK(name, value, "s");\
            }\
        }\
    }\
    while(0)

void write_ini_options(struct config *config, FILE *ini_file) {
    WRITE_NAME_VALUE("main_function", config->main_function_name, "s");
    WRITE_NAME_VALUE("init_function", config->init_function_name, "s");
    WRITE_NAME_VALUE("end_function", config->end_function_name, "s");
    WRITE_NAME_VALUE("library_file", config->library_file_path, "s");
    WRITE_EXTRA_CONFIG(config->config_data);
}

void options_to_ini_file(struct user_options *config, char *ini_file_path) {

    FILE *ini_file = fopen(ini_file_path, "w");

    if(!ini_file) {
        fprintf(stderr, "options_to_ini_file. Error creating %s\n",ini_file_path);
        return;
    }

    WRITE_INI_SECTION(MAIN_SECTION);
    WRITE_NAME_VALUE("num_threads", config->num_threads, "d");
    WRITE_NAME_VALUE("dt_pde" , config->dt_pde, "f");
    WRITE_NAME_VALUE("simulation_time", config->final_time, "f");
    WRITE_NAME_VALUE("beta", config->beta, "f");
    WRITE_NAME_VALUE("cm", config->cm, "f");
    WRITE_NAME_VALUE("start_adapting_at", config->start_adapting_at, "f");
    WRITE_NAME_VALUE("abort_on_no_activity", config->abort_no_activity, "d");
    WRITE_NAME_VALUE("vm_threshold", config->vm_threshold, "f");
    WRITE_NAME_VALUE("use_adaptivity", config->adaptive, "d");
    WRITE_NAME_VALUE("quiet", config->quiet, "d");
    WRITE_NAME_VALUE("refinement_bound", config->ref_bound, "f");
    WRITE_NAME_VALUE("derefinement_bound", config->deref_bound, "f");
    WRITE_NAME_VALUE("refine_each", config->refine_each, "d");
    WRITE_NAME_VALUE("derefine_each", config->derefine_each, "d");
    printf("\n");


    WRITE_INI_SECTION(ODE_SECTION);
    WRITE_NAME_VALUE("dt_ode", config->dt_ode, "f");
    WRITE_NAME_VALUE("use_gpu", config->gpu, "d");
    WRITE_NAME_VALUE("gpu_id", config->gpu_id, "d");
    WRITE_NAME_VALUE("library_file", config->model_file_path, "s");
    printf("\n");

    for(long i = 0; i < hmlen(config->stim_configs); i++) {
        struct string_voidp_hash_entry e = config->stim_configs[i];
        WRITE_INI_SECTION(e.key);
        struct config *tmp = (struct config*) e.value;
        write_ini_options(tmp, ini_file);
        printf("\n");
    }

    if(config->domain_config) {
        WRITE_INI_SECTION(DOMAIN_SECTION);
        write_ini_options(config->domain_config, ini_file);
        printf("\n");
    }

    if(config->purkinje_config) {
        WRITE_INI_SECTION(PURKINJE_SECTION);
        write_ini_options(config->purkinje_config, ini_file);
        printf("\n");

    }

    if(config->assembly_matrix_config) {
        WRITE_INI_SECTION(MATRIX_ASSEMBLY_SECTION);
        write_ini_options(config->assembly_matrix_config, ini_file);
        printf("\n");
    }

    if(config->update_monodomain_config) {
        WRITE_INI_SECTION(UPDATE_MONODOMAIN_SECTION);
        write_ini_options(config->update_monodomain_config, ini_file);
        printf("\n");
    }

    if(config->linear_system_solver_config) {
        WRITE_INI_SECTION(LINEAR_SYSTEM_SOLVER_SECTION);
        write_ini_options(config->linear_system_solver_config, ini_file);
        printf("\n");
    }

    if(config->extra_data_config) {
        WRITE_INI_SECTION(EXTRA_DATA_SECTION);
        write_ini_options(config->extra_data_config, ini_file);
        printf("\n");
    }

    if(config->save_mesh_config) {
        WRITE_INI_SECTION(SAVE_RESULT_SECTION);
        write_ini_options(config->save_mesh_config, ini_file);
        printf("\n");
    }

    if(config->save_state_config) {
        WRITE_INI_SECTION(SAVE_STATE_SECTION);
        write_ini_options(config->save_state_config, ini_file);
        printf("\n");
    }

    if(config->restore_state_config) {
        WRITE_INI_SECTION(RESTORE_STATE_SECTION);
        write_ini_options(config->restore_state_config, ini_file);
        printf("\n");
    }

    fclose(ini_file);

}
#undef WRITE_INI_SECTION
#undef WRITE_NAME_VALUE
#undef WRITE_NAME_VALUE_WITHOUT_CHECK
#undef WRITE_EXTRA_CONFIG

void configure_grid_from_options(struct grid *grid, struct user_options *options) {

    assert(grid);
    assert(options);

    grid->adaptive = options->adaptive;
}

void free_user_options(struct user_options *s) {

    free(s->model_file_path);

    if(s->stim_configs) {
        STIM_CONFIG_HASH_FOR_EACH_KEY_APPLY_FN_IN_VALUE(s->stim_configs, free_config_data);
                shfree(s->stim_configs);
    }

    if(s->extra_data_config)
        free_config_data(s->extra_data_config);

    if(s->domain_config)
        free_config_data(s->domain_config);

    if(s->save_mesh_config)
        free_config_data(s->save_mesh_config);

    if(s->assembly_matrix_config)
        free_config_data(s->assembly_matrix_config);

    if(s->linear_system_solver_config)
        free_config_data(s->linear_system_solver_config);

    if(s->restore_state_config)
        free_config_data(s->restore_state_config);

    if(s->save_state_config)
        free_config_data(s->save_state_config);

    if(s->update_monodomain_config)
        free_config_data(s->update_monodomain_config);

    free(s);
}
