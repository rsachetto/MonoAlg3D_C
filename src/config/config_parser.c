#include <assert.h>
#include <string.h>

#include "../ini_parser/ini_file_sections.h"
#include "../string/sds.h"
#include "../utils/file_utils.h"
#include "config_parser.h"
#include "stim_config.h"

#include "../common_types/common_types.h"
#include "../single_file_libraries/stb_ds.h"

static const struct option long_batch_options[] = {{"config_file", required_argument, NULL, 'c'}};

static const char *batch_opt_string = "c:h";

//TODO: we need to document the complex options. See comments below
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

static const char *opt_string = "c:abn:g:m:t:r:d:z:e:f:jR:D:G:k:v:qh";

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
    printf("--use_gpu | -g [yes|no|true|false]. Solve ODEs using GPU. Default: No \n");
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

void issue_overwrite_warning(const char *var, const char *section, const char *old_value, const char *new_value, const char *config_file) {
    fprintf(stderr,
            "WARNING: option %s in %s was set in the file %s to %s and is being overwritten "
            "by the command line flag to %s!\n",
            var, section, config_file, old_value, new_value);
}

struct batch_options *new_batch_options() {
    struct batch_options *user_args = (struct batch_options *)malloc(sizeof(struct batch_options));
    user_args->batch_config_file = NULL;
    user_args->initial_config = NULL;
    user_args->config_to_change = NULL;
    user_args->num_simulations = 0;

    return user_args;
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

    user_args->stim_configs = NULL;

    sh_new_arena(user_args->stim_configs);
    shdefault(user_args->stim_configs, NULL);

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

void set_stim_config(const char *args, struct string_voidp_hash_entry *stim_configs, const char *config_file) {

    sds extra_config;
    sds *extra_config_tokens;
    int tokens_count;
    extra_config = sdsnew(args);
    extra_config_tokens = sdssplit(extra_config, ",", &tokens_count);
    char *opt_value;
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

    struct stim_config *sc = (struct stim_config*) shget(stim_configs, stim_name);

    if(sc == NULL) {
        sc = new_stim_config();
        print_to_stdout_and_file("Creating new stimulus name %s from command line options!\n", stim_name);
        shput(stim_configs, stim_name, sc);
    }

    struct string_hash_entry *sh = sc->config_data.config;

    for(int i = 0; i < tokens_count; i++) {

        int values_count;
        sds *key_value = sdssplit(extra_config_tokens[i], "=", &values_count);

        key_value[0] = sdstrim(key_value[0], " ");
        key_value[1] = sdstrim(key_value[1], " ");

        key = key_value[0];
        value = key_value[1];

        if(strcmp(key, "start") == 0) {
            if(sc->stim_start_was_set) {
                sprintf(old_value, "%lf", sc->stim_start);
                issue_overwrite_warning("start", stim_name, old_value, value, config_file);
            }
            sc->stim_start = (real)strtod(value, NULL);
        } else if(strcmp(key, "duration") == 0) {
            if(sc->stim_duration_was_set) {
                sprintf(old_value, "%lf", sc->stim_duration);
                issue_overwrite_warning("duration", stim_name, old_value, value, config_file);
            }

            sc->stim_duration = (real)strtod(value, NULL);
        } else if(strcmp(key, "current") == 0) {
            if(sc->stim_current_was_set) {
                sprintf(old_value, "%lf", sc->stim_current);
                issue_overwrite_warning("current", stim_name, old_value, value, config_file);
            }
            sc->stim_current = (real)strtod(value, NULL);
        } else if(strcmp(key, "function") == 0) {
            if(sc->config_data.function_name_was_set) {
                issue_overwrite_warning("function", stim_name, sc->config_data.function_name, value, config_file);
            }
            free(sc->config_data.function_name);
            sc->config_data.function_name = strdup(value);
        } else if(strcmp(key, "period") == 0) {
            if(sc->stim_period_was_set) {
                sprintf(old_value, "%lf", sc->stim_period);
                issue_overwrite_warning("period", stim_name, old_value, value, config_file);
            }
            sc->stim_period = (real)strtod(value, NULL);
        } else if(strcmp(key, "library_file") == 0) {
            if(sc->config_data.library_file_path_was_set) {
                issue_overwrite_warning("library_file", stim_name, sc->config_data.library_file_path, value, config_file);
            }
            free(sc->config_data.library_file_path);
            sc->config_data.library_file_path = strdup(value);
        } else {
            opt_value = shget(sh, key);
            if(opt_value) {
                issue_overwrite_warning(key, stim_name, opt_value, value, config_file);
            }
            shput(sh, key, strdup(value));
        }
        sdsfreesplitres(key_value, values_count);
    }

    sdsfreesplitres(extra_config_tokens, tokens_count);
    free(stim_name);
}

void set_domain_config(const char *args, struct domain_config *dc, const char *config_file) {

    sds extra_config;
    sds *extra_config_tokens;
    int tokens_count;
    extra_config = sdsnew(args);
    extra_config_tokens = sdssplit(extra_config, ",", &tokens_count);
    char *opt_value;
    char old_value[32];
    char *key, *value;

    assert(dc);

    struct string_hash_entry *sh = dc->config_data.config;

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
            if(dc->domain_name_was_set) {
                issue_overwrite_warning("name", "domain", dc->domain_name, value, config_file);
            }
            free(dc->domain_name);
            dc->domain_name = strdup(value);
        } else if(strcmp(key, "start_dx") == 0) {
            if(dc->start_dx_was_set) {
                sprintf(old_value, "%lf", dc->start_dx);
                issue_overwrite_warning("start_dx", "domain", old_value, value, config_file);
            }
            dc->start_dx = (real)strtod(value, NULL);
        } else if(strcmp(key, "maximum_dx") == 0) {
            if(dc->max_dx_was_set) {
                sprintf(old_value, "%lf", dc->max_dx);
                issue_overwrite_warning("maximum_dx", "domain", old_value, value, config_file);
            }
            dc->max_dx = strtod(value, NULL);
        }

        else if(strcmp(key, "start_dy") == 0) {
            if(dc->start_dy_was_set) {
                sprintf(old_value, "%lf", dc->start_dy);
                issue_overwrite_warning("start_dy", "domain", old_value, value, config_file);
            }
            dc->start_dy = strtod(value, NULL);
        } else if(strcmp(key, "maximum_dy") == 0) {
            if(dc->max_dy_was_set) {
                sprintf(old_value, "%lf", dc->max_dy);
                issue_overwrite_warning("maximum_dy", "domain", old_value, value, config_file);
            }
            dc->max_dy = (real)strtod(value, NULL);
        }

        else if(strcmp(key, "start_dz") == 0) {
            if(dc->start_dz_was_set) {
                sprintf(old_value, "%lf", dc->start_dz);
                issue_overwrite_warning("start_dz", "domain", old_value, value, config_file);
            }
            dc->start_dz = (real)strtod(value, NULL);
        } else if(strcmp(key, "maximum_dz") == 0) {
            if(dc->max_dz_was_set) {
                sprintf(old_value, "%lf", dc->max_dz);
                issue_overwrite_warning("maximum_dz", "domain", old_value, value, config_file);
            }
            dc->max_dz = (real)strtod(value, NULL);
        }
        else if(strcmp(key, "function") == 0) {
            if(dc->config_data.function_name_was_set) {
                issue_overwrite_warning("function", "domain", dc->config_data.function_name, value, config_file);
            }
            free(dc->config_data.function_name);
            dc->config_data.function_name = strdup(value);
        } else if(strcmp(key, "library_file") == 0) {
            if(dc->config_data.library_file_path_was_set) {
                issue_overwrite_warning("library_file", "domain", dc->config_data.library_file_path, value, config_file);
            }
            free(dc->config_data.library_file_path);
            dc->config_data.library_file_path = strdup(value);
        } else {
            opt_value = shget(sh, key);
            if(opt_value) {
                issue_overwrite_warning(key, "domain", opt_value, value, config_file);
            }

            shput(sh, key, strdup(value));
        }
        sdsfreesplitres(key_value, values_count);
    }

    sdsfreesplitres(extra_config_tokens, tokens_count);
}

void set_save_mesh_config(const char *args, struct save_mesh_config *sm, const char *config_file) {

    sds extra_config;
    sds *extra_config_tokens;
    int tokens_count;
    extra_config = sdsnew(args);
    extra_config_tokens = sdssplit(extra_config, ",", &tokens_count);
    char *opt_value;
    char old_value[32];
    char *key, *value;

    assert(sm);

    struct string_hash_entry *sh = sm->config_data.config;

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
            if(sm->print_rate_was_set) {
                sprintf(old_value, "%d", sm->print_rate);
                issue_overwrite_warning("print_rate", "save_mesh", old_value, value, config_file);
            }
            sm->print_rate = (int)strtol(value, NULL, 10);
        } else if(strcmp(key, "output_dir") == 0) {
            if(sm->out_dir_name_was_set) {
                issue_overwrite_warning("output_dir", "save_mesh", sm->out_dir_name, value, config_file);
            }
            free(sm->out_dir_name);
            sm->out_dir_name = strdup(value);
        } else if(strcmp(key, "function") == 0) {
            if(sm->config_data.function_name_was_set) {
                issue_overwrite_warning("function", "save_mesh", sm->config_data.function_name, value, config_file);
            }
            free(sm->config_data.function_name);
            sm->config_data.function_name = strdup(value);
        } else if(strcmp(key, "library_file") == 0) {
            if(sm->config_data.library_file_path_was_set) {
                issue_overwrite_warning("library_file", "save_mesh", sm->config_data.library_file_path, value, config_file);
            }
            free(sm->config_data.library_file_path);
            sm->config_data.library_file_path = strdup(value);
        } else {
            opt_value = shget(sh, key);
            if(opt_value) {
                issue_overwrite_warning(key, "save_mesh", opt_value, value, config_file);
            }

            shput(sh, key, strdup(value));
        }
        sdsfreesplitres(key_value, values_count);
    }

    sdsfreesplitres(extra_config_tokens, tokens_count);
}

void set_config(const char *args, void *some_config, const char *config_file, char *config_type) {
    sds extra_config;
    sds *extra_config_tokens;
    int tokens_count;
    extra_config = sdsnew(args);
    extra_config_tokens = sdssplit(extra_config, ",", &tokens_count);
    char *opt_value;
    char *key, *value;

    assert(some_config);

    struct generic_config *config = (struct generic_config *)some_config;

    struct string_hash_entry *sh = config->config_data.config;

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

        if(strcmp(key, "function") == 0) {
            if(config->config_data.function_name_was_set) {
                issue_overwrite_warning("function", config_type, config->config_data.function_name, value, config_file);
            }
            free(config->config_data.function_name);
            config->config_data.function_name = strdup(value);
        } else if(strcmp(key, "library_file") == 0) {
            if(config->config_data.library_file_path_was_set) {
                issue_overwrite_warning("library_file", config_type, config->config_data.library_file_path, value, config_file);
            }
            free(config->config_data.library_file_path);
            config->config_data.library_file_path = strdup(value);
        } else {

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

void get_config_file(int argc, char **argv, struct user_options *user_args) {

    optind = 0;

    int opt = 0;

    int option_index;

    opt = getopt_long(argc, argv, opt_string, long_options, &option_index);

    while(opt != -1) {
        switch(opt) {
        case 'c':
            user_args->config_file = optarg;
            return;
        default:
            break;
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
            if(strcmp(optarg, "true") == 0 || strcmp(optarg, "yes") == 0) {
                user_args->gpu = true;
            } else if(strcmp(optarg, "false") == 0 || strcmp(optarg, "no") == 0) {
                user_args->gpu = false;
            } else {
                fprintf(stderr,
                        "Warning: Invalid value for use_gpu option: %s! Valid options are: true, yes, false, no. "
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
                user_args->domain_config = new_domain_config();
            }
            set_domain_config(optarg, user_args->domain_config, user_args->config_file);
            break;
        case SAVE_OPT:
            if(user_args->save_mesh_config == NULL) {
                print_to_stdout_and_file("Creating new save config from command line!\n");
                user_args->save_mesh_config = new_save_mesh_config();
            }
            set_save_mesh_config(optarg, user_args->save_mesh_config, user_args->config_file);
            break;
        case ASSEMBLY_MATRIX_OPT:
            if(user_args->assembly_matrix_config == NULL) {
                print_to_stdout_and_file("Creating new assembly_matrix config from command line!\n");
                user_args->assembly_matrix_config = new_assembly_matrix_config();
            }
            set_config(optarg, user_args->assembly_matrix_config, user_args->config_file, "assembly_matrix");
            break;
        case LINEAR_SYSTEM_SOLVER_OPT:
            if(user_args->linear_system_solver_config == NULL) {
                print_to_stdout_and_file("Creating new linear_system_solver config from command line!\n");
                user_args->linear_system_solver_config = new_linear_system_solver_config();
            }
            set_config(optarg, user_args->linear_system_solver_config, user_args->config_file, "linear_system_solver");
            break;
        case UPDATE_MONODOMAIN_SOLVER_OPT:
            if(user_args->update_monodomain_config == NULL) {
                print_to_stdout_and_file("Creating new update_monodomain config from command line!\n");
                user_args->update_monodomain_config = new_update_monodomain_config();
            }
            set_config(optarg, user_args->update_monodomain_config, user_args->config_file, "update_monodomain");
            break;
        case EXTRA_DATA_OPT:
            if(user_args->extra_data_config == NULL) {
                print_to_stdout_and_file("Creating new extra data config from command line!\n");
                user_args->extra_data_config = new_extra_data_config();
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
                user_args->save_state_config = new_save_state_config();
            }
            set_config(optarg, user_args->save_state_config, user_args->config_file, "save_state");
            break;
        case RESTORE_STATE_OPT:
            if(user_args->restore_state_config == NULL) {
                print_to_stdout_and_file("Creating new restore state config from command line!\n");
                user_args->restore_state_config = new_restore_state_config();
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
        } else if(MATCH_NAME("num_parameter_change")) {
            pconfig->num_par_change = (int)strtol(value, NULL, 10);
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
        if(strcmp(value, "true") == 0 || strcmp(value, "yes") == 0) {
            pconfig->abort_no_activity = true;
        } else {
            pconfig->abort_no_activity = false;
        }
        pconfig->abort_no_activity_was_set = true;
    } else if(MATCH_SECTION_AND_NAME(MAIN_SECTION, "vm_threshold")) {
        pconfig->vm_threshold = strtof(value, NULL);
        pconfig->vm_threshold_was_set = true;
    } else if(MATCH_SECTION_AND_NAME(MAIN_SECTION, "use_adaptivity")) {
        if(strcmp(value, "true") == 0 || strcmp(value, "yes") == 0) {
            pconfig->adaptive = true;
        } else {
            pconfig->adaptive = false;
        }
        pconfig->adaptive_was_set = true;
    } else if(MATCH_SECTION_AND_NAME(MAIN_SECTION, "quiet")) {
        if(strcmp(value, "true") == 0 || strcmp(value, "yes") == 0) {
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
    } else if(MATCH_SECTION_AND_NAME(ODE_SECTION, "dt_ode")) {
        pconfig->dt_ode = strtof(value, NULL);
        pconfig->dt_ode_was_set = true;
    } else if(MATCH_SECTION(ODE_SECTION)) {
        if(MATCH_NAME("use_gpu")) {
            if(strcmp(value, "true") == 0 || strcmp(value, "yes") == 0) {
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
    } else if(SECTION_STARTS_WITH(STIM_SECTION)) {

        struct stim_config *tmp = (struct stim_config *)shget(pconfig->stim_configs, section);

        if(tmp == NULL) {
            tmp = new_stim_config();
            shput(pconfig->stim_configs, section, tmp);
        }

        if(MATCH_NAME("start")) {
            tmp->stim_start = (real)strtod(value, NULL);
            tmp->stim_start_was_set = true;
        } else if(MATCH_NAME("duration")) {
            tmp->stim_duration = (real)strtod(value, NULL);
            tmp->stim_duration_was_set = true;
        } else if(MATCH_NAME("current")) {
            tmp->stim_current = (real)strtod(value, NULL);
            tmp->stim_current_was_set = true;
        } else if(MATCH_NAME("period")) {
            tmp->stim_period = (real)strtod(value, NULL);
            tmp->stim_period_was_set = true;
        } else if(MATCH_NAME("function")) {
            tmp->config_data.function_name = strdup(value);
            tmp->config_data.function_name_was_set = true;
        } else if(MATCH_NAME("library_file")) {
            tmp->config_data.library_file_path = strdup(value);
            tmp->config_data.library_file_path_was_set = true;
        } else {
            // name is a reserved word in stim config
            if(MATCH_NAME("name")) {
                fprintf(stderr,
                        "name is a reserved word and should not be used inside a stimulus config section. Found in %s. "
                        "Exiting!\n",
                        section);
                exit(EXIT_FAILURE);
            } else {
                shput(tmp->config_data.config, name, strdup(value));
            }
        }
    } else if(MATCH_SECTION(DOMAIN_SECTION)) {

        if(pconfig->domain_config == NULL) {
            pconfig->domain_config = new_domain_config();
        }

        if(MATCH_NAME("start_dx")) {
            pconfig->domain_config->start_dx = strtod(value, NULL);
            pconfig->domain_config->start_dx_was_set = true;
        } else if(MATCH_NAME("maximum_dx")) {
            pconfig->domain_config->max_dx = strtod(value, NULL);
            pconfig->domain_config->max_dx_was_set = true;
        } else if(MATCH_NAME("start_dy")) {
            pconfig->domain_config->start_dy = strtod(value, NULL);
            pconfig->domain_config->start_dy_was_set = true;
        } else if(MATCH_NAME("maximum_dy")) {
            pconfig->domain_config->max_dy = strtod(value, NULL);
            pconfig->domain_config->max_dy_was_set = true;
        } else if(MATCH_NAME("start_dz")) {
            pconfig->domain_config->start_dz = strtod(value, NULL);
            pconfig->domain_config->start_dz_was_set = true;
        } else if(MATCH_NAME("maximum_dz")) {
            pconfig->domain_config->max_dz = strtod(value, NULL);
            pconfig->domain_config->max_dz_was_set = true;
        } else if(MATCH_NAME("name")) {
            pconfig->domain_config->domain_name = strdup(value);
            pconfig->domain_config->domain_name_was_set = true;
        } else if(MATCH_NAME("function")) {
            pconfig->domain_config->config_data.function_name = strdup(value);
            pconfig->domain_config->config_data.function_name_was_set = true;

        } else if(MATCH_NAME("library_file")) {
            pconfig->domain_config->config_data.library_file_path = strdup(value);
            pconfig->domain_config->config_data.library_file_path_was_set = true;
        } else {
            shput(pconfig->domain_config->config_data.config, name, strdup(value));
        }
    }
    else if(MATCH_SECTION(PURKINJE_SECTION)) {

        if(pconfig->purkinje_config == NULL) {
            pconfig->purkinje_config = new_purkinje_config();
        }

        if (MATCH_NAME ("start_discretization")) {
            pconfig->purkinje_config->start_h = strtod (value, NULL);
            pconfig->purkinje_config->start_h_was_set = true;

        }
        else if(MATCH_NAME("name")) {
            pconfig->purkinje_config->domain_name = strdup(value);
            pconfig->purkinje_config->domain_name_was_set = true;

        }
        else if(MATCH_NAME("function")) {
            pconfig->purkinje_config->config_data.function_name = strdup(value);
            pconfig->purkinje_config->config_data.function_name_was_set = true;

        }
        else if(MATCH_NAME("library_file")) {
            pconfig->purkinje_config->config_data.library_file_path = strdup(value);
            pconfig->purkinje_config->config_data.library_file_path_was_set = true;

        }
        else {
            shput(pconfig->purkinje_config->config_data.config, name, strdup(value));
        }
    } else if(MATCH_SECTION(MATRIX_ASSEMBLY_SECTION)) {

        if(pconfig->assembly_matrix_config == NULL) {
            pconfig->assembly_matrix_config = new_assembly_matrix_config();
        }

        if(MATCH_NAME("function")) {
            pconfig->assembly_matrix_config->config_data.function_name = strdup(value);
            pconfig->assembly_matrix_config->config_data.function_name_was_set = true;

        } else if(MATCH_NAME("library_file")) {
            pconfig->assembly_matrix_config->config_data.library_file_path = strdup(value);
            pconfig->assembly_matrix_config->config_data.library_file_path_was_set = true;
        } else {
            shput(pconfig->assembly_matrix_config->config_data.config, name, strdup(value));
        }
    }
    else if(MATCH_SECTION(UPDATE_MONODOMAIN_SECTION)) {

        if(pconfig->update_monodomain_config == NULL) {
            pconfig->update_monodomain_config = new_update_monodomain_config();
        }

        if(MATCH_NAME("function")) {
            pconfig->update_monodomain_config->config_data.function_name = strdup(value);
            pconfig->update_monodomain_config->config_data.function_name_was_set = true;

        } else if(MATCH_NAME("library_file")) {
            pconfig->update_monodomain_config->config_data.library_file_path = strdup(value);
            pconfig->update_monodomain_config->config_data.library_file_path_was_set = true;
        } else {
            shput(pconfig->update_monodomain_config->config_data.config, name, strdup(value));
        }
    }
    else if(MATCH_SECTION(LINEAR_SYSTEM_SOLVER_SECTION)) {

        if(pconfig->linear_system_solver_config == NULL) {
            pconfig->linear_system_solver_config = new_linear_system_solver_config();
        }

        if(MATCH_NAME("function")) {
            pconfig->linear_system_solver_config->config_data.function_name = strdup(value);
            pconfig->linear_system_solver_config->config_data.function_name_was_set = true;

        } else if(MATCH_NAME("library_file")) {
            pconfig->linear_system_solver_config->config_data.library_file_path = strdup(value);
            pconfig->linear_system_solver_config->config_data.library_file_path_was_set = true;
        } else {
            shput(pconfig->linear_system_solver_config->config_data.config, name, strdup(value));
        }
    } else if(MATCH_SECTION(EXTRA_DATA_SECTION)) {

        if(pconfig->extra_data_config == NULL) {
            pconfig->extra_data_config = new_extra_data_config();
        }

        if(MATCH_NAME("function")) {
            pconfig->extra_data_config->config_data.function_name = strdup(value);
            pconfig->extra_data_config->config_data.function_name_was_set = true;
        } else if(MATCH_NAME("library_file")) {
            pconfig->extra_data_config->config_data.library_file_path = strdup(value);
            pconfig->extra_data_config->config_data.library_file_path_was_set = true;
        } else {
            shput(pconfig->extra_data_config->config_data.config, name, strdup(value));
        }
    } else if(MATCH_SECTION(SAVE_RESULT_SECTION)) {

        if(pconfig->save_mesh_config == NULL) {
            pconfig->save_mesh_config = new_save_mesh_config();
        }

        if(MATCH_NAME("print_rate")) {
            pconfig->save_mesh_config->print_rate = (int)strtol(value, NULL, 10);
            pconfig->save_mesh_config->print_rate_was_set = true;
        } else if(MATCH_NAME("output_dir")) {
            pconfig->save_mesh_config->out_dir_name = strdup(value);
            pconfig->save_mesh_config->out_dir_name_was_set = true;
        } else if(MATCH_NAME("function")) {
            pconfig->save_mesh_config->config_data.function_name = strdup(value);
            pconfig->save_mesh_config->config_data.function_name_was_set = true;
        } else if(MATCH_NAME("library_file")) {
            pconfig->save_mesh_config->config_data.library_file_path = strdup(value);
            pconfig->save_mesh_config->config_data.library_file_path_was_set = true;
        } else {
            shput(pconfig->save_mesh_config->config_data.config, name, strdup(value));
        }
    } else if(MATCH_SECTION(SAVE_STATE_SECTION)) {

        if(pconfig->save_state_config == NULL) {
            pconfig->save_state_config = new_save_state_config();
        }

        if(MATCH_NAME("save_rate")) {
            pconfig->save_state_config->save_rate = (int)strtol(value, NULL, 10);
            pconfig->save_state_config->save_rate_was_set = true;
        } else if(MATCH_NAME("function")) {
            pconfig->save_state_config->config_data.function_name = strdup(value);
            pconfig->save_state_config->config_data.function_name_was_set = true;
        } else if(MATCH_NAME("library_file")) {
            pconfig->save_state_config->config_data.library_file_path = strdup(value);
            pconfig->save_state_config->config_data.library_file_path_was_set = true;
        } else {
            shput(pconfig->save_state_config->config_data.config, name, strdup(value));
        }
    } else if(MATCH_SECTION(RESTORE_STATE_SECTION)) {

        if(pconfig->restore_state_config == NULL) {
            pconfig->restore_state_config = new_restore_state_config();
        }

        if(MATCH_NAME("function")) {
            pconfig->restore_state_config->config_data.function_name = strdup(value);
            pconfig->restore_state_config->config_data.function_name_was_set = true;
        } else if(MATCH_NAME("library_file")) {
            pconfig->restore_state_config->config_data.library_file_path = strdup(value);
            pconfig->restore_state_config->config_data.library_file_path_was_set = true;
        } else {
            shput(pconfig->restore_state_config->config_data.config, name, strdup(value));
        }
    } else {

        fprintf(stderr, "\033[33;5;7mInvalid name %s in section %s on the config file!\033[0m\n", name, section);
        return 0;
    }

    return 1;
}

void configure_grid_from_options(struct grid *grid, struct user_options *options) {

    assert(grid);
    assert(options);

    grid->adaptive = options->adaptive;
}

void free_user_options(struct user_options *s) {

    free(s->model_file_path);

    if(s->stim_configs) {
        STIM_CONFIG_HASH_FOR_EACH_KEY_APPLY_FN_IN_VALUE(s->stim_configs, free_stim_config);
        shfree(s->stim_configs);
    }

    if(s->extra_data_config)
        free_extra_data_config(s->extra_data_config);

    if(s->domain_config)
        free_domain_config(s->domain_config);

    if(s->save_mesh_config)
        free_save_mesh_config(s->save_mesh_config);

    if(s->assembly_matrix_config)
        free_assembly_matrix_config(s->assembly_matrix_config);

    if(s->linear_system_solver_config)
        free_linear_system_solver_config(s->linear_system_solver_config);

    if(s->restore_state_config)
        free_restore_state_config(s->restore_state_config);

    if(s->save_state_config)
        free_save_state_config(s->save_state_config);

    if(s->update_monodomain_config)
        free_update_monodomain_config(s->update_monodomain_config);

    free(s);
}
