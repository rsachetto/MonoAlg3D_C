#include "config_parser.h"
#include "../../ini_parser/ini_file_sections.h"
#include "../../utils/logfile_utils.h"

#include <string.h>
#include <assert.h>

static const struct option long_options[] = {
        { "config_file", required_argument, NULL, 'c' },
        { "output_dir", required_argument , NULL, 'o' },
        { "use_adaptivity", no_argument, NULL, 'a' },
        { "abort_on_no_activity", no_argument, NULL, 'b' },
        { "num_threads", required_argument, NULL, 'n' },
        { "use_gpu", required_argument, NULL, 'g' },
        { "print_rate",required_argument , NULL, 'p' },
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
        { "binary_output", no_argument, NULL, 'y'},
        { "sigma_x", required_argument, NULL, SIGMA_X},
        { "sigma_y", required_argument, NULL, SIGMA_Y},
        { "sigma_z", required_argument, NULL, SIGMA_Z},
        { "beta", required_argument, NULL, BETA},
        { "cm", required_argument, NULL, CM},
        { "start_adapting_at", required_argument, NULL, START_REFINING},
        { "domain", required_argument, NULL, DOMAIN_OPT},
        { "extra_data", required_argument, NULL, EXTRA_DATA_OPT},
        { "stimulus", required_argument, NULL, STIM_OPT},
        { "draw_gl_output", no_argument, NULL, DRAW_OPT},
        { "help", no_argument, NULL, 'h' },
        { NULL, no_argument, NULL, 0 }
};

static const char *opt_string =   "c:o:abn:g:p:m:t:r:d:z:e:f:jR:D:G:k:yh";


/* Display program usage, and exit.
 */
void display_usage (char **argv) {

    //TODO: help for domain, extra_data and stimuls flags
    printf ("Usage: %s [options]\n\n", argv[0]);
    printf ("Options:\n");
    printf ("--config_file | -c [configuration_file_path]. Simulation configuration file. Default NULL.\n");
    printf ("--simulation_time | -f [simulation final time]. Simulation final time. Default 10.\n");
    printf ("--output_dir | -o [output-dir]. Simulation output directory. If not set, the simulation will not be saved "
                    "to disk. \n");
    printf ("--use_adaptivity | -a. No argument required. Use adaptivity. Default No use.\n");
    printf ("--abort_on_no_activity | -b. No argument required. The simulation will be aborted if no activity is "
                    "verified after print_rate time steps. Default false.\n");
    printf ("--print_rate | -p [output-print-rate]. Output print rate (in number of iterations). Default: 1 \n");
    printf ("--max_cg_its | -m [max-its]. Maximum number of CG iterations. Default: number of volumes \n");
    printf ("--cg_tolerance | -t [tolerance]. Conjugate Gradiente tolerance. Default: 1e-16 \n");
    printf ("--refinement_bound | -r [ref-bound]. ALG refinement bound (Vm variantion between time steps). Default: "
                    "5.0 \n");
    printf ("--derefinement_bound | -d [deref-bound]. ALG derefinement bound (Vm variantion between time steps). "
                    "Default: 1.0 \n");
    printf ("--dt_edp | -z [dt]. Simulation time discretization (PDE). Default: 0.01 \n");
    printf ("--dt_edo | -e [dt]. Minimum ODE time discretization (using time adaptivity. Default: 0.01 \n");
    printf ("--beta. Value of beta for simulation (Default: 0.14 \n");
    printf ("--cm. Value cm (Default: 1.0 \n");
    printf ("--num_threads | -n [num-threads]. Solve using OpenMP. Default: 1 \n");
    printf ("--use_gpu | -g [yes|no|true|false]. Solve ODEs using GPU. Default: No \n");
    printf ("--binary_output | -y. Save output files in binary format. Default: No \n");
    printf ("--use_preconditioner | -j Use Jacobi Preconditioner. Default: No \n");
    printf ("--refine_each | -R [ts], Refine each ts timesteps. Default: 1 \n");
    printf ("--derefine_each | -D [ts], Derefine each ts timesteps. Default: 1 \n");
    printf ("--gpu_id | -G [id], ID of the GPU to be used. Default: 0 \n");
    printf ("--model_file_path | -k [.so file path], Path of the .so representing the cell model and the ode solver. "
                    "Default: NULL \n");
    printf ("--draw_gl_output, Draw a iterative 3D output of the simulation. Not recommended for big meshes. Default: not draw\n");
    printf ("--help | -h. Shows this help and exit \n");
    exit (EXIT_FAILURE);
}


void issue_overwrite_warning (const char *var, const char *old_value, const char *new_value, const char *config_file) {
    fprintf (stderr,
             "WARNING: option %s was set in the file %s to %s and is being overwritten "
                     "by the command line flag to %s!\n",
             var, config_file, old_value, new_value);
}

struct user_options *new_user_options () {

    struct user_options *user_args = (struct user_options *)malloc (sizeof (struct user_options));

    user_args->out_dir_name = NULL;
    user_args->out_dir_name_was_set = false;

    user_args->adaptive = false;
    user_args->adaptive_was_set = false;

    user_args->num_threads = 1;
    user_args->num_threads_was_set = false;

    user_args->gpu = false;
    user_args->gpu_was_set = false;

    user_args->final_time = 10.0;
    user_args->final_time_was_set = false;

    user_args->print_rate = 1;
    user_args->print_rate_was_set = false;

    user_args->max_its = 50;
    user_args->max_its_was_set = false;

    user_args->ref_bound = 0.11;
    user_args->ref_bound_was_set = false;

    user_args->deref_bound = 0.10;
    user_args->deref_bound_was_set = false;

    user_args->dt_edp = 0.02;
    user_args->dt_edp_was_set = false;

    user_args->dt_edo = 0.01;
    user_args->dt_edo_was_set = false;

    user_args->use_jacobi = false;
    user_args->use_jacobi_was_set = false;

    user_args->refine_each = 1;
    user_args->refine_each_was_set = false;

    user_args->derefine_each = 1;
    user_args->derefine_each_was_set = false;

    user_args->gpu_id = 0;
    user_args->gpu_id_was_set = false;

    user_args->abort_no_activity = false;
    user_args->abort_no_activity_was_set = false;

    user_args->cg_tol = 1e-16;
    user_args->cg_tol_was_set = false;

    user_args->model_file_path = NULL;
    user_args->model_file_path_was_set = false;

    user_args->config_file = NULL;

    user_args->sigma_x = 0.0000176;
    user_args->sigma_x_was_set = false;

    user_args->sigma_y = 0.0001334;
    user_args->sigma_y_was_set = false;

    user_args->sigma_z = 0.0000176;
    user_args->sigma_z_was_set = false;

    user_args->beta = 0.14;
    user_args->beta_was_set = false;

    user_args->cm = 1.0;
    user_args->cm_was_set = false;

    user_args->start_adapting_at = 1.0;
    user_args->start_adapting_at_was_set = false;

    user_args->binary = false;
    user_args->binary_was_set = false;


    user_args->stim_configs = NULL;
    user_args->domain_config = NULL;
    user_args->extra_data_config = NULL;

    user_args->draw = false;

    return user_args;
}

void set_stim_config(const char *args, struct stim_config_hash *stim_configs, const char *config_file) {

    sds extra_config;
    sds *extra_config_tokens;
    int tokens_count;
    extra_config = sdsnew(args);
    extra_config_tokens = sdssplit(extra_config, ",", &tokens_count);
    char * opt_value;
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

    struct stim_config *sc = stim_config_hash_search(stim_configs, stim_name);

    if(sc == NULL) {
        sc = new_stim_config();
        print_to_stdout_and_file("Creating new stimulus name %s from command line options!\n", stim_name);
        stim_config_hash_insert(stim_configs, stim_name, sc);
    }

    struct string_hash *sh = sc->config_data.config;

    for(int i = 0; i < tokens_count; i++) {

        int values_count;
        sds *key_value = sdssplit(extra_config_tokens[i], "=", &values_count);

        key_value[0] = sdstrim(key_value[0], " ");
        key_value[1] = sdstrim(key_value[1], " ");

        key   = key_value[0];
        value = key_value[1];

        if(strcmp(key, "start") == 0) {
            if(sc->stim_start_was_set) {
                sprintf (old_value, "%lf", sc->stim_start);
                print_to_stdout_and_file("WARNING: For stimulus %s:\n", stim_name);
                issue_overwrite_warning ("start", old_value, value, config_file);
            }
            sc->stim_start = (real)strtod(value, NULL);
        }
        else if (strcmp(key, "duration") == 0) {
            if(sc->stim_duration_was_set) {
                sprintf (old_value, "%lf", sc->stim_duration);
                print_to_stdout_and_file("WARNING: For stimulus %s:\n", stim_name);
                issue_overwrite_warning ("duration", old_value, value, config_file);
            }

            sc->stim_duration = (real)strtod(value, NULL);
        }
        else if (strcmp(key, "current") == 0) {
            if(sc->stim_current_was_set) {
                sprintf (old_value, "%lf", sc->stim_current);
                print_to_stdout_and_file("WARNING: For stimulus %s:\n", stim_name);
                issue_overwrite_warning ("current", old_value, value, config_file);
            }
            sc->stim_current = (real)strtod(value, NULL);
        }
        else if (strcmp(key, "function") == 0) {
            if(sc->config_data.function_name_was_set) {
                print_to_stdout_and_file("WARNING: For stimulus %s:\n", stim_name);
                issue_overwrite_warning ("function", sc->config_data.function_name, value, config_file);
            }
            free(sc->config_data.function_name);
            sc->config_data.function_name = strdup(value);
        }
        else if (strcmp(key, "library_file") == 0) {
            if(sc->config_data.library_file_path_was_set) {
                print_to_stdout_and_file("WARNING: For stimulus %s:\n", stim_name);
                issue_overwrite_warning ("library_file", sc->config_data.library_file_path, value, config_file);
            }
            free(sc->config_data.library_file_path);
            sc->config_data.library_file_path = strdup(value);
        }
        else {
            opt_value = string_hash_search(sh, key);
            if(opt_value) {
                issue_overwrite_warning(key, opt_value, value, config_file);
            }
            string_hash_insert_or_overwrite(sh, key, value);
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
    char * opt_value;
    char old_value[32];
    char *key, *value;

    assert(dc);

    struct string_hash *sh = dc->config_data.config;

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

        key   = key_value[0];
        value = key_value[1];

        if(strcmp(key, "name") == 0) {
            if(dc->domain_name_was_set) {
                print_to_stdout_and_file("WARNING: For domain configuration: \n");
                issue_overwrite_warning ("name", dc->domain_name, value, config_file);
            }
            free(dc->domain_name);
            dc->domain_name = strdup(value);
        }
        else if (strcmp(key, "start_discretization") == 0) {
            if(dc->start_h_was_set) {
                sprintf (old_value, "%lf", dc->start_h);
                print_to_stdout_and_file("WARNING: For domain configuration: \n");
                issue_overwrite_warning ("start_discretization", old_value, value, config_file);
            }
            dc->start_h = (real)strtod(value, NULL);
        }
        else if (strcmp(key, "maximum_discretization") == 0) {
            if(dc->max_h_was_set) {
                sprintf (old_value, "%lf", dc->max_h);
                print_to_stdout_and_file("WARNING: For domain configuration: \n");
                issue_overwrite_warning ("maximum_discretization", old_value, value, config_file);
            }
            dc->max_h = (real)strtod(value, NULL);
        }
        else if (strcmp(key, "function") == 0) {
            if(dc->config_data.function_name_was_set) {
                print_to_stdout_and_file("WARNING: For domain configuration: \n");
                issue_overwrite_warning ("function", dc->config_data.function_name, value, config_file);
            }
            free(dc->config_data.function_name);
            dc->config_data.function_name = strdup(value);
        }
        else if (strcmp(key, "library_file") == 0) {
            if(dc->config_data.library_file_path_was_set) {
                print_to_stdout_and_file("WARNING: For domain configuration: \n");
                issue_overwrite_warning ("library_file", dc->config_data.library_file_path, value, config_file);
            }
            free(dc->config_data.library_file_path);
            dc->config_data.library_file_path = strdup(value);
        }
        else {
            opt_value = string_hash_search(sh, key);
            if(opt_value) {
                issue_overwrite_warning(key, opt_value, value, config_file);
            }

            string_hash_insert_or_overwrite(sh, key, value);

        }
        sdsfreesplitres(key_value, values_count);

    }

    sdsfreesplitres(extra_config_tokens, tokens_count);

}

void set_extra_data_config(const char *args, struct extra_data_config *dc, const char *config_file) {

    sds extra_config;
    sds *extra_config_tokens;
    int tokens_count;
    extra_config = sdsnew(args);
    extra_config_tokens = sdssplit(extra_config, ",", &tokens_count);
    char * opt_value;
    char *key, *value;

    assert(dc);

    struct string_hash *sh = dc->config_data.config;

    for(int i = 0; i < tokens_count; i++) {
        extra_config_tokens[i] = sdstrim(extra_config_tokens[i], " ");

        int values_count;
        sds *key_value = sdssplit(extra_config_tokens[i], "=", &values_count);

        if(values_count != 2) {
            fprintf(stderr, "Invalid format for optios %s. Exiting!\n", args);
            exit(EXIT_FAILURE);
        }

        key_value[0] = sdstrim(key_value[0], " ");
        key_value[1] = sdstrim(key_value[1], " ");

        key   = key_value[0];
        value = key_value[1];

        if (strcmp(key, "function") == 0) {
            if(dc->config_data.function_name_was_set) {
                print_to_stdout_and_file("WARNING: For extra_data configuration: \n");
                issue_overwrite_warning ("function", dc->config_data.function_name, value, config_file);
            }
            free(dc->config_data.function_name);
            dc->config_data.function_name = strdup(value);
        }
        else if (strcmp(key, "library_file") == 0) {
            if(dc->config_data.library_file_path_was_set) {
                print_to_stdout_and_file("WARNING: For extra_data configuration: \n");
                issue_overwrite_warning ("library_file", dc->config_data.library_file_path, value, config_file);
            }
            free(dc->config_data.library_file_path);
            dc->config_data.library_file_path = strdup(value);
        }
        else {
            opt_value = string_hash_search(sh, key);
            if(opt_value) {
                issue_overwrite_warning(key, opt_value, value, config_file);
            }

            string_hash_insert_or_overwrite(sh, key, value);

        }
        sdsfreesplitres(key_value, values_count);

    }

    sdsfreesplitres(extra_config_tokens, tokens_count);

}

void get_config_file (int argc, char **argv, struct user_options *user_args) {
    int opt = 0;

    int option_index;

    opt = getopt_long (argc, argv, opt_string, long_options, &option_index);

    while (opt != -1) {
        switch (opt) {
            case 'c':
                user_args->config_file = optarg;
                return;
            default:
                break;
        }
        opt = getopt_long (argc, argv, opt_string, long_options, &option_index);
    }

    // We reset the index after parsing the config_file
    optind = 1;
}


void parse_options (int argc, char **argv, struct user_options *user_args) {

    int opt = 0;
    int option_index;

    opt =  getopt_long_only (argc, argv, opt_string, long_options, &option_index);
    char old_value[32];

    while (opt != -1) {
        switch (opt) {
            case 'R':
                if (user_args->refine_each_was_set) {
                    sprintf (old_value, "%d", user_args->refine_each);
                    issue_overwrite_warning ("refine_each", old_value, optarg, user_args->config_file);
                }
                user_args->refine_each = (int)strtol (optarg, NULL, 10);

                break;
            case 'D':
                if (user_args->derefine_each_was_set) {
                    sprintf (old_value, "%d", user_args->derefine_each);
                    issue_overwrite_warning ("derefine_each", old_value, optarg, user_args->config_file);
                }
                user_args->derefine_each =  (int)strtol (optarg, NULL, 10);

                break;
            case 'e':
                if (user_args->dt_edo_was_set) {
                    sprintf (old_value, "%lf", user_args->dt_edo);
                    issue_overwrite_warning ("dt_edo", old_value, optarg, user_args->config_file);
                }
                user_args->dt_edo = strtod (optarg, NULL);

                break;
            case 'm':
                if (user_args->max_its_was_set) {
                    sprintf (old_value, "%d", user_args->max_its);
                    issue_overwrite_warning ("max_cg_its", old_value, optarg, user_args->config_file);
                }
                user_args->max_its = (int)strtol (optarg, NULL, 10);

                break;
            case 'p':
                if (user_args->print_rate_was_set) {
                    sprintf (old_value, "%d", user_args->print_rate);
                    issue_overwrite_warning ("print_rate", old_value, optarg, user_args->config_file);
                }
                user_args->print_rate = (int)strtol (optarg, NULL, 10);
                break;
            case 'o':
                if (user_args->out_dir_name_was_set) {
                    if (user_args->out_dir_name) {
                        issue_overwrite_warning ("output_dir", user_args->out_dir_name, optarg, user_args->config_file);
                    } else {
                        issue_overwrite_warning ("output_dir", "No Save", optarg, user_args->config_file);
                    }
                }
                free(user_args->out_dir_name);
                user_args->out_dir_name = strdup(optarg);

                break;
            case 'k':
                if (user_args->model_file_path_was_set) {
                    if (user_args->model_file_path) {
                        issue_overwrite_warning ("model_file_path", user_args->model_file_path, optarg,
                                                 user_args->config_file);
                    } else {
                        issue_overwrite_warning ("model_file_path", "No Save", optarg, user_args->config_file);
                    }
                }
                free(user_args->model_file_path);
                user_args->model_file_path = strdup(optarg);

                break;
            case 'f':
                if (user_args->final_time_was_set) {
                    sprintf (old_value, "%lf", user_args->final_time);
                    issue_overwrite_warning ("simulation_time", old_value, optarg, user_args->config_file);
                }
                user_args->final_time = strtod (optarg, NULL);

                break;
            case 'r':
                if (user_args->ref_bound_was_set) {
                    sprintf (old_value, "%lf", user_args->ref_bound);
                    issue_overwrite_warning ("refinement_bound", old_value, optarg, user_args->config_file);
                }
                user_args->ref_bound = strtod (optarg, NULL);
                break;
            case 'd':
                if (user_args->deref_bound_was_set) {
                    sprintf (old_value, "%lf", user_args->deref_bound);
                    issue_overwrite_warning ("derefinement_bound", old_value, optarg, user_args->config_file);
                }
                user_args->deref_bound = strtod (optarg, NULL);
                break;
            case 'a':
                if (user_args->adaptive_was_set) {
                    sprintf (old_value, "%d", user_args->adaptive);
                    issue_overwrite_warning ("use_adaptivity", old_value, optarg, user_args->config_file);
                }
                user_args->adaptive = true;
                break;
            case 'n':
                if (((int)strtol (optarg, NULL, 10)) > 0) {
                    if (user_args->num_threads_was_set) {
                        sprintf (old_value, "%d", user_args->num_threads);
                        issue_overwrite_warning ("nu"
                                                         "m_threads", old_value, optarg, user_args->config_file);
                    }
                    user_args->num_threads = (int)strtol (optarg, NULL, 10);
                }
                break;
            case 'g':
                if (user_args->gpu_was_set) {

                    if(user_args->gpu) {
                        sprintf(old_value, "yes");
                    }
                    else {
                        sprintf(old_value, "no");
                    }

                    issue_overwrite_warning ("use_gpu", old_value, optarg, user_args->config_file);
                }
                if (strcmp(optarg, "true") == 0 || strcmp(optarg, "yes") == 0) {
                    user_args->gpu = true;
                } else if (strcmp(optarg, "false") == 0 || strcmp(optarg, "no") == 0) {
                    user_args->gpu = false;
                }
                else {
                    fprintf(stderr, "Warning: Invalid value for use_gpu option: %s! Valid options are: true, yes, false, no. Setting the value to false\n", optarg);
                    user_args->gpu = false;
                }
                break;
            case 'z':
                if (user_args->dt_edp_was_set) {
                    sprintf (old_value, "%lf", user_args->dt_edp);
                    issue_overwrite_warning ("dt_edp", old_value, optarg, user_args->config_file);
                }
                user_args->dt_edp = strtod (optarg, NULL);
                break;
            case 't':
                if (user_args->cg_tol_was_set) {
                    sprintf (old_value, "%lf", user_args->cg_tol);
                    issue_overwrite_warning ("cg_tolerance", old_value, optarg, user_args->config_file);
                }
                user_args->cg_tol = strtod (optarg, NULL);
                break;
            case SIGMA_X:
                if (user_args->sigma_x_was_set) {
                    sprintf (old_value, "%lf", user_args->sigma_x);
                    issue_overwrite_warning ("sigma_x", old_value, optarg, user_args->config_file);
                }
                user_args->sigma_x = strtod (optarg, NULL);
                break;
            case SIGMA_Y:
                if (user_args->sigma_y_was_set) {
                    sprintf (old_value, "%lf", user_args->sigma_y);
                    issue_overwrite_warning ("sigma_y", old_value, optarg, user_args->config_file);
                }
                user_args->sigma_y = strtod (optarg, NULL);
                break;
            case SIGMA_Z:
                if (user_args->sigma_z_was_set) {
                    sprintf (old_value, "%lf", user_args->sigma_z);
                    issue_overwrite_warning ("sigma_z", old_value, optarg, user_args->config_file);
                }
                user_args->sigma_z = strtod (optarg, NULL);
                break;
            case BETA:
                if (user_args->beta_was_set) {
                    sprintf (old_value, "%lf", user_args->beta);
                    issue_overwrite_warning ("beta", old_value, optarg, user_args->config_file);
                }
                user_args->beta = strtod (optarg, NULL);
                break;
            case CM:
                if (user_args->cm) {
                    sprintf (old_value, "%lf", user_args->cm);
                    issue_overwrite_warning ("cm", old_value, optarg, user_args->config_file);
                }
                user_args->cm = strtod (optarg, NULL);
                break;
            case START_REFINING:
                if (user_args->start_adapting_at_was_set) {
                    sprintf (old_value, "%lf", user_args->start_adapting_at);
                    issue_overwrite_warning ("start_adapting_at", old_value, optarg, user_args->config_file);
                }
                user_args->start_adapting_at = strtod (optarg, NULL);
                break;
            case 'G':
                if (user_args->gpu_id_was_set) {
                    sprintf (old_value, "%d", user_args->gpu_id);
                    issue_overwrite_warning ("gpu_id", old_value, optarg, user_args->config_file);
                }
                user_args->gpu_id = (int)strtol (optarg, NULL, 10);
                break;
            case 'b':
                if (user_args->abort_no_activity_was_set) {
                    sprintf (old_value, "%d", user_args->abort_no_activity);
                    issue_overwrite_warning ("abort_on_no_activity", old_value, optarg, user_args->config_file);
                }
                user_args->abort_no_activity = true;
                break;
            case 'j':
                if (user_args->use_jacobi_was_set) {
                    sprintf (old_value, "%d", user_args->use_jacobi);
                    issue_overwrite_warning ("use_preconditioner", old_value, optarg, user_args->config_file);
                }
                user_args->use_jacobi = true;
                break;
            case 'y':
                if (user_args->binary_was_set) {
                    sprintf (old_value, "%d", user_args->binary);
                    issue_overwrite_warning ("binary_output", old_value, optarg, user_args->config_file);
                }
                user_args->binary = true;
                break;
            case DOMAIN_OPT:
                if(user_args->domain_config == NULL) {
                    print_to_stdout_and_file("Creating new domain config from command line!\n");
                    user_args->domain_config = new_domain_config();
                }
                set_domain_config(optarg, user_args->domain_config, user_args->config_file );
                break;
            case EXTRA_DATA_OPT:
                if(user_args->extra_data_config == NULL) {
                    print_to_stdout_and_file("Creating new extra data config from command line!\n");
                    user_args->extra_data_config = new_extra_data_config();
                }
                set_extra_data_config(optarg, user_args->extra_data_config, user_args->config_file );
                break;
            case STIM_OPT:
                if(user_args->stim_configs == NULL) {
                    print_to_stdout_and_file("Creating new stim config from command line!\n");
                    user_args->stim_configs = stim_config_hash_create();
                }
                set_stim_config(optarg, user_args->stim_configs, user_args->config_file );
                break;
            case DRAW_OPT:
                user_args->draw = true;
                break;
            case 'h': /* fall-through is intentional */
            case '?':
                display_usage(argv);
                break;
            default:
                /* You won't actually get here. */
                break;
        }

        opt = getopt_long (argc, argv, opt_string, long_options, &option_index);
    }
}

int parse_config_file (void *user, const char *section, const char *name, const char *value) {
    struct user_options *pconfig = (struct user_options *) user;

    if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "cg_tolerance")) {
        pconfig->cg_tol = strtod(value, NULL);
        pconfig->cg_tol_was_set = true;
    } else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "max_cg_its")) {
        pconfig->max_its = (int)strtol (value, NULL, 10);
        pconfig->max_its_was_set = true;

    } else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "use_preconditioner")) {
        if (strcmp(value, "true") == 0 || strcmp(value, "yes") == 0) {
            pconfig->use_jacobi = true;
        } else {
            pconfig->use_jacobi = false;
        }
        pconfig->use_jacobi_was_set = true;
    } else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "num_threads")) {
        pconfig->num_threads = (int)strtol (value, NULL, 10);
        pconfig->num_threads_was_set = true;
    } else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "dt_edp")) {
        pconfig->dt_edp = strtod(value, NULL);
        pconfig->dt_edp_was_set = true;
    } else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "simulation_time")) {
        pconfig->final_time = strtod(value, NULL);
        pconfig->final_time_was_set = true;
    } else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "sigma_x")) {
        pconfig->sigma_x = strtod(value, NULL);
        pconfig->sigma_x_was_set = true;
    } else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "sigma_y")) {
        pconfig->sigma_y = strtod(value, NULL);
        pconfig->sigma_y_was_set = true;
    } else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "sigma_z")) {
        pconfig->sigma_z = strtod(value, NULL);
        pconfig->sigma_z_was_set = true;
    } else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "beta")) {
        pconfig->beta = strtod(value, NULL);
        pconfig->beta_was_set = true;
    } else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "cm")) {
        pconfig->cm = strtod(value, NULL);
        pconfig->cm_was_set = true;
    } else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "start_adapting_at")) {
        pconfig->start_adapting_at = strtod(value, NULL);
        pconfig->start_adapting_at_was_set = true;
    } else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "abort_on_no_activity")) {
        if (strcmp(value, "true") == 0 || strcmp(value, "yes") == 0) {
            pconfig->abort_no_activity = true;
        } else {
            pconfig->abort_no_activity = false;
        }
        pconfig->abort_no_activity_was_set = true;
    } else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "use_adaptivity")) {
        if (strcmp(value, "true") == 0 || strcmp(value, "yes") == 0) {
            pconfig->adaptive = true;
        } else {
            pconfig->adaptive = false;
        }
        pconfig->adaptive_was_set = true;
    } else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "print_rate")) {
        pconfig->print_rate = (int)strtol (value, NULL, 10);
        pconfig->print_rate_was_set = true;
    } else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "output_dir")) {
        pconfig->out_dir_name = strdup(value);
        pconfig->out_dir_name_was_set = true;
    } else if (MATCH_SECTION_AND_NAME (MAIN_SECTION, "binary_output")) {
        if (strcmp(value, "true") == 0 || strcmp(value, "yes") == 0) {
            pconfig->binary = true;
        } else {
            pconfig->binary = false;
        }
        pconfig->binary_was_set = true;
    } else if (MATCH_SECTION_AND_NAME (ALG_SECTION, "refinement_bound")) {
        pconfig->ref_bound = strtod(value, NULL);
        pconfig->ref_bound_was_set = true;
    } else if (MATCH_SECTION_AND_NAME (ALG_SECTION, "derefinement_bound")) {
        pconfig->deref_bound = strtod(value, NULL);
        pconfig->deref_bound_was_set = true;
    } else if (MATCH_SECTION_AND_NAME (ALG_SECTION, "refine_each")) {
        pconfig->refine_each = (int)strtol (value, NULL, 10);
        pconfig->refine_each_was_set = true;
    } else if (MATCH_SECTION_AND_NAME (ALG_SECTION, "derefine_each")) {
        pconfig->derefine_each = (int)strtol (value, NULL, 10);
        pconfig->derefine_each_was_set = true;
    } else if (MATCH_SECTION_AND_NAME (ODE_SECTION, "dt_edo")) {
        pconfig->dt_edo = strtod(value, NULL);
        pconfig->dt_edo_was_set = true;
    } else if (MATCH_SECTION(ODE_SECTION)) {
        if(MATCH_NAME("use_gpu")) {
            if (strcmp(value, "true") == 0 || strcmp(value, "yes") == 0) {
                pconfig->gpu = true;
            } else {
                pconfig->gpu = false;
            }
            pconfig->gpu_was_set = true;
        } else if (MATCH_NAME ("gpu_id")) {
            pconfig->gpu_id = (int)strtol (value, NULL, 10);
            pconfig->gpu_id_was_set = true;
        } else if (MATCH_NAME ("library_file")) {
            pconfig->model_file_path = strdup(value);
            pconfig->model_file_path_was_set = true;
        }
    }
    else if(SECTION_STARTS_WITH(STIM_SECTION)) {

        if(pconfig->stim_configs == NULL) {
            pconfig->stim_configs = stim_config_hash_create();
        }

        struct stim_config *tmp = stim_config_hash_search(pconfig->stim_configs, section);

        if( tmp == NULL) {
            tmp = new_stim_config();
            stim_config_hash_insert(pconfig->stim_configs, section, tmp);
        }

        if(MATCH_NAME("start")) {
            tmp->stim_start = (real)strtod(value, NULL);
            tmp->stim_start_was_set = true;
        } else if(MATCH_NAME("duration")) {
            tmp->stim_duration = (real)strtod(value, NULL);
            tmp->stim_duration_was_set = true;
        }else if(MATCH_NAME("current")) {
            tmp->stim_current = (real)strtod(value, NULL);
            tmp->stim_current_was_set = true;
        } else if(MATCH_NAME("function")) {
            tmp->config_data.function_name = strdup(value);
            tmp->config_data.function_name_was_set = true;
        } else if(MATCH_NAME("library_file")) {
            tmp->config_data.library_file_path = strdup(value);
            tmp->config_data.library_file_path_was_set = true;
        }
        else {
            //name is a reserved word in stim config
            if(MATCH_NAME("name")) {
                fprintf(stderr, "name is a reserved word and should not be used inside a stimulus config section. Found in %s. Exiting!\n", section);
                exit(EXIT_FAILURE);
            }
            else {
                string_hash_insert(tmp->config_data.config, name, value);
            }
        }
    }
    else if(MATCH_SECTION(DOMAIN_SECTION)) {

        if(pconfig->domain_config == NULL) {
            pconfig->domain_config = new_domain_config();
        }

        if (MATCH_NAME ("start_discretization")) {
            pconfig->domain_config->start_h = strtod (value, NULL);
            pconfig->domain_config->start_h_was_set = true;
        } else if (MATCH_NAME ("maximum_discretization")) {
            pconfig->domain_config->max_h = strtod (value, NULL);
            pconfig->domain_config->max_h_was_set = true;
        }
        else if(MATCH_NAME("name")) {
            pconfig->domain_config->domain_name = strdup(value);
            pconfig->domain_config->domain_name_was_set = true;
        }
        else if(MATCH_NAME("function")) {
            pconfig->domain_config->config_data.function_name = strdup(value);
            pconfig->domain_config->config_data.function_name_was_set = true;

        }
        else if(MATCH_NAME("library_file")) {
            pconfig->domain_config->config_data.library_file_path = strdup(value);
            pconfig->domain_config->config_data.library_file_path_was_set = true;
        }
        else {
            string_hash_insert(pconfig->domain_config->config_data.config, name, value);
        }
    }
    else if(MATCH_SECTION(EXTRA_DATA_SECTION)) {

        if(pconfig->extra_data_config == NULL) {
            pconfig->extra_data_config = new_extra_data_config();
        }

        if(MATCH_NAME("function")) {
            pconfig->extra_data_config->config_data.function_name = strdup(value);
            pconfig->extra_data_config->config_data.function_name_was_set = true;
        }
        else if(MATCH_NAME("library_file")) {
            pconfig->extra_data_config->config_data.library_file_path = strdup(value);
            pconfig->extra_data_config->config_data.library_file_path_was_set = true;
        }
        else {
            string_hash_insert(pconfig->extra_data_config->config_data.config, name, value);
        }
    }
    else {
        fprintf(stderr, "Invalid name %s in section %s on the config file!\n", name, section);
        return 0;
    }

    return 1;
}

void configure_grid_from_options(struct grid* grid, struct user_options *options) {
    assert(grid);
    assert(options);

    grid->adaptive = options->adaptive;
}


void free_user_options(struct user_options *s) {
    free(s->model_file_path);
    free(s->out_dir_name);
    if(s->stim_configs) {
        STIM_CONFIG_HASH_FOR_EACH_KEY_APPLY_FN_IN_VALUE(s->stim_configs, free_stim_config);
        stim_config_hash_destroy(s->stim_configs);
    }
    if(s->extra_data_config)
        free_extra_data_config(s->extra_data_config);
    if(s->domain_config)
        free_domain_config(s->domain_config);
    free(s);
}