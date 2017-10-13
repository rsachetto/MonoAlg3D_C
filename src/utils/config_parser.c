#include <getopt.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "ini_parser/ini_file_sections.h"
#include "config_parser.h"
#include <string.h>

static const char *optString = "c:l:k:p:abn:gp:s:x:m:t:r:d:z:e:f:l:R:D:G:h?";

/* Display program usage, and exit.
 */
void display_usage (char **argv) {

    printf ("Usage: %s [options]\n\n", argv[0]);
    printf ("Options:\n");
    printf ("--config_file | -c [configuration_file_path]. Simulation final time. Default NULL.\n");
    printf ("--simulation_time | -f [simulation final time]. Simulation final time. Default 10.\n");
    printf ("--side_lenght | -l [length]. Tissue side lenght (micro m). Default: 12800 micro m \n");
    printf ("--output_dir | -o [output-dir]. Simulation output directory. If not set, the simulation will not be saved "
            "to disk. \n");
    printf ("--use_adaptivity | -a. No argument required. Use adaptivity. Default No use.\n");
    printf ("--abort_on_no_activity | -b. No argument required. The simulation will be aborted if no activity is "
            "verified after print_rate time steps. Default false.\n");
    printf ("--print_rate | -p [output-print-rate]. Output print rate (in number of iterations). Default: 1 \n");
    printf ("--start_discretization | -s [start-h]. Initial space discretization (micro m). Default: 12.5 micro m \n");
    printf ("--maximum_discretization | -x [max-h].Maximum space discretization (micro m). Default: 200.0 micro m \n");
    printf ("--max_cg_its | -m [max-its]. Maximum number of CG iterations. Default: number of volumes \n");
    printf ("--cg_tolerance | -t [tolerance]. Conjugate Gradiente tolerance. Default: 1e-16 \n");
    printf ("--refinement_bound | -r [ref-bound]. ALG refinement bound (Vm variantion between time steps). Default: "
            "5.0 \n");
    printf ("--derefinement_bound | -d [deref-bound]. ALG derefinement bound (Vm variantion between time steps). "
            "Default: 1.0 \n");
    printf ("--dt_edp | -z [dt]. Simulation time discretization (PDE). Default: 0.01 \n");
    printf ("--dt_edo | -e [dt]. Minimum ODE time discretization (using time adaptivity. Default: 0.01 \n");
    printf ("--num_threads | -n [num-threads]. Solve using OpenMP. Default: 1 \n");
    printf ("--use_gpu | -g. Solve ODEs using GPU. Default: No \n");
    printf ("--use_preconditioner | -j Use Jacobi Preconditioner. Default: No \n");
    printf ("--refine_each | -R [ts], Refine each ts timesteps. Default: 1 \n");
    printf ("--derefine_each | -D [ts], Derefine each ts timesteps. Default: 1 \n");
    printf ("--gpu_id | -G [id], ID of the GPU to be used. Default: 0 \n");
    printf ("--model_file_path | -k [.so file path], Path of the .so representing the cell model and the ode solver. "
            "Default: NULL \n");
    printf ("--help | -h. Shows this help and exit \n");
    exit (EXIT_FAILURE);
}

struct user_options *new_user_options () {

    struct user_options *user_args = (struct user_options *)malloc (sizeof (struct user_options));

    user_args->out_dir_name = NULL;
    user_args->out_dir_name_was_set = false;

    user_args->side_lenght = 12800.0;
    user_args->side_lenght_was_set = false;

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

    user_args->start_h = 100.0;
    user_args->start_h_was_set = false;

    user_args->max_h = 200.0;
    user_args->max_h_was_set = false;

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

    return user_args;
}

void get_config_file (int argc, char **argv, struct user_options *user_args) {
    int opt = 0;

    int longIndex;

    opt = getopt_long (argc, argv, optString, longOpts, &longIndex);

    while (opt != -1) {
        switch (opt) {
        case 'c':
            user_args->config_file = optarg;
            return;
        }
        opt = getopt_long (argc, argv, optString, longOpts, &longIndex);
    }
}

void issue_overwrite_warning (const char *var, const char *old_value, const char *new_value, const char *config_file) {
    fprintf (stderr,
             "WARNING: option %s was set in the file %s to %s and is being overwritten "
             "by the command line flag to %s!\n",
             var, config_file, old_value, new_value);
}

void parse_options (int argc, char **argv, struct user_options *user_args) {
    int opt = 0;

    // We reset the index after parsing the config_file
    optind = 1;

    int longIndex;

    opt = getopt_long (argc, argv, optString, longOpts, &longIndex);
    char old_value[32];

    while (opt != -1) {
        switch (opt) {

        case 'R':
            if (user_args->refine_each_was_set) {
                sprintf (old_value, "%d", user_args->refine_each);
                issue_overwrite_warning ("refine_each", old_value, optarg, user_args->config_file);
            }
            user_args->refine_each = atoi (optarg);

            break;
        case 'D':
            if (user_args->derefine_each_was_set) {
                sprintf (old_value, "%d", user_args->derefine_each);
                issue_overwrite_warning ("derefine_each", old_value, optarg, user_args->config_file);
            }
            user_args->derefine_each = atoi (optarg);

            break;
        case 'e':
            if (user_args->dt_edo_was_set) {
                sprintf (old_value, "%lf", user_args->dt_edo);
                issue_overwrite_warning ("dt_edo", old_value, optarg, user_args->config_file);
            }
            user_args->dt_edo = atof (optarg);

            break;
        case 'm':
            if (user_args->max_its_was_set) {
                sprintf (old_value, "%d", user_args->max_its);
                issue_overwrite_warning ("max_cg_its", old_value, optarg, user_args->config_file);
            }
            user_args->max_its = atoi (optarg);

            break;
        case 'l':
            if (user_args->side_lenght_was_set) {
                sprintf (old_value, "%lf", user_args->side_lenght);
                issue_overwrite_warning ("side_lenght", old_value, optarg, user_args->config_file);
            }
            user_args->side_lenght = atof (optarg);

            break;
        case 'p':
            if (user_args->print_rate_was_set) {
                sprintf (old_value, "%d", user_args->print_rate);
                issue_overwrite_warning ("print_rate", old_value, optarg, user_args->config_file);
            }
            user_args->print_rate = atoi (optarg);

            break;
        case 'o':
            if (user_args->out_dir_name_was_set) {
                if (user_args->out_dir_name) {
                    issue_overwrite_warning ("output_dir", user_args->out_dir_name, optarg, user_args->config_file);
                } else {
                    issue_overwrite_warning ("output_dir", "No Save", optarg, user_args->config_file);
                }
            }
            user_args->out_dir_name = optarg;

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
            user_args->model_file_path = optarg;

            break;
        case 'f':
            if (user_args->final_time_was_set) {
                sprintf (old_value, "%lf", user_args->final_time);
                issue_overwrite_warning ("simulation_time", old_value, optarg, user_args->config_file);
            }
            user_args->final_time = atof (optarg);

            break;
        case 's':
            if (user_args->start_h_was_set) {
                sprintf (old_value, "%lf", user_args->start_h);
                issue_overwrite_warning ("start_discretization", old_value, optarg, user_args->config_file);
            }
            user_args->start_h = atof (optarg);
            break;
        case 'x':
            if (user_args->max_h_was_set) {
                sprintf (old_value, "%lf", user_args->max_h);
                issue_overwrite_warning ("maximum_discretization", old_value, optarg, user_args->config_file);
            }
            user_args->max_h = atof (optarg);
            break;
        case 'r':
            if (user_args->ref_bound_was_set) {
                sprintf (old_value, "%lf", user_args->ref_bound);
                issue_overwrite_warning ("refinement_bound", old_value, optarg, user_args->config_file);
            }
            user_args->ref_bound = atof (optarg);
            break;
        case 'd':
            if (user_args->deref_bound_was_set) {
                sprintf (old_value, "%lf", user_args->deref_bound);
                issue_overwrite_warning ("derefinement_bound", old_value, optarg, user_args->config_file);
            }
            user_args->deref_bound = atof (optarg);
            break;
        case 'a':
            if (user_args->adaptive_was_set) {
                sprintf (old_value, "%d", user_args->adaptive);
                issue_overwrite_warning ("use_adaptivity", old_value, optarg, user_args->config_file);
            }
            user_args->adaptive = true;
            break;
        case 'n':
            if (atoi (optarg) > 0) {
                if (user_args->num_threads_was_set) {
                    sprintf (old_value, "%d", user_args->num_threads);
                    issue_overwrite_warning ("num_threads", old_value, optarg, user_args->config_file);
                }
                user_args->num_threads = atoi (optarg);
            }
            break;
        case 'g':
            if (user_args->gpu_was_set) {
                sprintf (old_value, "%d", user_args->gpu);
                issue_overwrite_warning ("use_gpu", old_value, optarg, user_args->config_file);
            }
            user_args->gpu = true;
            break;
        case 'z':
            if (user_args->dt_edp_was_set) {
                sprintf (old_value, "%lf", user_args->dt_edp);
                issue_overwrite_warning ("dt_edp", old_value, optarg, user_args->config_file);
            }
            user_args->dt_edp = atof (optarg);
            break;
        case 't':
            if (user_args->cg_tol_was_set) {
                sprintf (old_value, "%lf", user_args->cg_tol);
                issue_overwrite_warning ("cg_tolerance", old_value, optarg, user_args->config_file);
            }
            user_args->cg_tol = atof (optarg);
            break;
        case 'G':
            if (user_args->gpu_id_was_set) {
                sprintf (old_value, "%d", user_args->gpu_id);
                issue_overwrite_warning ("gpu_id", old_value, optarg, user_args->config_file);
            }
            user_args->gpu_id = atoi (optarg);
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
        case 'h': /* fall-through is intentional */
        case '?':
            display_usage (argv);
            break;

        default:
            /* You won't actually get here. */
            break;
        }

        opt = getopt_long (argc, argv, optString, longOpts, &longIndex);
    }
}

int parse_config_file (void *user, const char *section, const char *name, const char *value) {
    struct user_options *pconfig = (struct user_options *)user;

    if (MATCH (MAIN_SECTION, "cg_tolerance")) {
        pconfig->cg_tol = atof (value);
        pconfig->cg_tol_was_set = true;
    } else if (MATCH (MAIN_SECTION, "max_cg_its")) {
        pconfig->max_its = atoi (value);
        pconfig->max_its_was_set = true;

    } else if (MATCH (MAIN_SECTION, "use_preconditioner")) {
        if (strcmp (value, "true") == 0 || strcmp (value, "yes") == 0) {
            pconfig->use_jacobi = true;
        } else {
            pconfig->use_jacobi = false;
        }
        pconfig->use_jacobi_was_set = true;
    } else if (MATCH (MAIN_SECTION, "num_threads")) {
        pconfig->num_threads = atoi (value);
        pconfig->num_threads_was_set = true;
    } else if (MATCH (MAIN_SECTION, "dt_edp")) {
        pconfig->dt_edp = atof (value);
        pconfig->dt_edp_was_set = true;
    } else if (MATCH (MAIN_SECTION, "simulation_time")) {
        pconfig->final_time = atof (value);
        pconfig->final_time_was_set = true;
    } else if (MATCH (MAIN_SECTION, "abort_on_no_activity")) {
        if (strcmp (value, "true") == 0 || strcmp (value, "yes") == 0) {
            pconfig->abort_no_activity = true;
        } else {
            pconfig->abort_no_activity = false;
        }

        pconfig->abort_no_activity_was_set = true;
    } else if (MATCH (MAIN_SECTION, "use_adaptivity")) {
        if (strcmp (value, "true") == 0 || strcmp (value, "yes") == 0) {
            pconfig->adaptive = true;
        } else {
            pconfig->adaptive = false;
        }
        pconfig->adaptive_was_set = true;
    } else if (MATCH (MAIN_SECTION, "number_of_stims")) {
        // TODO: we have to decide how to handle the stims
    } else if (MATCH (MAIN_SECTION, "print_rate")) {
        pconfig->print_rate = atoi (value);
        pconfig->print_rate_was_set = true;
    } else if (MATCH (MAIN_SECTION, "output_dir")) {
        pconfig->out_dir_name = strdup (value);
        pconfig->out_dir_name_was_set = true;
    }

    else if (MATCH (ALG_SECTION, "start_discretization")) {
        pconfig->start_h = atof (value);
        pconfig->start_h_was_set = true;
    } else if (MATCH (ALG_SECTION, "maximum_discretization")) {
        pconfig->max_h = atof (value);
        pconfig->max_h_was_set = true;
    } else if (MATCH (ALG_SECTION, "refinement_bound")) {
        pconfig->ref_bound = atof (value);
        pconfig->ref_bound_was_set = true;
    } else if (MATCH (ALG_SECTION, "derefinement_bound")) {
        pconfig->deref_bound = atof (value);
        pconfig->deref_bound_was_set = true;
    } else if (MATCH (ALG_SECTION, "refine_each")) {
        pconfig->refine_each = atoi (value);
        pconfig->refine_each_was_set = true;
    } else if (MATCH (ALG_SECTION, "derefine_each")) {
        pconfig->derefine_each = atoi (value);
        pconfig->derefine_each_was_set = true;
    }

    else if (MATCH (ODE_SECTION, "dt_edo")) {
        pconfig->dt_edo = atof (value);
        pconfig->dt_edo_was_set = true;
    } else if (MATCH (ODE_SECTION, "use_gpu")) {
        if (strcmp (value, "true") == 0 || strcmp (value, "yes") == 0) {
            pconfig->gpu = true;
        } else {
            pconfig->gpu = false;
        }
        pconfig->gpu_was_set = true;
    } else if (MATCH (ODE_SECTION, "gpu_id")) {
        pconfig->gpu_id = atoi (value);
        pconfig->gpu_id_was_set = true;
    } else if (MATCH (ODE_SECTION, "model_file_path")) {
        pconfig->model_file_path = strdup (value);
        pconfig->model_file_path_was_set = true;
    }

    else {
        return 0; /* unknown section/name, error */
    }
    return 1;
}