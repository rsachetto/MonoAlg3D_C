#include "alg/grid/grid.h"
#include "ini_parser/ini.h"
#include "monodomain/monodomain_solver.h"
#include "monodomain/ode_solver.h"
#include "string/sds.h"
#include "utils/file_utils.h"
#include <mpi.h>
#include <string.h>


#include "single_file_libraries/stb_ds.h"

struct changed_parameters {
    char *section;
    char *name;
    union {
        f32_array values;
        float value;
    };
};

struct simulation {
    uint32_t run_number;
    struct changed_parameters *parameters;
};

static void configure_new_parameters(struct changed_parameters *changed, struct user_options *options) {

    int len = arrlen(changed);

    for(int n = 0; n < len; n++) {
        sds char_value = sdscatprintf(sdsempty(), "%lf", changed[n].value);

        if (strcmp("extra_data", changed[n].section) == 0) {
            if (options->extra_data_config) {
                shput(options->extra_data_config->config_data.config, changed[n].name,
                      strdup(char_value));
            }
        } else if (strcmp("domain", changed[n].section) == 0) {

            if (options->domain_config) {

                if (strcmp("start_dx", changed[n].name) == 0) {
                    options->domain_config->start_dx = changed[n].value;
                } else if (strcmp("start_dy", changed[n].name) == 0) {
                    options->domain_config->start_dy = changed[n].value;
                } else if (strcmp("start_dz", changed[n].name) == 0) {
                    options->domain_config->start_dz = changed[n].value;
                } else if (strcmp("maximum_dx", changed[n].name) == 0) {
                    options->domain_config->max_dx = changed[n].value;
                } else if (strcmp("maximum_dy", changed[n].name) == 0) {
                    options->domain_config->max_dy = changed[n].value;
                } else if (strcmp("maximum_dz", changed[n].name) == 0) {
                    options->domain_config->max_dz = changed[n].value;
                } else {
                    shput(options->domain_config->config_data.config, changed[n].name,
                          strdup(char_value));
                }
            }
        } else if (strcmp("assembly_matrix", changed[n].section) == 0) {
            if (options->assembly_matrix_config) {
                shput(options->assembly_matrix_config->config_data.config,
                      changed[n].name, strdup(char_value));
            }
        } else if (strcmp("linear_system_solver", changed[n].section) == 0) {
            if (options->linear_system_solver_config) {
                shput(options->linear_system_solver_config->config_data.config,
                      changed[n].name, strdup(char_value));
            }
        }


        sdsfree(char_value);
    }

}

static void free_current_simulation_resources(struct monodomain_solver *monodomain_solver,  struct ode_solver *ode_solver, struct grid *the_grid) {
    clean_and_free_grid(the_grid);
    free_ode_solver(ode_solver);
    free(monodomain_solver);
    close_logfile();
}

static struct changed_parameters parse_range_or_list_values(char *directive_rhs, char *directive_lhs) {

    int count_sec_name;
    int count_value_directive;

    sds *section_name = sdssplit(directive_rhs, "|", &count_sec_name);

    if(count_sec_name != 2 ) {
        fprintf(stderr, "Error parsing section|name directive %s\n", directive_rhs);
        MPI_Abort(MPI_COMM_WORLD, 0);
    }

    sds *value_directive = sdssplit(directive_lhs, "|", &count_value_directive);

    f32_array values = NULL;

    if(strcmp(value_directive[0], "range") == 0) {

        if(count_value_directive != 4 ) {
            fprintf(stderr, "Error parsing range directive. Valid range is range|start|end|increment, given %s\n", directive_lhs);
            MPI_Abort(MPI_COMM_WORLD, 0);
        }


        double value_start = strtod(value_directive[1], NULL);
        double value_end = strtod(value_directive[2], NULL);
        double inc = strtod(value_directive[3], NULL);

        while (value_start <= value_end) {
            arrpush(values, value_start);
            value_start += inc;
        }
    }
    else if(strcmp(value_directive[0], "list") == 0) {
        for(int i = 1; i < count_value_directive; i++) {
            float val = strtof(value_directive[i], NULL);
                    arrpush(values, val);
        }
    }
    else {
        fprintf(stderr, "Error parsing (list or range)|start|end|increment directive %s\n", directive_lhs);
        MPI_Abort(MPI_COMM_WORLD, 0);
    }


    struct changed_parameters c;
    c.section = strdup(section_name[0]);
    c.name = strdup(section_name[1]);
    c.values = values;

    sdsfreesplitres(section_name, count_sec_name);
    sdsfreesplitres(value_directive, count_value_directive);

    return c;
}

static f32_array get_combination(ui32_array counters, f32_array *sets);
static bool increment(ui32_array counters, f32_array *sets);

static f32_array * get_combinations(f32_array *sets) {

    f32_array *combinations = NULL;
    ui32_array counters = NULL;

    for(int i = 0; i < arrlen(sets); i++) {
        arrput(counters, 0);
    }

    do {
        arrput(combinations, get_combination(counters, sets));
    } while(increment(counters, sets));

            arrfree(counters);

    return combinations;
}

static f32_array get_combination(ui32_array counters, f32_array *sets) {

    f32_array o = NULL;
    for(int i = 0; i < arrlen(counters); i++) {
        arrput(o, sets[i][counters[i]]);
    }
    return o;
}

static bool increment(ui32_array counters, f32_array *sets) {
    for(int i = arrlen(counters) - 1; i >= 0; i--) {
        if(counters[i] < arrlen(sets[i]) - 1) {
            counters[i]++;
            return true;
        } else {
            counters[i] = 0;
        }
    }
    return false;
}

struct simulation * generate_all_simulations(struct string_hash_entry* modify_directives, int num_sims) {

    if(!modify_directives) return NULL;

    int n = (int)shlen(modify_directives);

    struct simulation *all_simulations = NULL;
    f32_array *sets = NULL;
    struct changed_parameters *section_names = NULL;

    for(int i = 0; i < n; i++) {
        struct changed_parameters c = parse_range_or_list_values(modify_directives[i].key, modify_directives[i].value);
        arrput(section_names, c);
        arrput(sets, c.values);
    }

    f32_array *a = get_combinations(sets);
    uint32_t sim_count = 1;

    for(int s = 0; s < num_sims; s++) {

        for(int i = 0; i < arrlen(a); i++) {

            struct simulation s;
            s.run_number = sim_count;
            s.parameters = NULL;

            for(int k = 0; k < arrlen(a[i]); k++) {
                struct changed_parameters comb;
                comb.section = strdup(section_names[k].section);
                comb.name = strdup(section_names[k].name);
                comb.value = a[i][k];
                arrput(s.parameters, comb);
            }

            arrput(all_simulations, s);
            sim_count++;
        }
    }

    for(int i = 0; i < arrlen(sets); i++) {
        arrfree(sets[i]);
    }

    arrfree(sets);
    arrfree(a);
    for(int i = 0; i < arrlen(section_names); i++) {
        free(section_names[i].section);
        free(section_names[i].name);
    }
    arrfree(section_names);

    return all_simulations;


}

void print_simulations( struct simulation *all_simulations) {

    if(!all_simulations) return;

    int n = arrlen(all_simulations);
    for (int i = 0; i < n; i++) {
        struct changed_parameters *c = all_simulations[i].parameters;
        if(!c) continue;

        int k = arrlen(c);

        printf("----------------SIMULATION %d--------------\n", all_simulations[i].run_number);

        for(int j = 0; j < k; j++) {
            printf("SECTION: %s | NAME: %s | VALUE %f\n", c[j].section, c[j].name, c[j].value);
        }

        printf("-------------------------------------------\n");

    }

}

int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);

    int rank, num_proccess, num_max_proc;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proccess);

    int simulation_number_start;
    int num_simulations = 0;

    char *output_folder;

    bool skip_existing;

    struct batch_options *batch_options;
    batch_options = new_batch_options();

    parse_batch_options(argc, argv, batch_options);

    if(ini_parse(batch_options->batch_config_file, parse_batch_config_file, batch_options) < 0) {
        fprintf(stderr, "Error: Can't load the config file %s\n", batch_options->batch_config_file);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    skip_existing = batch_options->skip_existing_run;

    if(rank == 0) {
        create_dir(batch_options->output_folder);
    }

    struct simulation *all_simulations = generate_all_simulations(batch_options->config_to_change, batch_options->num_simulations);
    //print_simulations(all_simulations);

    int total_simulations = arrlen(all_simulations);

    if(num_proccess > total_simulations) {
        num_max_proc = total_simulations;
    } else {
        num_max_proc = num_proccess;
    }

    if(rank < num_max_proc)
        num_simulations = total_simulations / num_max_proc;

    int last_rank_extra = total_simulations % num_max_proc;

    simulation_number_start = rank*num_simulations;

    if(rank == num_max_proc - 1) {
        num_simulations += last_rank_extra;
    }

    if(num_simulations == 0) {
        MPI_Finalize();
        return EXIT_SUCCESS;
    }

    struct user_options *options;
    options = new_user_options();

    struct grid *the_grid;
    struct monodomain_solver *monodomain_solver;
    struct ode_solver *ode_solver;

    if(ini_parse(batch_options->initial_config, parse_config_file, options) < 0) {
        fprintf(stderr, "Error parsing config file %s\n", batch_options->initial_config);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if(!options->save_mesh_config) {
        options->save_mesh_config->out_dir_name = strdup("batch_run");
    }

    char *initial_out_dir_name  = strdup(options->save_mesh_config->out_dir_name);

    no_stdout = options->quiet;

    output_folder = batch_options->output_folder;
    options->draw = false;

    for(int s = simulation_number_start; s < simulation_number_start+num_simulations; s++) {

            the_grid = new_grid();
            monodomain_solver = new_monodomain_solver();
            ode_solver = new_ode_solver();

            struct simulation simulation = all_simulations[s];
            int config_n = arrlen(simulation.parameters);

            sds new_out_dir_name = sdscatprintf(sdsempty(), "%s/%s_run_%d", output_folder, initial_out_dir_name, simulation.run_number);

            for(int n = 0; n < config_n; n++) {
                sdscatprintf(new_out_dir_name, "_%s_%f", simulation.parameters[n].name, simulation.parameters[n].value);
            }

            if(skip_existing) {
                if(check_simulation_completed(new_out_dir_name)) {
                    printf("Rank %d skipping existing simulation on %s\n", rank, new_out_dir_name);
                    sdsfree(new_out_dir_name);
                    continue;
                }
            }

            configure_new_parameters(simulation.parameters, options);

            free(options->save_mesh_config->out_dir_name);
            options->save_mesh_config->out_dir_name = strdup(new_out_dir_name);
            sdsfree(new_out_dir_name);

            printf("Rank %d, performing simulation %d and saving in %s\n", rank, simulation_number_start, options->save_mesh_config->out_dir_name);

            // Create the output dir and the logfile
            sds buffer_log = sdsnew("");
            sds buffer_ini = sdsnew("");

            remove_directory(options->save_mesh_config->out_dir_name);
            create_dir(options->save_mesh_config->out_dir_name);

            buffer_log = sdscatfmt(buffer_log, "%s/outputlog.txt", options->save_mesh_config->out_dir_name);
            open_logfile(buffer_log);

            print_to_stdout_and_file("Command line to reproduce this simulation:\n");
            for (int i = 0; i < argc; i++) {
                print_to_stdout_and_file("%s ", argv[i]);
            }

            print_to_stdout_and_file("\n");

            buffer_ini =
                    sdscatfmt(buffer_ini, "%s/original_configuration.ini", options->save_mesh_config->out_dir_name);

            print_to_stdout_and_file("For reproducibility purposes the configuration file was copied to file: %s\n",
                                     buffer_ini);

            options_to_ini_file(options, buffer_ini);

            sdsfree(buffer_log);
            sdsfree(buffer_ini);

            configure_ode_solver_from_options(ode_solver, options);
            configure_monodomain_solver_from_options(monodomain_solver, options);
            configure_grid_from_options(the_grid, options);

#ifndef COMPILE_CUDA
            if(ode_solver->gpu) {
                print_to_stdout_and_file("Cuda runtime not found in this system. Fallbacking to CPU solver!!\n");
                ode_solver->gpu = false;
            }
#endif

            int nt = monodomain_solver->num_threads;

            if (nt == 0)
                nt = 1;

#if defined(_OPENMP)
            omp_set_num_threads(nt);
#endif
            solve_monodomain(monodomain_solver, ode_solver, the_grid, options);
            free_current_simulation_resources(monodomain_solver, ode_solver, the_grid);

    }

    free_batch_options(batch_options);
    free_user_options(options);
    free(initial_out_dir_name);

    MPI_Finalize();

    return EXIT_SUCCESS;
}
