#include "../3dparty/ini_parser/ini_file_sections.h"
#include "../3dparty/sds/sds.h"
#include "../3dparty/stb_ds.h"
#include "../config/config_parser.h"
#include "batch_utils.h"
#include <stdio.h>
#include <string.h>

void configure_new_parameters(struct changed_parameters *changed, struct user_options *options) {

    int len = arrlen(changed);

    for (int n = 0; n < len; n++) {
        float float_value = strtof(changed[n].value, NULL);

        if (strcmp(MAIN_SECTION, changed[n].section) == 0) {

            if (strcmp("num_threads", changed[n].name) == 0) {
                options->num_threads = (int) float_value;
            } else if (strcmp("dt_pde", changed[n].name) == 0) {
                options->dt_pde = float_value;
            } else if (strcmp("simulation_time", changed[n].name) == 0) {
                options->final_time = float_value;
            } else if (strcmp("vm_threshold", changed[n].name) == 0) {
                options->vm_threshold = float_value;
            } else if (strcmp("use_adaptivity", changed[n].name) == 0) {
                options->adaptive = IS_TRUE(changed[n].value);
            } else if (strcmp("abort_on_no_activity", changed[n].name) == 0) {
                options->abort_no_activity = IS_TRUE(changed[n].value);
            } else if (strcmp("gpu_id", changed[n].name) == 0) {
                options->gpu_id = (int) float_value;
            } else {
                fprintf(stderr, "Parameter %s in section %s cannot be changed\n", changed[n].name, changed[n].section);
            }

        } else if (strcmp(ODE_SECTION, changed[n].section) == 0) {
            if (strcmp("dt", changed[n].name) == 0) {
                options->dt_ode = float_value;
                //TODO: this is deprecated
            } else if (strcmp("gpu_id", changed[n].name) == 0) {
                options->gpu_id = (int) float_value;
            } else if (strcmp("library_file", changed[n].name) == 0) {
                free(options->model_file_path);
                options->model_file_path = strdup(changed[n].value);
            }
        } else if (strcmp(ODE_PURKINJE_SECTION, changed[n].section) == 0) {
            if (strcmp("dt", changed[n].name) == 0) {
                options->dt_ode = float_value;
                //TODO: this is deprecated
            } else if (strcmp("gpu_id", changed[n].name) == 0) {
                options->gpu_id = (int) float_value;
            } else if (strcmp("library_file", changed[n].name) == 0) {
                free(options->model_file_path);
                options->model_file_path = strdup(changed[n].value);
            }
        }

        if (strcmp(EXTRA_DATA_SECTION, changed[n].section) == 0) {
            if (options->extra_data_config) {
                shput(options->extra_data_config->config_data, changed[n].name,
                        strdup(changed[n].value));
            }
        } else if (strcmp(DOMAIN_SECTION, changed[n].section) == 0) {

            if (options->domain_config) {
                set_or_overwrite_common_data(options->domain_config, changed[n].name, changed[n].value, NULL, NULL);
            }
        } else if (strcmp(PURKINJE_SECTION, changed[n].section) == 0) {

            if (options->purkinje_config) {
                set_or_overwrite_common_data(options->purkinje_config, changed[n].name, changed[n].value, NULL, NULL);
            }
        }
          else if (strcmp(MATRIX_ASSEMBLY_SECTION, changed[n].section) == 0) {
            if (options->assembly_matrix_config) {
                set_or_overwrite_common_data(options->assembly_matrix_config, changed[n].name, changed[n].value, NULL, NULL);
            }

        } else if (strcmp(LINEAR_SYSTEM_SOLVER_SECTION, changed[n].section) == 0) {
            if (options->linear_system_solver_config) {
                set_or_overwrite_common_data(options->linear_system_solver_config, changed[n].name, changed[n].value, NULL, NULL);

            }
        } else if (strcmp(LINEAR_SYSTEM_SOLVER_PURKINJE_SECTION, changed[n].section) == 0) {
            if (options->purkinje_linear_system_solver_config) {
                set_or_overwrite_common_data(options->purkinje_linear_system_solver_config, changed[n].name, changed[n].value, NULL, NULL);

            }
        } else if (strncmp(changed[n].section, STIM_SECTION, strlen(STIM_SECTION)) == 0) {

            struct config *sc = (struct config *) shget(options->stim_configs, changed[n].section);
            if (sc) {
                set_or_overwrite_common_data(sc, changed[n].name, changed[n].value, NULL, NULL);
                shput(options->stim_configs, changed[n].section, sc);

            }

        } else if (strncmp(changed[n].section, STIM_PURKINJE_SECTION, strlen(STIM_PURKINJE_SECTION)) == 0) {

            struct config *sc = (struct config *) shget(options->stim_configs, changed[n].section);
            if (sc) {
                set_or_overwrite_common_data(sc, changed[n].name, changed[n].value, NULL, NULL);
                shput(options->stim_configs, changed[n].section, sc);

            }

        }

    }

}

void free_current_simulation_resources(struct monodomain_solver *monodomain_solver, struct ode_solver *ode_solver,
        struct grid *the_grid) {
    clean_and_free_grid(the_grid);
    free_ode_solver(ode_solver);
    free(monodomain_solver);
    close_logfile();
}

struct changed_parameters parse_range_or_list_values(char *directive_rhs, char *directive_lhs) {

    int count_sec_name;
    int count_value_directive;

    sds *section_name = sdssplit(directive_rhs, "|", &count_sec_name);

    if (count_sec_name != 2) {
        fprintf(stderr, "Error parsing section|name directive %s\n", directive_rhs);
        return (struct changed_parameters){0};
    }

    sds *value_directive = sdssplit(directive_lhs, "|", &count_value_directive);

    string_array values = NULL;

    if (strcmp(value_directive[0], "range") == 0) {

        if (count_value_directive != 4) {
            fprintf(stderr, "Error parsing range directive. Valid range is range|start|end|increment, given %s\n",
                    directive_lhs);
            return (struct changed_parameters){0};
        }


        double value_start = strtod(value_directive[1], NULL);
        double value_end = strtod(value_directive[2], NULL);
        double inc = strtod(value_directive[3], NULL);

        while (value_start <= value_end) {
            char tmp[256];
            sprintf(tmp, "%f", value_start);
            arrpush(values, strdup(tmp));
            value_start += inc;
        }
    } else if (strcmp(value_directive[0], "list") == 0) {
        for (int i = 1; i < count_value_directive; i++) {
            arrpush(values, strdup(value_directive[i]));
        }
    } else {
        fprintf(stderr, "Error parsing (list or range)|start|end|increment directive %s\n", directive_lhs);
        return (struct changed_parameters){0};
    }

    struct changed_parameters c;
    c.section = strdup(section_name[0]);
    c.name = strdup(section_name[1]);
    c.values = values;

    sdsfreesplitres(section_name, count_sec_name);
    sdsfreesplitres(value_directive, count_value_directive);

    return c;
}

static bool increment(ui32_array counters, string_array *sets) {
    for (int i = arrlen(counters) - 1; i >= 0; i--) {
        if (counters[i] < arrlen(sets[i]) - 1) {
            counters[i]++;
            return true;
        } else {
            counters[i] = 0;
        }
    }
    return false;
}

static string_array get_combination(const uint32_t *counters, string_array *sets) {

    string_array o = NULL;
    for (int i = 0; i < arrlen(counters); i++) {
        arrput(o, strdup(sets[i][counters[i]]));
    }
    return o;
}

static string_array *get_combinations(string_array *sets) {

    string_array *combinations = NULL;
    ui32_array counters = NULL;

    for (int i = 0; i < arrlen(sets); i++) {
        arrput(counters, 0);
    }

    do {
        arrput(combinations, get_combination(counters, sets));
    } while (increment(counters, sets));

    arrfree(counters);

    return combinations;
}

struct simulation *generate_all_simulations(struct string_hash_entry *modify_directives, int num_sims) {

    if (!modify_directives) return NULL;

    int n = (int) shlen(modify_directives);

    struct simulation *all_simulations = NULL;
    string_array *sets = NULL;
    struct changed_parameters *section_names = NULL;

    for (int i = 0; i < n; i++) {
        struct changed_parameters c = parse_range_or_list_values(modify_directives[i].key, modify_directives[i].value);
        if(c.name == NULL) {
            return NULL;
        }
        arrput(section_names, c);
        arrput(sets, c.values);
    }

    string_array *combinations = get_combinations(sets);
    uint32_t sim_count = 1;

    long num_combinations = arrlen(combinations);

    for (int s = 0; s < num_sims; s++) {

        for (long i = 0; i < num_combinations; i++) {

            struct simulation sim;
            sim.run_number = sim_count;
            sim.parameters = NULL;

            long conbinations_i_size = arrlen(combinations[i]);

            for (long k = 0; k < conbinations_i_size; k++) {
                struct changed_parameters comb;
                comb.section = strdup(section_names[k].section);
                comb.name = strdup(section_names[k].name);
                comb.value = combinations[i][k];
                arrput(sim.parameters, comb);
            }

            arrput(all_simulations, sim);
            sim_count++;
        }
    }

    for (long i = 0; i < arrlen(sets); i++) {
        arrfree(sets[i]);
    }

    arrfree(sets);
    arrfree(combinations);
    for (long i = 0; i < arrlen(section_names); i++) {
        free(section_names[i].section);
        free(section_names[i].name);
    }

    arrfree(section_names);

    return all_simulations;

}

void print_simulations(struct simulation *all_simulations) {

    if (!all_simulations) return;

    int n = arrlen(all_simulations);
    for (int i = 0; i < n; i++) {
        struct changed_parameters *c = all_simulations[i].parameters;
        if (!c) continue;

        int k = arrlen(c);

        printf("----------------SIMULATION %u--------------\n", all_simulations[i].run_number);

        for (int j = 0; j < k; j++) {
            printf("SECTION: %s | NAME: %s | VALUE: %s\n", c[j].section, c[j].name, c[j].value);
        }

        printf("-------------------------------------------\n");

    }

}
