#include <criterion/alloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include "3dparty/stb_ds.h"
#include "common_types/common_types.h"
#include "utils/file_utils.h"
#include "config/config_common.h"
#include "config/domain_config.h"
#include "domains_library/mesh_info_data.h"

static const char *expansion_opt_string = "i:o:n:v:h?";
static const struct option long_expansion_options[] = {{"input", required_argument, NULL, 'i'},
                                                        {"output", required_argument, NULL, 'o'},
                                                        {"n_rows", required_argument, NULL, 'n'},
                                                        {"n_volumes", required_argument, NULL, 'v'}};
struct expansion_options {
    char *input;
    char *output;
    int n_rows;
    char* n_volumes;
};

static void display_expansion_usage(char **argv) {

    printf("Usage: %s [options]\n\n", argv[0]);
    printf("Options:\n");
    printf("--input  | -i [input]. File. Default NULL.\n");
    printf("--output | -o [output]. Output file. Default NULL.\n");
    printf("--n_rows | -n [n_rows]. Number of rows to expand.\n");
    printf("--n_volumes | -v [n_volumes]. Number of volumes in the mesh\n");
    printf("--help | -h. Shows this help and exit \n");
    exit(EXIT_FAILURE);
}


static void parse_expansion_options(int argc, char **argv, struct expansion_options *user_args) {

    int opt = 0;
    int option_index;

    opt = getopt_long_only(argc, argv, expansion_opt_string, long_expansion_options, &option_index);

    while(opt != -1) {
        switch(opt) {
        case 'i':
            user_args->input = strdup(optarg);
            break;
        case 'o':
            user_args->output = strdup(optarg);
            break;
        case 'n':
            user_args->n_rows = atoi(optarg);
            break;
        case 'v':
            user_args->n_volumes = strdup(optarg);
            break;
        case 'h': /* fall-through is intentional */
        case '?':
            display_expansion_usage(argv);
            break;
        default:
            /* You won't actually get here. */
            break;
        }

        opt = getopt_long(argc, argv, expansion_opt_string, long_expansion_options, &option_index);
    }

    if(!user_args->input) {
        display_expansion_usage(argv);
    }
}

static void expand_scar(const char *input, const char *output, int n_rows, char* n_volumes) {

    //set_no_stdout(true);

    struct grid *grid = new_grid();
    struct config *domain_config;

    domain_config = alloc_and_init_config_data();
    char *discretization = "200";

    shput_dup_value(domain_config->config_data, "start_discretization", strdup(discretization));
    shput_dup_value(domain_config->config_data, "maximum_discretization", strdup(discretization));
    shput_dup_value(domain_config->config_data, "mesh_file", strdup(input));

    domain_config->main_function_name = strdup("initialize_grid_with_custom_mesh");
    shput_dup_value(domain_config->config_data, "name", "Test custom mesh");

    shput(domain_config->config_data, "side_length_x", strdup(discretization));
    shput(domain_config->config_data, "side_length_y", strdup(discretization));
    shput(domain_config->config_data, "side_length_z", strdup(discretization));
    shput(domain_config->config_data, "number_of_points", strdup(n_volumes));
    shput(domain_config->config_data, "num_extra_fields", "5");

    init_config_functions(domain_config, "./shared_libs/libdefault_domains.so", "domain");

    int success = ((set_spatial_domain_fn*)domain_config->main_function)(domain_config, grid);

    if(success == 0) {
        printf("Error loading mesh in %s. Exiting!\n", input);
        exit(EXIT_FAILURE);
    }

    struct cell_node *neighbour;
    real *extra_info_n;
    real *extra_info;

    bool evaluate;

    for (int i = 0; i < n_rows; i++) {

        FOR_EACH_CELL(grid) {
            cell->visited = false;
        }

        FOR_EACH_CELL(grid) {

            extra_info = (real*)cell->mesh_extra_info;

            if(i == 0) {
                evaluate = cell->active && extra_info[4] == 1;
            } else {
                evaluate = cell->active && extra_info[4] == 2 && !cell->visited;
            }

            if(evaluate) {

                cell->visited = true;

                neighbour = get_cell_neighbour(cell, cell->neighbours[FRONT]);

                if(neighbour) {
                    extra_info_n = (real*)neighbour->mesh_extra_info;
                    if(extra_info_n[4] != 1 && extra_info_n[4] != 2) {
                        extra_info_n[4] = 2;
                        neighbour->visited = true;
                    }
                }

                neighbour = get_cell_neighbour(cell, cell->neighbours[BACK]);

                if(neighbour) {
                    extra_info_n = (real*)neighbour->mesh_extra_info;
                    if(extra_info_n[4] != 1 && extra_info_n[4] != 2) {
                        extra_info_n[4] = 2;
                        neighbour->visited = true;
                    }
                }

                neighbour = get_cell_neighbour(cell, cell->neighbours[DOWN]);
                if(neighbour) {
                    extra_info_n = (real*)neighbour->mesh_extra_info;
                    if(extra_info_n[4] != 1 && extra_info_n[4] != 2) {
                        extra_info_n[4] = 2;
                        neighbour->visited = true;
                    }
                }


                neighbour = get_cell_neighbour(cell, cell->neighbours[TOP]);
                if(neighbour) {
                    extra_info_n = (real*)neighbour->mesh_extra_info;
                    if(extra_info_n[4] != 1 && extra_info_n[4] != 2) {
                        extra_info_n[4] = 2;
                        neighbour->visited = true;
                    }
                }

                neighbour = get_cell_neighbour(cell, cell->neighbours[RIGHT]);
                if(neighbour) {
                    extra_info_n = (real*)neighbour->mesh_extra_info;
                    if(extra_info_n[4] != 1 && extra_info_n[4] != 2) {
                        extra_info_n[4] = 2;
                        neighbour->visited = true;
                    }
                }

                neighbour = get_cell_neighbour(cell, cell->neighbours[LEFT]);
                if(neighbour) {
                    extra_info_n = (real*)neighbour->mesh_extra_info;
                    if(extra_info_n[4] != 1 && extra_info_n[4] != 2) {
                        extra_info_n[4] = 2;
                        neighbour->visited = true;
                    }
                }
            }
        }
    }

    FILE *out = fopen(output, "w");
    real_cpu dx, dy, dz;

    FOR_EACH_CELL(grid) {
        if(cell->active) {

            extra_info = (real*)cell->mesh_extra_info;

            dx = cell->discretization.x / 2.0;
            dy = cell->discretization.y / 2.0;
            dz = cell->discretization.z / 2.0;

            fprintf(out, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d\n", cell->center.x, cell->center.y, cell->center.z, dx, dy, dz, extra_info[0], extra_info[1], extra_info[2], extra_info[3], (int)extra_info[4]);
        }
    }


}

int main(int argc, char **argv) {

    struct expansion_options *options  = CALLOC_ONE_TYPE(struct expansion_options);

    parse_expansion_options(argc, argv, options);

    char *input = options->input;
    char *output = options->output;

    struct path_information input_info;

    get_path_information(input, &input_info);

    if(!input_info.exists) {
        fprintf(stderr, "Invalid file (%s)! The input parameter should be an existing alg file!\n", input);
        return EXIT_FAILURE;
    }

    if(!output) {
        output = "expanded_mesh.alg";
    }

    expand_scar(input, output, options->n_rows, options->n_volumes);


    return EXIT_SUCCESS;
}
