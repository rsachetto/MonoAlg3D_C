// Author: Lucas Berg (@bergolho)
// Script used to clip a section of a mesh and extract the extra data information
// Version: 29/01/2024
// Last change: 29/01/2024

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
#include "extra_data_library/helper_functions.h"

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
    char *discretization = "300";
    //char *discretization = "500";

    shput_dup_value(domain_config->config_data, "start_discretization", strdup(discretization));
    shput_dup_value(domain_config->config_data, "maximum_discretization", strdup(discretization));
    shput_dup_value(domain_config->config_data, "mesh_file", strdup(input));

    domain_config->main_function_name = strdup("initialize_grid_with_dti_mesh");
    shput_dup_value(domain_config->config_data, "name", "Oxford DTI004 with Transmurality and Fiber orientation");

    //shput(domain_config->config_data, "side_length_x", strdup(discretization));
    //shput(domain_config->config_data, "side_length_y", strdup(discretization));
    //shput(domain_config->config_data, "side_length_z", strdup(discretization));
    shput(domain_config->config_data, "num_volumes", strdup(n_volumes));
    //shput(domain_config->config_data, "num_extra_fields", "5");

    init_config_functions(domain_config, "./shared_libs/libdefault_domains.so", "domain");

    int success = ((set_spatial_domain_fn*)domain_config->main_function)(domain_config, grid);

    if(success == 0) {
        printf("Error loading mesh in %s. Exiting!\n", input);
        exit(EXIT_FAILURE);
    }

    struct dti_mesh_info *extra_data = NULL;
    bool ignore_cell;
    real_cpu dx, dy, dz;
    real_cpu min_x, max_x, min_y, max_y, min_z, max_z;
    real *f, *s, *n;
    real_cpu transmurality, base_apex_heterogeneity, apicobasal;
    uint32_t transmurality_labels;

    min_x=76803.6;
    min_y=51176.3;
    min_z=8740.9;
    max_x=96803.6;
    max_y=71176.3;
    max_z=28740.9;

    FILE *out = fopen(output, "w");
    FOR_EACH_CELL(grid) {
        if(cell->active) {
            
            f = cell->sigma.fibers.f;
            s = cell->sigma.fibers.s;
            n = cell->sigma.fibers.n;

            dx = cell->discretization.x / 2.0;
            dy = cell->discretization.x / 2.0;
            dz = cell->discretization.x / 2.0;

            extra_data = (struct dti_mesh_info *)cell->mesh_extra_info;
            transmurality = extra_data->transmurality;
            base_apex_heterogeneity = extra_data->base_apex_heterogeneity;
            apicobasal = extra_data->apicobasal;
            transmurality_labels = extra_data->dti_transmurality_labels;
            
            // Change the FAST_ENDO tag to ENDO
            if (transmurality_labels == 3) {
                transmurality_labels = 0;
            }

            ignore_cell = cell->center.x < min_x || cell->center.x > max_x || cell->center.y < min_y || cell->center.y > max_y || cell->center.z < min_z || cell->center.z > max_z;

            if(!ignore_cell) {
                fprintf(out, "%g,%g,%g,%g,%g,%g,%g,%g,%g,%u,%g,%g,%g,%g,%g,%g,%g,%g,%g\n", cell->center.x, cell->center.y, cell->center.z, dx, dy, dz, \
                                                    transmurality, base_apex_heterogeneity, apicobasal, \
                                                    transmurality_labels, \
                                                    f[0], f[1], f[2], s[0], s[1], s[2], n[0], n[1], n[2]);
            }

        }
    }
    fclose(out);
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