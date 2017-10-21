//
// Created by sachetto on 13/10/17.
//

#include <unitypes.h>
#include <stdbool.h>
#include "../utils/erros_helpers.h"
#include "../utils/logfile_utils.h"
#include "../utils/utils.h"
#include "../main/constants.h"
#include "../hash/string_hash.h"
#include "stdlib.h"
#include "../alg/grid/grid.h"
#include "../main/config/stim_config.h"

SET_SPATIAL_STIM(set_benchmark_spatial_stim) {

    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    bool stim;
    Real stim_current = config->stim_current;
    Real stim_value;

#pragma omp parallel for private(stim, stim_value)
    for (int i = 0; i < n_active; i++) {

        stim = ac[i]->center_x > 5500.0;
        stim &= ac[i]->center_x < 7000.0;
        stim &= ac[i]->center_y < 1500.0;
        stim &= ac[i]->center_z < 1500.0;

        if (stim) {
            stim_value = stim_current;
        } else {
            stim_value = 0.0;
        }

        config->spatial_stim_currents[i] = stim_value;
    }
}

SET_SPATIAL_STIM(stim_if_x_less_than) {

    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    char *config_char = string_hash_search(config->config_data.config, "x_limit");
    if(config_char == NULL) {
        report_parameter_error_on_function("stim_if_x_less_than", "x_limit");
    }
    double x_limit = atof(config_char);
    free(config_char);

    bool stim;
    Real stim_current = config->stim_current;
    Real stim_value;

    #pragma omp parallel for private(stim, stim_value)
    for (int i = 0; i < n_active; i++) {
        stim = ac[i]->center_x < x_limit;

        if (stim) {
            stim_value = stim_current;
        } else {
            stim_value = 0.0;
        }

        config->spatial_stim_currents[i] = stim_value;

    }
}

SET_SPATIAL_STIM(set_stim_from_file) {

    char *stim_file = string_hash_search(config->config_data.config, "stim_file");
    if(stim_file == NULL) {
        report_parameter_error_on_function("set_stim_from_file", "stim_file");
    }

    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;
    int s_size;

    bool stim;
    Real stim_current = config->stim_current;
    Real stim_value;


    FILE *s_file = fopen(stim_file,"r");

    if(!s_file) {
        fprintf(stderr, "Error opening stim file %s! Exiting!\n", stim_file);
        exit(EXIT_FAILURE);
    }

    fscanf(s_file, "%d\n", &s_size);

    double **cell_stims = (double**) malloc(sizeof(double*)*s_size);
    for(int i=0; i< s_size; i++){
        cell_stims[i] = (double*) malloc(sizeof(double) * 3);
        if(cell_stims[i] == NULL) {
            fprintf(stderr, "Failed to allocate memory for the stim file\n");
            exit(0);
        }

        fscanf(s_file, "%lf %lf %lf\n",&cell_stims[i][0],&cell_stims[i][1],&cell_stims[i][2]);
    }

    sort_vector(cell_stims, s_size);

    fclose(s_file);

    #pragma omp parallel for private(stim, stim_value)
    for (int i = 0; i < n_active; i++) {

        double center_x = ac[i]->center_x;
        double center_y = ac[i]->center_y;
        double center_z = ac[i]->center_z;

        int index = inside_mesh(cell_stims, center_x, center_y, center_z, 0, s_size - 1);
        stim = (index != -1);

        if (stim) {
            stim_value = stim_current;
        } else {
            stim_value = 0.0;
        }

        config->spatial_stim_currents[i] = stim_value;

    }

}

SET_SPATIAL_STIM(stim_if_x_greater_equal_than) {

    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    char *config_char = string_hash_search(config->config_data.config, "x_limit");
    if(config_char == NULL) {
        report_parameter_error_on_function("stim_if_x_more_than", "x_limit");
    }
    double x_limit = atof(config_char);
    free(config_char);

    bool stim;
    Real stim_current = config->stim_current;
    Real stim_value;

#pragma omp parallel for private(stim)
    for (int i = 0; i < n_active; i++) {
        stim = ac[i]->center_x >= x_limit;

        if (stim) {
            stim_value = stim_current;
        } else {
            stim_value = 0.0;
        }

        config->spatial_stim_currents[i] = stim_value;

    }
}