//
// Created by sachetto on 13/10/17.
//

#include <unitypes.h>
#include <stdbool.h>
#include <stdlib.h>

#include "../utils/erros_helpers.h"
#include "../utils/logfile_utils.h"
#include "../utils/utils.h"
#include "../main/constants.h"
#include "../hash/string_hash.h"
#include "../alg/grid/grid.h"
#include "../main/config/stim_config.h"
#include "../libraries_common/config_helpers.h"

SET_SPATIAL_STIM(set_benchmark_spatial_stim) {

    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    bool stim;
    real stim_current = config->stim_current;
    real stim_value;

    if(config->spatial_stim_currents) {
        free(config->spatial_stim_currents);
    }

    config->spatial_stim_currents = (real *)malloc(n_active*sizeof(real));

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

    bool stim;
    real stim_current = config->stim_current;
    real stim_value;

    double x_limit = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(double, x_limit, config->config_data.config, "x_limit");

    if(config->spatial_stim_currents) {
        free(config->spatial_stim_currents);
    }

    config->spatial_stim_currents = (real *)malloc(n_active*sizeof(real));
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

    char *stim_file = NULL;

    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(stim_file, config->config_data.config, "stim_file");

    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;
    int s_size;

    bool stim;
    real stim_current = config->stim_current;
    real stim_value;

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

    if(config->spatial_stim_currents) {
        free(config->spatial_stim_currents);
    }

    config->spatial_stim_currents = (real *)malloc(n_active*sizeof(real));

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

    double x_limit = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(double, x_limit, config->config_data.config, "x_limit");

    bool stim;
    real stim_current = config->stim_current;
    real stim_value;

    if(config->spatial_stim_currents) {
        free(config->spatial_stim_currents);
    }

    config->spatial_stim_currents = (real *)malloc(n_active*sizeof(real));

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