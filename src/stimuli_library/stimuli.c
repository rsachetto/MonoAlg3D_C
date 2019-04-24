//
// Created by sachetto on 13/10/17.
//

#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>

#include "../utils/utils.h"
#include "../monodomain/constants.h"
#include "../alg/grid/grid.h"
#include "../config/stim_config.h"
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

	int i;

    #pragma omp parallel for private(stim, stim_value)
    for (i = 0; i < n_active; i++) {

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

    real_cpu x_limit = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, x_limit, config->config_data.config, "x_limit");

    if(config->spatial_stim_currents) {
        free(config->spatial_stim_currents);
    }

    config->spatial_stim_currents = (real *)malloc(n_active*sizeof(real));

	int i;

    #pragma omp parallel for private(stim, stim_value)
    for (i = 0; i < n_active; i++) {
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
    size_t s_size;

    bool stim;
    real stim_current = config->stim_current;
    real stim_value;

    FILE *s_file = fopen(stim_file,"r");

    if(!s_file) {
        fprintf(stderr, "Error opening stim file %s! Exiting!\n", stim_file);
        exit(EXIT_FAILURE);
    }

    fscanf(s_file, "%zu\n", &s_size);

    real_cpu **cell_stims = (real_cpu**) malloc(sizeof(real_cpu*)*s_size);
    for(int i=0; i< s_size; i++){
        cell_stims[i] = (real_cpu*) malloc(sizeof(real_cpu) * 3);
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
	
	int i;
	
	#pragma omp parallel for private(stim, stim_value)
    for (i = 0; i < n_active; i++) {

        real_cpu center_x = ac[i]->center_x;
        real_cpu center_y = ac[i]->center_y;
        real_cpu center_z = ac[i]->center_z;

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

    real_cpu x_limit = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, x_limit, config->config_data.config, "x_limit");

    real stim_current = config->stim_current;
    real stim_value;

    if(config->spatial_stim_currents) {
        free(config->spatial_stim_currents);
    }

    config->spatial_stim_currents = (real *)malloc(n_active*sizeof(real));

	int i;
	
	#pragma omp parallel for private(stim_value)
    for (i = 0; i < n_active; i++) {
        bool stim = (ac[i]->center_x >= x_limit);

        if (stim) {
            stim_value = stim_current;
        } else {
            stim_value = 0.0;
        }

        config->spatial_stim_currents[i] = stim_value;

    }
}

SET_SPATIAL_STIM(stim_base_mouse) {

    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    real_cpu stim_size = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, stim_size, config->config_data.config, "stim_size");

    real stim_current = config->stim_current;
    real stim_value;

    if(config->spatial_stim_currents) {
        free(config->spatial_stim_currents);
    }

    config->spatial_stim_currents = (real *)malloc(n_active*sizeof(real));

	int i;

    #pragma omp parallel for private(stim_value)
    for (i = 0; i < n_active; i++) {

        bool stim;
        stim  = (ac[i]->center_x >= 3000.0 - stim_size) && (ac[i]->center_x <= 3000.0 + stim_size);
        stim &= (ac[i]->center_y >= 2400.0 - stim_size) && (ac[i]->center_y <= 2400.0 + stim_size);
        stim &= (ac[i]->center_z >= 300 - stim_size) && (ac[i]->center_z <= 300 + stim_size);

        if (stim) {
            stim_value = stim_current;
        } else {
            stim_value = 0.0;
        }

        config->spatial_stim_currents[i] = stim_value;

    }
}

SET_SPATIAL_STIM(stim_mouse_spiral) {

    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    real stim_current = config->stim_current;
    real stim_value;

    if(config->spatial_stim_currents) {
        free(config->spatial_stim_currents);
    }

    config->spatial_stim_currents = (real *)malloc(n_active*sizeof(real));

	int i;

    #pragma omp parallel for private(stim_value)
    for (i = 0; i < n_active; i++) {

        bool stim;        

        stim  = (ac[i]->center_x >= 3000.0) && (ac[i]->center_x <= 6000.0);
        stim &= (ac[i]->center_y >= 1940.0) && (ac[i]->center_y <= 6100.0);
        stim &= (ac[i]->center_z >=  2230.0) && (ac[i]->center_z <= 5800.0);

        if (stim) {
            stim_value = stim_current;
        } else {
            stim_value = 0.0;
        }

        config->spatial_stim_currents[i] = stim_value;

    }
}

SET_SPATIAL_STIM(stim_x_y_limits) {

    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    real_cpu max_x = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, max_x, config->config_data.config, "max_x");

    real_cpu min_x = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, min_x, config->config_data.config, "min_x");

    real_cpu max_y = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, max_y, config->config_data.config, "max_y");

    real_cpu min_y = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, min_y, config->config_data.config, "min_y");

    real stim_current = config->stim_current;
    real stim_value;

    if(config->spatial_stim_currents) {
        free(config->spatial_stim_currents);
    }

    config->spatial_stim_currents = (real *)malloc(n_active*sizeof(real));

    int i;

#pragma omp parallel for private(stim_value)
    for (i = 0; i < n_active; i++) {

        bool stim;
        stim  = (ac[i]->center_x >= min_x) && (ac[i]->center_x <= max_x);
        stim &= (ac[i]->center_y >= min_y) && (ac[i]->center_y <= max_y);

        if (stim) {
            stim_value = stim_current;
        } else {
            stim_value = 0.0;
        }

        config->spatial_stim_currents[i] = stim_value;

    }
}

SET_SPATIAL_STIM(stim_x_y_z_limits) {

    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    real_cpu max_x = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, max_x, config->config_data.config, "max_x");

    real_cpu min_x = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, min_x, config->config_data.config, "min_x");

    real_cpu max_y = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, max_y, config->config_data.config, "max_y");

    real_cpu min_y = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, min_y, config->config_data.config, "min_y");

    real_cpu max_z = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, max_z, config->config_data.config, "max_z");

    real_cpu min_z = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, min_z, config->config_data.config, "min_z");

    real stim_current = config->stim_current;
    real stim_value;

    if(config->spatial_stim_currents) {
        free(config->spatial_stim_currents);
    }

    config->spatial_stim_currents = (real *)malloc(n_active*sizeof(real));

    int i;

    #pragma omp parallel for private(stim_value)
    for (i = 0; i < n_active; i++) {

        bool stim;
        stim  = (ac[i]->center_x >= min_x) && (ac[i]->center_x <= max_x);
        stim &= (ac[i]->center_y >= min_y) && (ac[i]->center_y <= max_y);
        stim &= (ac[i]->center_z >= min_z) && (ac[i]->center_z <= max_z);

        if (stim) {
            stim_value = stim_current;
        } else {
            stim_value = 0.0;
        }

        config->spatial_stim_currents[i] = stim_value;

    }
}

// ***********************************************************************************************************
// New Berg's stimulus
SET_SPATIAL_STIM(stim_if_inside_circle_than) {

    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    bool stim;
    real stim_current = config->stim_current;
    real stim_value;

    real_cpu center_x = 0.0;
    real_cpu center_y = 0.0;
    real_cpu center_z = 0.0;
    real_cpu radius = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, center_x, config->config_data.config, "center_x");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, center_y, config->config_data.config, "center_y");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, center_z, config->config_data.config, "center_z");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, radius, config->config_data.config, "radius");

    if(config->spatial_stim_currents) {
        free(config->spatial_stim_currents);
    }

    config->spatial_stim_currents = (real *)malloc(n_active*sizeof(real));

	int i;

    #pragma omp parallel for private(stim, stim_value)
    for (i = 0; i < n_active; i++) 
    {
        real_cpu dist = sqrt(pow(ac[i]->center_x-center_x,2)+pow(ac[i]->center_y-center_y,2)+pow(ac[i]->center_z-center_z,2));
        stim = dist <= radius;

        if (stim) {
            stim_value = stim_current;
        } else {
            stim_value = 0.0;
        }

        config->spatial_stim_currents[i] = stim_value;

    }
}

SET_SPATIAL_STIM(stim_if_id_less_than) {

    uint32_t n_active = the_grid->num_active_cells;
    //struct cell_node **ac = the_grid->active_cells;

    bool stim;
    real stim_current = config->stim_current;
    real stim_value;

    int id = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(int, id, config->config_data.config, "id_limit");

    if(config->spatial_stim_currents) 
    {
        free(config->spatial_stim_currents);
    }

    config->spatial_stim_currents = (real *)malloc(n_active*sizeof(real));

	int i;

    #pragma omp parallel for private(stim, stim_value)
    for (i = 0; i < n_active; i++) 
    {
        stim = i <= id;

        if (stim) 
        {
            stim_value = stim_current;
        } 
        else 
        {
            stim_value = 0.0;
        }

        config->spatial_stim_currents[i] = stim_value;

    }
}

SET_SPATIAL_STIM(stim_if_id_greater_than) {

    uint32_t n_active = the_grid->num_active_cells;
    //struct cell_node **ac = the_grid->active_cells;

    bool stim;
    real stim_current = config->stim_current;
    real stim_value;

    int id = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(int, id, config->config_data.config, "id_limit");

    if(config->spatial_stim_currents) 
    {
        free(config->spatial_stim_currents);
    }

    config->spatial_stim_currents = (real *)malloc(n_active*sizeof(real));

	int i;

    #pragma omp parallel for private(stim, stim_value)
    for (i = 0; i < n_active; i++) 
    {
        stim = i >= id;

        if (stim) 
        {
            stim_value = stim_current;
        } 
        else 
        {
            stim_value = 0.0;
        }

        config->spatial_stim_currents[i] = stim_value;

    }
}

SET_SPATIAL_STIM(stim_concave) 
{

    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    real_cpu max_x_1 = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, max_x_1, config->config_data.config, "max_x_1");

    real_cpu min_x_1 = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, min_x_1, config->config_data.config, "min_x_1");

    real_cpu max_y_1 = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, max_y_1, config->config_data.config, "max_y_1");

    real_cpu min_y_1 = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, min_y_1, config->config_data.config, "min_y_1");

    real_cpu max_x_2 = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, max_x_2, config->config_data.config, "max_x_2");

    real_cpu min_x_2 = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, min_x_2, config->config_data.config, "min_x_2");

    real_cpu max_y_2 = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, max_y_2, config->config_data.config, "max_y_2");

    real_cpu min_y_2 = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, min_y_2, config->config_data.config, "min_y_2");

    real stim_current = config->stim_current;
    real stim_value;

    if(config->spatial_stim_currents) 
    {
        free(config->spatial_stim_currents);
    }

    config->spatial_stim_currents = (real *)malloc(n_active*sizeof(real));

    int i;

#pragma omp parallel for private(stim_value)
    for (i = 0; i < n_active; i++) 
    {

        bool stim_1, stim_2;
        // First corner
        stim_1  = (ac[i]->center_x >= min_x_1) && (ac[i]->center_x <= max_x_1);
        stim_1 &= (ac[i]->center_y >= min_y_1) && (ac[i]->center_y <= max_y_1);
        // Second corner
        stim_2 = (ac[i]->center_x >= min_x_2) && (ac[i]->center_x <= max_x_2);
        stim_2 &= (ac[i]->center_y >= min_y_2) && (ac[i]->center_y <= max_y_2);

        if (stim_1 || stim_2) 
        {
            stim_value = stim_current;
        } 
        else 
        {
            stim_value = 0.0;
        }

        config->spatial_stim_currents[i] = stim_value;

    }
}