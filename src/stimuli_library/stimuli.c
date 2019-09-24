//
// Created by sachetto on 13/10/17.
//

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#include "../alg/grid/grid.h"
#include "../config/stim_config.h"
#include "../config_helpers/config_helpers.h"
#include "../utils/utils.h"

#define SET_STIM_VALUE(i, stim_value) ((real *)(config->persistent_data))[i] = stim_value

#define ALLOCATE_STIMS()                                                                                                \
	if(adaptive) {                                                                                            \
		if(config->persistent_data) {                                                                                   \
			free(config->persistent_data);                                                                              \
		}                                                                                                               \
		config->persistent_data = (real *)malloc(n_active * sizeof(real));                                              \
	}                                                                                                                   \
	else {  																											\
		if(!config->persistent_data) {                                                                                  \
			config->persistent_data = (real *)malloc(n_active * sizeof(real));                                          \
		}                                                                                                               \
	}

SET_SPATIAL_STIM(set_benchmark_spatial_stim) {

    bool stim;
    real stim_current = 0.0;

    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_current, config->config_data, "current");
    real stim_value;

    ALLOCATE_STIMS();

    #pragma omp parallel for private(stim, stim_value)
    for(uint32_t i = 0; i < n_active; i++) {

        stim = ac[i]->center.x > 5500.0;
        stim &= ac[i]->center.x < 7000.0;
        stim &= ac[i]->center.y < 1500.0;
        stim &= ac[i]->center.z < 1500.0;

        if(stim) {
            stim_value = stim_current;
        } else {
            stim_value = 0.0;
        }

        SET_STIM_VALUE(i, stim_value);
        //((real*)(config->persistent_data))[i] = stim_value;
    }
}

SET_SPATIAL_STIM(stim_if_x_less_than) {

    real stim_current = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_current, config->config_data, "current");

    bool stim;
    real stim_value;

    real_cpu x_limit = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, x_limit, config->config_data, "x_limit");

    uint32_t i;

    ALLOCATE_STIMS();

    #pragma omp parallel for private(stim, stim_value)
    for(i = 0; i < n_active; i++) {
        stim = ac[i]->center.x < x_limit;

        if(stim) {
            stim_value = stim_current;
        } else {
            stim_value = 0.0;
        }

        SET_STIM_VALUE(i, stim_value);
    }
}

SET_SPATIAL_STIM(stim_if_y_less_than) {

    real stim_current = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_current, config->config_data, "current");

    bool stim;
    real stim_value;

    real_cpu y_limit = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, y_limit, config->config_data, "y_limit");

    uint32_t i;

    ALLOCATE_STIMS();

    #pragma omp parallel for private(stim, stim_value)
    for (i = 0; i < n_active; i++) {
        stim = ac[i]->center.y < y_limit;

        if (stim) {
            stim_value = stim_current;
        } else {
            stim_value = 0.0;
        }

        SET_STIM_VALUE(i, stim_value);

    }
}

SET_SPATIAL_STIM(stim_if_z_less_than) {

    real stim_current = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_current, config->config_data, "current");

    bool stim;
    real stim_value;

    real_cpu z_limit = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, z_limit, config->config_data, "z_limit");

    uint32_t i;

    ALLOCATE_STIMS();

    #pragma omp parallel for private(stim, stim_value)
    for (i = 0; i < n_active; i++) {
        stim = ac[i]->center.z < z_limit;

        if (stim) {
            stim_value = stim_current;
        } else {
            stim_value = 0.0;
        }

        SET_STIM_VALUE(i, stim_value);

    }
}

SET_SPATIAL_STIM(stim_if_y_greater_or_equal_than) {

    real stim_current = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_current, config->config_data, "current");

    bool stim;
    real stim_value;

    real_cpu y_limit = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, y_limit, config->config_data, "y_limit");

    uint32_t i;

    ALLOCATE_STIMS();

    #pragma omp parallel for private(stim, stim_value)
    for (i = 0; i < n_active; i++) {
        stim = ac[i]->center.y >= y_limit;

        if (stim) {
            stim_value = stim_current;
        } else {
            stim_value = 0.0;
        }

        SET_STIM_VALUE(i, stim_value);

    }
}

SET_SPATIAL_STIM(set_stim_from_file) {

    char *stim_file = NULL;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(stim_file, config->config_data, "stim_file");

    size_t s_size;

    bool stim;
    real stim_current = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_current, config->config_data, "current");

    real stim_value;

    FILE *s_file = fopen(stim_file, "r");

    if(!s_file) {
        fprintf(stderr, "Error opening stim file %s! Exiting!\n", stim_file);
        exit(EXIT_FAILURE);
    }

    fscanf(s_file, "%zu\n", &s_size);

    real_cpu **cell_stims = (real_cpu **)malloc(sizeof(real_cpu *) * s_size);
    for(int i = 0; i < s_size; i++) {
        cell_stims[i] = (real_cpu *)malloc(sizeof(real_cpu) * 3);
        if(cell_stims[i] == NULL) {
            fprintf(stderr, "Failed to allocate memory for the stim file\n");
            exit(0);
        }

        fscanf(s_file, "%lf %lf %lf\n", &cell_stims[i][0], &cell_stims[i][1], &cell_stims[i][2]);
    }

    sort_vector(cell_stims, s_size);

    fclose(s_file);

    uint32_t i;

    ALLOCATE_STIMS();

    #pragma omp parallel for private(stim, stim_value)
    for(i = 0; i < n_active; i++) {

        real_cpu center_x = ac[i]->center.x;
        real_cpu center_y = ac[i]->center.y;
        real_cpu center_z = ac[i]->center.z;

        int index = inside_mesh(cell_stims, center_x, center_y, center_z, 0, s_size - 1);
        stim = (index != -1);

        if(stim) {
            stim_value = stim_current;
        } else {
            stim_value = 0.0;
        }

        SET_STIM_VALUE(i, stim_value);
    }
}

SET_SPATIAL_STIM(stim_if_x_greater_equal_than) {

    real_cpu x_limit = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, x_limit, config->config_data, "x_limit");

    real stim_current = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_current, config->config_data, "current");

    real stim_value;

    uint32_t i;

    ALLOCATE_STIMS();

    #pragma omp parallel for private(stim_value)
    for(i = 0; i < n_active; i++) {
        bool stim = (ac[i]->center.x >= x_limit);

        if(stim) {
            stim_value = stim_current;
        } else {
            stim_value = 0.0;
        }

        SET_STIM_VALUE(i, stim_value);
    }
}

SET_SPATIAL_STIM(stim_base_mouse) {

    real_cpu stim_size = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, stim_size, config->config_data, "stim_size");

    real stim_current = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_current, config->config_data, "current");

    real stim_value;

    uint32_t i;

    ALLOCATE_STIMS();

    #pragma omp parallel for private(stim_value)
    for(i = 0; i < n_active; i++) {

        bool stim;
        stim = (ac[i]->center.x >= 3000.0 - stim_size) && (ac[i]->center.x <= 3000.0 + stim_size);
        stim &= (ac[i]->center.y >= 2400.0 - stim_size) && (ac[i]->center.y <= 2400.0 + stim_size);
        stim &= (ac[i]->center.z >= 300 - stim_size) && (ac[i]->center.z <= 300 + stim_size);

        if(stim) {
            stim_value = stim_current;
        } else {
            stim_value = 0.0;
        }

        SET_STIM_VALUE(i, stim_value);
    }
}

SET_SPATIAL_STIM(stim_mouse_spiral) {

    real stim_current = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_current, config->config_data, "current");

    real stim_value;

    int i;

    ALLOCATE_STIMS();

    #pragma omp parallel for private(stim_value)
    for(i = 0; i < n_active; i++) {

        bool stim;

        stim = (ac[i]->center.x >= 3000.0) && (ac[i]->center.x <= 6000.0);
        stim &= (ac[i]->center.y >= 1940.0) && (ac[i]->center.y <= 6100.0);
        stim &= (ac[i]->center.z >= 2230.0) && (ac[i]->center.z <= 5800.0);

        if(stim) {
            stim_value = stim_current;
        } else {
            stim_value = 0.0;
        }

        SET_STIM_VALUE(i, stim_value);
    }
}

SET_SPATIAL_STIM(stim_x_y_limits) {

    real_cpu max_x = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, max_x, config->config_data, "max_x");

    real_cpu min_x = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, min_x, config->config_data, "min_x");

    real_cpu max_y = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, max_y, config->config_data, "max_y");

    real_cpu min_y = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, min_y, config->config_data, "min_y");

    real stim_current = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_current, config->config_data, "current");

    real stim_value;

    ALLOCATE_STIMS();

    int i;

#pragma omp parallel for private(stim_value)
    for(i = 0; i < n_active; i++) {

        bool stim;
        stim = (ac[i]->center.x >= min_x) && (ac[i]->center.x <= max_x);
        stim &= (ac[i]->center.y >= min_y) && (ac[i]->center.y <= max_y);

        if(stim) {
            stim_value = stim_current;
        } else {
            stim_value = 0.0;
        }

        SET_STIM_VALUE(i, stim_value);
    }
}

SET_SPATIAL_STIM(stim_x_y_z_limits) {

    real_cpu max_x = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, max_x, config->config_data, "max_x");

    real_cpu min_x = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, min_x, config->config_data, "min_x");

    real_cpu max_y = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, max_y, config->config_data, "max_y");

    real_cpu min_y = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, min_y, config->config_data, "min_y");

    real_cpu max_z = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, max_z, config->config_data, "max_z");

    real_cpu min_z = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, min_z, config->config_data, "min_z");

    real stim_current = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_current, config->config_data, "current");

    real stim_value;

    uint32_t i;

    ALLOCATE_STIMS();

    #pragma omp parallel for private(stim_value)
    for(i = 0; i < n_active; i++) {

        bool stim;
        stim = (ac[i]->center.x >= min_x) && (ac[i]->center.x <= max_x);
        stim &= (ac[i]->center.y >= min_y) && (ac[i]->center.y <= max_y);
        stim &= (ac[i]->center.z >= min_z) && (ac[i]->center.z <= max_z);

        if(stim) {
            stim_value = stim_current;
        } else {
            stim_value = 0.0;
        }

        SET_STIM_VALUE(i, stim_value);
    }
}

// ***********************************************************************************************************
// New Berg's stimulus
SET_SPATIAL_STIM(stim_if_inside_circle_than) {

    bool stim;
    real stim_value;

    real stim_current = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_current, config->config_data, "current");

    real_cpu center_x = 0.0;
    real_cpu center_y = 0.0;
    real_cpu center_z = 0.0;
    real_cpu radius = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, center_x, config->config_data, "center_x");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, center_y, config->config_data, "center_y");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, center_z, config->config_data, "center_z");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, radius, config->config_data, "radius");

    uint32_t i;

    ALLOCATE_STIMS();

    #pragma omp parallel for private(stim, stim_value)
    for(i = 0; i < n_active; i++) {
        real_cpu dist = sqrt(pow(ac[i]->center.x - center_x, 2) + pow(ac[i]->center.y - center_y, 2) +
                             pow(ac[i]->center.z - center_z, 2));
        stim = dist <= radius;

        if(stim) {
            stim_value = stim_current;
        } else {
            stim_value = 0.0;
        }

        SET_STIM_VALUE(i, stim_value);
    }
}

SET_SPATIAL_STIM(stim_if_id_less_than) {

    bool stim;
    real stim_value;

    real stim_current = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_current, config->config_data, "current");

    int id = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(int, id, config->config_data, "id_limit");

    int i;

    ALLOCATE_STIMS();

    #pragma omp parallel for private(stim, stim_value)
    for(i = 0; i < n_active; i++) {
        stim = i <= id;

        if(stim) {
            stim_value = stim_current;
        } else {
            stim_value = 0.0;
        }

        SET_STIM_VALUE(i, stim_value);
    }
}

SET_SPATIAL_STIM(stim_if_id_greater_than) {

    bool stim;
    real stim_value;

    real stim_current = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_current, config->config_data, "current");

    int id = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(int, id, config->config_data, "id_limit");

    int i;

    ALLOCATE_STIMS();

    #pragma omp parallel for private(stim, stim_value)
    for(i = 0; i < n_active; i++) {
        stim = i >= id;

        if(stim) {
            stim_value = stim_current;
        } else {
            stim_value = 0.0;
        }

        SET_STIM_VALUE(i, stim_value);
    }
}

SET_SPATIAL_STIM(stim_concave) {

    real_cpu max_x_1 = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, max_x_1, config->config_data, "max_x_1");

    real_cpu min_x_1 = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, min_x_1, config->config_data, "min_x_1");

    real_cpu max_y_1 = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, max_y_1, config->config_data, "max_y_1");

    real_cpu min_y_1 = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, min_y_1, config->config_data, "min_y_1");

    real_cpu max_x_2 = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, max_x_2, config->config_data, "max_x_2");

    real_cpu min_x_2 = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, min_x_2, config->config_data, "min_x_2");

    real_cpu max_y_2 = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, max_y_2, config->config_data, "max_y_2");

    real_cpu min_y_2 = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, min_y_2, config->config_data, "min_y_2");

    real stim_current = 0.0;
    real stim_value;

    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_current, config->config_data, "current");

    uint32_t i;

    ALLOCATE_STIMS();

    #pragma omp parallel for private(stim_value)
    for(i = 0; i < n_active; i++) {

        bool stim_1, stim_2;
        // First corner
        stim_1 = (ac[i]->center.x >= min_x_1) && (ac[i]->center.x <= max_x_1);
        stim_1 &= (ac[i]->center.y >= min_y_1) && (ac[i]->center.y <= max_y_1);
        // Second corner
        stim_2 = (ac[i]->center.x >= min_x_2) && (ac[i]->center.x <= max_x_2);
        stim_2 &= (ac[i]->center.y >= min_y_2) && (ac[i]->center.y <= max_y_2);

        if(stim_1 || stim_2) {
            stim_value = stim_current;
        } else {
            stim_value = 0.0;
        }

        SET_STIM_VALUE(i, stim_value);
    }
}

/*
SET_SPATIAL_STIM(stim_purkinje_if_id_less_than) 
{

    bool stim;
    real stim_value;

    real stim_current = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_current, config->config_data, "current");

    int id = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(int, id, config->config_data, "id_limit");

    int i;

    ALLOCATE_STIMS();

    #pragma omp parallel for private(stim, stim_value)
    for(i = 0; i < n_active; i++) 
    {
        stim = i <= id;

        if(stim) 
        {
            stim_value = stim_current;
        } 
        else 
        {
            stim_value = 0.0;
        }

        SET_STIM_VALUE(i, stim_value);
    }
}

SET_SPATIAL_STIM(stim_purkinje_if_x_less_than) 
{

    // Total number of cells on the grid Purkinje + Tissue
    uint32_t n_active_tissue = the_grid->num_active_cells;
    uint32_t n_active_purkinje = the_grid->the_purkinje->num_active_purkinje_cells;
    uint32_t n_active = n_active_purkinje + n_active_tissue;

    struct cell_node **ac_purkinje = the_grid->the_purkinje->purkinje_cells;
    struct cell_node **ac_tissue = the_grid->active_cells;

    real stim_current = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_current, config->config_data, "current");

    bool stim;
    real stim_value;

    real_cpu x_limit = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, x_limit, config->config_data, "x_limit");

    uint32_t i, j;

    ALLOCATE_STIMS();

    #pragma omp parallel for private(stim, stim_value)
    for(i = 0; i < n_active_purkinje; i++) 
    {
        stim = ac_purkinje[i]->center.x < x_limit;

        if(stim) 
        {
            stim_value = stim_current;
        } 
        else 
        {
            stim_value = 0.0;
        }

        SET_STIM_VALUE(i, stim_value);
    }

    i = n_active_purkinje;

    for(j = 0; j < n_active_tissue; j++) 
    {
        stim = ac_tissue[j]->center.x < x_limit;

        if(stim) 
        {
            stim_value = stim_current;
        } 
        else 
        {
            stim_value = 0.0;
        }

        SET_STIM_VALUE(i, stim_value);
        
        i++;
    }

}

SET_SPATIAL_STIM(stim_purkinje_if_x_greater_equal_than) {

    // Total number of cells on the grid Purkinje + Tissue
    uint32_t n_active_tissue = the_grid->num_active_cells;
    uint32_t n_active_purkinje = the_grid->the_purkinje->num_active_purkinje_cells;
    uint32_t n_active = n_active_purkinje + n_active_tissue;

    struct cell_node **ac_purkinje = the_grid->the_purkinje->purkinje_cells;
    struct cell_node **ac_tissue = the_grid->active_cells;

    real stim_current = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_current, config->config_data, "current");

    bool stim;
    real stim_value;

    real_cpu x_limit = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, x_limit, config->config_data, "x_limit");

    uint32_t i, j;

    ALLOCATE_STIMS();

    #pragma omp parallel for private(stim, stim_value)
    for(i = 0; i < n_active_tissue; i++) 
    {
        stim = ac_tissue[i]->center.x >= x_limit;

        if(stim) 
        {
            stim_value = stim_current;
        } 
        else 
        {
            stim_value = 0.0;
        }

        SET_STIM_VALUE(i, stim_value);
    }

    i = n_active_purkinje;

    for(j = 0; j < n_active_purkinje; j++) 
    {
        stim = ac_purkinje[j]->center.x >= x_limit;

        if(stim) 
        {
            stim_value = stim_current;
        } 
        else 
        {
            stim_value = 0.0;
        }

        SET_STIM_VALUE(i, stim_value);
        
        i++;
    }
}
*/