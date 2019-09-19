//
// Created by berg on 03/04/19.
//

#include "../config/extra_data_config.h"
#include "../config_helpers/config_helpers.h"
#include "../libraries_common/common_data_structures.h"

SET_EXTRA_DATA (set_mixed_model_if_x_less_than)
{
    uint32_t num_active_cells = the_grid->num_active_cells;

    *extra_data_size = sizeof(uint32_t)*(num_active_cells + 1);

    uint32_t *mapping = (uint32_t*)malloc(*extra_data_size);

    struct cell_node ** ac = the_grid->active_cells;

    real x_limit = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, x_limit, config->config_data, "x_limit");

    int i;
    bool inside;

    #pragma omp parallel for
    for (i = 0; i < num_active_cells; i++)
    {
        real center_x = ac[i]->center.x;
        real center_y = ac[i]->center.y;
        real center_z = ac[i]->center.z;

        inside = (center_x <= x_limit);

        if (inside)
            mapping[i] = 0;
        else
            mapping[i] = 1;        
    }

    return (void*)mapping;
}

SET_EXTRA_DATA (set_mixed_model_purkinje_and_tissue)
{
    uint32_t num_active_tissue_cells = the_grid->num_active_cells;
    uint32_t num_active_purkinje_cells = the_grid->the_purkinje->num_active_purkinje_cells; 
    uint32_t num_active_cells = num_active_tissue_cells + num_active_purkinje_cells;

    *extra_data_size = sizeof(uint32_t)*(num_active_cells + 1);

    uint32_t *mapping = (uint32_t*)malloc(*extra_data_size);

    struct cell_node ** ac = the_grid->active_cells;
    struct cell_node ** ac_purkinje = the_grid->the_purkinje->purkinje_cells;

    int i;
    bool inside;

    // Purkinje section
    #pragma omp parallel for
    for (i = 0; i < num_active_purkinje_cells; i++)
    {
        mapping[i] = 0;   
    }

    // Tissue section
    #pragma omp parallel for
    for (i = num_active_purkinje_cells; i < num_active_cells; i++)
    {
        mapping[i] = 1;        
    }

    return (void*)mapping;
}