//
// Created by sachetto on 01/10/17.
//

#include <unistd.h>

#include "../config/extra_data_config.h"
#include "../config_helpers/config_helpers.h"
#include "../libraries_common/common_data_structures.h"
#include "../domains_library/mesh_info_data.h"
#include "helper_functions.h"


SET_EXTRA_DATA(set_extra_data_for_fibrosis_sphere) {

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node ** ac = the_grid->active_cells;

    struct extra_data_for_fibrosis *extra_data = NULL;

    real plain_center = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_center, config, "plain_center");

    real border_zone_size = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, border_zone_size, config, "border_zone_size");

    real sphere_radius = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sphere_radius, config, "sphere_radius");

    extra_data = set_common_schemia_data(config, num_active_cells);


    OMP(parallel for)
    for (uint32_t i = 0; i < num_active_cells; i++) {

        if(FIBROTIC(ac[i])) {
            extra_data->fibrosis[i] = 0.0;
        }
        else if(BORDER_ZONE(ac[i])) {

            real center_x = (real)ac[i]->center.x;
            real center_y = (real)ac[i]->center.y;
            //TODO: Maybe we want the distance from the Z as well
            //real center_z = (real)ac[i]->center_z;

            real distanceFromCenter = sqrtf((center_x - plain_center)*(center_x - plain_center) + (center_y - plain_center)*(center_y - plain_center));
            distanceFromCenter = (distanceFromCenter - sphere_radius)/border_zone_size;
            extra_data->fibrosis[i] = distanceFromCenter;

        }
        else {
            extra_data->fibrosis[i] = 1.0;
        }

    }

    SET_EXTRA_DATA_SIZE(sizeof(struct extra_data_for_fibrosis));

    return (void*)extra_data;
}

SET_EXTRA_DATA(set_extra_data_for_fibrosis_plain) {

    uint32_t num_active_cells = the_grid->num_active_cells;

    struct extra_data_for_fibrosis *extra_data = NULL;

    extra_data = set_common_schemia_data(config, num_active_cells);

    OMP(parallel for)
    for(uint32_t i = 0; i < num_active_cells; i++) {
        extra_data->fibrosis[i] = 0.0;
    }

   SET_EXTRA_DATA_SIZE(sizeof(struct extra_data_for_fibrosis));

    return (void*)extra_data;
}

SET_EXTRA_DATA(set_extra_data_for_no_fibrosis) {

    uint32_t num_active_cells = the_grid->num_active_cells;

    struct extra_data_for_fibrosis *extra_data = NULL;

    extra_data = set_common_schemia_data(config, num_active_cells);

    OMP(parallel for)
    for(uint32_t i = 0; i < num_active_cells; i++) {
        extra_data->fibrosis[i] = 1.0;
    }

    return (void*)extra_data;
}

SET_EXTRA_DATA(set_extra_data_for_benchmark) {

    *extra_data_size = sizeof(real)*19;

    real *initial_conditions = (real*)malloc(*extra_data_size);

    // Initial conditions  // Var      Units          Initial value
    initial_conditions[ 0] = -85.423f;  // V;       millivolt;     -85.423
    initial_conditions[ 1] = 0.0165;   // Xr1;     dimensionless; 0.0165
    initial_conditions[ 2] = 0.473;    // Xr2;     dimensionless; 0.473
    initial_conditions[ 3] = 0.0174;   // Xs;      dimensionless; 0.0174
    initial_conditions[ 4] = 0.00165;  // m;       dimensionless; 0.00165
    initial_conditions[ 5] = 0.749;    // h;       dimensionless; 0.749
    initial_conditions[ 6] = 0.6788;   // j;       dimensionless; 0.6788
    initial_conditions[ 7] = 3.288e-5; // d;       dimensionless; 3.288e-5
    initial_conditions[ 8] = 0.7026;   // f;       dimensionless; 0.7026
    initial_conditions[ 9] = 0.9526;   // f2;      dimensionless; 0.9526
    initial_conditions[10] = 0.9942;   // fCass;   dimensionless; 0.9942
    initial_conditions[11] = 0.999998; // s;       dimensionless; 0.999998
    initial_conditions[12] = 2.347e-8; // r;       dimensionless; 2.347e-8
    initial_conditions[13] = 0.000153; // Ca_i;    millimolar;    0.000153
    initial_conditions[14] = 4.272;    // Ca_SR;   millimolar;    4.272
    initial_conditions[15] = 0.00042;  // Ca_ss;   millimolar;    0.00042
    initial_conditions[16] = 0.8978;   // R_prime; dimensionless; 0.8978
    initial_conditions[17] = 10.132;   // Na_i;    millimolar;    10.132
    initial_conditions[18] = 138.52;   // K_i;     millimolar;    138.52

    return (void*)initial_conditions;
}

SET_EXTRA_DATA (set_mixed_model_if_x_less_than)
{
    uint32_t num_active_cells = the_grid->num_active_cells;

    *extra_data_size = sizeof(uint32_t)*(num_active_cells);

    uint32_t *mapping = (uint32_t*)malloc(*extra_data_size);

    struct cell_node ** ac = the_grid->active_cells;

    real x_limit = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, x_limit, config, "x_limit");

    int i;
    bool inside;

    OMP(parallel for)
    for (i = 0; i < num_active_cells; i++)
    {
        real center_x = ac[i]->center.x;

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
    uint32_t num_active_purkinje_cells = the_grid->purkinje->num_active_purkinje_cells; 
    uint32_t num_active_cells = num_active_tissue_cells + num_active_purkinje_cells;

    *extra_data_size = sizeof(uint32_t)*(num_active_cells + 1);

    uint32_t *mapping = (uint32_t*)malloc(*extra_data_size);

    int i;

    // Purkinje section
    OMP(parallel for)
    for (i = 0; i < num_active_purkinje_cells; i++) {
        mapping[i] = 0;   
    }

    // Tissue section
    OMP(parallel for)
    for (i = num_active_purkinje_cells; i < num_active_cells; i++) {
        mapping[i] = 1;        
    }

    return (void*)mapping;
}