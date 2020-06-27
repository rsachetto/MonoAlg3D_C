//
// Created by sachetto on 01/10/17.
//

#include <unistd.h>

#include "../config/extra_data_config.h"
#include "../config_helpers/config_helpers.h"
#include "../libraries_common/common_data_structures.h"
#include "../utils/file_utils.h"


real* set_commom_schemia_data(struct config *config, uint32_t num_cells, int num_par, size_t *extra_data_size) {

    *extra_data_size = sizeof(real)*(num_cells + num_par);

    real *extra_data = (real*)malloc(*extra_data_size);

    real atpi = 6.8;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, atpi, config->config_data, "atpi");

    real Ko = 5.4;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Ko, config->config_data, "Ko");

    real Ki = 138.3;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Ki, config->config_data, "Ki");

    real GNa_multiplicator = 1.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, GNa_multiplicator, config->config_data, "GNa_multiplicator");

    real GCaL_multiplicator = 1.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, GCaL_multiplicator, config->config_data, "GCaL_multiplicator");

    real INaCa_multiplicator = 1.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, INaCa_multiplicator, config->config_data, "INaCa_multiplicator");

    real Vm_modifier = 0.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Vm_modifier, config->config_data, "Vm_modifier");

    extra_data[0] = atpi;
    extra_data[1] = Ko;
    extra_data[2] = Ki;
    extra_data[3] = Vm_modifier;
    extra_data[4] = GNa_multiplicator;
    extra_data[5] = GCaL_multiplicator;
    extra_data[6] = INaCa_multiplicator;

    return extra_data;

}

SET_EXTRA_DATA(set_extra_data_for_fibrosis_sphere) {

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node ** ac = the_grid->active_cells;

    real *fibs = NULL;

    real plain_center = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_center, config->config_data, "plain_center");

    real border_zone_size = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, border_zone_size, config->config_data, "border_zone_size");

    real sphere_radius = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sphere_radius, config->config_data, "sphere_radius");

    int num_par = 7;
    fibs = set_commom_schemia_data(config, num_active_cells, num_par, extra_data_size);

    OMP(parallel for)
    for (uint32_t i = 0; i < num_active_cells; i++) {

        if(FIBROTIC(ac[i])) {
            fibs[i+num_par] = 0.0;
        }
        else if(BORDER_ZONE(ac[i])) {

            real center_x = (real)ac[i]->center.x;
            real center_y = (real)ac[i]->center.y;
            //TODO: Maybe we want the distance from the Z as well
            //real center_z = (real)ac[i]->center_z;

            real distanceFromCenter = sqrtf((center_x - plain_center)*(center_x - plain_center) + (center_y - plain_center)*(center_y - plain_center));
            distanceFromCenter = (distanceFromCenter - sphere_radius)/border_zone_size;
            fibs[i+num_par] = distanceFromCenter;

        }
        else {
            fibs[i+num_par] = 1.0;
        }

    }

    return (void*)fibs;
}

SET_EXTRA_DATA(set_extra_data_for_fibrosis_plain) {

    uint32_t num_active_cells = the_grid->num_active_cells;
    int num_par = 7;

    real *fibs = NULL;

    fibs = set_commom_schemia_data(config, num_active_cells, num_par, extra_data_size);

    for(uint32_t i = num_par; i < num_active_cells + num_par; i++) {
        fibs[i] = 0.0;
    }

    return (void*)fibs;
}

SET_EXTRA_DATA(set_extra_data_for_no_fibrosis) {

    uint32_t num_active_cells = the_grid->num_active_cells;

    int num_par = 7;
    real *fibs = NULL;

    fibs = set_commom_schemia_data(config, num_active_cells, num_par, extra_data_size);

    for(uint32_t i = num_par; i < num_active_cells + num_par; i++) {
        fibs[i] = 1.0;
    }

    return (void*)fibs;
}

SET_EXTRA_DATA(set_extra_data_for_human_full_mesh) {

    uint32_t num_active_cells = the_grid->num_active_cells;

     int num_par = 7;
    real *fibs = NULL;
    fibs = set_commom_schemia_data(config, num_active_cells, num_par, extra_data_size);

    struct cell_node ** ac = the_grid->active_cells;

    real_cpu small_scar_center_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, small_scar_center_x, config->config_data, "small_scar_center_x");

    real_cpu small_scar_center_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, small_scar_center_y, config->config_data, "small_scar_center_y");

    real_cpu small_scar_center_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, small_scar_center_z, config->config_data, "small_scar_center_z");

    real_cpu big_scar_center_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, big_scar_center_x, config->config_data, "big_scar_center_x");

    real_cpu big_scar_center_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, big_scar_center_y, config->config_data, "big_scar_center_y");

    real_cpu big_scar_center_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, big_scar_center_z, config->config_data, "big_scar_center_z");

    real_cpu bz_size_big = 0;
    real_cpu bz_size_small = 0;
    real_cpu dist_big = 0;
    real_cpu dist_small = 0;

	uint32_t i;
	bool fibrotic, border_zone;
	char scar_type;

	OMP(parallel for private(dist_big, dist_small))
    for (i = 0; i < num_active_cells; i++) {

        border_zone = BORDER_ZONE(ac[i]);
        scar_type = SCAR_TYPE(ac[i]);

        if (ac[i]->active && border_zone) {
            real_cpu center_x = ac[i]->center.x;
            real_cpu center_y = ac[i]->center.y;
            real_cpu center_z = ac[i]->center.z;
            if(scar_type == 'b') {
                dist_big = sqrt((center_x - big_scar_center_x) * (center_x - big_scar_center_x) +
                                (center_y - big_scar_center_y) * (center_y - big_scar_center_y) +
                                (center_z - big_scar_center_z) * (center_z - big_scar_center_z));
                OMP(critical(big))
                if (dist_big > bz_size_big) {
                    bz_size_big = dist_big;
                }
            }
            else if(scar_type == 's') {
                dist_small = sqrt((center_x - small_scar_center_x) * (center_x - small_scar_center_x) +
                                  (center_y - small_scar_center_y) * (center_y - small_scar_center_y) +
                                  (center_z - small_scar_center_z) * (center_z - small_scar_center_z));
                OMP(critical(small))
                if (dist_small > bz_size_small) {
                    bz_size_small = dist_small;
                }
            }
        }
    }

    OMP(parallel for private(dist_big, dist_small))
    for (i = 0; i < num_active_cells; i++) {

        if (ac[i]->active) {
            fibrotic = FIBROTIC(ac[i]);
            border_zone = BORDER_ZONE(ac[i]);
            scar_type = SCAR_TYPE(ac[i]);

            if(fibrotic) {
                fibs[i+num_par] = 0.0f;
            }
            else if (border_zone) {
                real_cpu center_x = ac[i]->center.x;
                real_cpu center_y = ac[i]->center.y;
                real_cpu center_z = ac[i]->center.z;
                if(scar_type == 'b') {
                    dist_big = sqrt((center_x - big_scar_center_x) * (center_x - big_scar_center_x) +
                                    (center_y - big_scar_center_y) * (center_y - big_scar_center_y) +
                                    (center_z - big_scar_center_z) * (center_z - big_scar_center_z));
                    fibs[i+num_par] = (real)(dist_big / bz_size_big);

                }
                else if(scar_type == 's') {
                    dist_small = sqrt((center_x - small_scar_center_x) * (center_x - small_scar_center_x) +
                                      (center_y - small_scar_center_y) * (center_y - small_scar_center_y) +
                                      (center_z - small_scar_center_z) * (center_z - small_scar_center_z));
                    fibs[i+num_par] = (real)(dist_small / bz_size_small);
                }
                else {
                    fibs[i+num_par] = 1.0f;
                }
            }
        }
    }

    return (void*)fibs;
}

SET_EXTRA_DATA(set_extra_data_for_scar_wedge) {

    uint32_t num_active_cells = the_grid->num_active_cells;
    real *fibs = NULL;

    int num_par = 7;
    fibs = set_commom_schemia_data(config, num_active_cells, num_par, extra_data_size);

    struct cell_node ** ac = the_grid->active_cells;

    char *scar_size;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR (scar_size, config->config_data, "scar_size");

    uint8_t size_code;

    if(strcmp(scar_size, "big") == 0) {
        size_code = 0;
    }
    else if(strcmp(scar_size, "small") == 0) {
        size_code = 1;
    }
    else {
        printf("Function: set_extra_data_for_scar_edge, invalid scar size %s. Valid sizes are big or small. Exiting!\n", scar_size);
        exit(EXIT_FAILURE);
    }

    real_cpu scar_center_x;
    real_cpu scar_center_y;
    real_cpu scar_center_z;

    ////Fibrosis configuration
    //BIG SCAR
    if(size_code == 0) {
        scar_center_x = 95300;
        scar_center_y = 81600;
        scar_center_z = 36800;
    }
    else {
        scar_center_x = 52469;
        scar_center_y = 83225;
        scar_center_z = 24791;
    }


    real_cpu bz_size = 0.0;
    real_cpu dist;

	uint32_t i;
	bool border_zone, fibrotic;

	OMP(parallel for private(dist))
    for (i = 0; i < num_active_cells; i++) {
        if(ac[i]->active) {
            border_zone = BORDER_ZONE(ac[i]);
            if(border_zone) {
                real_cpu center_x = ac[i]->center.x;
                real_cpu center_y = ac[i]->center.y;
                real_cpu center_z = ac[i]->center.z;
                dist =  sqrt((center_x - scar_center_x)*(center_x - scar_center_x) + (center_y - scar_center_y)*(center_y - scar_center_y)  + (center_z - scar_center_z)*(center_z - scar_center_z)  );
                OMP(critical)
                if(dist > bz_size) {
                    bz_size = dist;
                }
            }

        }
    }

    OMP(parallel for private(dist))
    for (i = 0; i < num_active_cells; i++) {

        if(ac[i]->active) {

            border_zone = BORDER_ZONE(ac[i]);
            fibrotic = FIBROTIC(ac[i]);

            if(fibrotic) {
                fibs[i+num_par] = 0.0;
            }
            else if(border_zone) {
                real_cpu center_x = ac[i]->center.x;
                real_cpu center_y = ac[i]->center.y;
                real_cpu center_z = ac[i]->center.z;
                dist =  sqrt((center_x - scar_center_x)*(center_x - scar_center_x) + (center_y - scar_center_y)*(center_y - scar_center_y)  + (center_z - scar_center_z)*(center_z - scar_center_z)  );
                dist = dist/bz_size;

                fibs[i + num_par] = (real)dist;

            }
            else {
                fibs[i + num_par] = 1.0f;
            }

        }
    }

    return (void*)fibs;
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

SET_EXTRA_DATA(set_extra_data_for_fibrosis_sphere_atpi_changed) {


    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node ** ac = the_grid->active_cells;

    real plain_center = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_center, config->config_data, "plain_center");

    real border_zone_size = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, border_zone_size, config->config_data, "border_zone_size");

    real sphere_radius = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sphere_radius, config->config_data, "sphere_radius");
    
    int num_par = 6;
    int num_tt_par = 12;

//num_tt_par = 12 initial conditions of tt3, num_par 6, extra data

    *extra_data_size = sizeof(real)*(num_par+num_tt_par+num_active_cells);

    real *extra_data = (real*)malloc(*extra_data_size);

    real atpi = 6.8;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, atpi, config->config_data, "atpi");

    real Ko = 5.4;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Ko, config->config_data, "Ko");

    real Ki = 138.3;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Ki, config->config_data, "Ki");

    real GNa_multiplicator = 1.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, GNa_multiplicator, config->config_data, "GNa_multiplicator");

    real GCa_multiplicator = 1.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, GCa_multiplicator, config->config_data, "GCa_multiplicator");

    real Vm_modifier = 0.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Vm_modifier, config->config_data, "Vm_modifier");

    // Extra parameters section
    extra_data[0] = atpi;
    extra_data[1] = Ko;
    extra_data[2] = Ki;
    extra_data[3] = Vm_modifier;
    extra_data[4] = GNa_multiplicator;
    extra_data[5] = GCa_multiplicator;
  
    // Extra initial conditions section
    extra_data[6] = -79.5089;
    extra_data[7] = 0.0056681;
    extra_data[8] = 0.554756;
    extra_data[9] = 0.54673;
    extra_data[10] = 0.000565801;
    extra_data[11] = 0.00486328;
    extra_data[12] = 0.787571;
    extra_data[13] = 0.998604;
    extra_data[14] = 0.998842;
    extra_data[15] = 7.23073e-05;
    extra_data[16] = 6.2706e-08;
    extra_data[17] = 0.412462;

    // Fibrotic cells configuration
	#pragma omp parallel for
    for (uint32_t i = 0; i < num_active_cells; i++) {

        if(FIBROTIC(ac[i])) {
            extra_data[i+num_par+num_tt_par] = 0.0;
        }
        else if(BORDER_ZONE(ac[i])) {

            real center_x = (real)ac[i]->center.x;
            real center_y = (real)ac[i]->center.y;
            //TODO: Maybe we want the distance from the Z as well
            //real center_z = (real)ac[i]->center_z;

            real distanceFromCenter = sqrtf((center_x - plain_center)*(center_x - plain_center) + (center_y - plain_center)*(center_y - plain_center));
            distanceFromCenter = (distanceFromCenter - sphere_radius)/border_zone_size;
            extra_data[i+num_par+num_tt_par] = distanceFromCenter;

        }
        else {
            extra_data[i+num_par+num_tt_par] = 1.0;
        }

    }

    return (void*)extra_data;
}
