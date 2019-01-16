//
// Created by sachetto on 01/10/17.
//

#include "../config/extra_data_config.h"
#include "../libraries_common/config_helpers.h"
#include "../libraries_common/common_data_structures.h"


SET_EXTRA_DATA(set_extra_data_for_fibrosis_sphere) {

    uint32_t num_active_cells = the_grid->num_active_cells;

    *extra_data_size = sizeof(float)*(num_active_cells+1);

    float *fibs = (float*)malloc(*extra_data_size);

    struct cell_node ** ac = the_grid->active_cells;

    real atpi = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, atpi, config, "atpi");

    float plain_center = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_center, config, "plain_center");

    float border_zone_size = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, border_zone_size, config, "border_zone_size");

    float sphere_radius = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sphere_radius, config, "sphere_radius");

    fibs[0] = atpi;    

	int i;

	#pragma omp parallel for
    for (i = 0; i < num_active_cells; i++) {

        if(FIBROTIC(ac[i])) {
            fibs[i+1] = 0.0;
        }
        else if(BORDER_ZONE(ac[i])) {

            float center_x = (float)ac[i]->center_x;
            float center_y = (float)ac[i]->center_y;
            //TODO: Maybe we want the distance from the Z as well
            //float center_z = (float)ac[i]->center_z;

            float distanceFromCenter = sqrtf((center_x - plain_center)*(center_x - plain_center) + (center_y - plain_center)*(center_y - plain_center));
            distanceFromCenter = (distanceFromCenter - sphere_radius)/border_zone_size;
            fibs[i+1] = distanceFromCenter;

        }
        else {
            fibs[i+1] = 1.0;
        }

    }


    return (void*)fibs;
}

SET_EXTRA_DATA(set_extra_data_for_fibrosis_plain) {

    uint32_t num_active_cells = the_grid->num_active_cells;

    *extra_data_size = sizeof(float)*(num_active_cells+5);

    float *fibs = (float*)calloc(*extra_data_size, sizeof(float));

    real atpi = 6.8;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, atpi, config, "atpi");

    real Ko = 5.4;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Ko, config, "Ko");

    real Ki_multiplicator = 1.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Ki_multiplicator, config, "Ki_multiplicator");

    real K1_multiplicator = 1.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, K1_multiplicator, config, "K1_multiplicator");

    real acidosis = false;
    GET_PARAMETER_BINARY_VALUE_OR_USE_DEFAULT(acidosis, config, "acidosis");

    fibs[0] = atpi;
    fibs[1] = Ko;
    fibs[2] = Ki_multiplicator;
    fibs[3] = K1_multiplicator;
    fibs[4] = (real)acidosis;


    return (void*)fibs;
}



SET_EXTRA_DATA(set_extra_data_for_no_fibrosis) {

    uint32_t num_active_cells = the_grid->num_active_cells;

    *extra_data_size = sizeof(float)*(num_active_cells+5);

    float *fibs = (float*)calloc(*extra_data_size, sizeof(float));

    real atpi = 6.8;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, atpi, config, "atpi");

    real Ko = 5.4;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Ko, config, "Ko");

    real Ki_multiplicator = 1.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Ki_multiplicator, config, "Ki_multiplicator");

    real K1_multiplicator = 1.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, K1_multiplicator, config, "K1_multiplicator");

    real acidosis = false;
    GET_PARAMETER_BINARY_VALUE_OR_USE_DEFAULT(acidosis, config, "acidosis");


    fibs[0] = atpi;
    fibs[1] = Ko;
    fibs[2] = Ki_multiplicator;
    fibs[3] = K1_multiplicator;
    fibs[4] = (real)acidosis;

    for(int i = 5; i < num_active_cells+5; i++) {
        fibs[i] = 0.0;
    }

    return (void*)fibs;
}

SET_EXTRA_DATA(set_extra_data_for_fibrosis) {

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node ** ac = the_grid->active_cells;

    *extra_data_size = sizeof(real)*(num_active_cells+1);

    real *fibs = (real*)malloc(*extra_data_size);

    fibs[0] = 6.8;

    for(int i = 0; i < num_active_cells; i++) {
        if(FIBROTIC(ac[i])) {
            fibs[i+1] = 0.0;
        }
        else {
            fibs[i+1] = 1.0;
        }
    }

    return (void*)fibs;
}

SET_EXTRA_DATA(set_extra_data_for_human_full_mesh) {

    uint32_t num_active_cells = the_grid->num_active_cells;

    *extra_data_size = sizeof(float)*(num_active_cells+4);

    float *fibs = (float*)calloc(*extra_data_size, sizeof(float));

    real atpi = 6.8;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, atpi, config, "atpi");

    real Ko = 5.4;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Ko, config, "Ko");

    real Ki_multiplicator = 1.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Ki_multiplicator, config, "Ki_multiplicator");

    real acidosis;
    GET_PARAMETER_BINARY_VALUE_OR_USE_DEFAULT(acidosis, config, "acidosis");

    fibs[0] = atpi;
    fibs[1] = Ko;
    fibs[2] = Ki_multiplicator;
    fibs[3] = (real)acidosis;

    struct cell_node ** ac = the_grid->active_cells;

    double small_scar_center_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, small_scar_center_x, config, "small_scar_center_x");

    double small_scar_center_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, small_scar_center_y, config, "small_scar_center_y");

    double small_scar_center_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, small_scar_center_z, config, "small_scar_center_z");

    double big_scar_center_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, big_scar_center_x, config, "big_scar_center_x");

    double big_scar_center_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, big_scar_center_y, config, "big_scar_center_y");

    double big_scar_center_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, big_scar_center_z, config, "big_scar_center_z");

    double bz_size_big = 0;
    double bz_size_small = 0;
    double dist_big = 0;
    double dist_small = 0;

	int i;
	bool fibrotic, border_zone;
	char scar_type;
    
    //#pragma omp parallel for private(dist_big, dist_small) reduction(max: bz_size_big, bz_size_small)
	#pragma omp parallel for private(dist_big, dist_small)
    for (i = 0; i < num_active_cells; i++) {

        border_zone = BORDER_ZONE(ac[i]);
        scar_type = SCAR_TYPE(ac[i]);

        if (ac[i]->active && border_zone) {
            double center_x = ac[i]->center_x;
            double center_y = ac[i]->center_y;
            double center_z = ac[i]->center_z;
            if(scar_type == 'b') {
                dist_big = sqrt((center_x - big_scar_center_x) * (center_x - big_scar_center_x) +
                                (center_y - big_scar_center_y) * (center_y - big_scar_center_y) +
                                (center_z - big_scar_center_z) * (center_z - big_scar_center_z));
				#pragma omp critical(big)
                if (dist_big > bz_size_big) {
                    bz_size_big = dist_big;
                }
            }
            else if(scar_type == 's') {
                dist_small = sqrt((center_x - small_scar_center_x) * (center_x - small_scar_center_x) +
                                  (center_y - small_scar_center_y) * (center_y - small_scar_center_y) +
                                  (center_z - small_scar_center_z) * (center_z - small_scar_center_z));
				#pragma omp critical(small)
                if (dist_small > bz_size_small) {
                    bz_size_small = dist_small;
                }
            }
        }
    }

    #pragma omp parallel for private(dist_big, dist_small)
    for (i = 0; i < num_active_cells; i++) {

        if (ac[i]->active) {
            fibrotic = FIBROTIC(ac[i]);
            border_zone = BORDER_ZONE(ac[i]);
            scar_type = SCAR_TYPE(ac[i]);

            if(fibrotic) {
                fibs[i+1] = 0.0f;
            }
            else if (border_zone) {
                double center_x = ac[i]->center_x;
                double center_y = ac[i]->center_y;
                double center_z = ac[i]->center_z;
                if(scar_type == 'b') {
                    dist_big = sqrt((center_x - big_scar_center_x) * (center_x - big_scar_center_x) +
                                    (center_y - big_scar_center_y) * (center_y - big_scar_center_y) +
                                    (center_z - big_scar_center_z) * (center_z - big_scar_center_z));
                    fibs[i+1] = (real)(dist_big / bz_size_big);
                    
                }
                else if(scar_type == 's') {
                    dist_small = sqrt((center_x - small_scar_center_x) * (center_x - small_scar_center_x) +
                                      (center_y - small_scar_center_y) * (center_y - small_scar_center_y) +
                                      (center_z - small_scar_center_z) * (center_z - small_scar_center_z));
                    fibs[i+1] = (real)(dist_small / bz_size_small);
                }
                else {
                    fibs[i+1] = 1.0f;
                }
            }
        }
    }

    return (void*)fibs;
}

SET_EXTRA_DATA(set_extra_data_for_scar_wedge) {

    uint32_t num_active_cells = the_grid->num_active_cells;

    *extra_data_size = sizeof(float)*(num_active_cells+1);

    float *fibs = (float*)malloc(*extra_data_size);

    struct cell_node ** ac = the_grid->active_cells;

    real atpi = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real,atpi,  config, "atpi");
    fibs[0] = atpi;

    char *scar_size;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR (scar_size, config, "scar_size");

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

    double scar_center_x;
    double scar_center_y;
    double scar_center_z;

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


    double bz_size = 0.0;
    double dist;

	int i;
	bool border_zone, fibrotic;

//    #pragma omp parallel for private(dist) reduction(max: bz_size)
	#pragma omp parallel for private(dist)
    for (i = 0; i < num_active_cells; i++) {
        if(ac[i]->active) {
            border_zone = BORDER_ZONE(ac[i]);
            if(border_zone) {
                double center_x = ac[i]->center_x;
                double center_y = ac[i]->center_y;
                double center_z = ac[i]->center_z;
                dist =  sqrt((center_x - scar_center_x)*(center_x - scar_center_x) + (center_y - scar_center_y)*(center_y - scar_center_y)  + (center_z - scar_center_z)*(center_z - scar_center_z)  );
				#pragma omp critical
                if(dist > bz_size) {
                    bz_size = dist;
                }
            }

        }
    }	

    #pragma omp parallel for private(dist)
    for (i = 0; i < num_active_cells; i++) {

        if(ac[i]->active) {

            border_zone = BORDER_ZONE(ac[i]);
            fibrotic = FIBROTIC(ac[i]);

            if(fibrotic) {
                fibs[i+1] = 0.0;
            }
            else if(border_zone) {
                double center_x = ac[i]->center_x;
                double center_y = ac[i]->center_y;
                double center_z = ac[i]->center_z;
                dist =  sqrt((center_x - scar_center_x)*(center_x - scar_center_x) + (center_y - scar_center_y)*(center_y - scar_center_y)  + (center_z - scar_center_z)*(center_z - scar_center_z)  );
                dist = dist/bz_size;

                fibs[i+1] = (real)dist;

            }
            else {
                fibs[i+1] = 1.0f;
            }

        }
    }

    return (void*)fibs;
}

SET_EXTRA_DATA(set_extra_data_for_benchmark) {

    *extra_data_size = sizeof(float)*19;

    float *initial_conditions = (float*)malloc(*extra_data_size);

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
