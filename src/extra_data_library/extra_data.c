//
// Created by sachetto on 01/10/17.
//

#include "../monodomain/config/extra_data_config.h"
#include "../libraries_common/config_helpers.h"


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

        if(ac[i]->fibrotic) {
            fibs[i+1] = 0.0;
        }
        else if(ac[i]->border_zone) {

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

    *extra_data_size = sizeof(float)*(num_active_cells+1);

    float *fibs = (float*)calloc(num_active_cells+1, sizeof(float));

    real atpi = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, atpi, config, "atpi");

    fibs[0] = atpi;

    return (void*)fibs;
}

SET_EXTRA_DATA(set_extra_data_for_no_fibrosis) {

    uint32_t num_active_cells = the_grid->num_active_cells;

    *extra_data_size = sizeof(float)*(num_active_cells+1);

    float *fibs = (float*)malloc(*extra_data_size);

    real atpi = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real,atpi,  config, "atpi");
    fibs[0] = atpi;

    for(int i = 1; i < num_active_cells+1; i++) {
        fibs[i] = 1.0;
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
        if(ac[i]->fibrotic) {
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

    *extra_data_size = sizeof(float)*(num_active_cells+1);

    float *fibs = (float*)malloc(*extra_data_size);

    struct cell_node ** ac = the_grid->active_cells;

    real atpi = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real,atpi,  config, "atpi");
    fibs[0] = atpi;

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
    
    //#pragma omp parallel for private(dist_big, dist_small) reduction(max: bz_size_big, bz_size_small)
	#pragma omp parallel for private(dist_big, dist_small)
    for (i = 0; i < num_active_cells; i++) {

        if (ac[i]->active && ac[i]->border_zone) {
            double center_x = ac[i]->center_x;
            double center_y = ac[i]->center_y;
            double center_z = ac[i]->center_z;
            if(ac[i]->scar_type == 'b') {
                dist_big = sqrt((center_x - big_scar_center_x) * (center_x - big_scar_center_x) +
                                (center_y - big_scar_center_y) * (center_y - big_scar_center_y) +
                                (center_z - big_scar_center_z) * (center_z - big_scar_center_z));
				#pragma omp critical(big)
                if (dist_big > bz_size_big) {
                    bz_size_big = dist_big;
                }
            }
            else if(ac[i]->scar_type == 's') {
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
            if(ac[i]->fibrotic) {
                fibs[i+1] = 0.0f;
            }
            else if (ac[i]->border_zone) {
                double center_x = ac[i]->center_x;
                double center_y = ac[i]->center_y;
                double center_z = ac[i]->center_z;
                if(ac[i]->scar_type == 'b') {
                    dist_big = sqrt((center_x - big_scar_center_x) * (center_x - big_scar_center_x) +
                                    (center_y - big_scar_center_y) * (center_y - big_scar_center_y) +
                                    (center_z - big_scar_center_z) * (center_z - big_scar_center_z));
                    fibs[i+1] = (real)(dist_big / bz_size_big);
                    
                }
                else if(ac[i]->scar_type == 's') {
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

//    #pragma omp parallel for private(dist) reduction(max: bz_size)
	#pragma omp parallel for private(dist)
    for (i = 0; i < num_active_cells; i++) {
        if(ac[i]->active) {
            if(ac[i]->border_zone) {
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
            if(ac[i]->fibrotic) {
                fibs[i+1] = 0.0;
            }
            else if(ac[i]->border_zone) {
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
