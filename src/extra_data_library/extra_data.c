//
// Created by sachetto on 01/10/17.
//

#include "extra_data.h"

#include "../hash/point_hash.h"
#include "../utils/utils.h"
#include <math.h>

void * set_extra_data_for_fibrosis_sphere(struct grid *the_grid, struct string_hash *extra_data_config, size_t *extra_data_bytes) {

    uint32_t num_active_cells = the_grid->num_active_cells;

    *extra_data_bytes = sizeof(float)*(num_active_cells+1);

    float *fibs = (float*)malloc(*extra_data_bytes);

    struct cell_node ** ac = the_grid->active_cells;

    char *config_char = string_hash_search(extra_data_config, "atpi");
    float atpi = (float)atof(config_char);
    free(config_char);

    fibs[0] = atpi;

    config_char = string_hash_search(extra_data_config, "plain_center");
    float plain_center = (float)atof(config_char);
    free(config_char);

    config_char = string_hash_search(extra_data_config, "bz_size");
    float bz_size = (float)atof(config_char);
    free(config_char);

    config_char = string_hash_search(extra_data_config, "fib_radius");
    float fib_radius = (float)atof(config_char);
    free(config_char);


    #pragma omp parallel for
    for (int i = 0; i < num_active_cells; i++) {

        bool stim;

        if(ac[i]->fibrotic) {
            fibs[i] = 0.0;
        }
        else if(ac[i]->border_zone) {

            float center_x = (float)ac[i]->center_x;
            float center_y = (float)ac[i]->center_y;
            float center_z = (float)ac[i]->center_z;

                float distanceFromCenter = sqrtf((center_x - plain_center)*(center_x - plain_center) + (center_y - plain_center)*(center_y - plain_center));
                distanceFromCenter = (distanceFromCenter - fib_radius)/bz_size;
                fibs[i] = distanceFromCenter;

        }
        else {
            fibs[i] = 1.0;
        }

    }


    return (void*)fibs;
}

void * set_extra_data_for_fibrosis_plain(struct grid *the_grid, struct string_hash *extra_data_config, size_t *extra_data_bytes) {

    uint32_t num_active_cells = the_grid->num_active_cells;

    *extra_data_bytes = sizeof(float)*(num_active_cells+1);

    float *fibs = (float*)calloc(num_active_cells+1, sizeof(float));

    struct cell_node ** ac = the_grid->active_cells;
    char *config_char = string_hash_search(extra_data_config, "atpi");
    float atpi = (float)atof(config_char);
    free(config_char);

    fibs[0] = atpi;

    return (void*)fibs;
}
