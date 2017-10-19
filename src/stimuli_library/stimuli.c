//
// Created by sachetto on 13/10/17.
//

#include <unitypes.h>
#include "stimuli.h"
#include "../utils/erros_helpers.h"

void set_benchmark_spatial_stim (struct grid *the_grid, Real stim_current, Real *spatial_currents, struct string_hash *config) {

    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    bool stim;

#pragma omp parallel for private(stim)
    for (int i = 0; i < n_active; i++) {

        stim = ac[i]->center_x > 5500.0;
        stim &= ac[i]->center_x < 7000.0;
        stim &= ac[i]->center_y < 1500.0;
        stim &= ac[i]->center_z < 1500.0;

        if (stim) {
            spatial_currents[i] = stim_current;
        } else {
            spatial_currents[i] = 0.0;
        }
    }
}

void stim_if_x_less_than (struct grid *the_grid, Real stim_current, Real *spatial_currents, struct string_hash *config) {

    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    char *config_char = string_hash_search(config, "x_limit");
    if(config_char == NULL) {
        report_parameter_error_on_function("stim_if_x_less_than", "x_limit");
    }
    double x_limit = atof(config_char);
    free(config_char);

    bool stim;

    #pragma omp parallel for private(stim)
    for (int i = 0; i < n_active; i++) {
        stim = ac[i]->center_x < x_limit;

        if (stim) {
            spatial_currents[i] = stim_current;
        } else {
            spatial_currents[i] = 0.0;
        }

    }
}

void set_stim_from_file(struct grid *the_grid, Real stim_current, Real *spatial_currents, struct string_hash *config) {

    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    bool stim = false;

    #pragma omp parallel for private(stim)
    for (int i = 0; i < n_active; i++) {

        double center_x = ac[i]->center_x;
        double center_y = ac[i]->center_y;
        double center_z = ac[i]->center_z;

        //TODO: impelement this
        //if (globalArgs.stim_file) {
        //    int index = inMesh(cell_stims, center_x, center_y, center_z, 0, s_size - 1);
        //    return index != -1;
       // }
    }

}

void stim_if_x_greater_equal_than (struct grid *the_grid, Real stim_current, Real *spatial_currents, struct string_hash *config) {

    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    char *config_char = string_hash_search(config, "x_limit");
    if(config_char == NULL) {
        report_parameter_error_on_function("stim_if_x_more_than", "x_limit");
    }
    double x_limit = atof(config_char);
    free(config_char);

    bool stim;

#pragma omp parallel for private(stim)
    for (int i = 0; i < n_active; i++) {
        stim = ac[i]->center_x >= x_limit;

        if (stim) {
            spatial_currents[i] = stim_current;
        } else {
            spatial_currents[i] = 0.0;
        }

    }
}