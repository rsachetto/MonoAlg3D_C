//
// Created by sachetto on 13/10/17.
//

#include <unitypes.h>
#include "stimuli.h"
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