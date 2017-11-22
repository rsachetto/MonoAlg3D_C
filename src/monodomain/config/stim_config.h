//
// Created by sachetto on 13/10/17.
//

#ifndef MONOALG3D_STIM_CONFIG_H
#define MONOALG3D_STIM_CONFIG_H

#include "../../alg/grid/grid.h"
#include "../constants.h"
#include "config_common.h"

struct stim_config;

#define SET_SPATIAL_STIM(name) EXPORT_FN void name(struct stim_config *config, struct grid *the_grid)
typedef SET_SPATIAL_STIM(set_spatial_stim_fn);

struct stim_config {

    struct config_common config_data;

    real stim_start;
    bool stim_start_was_set;
    real stim_duration;
    bool stim_duration_was_set;
    real stim_current;
    bool stim_current_was_set;

    real *spatial_stim_currents;
    set_spatial_stim_fn *set_spatial_stim;
};

void init_stim_functions(struct stim_config *config, char* stim_name);
struct stim_config* new_stim_config();
void print_stim_config_values(struct stim_config* s);
void free_stim_config(struct stim_config* s);

#endif //MONOALG3D_STIM_CONFIG_H
