//
// Created by sachetto on 13/10/17.
//

#ifndef MONOALG3D_STIM_CONFIG_H
#define MONOALG3D_STIM_CONFIG_H

#include "../alg/grid/grid.h"
#include "../main/constants.h"

typedef void (*set_spatial_stim_fn_pt)(struct grid *, Real, Real *);

struct stim_config {
    void *handle;
    Real stim_start;
    Real stim_duration;
    Real stim_current;
    char *stim_function;
    Real *spatial_stim_currents;
    char *stim_library_file;
    set_spatial_stim_fn_pt set_spatial_stim_fn;
};

void init_stim_functions(struct stim_config *config);
struct stim_config* new_stim_config();
void print_stim_config_values(struct stim_config* s);
void free_stim_config(struct stim_config* s);

#endif //MONOALG3D_STIM_CONFIG_H
