//
// Created by bergolho on 19/07/18.
//

#ifndef MONOALG3D_PURKINJE_CONFIG_H
#define MONOALG3D_PURKINJE_CONFIG_H

#include "../alg/grid/grid.h"
#include "config_common.h"
#include "../monodomain/constants.h"

struct purkinje_config;

#define SET_SPATIAL_PURKINJE(name) EXPORT_FN int name(struct purkinje_config *config, struct grid *the_grid)
typedef SET_SPATIAL_PURKINJE(set_spatial_purkinje_fn);

struct purkinje_config 
{
    struct config_common config_data;
    char *domain_name;
    bool domain_name_was_set;
    real_cpu start_h;
    bool start_h_was_set;
    set_spatial_purkinje_fn *set_spatial_purkinje;
};

void init_purkinje_functions(struct purkinje_config *config);
struct purkinje_config* new_purkinje_config();
void free_purkinje_config(struct purkinje_config* s);
void print_purkinje_config_values(struct purkinje_config* s);



#endif //MONOALG3D_STIM_CONFIG_H
