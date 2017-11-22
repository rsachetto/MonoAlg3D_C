//
// Created by sachetto on 13/10/17.
//

#ifndef MONOALG3D_DOMAIN_CONFIG_H
#define MONOALG3D_DOMAIN_CONFIG_H

#include "../../alg/grid/grid.h"
#include "config_common.h"
#include "../constants.h"

struct domain_config;

#define SET_SPATIAL_DOMAIN(name) EXPORT_FN void name(struct domain_config *config, struct grid *the_grid)
typedef SET_SPATIAL_DOMAIN(set_spatial_domain_fn);

struct domain_config {
    struct config_common config_data;
    char *domain_name;
    bool domain_name_was_set;
    double start_h, max_h;
    bool start_h_was_set, max_h_was_set;
    set_spatial_domain_fn *set_spatial_domain;
};

void init_domain_functions(struct domain_config *config);
struct domain_config* new_domain_config();
void free_domain_config(struct domain_config* s);
void print_domain_config_values(struct domain_config* s);



#endif //MONOALG3D_STIM_CONFIG_H
