//
// Created by sachetto on 13/10/17.
//

#ifndef MONOALG3D_DOMAIN_CONFIG_H
#define MONOALG3D_DOMAIN_CONFIG_H

#include "../alg/grid/grid.h"
#include "../hash/string_hash.h"

typedef void (*set_spatial_domain_fn_pt)(struct grid *, struct string_hash *);

struct domain_config {
    void *handle;
    char *domain_function;
    char *domain_library_file;
    char *domain_name;
    struct string_hash *config;
    set_spatial_domain_fn_pt set_spatial_domain_fn;
    bool configured;
};

void init_domain_functions(struct domain_config *config);
struct domain_config* new_domain_config();
void free_domain_config(struct domain_config* s);
void print_domain_config_values(struct domain_config* s);

#endif //MONOALG3D_STIM_CONFIG_H
