//
// Created by sachetto on 13/10/17.
//

#ifndef MONOALG3D_MODIFY_DOMAIN_CONFIG_H
#define MONOALG3D_MODIFY_DOMAIN_CONFIG_H

#include "../alg/grid/grid.h"
#include "../monodomain/constants.h"
#include "config_common.h"

#define MODIFY_DOMAIN(name) void name(struct time_info *time_info, struct config *config, struct grid *the_grid)
typedef MODIFY_DOMAIN(modify_current_domain_fn);

void print_modify_domain_config_values(struct config *s);

#endif // MONOALG3D_MODIFY_DOMAIN_CONFIG_H
