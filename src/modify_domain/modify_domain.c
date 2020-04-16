//
// Created by sachetto on 06/09/2019.
//

#include "../config/modify_current_domain_config.h"
#include "../libraries_common/common_data_structures.h"
#include "../config_helpers/config_helpers.h"
#include <time.h>
#include <unistd.h>

MODIFY_DOMAIN(set_scar_fibrosis) {

    real_cpu phi = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, phi, config->config_data, "phi");

    real_cpu plain_center = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, plain_center, config->config_data, "plain_center");

    real_cpu sphere_radius = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sphere_radius, config->config_data, "sphere_radius");

    real_cpu border_zone_size = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, border_zone_size, config->config_data,
                                                "border_zone_size");

    unsigned fib_seed = 0;
    bool success;
    GET_PARAMETER_NUMERIC_VALUE(unsigned, fib_seed, config->config_data, "seed", success);

    if(!success)
        fib_seed = (unsigned)time(NULL) + getpid();

    char tmp[256];
    sprintf(tmp, "%u", fib_seed);
    shput_dup_value(config->config_data, "seed", "tmp");

    srand(fib_seed);

    struct cell_node *grid_cell = the_grid->first_cell;

    while (grid_cell != 0) {

        if (grid_cell->active) {
            if (FIBROTIC(grid_cell)) {
                real_cpu p = (real_cpu) (rand()) / (RAND_MAX);
                if (p < phi)
                    grid_cell->active = false;
                grid_cell->can_change = false;
            } else if (BORDER_ZONE(grid_cell)) {
                real_cpu distance_from_center =
                        sqrt((grid_cell->center.x - plain_center) * (grid_cell->center.x - plain_center) +
                             (grid_cell->center.y - plain_center) * (grid_cell->center.y - plain_center));
                distance_from_center = (distance_from_center - sphere_radius) / border_zone_size;
                real_cpu phi_local = phi - phi * distance_from_center;
                real_cpu p = (real_cpu) (rand()) / (RAND_MAX);

                if (p < phi_local) {
                    grid_cell->active = false;
                }

                grid_cell->can_change = false;
            }
        }
        grid_cell = grid_cell->next;
    }


}
