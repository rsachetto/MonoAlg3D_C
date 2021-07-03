//
// Created by sachetto on 13/10/17.
//

#include "domain_config.h"
#include <dlfcn.h>
#include <string.h>
#include "../logger/logger.h"

#include "../3dparty/stb_ds.h"
#include "../config_helpers/config_helpers.h"

void print_domain_config_values(struct config* config) {

    if(config == NULL) {
        log_info("[domain] No Domain configuration.\n");
        return;
    }

    char *name;
    real_cpu start_dx = 0;
    real_cpu start_dy = 0;
    real_cpu start_dz = 0;

    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(name, config, "name");
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, start_dx, config, "start_dx");
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, start_dy, config, "start_dy");
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, start_dz, config, "start_dz");

    log_info("Domain configuration:\n");
    log_info("[domain] Domain name: %s\n", name);

    if(start_dx && start_dy && start_dz) {
        log_info("[domain] Domain initial Space Discretization: dx %lf um, dy %lf um, dz %lf um\n", start_dx, start_dy, start_dz);
    }

    if(shlen(config->config_data) == 1)
    {
        log_info("[domain] Domain extra parameter:\n");
    }
    else if(shlen(config->config_data) > 1)
    {
        log_info("[domain] Domain extra parameters:\n");
    }

    STRING_HASH_PRINT_KEY_VALUE_LOG(config->config_data);
}
