 
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
        log_to_stdout_and_file("[domain] No Domain configuration.\n");
        return;
    }

    char *name;
    real_cpu start_dx = 0;
    real_cpu start_dy = 0;
    real_cpu start_dz = 0;

    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(name, config->config_data, "name");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, start_dx, config->config_data, "start_dx");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, start_dy, config->config_data, "start_dy");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, start_dz, config->config_data, "start_dz");


    log_to_stdout_and_file("Domain configuration:\n");
    log_to_stdout_and_file("[domain] Domain name: %s\n", name);
    log_to_stdout_and_file("[domain] Domain initial Space Discretization: dx %lf um, dy %lf um, dz %lf um\n",
                             start_dx, start_dy, start_dz);

    if(shlen(config->config_data) == 1)
    {
        log_to_stdout_and_file("[domain] Domain extra parameter:\n");
    }
    else if(shlen(config->config_data) > 1)
    {
        log_to_stdout_and_file("[domain] Domain extra parameters:\n");
    }

    STRING_HASH_PRINT_KEY_VALUE_LOG(config->config_data);
}