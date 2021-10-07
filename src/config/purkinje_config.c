//
// Created by bergolho on 19/07/18.
//
#include "purkinje_config.h"

#include <dlfcn.h>
#include <string.h>
#include "../utils/file_utils.h"

#include "../3dparty/stb_ds.h"
#include "../config_helpers/config_helpers.h"

void print_purkinje_config_values (struct config* config) {

    if(config == NULL) {
        log_info("[purkinje] No Purkinje configuration.\n");
        return;
    }

    char *name = NULL;
    real_cpu dx = 0.0;
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(name, config, "name");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, dx, config, "dx");

    log_info("Purkinje configuration:\n");
    log_info("[purkinje] Purkinje network name: %s\n", name);
    log_info("[purkinje] Purkinje network initial Space Discretization: %lf um\n", dx);

    LOG_CONFIG_INFO("[purkinje] Purkinje extra parameter", config->config_data);

    free(name);
}
