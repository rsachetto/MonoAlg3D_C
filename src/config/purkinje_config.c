//
// Created by bergolho on 19/07/18.
//
#include "purkinje_config.h"

#include <dlfcn.h>
#include <string.h>
#include "../utils/file_utils.h"

#include "../single_file_libraries/stb_ds.h"
#include "../config_helpers/config_helpers.h"

void print_purkinje_config_values (struct config* config) {

    if(config == NULL) {
        print_to_stdout_and_file("[purkinje] No Purkinje configuration.\n");
        return;
    }

    char *name = NULL;
    real_cpu start_h = 0.0;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(name, config->config_data, "name");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, start_h, config->config_data, "start_discretization");

    print_to_stdout_and_file("Purkinje configuration:\n");
    print_to_stdout_and_file("[purkinje] Purkinje network name: %s\n", name);
    print_to_stdout_and_file ("[purkinje] Purkinje network initial Space Discretization: %lf um\n", start_h);

    if (shlen(config->config_data) == 1)
    {
        print_to_stdout_and_file ("[purkinje] Purkinje extra parameter:\n");
    }
    else if (shlen(config->config_data) > 1) {
        print_to_stdout_and_file ("[purkinje] Purkinje extra parameters:\n");
    }

    STRING_HASH_PRINT_KEY_VALUE_LOG (config->config_data);

    free(name);
}
