//
// Created by sachetto on 13/10/17.
//

#include "save_mesh_config.h"

#include <dlfcn.h>
#include <string.h>
#include <stdlib.h>

#include "../utils/file_utils.h"

#include "../single_file_libraries/stb_ds.h"
#include "../config_helpers/config_helpers.h"

void print_save_mesh_config_values(struct config* s) {

    if(s == NULL) {
        print_to_stdout_and_file("No Save results configuration.\n");
        return;
    }

    real_cpu print_rate = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, print_rate, s->config_data, "print_rate");

    char *out_dir_name = NULL;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(out_dir_name, s->config_data, "output_dir");

    print_to_stdout_and_file("Save results configuration:\n");
    print_to_stdout_and_file("Print Rate = %d\n", print_rate);

    if (out_dir_name != NULL) {
        print_to_stdout_and_file("Saving simulation results to: %s\n", out_dir_name);
    }

    if (shlen(s->config_data) == 1) {
        print_to_stdout_and_file("Save mesh extra parameter:\n");
    } else if (shlen(s->config_data) > 1) {
        print_to_stdout_and_file("Save mesh extra parameters:\n");
    }

    STRING_HASH_PRINT_KEY_VALUE_LOG(s->config_data);

    free(out_dir_name);
}