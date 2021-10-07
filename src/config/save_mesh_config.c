//
// Created by sachetto on 13/10/17.
//

#include "save_mesh_config.h"

#include <string.h>
#include <stdlib.h>

#include "../logger/logger.h"

#include "../3dparty/stb_ds.h"
#include "../config_helpers/config_helpers.h"

void print_save_mesh_config_values(struct config* s) {

    if(s == NULL) {
        log_info("[save_mesh] No Save results configuration.\n");
        return;
    }

    int print_rate = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(int, print_rate, s, "print_rate");

    char *out_dir_name = NULL;
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(out_dir_name, s, "output_dir");

    log_info("Save results configuration:\n");
    log_info("[save_mesh] Print Rate = %d\n", print_rate);

    if (out_dir_name != NULL) {
        log_info("[save_mesh] Saving simulation results to: %s\n", out_dir_name);
    }

    LOG_CONFIG_INFO("[save_mesh] Save mesh extra parameter", s->config_data);    

    free(out_dir_name);
}
