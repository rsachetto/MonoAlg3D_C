//
// Created by sachetto on 13/10/17.
//

#include "stim_config.h"
#include <dlfcn.h>
#include "../logger/logger.h"

#include "../3dparty/stb_ds.h"
#include "../config_helpers/config_helpers.h"

void print_stim_config_values(struct config* s) {

    if(s == NULL) {
        log_info("No Stimulus configuration.\n");
        return;
    }

    real stim_start = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_start, s, "start");

    real stim_duration = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_duration, s, "duration");

    real stim_current = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_current, s, "current");

    log_info("[stim] Stimulus start: %lf\n",    stim_start);
    log_info("[stim] Stimulus duration: %lf\n", stim_duration);
    log_info("[stim] Stimulus current: %lf\n",  stim_current);
    log_info("[stim] Stimulus library: %s\n",   s->library_file_path);
    log_info("[stim] Main Stimulus function: %s\n",  s->main_function_name);

    if(s->init_function_name)
        log_info("[stim] Init Stimulus function: %s\n",  s->init_function_name);

    if(s->end_function_name)
        log_info("[stim] End Stimulus function: %s\n",  s->end_function_name);

    struct string_hash_entry *tmp = s->config_data;

    LOG_CONFIG_INFO("[stim] Stimulus extra parameter", tmp);

}
