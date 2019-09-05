//
// Created by sachetto on 13/10/17.
//

#include "stim_config.h"
#include <dlfcn.h>
#include "../utils/file_utils.h"

#include "../single_file_libraries/stb_ds.h"
#include "../config_helpers/config_helpers.h"

void print_stim_config_values(struct config* s) {

    if(s == NULL) {
        print_to_stdout_and_file("No Stimulus configuration.\n");
        return;
    }
    
    real stim_start = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_start, s->config_data, "start");

    real stim_duration = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_duration, s->config_data, "duration");

    real stim_current = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_current, s->config_data, "current");

    print_to_stdout_and_file("Stimulus start: %lf\n",    stim_start);
    print_to_stdout_and_file("Stimulus duration: %lf\n", stim_duration);
    print_to_stdout_and_file("Stimulus current: %lf\n",  stim_current);
    print_to_stdout_and_file("Stimulus library: %s\n",   s->library_file_path);
    print_to_stdout_and_file("Main Stimulus function: %s\n",  s->main_function_name);

    if(s->init_function_name)
        print_to_stdout_and_file("Init Stimulus function: %s\n",  s->init_function_name);

    if(s->end_function_name)
        print_to_stdout_and_file("End Stimulus function: %s\n",  s->end_function_name);

    struct string_hash_entry *tmp = s->config_data;

    if(shlen(tmp) == 1)
    {
        print_to_stdout_and_file("Stimulus extra parameter:\n");
    }
    else if(shlen(tmp)  > 1)
    {
        print_to_stdout_and_file("Stimulus extra parameters:\n");
    }

    STRING_HASH_PRINT_KEY_VALUE_LOG(tmp);

}