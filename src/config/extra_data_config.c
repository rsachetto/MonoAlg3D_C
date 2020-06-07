//
// Created by sachetto on 13/10/17.
//

#include "extra_data_config.h"
#include "../logger/logger.h"
#include "../3dparty/stb_ds.h"

void print_extra_data_config_values(struct config* s) {

    if(s == NULL) {
        log_to_stdout_and_file("[extra_data] No extra data configuration.\n");
        return;
    }
    log_to_stdout_and_file("[extra_data] Extra data ODE function configuration:\n");

    log_to_stdout_and_file("[extra_data] Extra data library: %s\n", s->library_file_path);
    log_to_stdout_and_file("[extra_data] Extra data main function: %s\n", s->main_function_name);

    if(s->init_function_name) {
        log_to_stdout_and_file("[extra_data] Extra data init function: %s\n", s->init_function_name);
    }

    if(s->end_function_name) {
        log_to_stdout_and_file("[extra_data] Extra data end function: %s\n", s->end_function_name);
    }

    if(shlen(s->config_data) == 1) {
        log_to_stdout_and_file("[extra_data] Extra data parameter:\n");
    } else if(shlen(s->config_data) > 1) {
        log_to_stdout_and_file("[extra_data] Extra data parameters:\n");
    }

    STRING_HASH_PRINT_KEY_VALUE_LOG(s->config_data);
}