//
// Created by sachetto on 13/10/17.
//

#include "linear_system_solver_config.h"

#include <string.h>
#include "../logger/logger.h"
#include "../3dparty/stb_ds.h"

void print_linear_system_solver_config_values(struct config* s) {

    if(s == NULL) {
        log_to_stdout_and_file("[linear_system_solver] No Linear system solver configuration.\n");
        return;
    }

    log_to_stdout_and_file("[linear_system_solver] Linear system solver configuration:\n");

    log_to_stdout_and_file("[linear_system_solver] Linear system solver library: %s\n", s->library_file_path);
    log_to_stdout_and_file("[linear_system_solver] Linear system solver main function: %s\n", s->main_function_name);

    if(s->init_function_name)
        log_to_stdout_and_file("[linear_system_solver] Linear system solver init function: %s\n", s->init_function_name);

    if(s->end_function_name)
        log_to_stdout_and_file("[linear_system_solver] Liner system solver end function: %s\n", s->end_function_name);

    STRING_HASH_PRINT_KEY_VALUE_LOG(s->config_data);
}