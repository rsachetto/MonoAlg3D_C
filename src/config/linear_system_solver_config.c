//
// Created by sachetto on 13/10/17.
//

#include "linear_system_solver_config.h"

#include <dlfcn.h>
#include <string.h>
#include "../utils/file_utils.h"
#include "../common_types/common_types.h"
#include "../single_file_libraries/stb_ds.h"

void print_linear_system_solver_config_values(struct config* s) {

    if(s == NULL) {
        print_to_stdout_and_file("[linear_system_solver] No Linear system solver configuration.\n");
        return;
    }

    print_to_stdout_and_file("[linear_system_solver] Linear system solver configuration:\n");

    print_to_stdout_and_file("[linear_system_solver] Linear system solver library: %s\n", s->library_file_path);
    print_to_stdout_and_file("[linear_system_solver] Linear system solver main function: %s\n", s->main_function_name);

    if(s->init_function_name)
        print_to_stdout_and_file("[linear_system_solver] Linear system solver init function: %s\n", s->init_function_name);

    if(s->end_function_name)
        print_to_stdout_and_file("[linear_system_solver] Liner system solver end function: %s\n", s->end_function_name);

    STRING_HASH_PRINT_KEY_VALUE_LOG(s->config_data);
}

