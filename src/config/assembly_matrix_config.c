//
// Created by sachetto on 13/10/17.
//

#include "assembly_matrix_config.h"

#include <dlfcn.h>
#include <string.h>
#include "../utils/file_utils.h"

#include "../common_types/common_types.h"
#include "../single_file_libraries/stb_ds.h"

void print_assembly_matrix_config_values(struct config* s) {

    if(s == NULL) {
        print_to_stdout_and_file("[assembly_matrix] No Assembly Matrix configuration.\n");
        return;
    }

    print_to_stdout_and_file("[assembly_matrix] Assembly Matrix function configuration:\n");

    print_to_stdout_and_file("[assembly_matrix] Assembly Matrix library: %s\n", s->library_file_path);
//    print_to_stdout_and_file("Assembly Matrix main function: %s\n", s->main_function_name);
//
//    if(s->init_function_name) {
//        print_to_stdout_and_file("Set initial condition function: %s\n", s->init_function_name);
//    }
//
//    if(s->init_function_name) {
//        print_to_stdout_and_file("Assembly Matrix end function: %s\n", s->end_function_name);
//    }

    if(shlen(s->config_data) == 1) {
        print_to_stdout_and_file("[assembly_matrix] Assembly Matrix parameter:\n");
    } else if(shlen(s->config_data) > 1) {
        print_to_stdout_and_file("[assembly_matrix] Assembly Matrix parameters:\n");
    }

    STRING_HASH_PRINT_KEY_VALUE_LOG(s->config_data);
}