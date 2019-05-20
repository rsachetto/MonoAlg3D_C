//
// Created by sachetto on 13/10/17.
//

#include "linear_system_solver_config.h"

#include <dlfcn.h>
#include <string.h>
#include "../utils/file_utils.h"
#include "../common_types/common_types.h"
#include "../single_file_libraries/stb_ds.h"

void init_linear_system_solver_functions(struct linear_system_solver_config *config) {

    char *function_name = config->config_data.function_name;

	char *default_function = "./shared_libs/libdefault_linear_system_solver.so";

    if(config->config_data.library_file_path == NULL) {
        config->config_data.library_file_path = strdup(default_function);
        config->config_data.library_file_path_was_set = true;
    }

    config->config_data.handle = dlopen (config->config_data.library_file_path, RTLD_LAZY);
    if (!config->config_data.handle) {
        fputs (dlerror(), stderr);
        fprintf(stderr, "\n");
        exit(1);
    }

    if(function_name){
        config->solve_linear_system = dlsym(config->config_data.handle, function_name);
        if (dlerror() != NULL)  {
            print_to_stderr_and_file_and_exit("\n%s function not found in the provided linear_system_solver library\n", function_name);
        }
    }
    else {
        print_to_stderr_and_file_and_exit("No function name for matrix assembly library provided. Exiting!\n");
    }

}

struct linear_system_solver_config* new_linear_system_solver_config() {
    struct linear_system_solver_config *result = (struct linear_system_solver_config*) malloc(sizeof(struct linear_system_solver_config));

    init_config_common_data(&(result->config_data));
    result->solve_linear_system = NULL;
    return result;
}

void print_linear_system_solver_config_values(struct linear_system_solver_config* s) {

    if(s == NULL) {
        print_to_stdout_and_file("No Linear system solver configuration.\n");
        return;
    }

    print_to_stdout_and_file("Linear system solver configuration:\n");

    print_to_stdout_and_file("Linear system solver library: %s\n", s->config_data.library_file_path);
    print_to_stdout_and_file("Linear system solver function: %s\n", s->config_data.function_name);

    STRING_HASH_PRINT_KEY_VALUE_LOG(s->config_data.config);
}

void free_linear_system_solver_config(struct linear_system_solver_config* s) {

    free_config_common_data(&(s->config_data));
    free(s);
}
