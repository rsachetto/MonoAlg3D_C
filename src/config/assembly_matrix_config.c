//
// Created by sachetto on 13/10/17.
//

#include "assembly_matrix_config.h"

#include <dlfcn.h>
#include <string.h>
#include "../utils/file_utils.h"

#include "../common_types/common_types.h"
#include "../single_file_libraries/stb_ds.h"

void init_assembly_matrix_functions(struct assembly_matrix_config *config) {

    char *function_name = config->config_data.function_name;

    char *default_function = "./shared_libs/libdefault_matrix_assembly.so";

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
        config->assembly_matrix = dlsym(config->config_data.handle, function_name);
        if (dlerror() != NULL)  {
            fprintf(stderr, "\n%s function not found in the provided assembly_matrix library\n", function_name);
            exit(EXIT_FAILURE);
        }
    }
    else {
        fprintf(stderr, "No function name for matrix assembly library provided. Exiting!\n");
        exit(EXIT_FAILURE);
    }

}

struct assembly_matrix_config* new_assembly_matrix_config() {
    struct assembly_matrix_config *result = (struct assembly_matrix_config*) malloc(sizeof(struct assembly_matrix_config));

    init_config_common_data(&(result->config_data));
    result->assembly_matrix = NULL;
    return result;
}

void print_assembly_matrix_config_values(struct assembly_matrix_config* s) {

    if(s == NULL) {
        print_to_stdout_and_file("No Assembly Matrix configuration.\n");
        return;
    }

    print_to_stdout_and_file("Assembly Matrix function configuration:\n");

    print_to_stdout_and_file("Assembly Matrix library: %s\n", s->config_data.library_file_path);
    print_to_stdout_and_file("Assembly Matrix function: %s\n", s->config_data.function_name);

    if(shlen(s->config_data.config) == 1) {
        print_to_stdout_and_file("Assembly Matrix parameter:\n");
    } else if(shlen(s->config_data.config) > 1) {
        print_to_stdout_and_file("Assembly Matrix parameters:\n");
    }

    STRING_HASH_PRINT_KEY_VALUE_LOG(s->config_data.config);
}

void free_assembly_matrix_config(struct assembly_matrix_config* s) {
    free_config_common_data(&(s->config_data));
    free(s);
}
