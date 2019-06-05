//
// Created by sachetto on 13/10/17.
//

#include "extra_data_config.h"
#include "../utils/file_utils.h"
#include "../single_file_libraries/stb_ds.h"

#include <dlfcn.h>

#include <string.h>

void init_extra_data_functions(struct extra_data_config *config) {

    char *function_name = config->config_data.function_name;

	char *default_function = "./shared_libs/libdefault_extra_data.so";

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

    if(function_name) {
        config->set_extra_data = dlsym(config->config_data.handle, function_name);
        if (dlerror() != NULL) {
            fprintf(stderr, "\n%s function not found in the provided extra data library\n", function_name);
            exit(EXIT_FAILURE);
        }
    }
    else {
        fprintf(stderr, "No function name for extra data library provided. Exiting!\n");
        exit(EXIT_FAILURE);
    }

}

struct extra_data_config* new_extra_data_config() {
    struct extra_data_config *result = (struct extra_data_config*) malloc(sizeof(struct extra_data_config));
    init_config_common_data(&(result->config_data));

    result->set_extra_data = NULL;
    return result;
}

void print_extra_data_config_values(struct extra_data_config* s) {

    if(s == NULL) {
        print_to_stdout_and_file("No extra data configuration.\n");
        return;
    }


    print_to_stdout_and_file("Extra data ODE function configuration:\n");

    print_to_stdout_and_file("Extra data library: %s\n", s->config_data.library_file_path);
    print_to_stdout_and_file("Extra data function: %s\n", s->config_data.function_name);

    if(shlen(s->config_data.config) == 1) {
        print_to_stdout_and_file("Extra data parameter:\n");
    } else if(shlen(s->config_data.config) > 1) {
        print_to_stdout_and_file("Extra data parameters:\n");
    }

    STRING_HASH_PRINT_KEY_VALUE_LOG(s->config_data.config);
}

void free_extra_data_config(struct extra_data_config* s) {
    free_config_common_data(&(s->config_data));
    free(s);
}
