//
// Created by sachetto on 13/10/17.
//

#include "extra_data_config.h"
#include <dlfcn.h>
#include <string.h>

void init_extra_data_functions(struct extra_data_config *config) {

    if(!config->config_data.configured) {
        printf("No extra data function configured!\n");
        return;
    }

    char *error;
    char *function_name = config->config_data.function_name;

    char *default_function = "./shared_libs/libdefault_extra_data.so";

    if(config->config_data.library_file_path == NULL) {
        printf("Using the default library for extra data functions\n");
        config->config_data.library_file_path = strdup(default_function);
    }
    else {
        printf("Opening %s as ODE extra data lib\n", config->config_data.library_file_path);

    }

    config->config_data.handle = dlopen (config->config_data.library_file_path, RTLD_LAZY);
    if (!config->config_data.handle) {
        fputs (dlerror(), stderr);
        fprintf(stderr, "\n");
        exit(1);
    }

    if(function_name) {
        config->set_extra_data_fn = dlsym(config->config_data.handle, function_name);
        if ((error = dlerror()) != NULL) {
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

    result->set_extra_data_fn = NULL;
    return result;
}

void print_extra_data_config_values(struct extra_data_config* s) {

    printf("extra_data_function: %s\n",s->config_data.function_name);
    printf("extra_data_library_file: %s\n",s->config_data.library_file_path);
    printf("extra_data_config:\n");

    STRING_HASH_PRINT_KEY_VALUE(s->config_data.config);
}

void free_extra_data_config(struct extra_data_config* s) {
    free(s->config_data.library_file_path);
    free(s->config_data.function_name);
    string_hash_destroy(s->config_data.config);
    if(s->config_data.handle)
        dlclose(s->config_data.handle);
    free(s);
}
