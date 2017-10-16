//
// Created by sachetto on 13/10/17.
//

#include "domain_config.h"
#include <dlfcn.h>
#include <string.h>

void init_domain_functions(struct domain_config *config) {

    if(!config->config_data.configured) {
        printf("No domain function provided! Exiting!\n");
        exit(EXIT_FAILURE);
    }

    char *error;
    char *function_name = config->config_data.function_name;
    char *default_function = "./shared_libs/libdefault_domains.so";

    if(config->config_data.library_file_path == NULL) {
        printf("Using the default library for domain functions for %s\n", config->domain_name);
        config->config_data.library_file_path = strdup(default_function);
    }
    else {
        printf("Opening %s as stimuli lib\n", config->config_data.library_file_path);

    }

    config->config_data.handle = dlopen (config->config_data.library_file_path, RTLD_LAZY);
    if (!config->config_data.handle) {
        fputs (dlerror(), stderr);
        fprintf(stderr, "\n");
        exit(1);
    }

    if(function_name){
        config->set_spatial_domain_fn = dlsym(config->config_data.handle, function_name);
        if ((error = dlerror()) != NULL)  {
            fprintf(stderr, "\n%s function not found in the provided domain library\n", function_name);
            exit(EXIT_FAILURE);
        }
    }
    else {
        fprintf(stderr, "No function name for domain library provided. Exiting!\n");
        exit(EXIT_FAILURE);
    }

}

struct domain_config* new_domain_config() {
    struct domain_config *result = (struct domain_config*) malloc(sizeof(struct domain_config));

    init_config_common_data(&(result->config_data));

    result->set_spatial_domain_fn = NULL;
    result->domain_name = NULL;
    return result;
}

void print_domain_config_values(struct domain_config* s) {

    printf("domain_name: %s\n",s->domain_name);
    printf("domain_function: %s\n",s->config_data.function_name);
    printf("domain_library_file: %s\n",s->config_data.library_file_path);
    printf("start_discretization: %lf\n",s->start_h);
    printf("maximum_discretization: %lf\n",s->max_h);
    printf("domain_config:\n");
    STRING_HASH_PRINT_KEY_VALUE(s->config_data.config);
}

void free_domain_config(struct domain_config* s) {
    free(s->config_data.library_file_path);
    free(s->config_data.function_name);
    string_hash_destroy(s->config_data.config);
    free(s->domain_name);
    if(s->config_data.handle)
        dlclose(s->config_data.handle);
    free(s);
}
