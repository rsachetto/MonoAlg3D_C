//
// Created by sachetto on 13/10/17.
//

#include "domain_config.h"
#include <dlfcn.h>
#include <string.h>

void init_domain_functions(struct domain_config *config) {

    char *error;

    if(config->domain_library_file == NULL) {
        printf("Using the default library for domain functions for %s\n", config->domain_name);
        config->domain_library_file = strdup("./shared_libs/libdefault_domains.so");
    }
    else {
        printf("Opening %s as stimuli lib\n", config->domain_library_file);

    }

    config->handle = dlopen (config->domain_library_file, RTLD_LAZY);
    if (!config->handle) {
        fputs (dlerror(), stderr);
        fprintf(stderr, "\n");
        exit(1);
    }

    config->set_spatial_domain_fn = dlsym(config->handle, config->domain_function);
    if ((error = dlerror()) != NULL)  {
        fprintf(stderr, "\n%s function not found in the provided domain library\n", config->domain_function);
        exit(EXIT_FAILURE);
    }

}

struct domain_config* new_domain_config() {
    struct domain_config *result = (struct domain_config*) malloc(sizeof(struct domain_config));
    result->domain_function = NULL;
    result->domain_library_file = NULL;
    result->set_spatial_domain_fn = NULL;
    result->domain_name = NULL;
    result->config = string_hash_create();
    return result;
}

void print_domain_config_values(struct domain_config* s) {

    printf("domain_name: %s\n",s->domain_name);
    printf("domain_function: %s\n",s->domain_function);
    printf("domain_library_file: %s\n",s->domain_library_file);
    printf("domain_config:\n");

    STRING_HASH_PRINT_KEY_VALUE(s->config);
}

void free_domain_config(struct domain_config* s) {
    free(s->domain_library_file);
    free(s->domain_function);
    string_hash_destroy(s->config);
    free(s->domain_name);
    dlclose(s->handle);
    free(s);
}
