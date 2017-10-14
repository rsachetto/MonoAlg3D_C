//
// Created by sachetto on 13/10/17.
//

#include "extra_data_config.h"
#include <dlfcn.h>

void init_extra_data_functions(struct extra_data_config *config) {

    if(!config->configured) {
        printf("No extra data function configured!\n");
        return;
    }

    char *error;

    if(config->extra_data_library_file == NULL) {
        fprintf(stderr, "extra_data_library_path not provided. Exiting!\n");
        exit(1);
    }

    printf("Opening %s as ODE extra data lib\n", config->extra_data_library_file);


    config->handle = dlopen (config->extra_data_library_file, RTLD_LAZY);
    if (!config->handle) {
        fputs (dlerror(), stderr);
        fprintf(stderr, "\n");
        exit(1);
    }

    config->set_extra_data_fn = dlsym(config->handle, config->extra_data_function);
    if ((error = dlerror()) != NULL)  {
        fprintf(stderr, "\n%s function not found in the provided extra data library\n", config->extra_data_function);
        exit(EXIT_FAILURE);
    }

}

struct extra_data_config* new_extra_data_config() {
    struct extra_data_config *result = (struct extra_data_config*) malloc(sizeof(struct extra_data_config));
    result->extra_data_function = NULL;
    result->extra_data_library_file = NULL;
    result->set_extra_data_fn = NULL;
    result->config = string_hash_create();
    result->configured = false;
    return result;
}

void print_extra_data_config_values(struct extra_data_config* s) {

    printf("extra_data_function: %s\n",s->extra_data_function);
    printf("extra_data_library_file: %s\n",s->extra_data_library_file);
    printf("extra_data_config:\n");

    STRING_HASH_PRINT_KEY_VALUE(s->config);
}

void free_extra_data_config(struct extra_data_config* s) {
    free(s->extra_data_library_file);
    free(s->extra_data_function);
    string_hash_destroy(s->config);
    dlclose(s->handle);
    free(s);
}
