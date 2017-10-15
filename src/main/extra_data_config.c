//
// Created by sachetto on 13/10/17.
//

#include "extra_data_config.h"
#include <dlfcn.h>

void init_extra_data_functions(struct extra_data_config *config) {

    if(!config->config_data.configured) {
        printf("No extra data function configured!\n");
        return;
    }

    char *error;
    char *library_path = config->config_data.library_file_path;
    char *function_name = config->config_data.function_name;

    if(library_path == NULL) {
        fprintf(stderr, "extra_data_library_path not provided. Ignoring!\n");
        config->config_data.configured = false;
        return;
    }

    printf("Opening %s as ODE extra data lib\n", library_path);


    config->config_data.handle = dlopen (library_path, RTLD_LAZY);
    if (!config->config_data.handle) {
        fputs (dlerror(), stderr);
        fprintf(stderr, "\n");
        exit(1);
    }

    config->set_extra_data_fn = dlsym(config->config_data.handle, function_name);
    if ((error = dlerror()) != NULL)  {
        fprintf(stderr, "\n%s function not found in the provided extra data library\n", function_name);
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
    dlclose(s->config_data.handle);
    free(s);
}
