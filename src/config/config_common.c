//
// Created by sachetto on 14/10/17.
//

#include <stddef.h>
#include <dlfcn.h>
#include <stdio.h>

#include "config_common.h"


#include "../single_file_libraries/stb_ds.h"

struct config *alloc_and_init_config_data() {
    struct config *cm = (struct config *)calloc(1, sizeof(struct config));
    sh_new_arena(cm->config_data);
    shdefault(cm->config_data, NULL);
    cm->persistent_data = NULL;

    return cm;
}

void init_config_functions(struct config *config, char *default_lib, char *config_type) {

    char *main_function_name = config->main_function_name;
    char *init_function_name = config->init_function_name;
    char *end_function_name = config->end_function_name;

    if(config->library_file_path == NULL) {
        config->library_file_path = strdup(default_lib);
        config->library_file_path_was_set = true;
    }

    config->handle = dlopen (config->library_file_path, RTLD_LAZY);
    if (!config->handle) {
        fputs (dlerror(), stderr);
        fprintf(stderr, "\n");
        exit(1);
    }

    if(main_function_name){
        config->main_function = dlsym(config->handle, main_function_name);
        if (dlerror() != NULL)  {
            fprintf(stderr, "\n%s function not found in the provided in library %s\n", main_function_name, config->library_file_path);
            exit(EXIT_FAILURE);
        }
    }
    else {
        fprintf(stderr, "No function name for [%s] provided. Exiting!\n", config_type);
        exit(EXIT_FAILURE);
    }

    if(init_function_name){
        config->init_function = dlsym(config->handle, init_function_name);
        if (dlerror() != NULL)  {
            fprintf(stderr, "\n%s function not found in the provided in library %s\n", init_function_name, config->library_file_path);
            exit(EXIT_FAILURE);
        }
    }

    if(end_function_name){
        config->end_function = dlsym(config->handle, end_function_name);
        if (dlerror() != NULL)  {
            fprintf(stderr, "\n%s function not found in the provided in library %s\n", end_function_name, config->library_file_path);
            exit(EXIT_FAILURE);
        }
    }
}


void free_config_data(struct config *cm) {

    ptrdiff_t n = shlen(cm->config_data);

    for(long i = 0; i < n; i++) {
        free(cm->config_data[i].value);
    }

    shfree(cm->config_data);

    free(cm->library_file_path);
    free(cm->main_function_name);
    free(cm->init_function_name);
    free(cm->end_function_name);

    if(cm->handle)
        dlclose(cm->handle);

    free(cm);

}
