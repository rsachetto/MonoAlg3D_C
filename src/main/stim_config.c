//
// Created by sachetto on 13/10/17.
//

#include "stim_config.h"
#include <dlfcn.h>
#include <string.h>

void init_stim_functions(struct stim_config *config, char* stim_name) {

    if(!config->configured) {
        printf("No stimuli configured! None will be applied!\n");
        return;
    }

    char *error;

    if(config->stim_library_file == NULL) {
        printf("Using the default library for stimuli functions for %s\n", stim_name);
        config->stim_library_file = strdup("./shared_libs/libdefault_stimuli.so");
    }
    else {
        printf("Opening %s as stimuli lib for %s\n", config->stim_library_file, stim_name);

    }

    config->handle = dlopen (config->stim_library_file, RTLD_LAZY);
    if (!config->handle) {
        fputs (dlerror(), stderr);
        fprintf(stderr, "\n");
        exit(1);
    }

    config->set_spatial_stim_fn = dlsym(config->handle, config->stim_function);
    if ((error = dlerror()) != NULL)  {
        fprintf(stderr, "\n%s function not found in the provided stimuli library\n", config->stim_function);
        exit(EXIT_FAILURE);
    }

}

struct stim_config* new_stim_config() {
    struct stim_config *result = (struct stim_config*) malloc(sizeof(struct stim_config));
    result->stim_function = NULL;
    result->stim_library_file = NULL;
    result->set_spatial_stim_fn = NULL;
    result->spatial_stim_currents = NULL;
    result->configured = false;
    return result;
}

void print_stim_config_values(struct stim_config* s) {
    printf("stim_start %lf\n", s->stim_start);
    printf("stim_duration %lf\n", s->stim_duration);
    printf("stim_current %lf\n", s->stim_current);
    printf("stim_function %s\n", s->stim_function);
}

void free_stim_config(struct stim_config* s) {
    free(s->stim_function);
    free(s->stim_library_file);
    dlclose(s->handle);
    free(s);
}
