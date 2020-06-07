//
// Created by sachetto on 14/10/17.
//

#include <stddef.h>
#include <dlfcn.h>

#include "config_common.h"
#include "../3dparty/stb_ds.h"
#include "../logger/logger.h"

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
        log_to_stderr_and_file_and_exit("%s\n", dlerror());
    }

    if(main_function_name){
        config->main_function = dlsym(config->handle, main_function_name);
        char *error = dlerror();
        if (error != NULL)  {
            log_to_stderr_and_file_and_exit("\n%s function not found in the provided in library %s. Error from dlsym %s\n", main_function_name, config->library_file_path, error);
        }
    }
    else {
        log_to_stderr_and_file_and_exit("No function name for [%s] provided. Exiting!\n", config_type);
    }

    if(init_function_name){
        config->init_function = dlsym(config->handle, init_function_name);
        if (dlerror() != NULL)  {
            log_to_stderr_and_file_and_exit("\n%s function not found in the provided in library %s\n", init_function_name, config->library_file_path);
        }
    }

    if(end_function_name){
        config->end_function = dlsym(config->handle, end_function_name);
        if (dlerror() != NULL)  {
            log_to_stderr_and_file_and_exit("\n%s function not found in the provided in library %s\n", end_function_name, config->library_file_path);
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
