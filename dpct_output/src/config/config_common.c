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
        log_error_and_exit("%s\n", dlerror());
    }

    if(main_function_name) {
        config->main_function = dlsym(config->handle, main_function_name);
        char *error = dlerror();
        if (error != NULL)  {
            log_error_and_exit("'%s' function not found in the provided in library %s. Error from dlsym %s\n", main_function_name, config->library_file_path, error);
        }
    } else {
        log_error_and_exit("No function name for [%s] provided. Exiting!\n", config_type);
    }

    if(init_function_name) {
        config->init_function = dlsym(config->handle, init_function_name);
        if (dlerror() != NULL)  {
            log_error_and_exit("'%s' function not found in the provided in library %s\n", init_function_name, config->library_file_path);
        }
    }
    else {
        sds new_init_function_name = sdscatfmt(sdsempty(), "init_%s", main_function_name);
        config->init_function = dlsym(config->handle, new_init_function_name);
        if (dlerror() != NULL)  {
            config->init_function = NULL;
        }
        else {
            config->init_function_name = strdup(new_init_function_name);
        }

        sdsfree(new_init_function_name);

    }

    if(end_function_name) {
        config->end_function = dlsym(config->handle, end_function_name);
        if (dlerror() != NULL)  {
            log_error_and_exit("'%s' function not found in the provided in library %s\n", end_function_name, config->library_file_path);
        }
    } else {
        sds new_end_function_name = sdscatfmt(sdsempty(), "end_%s", main_function_name);
        config->end_function = dlsym(config->handle, new_end_function_name);
        if (dlerror() != NULL)  {
            config->end_function = NULL;
        }
        else {
            config->end_function_name = strdup(new_end_function_name);
        }

        sdsfree(new_end_function_name);
    }


    for(int i = 0; i < arrlen(config->extra_function_names); i++) {
         arrput(config->extra_functions, dlsym(config->handle, config->extra_function_names[i]));
        if (dlerror() != NULL)  {
            log_error_and_exit("'%s' function not found in the provided in library %s\n", config->extra_function_names[i], config->library_file_path);
        }
    }

}

void free_config_data(struct config *cm) {

    size_t n = shlen(cm->config_data);

    for(size_t i = 0; i < n; i++) {
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
