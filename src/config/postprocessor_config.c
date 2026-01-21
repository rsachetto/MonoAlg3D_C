//
// Created by sachetto on 13/10/17.
//

#include "postprocessor_config.h"
#include "../3dparty/stb_ds.h"
#include "../logger/logger.h"
#include <stdlib.h>

#include <dlfcn.h>

struct postprocess_function *new_postprocess_function() {

    struct postprocess_function *ret = (struct postprocess_function *)malloc(sizeof(struct postprocess_function));
    ret->function_params = alloc_and_init_config_data();
    ret->function = NULL;

    return ret;
}

void init_postprocess_function(struct postprocess_function *function, const char *f_name) {

    if(!f_name)
        log_error_and_exit("No function name for [%s] provided. Exiting!\n", f_name);

    void *handle = dlopen(DEFAULT_POSTPROCESSOR_LIB, RTLD_LAZY);
    if(!handle) {
        log_error_and_exit("%s\n", dlerror());
    }

    function->function = (postprocess_fn *)dlsym(handle, f_name);
    char *error = dlerror();

    if(error != NULL) {
        log_error_and_exit("\n%s function not found in the provided in library %s. Error from dlsym %s\n", f_name, DEFAULT_POSTPROCESSOR_LIB, error);
    }
}
