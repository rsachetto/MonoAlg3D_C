//
// Created by sachetto on 13/10/17.
//

#include "postprocessor_config.h"
#include <stdlib.h>
#include "../3dparty/stb_ds.h"
#include "../logger/logger.h"

#include <dlfcn.h>

struct postprocess_function *new_postprocess_function() {
	struct postprocess_function *ret = malloc(sizeof(struct postprocess_function));

    sh_new_arena(ret->function_params);
    shdefault(ret->function_params, NULL);
    ret->function = NULL;

	return ret;
}

void init_postprocess_function(struct postprocess_function *function, const char *f_name) {
	
	if(!f_name)
		log_to_stderr_and_file_and_exit("No function name for [%s] provided. Exiting!\n", f_name);
	
	void *handle = dlopen (DEFAULT_POSTPROCESSOR_LIB, RTLD_LAZY);
	if (!handle) {
		log_to_stderr_and_file_and_exit("%s\n", dlerror());
	}

	function->function = dlsym(handle, f_name);
	char *error = dlerror();
	
	if (error != NULL)  {
		log_to_stderr_and_file_and_exit("\n%s function not found in the provided in library %s. Error from dlsym %s\n", f_name, DEFAULT_POSTPROCESSOR_LIB, error);
	}

}
