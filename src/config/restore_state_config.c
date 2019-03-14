//
// Created by sachetto on 13/10/17.
//

#include "restore_state_config.h"

#ifdef _MSC_VER
#include "../dlfcn-win32/dlfcn.h"
#else
#include <dlfcn.h>
#endif
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "../utils/file_utils.h"

#include "../common_types/common_types.h"
#include "../single_file_libraries/stb_ds.h"

void init_restore_state_functions(struct restore_state_config *config) {

    char *function_name = config->config_data.function_name;

#ifdef _MSC_VER
	char *default_library = "./shared_libs/libdefault_restore_state.dll";
#else
	char *default_library = "./shared_libs/libdefault_restore_state.so";
#endif

    if(config->config_data.library_file_path == NULL) {
        print_to_stdout_and_file("Using the default library for restoring simulation state\n");
        config->config_data.library_file_path = strdup(default_library);
        config->config_data.library_file_path_was_set = true;
    }
    else {
        print_to_stdout_and_file("Using %s as restore state lib\n", config->config_data.library_file_path);

    }

    config->config_data.handle = dlopen (config->config_data.library_file_path, RTLD_LAZY);
    if (!config->config_data.handle) {
        fputs (dlerror(), stderr);
        fprintf(stderr, "\n");
        exit(1);
    }

    if(function_name) {
        config->restore_state = dlsym(config->config_data.handle, function_name);
        if (dlerror() != NULL)  {
            fprintf(stderr, "\n%s function not found in the provided restore state library\n", function_name);
            exit(EXIT_FAILURE);
        }
    }
    else {
        fprintf(stderr, "No function name for restore state library provided. Exiting!\n");
        exit(EXIT_FAILURE);
    }

}

struct restore_state_config* new_restore_state_config() {
    struct restore_state_config *result = (struct restore_state_config*) malloc(sizeof(struct restore_state_config));

    init_config_common_data(&(result->config_data));
    result->restore_state = NULL;
    return result;
}

void print_restore_state_config_values(struct restore_state_config* s) {
    printf("restore_state_function: %s\n",s->config_data.function_name);
    printf("restore_stateh_library_file: %s\n",s->config_data.library_file_path);
    printf("restore_state_config:\n");
    //TODO: chan
    STRING_HASH_PRINT_KEY_VALUE(s->config_data.config);
}

void free_restore_state_config(struct restore_state_config* s) {
    free_config_common_data(&(s->config_data));
    free(s);

}
