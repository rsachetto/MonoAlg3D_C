//
// Created by sachetto on 13/10/17.
//

#include "save_state_config.h"

#ifdef _MSC_VER
#include "../dlfcn-win32/dlfcn.h"
#else
#include <dlfcn.h>
#endif
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "../utils/file_utils.h"
#include "../single_file_libraries/stb_ds.h"

void init_save_state_functions(struct save_state_config *config) {

    char *function_name = config->config_data.function_name;

#ifdef _MSC_VER
	char *default_library = "./shared_libs/libdefault_save_state.dll";
#else
	char *default_library = "./shared_libs/libdefault_save_state.so";
#endif

    if(config->config_data.library_file_path == NULL) {
        print_to_stdout_and_file("Using the default library for saving simulation state\n");
        config->config_data.library_file_path = strdup(default_library);
        config->config_data.library_file_path_was_set = true;
    }
    else {
        print_to_stdout_and_file("Using %s as save state lib\n", config->config_data.library_file_path);

    }

    config->config_data.handle = dlopen (config->config_data.library_file_path, RTLD_LAZY);
    if (!config->config_data.handle) {
        fputs (dlerror(), stderr);
        fprintf(stderr, "\n");
        exit(1);
    }

    if(function_name) {
        config->save_state = dlsym(config->config_data.handle, function_name);
        if (dlerror() != NULL)  {
            fprintf(stderr, "\n%s function not found in the provided save state library\n", function_name);
            exit(EXIT_FAILURE);
        }
    }
    else {
        fprintf(stderr, "No function name for save state library provided. Exiting!\n");
        exit(EXIT_FAILURE);
    }

}

struct save_state_config* new_save_state_config() {
    struct save_state_config *result = (struct save_state_config*) malloc(sizeof(struct save_state_config));

    init_config_common_data(&(result->config_data));
    result->save_state = NULL;    
    result->save_rate = 0;   

    return result;
}

void print_save_state_config_values(struct save_state_config* s) {        
    printf("save_state_save_rate: %d\n",s->save_rate);
    printf("save_state_function: %s\n",s->config_data.function_name);
    printf("save_stateh_library_file: %s\n",s->config_data.library_file_path);
    printf("save_state_config:\n");
    STRING_HASH_PRINT_KEY_VALUE(s->config_data.config);
}

void free_save_state_config(struct save_state_config* s) {
    free_config_common_data(&(s->config_data));
    free(s);

}
