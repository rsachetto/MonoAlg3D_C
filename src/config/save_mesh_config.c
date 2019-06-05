//
// Created by sachetto on 13/10/17.
//

#include "save_mesh_config.h"

#include <dlfcn.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "../utils/file_utils.h"

#include "../common_types/common_types.h"
#include "../single_file_libraries/stb_ds.h"

void init_save_mesh_functions(struct save_mesh_config *config) {

    char *function_name = config->config_data.function_name;

	char *default_function = "./shared_libs/libdefault_save_mesh.so";

    if(config->config_data.library_file_path == NULL) {
        config->config_data.library_file_path = strdup(default_function);
        config->config_data.library_file_path_was_set = true;
    }

    config->config_data.handle = dlopen (config->config_data.library_file_path, RTLD_LAZY);
    if (!config->config_data.handle) {
        fputs (dlerror(), stderr);
        fprintf(stderr, "\n");
        exit(1);
    }

    if(function_name){
        config->save_mesh = dlsym(config->config_data.handle, function_name);
        if (dlerror() != NULL)  {
            fprintf(stderr, "\n%s function not found in the provided save_mesh library\n", function_name);
            exit(EXIT_FAILURE);
        }
    }
    else {
        fprintf(stderr, "No function name for save_mesh library provided. Exiting!\n");
        exit(EXIT_FAILURE);
    }

}

struct save_mesh_config* new_save_mesh_config() {
    struct save_mesh_config *result = (struct save_mesh_config*) malloc(sizeof(struct save_mesh_config));

    init_config_common_data(&(result->config_data));
    result->save_mesh = NULL;
    result->out_dir_name = NULL;
    result->out_dir_name_was_set = false;
    result->print_rate_was_set = false;
    result->print_rate = 0;
    result->remove_older_simulation_dir = false;
    result->remove_older_simulation_dir_was_set = false;
    return result;
}

void print_save_mesh_config_values(struct save_mesh_config* s) {

    if(s == NULL) {
        print_to_stdout_and_file("No Save results configuration.\n");
        return;
    }

    print_to_stdout_and_file("Save results configuration:\n");
    print_to_stdout_and_file("Print Rate = %d\n", s->print_rate);

    if (s->out_dir_name != NULL) {
        print_to_stdout_and_file("Saving simulation results to: %s\n", s->out_dir_name);
    }

    if (shlen(s->config_data.config) == 1) {
        print_to_stdout_and_file("Save mesh extra parameter:\n");
    } else if (shlen(s->config_data.config) > 1) {
        print_to_stdout_and_file("Save mesh extra parameters:\n");
    }

    STRING_HASH_PRINT_KEY_VALUE_LOG(s->config_data.config);
}

void free_save_mesh_config(struct save_mesh_config* s) {

    free_config_common_data(&(s->config_data));
    free(s->out_dir_name);
    free(s);
}
