//
// Created by sachetto on 13/10/17.
//

#include "save_mesh_config.h"

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

void init_save_mesh_functions(struct save_mesh_config *config) {

    char *function_name = config->config_data.function_name;

#ifdef _MSC_VER
	char *default_function = "./shared_libs/libdefault_save_mesh.dll";
#else
	char *default_function = "./shared_libs/libdefault_save_mesh.so";
#endif

    if(config->config_data.library_file_path == NULL) {
        print_to_stdout_and_file("Using the default library for saving the mesh\n");
        config->config_data.library_file_path = strdup(default_function);
        config->config_data.library_file_path_was_set = true;
    }
    else {
        print_to_stdout_and_file("Using %s as save mesh lib\n", config->config_data.library_file_path);

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
    result->last_count = 0;
    return result;
}

void print_save_mesh_config_values(struct save_mesh_config* s) {

    printf("save_mesh_print_rate: %d\n",s->print_rate);
    printf("save_mesh_output_dir: %s\n",s->out_dir_name);
    printf("save_mesh_function: %s\n",s->config_data.function_name);
    printf("save_mesh_library_file: %s\n",s->config_data.library_file_path);
    printf("save_mesh_config:\n");
    STRING_HASH_PRINT_KEY_VALUE(s->config_data.config);
}

void free_save_mesh_config(struct save_mesh_config* s) {

    free_config_common_data(&(s->config_data));
    free(s->out_dir_name);
    free(s);
}
