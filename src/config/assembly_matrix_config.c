//
// Created by sachetto on 13/10/17.
//

#include "assembly_matrix_config.h"

#ifdef _MSC_VER
#include "../dlfcn-win32/dlfcn.h"
#else
#include <dlfcn.h>
#endif
#include <string.h>
#include "../utils/file_utils.h"

void init_assembly_matrix_functions(struct assembly_matrix_config *config) {

    char *function_name = config->config_data.function_name;

#ifdef _MSC_VER
	char *default_function = "./shared_libs/default_matrix_assembly.dll";
#else
	char *default_function = "./shared_libs/libdefault_matrix_assembly.so";
#endif

    if(config->config_data.library_file_path == NULL) {
        print_to_stdout_and_file("Using the default library for matrix assembly\n");
        config->config_data.library_file_path = strdup(default_function);
        config->config_data.library_file_path_was_set = true;
    }
    else {
        print_to_stdout_and_file("Using %s as matrix assembly lib\n", config->config_data.library_file_path);

    }

    config->config_data.handle = dlopen (config->config_data.library_file_path, RTLD_LAZY);
    if (!config->config_data.handle) {
        fputs (dlerror(), stderr);
        fprintf(stderr, "\n");
        exit(1);
    }

    if(function_name){
        config->assembly_matrix = dlsym(config->config_data.handle, function_name);
        if (dlerror() != NULL)  {
            fprintf(stderr, "\n%s function not found in the provided assembly_matrix library\n", function_name);
            exit(EXIT_FAILURE);
        }
    }
    else {
        fprintf(stderr, "No function name for matrix assembly library provided. Exiting!\n");
        exit(EXIT_FAILURE);
    }

}

struct assembly_matrix_config* new_assembly_matrix_config() {
    struct assembly_matrix_config *result = (struct assembly_matrix_config*) malloc(sizeof(struct assembly_matrix_config));

    init_config_common_data(&(result->config_data));
    result->assembly_matrix = NULL;
    return result;
}

void print_assembly_matrix_config(struct assembly_matrix_config* s) {

    printf("assembly_matrix_function: %s\n",s->config_data.function_name);
    printf("assembly_matrix_library_file: %s\n",s->config_data.library_file_path);
    printf("assembly_matrix_config:\n");
    STRING_HASH_PRINT_KEY_VALUE(s->config_data.config);
}

void free_assembly_matrix_config(struct assembly_matrix_config* s) {
    free(s->config_data.library_file_path);
    free(s->config_data.function_name);
    string_hash_destroy(s->config_data.config);

    if(s->config_data.handle)
        dlclose(s->config_data.handle);
    free(s);
}
