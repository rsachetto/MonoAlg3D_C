//
// Created by sachetto on 13/10/17.
//

#include "update_monodomain_config.h"

#ifdef _MSC_VER
#include "../dlfcn-win32/dlfcn.h"
#else
#include <dlfcn.h>
#endif
#include <string.h>
#include "../utils/file_utils.h"

#include "../common_types/common_types.h"
#include "../single_file_libraries/stb_ds.h"

void init_update_monodomain_functions(struct update_monodomain_config *config) {

    char *function_name = config->config_data.function_name;

#ifdef _MSC_VER
	char *default_library = "./shared_libs/default_update_monodomain.dll";
#else
	char *default_library = "./shared_libs/libdefault_update_monodomain.so";
#endif

    if(config->config_data.library_file_path == NULL) {
        print_to_stdout_and_file("Using the default library for update monodomain\n");
        config->config_data.library_file_path = strdup(default_library);
        config->config_data.library_file_path_was_set = true;
    }
    else {
        print_to_stdout_and_file("Using %s as update monodomain lib\n", config->config_data.library_file_path);

    }

    config->config_data.handle = dlopen (config->config_data.library_file_path, RTLD_LAZY);
    if (!config->config_data.handle) {
        fputs (dlerror(), stderr);
        fprintf(stderr, "\n");
        exit(1);
    }

    if(function_name) {
        config->update_monodomain = dlsym(config->config_data.handle, function_name);
        if (dlerror() != NULL) {
            fprintf(stderr, "\n%s function not found in the provided update_monodomain library\n", function_name);
            exit(EXIT_FAILURE);
        }
    }
    else {
        fprintf(stderr, "No function name for update_monodomain library provided. Exiting!\n");
        exit(EXIT_FAILURE);
    }

}

struct update_monodomain_config* new_update_monodomain_config() {
    struct update_monodomain_config *result = (struct update_monodomain_config*) malloc(sizeof(struct update_monodomain_config));

    init_config_common_data(&(result->config_data));
    result->update_monodomain = NULL;
    return result;
}

void print_update_monodomain_config(struct update_monodomain_config* s) {

    printf("update_monodomain_function: %s\n",s->config_data.function_name);
    printf("update_monodomain_library_file: %s\n",s->config_data.library_file_path);
    printf("update_monodomain_config:\n");
    STRING_HASH_PRINT_KEY_VALUE(s->config_data.config);
}

void free_update_monodomain_config(struct update_monodomain_config* s) {
    free_config_common_data(&(s->config_data));
    free(s);
}
