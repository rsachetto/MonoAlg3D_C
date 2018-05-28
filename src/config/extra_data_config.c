//
// Created by sachetto on 13/10/17.
//

#include "extra_data_config.h"
#include "../../utils/logfile_utils.h"

#ifdef _MSC_VER
#include "../../dlfcn-win32/dlfcn.h"
#else
#include <dlfcn.h>
#endif

#include <string.h>

void init_extra_data_functions(struct extra_data_config *config) {

    char *function_name = config->config_data.function_name;

#ifdef _MSC_VER
	char *default_function = "./shared_libs/default_extra_data.dll";
#else
	char *default_function = "./shared_libs/libdefault_extra_data.so";
#endif

    if(config->config_data.library_file_path == NULL) {
        print_to_stdout_and_file("Using the default library for extra data functions\n");
        config->config_data.library_file_path = strdup(default_function);
        config->config_data.library_file_path_was_set = true;
    }
    else {
        print_to_stdout_and_file("Opening %s as ODE extra data lib\n", config->config_data.library_file_path);

    }

    config->config_data.handle = dlopen (config->config_data.library_file_path, RTLD_LAZY);
    if (!config->config_data.handle) {
        fputs (dlerror(), stderr);
        fprintf(stderr, "\n");
        exit(1);
    }

    if(function_name) {
        config->set_extra_data = dlsym(config->config_data.handle, function_name);
        if (dlerror() != NULL) {
            fprintf(stderr, "\n%s function not found in the provided extra data library\n", function_name);
            exit(EXIT_FAILURE);
        }
    }
    else {
        fprintf(stderr, "No function name for extra data library provided. Exiting!\n");
        exit(EXIT_FAILURE);
    }

}

struct extra_data_config* new_extra_data_config() {
    struct extra_data_config *result = (struct extra_data_config*) malloc(sizeof(struct extra_data_config));
    init_config_common_data(&(result->config_data));

    result->set_extra_data = NULL;
    return result;
}

void print_extra_data_config_values(struct extra_data_config* s) {

    printf("extra_data_function: %s\n",s->config_data.function_name);
    printf("extra_data_library_file: %s\n",s->config_data.library_file_path);
    printf("extra_data_config:\n");

    STRING_HASH_PRINT_KEY_VALUE(s->config_data.config);
}

void free_extra_data_config(struct extra_data_config* s) {
    free(s->config_data.library_file_path);
    free(s->config_data.function_name);
    string_hash_destroy(s->config_data.config);
    if(s->config_data.handle)
        dlclose(s->config_data.handle);
    free(s);
}
