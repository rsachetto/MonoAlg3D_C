//
// Created by sachetto on 13/10/17.
//

#include "domain_config.h"
#ifdef _MSC_VER
#include "../dlfcn-win32/dlfcn.h"
#else
#include <dlfcn.h>
#endif
#include <string.h>
#include "../utils/file_utils.h"

#include "../common_types/common_types.h"
#include "../single_file_libraries/stb_ds.h"

void init_domain_functions(struct domain_config *config) {

    char *function_name = config->config_data.function_name;

#ifdef _MSC_VER
	char *default_function = "./shared_libs/default_domains.dll";
#else
	char *default_function = "./shared_libs/libdefault_domains.so";
#endif

    if(config->config_data.library_file_path == NULL) {
        print_to_stdout_and_file("Using the default library for domain functions for %s\n", config->domain_name);
        config->config_data.library_file_path = strdup(default_function);
        config->config_data.library_file_path_was_set = true;
    }
    else {
        print_to_stdout_and_file("Using %s as domain lib\n", config->config_data.library_file_path);

    }

    config->config_data.handle = dlopen (config->config_data.library_file_path, RTLD_LAZY);
    if (!config->config_data.handle) {
        fputs (dlerror(), stderr);
        fprintf(stderr, "\n");
        exit(1);
    }

    if(function_name){
        config->set_spatial_domain = dlsym(config->config_data.handle, function_name);
        if (dlerror() != NULL)  {
            fprintf(stderr, "\n%s function not found in the provided domain library\n", function_name);
            exit(EXIT_FAILURE);
        }
    }
    else {
        fprintf(stderr, "No function name for domain library provided. Exiting!\n");
        exit(EXIT_FAILURE);
    }

}

struct domain_config* new_domain_config() {
    struct domain_config *result = (struct domain_config*) malloc(sizeof(struct domain_config));

    init_config_common_data(&(result->config_data));

    result->start_dx_was_set = false;
    result->max_dx_was_set = false;

    result->start_dy_was_set = false;
    result->max_dy_was_set = false;

    result->start_dz_was_set = false;
    result->max_dz_was_set = false;

    result->domain_name_was_set = false;

    result->set_spatial_domain = NULL;
    result->domain_name = NULL;
    return result;
}

void print_domain_config_values(struct domain_config* s) {

    printf("domain_name: %s\n",s->domain_name);
    printf("domain_function: %s\n",s->config_data.function_name);
    printf("domain_library_file: %s\n",s->config_data.library_file_path);
    printf("start_discretization: dx %lf, dy %lf, dz %lf\n", s->start_dx, s->start_dy, s->start_dz);
    printf("maximum_discretization: dx %lf, dy %lf, dz %lf\n", s->max_dx, s->max_dy, s->max_dz);
    printf("domain_config:\n");
    STRING_HASH_PRINT_KEY_VALUE(s->config_data.config);
}

void free_domain_config(struct domain_config* s) {

    free_config_common_data(&(s->config_data));
    free(s->domain_name);

    free(s);
}
