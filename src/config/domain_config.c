//
// Created by sachetto on 13/10/17.
//

#include "domain_config.h"
#include <dlfcn.h>
#include <string.h>
#include "../utils/file_utils.h"

#include "../common_types/common_types.h"
#include "../single_file_libraries/stb_ds.h"

void init_domain_functions(struct domain_config *config) {

    char *function_name = config->config_data.function_name;

	char *default_function = "./shared_libs/libdefault_domains.so";

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

    if(s == NULL) {
        print_to_stdout_and_file(" No Domain configuration.\n");
        return;
    }

    print_to_stdout_and_file("Domain configuration:\n");
    print_to_stdout_and_file("Domain name: %s\n", s->domain_name);
    print_to_stdout_and_file("Domain initial Space Discretization: dx %lf um, dy %lf um, dz %lf um\n",
                             s->start_dx, s->start_dy, s->start_dz);

    if(shlen(s->config_data.config) == 1)
    {
        print_to_stdout_and_file("Domain extra parameter:\n");
    }
    else if(shlen(s->config_data.config) > 1)
    {
        print_to_stdout_and_file("Domain extra parameters:\n");
    }

    STRING_HASH_PRINT_KEY_VALUE_LOG(s->config_data.config);
}

void free_domain_config(struct domain_config* s) {

    free_config_common_data(&(s->config_data));
    free(s->domain_name);

    free(s);
}
