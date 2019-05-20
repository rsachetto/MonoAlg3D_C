//
// Created by sachetto on 13/10/17.
//

#include "stim_config.h"
#include <dlfcn.h>
#include <string.h>
#include "../utils/file_utils.h"
#include "../common_types/common_types.h"

#include "../single_file_libraries/stb_ds.h"

void init_stim_functions(struct stim_config *config, char* stim_name) {

	char *default_library_name = "./shared_libs/libdefault_stimuli.so";

    char *function_name = config->config_data.function_name;

    if(config->config_data.library_file_path == NULL) {
        config->config_data.library_file_path = strdup(default_library_name);
        config->config_data.library_file_path_was_set = true;
    }

    config->config_data.handle = dlopen (config->config_data.library_file_path, RTLD_LAZY);
    if (!config->config_data.handle) {
        fputs (dlerror(), stderr);
        fprintf(stderr, "\n");
        exit(1);
    }

    if(function_name) {
        config->set_spatial_stim = dlsym(config->config_data.handle, function_name);
        if (dlerror() != NULL) {
            fprintf(stderr, "\n%s function not found in the provided stimuli library\n", function_name);
            exit(EXIT_FAILURE);
        }
    }
    else {
        fprintf(stderr, "No function name for sitimuli library provided. Exiting!\n");
        exit(EXIT_FAILURE);
    }


}

struct stim_config* new_stim_config() {

    struct stim_config *result = (struct stim_config*) malloc(sizeof(struct stim_config));
    init_config_common_data(&(result->config_data));

    result->set_spatial_stim = NULL;
    result->spatial_stim_currents = NULL;
    result->stim_period = 0.0;

    result->stim_current_was_set = false;
    result->stim_duration_was_set = false;
    result->stim_start_was_set = false;
    result->stim_period_was_set = false;

    return result;
}

void print_stim_config_values(struct stim_config* s) {


    if(s == NULL) {
        print_to_stdout_and_file("No Stimulus configuration.\n");
        return;
    }

    print_to_stdout_and_file("Stimulus start: %lf\n",    s->stim_start);
    print_to_stdout_and_file("Stimulus duration: %lf\n", s->stim_duration);
    print_to_stdout_and_file("Stimulus current: %lf\n",  s->stim_current);
    print_to_stdout_and_file("Stimulus library: %s\n",   s->config_data.library_file_path);
    print_to_stdout_and_file("Stimulus function: %s\n",  s->config_data.function_name);

    struct string_hash_entry *tmp = s->config_data.config;

    if(shlen(tmp) == 1)
    {
        print_to_stdout_and_file("Stimulus extra parameter:\n");
    }
    else if(shlen(tmp)  > 1)
    {
        print_to_stdout_and_file("Stimulus extra parameters:\n");
    }

    STRING_HASH_PRINT_KEY_VALUE_LOG(tmp);

}

void free_stim_config(struct stim_config* s) {
    free_config_common_data(&(s->config_data));
    free(s);
}
