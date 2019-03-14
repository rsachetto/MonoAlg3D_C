//
// Created by bergolho on 19/07/18.
//

#include "purkinje_config.h"
#ifdef _MSC_VER
#include "../dlfcn-win32/dlfcn.h"
#else
#include <dlfcn.h>
#endif
#include <string.h>
#include "../utils/file_utils.h"
#include "../common_types/common_types.h"
#include "../single_file_libraries/stb_ds.h"

void init_purkinje_functions (struct purkinje_config *config) 
{

    char *function_name = config->config_data.function_name;

#ifdef _MSC_VER
	char *default_function = "./shared_libs/default_purkinje.dll";
#else
	char *default_function = "./shared_libs/libdefault_purkinje.so";
#endif

    if(config->config_data.library_file_path == NULL) 
    {
        print_to_stdout_and_file("Using the default library for Purkinje functions for %s\n", config->domain_name);
        config->config_data.library_file_path = strdup(default_function);
        config->config_data.library_file_path_was_set = true;
    }
    else 
    {
        print_to_stdout_and_file("Opening %s as Purkinje lib\n", config->config_data.library_file_path);
    }

    config->config_data.handle = dlopen (config->config_data.library_file_path, RTLD_LAZY);
    if (!config->config_data.handle) 
    {
        fputs (dlerror(), stderr);
        fprintf(stderr, "\n");
        exit(1);
    }

    if(function_name)
    {
        config->set_spatial_purkinje = dlsym(config->config_data.handle, function_name);
        if (dlerror() != NULL)  
        {
            fprintf(stderr, "\n%s function not found in the provided Purkinje library\n", function_name);
            exit(EXIT_FAILURE);
        }
    }
    else 
    {
        fprintf(stderr, "No function name for Purkinje library provided. Exiting!\n");
        exit(EXIT_FAILURE);
    }

}

struct purkinje_config* new_purkinje_config() 
{
    struct purkinje_config *result = (struct purkinje_config*) malloc(sizeof(struct purkinje_config));

    init_config_common_data(&(result->config_data));

    result->start_h_was_set = false;
    result->domain_name_was_set = false;

    result->set_spatial_purkinje = NULL;
    result->domain_name = NULL;
    return result;
}

void print_purkinje_config_values (struct purkinje_config* s) 
{

    printf("purkinje_name: %s\n",s->domain_name);
    printf("purkinje_function: %s\n",s->config_data.function_name);
    printf("purkinje_library_file: %s\n",s->config_data.library_file_path);
    printf("start_discretization: %lf\n",s->start_h);
    printf("purkinje_config:\n");
    STRING_HASH_PRINT_KEY_VALUE(s->config_data.config);
}

void free_purkinje_config(struct purkinje_config* s) 
{
    free_config_common_data(&(s->config_data));
    free(s->domain_name);
    free(s);
}
