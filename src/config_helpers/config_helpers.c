//
// Created by sachetto on 22/10/17.
//

#include "config_helpers.h"
#include "../single_file_libraries/stb_ds.h"
#include <stdio.h>
#include <stdlib.h>

void report_error_on_function(const char * function, const char *error) {
    fprintf(stderr, "Error on function %s: %s.\n", function, error);
}

void report_parameter_error_on_function(const char * function, const char *parameter) {
    fprintf(stderr, "Error: User provided function %s needs a parameter named %s. Exiting!\n", function, parameter);
    exit(EXIT_FAILURE);
}

char *get_char_parameter (struct string_hash_entry *config, const char *parameter) {
    char *config_value = shget (config, parameter);
    return config_value;
}