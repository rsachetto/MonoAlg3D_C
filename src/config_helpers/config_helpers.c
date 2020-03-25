//
// Created by sachetto on 22/10/17.
//

#include "config_helpers.h"
#include "../3dparty/stb_ds.h"
#include <stdio.h>
#include <stdlib.h>

void report_error_on_function(int line, const char * file, const char *error) {
    fprintf(stderr, "Error on line %d of file %s: %s.\n", line, file, error);
}

void report_parameter_error_on_function(int line, const char * file, const char *parameter) {
    fprintf(stderr, "Error: Call on line %d of file %s needs a parameter named %s. Exiting!\n", line, file, parameter);
    exit(EXIT_FAILURE);
}

char *get_char_parameter (struct string_hash_entry *config, const char *parameter) {
    char *config_value = shget (config, parameter);
    return config_value;
}