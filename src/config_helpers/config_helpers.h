//
// Created by sachetto on 22/10/17.
//

#ifndef MONOALG3D_CONFIG_HELPERS_H
#define MONOALG3D_CONFIG_HELPERS_H
#include "../common_types/common_types.h"
#include "../monodomain/constants.h"
#include <stdbool.h>
#include <string.h>

#include "../single_file_libraries/stb_ds.h"

void report_parameter_error_on_function(const char *function, const char *parameter);
void report_error_on_function(const char *function, const char *error);
char *get_char_parameter(struct string_hash_entry *config, const char *parameter);

#define IS_TRUE(str) ( strcmp((str), "true") == 0 || strcmp((str), "yes") == 0 || strcmp((str), "1") == 0 )
#define IS_FALSE(str) ( strcmp((str), "false") == 0 || strcmp((str), "no") == 0 || strcmp((str), "0") == 0 )

#define GET_PARAMETER_VALUE_CHAR_OR_USE_DEFAULT(value, config, parameter)                                              \
    do {                                                                                                               \
        char *config_char = get_char_parameter(config, parameter);                                                     \
        if(config_char) {                                                                                              \
            (value) = strdup(config_char);                                                                             \
        }                                                                                                              \
    } while(0)

#define GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(value, config, parameter)                                             \
    do {                                                                                                               \
        (value) = NULL;                                                                                                \
        GET_PARAMETER_VALUE_CHAR_OR_USE_DEFAULT(value, config, parameter);                                             \
        if(!(value)) {                                                                                                 \
            report_parameter_error_on_function(__func__, parameter);                                                   \
        }                                                                                                              \
    } while(0)

#define GET_PARAMETER_NUMERIC_VALUE(type, value, config, parameter, success)                                           \
    do {                                                                                                               \
        char *config_char = get_char_parameter(config, parameter);                                                     \
        (success) = false;                                                                                             \
        if(config_char) {                                                                                              \
            (value) = (type)strtod(config_char, NULL);                                                                 \
            (success) = true;                                                                                          \
        }                                                                                                              \
    } while(0)

#define GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(value, config, parameter)                                           \
    do {                                                                                                               \
        char *value_char = NULL;                                                                                       \
        GET_PARAMETER_VALUE_CHAR_OR_USE_DEFAULT(value_char, config, parameter);                                        \
        if(value_char != NULL) {                                                                                       \
            (value) = IS_TRUE(value_char);                                                                             \
        }                                                                                                              \
        free(value_char);                                                                                              \
    } while(0)

#define GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(type, value, config, parameter)                                     \
    do {                                                                                                               \
        char *config_char = get_char_parameter(config, parameter);                                                     \
        if(config_char) {                                                                                              \
            (value) = (type)strtod(config_char, NULL);                                                                 \
        }                                                                                                              \
    } while(0)

#define GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(type, value, config, parameter)                                    \
    do {                                                                                                               \
        bool success;                                                                                                  \
        GET_PARAMETER_NUMERIC_VALUE(type, value, config, parameter, success);                                          \
        if(!success) {                                                                                                 \
            report_parameter_error_on_function(__func__, parameter);                                                   \
        }                                                                                                              \
    } while(0)

#define REPORT_ERROR_ON_FUNCTION(error)                                                                                \
    do {                                                                                                               \
        report_error_on_function(__func__, error);                                                                     \
    } while(0)

#endif // MONOALG3D_CONFIG_HELPERS_H