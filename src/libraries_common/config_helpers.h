//
// Created by sachetto on 22/10/17.
//

#ifndef MONOALG3D_CONFIG_HELPERS_H
#define MONOALG3D_CONFIG_HELPERS_H
#include "../hash/string_hash.h"
#include "../monodomain/constants.h"
#include <stdbool.h>
#include <string.h>


void report_parameter_error_on_function(const char * function, const char *parameter);
char *get_char_parameter (struct string_hash *config, const char *parameter);

#define GET_PARAMETER_VALUE_CHAR(value, config, parameter)                                                             \
    do {                                                                                                               \
        (value) = NULL;                                                                                                \
        char *config_char = get_char_parameter (config, parameter);                                                    \
        if (config_char) {                                                                                             \
            (value) = strdup (config_char);                                                                            \
            free (config_char);                                                                                        \
        }                                                                                                              \
    } while (0)

#define GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(value, config, parameter)                                             \
    do {                                                                                                               \
        GET_PARAMETER_VALUE_CHAR (value, config, parameter);                                                           \
        if (!(value)) {                                                                                                \
            report_parameter_error_on_function (__func__, parameter);                                                  \
        }                                                                                                              \
    } while (0)

#define GET_PARAMETER_NUMERIC_VALUE(type, value, config, parameter, success)                                           \
    do {                                                                                                               \
        char *config_char = get_char_parameter (config, parameter);                                                    \
        (success) = false;                                                                                             \
        if (config_char) {                                                                                             \
            (value) = (type)strtod (config_char, NULL);                                                                \
            (success) = true;                                                                                          \
            free (config_char);                                                                                        \
        }                                                                                                              \
    } while (0)

#define GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(type, value, config, parameter)                                     \
    do {                                                                                                               \
        char *config_char = get_char_parameter (config, parameter);                                                    \
        if (config_char) {                                                                                             \
            (value) = (type)strtod (config_char, NULL);                                                                \
            free (config_char);                                                                                        \
        }                                                                                                              \
    } while (0)

#define GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(type, value, config, parameter)                                    \
    do {                                                                                                               \
        bool s;                                                                                                        \
        GET_PARAMETER_NUMERIC_VALUE (type, value, config, parameter, s);                                               \
        if (!s) {                                                                                                      \
            report_parameter_error_on_function (__func__, parameter);                                                  \
        }                                                                                                              \
    } while (0)

#endif // MONOALG3D_CONFIG_HELPERS_H
