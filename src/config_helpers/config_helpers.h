//
// Created by sachetto on 22/10/17.
//

#ifndef MONOALG3D_CONFIG_HELPERS_H
#define MONOALG3D_CONFIG_HELPERS_H
#include "../3dparty/tinyexpr/tinyexpr.h"
#include "../common_types/common_types.h"
#include "../monodomain/constants.h"
#include <stdbool.h>
#include <string.h>

#include "../3dparty/stb_ds.h"

void report_parameter_error_on_function(int line, const char *file, const char *parameter);
void report_error_on_function(int line, const char *file, const char *error);
char *get_string_parameter(struct string_hash_entry *config, const char *parameter);
bool get_vector_parameter(real_cpu **v, struct string_hash_entry *config, const char *parameter, int n);

#define STRINGS_EQUAL(str1, str2) (strcmp((str1), (str2)) == 0)

#define IS_TRUE(str) (strcmp((str), "true") == 0 || strcmp((str), "yes") == 0 || strcmp((str), "1") == 0)
#define IS_FALSE(str) (strcmp((str), "false") == 0 || strcmp((str), "no") == 0 || strcmp((str), "0") == 0)

#define ARRAY_FOR_EACH(array) for(int i = 0; i < arrlen((array)); i++)

#define REPORT_ERROR_ON_FUNCTION(error)                                                                                                                        \
    do {                                                                                                                                                       \
        report_error_on_function(__LINE__, __FILE__, error);                                                                                                   \
    } while(0)

#define GET_PARAMETER_STRING_VALUE_OR_USE_DEFAULT(value, config, parameter)                                                                                    \
    do {                                                                                                                                                       \
        char *__config_char = get_string_parameter(config, parameter);                                                                                           \
        if(__config_char) {                                                                                                                                    \
            (value) = strdup(__config_char);                                                                                                                   \
        }                                                                                                                                                      \
    } while(0)

#define GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(value, config, parameter)                                                                                   \
    do {                                                                                                                                                       \
        (value) = NULL;                                                                                                                                        \
        GET_PARAMETER_STRING_VALUE_OR_USE_DEFAULT(value, config, parameter);                                                                                   \
        if(!(value)) {                                                                                                                                         \
            report_parameter_error_on_function(__LINE__, __FILE__, parameter);                                                                                 \
        }                                                                                                                                                      \
    } while(0)

#define GET_PARAMETER_NUMERIC_VALUE(type, value, config, parameter, __success)                                                                                 \
    do {                                                                                                                                                       \
        char *__config_char = get_string_parameter(config, parameter);                                                                                           \
        (__success) = false;                                                                                                                                   \
        if(__config_char) {                                                                                                                                    \
            int expr_parse_error__;                                                                                                                            \
            real_cpu expr_parse_result__ = (real_cpu)te_interp(__config_char, &expr_parse_error__);                                                            \
            (value) = (type)expr_parse_result__;                                                                                                               \
            if(expr_parse_error__ == 0)                                                                                                                        \
                (__success) = true;                                                                                                                            \
        }                                                                                                                                                      \
    } while(0)

#define GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(value, config, parameter)                                                                                   \
    do {                                                                                                                                                       \
        char *value_char = NULL;                                                                                                                               \
        GET_PARAMETER_STRING_VALUE_OR_USE_DEFAULT(value_char, config, parameter);                                                                              \
        if(value_char != NULL) {                                                                                                                               \
            (value) = IS_TRUE(value_char);                                                                                                                     \
        }                                                                                                                                                      \
        free(value_char);                                                                                                                                      \
    } while(0)

#define GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(type, value, config, parameter)                                                                             \
    do {                                                                                                                                                       \
        char *__config_char = get_string_parameter(config, parameter);                                                                                           \
        if(__config_char) {                                                                                                                                    \
            int expr_parse_error__;                                                                                                                            \
            real_cpu expr_parse_result__ = (real_cpu)te_interp(__config_char, &expr_parse_error__);                                                            \
            (value) = (type)expr_parse_result__;                                                                                                               \
        }                                                                                                                                                      \
    } while(0)

#define GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(type, value, config, parameter)                                                                            \
    do {                                                                                                                                                       \
        bool __success;                                                                                                                                        \
        GET_PARAMETER_NUMERIC_VALUE(type, value, config, parameter, __success);                                                                                \
        if(!__success) {                                                                                                                                       \
            report_parameter_error_on_function(__LINE__, __FILE__, parameter);                                                                                 \
        }                                                                                                                                                      \
    } while(0)

#define GET_PARAMETER_VECTOR_VALUE_OR_USE_DEFAULT(value, config, parameter, n)                                                                                 \
    do {                                                                                                                                                       \
        char *__config_char = get_string_parameter(config, parameter);                                                                                           \
        if(__config_char) {                                                                                                                                    \
            bool __success = get_vector_parameter(&value, config, __config_char, n);                                                                           \
            if(!__success) {                                                                                                                                   \
                REPORT_ERROR_ON_FUNCTION("Error parsing vector parameter!\n");                                                                                 \
            }                                                                                                                                                  \
        }                                                                                                                                                      \
    } while(0)

#endif // MONOALG3D_CONFIG_HELPERS_H
