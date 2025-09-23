//
// Created by sachetto on 22/10/17.
//

#include "config_helpers.h"
#include "../3dparty/stb_ds.h"
#include "../logger/logger.h"
#include <stdlib.h>

void report_error_on_function(int line, const char *file, const char *error) {
    log_error_and_exit("Error on line %d of file %s: %s.\n", line, file, error);
}

void report_error_on_function_and_continue(int line, const char *file, const char *error) {
    log_error("Error on line %d of file %s: %s.\n", line, file, error);
}

void report_parameter_error_on_function(int line, const char *file, const char *parameter) {
    log_error_and_exit("Error: Call on line %d of file %s needs a parameter named %s. Exiting!\n", line, file, parameter);
}

char *get_string_parameter(struct string_hash_entry *config, const char *parameter) {
    char *config_value = shget(config, parameter);
    return config_value;
}

bool get_vector_parameter(real_cpu **v, const char *parameter, int n) {

    *v = (real_cpu *)malloc(sizeof(real_cpu) * n);

    int c;

    sds config_value = sdsnew(parameter);

    config_value = sdstrim(config_value, "[ ]");
    sds *components = sdssplit(config_value, ",", &c);

    for(int i = 0; i < n; i++) {
        int expr_parse_error = 0;
        real_cpu expr_parse_result = (real_cpu)te_interp(components[i], &expr_parse_error);
        if(expr_parse_error) {
            return false;
        } else {
            (*v)[i] = expr_parse_result;
        }
    }
    sdsfreesplitres(components, c);

    return true;
}

bool get_vector3_parameter(real_cpu v[3], const char *parameter) {

    int c;

    sds config_value = sdsnew(parameter);

    config_value = sdstrim(config_value, "[ ]");
    sds *components = sdssplit(config_value, ",", &c);

    if(c != 3) {
        sdsfreesplitres(components, c);
        return false;
    }

    for(int i = 0; i < 3; i++) {
        int expr_parse_error = 0;
        real_cpu expr_parse_result = (real_cpu)te_interp(components[i], &expr_parse_error);
        if(expr_parse_error) {
            return false;
        } else {
            v[i] = expr_parse_result;
        }
    }
    sdsfreesplitres(components, c);

    return true;
}

bool get_matrix_parameter(real_cpu **v, const char *parameter, int nlin, int ncol) {

    *v = (real_cpu *)malloc(sizeof(real_cpu) * nlin * ncol);

    int c;

    sds config_value = sdsnew(parameter);

    config_value = sdstrim(config_value, "{ }");
    sds *components = sdssplit(config_value, "_", &c);

    for(int i = 0; i < nlin; i++) {

        sds config_value_2 = sdsnew(components[i]);

        config_value_2 = sdstrim(config_value_2, "[ ]");
        sds *components_2 = sdssplit(config_value_2, ",", &c);

        for(int j = 0; j < ncol; j++) {
            int expr_parse_error = 0;
            real_cpu expr_parse_result = (real_cpu)te_interp(components_2[j], &expr_parse_error);
            if(expr_parse_error) {
                return false;
            } else {
                (*v)[i * ncol + j] = expr_parse_result;
            }
        }

        sdsfreesplitres(components_2, c);
    }
    sdsfreesplitres(components, c);

    return true;
}
