//
// Created by sachetto on 22/10/17.
//

#include "config_helpers.h"
#include "../3dparty/sds/sds.h"
#include "../3dparty/stb_ds.h"
#include "../logger/logger.h"
#include <stdio.h>
#include <stdlib.h>

void report_error_on_function(int line, const char * file, const char *error) {
    log_error_and_exit("Error on line %d of file %s: %s.\n", line, file, error);
}

void report_parameter_error_on_function(int line, const char * file, const char *parameter) {
    log_error_and_exit("Error: Call on line %d of file %s needs a parameter named %s. Exiting!\n", line, file, parameter);
}

char *get_string_parameter (struct string_hash_entry *config, const char *parameter) {
    char *config_value = shget (config, parameter);
    return config_value;
}

bool get_vector_parameter(real_cpu **v, struct string_hash_entry *config, const char *parameter, int n) {

    *v = malloc(sizeof(real_cpu) * n);

    int c;

    sds config_value = sdsnew(parameter);

    config_value = sdstrim(config_value, "[ ]");
    sds *components = sdssplit(config_value, ",", &c);

    for(int i =0; i < n; i++) {
        int expr_parse_error = 0;
        real_cpu expr_parse_result = (real_cpu) te_interp(components[i], &expr_parse_error);
        if(expr_parse_error) {
            return false;
        }
        else {
            (*v)[i] = expr_parse_result;
        }
    }
    sdsfreesplitres(components, c);

    return true;
}

bool get_vector3_parameter(real_cpu v[3], struct string_hash_entry *config, const char *parameter) {

    int c;

    sds config_value = sdsnew(parameter);

    config_value = sdstrim(config_value, "[ ]");
    sds *components = sdssplit(config_value, ",", &c);

    if (c != 3) {
        sdsfreesplitres(components, c);
        return false;
    }

    for(int i =0; i < 3; i++) {
        int expr_parse_error = 0;
        real_cpu expr_parse_result = (real_cpu) te_interp(components[i], &expr_parse_error);
        if(expr_parse_error) {
            return false;
        }
        else {
            v[i] = expr_parse_result;
        }
    }
    sdsfreesplitres(components, c);

    return true;
}
