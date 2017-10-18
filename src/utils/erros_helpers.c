//
// Created by sachetto on 18/10/17.
//

#include <stdio.h>
#include <stdlib.h>
#include "erros_helpers.h"

void report_parameter_error_on_function(const char * function, const char *parameter) {
    fprintf(stderr, "Error: User provided function %s needs a parameter named %s. Exiting!\n", function, parameter);
    exit(EXIT_FAILURE);
}
