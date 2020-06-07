//
// Created by sachetto on 13/10/17.
//

#include "save_state_config.h"

#include <dlfcn.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "../utils/file_utils.h"
#include "../3dparty/stb_ds.h"
#include "../config_helpers/config_helpers.h"

void print_save_state_config_values(struct config* s) {

    int save_rate = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(int, save_rate, s->config_data, "save_rate");

    printf("Save state save_rate: %d\n",save_rate);
    printf("Save state library_file: %s\n", s->library_file_path);
    printf("Save state main function: %s\n",s->main_function_name);
    printf("Save state init function: %s\n",s->init_function_name);
    printf("Save state end function: %s\n", s->end_function_name);

    printf("Save state config:\n");
    STRING_HASH_PRINT_KEY_VALUE(s->config_data);
}


