//
// Created by sachetto on 13/10/17.
//

#include "restore_state_config.h"

#include <dlfcn.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "../utils/file_utils.h"

#include "../common_types/common_types.h"
#include "../3dparty/stb_ds.h"

void print_restore_state_config_values(struct config* s) {
    printf("Restore state library_file: %s\n", s->library_file_path);
    printf("Restore state main function: %s\n", s->main_function_name);

    printf("Restore state init function: %s\n", s->init_function_name);
    printf("Restore state end function: %s\n", s->end_function_name);

    printf("restore state_config:\n");
    STRING_HASH_PRINT_KEY_VALUE(s->config_data);
}