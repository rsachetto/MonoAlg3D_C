//
// Created by sachetto on 14/10/17.
//

#include <stddef.h>
#include "config_common.h"

void init_config_common_data(struct config_common *cm) {
    cm->handle = NULL;
    cm->function_name = NULL;
    cm->function_name_was_set = false;
    cm->library_file_path = NULL;
    cm->library_file_path_was_set = false;
    cm->config = string_hash_create();
}
