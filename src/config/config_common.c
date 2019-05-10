//
// Created by sachetto on 14/10/17.
//

#include <stddef.h>
#include "config_common.h"
#include <dlfcn.h>


#include "../single_file_libraries/stb_ds.h"

void init_config_common_data(struct config_common *cm) {
    cm->handle = NULL;
    cm->function_name = NULL;
    cm->function_name_was_set = false;
    cm->library_file_path = NULL;
    cm->library_file_path_was_set = false;
    cm->config = NULL;
    sh_new_arena(cm->config);
    shdefault(cm->config, NULL);
}

void free_config_common_data(struct config_common *cm) {

    ptrdiff_t n = shlen(cm->config);

    for(int i = 0; i < n; i++) {
        free(cm->config[i].value);
    }

    free(cm->library_file_path);
    free(cm->function_name);
    if(cm->handle)
        dlclose(cm->handle);

    shfree(cm->config);
}
