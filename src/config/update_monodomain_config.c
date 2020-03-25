//
// Created by sachetto on 13/10/17.
//

#include "update_monodomain_config.h"

#include <dlfcn.h>
#include <string.h>
#include "../utils/file_utils.h"

#include "../common_types/common_types.h"
#include "../3dparty/stb_ds.h"

void print_update_monodomain_config_values(struct config* s) {

    if (s == NULL) {
        log_to_stdout_and_file("[update_monodomain] No Update Monodomain configuration.\n");
        return;
    }

    log_to_stdout_and_file("[update_monodomain] Update Monodomain function configuration:\n");

    log_to_stdout_and_file("[update_monodomain] Update Monodomain library: %s\n", s->library_file_path);
    log_to_stdout_and_file("[update_monodomain] Update Monodomain main function: %s\n", s->main_function_name);

    if(s->init_function_name)
        log_to_stdout_and_file("[update_monodomain] Update Monodomain init function: %s\n", s->init_function_name);

    if(s->end_function_name)
        log_to_stdout_and_file("[update_monodomain] Update Monodomain end function: %s\n", s->end_function_name);

    if (shlen(s->config_data) == 1) {
        log_to_stdout_and_file("[update_monodomain] Update Monodomain parameter:\n");
    } else if (shlen(s->config_data) > 1) {
        log_to_stdout_and_file("[update_monodomain] Update Monodomain parameters:\n");
    }

    STRING_HASH_PRINT_KEY_VALUE_LOG(s->config_data);

}