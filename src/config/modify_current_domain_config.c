//
// Created by sachetto on 13/10/17.
//

#include "stim_config.h"
#include <dlfcn.h>
#include "../utils/file_utils.h"

#include "../3dparty/stb_ds.h"
#include "../config_helpers/config_helpers.h"

void print_modify_domain_config_values(struct config* s) {

    if(s == NULL) {
        log_to_stdout_and_file("No Modify domain configuration.\n");
        return;
    }

    //TODO: implement...


}