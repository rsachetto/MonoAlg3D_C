//
// Created by sachetto on 13/10/17.
//

#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>

#include "../alg/grid/grid.h"
#include "../config/restore_state_config.h"
#include "../libraries_common/config_helpers.h"
#include "../string/sds.h"

RESTORE_STATE(restore_simulation_state) {   

    sds tmp = sdsnew (output_dir);    
    tmp = sdscat (tmp, "/simulation_state.dat");    

    FILE *input_file = fopen (tmp, "r");    

    sdsfree (tmp);

    if(!input_file) {
        fprintf(stderr, "Error opening %s file for restoring the simulation state\n", tmp);
        exit(10);
    }   

    //TODO: we need to restore the mash and ode solver states


    fclose(input_file);


}