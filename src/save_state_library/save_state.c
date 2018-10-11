//
// Created by sachetto on 13/10/17.
//

#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>

#include "../alg/grid/grid.h"
#include "../config/save_state_config.h"
#include "../libraries_common/config_helpers.h"
#include "../string/sds.h"

SAVE_STATE(save_simulation_state) {   

    sds tmp = sdsnew (output_dir);    
    tmp = sdscat (tmp, "/simulation_state.dat");    

    FILE *output_file = fopen (tmp, "w");    

    sdsfree (tmp);

    if(!output_file) {
        fprintf(stderr, "Error opening %s file for saving the simulation state\n", tmp);
        return;
    }

    //FIXME: maybe we should write as binary here. Using text for debugging

    fprintf(output_file, "%d\n", the_grid->number_of_cells);
    fprintf(output_file, "%d\n", the_grid->num_active_cells);

    struct cell_node *grid_cell = the_grid->first_cell;    

    while (grid_cell != 0) {
        if (grid_cell->active) {
            fprintf (output_file, "%lf,%lf,%lf,%lf,%lf\n", grid_cell->center_x, grid_cell->center_y,
                     grid_cell->center_z, grid_cell->half_face_length, grid_cell->v);
        }
        grid_cell = grid_cell->next;
    }   

    fclose(output_file);
    

    //TODO: we need to restore the ode solver state as well

}