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

    //Here we save the domain state
    if(the_grid){
        sds tmp = sdsnew(output_dir);
        tmp = sdscat(tmp, "/grid_checkpoint.dat");

        FILE *output_file = fopen(tmp, "wb");

        sdsfree(tmp);

        if (!output_file) {
            fprintf(stderr, "Error opening %s file for saving the simulation state\n", tmp);
            return;
        }

        fwrite(&(the_grid->side_length), sizeof(the_grid->side_length), 1, output_file);
        fwrite(&(the_grid->number_of_cells), sizeof(the_grid->number_of_cells), 1, output_file);
        fwrite(&(the_grid->num_active_cells), sizeof(the_grid->num_active_cells), 1, output_file);

        struct cell_node *grid_cell = the_grid->first_cell;

        while (grid_cell != 0) {

            fwrite(&(grid_cell->center_x), sizeof(grid_cell->center_x), 1, output_file);
            fwrite(&(grid_cell->center_y), sizeof(grid_cell->center_y), 1, output_file);
            fwrite(&(grid_cell->center_z), sizeof(grid_cell->center_z), 1, output_file);
            fwrite(&(grid_cell->v), sizeof(grid_cell->v), 1, output_file);
            fwrite(&(grid_cell->north_flux), sizeof(grid_cell->north_flux), 1, output_file);
            fwrite(&(grid_cell->south_flux), sizeof(grid_cell->south_flux), 1, output_file);
            fwrite(&(grid_cell->east_flux), sizeof(grid_cell->east_flux), 1, output_file);
            fwrite(&(grid_cell->west_flux), sizeof(grid_cell->west_flux), 1, output_file);
            fwrite(&(grid_cell->front_flux), sizeof(grid_cell->front_flux), 1, output_file);
            fwrite(&(grid_cell->back_flux), sizeof(grid_cell->back_flux), 1, output_file);
            fwrite(&(grid_cell->b), sizeof(grid_cell->b), 1, output_file);
            fwrite(&(grid_cell->can_change), sizeof(grid_cell->can_change), 1, output_file);
            fwrite(&(grid_cell->active), sizeof(grid_cell->active), 1, output_file);

            grid_cell = grid_cell->next;
        }

        fclose(output_file);

    }

    //Here we save the monodomain solver state
    if(the_monodomain_solver) {
        sds tmp = sdsnew(output_dir);
        tmp = sdscat(tmp, "/monodomain_solver_checkpoint.dat");

        FILE *output_file = fopen(tmp, "wb");

        sdsfree(tmp);

        if (!output_file) {
            fprintf(stderr, "Error opening %s file for saving the simulation state\n", tmp);
            return;
        }


        fwrite(the_monodomain_solver, sizeof(struct monodomain_solver), 1, output_file);
        fclose(output_file);

    }

    if(the_ode_solver) {

        sds tmp = sdsnew(output_dir);
        tmp = sdscat(tmp, "/ode_solver_checkpoint.dat");

        FILE *output_file = fopen(tmp, "wb");

        sdsfree(tmp);

        if (!output_file) {
            fprintf(stderr, "Error opening %s file for saving the simulation state\n", tmp);
            return;
        }

        //TODO: here we have to save the ode_solver state....


        fclose(output_file);

    }

}