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
#include "../hash/point_voidp_hash.h"

RESTORE_STATE(restore_simulation_state) {



    //Here we restore the saved domain
    if(the_grid) {

        sds tmp = sdsnew (input_dir);
        tmp = sdscat (tmp, "/grid_checkpoint.dat");

        FILE *input_file = fopen (tmp, "rb");

        sdsfree (tmp);

        if(!input_file) {
            fprintf(stderr, "Error opening %s file for restoring the simulation state\n", tmp);
            exit(10);
        }

        int number_of_cells = 0;
        int num_active_cells = 0;
        float side_length;
        struct point_voidp_hash *mesh_hash = point_voidp_hash_create();

        fread(&side_length, sizeof(the_grid->side_length), 1, input_file);
        fread(&number_of_cells, sizeof(the_grid->number_of_cells), 1, input_file);
        fread(&num_active_cells, sizeof(the_grid->num_active_cells), 1, input_file);

        struct point_3d mp;

        initialize_and_construct_grid (the_grid, side_length);
        struct cell_node *grid_cell = the_grid->first_cell;

        int num_data = 11;

        //Read the mesh to a point hash
        for(int i = 0; i < number_of_cells; i++ ) {

            double *mesh_data = (double*) calloc(num_data, sizeof(double));

            //Read center_x, center_y, center_z
            fread(&mp, sizeof(mp), 1, input_file);

            //Read v, north_flux, south_flux, east_flux, west_flux, front_flux,
            //back_flux, b
            fread(&mesh_data[0], sizeof(double), 8, input_file);

            //read can_change
            fread(&mesh_data[8], sizeof(grid_cell->can_change), 1, input_file);

            //read active
            fread(&mesh_data[9], sizeof(grid_cell->active), 1, input_file);

            point_voidp_hash_insert(mesh_hash, mp, mesh_data);
        }

        printf("Restoring grid state...\n");

        double *cell_data;

        while(grid_cell) {

            if(grid_cell->visited) {
                //This cell is already restored
                grid_cell = grid_cell->next;
            }
            else {

                #ifndef NDEBUG
                if( (the_grid->number_of_cells % 1000 == 0) || (the_grid->number_of_cells == number_of_cells) ) {
                    printf("Target mesh size: %d, current mesh size: %d\n", number_of_cells, the_grid->number_of_cells);
                }
                #endif

                struct point_3d mesh_point;
                mesh_point.x = grid_cell->center_x;
                mesh_point.y = grid_cell->center_y;
                mesh_point.z = grid_cell->center_z;

                if((cell_data = point_voidp_hash_search(mesh_hash, mesh_point)) != (void*)-1 ) {
                    //This grid_cell is already in the mesh. We only need to restore the data associated to it...

                    //If the cell is not active we don't need to recover its state
                    if(cell_data[9] == 0.0) {
                        grid_cell->active = false;
                    }
                    else {
                        grid_cell->v = cell_data[0];
                        grid_cell->north_flux = cell_data[1];
                        grid_cell->south_flux = cell_data[2];
                        grid_cell->east_flux = cell_data[3];
                        grid_cell->west_flux = cell_data[4];
                        grid_cell->front_flux = cell_data[5];
                        grid_cell->back_flux = cell_data[6];
                        grid_cell->b = cell_data[7];

                        if (cell_data[8] == 0.0) {
                            grid_cell->can_change = false;
                        }
                    }

                    grid_cell->visited = true;
                    grid_cell = grid_cell->next;
                }
                else {
                    //This grid_cell is not in the  mesh. We need to refine it and check again and again....
                    refine_grid_cell(the_grid, grid_cell);
                    grid_cell = the_grid->first_cell;
                }
            }
        }

        assert(number_of_cells == the_grid->number_of_cells);

        fclose(input_file);
        point_voidp_hash_destroy(mesh_hash);
    }

    if(the_monodomain_solver) {
        sds tmp = sdsnew (input_dir);
        tmp = sdscat (tmp, "/monodomain_solver_checkpoint.dat");

        FILE *input_file = fopen (tmp, "rb");

        sdsfree (tmp);

        if(!input_file) {
            fprintf(stderr, "Error opening %s file for restoring the simulation state\n", tmp);
            exit(10);
        }

        fread(the_monodomain_solver, sizeof(struct monodomain_solver), 1, input_file);
        fclose(input_file);
    }


    if(the_ode_solver) {

        sds tmp = sdsnew (input_dir);
        tmp = sdscat (tmp, "/ode_solver_checkpoint.dat");

        FILE *input_file = fopen (tmp, "rb");

        sdsfree (tmp);

        //TODO: here we have to restore the ode_solver state. We need to be careful about the SV GPU data
        //

        fclose(input_file);
    }

}