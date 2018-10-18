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

    sds tmp = sdsnew (input_dir);
    tmp = sdscat (tmp, "/simulation_state.dat");

    FILE *input_file = fopen (tmp, "r");

    sdsfree (tmp);

    if(!input_file) {
        fprintf(stderr, "Error opening %s file for restoring the simulation state\n", tmp);
        exit(10);
    }


    //Here we restore the saved domain
    if(the_grid) {
        int number_of_cells = 0;
        int num_active_cells = 0;
        double side_length;
        struct point_voidp_hash *mesh_hash = point_voidp_hash_create();

        fscanf(input_file, "%lf\n", &side_length);
        fscanf(input_file, "%d\n", &number_of_cells);
        fscanf(input_file, "%d\n", &num_active_cells);

        struct point_3d mp;
        size_t mesh_data_size = 10*sizeof(double); //TODO: we will probably want to change this size
        //Read the mesh to a point hash
        for(int i = 0; i < number_of_cells; i++ ) {

            //This is how we save the grid state:
            //fprintf (output_file, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
            //         grid_cell->center_x, grid_cell->center_y, grid_cell->center_z,
            //         grid_cell->v,
            //         grid_cell->north_flux, grid_cell->south_flux, grid_cell->east_flux,
            //         grid_cell->west_flux, grid_cell->front_flux,  grid_cell->back_flux,
            //         grid_cell->b,
            //         grid_cell->can_change, grid_cell->active);
            //We have to restore it the same way...

            double *mesh_data = (double*) malloc(mesh_data_size);
            fscanf (input_file, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
                    &(mp.x), &(mp.y), &(mp.z),
                    &mesh_data[0], &mesh_data[1], &mesh_data[2],
                    &mesh_data[3], &mesh_data[4], &mesh_data[5],
                    &mesh_data[6], &mesh_data[7], &mesh_data[8],
                    &mesh_data[9]);

            point_voidp_hash_insert(mesh_hash, mp, mesh_data);
        }

        initialize_and_construct_grid (the_grid, side_length);
        struct cell_node *grid_cell = the_grid->first_cell;

        printf("Restoring grid state...\n");

        double *cell_data;

        while(grid_cell) {

            if(grid_cell->visited) {
                //we have already seen this cell before
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

                    grid_cell->v          = cell_data[0];
                    grid_cell->north_flux = cell_data[1];
                    grid_cell->south_flux = cell_data[2];
                    grid_cell->east_flux  = cell_data[3];
                    grid_cell->west_flux  = cell_data[4];
                    grid_cell->front_flux = cell_data[5];
                    grid_cell->back_flux  = cell_data[6];
                    grid_cell->b          = cell_data[7];


                    if(cell_data[8] == 0.0) {
                        grid_cell->can_change = false;
                    }

                    if(cell_data[9] == 0.0) {
                        grid_cell->active = false;
                    }

                    grid_cell->visited = true;
                    grid_cell = grid_cell->next;
                }
                else {
                    //This grid_cell is not in the mesh. We need to refine it and check again and again....
                    refine_grid_cell(the_grid, grid_cell);
                    grid_cell = the_grid->first_cell;
                }
            }
        }

        assert(number_of_cells == the_grid->number_of_cells);

        fclose(input_file);
        point_voidp_hash_destroy(mesh_hash);
    }

}