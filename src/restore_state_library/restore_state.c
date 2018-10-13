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

    //We are going to start recovering the domain mesh first. After that,
    //we need to restore the ode solver and the monodomain solver states

    int number_of_cells = 0;
    struct point_voidp_hash *mesh_hash = point_voidp_hash_create();

    //Read the mesh to a point hash
    fscanf(input_file, "%d\n", &number_of_cells);
    struct point_3d mesh_point;
    size_t mesh_data_size = 3*sizeof(double); //TODO: we will probably want to change this size
    double *mesh_data = (double*) malloc(mesh_data);

    for(int i = 0; i < number_of_cells; i++ ) {
        fscanf (input_file, "%lf,%lf,%lf,%lf,%lf\n", &(mesh_point.x), &(mesh_point.y),
                &(mesh_point.z), &mesh_data[0], &mesh_data[1]);

        point_voidp_hash_insert(mesh_hash, mesh_point, mesh_data);
    }

    //TODO: how do we restore the exact same mesh??????


    fclose(input_file);


}