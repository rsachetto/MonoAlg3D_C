//
// Created by sachetto on 13/10/17.
//

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#include "../alg/grid/grid.h"
#include "../config/restore_state_config.h"

#include "../common_types/common_types.h"
#include "../single_file_libraries/stb_ds.h"

#include "../libraries_common/config_helpers.h"

#ifdef COMPILE_CUDA
#include "../models_library/model_gpu_utils.h"
#endif

#include "../string/sds.h"

RESTORE_STATE (restore_simulation_state) {

    // Here we restore the saved domain
    if (the_grid) {

        sds tmp = sdsnew (input_dir);
        tmp = sdscat (tmp, "/grid_checkpoint.dat");

        FILE *input_file = fopen (tmp, "rb");

        sdsfree (tmp);

        if (!input_file) {
            fprintf (stderr, "Error opening %s file for restoring the simulation state\n", tmp);
            exit (10);
        }

        int number_of_cells = 0;
        int num_active_cells = 0;
        real_cpu side_length_x, side_length_y, side_length_z;
        struct point_voidp_hash_entry *mesh_hash = NULL;
        hmdefault(mesh_hash, NULL);

        fread (&side_length_x, sizeof (the_grid->side_length_x), 1, input_file);
        fread (&side_length_y, sizeof (the_grid->side_length_y), 1, input_file);
        fread (&side_length_z, sizeof (the_grid->side_length_z), 1, input_file);
        fread (&number_of_cells, sizeof (the_grid->number_of_cells), 1, input_file);
        fread (&num_active_cells, sizeof (the_grid->num_active_cells), 1, input_file);

        struct point_3d mp;

        initialize_and_construct_grid (the_grid, side_length_x, side_length_y, side_length_z);
        struct cell_node *grid_cell = the_grid->first_cell;

        int num_data = 11;

        // Read the mesh to a point hash
        for (int i = 0; i < number_of_cells; i++) {

            real_cpu *mesh_data = (real_cpu *)calloc (num_data, sizeof (real_cpu));

            // Read center_x, center_y, center_z
            fread (&mp, sizeof (mp), 1, input_file);

            // Read v, north_flux, south_flux, east_flux, west_flux, front_flux,
            // back_flux, b
            fread (&mesh_data[0], sizeof (real_cpu), 8, input_file);

            // read can_change
            fread (&mesh_data[8], sizeof (grid_cell->can_change), 1, input_file);

            // read active
            fread (&mesh_data[9], sizeof (grid_cell->active), 1, input_file);

            hmput(mesh_hash, mp, mesh_data);
        }

        printf ("Restoring grid state...\n");

        real_cpu *cell_data;

        while (grid_cell) {

            if (grid_cell->visited) {
                // This cell is already restored
                grid_cell = grid_cell->next;
            } else {

                struct point_3d mesh_point;
                mesh_point.x = grid_cell->center_x;
                mesh_point.y = grid_cell->center_y;
                mesh_point.z = grid_cell->center_z;

                if ((cell_data = hmget (mesh_hash, mesh_point)) != NULL) {

                    // This grid_cell is already in the mesh. We only need to restore the data associated to it...
                    // If the cell is not active we don't need to recover its state
                    if (cell_data[9] == 0.0) {
                        grid_cell->active = false;
                    } else {
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

                } else {
                    // This grid_cell is not in the  mesh. We need to refine it and check again and again....
                    refine_grid_cell (the_grid, grid_cell);
                    grid_cell = the_grid->first_cell;
                }
            }
        }

        assert (number_of_cells == the_grid->number_of_cells);

        fclose (input_file);
        hmfree(mesh_hash);
    }

    if (the_monodomain_solver) {

        printf ("Restoring monodomain solver state...\n");

        sds tmp = sdsnew (input_dir);
        tmp = sdscat (tmp, "/monodomain_solver_checkpoint.dat");

        FILE *input_file = fopen (tmp, "rb");

        if (!input_file) {
            fprintf (stderr, "Error opening %s file for restoring the simulation state\n", tmp);
            exit (10);
        }

        sdsfree (tmp);

        fread (the_monodomain_solver, sizeof (struct monodomain_solver), 1, input_file);
        fclose (input_file);
    }

    if (the_ode_solver) {

        printf ("Restoring ode solver state...\n");

        sds tmp = sdsnew (input_dir);
        tmp = sdscat (tmp, "/ode_solver_checkpoint.dat");

        FILE *input_file = fopen (tmp, "rb");

        sdsfree (tmp);

        fread (&(the_ode_solver->max_dt), sizeof (the_ode_solver->max_dt), 1, input_file);
        fread (&(the_ode_solver->min_dt), sizeof (the_ode_solver->min_dt), 1, input_file);
        fread (&(the_ode_solver->rel_tol), sizeof (the_ode_solver->rel_tol), 1, input_file);
        fread (&(the_ode_solver->abs_tol), sizeof (the_ode_solver->abs_tol), 1, input_file);

        fread (&(the_ode_solver->previous_dt), sizeof (the_ode_solver->previous_dt), 1, input_file);
        fread (&(the_ode_solver->time_new), sizeof (the_ode_solver->time_new), 1, input_file);

        fread (&(the_ode_solver->num_cells_to_solve), sizeof (the_ode_solver->num_cells_to_solve), 1, input_file);

        size_t num_cells_to_solve = the_ode_solver->num_cells_to_solve;

        the_ode_solver->cells_to_solve = NULL;

        if(num_cells_to_solve) {
            free(the_ode_solver->cells_to_solve);
            the_ode_solver->cells_to_solve =
                    malloc(sizeof(the_ode_solver->cells_to_solve[0]) * the_ode_solver->num_cells_to_solve);
            fread(the_ode_solver->cells_to_solve, sizeof(the_ode_solver->cells_to_solve[0]), num_cells_to_solve,
                  input_file);
        }

        fread (&(the_ode_solver->gpu), sizeof (the_ode_solver->gpu), 1, input_file);
        fread (&(the_ode_solver->gpu_id), sizeof (the_ode_solver->gpu_id), 1, input_file);

        fread (&(the_ode_solver->model_data), sizeof (the_ode_solver->model_data), 1, input_file);
        fread (&(the_ode_solver->pitch), sizeof (the_ode_solver->pitch), 1, input_file);

        fread (&(the_ode_solver->original_num_cells), sizeof (the_ode_solver->original_num_cells), 1, input_file);
        if (the_ode_solver->gpu) {
#ifdef COMPILE_CUDA
            real *sv_cpu;
            sv_cpu = (real *)malloc (the_ode_solver->original_num_cells * the_ode_solver->model_data.number_of_ode_equations *
                                     sizeof (real));
            fread (sv_cpu, sizeof (real),
                   the_ode_solver->original_num_cells * the_ode_solver->model_data.number_of_ode_equations, input_file);

            check_cuda_error(cudaMemcpy2D (the_ode_solver->sv, the_ode_solver->pitch, sv_cpu, the_ode_solver->original_num_cells * sizeof (real),
                          the_ode_solver->original_num_cells * sizeof (real),
                          (size_t)the_ode_solver->model_data.number_of_ode_equations, cudaMemcpyHostToDevice));
#endif
        } else {
            fread (the_ode_solver->sv, sizeof (real),
                   the_ode_solver->original_num_cells * the_ode_solver->model_data.number_of_ode_equations, input_file);
        }

        fread(&(the_ode_solver->extra_data_size), sizeof(the_ode_solver->extra_data_size), 1, input_file);
        fread (the_ode_solver->ode_extra_data, the_ode_solver->extra_data_size, 1, input_file);

        fclose (input_file);
    }
}
