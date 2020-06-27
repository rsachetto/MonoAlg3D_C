//
// Created by sachetto on 13/10/17.
//

#include <stdlib.h>
#include <string.h>

#include "../alg/grid/grid.h"
#include "../config/save_state_config.h"
#include "../3dparty/sds/sds.h"


#ifdef COMPILE_CUDA
#include "../models_library/model_gpu_utils.h"
#endif

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

        fwrite(&(the_grid->cube_side_length.x), sizeof(the_grid->cube_side_length.x), 1, output_file);
        fwrite(&(the_grid->cube_side_length.y), sizeof(the_grid->cube_side_length.y), 1, output_file);
        fwrite(&(the_grid->cube_side_length.z), sizeof(the_grid->cube_side_length.z), 1, output_file);
        fwrite(&(the_grid->mesh_side_length.x), sizeof(the_grid->mesh_side_length.x), 1, output_file);
        fwrite(&(the_grid->mesh_side_length.y), sizeof(the_grid->mesh_side_length.y), 1, output_file);
        fwrite(&(the_grid->mesh_side_length.z), sizeof(the_grid->mesh_side_length.z), 1, output_file);
        fwrite(&(the_grid->number_of_cells), sizeof(the_grid->number_of_cells), 1, output_file);
        fwrite(&(the_grid->num_active_cells), sizeof(the_grid->num_active_cells), 1, output_file);

        struct cell_node *grid_cell = the_grid->first_cell;

        while (grid_cell != 0) {

            fwrite(&(grid_cell->center.x), sizeof(grid_cell->center.x), 1, output_file);
            fwrite(&(grid_cell->center.y), sizeof(grid_cell->center.y), 1, output_file);
            fwrite(&(grid_cell->center.z), sizeof(grid_cell->center.z), 1, output_file);
            fwrite(&(grid_cell->v), sizeof(grid_cell->v), 1, output_file);
            fwrite(&(grid_cell->z_front_flux), sizeof(grid_cell->z_front_flux), 1, output_file);
            fwrite(&(grid_cell->z_back_flux), sizeof(grid_cell->z_back_flux), 1, output_file);
            fwrite(&(grid_cell->y_top_flux), sizeof(grid_cell->y_top_flux), 1, output_file);
            fwrite(&(grid_cell->y_down_flux), sizeof(grid_cell->y_down_flux), 1, output_file);
            fwrite(&(grid_cell->x_right_flux), sizeof(grid_cell->x_right_flux), 1, output_file);
            fwrite(&(grid_cell->x_left_flux), sizeof(grid_cell->x_left_flux), 1, output_file);
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
        fwrite(time_info, sizeof(struct time_info), 1, output_file);

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

        fwrite(&(the_ode_solver->max_dt), sizeof(the_ode_solver->max_dt), 1, output_file);
        fwrite(&(the_ode_solver->min_dt), sizeof(the_ode_solver->min_dt), 1, output_file);
        fwrite(&(the_ode_solver->rel_tol), sizeof(the_ode_solver->rel_tol), 1, output_file);
        fwrite(&(the_ode_solver->abs_tol), sizeof(the_ode_solver->abs_tol), 1, output_file);

//        fwrite(&(the_ode_solver->previous_dt), sizeof(the_ode_solver->previous_dt), 1, output_file);
//        fwrite(&(the_ode_solver->time_new), sizeof(the_ode_solver->time_new), 1, output_file);


        size_t num_cells_to_solve = 0;

        if(the_ode_solver->cells_to_solve) {
            num_cells_to_solve = the_ode_solver->num_cells_to_solve;
            fwrite(&num_cells_to_solve, sizeof(the_ode_solver->num_cells_to_solve), 1, output_file);
            fwrite(the_ode_solver->cells_to_solve, sizeof(the_ode_solver->cells_to_solve[0]), num_cells_to_solve,
                   output_file);
        }

        else {
            fwrite(&num_cells_to_solve, sizeof(the_ode_solver->num_cells_to_solve), 1, output_file);
        }

        fwrite(&(the_ode_solver->gpu), sizeof(the_ode_solver->gpu), 1, output_file);
        fwrite(&(the_ode_solver->gpu_id), sizeof(the_ode_solver->gpu_id), 1, output_file);

//        fwrite(&(the_ode_solver->model_data), sizeof(the_ode_solver->model_data), 1, output_file);
//        unsigned long data_size = strlen(the_ode_solver->model_data.model_library_path);
//        fwrite(&(data_size), sizeof(data_size), 1, output_file);
//        fwrite(the_ode_solver->model_data.model_library_path, data_size, 1, output_file);

        fwrite(&(the_ode_solver->pitch), sizeof(the_ode_solver->pitch), 1, output_file);
        fwrite(&(the_ode_solver->original_num_cells), sizeof(the_ode_solver->original_num_cells), 1, output_file);
        
        if(the_ode_solver->gpu) {

        #ifdef COMPILE_CUDA
            real *sv_cpu;
            sv_cpu = (real*) malloc(the_ode_solver->original_num_cells * the_ode_solver->model_data.number_of_ode_equations * sizeof(real));

            check_cuda_error(cudaMemcpy2D(sv_cpu, the_ode_solver->original_num_cells * sizeof(real), the_ode_solver->sv, the_ode_solver->pitch,
                         the_ode_solver->original_num_cells * sizeof(real), (size_t)the_ode_solver->model_data.number_of_ode_equations,
                         cudaMemcpyDeviceToHost));

            fwrite(sv_cpu, sizeof(real), the_ode_solver->original_num_cells * the_ode_solver->model_data.number_of_ode_equations,
                    output_file);
        #endif

        }
        else {
            fwrite(the_ode_solver->sv, sizeof(real), the_ode_solver->original_num_cells * the_ode_solver->model_data.number_of_ode_equations,
                   output_file);

        }

        fwrite(&(the_ode_solver->extra_data_size), sizeof(the_ode_solver->extra_data_size), 1, output_file);
        fwrite(the_ode_solver->ode_extra_data, the_ode_solver->extra_data_size, 1, output_file);

        fclose(output_file);

    }

}
