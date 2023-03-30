//
// Created by sachetto on 13/10/17.
//

#include <stdlib.h>
#include <string.h>

#include "../alg/grid/grid.h"
#include "../config/save_state_config.h"
#include "../save_mesh_library/save_mesh_helper.h"

#ifdef COMPILE_CUDA
#include "../gpu_utils/gpu_utils.h"
#endif

SAVE_STATE(save_simulation_state) {
    //Here we save the domain state
    if(the_grid) {
        sds tmp = sdsnew(output_dir);
        tmp = sdscat(tmp, "/grid_checkpoint.dat");

        FILE *output_file = fopen(tmp, "wb");

        sdsfree(tmp);

        if (!output_file) {
            fprintf(stderr, "Error opening %s file for saving the simulation state\n", tmp);
            return;
        }

        fwrite(&(the_grid->cube_side_length), sizeof(the_grid->cube_side_length), 1, output_file);
        fwrite(&(the_grid->mesh_side_length), sizeof(the_grid->mesh_side_length), 1, output_file);
        fwrite(&(the_grid->number_of_cells),  sizeof(the_grid->number_of_cells),  1, output_file);
        fwrite(&(the_grid->num_active_cells), sizeof(the_grid->num_active_cells), 1, output_file);

        struct cell_node *grid_cell = the_grid->first_cell;

        while (grid_cell != 0) {

            if(!grid_cell->mesh_extra_info) grid_cell->mesh_extra_info_size = 0;

            fwrite(&(grid_cell->center),               sizeof(grid_cell->center),                   1, output_file);
            fwrite(&(grid_cell->v),                    sizeof(grid_cell->v),                        1, output_file);
            fwrite(&(grid_cell->front_flux),           sizeof(grid_cell->front_flux),               1, output_file);
            fwrite(&(grid_cell->back_flux),            sizeof(grid_cell->back_flux),                1, output_file);
            fwrite(&(grid_cell->top_flux),             sizeof(grid_cell->top_flux),                 1, output_file);
            fwrite(&(grid_cell->down_flux),            sizeof(grid_cell->down_flux),                1, output_file);
            fwrite(&(grid_cell->right_flux),           sizeof(grid_cell->right_flux),               1, output_file);
            fwrite(&(grid_cell->left_flux),            sizeof(grid_cell->left_flux),                1, output_file);
            fwrite(&(grid_cell->b),                    sizeof(grid_cell->b),                        1, output_file);
            fwrite(&(grid_cell->can_change),           sizeof(grid_cell->can_change),               1, output_file);
            fwrite(&(grid_cell->active),               sizeof(grid_cell->active),                   1, output_file);
            fwrite(&(grid_cell->mesh_extra_info_size), sizeof(grid_cell->mesh_extra_info_size),     1, output_file);

            if(grid_cell->mesh_extra_info_size)
                fwrite(grid_cell->mesh_extra_info,    grid_cell->mesh_extra_info_size,             1, output_file);

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

        fwrite(&(the_ode_solver->adaptive), sizeof(the_ode_solver->adaptive), 1, output_file);
        fwrite(&(the_ode_solver->max_dt), sizeof(the_ode_solver->max_dt), 1, output_file);
        fwrite(&(the_ode_solver->min_dt), sizeof(the_ode_solver->min_dt), 1, output_file);
        fwrite(&(the_ode_solver->rel_tol), sizeof(the_ode_solver->rel_tol), 1, output_file);
        fwrite(&(the_ode_solver->abs_tol), sizeof(the_ode_solver->abs_tol), 1, output_file);

        bool read_cells_to_solve = false;

        fwrite(&(the_ode_solver->num_cells_to_solve), sizeof(the_ode_solver->num_cells_to_solve), 1, output_file);

        if(the_ode_solver->cells_to_solve) {
            read_cells_to_solve = true;
            fwrite(the_ode_solver->cells_to_solve, sizeof(the_ode_solver->cells_to_solve[0]), the_ode_solver->num_cells_to_solve, output_file);
        }

        fwrite(&(read_cells_to_solve), sizeof(read_cells_to_solve), 1, output_file);

        fwrite(&(the_ode_solver->gpu), sizeof(the_ode_solver->gpu), 1, output_file);
        fwrite(&(the_ode_solver->gpu_id), sizeof(the_ode_solver->gpu_id), 1, output_file);

        fwrite(&(the_ode_solver->pitch), sizeof(the_ode_solver->pitch), 1, output_file);
        fwrite(&(the_ode_solver->original_num_cells), sizeof(the_ode_solver->original_num_cells), 1, output_file);

        size_t num_sv_entries = the_ode_solver->model_data.number_of_ode_equations;

        if(the_ode_solver->gpu) {

            #ifdef COMPILE_CUDA
            if(the_ode_solver->adaptive) {
                num_sv_entries = num_sv_entries + 3;
            }

            real *sv_cpu;
            sv_cpu = MALLOC_ARRAY_OF_TYPE(real, the_ode_solver->original_num_cells * num_sv_entries);

            check_cuda_error(cudaMemcpy2D(sv_cpu, the_ode_solver->original_num_cells * sizeof(real), the_ode_solver->sv, the_ode_solver->pitch,
                         the_ode_solver->original_num_cells * sizeof(real), num_sv_entries, cudaMemcpyDeviceToHost));

            fwrite(sv_cpu, sizeof(real), the_ode_solver->original_num_cells * num_sv_entries, output_file);
            #endif

        }
        else {
            fwrite(the_ode_solver->sv, sizeof(real), the_ode_solver->original_num_cells * num_sv_entries,  output_file);
        }

        //fwrite(&(the_ode_solver->extra_data_size), sizeof(the_ode_solver->extra_data_size), 1, output_file);
        //fwrite(the_ode_solver->ode_extra_data, the_ode_solver->extra_data_size, 1, output_file);

        fclose(output_file);

    }

    CALL_EXTRA_FUNCTIONS(save_state_fn, time_info, config, save_mesh_config, the_grid, the_monodomain_solver, the_ode_solver, the_purkinje_ode_solver, output_dir);
}

SAVE_STATE(save_purkinje_coupled_simulation_state) {
    
    // Firstly, we save the [domain] and [purkinje] state
    if (the_grid && the_grid->purkinje) {

        sds tmp = sdsnew(output_dir);
        tmp = sdscat(tmp, "/grid_checkpoint.dat");

        FILE *output_file = fopen(tmp, "wb");

        sdsfree(tmp);

        if (!output_file) {
            fprintf(stderr, "Error opening %s file for saving the simulation state\n", tmp);
            return;
        }

    // DOMAIN section
        fwrite(&(the_grid->cube_side_length), sizeof(the_grid->cube_side_length), 1, output_file);
        fwrite(&(the_grid->mesh_side_length), sizeof(the_grid->mesh_side_length), 1, output_file);
        fwrite(&(the_grid->number_of_cells),  sizeof(the_grid->number_of_cells),  1, output_file);
        fwrite(&(the_grid->num_active_cells), sizeof(the_grid->num_active_cells), 1, output_file);

        struct cell_node *grid_cell = the_grid->first_cell;

        while (grid_cell != 0) {

            if(!grid_cell->mesh_extra_info) grid_cell->mesh_extra_info_size = 0;

            fwrite(&(grid_cell->center),               sizeof(grid_cell->center),                   1, output_file);
            fwrite(&(grid_cell->v),                    sizeof(grid_cell->v),                        1, output_file);
            fwrite(&(grid_cell->front_flux),           sizeof(grid_cell->front_flux),               1, output_file);
            fwrite(&(grid_cell->back_flux),            sizeof(grid_cell->back_flux),                1, output_file);
            fwrite(&(grid_cell->top_flux),             sizeof(grid_cell->top_flux),                 1, output_file);
            fwrite(&(grid_cell->down_flux),            sizeof(grid_cell->down_flux),                1, output_file);
            fwrite(&(grid_cell->right_flux),           sizeof(grid_cell->right_flux),               1, output_file);
            fwrite(&(grid_cell->left_flux),            sizeof(grid_cell->left_flux),                1, output_file);
            fwrite(&(grid_cell->b),                    sizeof(grid_cell->b),                        1, output_file);
            fwrite(&(grid_cell->can_change),           sizeof(grid_cell->can_change),               1, output_file);
            fwrite(&(grid_cell->active),               sizeof(grid_cell->active),                   1, output_file);
            fwrite(&(grid_cell->mesh_extra_info_size), sizeof(grid_cell->mesh_extra_info_size),     1, output_file);

            if(grid_cell->mesh_extra_info_size)
                fwrite(grid_cell->mesh_extra_info,    grid_cell->mesh_extra_info_size,             1, output_file);

            grid_cell = grid_cell->next;
        }

    // PURKINJE section
        struct grid_purkinje *the_grid_purkinje = the_grid->purkinje;
        fwrite(&(the_grid_purkinje->cube_side_length),sizeof(the_grid_purkinje->cube_side_length), 1, output_file);
        fwrite(&(the_grid_purkinje->mesh_side_length),sizeof(the_grid_purkinje->mesh_side_length), 1, output_file);
        fwrite(&(the_grid_purkinje->number_of_purkinje_cells),sizeof(the_grid_purkinje->number_of_purkinje_cells), 1, output_file);
        fwrite(&(the_grid_purkinje->num_active_purkinje_cells),sizeof(the_grid_purkinje->num_active_purkinje_cells), 1, output_file);

        struct cell_node *pk_grid_cell = the_grid_purkinje->first_cell;
        while (pk_grid_cell != 0) {

            if (!pk_grid_cell->mesh_extra_info) pk_grid_cell->mesh_extra_info_size = 0;

            fwrite(&(pk_grid_cell->center),               sizeof(pk_grid_cell->center),                   1, output_file);
            fwrite(&(pk_grid_cell->v),                    sizeof(pk_grid_cell->v),                        1, output_file);
            fwrite(&(pk_grid_cell->front_flux),           sizeof(pk_grid_cell->front_flux),               1, output_file);
            fwrite(&(pk_grid_cell->back_flux),            sizeof(pk_grid_cell->back_flux),                1, output_file);
            fwrite(&(pk_grid_cell->top_flux),             sizeof(pk_grid_cell->top_flux),                 1, output_file);
            fwrite(&(pk_grid_cell->down_flux),            sizeof(pk_grid_cell->down_flux),                1, output_file);
            fwrite(&(pk_grid_cell->right_flux),           sizeof(pk_grid_cell->right_flux),               1, output_file);
            fwrite(&(pk_grid_cell->left_flux),            sizeof(pk_grid_cell->left_flux),                1, output_file);
            fwrite(&(pk_grid_cell->b),                    sizeof(pk_grid_cell->b),                        1, output_file);
            fwrite(&(pk_grid_cell->can_change),           sizeof(pk_grid_cell->can_change),               1, output_file);
            fwrite(&(pk_grid_cell->active),               sizeof(pk_grid_cell->active),                   1, output_file);
            fwrite(&(pk_grid_cell->mesh_extra_info_size), sizeof(pk_grid_cell->mesh_extra_info_size),     1, output_file);

            if(pk_grid_cell->mesh_extra_info_size)
                fwrite(pk_grid_cell->mesh_extra_info,    pk_grid_cell->mesh_extra_info_size,              1, output_file);

            pk_grid_cell = pk_grid_cell->next;
        }

        fclose(output_file);
    }

    // Secondly, we save the monodomain solver state
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

    // Finally, we save the ODE state variables from both [ode_solver] and [purkinje_ode_solver]
    if(the_ode_solver && the_purkinje_ode_solver) {

        sds tmp = sdsnew(output_dir);
        tmp = sdscat(tmp, "/ode_solver_checkpoint.dat");

        FILE *output_file = fopen(tmp, "wb");

        sdsfree(tmp);

        if (!output_file) {
            fprintf(stderr, "Error opening %s file for saving the simulation state\n", tmp);
            return;
        }

    // DOMAIN section
        fwrite(&(the_ode_solver->adaptive), sizeof(the_ode_solver->adaptive), 1, output_file);
        fwrite(&(the_ode_solver->max_dt), sizeof(the_ode_solver->max_dt), 1, output_file);
        fwrite(&(the_ode_solver->min_dt), sizeof(the_ode_solver->min_dt), 1, output_file);
        fwrite(&(the_ode_solver->rel_tol), sizeof(the_ode_solver->rel_tol), 1, output_file);
        fwrite(&(the_ode_solver->abs_tol), sizeof(the_ode_solver->abs_tol), 1, output_file);

        bool read_cells_to_solve = false;

        fwrite(&(the_ode_solver->num_cells_to_solve), sizeof(the_ode_solver->num_cells_to_solve), 1, output_file);

        if(the_ode_solver->cells_to_solve) {
            read_cells_to_solve = true;
            fwrite(the_ode_solver->cells_to_solve, sizeof(the_ode_solver->cells_to_solve[0]), the_ode_solver->num_cells_to_solve, output_file);
        }

        fwrite(&(read_cells_to_solve), sizeof(read_cells_to_solve), 1, output_file);

        fwrite(&(the_ode_solver->gpu), sizeof(the_ode_solver->gpu), 1, output_file);
        fwrite(&(the_ode_solver->gpu_id), sizeof(the_ode_solver->gpu_id), 1, output_file);

        fwrite(&(the_ode_solver->pitch), sizeof(the_ode_solver->pitch), 1, output_file);
        fwrite(&(the_ode_solver->original_num_cells), sizeof(the_ode_solver->original_num_cells), 1, output_file);

        size_t num_sv_entries = the_ode_solver->model_data.number_of_ode_equations;

        if(the_ode_solver->gpu) {

            #ifdef COMPILE_CUDA
            if(the_ode_solver->adaptive) {
                num_sv_entries = num_sv_entries + 3;
            }

            real *sv_cpu;
            sv_cpu = MALLOC_ARRAY_OF_TYPE(real, the_ode_solver->original_num_cells * num_sv_entries);

            check_cuda_error(cudaMemcpy2D(sv_cpu, the_ode_solver->original_num_cells * sizeof(real), the_ode_solver->sv, the_ode_solver->pitch,
                         the_ode_solver->original_num_cells * sizeof(real), num_sv_entries, cudaMemcpyDeviceToHost));

            fwrite(sv_cpu, sizeof(real), the_ode_solver->original_num_cells * num_sv_entries, output_file);
            #endif

        }
        else {
            fwrite(the_ode_solver->sv, sizeof(real), the_ode_solver->original_num_cells * num_sv_entries,  output_file);
        }

        //fwrite(&(the_ode_solver->extra_data_size), sizeof(the_ode_solver->extra_data_size), 1, output_file);
        //fwrite(the_ode_solver->ode_extra_data, the_ode_solver->extra_data_size, 1, output_file);

    // PURKINJE section
        fwrite(&(the_purkinje_ode_solver->adaptive), sizeof(the_purkinje_ode_solver->adaptive), 1, output_file);
        fwrite(&(the_purkinje_ode_solver->max_dt), sizeof(the_purkinje_ode_solver->max_dt), 1, output_file);
        fwrite(&(the_purkinje_ode_solver->min_dt), sizeof(the_purkinje_ode_solver->min_dt), 1, output_file);
        fwrite(&(the_purkinje_ode_solver->rel_tol), sizeof(the_purkinje_ode_solver->rel_tol), 1, output_file);
        fwrite(&(the_purkinje_ode_solver->abs_tol), sizeof(the_purkinje_ode_solver->abs_tol), 1, output_file);

        read_cells_to_solve = false;

        fwrite(&(the_purkinje_ode_solver->num_cells_to_solve), sizeof(the_purkinje_ode_solver->num_cells_to_solve), 1, output_file);

        if(the_purkinje_ode_solver->cells_to_solve) {
            read_cells_to_solve = true;
            fwrite(the_purkinje_ode_solver->cells_to_solve, sizeof(the_purkinje_ode_solver->cells_to_solve[0]), the_purkinje_ode_solver->num_cells_to_solve, output_file);
        }

        fwrite(&(read_cells_to_solve), sizeof(read_cells_to_solve), 1, output_file);

        fwrite(&(the_purkinje_ode_solver->gpu), sizeof(the_purkinje_ode_solver->gpu), 1, output_file);
        fwrite(&(the_purkinje_ode_solver->gpu_id), sizeof(the_purkinje_ode_solver->gpu_id), 1, output_file);

        fwrite(&(the_purkinje_ode_solver->pitch), sizeof(the_purkinje_ode_solver->pitch), 1, output_file);
        fwrite(&(the_purkinje_ode_solver->original_num_cells), sizeof(the_purkinje_ode_solver->original_num_cells), 1, output_file);

        num_sv_entries = the_purkinje_ode_solver->model_data.number_of_ode_equations;

        if(the_purkinje_ode_solver->gpu) {

            #ifdef COMPILE_CUDA
            if(the_purkinje_ode_solver->adaptive) {
                num_sv_entries = num_sv_entries + 3;
            }

            real *sv_cpu;
            sv_cpu = MALLOC_ARRAY_OF_TYPE(real, the_purkinje_ode_solver->original_num_cells * num_sv_entries);

            check_cuda_error(cudaMemcpy2D(sv_cpu, the_purkinje_ode_solver->original_num_cells * sizeof(real), the_purkinje_ode_solver->sv, the_purkinje_ode_solver->pitch,
                         the_purkinje_ode_solver->original_num_cells * sizeof(real), num_sv_entries, cudaMemcpyDeviceToHost));

            fwrite(sv_cpu, sizeof(real), the_purkinje_ode_solver->original_num_cells * num_sv_entries, output_file);
            #endif

        }
        else {
            fwrite(the_purkinje_ode_solver->sv, sizeof(real), the_purkinje_ode_solver->original_num_cells * num_sv_entries,  output_file);
        }

        //fwrite(&(the_purkinje_ode_solver->extra_data_size), sizeof(the_purkinje_ode_solver->extra_data_size), 1, output_file);
        //fwrite(the_purkinje_ode_solver->ode_extra_data, the_purkinje_ode_solver->extra_data_size, 1, output_file);

        fclose(output_file);
    }

    CALL_EXTRA_FUNCTIONS(save_state_fn, time_info, config, save_mesh_config, the_grid, the_monodomain_solver, the_ode_solver, the_purkinje_ode_solver, output_dir);
}

static void save_point_array_hash(struct point_voidp_hash_entry *p, FILE *output_file) {

        size_t n = hmlen(p);
        fwrite(&n, sizeof(n), 1, output_file);

        for(size_t i = 0; i < n; i++) {
            struct point_3d key = p[i].key;
            float  *value       = p[i].value;

            fwrite(&key, sizeof(key), 1, output_file);

            size_t n2 = arrlen(value);
            fwrite(&n2, sizeof(n2), 1, output_file);

            for(size_t j = 0; j < n2; j++) {
                fwrite(&value[j], sizeof(value[j]), 1, output_file);
            }

        }
}

static void save_point_float_hash(struct point_hash_entry *p, FILE *output_file) {

        size_t n = hmlen(p);
        fwrite(&n, sizeof(n), 1, output_file);

        for(size_t i = 0; i < n; i++) {
            struct point_3d key = p[i].key;
            float         value = p[i].value;
            fwrite(&key, sizeof(key), 1, output_file);
            fwrite(&value, sizeof(value), 1, output_file);
        }
}

static void save_activation_times_data(FILE *output_file, struct common_persistent_data* persistent_data) {
    save_point_float_hash(persistent_data->last_time_v, output_file);
    save_point_float_hash(persistent_data->num_activations, output_file);
    save_point_float_hash(persistent_data->cell_was_active, output_file);

    save_point_array_hash(persistent_data->activation_times, output_file);
    save_point_array_hash(persistent_data->apds, output_file);
}

static void save_ensight_data(FILE *output_file,  struct common_persistent_data* persistent_data ) {
    fwrite(&(persistent_data->file_count), sizeof(persistent_data->file_count), 1, output_file);
    fwrite(&(persistent_data->n_digits), sizeof(persistent_data->n_digits), 1, output_file);
}

SAVE_STATE(save_simulation_state_with_activation_times_extra_fn) {

    if(save_mesh_config->persistent_data) {
        sds tmp = sdsnew(output_dir);
        tmp = sdscat(tmp, "/persistent_data_checkpoint.dat");

        FILE *output_file = fopen(tmp, "wb");
        sdsfree(tmp);

        struct common_persistent_data* persistent_data = (struct common_persistent_data*)save_mesh_config->persistent_data;
        save_activation_times_data(output_file, persistent_data);
        fclose(output_file);

    }
}

SAVE_STATE(save_simulation_state_with_ensight_extra_fn) {

    if(save_mesh_config->persistent_data) {
        sds tmp = sdsnew(output_dir);
        tmp = sdscat(tmp, "/persistent_data_checkpoint.dat");

        FILE *output_file = fopen(tmp, "wb");
        sdsfree(tmp);

        struct common_persistent_data* persistent_data = (struct common_persistent_data*)save_mesh_config->persistent_data;
        save_ensight_data(output_file, persistent_data);
        fclose(output_file);

    }
}

SAVE_STATE(save_simulation_state_with_activation_times) {

    if(save_mesh_config->persistent_data) {
        sds tmp = sdsnew(output_dir);
        tmp = sdscat(tmp, "/persistent_data_checkpoint.dat");

        FILE *output_file = fopen(tmp, "wb");
        sdsfree(tmp);

        struct common_persistent_data* persistent_data = (struct common_persistent_data*)save_mesh_config->persistent_data;

        save_activation_times_data(output_file, persistent_data);

         char *mesh_format = NULL;
         GET_PARAMETER_STRING_VALUE_OR_USE_DEFAULT(mesh_format, save_mesh_config, "mesh_format");
         if(mesh_format != NULL && STRINGS_EQUAL(mesh_format, "ensight")) {
             save_ensight_data(output_file, persistent_data);
             free(mesh_format);
         }
         fclose(output_file);
    }

    save_simulation_state(time_info, config, save_mesh_config, the_grid, the_monodomain_solver, the_ode_solver, the_purkinje_ode_solver, output_dir);

}
