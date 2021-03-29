//
// Created by sachetto on 13/10/17.
//

#include <stdbool.h>
#include <stdlib.h>

#include "../3dparty/sds/sds.h"
#include "../alg/grid/grid.h"
#include "../config/save_mesh_config.h"
#include "../utils/utils.h"

#include "../vtk_utils/vtk_unstructured_grid.h"
#include "../domains_library/mesh_info_data.h"

#include "save_mesh_helper.h"

#ifdef COMPILE_CUDA
#include "../gpu_utils/gpu_utils.h"
#endif

static char *file_prefix;
static bool binary = false;
static bool clip_with_plain = false;
static bool clip_with_bounds = false;
static bool save_pvd = true;
static bool save_inactive = false;
static bool compress = false;
static bool save_f = false;
static int compression_level = 3;
char *output_dir;
bool save_visible_mask = true;

static bool initialized = false;

static void save_visibility_mask(sds output_dir_with_file, ui8_array visible_cells) {
        sds output_dir_with_new_file = sdsnew(output_dir_with_file);
        output_dir_with_new_file = sdscat(output_dir_with_new_file, ".vis");
		FILE *vis = fopen(output_dir_with_new_file, "wb");
		fwrite(visible_cells, sizeof(uint8_t), arrlen(visible_cells), vis);
        sdsfree(output_dir_with_new_file);
		fclose(vis);
}

SAVE_MESH(save_as_adjacency_list) {

    int iteration_count = time_info->iteration;

    if(!initialized) {
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config, "output_dir");
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(file_prefix, config, "file_prefix");
        initialized = true;
    }

    sds tmp = sdsnew(output_dir);
    tmp = sdscat(tmp, "/");

    sds base_name = NULL;
    if(binary) {
        base_name = create_base_name(file_prefix, iteration_count, "bin");
    }
    else {
        base_name = create_base_name(file_prefix, iteration_count, "txt");
    }

    tmp = sdscat(tmp, base_name);

    FILE *output_file = fopen(tmp, "w");

    sdsfree(base_name);
    sdsfree(tmp);

    struct cell_node *neighbour;

    FOR_EACH_CELL(the_grid) {

        if(cell->active) {

            fprintf(output_file, "%d ", cell->grid_position);

            neighbour = get_cell_neighbour(cell, cell->neighbours[FRONT]);
            if(neighbour) {
                fprintf(output_file, "%d ", neighbour->grid_position);
            }
            neighbour = get_cell_neighbour(cell, cell->neighbours[BACK]);
            if(neighbour) {
                fprintf(output_file, "%d ", neighbour->grid_position);
            }

            neighbour = get_cell_neighbour(cell, cell->neighbours[DOWN]);
            if(neighbour) {
                fprintf(output_file, "%d ", neighbour->grid_position);
            }

            neighbour = get_cell_neighbour(cell, cell->neighbours[TOP]);
            if(neighbour) {
                fprintf(output_file, "%d ", neighbour->grid_position);
            }

            neighbour = get_cell_neighbour(cell, cell->neighbours[RIGHT]);
            if(neighbour) {
                fprintf(output_file, "%d ", neighbour->grid_position);
            }

            neighbour = get_cell_neighbour(cell, cell->neighbours[LEFT]);
            if(neighbour) {
                fprintf(output_file, "%d", neighbour->grid_position);
            }

            fprintf(output_file, "\n");
        }
    }

    fclose(output_file);


}

INIT_SAVE_MESH(init_save_one_cell_state_variables) {
    config->persistent_data = malloc(sizeof(struct save_one_cell_state_variables_persistent_data));
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR( ((struct save_one_cell_state_variables_persistent_data *) config->persistent_data)->file_name, config, "file_name");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, ((struct save_one_cell_state_variables_persistent_data *) config->persistent_data)->cell_center_x, config, "cell_center_x");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, ((struct save_one_cell_state_variables_persistent_data *) config->persistent_data)->cell_center_y, config, "cell_center_y");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, ((struct save_one_cell_state_variables_persistent_data *) config->persistent_data)->cell_center_z, config, "cell_center_z");

    ((struct save_one_cell_state_variables_persistent_data *) config->persistent_data)->file = fopen(((struct save_one_cell_state_variables_persistent_data *) config->persistent_data)->file_name, "w");
    ((struct save_one_cell_state_variables_persistent_data *) config->persistent_data)->cell_sv_position = -1;

}

SAVE_MESH(save_one_cell_state_variables) {

    struct save_one_cell_state_variables_persistent_data *params = ((struct save_one_cell_state_variables_persistent_data *) config->persistent_data);

    if(params->cell_sv_position == -1) {
        if(!the_grid->adaptive) {
            FOR_EACH_CELL(the_grid) {
                if(cell->center.x == params->cell_center_x && cell->center.y == params->cell_center_y && cell->center.z == params->cell_center_z) {
                    params->cell_sv_position = cell->sv_position;
                    break;
                }
            }
        }
    }

    if(ode_solver->gpu) {
#ifdef COMPILE_CUDA
            real *cell_sv;

            cell_sv = MALLOC_ARRAY_OF_TYPE(real, ode_solver->model_data.number_of_ode_equations);

            check_cuda_error(cudaMemcpy2D(cell_sv, sizeof(real), ode_solver->sv + params->cell_sv_position,
                                           ode_solver->pitch, sizeof(real),
                                           ode_solver->model_data.number_of_ode_equations, cudaMemcpyDeviceToHost));

            fprintf(params->file, "%lf %lf ", time_info->current_t, cell_sv[0]);
            for(int i = 2; i < 11; i++) {
                fprintf(params->file, " %lf ", cell_sv[i]);
            }

            fprintf(params->file, " %lf ", cell_sv[13]);
            fprintf(params->file, " %lf ", cell_sv[11]);
            fprintf(params->file, " %lf ", cell_sv[12]);

            for(int i = 14; i < 30; i++) {
                fprintf(params->file, " %lf ", cell_sv[i]);
            }

            fprintf(params->file, " %lf ", cell_sv[32]);
            fprintf(params->file, " %lf ", cell_sv[33]);

            fprintf(params->file, " %lf ", cell_sv[30]);
            fprintf(params->file, " %lf ", cell_sv[31]);

            fprintf(params->file, " %lf ", cell_sv[39]);
            fprintf(params->file, " %lf ", cell_sv[40]);
            fprintf(params->file, " %lf ", cell_sv[41]);
            fprintf(params->file, " %lf ", cell_sv[1]);

            fprintf(params->file, " %lf ", cell_sv[34]);
            fprintf(params->file, " %lf ", cell_sv[35]);
            fprintf(params->file, " %lf ", cell_sv[36]);
            fprintf(params->file, " %lf ", cell_sv[37]);
            fprintf(params->file, " %lf ", cell_sv[38]);
            fprintf(params->file, " %lf\n", cell_sv[42]);
            free(cell_sv);
#endif
    }
    else {

        real *cell_sv =  &ode_solver->sv[params->cell_sv_position * ode_solver->model_data.number_of_ode_equations];

        fprintf(params->file, "%lf %lf ", time_info->current_t, cell_sv[0]);
        for(int i = 2; i < 11; i++) {
            fprintf(params->file, " %lf ", cell_sv[i]);
        }

        fprintf(params->file, " %lf ", cell_sv[13]);
        fprintf(params->file, " %lf ", cell_sv[11]);
        fprintf(params->file, " %lf ", cell_sv[12]);

        for(int i = 14; i < 30; i++) {
            fprintf(params->file, " %lf ", cell_sv[i]);
        }

        fprintf(params->file, " %lf ", cell_sv[32]);
        fprintf(params->file, " %lf ", cell_sv[33]);

        fprintf(params->file, " %lf ", cell_sv[30]);
        fprintf(params->file, " %lf ", cell_sv[31]);

        fprintf(params->file, " %lf ", cell_sv[39]);
        fprintf(params->file, " %lf ", cell_sv[40]);
        fprintf(params->file, " %lf ", cell_sv[41]);
        fprintf(params->file, " %lf ", cell_sv[1]);

        fprintf(params->file, " %lf ", cell_sv[34]);
        fprintf(params->file, " %lf ", cell_sv[35]);
        fprintf(params->file, " %lf ", cell_sv[36]);
        fprintf(params->file, " %lf ", cell_sv[37]);
        fprintf(params->file, " %lf ", cell_sv[38]);
        fprintf(params->file, " %lf\n", cell_sv[42]);


    }

}

END_SAVE_MESH(end_save_one_cell_state_variables) {
    free(((struct save_one_cell_state_variables_persistent_data *) config->persistent_data)->file_name);
    fclose(((struct save_one_cell_state_variables_persistent_data *) config->persistent_data)->file);
    free(config->persistent_data);
}

SAVE_MESH(save_as_text_or_binary) {

    int iteration_count = time_info->iteration;

//    if(!initialized) {
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config, "output_dir");
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(file_prefix, config, "file_prefix");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(binary, config, "binary");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_plain, config, "clip_with_plain");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_bounds, config, "clip_with_bounds");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_inactive, config, "save_inactive_cells");
		GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_visible_mask, config, "save_visible_mask");
        initialized = true;
//    }

    real_cpu min_x = 0.0;
    real_cpu min_y = 0.0;
    real_cpu min_z = 0.0;
    real_cpu max_x = 0.0;
    real_cpu max_y = 0.0;
    real_cpu max_z = 0.0;

    float p0[3] = {0, 0, 0};
    float n[3] = {0, 0, 0};

    if(clip_with_plain) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, n[0], config, "normal_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, n[1], config, "normal_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, n[2], config, "normal_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, p0[0], config, "origin_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, p0[1], config, "origin_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, p0[2], config, "origin_z");
    }

    if(clip_with_bounds) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, min_x, config, "min_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, min_y, config, "min_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, min_z, config, "min_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, max_x, config, "max_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, max_y, config, "max_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, max_z, config, "max_z");
    }

    real_cpu l = sqrtf(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    real_cpu A = n[0] / l;
    real_cpu B = n[1] / l;
    real_cpu C = n[2] / l;
    real_cpu D = -(n[0] * p0[0] + n[1] * p0[1] + n[2] * p0[2]);

    real_cpu side;

    sds tmp = sdsnew(output_dir);
    tmp = sdscat(tmp, "/");

    sds base_name = NULL;
    if(binary) {
        base_name = create_base_name(file_prefix, iteration_count, "bin");
    }
    else {
        base_name = create_base_name(file_prefix, iteration_count, "txt");
    }

    tmp = sdscat(tmp, base_name);

    FILE *output_file = fopen(tmp, "w");

       struct cell_node *grid_cell = the_grid->first_cell;

    real_cpu center_x, center_y, center_z, dx, dy, dz;
    real_cpu v;

	ui8_array cell_visibility = NULL;
	arrsetcap(cell_visibility, the_grid->num_active_cells);

    while(grid_cell != 0) {

        if(grid_cell->active || save_inactive) {

            center_x = grid_cell->center.x;
            center_y = grid_cell->center.y;
            center_z = grid_cell->center.z;

            if(clip_with_plain) {
                side = A * center_x + B * center_y + C * center_z + D;
                if(side < 0) {
                    grid_cell = grid_cell->next;
                    continue;
                }
            }

            if(clip_with_bounds) {
                bool ignore_cell = center_x < min_x || center_x > max_x || center_y < min_y || center_y > max_y ||
                                   center_z < min_z || center_z > max_z;

                if(ignore_cell) {
                    grid_cell = grid_cell->next;
                    continue;
                }
            }

            v = grid_cell->v;
            dx = grid_cell->discretization.x/2.0;
            dy = grid_cell->discretization.y/2.0;
            dz = grid_cell->discretization.z/2.0;

            if(binary) {
                fwrite(&center_x, sizeof(center_x), 1, output_file);
                fwrite(&center_y, sizeof(center_y), 1, output_file);
                fwrite(&center_z, sizeof(center_z), 1, output_file);
                fwrite(&dx, sizeof(dx), 1, output_file);
                fwrite(&dy, sizeof(dy), 1, output_file);
                fwrite(&dz, sizeof(dz), 1, output_file);
                fwrite(&v, sizeof(v), 1, output_file);
            } else {
                fprintf(output_file, "%g,%g,%g,%g,%g,%g,%g\n", center_x, center_y, center_z, dx, dy, dz, v);
            }
			arrput(cell_visibility, grid_cell->visible);
        }
        grid_cell = grid_cell->next;
    }

	if(save_visible_mask) {
		save_visibility_mask(tmp, cell_visibility);
	}

   	sdsfree(base_name);
    sdsfree(tmp);


    fclose(output_file);
}

struct save_as_vtk_or_vtu_persistent_data {
    struct vtk_unstructured_grid * grid;
    bool first_save_call;
};

INIT_SAVE_MESH(init_save_as_vtk_or_vtu) {
    config->persistent_data = malloc(sizeof(struct save_as_vtk_or_vtu_persistent_data));
    ((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->grid = NULL;
    ((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->first_save_call = true;
}

END_SAVE_MESH(end_save_as_vtk_or_vtu) {
    free_vtk_unstructured_grid(((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->grid);
    free(config->persistent_data);
}

SAVE_MESH(save_as_vtk) {

    int iteration_count = time_info->iteration;

    if(((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->first_save_call) {
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config, "output_dir");
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(file_prefix, config, "file_prefix");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_plain, config, "clip_with_plain");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_bounds, config, "clip_with_bounds");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(binary, config, "binary");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_f, config, "save_f");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_visible_mask, config, "save_visible_mask");

        ((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->first_save_call = false;

    }
    float plain_coords[6] = {0, 0, 0, 0, 0, 0};
    float bounds[6] = {0, 0, 0, 0, 0, 0};

    if(clip_with_plain) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[0], config, "origin_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[1], config, "origin_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[2], config, "origin_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[3], config, "normal_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[4], config, "normal_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[5], config, "normal_z");
    }

    if(clip_with_bounds) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[0], config, "min_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[1], config, "min_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[2], config, "min_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[3], config, "max_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[4], config, "max_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[5], config, "max_z");
    }

    sds output_dir_with_file = sdsnew(output_dir);
    output_dir_with_file = sdscat(output_dir_with_file, "/");
    sds base_name = create_base_name(file_prefix, iteration_count, "vtk");

    real_cpu current_t = time_info->current_t;

    //TODO: change this. We dont need the current_t here
    output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, current_t);

    bool read_only_data = ((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->grid != NULL;

	new_vtk_unstructured_grid_from_alg_grid(&(((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->grid), the_grid, clip_with_plain, plain_coords, clip_with_bounds, bounds, read_only_data, save_f);

    save_vtk_unstructured_grid_as_legacy_vtk(((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->grid, output_dir_with_file, binary, save_f);

	if(save_visible_mask) {
		save_visibility_mask(output_dir_with_file, (((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->grid)->cell_visibility);
	}	

    if(the_grid->adaptive) {
        free_vtk_unstructured_grid(((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->grid);
        ((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->grid = NULL;
    }

    sdsfree(output_dir_with_file);
    sdsfree(base_name);
}

SAVE_MESH(save_as_vtu) {

    int iteration_count = time_info->iteration;

    if(((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->first_save_call) {
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config, "output_dir");
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(file_prefix, config, "file_prefix");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_plain, config, "clip_with_plain");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_bounds, config, "clip_with_bounds");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(binary, config, "binary");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_pvd, config, "save_pvd");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(compress, config, "compress");
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, compression_level, config, "compression_level");
		GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_visible_mask, config, "save_visible_mask");

        if(compress) binary = true;

        if(!save_pvd) {
            ((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->first_save_call = false;
        }
    }

    float plain_coords[6] = {0, 0, 0, 0, 0, 0};
    float bounds[6] = {0, 0, 0, 0, 0, 0};

    if(clip_with_plain) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[0], config, "origin_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[1], config, "origin_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[2], config, "origin_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[3], config, "normal_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[4], config, "normal_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[5], config, "normal_z");
    }

    if(clip_with_bounds) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[0], config, "min_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[1], config, "min_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[2], config, "min_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[3], config, "max_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[4], config, "max_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[5], config, "max_z");
    }


    sds output_dir_with_file = sdsnew(output_dir);
    output_dir_with_file = sdscat(output_dir_with_file, "/");
    sds base_name = create_base_name(file_prefix, iteration_count, "vtu");

    real_cpu current_t = time_info->current_t;

    output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, current_t);

    if(save_pvd) {
        add_file_to_pvd(current_t, output_dir, base_name, ((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->first_save_call);
        ((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->first_save_call = false;
    }

    bool read_only_data = ((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->grid != NULL;
    new_vtk_unstructured_grid_from_alg_grid(&((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->grid, the_grid, clip_with_plain, plain_coords, clip_with_bounds, bounds, read_only_data, save_f);

    if(compress) {
        save_vtk_unstructured_grid_as_vtu_compressed(((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->grid, output_dir_with_file, compression_level);
    }
    else {
        save_vtk_unstructured_grid_as_vtu(((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->grid, output_dir_with_file, binary);
    }

	if(save_visible_mask) {
		save_visibility_mask(output_dir_with_file, (((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->grid)->cell_visibility);
	}

    sdsfree(output_dir_with_file);
    sdsfree(base_name);

    //TODO: I do not know if we should to this here or call the end and init save functions on the adaptivity step.....
    if(the_grid->adaptive) {
        free_vtk_unstructured_grid(((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->grid);
        ((struct save_as_vtk_or_vtu_persistent_data *) config->persistent_data)->grid = NULL;
    }

}

struct save_with_activation_times_persistent_data {
    struct point_hash_entry *last_time_v;
    struct point_hash_entry *num_activations;
    struct point_hash_entry *cell_was_active;
    struct point_voidp_hash_entry *activation_times;
    struct point_voidp_hash_entry *apds;
    bool first_save_call;

};

INIT_SAVE_MESH(init_save_with_activation_times) {

    config->persistent_data = calloc(1, sizeof(struct save_with_activation_times_persistent_data));
    hmdefault(((struct save_with_activation_times_persistent_data*)config->persistent_data)->cell_was_active, 0.0);
    hmdefault(((struct save_with_activation_times_persistent_data*)config->persistent_data)->last_time_v, -100.0);
    hmdefault(((struct save_with_activation_times_persistent_data*)config->persistent_data)->num_activations, 0);
    hmdefault(((struct save_with_activation_times_persistent_data*)config->persistent_data)->activation_times, NULL);
    hmdefault(((struct save_with_activation_times_persistent_data*)config->persistent_data)->apds, NULL);
    ((struct save_with_activation_times_persistent_data*)config->persistent_data)->first_save_call = true;

}

END_SAVE_MESH(end_save_with_activation_times) {
    free(config->persistent_data);
}

SAVE_MESH(save_with_activation_times) {

    int mesh_output_pr = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, mesh_output_pr, config, "mesh_print_rate");

    int iteration_count = time_info->iteration;

    if(mesh_output_pr) {
        if (iteration_count % mesh_output_pr == 0)
            save_as_text_or_binary(time_info, config, the_grid, ode_solver);
    }

    float time_threshold = 10.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, time_threshold, config, "time_threshold");

    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config, "output_dir");

    float activation_threshold = -30.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, activation_threshold, config, "activation_threshold");

    float apd_threshold = -83.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, apd_threshold, config, "apd_threshold");

    real_cpu current_t = time_info->current_t;
    real_cpu last_t = time_info->final_t;
    real_cpu dt = time_info->dt;

    sds output_dir_with_file = sdsnew(output_dir);
    output_dir_with_file = sdscat(output_dir_with_file, "/");
    sds base_name = create_base_name("activation_info", 0, "acm");
    output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, current_t);

    struct save_with_activation_times_persistent_data *persistent_data =
            (struct save_with_activation_times_persistent_data*)config->persistent_data;

    struct cell_node *grid_cell = the_grid->first_cell;

    real_cpu center_x, center_y, center_z, dx, dy, dz;
    real_cpu v;

    FILE *act_file = fopen(output_dir_with_file, "w");

    fprintf(act_file, "%d\n", (last_t-current_t) <= dt ); //rounding errors

    while(grid_cell != 0) {

        if( grid_cell->active || ( grid_cell->mesh_extra_info && ( FIBROTIC(grid_cell) || BORDER_ZONE(grid_cell) ) ) ) {

            center_x = grid_cell->center.x;
            center_y = grid_cell->center.y;
            center_z = grid_cell->center.z;

            v = grid_cell->v;

            struct point_3d cell_coordinates;
            cell_coordinates.x = center_x;
            cell_coordinates.y = center_y;
            cell_coordinates.z = center_z;

            dx = grid_cell->discretization.x / 2.0;
            dy = grid_cell->discretization.y / 2.0;
            dz = grid_cell->discretization.z / 2.0;

            if(grid_cell->mesh_extra_info) {
                fprintf(act_file, "%g,%g,%g,%g,%g,%g,%d,%d,%d ", center_x, center_y, center_z, dx, dy, dz, grid_cell->active, FIBROTIC(grid_cell),
                        BORDER_ZONE(grid_cell));
            }
            else {
                fprintf(act_file, "%g,%g,%g,%g,%g,%g ", center_x, center_y, center_z, dx, dy, dz);
            }

            int n_activations = 0;
            float *apds_array = NULL;
            float *activation_times_array = NULL;

            if(grid_cell->active) {
                float last_v = hmget(persistent_data->last_time_v, cell_coordinates);

                n_activations = (int) hmget(persistent_data->num_activations, cell_coordinates);
                activation_times_array = (float *) hmget(persistent_data->activation_times, cell_coordinates);
                apds_array = (float *) hmget(persistent_data->apds, cell_coordinates);

                int act_times_len = arrlen(activation_times_array);

                if (current_t == 0.0f) {
                    hmput(persistent_data->last_time_v, cell_coordinates, v);
                } else {
                    if ((last_v < activation_threshold) && (v >= activation_threshold)) {

                        if (act_times_len == 0) {
                            n_activations++;
                            hmput(persistent_data->num_activations, cell_coordinates, n_activations);
                                    arrput(activation_times_array, current_t);
                            float tmp = hmget(persistent_data->cell_was_active, cell_coordinates);
                            hmput(persistent_data->cell_was_active, cell_coordinates, tmp + 1);
                            hmput(persistent_data->activation_times, cell_coordinates, activation_times_array);
                        } else { //This is to avoid spikes in the middle of an Action Potential
                            float last_act_time = activation_times_array[act_times_len - 1];
                            if (current_t - last_act_time > time_threshold) {
                                n_activations++;
                                hmput(persistent_data->num_activations, cell_coordinates, n_activations);
                                        arrput(activation_times_array, current_t);
                                float tmp = hmget(persistent_data->cell_was_active, cell_coordinates);
                                hmput(persistent_data->cell_was_active, cell_coordinates, tmp + 1);
                                hmput(persistent_data->activation_times, cell_coordinates, activation_times_array);
                            }
                        }
                    }

                    //CHECK APD
                    bool was_active = (hmget(persistent_data->cell_was_active, cell_coordinates) != 0.0);
                    if (was_active) {
                        if (v <= apd_threshold || (hmget(persistent_data->cell_was_active, cell_coordinates) == 2.0) || (last_t-current_t) <= dt) {

                            int tmp = (int)hmget(persistent_data->cell_was_active, cell_coordinates);
                            int act_time_array_len = arrlen(activation_times_array);
                            //if this in being calculated because we had a new activation before the cell achieved the rest potential,
                            // we need to get the activation before this one
                            real_cpu last_act_time = activation_times_array[act_time_array_len  - tmp];
                            real_cpu apd = current_t - last_act_time;
                                    arrput(apds_array, apd);
                            hmput(persistent_data->apds, cell_coordinates, apds_array);
                            hmput(persistent_data->cell_was_active, cell_coordinates, tmp - 1);
                        }
                    }

                    hmput(persistent_data->last_time_v, cell_coordinates, v);
                }
            }

            fprintf(act_file, "%d [ ", n_activations);

            for (unsigned long i = 0; i < n_activations; i++) {
                fprintf(act_file, "%lf ", activation_times_array[i]);
            }
            fprintf(act_file, "] ");

            fprintf(act_file, "[ ");

            for (unsigned long i = 0; i < arrlen(apds_array); i++) {
                fprintf(act_file, "%lf ", apds_array[i]);
            }
            fprintf(act_file, "]\n");
        }

        grid_cell = grid_cell->next;
    }

    fclose(act_file);


}

SAVE_MESH(no_save) {
    //Nop
}
