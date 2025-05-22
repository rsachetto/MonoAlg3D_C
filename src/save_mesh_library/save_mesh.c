//
// Created by sachetto on 13/10/17.
//

#include <stdbool.h>
#include <stdlib.h>

#include "../3dparty/sds/sds.h"
#include "../alg/grid/grid.h"
#include "../config/save_mesh_config.h"
#include "../extra_data_library/helper_functions.h"
#include "../utils/utils.h"

#include "../domains_library/mesh_info_data.h"
#include "../ensight_utils/ensight_grid.h"
#include "../vtk_utils/vtk_unstructured_grid.h"

#include "save_mesh_helper.h"

#if defined(COMPILE_CUDA) || defined(COMPILE_SYCL)
#define COMPILE_GPU
#endif

#ifdef COMPILE_GPU
#include "../gpu_utils/accel_utils.h"
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
bool save_scar_cells = false;
bool save_purkinje = true;
static bool initialized = false;
static bool save_ode_state_variables = false;

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
    } else {
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
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(((struct save_one_cell_state_variables_persistent_data *)config->persistent_data)->file_name, config,
                                               "file_name");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, ((struct save_one_cell_state_variables_persistent_data *)config->persistent_data)->cell_center_x,
                                                config, "cell_center_x");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, ((struct save_one_cell_state_variables_persistent_data *)config->persistent_data)->cell_center_y,
                                                config, "cell_center_y");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, ((struct save_one_cell_state_variables_persistent_data *)config->persistent_data)->cell_center_z,
                                                config, "cell_center_z");

    ((struct save_one_cell_state_variables_persistent_data *)config->persistent_data)->file =
        fopen(((struct save_one_cell_state_variables_persistent_data *)config->persistent_data)->file_name, "w");
    ((struct save_one_cell_state_variables_persistent_data *)config->persistent_data)->cell_sv_position = -1;
}

SAVE_MESH(save_one_cell_state_variables) {

    struct save_one_cell_state_variables_persistent_data *params = ((struct save_one_cell_state_variables_persistent_data *)config->persistent_data);

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
#ifdef COMPILE_GPU
        real *cell_sv;

        cell_sv = MALLOC_ARRAY_OF_TYPE(real, ode_solver->model_data.number_of_ode_equations);

        //        check_cuda_error(cudaMemcpy2D(cell_sv, sizeof(real), ode_solver->sv + params->cell_sv_position, ode_solver->pitch, sizeof(real),
        //                                     ode_solver->model_data.number_of_ode_equations, cudaMemcpyDeviceToHost));

        memcpy2d_device(cell_sv, sizeof(real), ode_solver->sv + params->cell_sv_position, ode_solver->pitch, sizeof(real),
                        ode_solver->model_data.number_of_ode_equations, DEVICE_TO_HOST);

        // All state variables
        for(uint32_t i = 0; i < ode_solver->model_data.number_of_ode_equations; i++) {
            fprintf(params->file, "%g, ", cell_sv[i]);
        }
        fprintf(params->file, "\n");

        free(cell_sv);
#endif
    } else {

        real *cell_sv = &ode_solver->sv[params->cell_sv_position * ode_solver->model_data.number_of_ode_equations];

        // Time and transmembrane potential
        // fprintf(params->file, "%g %g\n", time_info->current_t, cell_sv[0]);

        // Time, Cai and Vm
        // fprintf(params->file, "%g %g %g\n", time_info->current_t, cell_sv[0], cell_sv[5]);

        // Only transmembrane potential
        // fprintf(params->file, "%g\n", cell_sv[0]);

        // All state variables
        for(uint32_t i = 0; i < ode_solver->model_data.number_of_ode_equations; i++) {
            fprintf(params->file, "%g, ", cell_sv[i]);
        }
        fprintf(params->file, "\n");

        // Time, Cai and Vm at certain timestep
        // if (time_info->current_t >= 15200) {
        //    fprintf(params->file, "%.3lf %g %g\n", time_info->current_t, cell_sv[0], cell_sv[5]);
        //}

        // All state-variables at certain timestep
        // if (time_info->current_t >= 15200) {
        //    for (uint32_t i = 0; i < ode_solver->model_data.number_of_ode_equations; i++) {
        //        fprintf(params->file, "%g, ", cell_sv[i]);
        //    }
        //    fprintf(params->file, "\n");
        //}
    }
}

END_SAVE_MESH(end_save_one_cell_state_variables) {
    free(((struct save_one_cell_state_variables_persistent_data *)config->persistent_data)->file_name);
    fclose(((struct save_one_cell_state_variables_persistent_data *)config->persistent_data)->file);
    free(config->persistent_data);
    config->persistent_data = NULL;
}

SAVE_MESH(save_as_text_or_binary) {

    int iteration_count = time_info->iteration;

    real_cpu min_x = 0.0;
    real_cpu min_y = 0.0;
    real_cpu min_z = 0.0;
    real_cpu max_x = 0.0;
    real_cpu max_y = 0.0;
    real_cpu max_z = 0.0;

    float p0[3] = {1, 1, 1};
    float n[3] = {1, 1, 1};

    //    if(!initialized) {
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config, "output_dir");
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(file_prefix, config, "file_prefix");
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(binary, config, "binary");
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_plain, config, "clip_with_plain");
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_bounds, config, "clip_with_bounds");
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_inactive, config, "save_inactive_cells");
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_visible_mask, config, "save_visible_mask");
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_ode_state_variables, config, "save_ode_state_variables");

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

    initialized = true;
    //   }

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
    } else {
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
                bool ignore_cell = center_x < min_x || center_x > max_x || center_y < min_y || center_y > max_y || center_z < min_z || center_z > max_z;

                if(ignore_cell) {
                    grid_cell = grid_cell->next;
                    continue;
                }
            }

            v = grid_cell->v;
            dx = grid_cell->discretization.x / 2.0;
            dy = grid_cell->discretization.y / 2.0;
            dz = grid_cell->discretization.z / 2.0;

            if(binary) {
                // TODO: maybe the size of the data should be always fixed (double as instance)
                fwrite(&center_x, sizeof(center_x), 1, output_file);
                fwrite(&center_y, sizeof(center_y), 1, output_file);
                fwrite(&center_z, sizeof(center_z), 1, output_file);
                fwrite(&dx, sizeof(dx), 1, output_file);
                fwrite(&dy, sizeof(dy), 1, output_file);
                fwrite(&dz, sizeof(dz), 1, output_file);
                fwrite(&v, sizeof(v), 1, output_file);
            } else {
                if(save_ode_state_variables) {

                    int n_state_vars = ode_solver->model_data.number_of_ode_equations - 1; // Vm is always saved
                    size_t num_sv_entries = n_state_vars + 1;
                    real *sv_cpu;

                    if(ode_solver->gpu) {

#ifdef COMPILE_GPU
                        sv_cpu = MALLOC_ARRAY_OF_TYPE(real, ode_solver->original_num_cells * num_sv_entries);
                        // check_cuda_error(cudaMemcpy2D(sv_cpu, ode_solver->original_num_cells * sizeof(real), ode_solver->sv, ode_solver->pitch,
                        //                              ode_solver->original_num_cells * sizeof(real), num_sv_entries, cudaMemcpyDeviceToHost));
                        memcpy2d_device(sv_cpu, ode_solver->original_num_cells * sizeof(real), ode_solver->sv, ode_solver->pitch,
                                        ode_solver->original_num_cells * sizeof(real), num_sv_entries, DEVICE_TO_HOST);
#endif
                    } else {
                        sv_cpu = ode_solver->sv;
                    }

                    fprintf(output_file, "%g,%g,%g,%g,%g,%g,%g", center_x, center_y, center_z, dx, dy, dz, v);

                    for(int i = 1; i <= n_state_vars; i++) {
                        float value;
                        if(ode_solver->gpu) {
                            value = (float)sv_cpu[i * ode_solver->original_num_cells];
                        } else {
                            value = sv_cpu[i];
                        }

                        fprintf(output_file, ",%g", value);
                    }

                    if(ode_solver->gpu) {
                        free(sv_cpu);
                    }

                    fprintf(output_file, "\n");
                } else {
                    fprintf(output_file, "%g,%g,%g,%g,%g,%g,%g\n", center_x, center_y, center_z, dx, dy, dz, v);
                }
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

    CALL_EXTRA_FUNCTIONS(save_mesh_fn, time_info, config, the_grid, ode_solver, purkinje_ode_solver);
}

INIT_SAVE_MESH(init_save_as_vtk_or_vtu) {
    if(config->persistent_data == NULL) {
        config->persistent_data = malloc(sizeof(struct common_persistent_data));
        ((struct common_persistent_data *)config->persistent_data)->grid = NULL;
        ((struct common_persistent_data *)config->persistent_data)->first_save_call = true;
    }
}

END_SAVE_MESH(end_save_as_vtk_or_vtu) {
    free_vtk_unstructured_grid(((struct common_persistent_data *)config->persistent_data)->grid);
    free(config->persistent_data);
    config->persistent_data = NULL;
}

SAVE_MESH(save_as_vtk) {

    int iteration_count = time_info->iteration;

    if(((struct common_persistent_data *)config->persistent_data)->first_save_call) {
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config, "output_dir");
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(file_prefix, config, "file_prefix");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_plain, config, "clip_with_plain");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_bounds, config, "clip_with_bounds");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(binary, config, "binary");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_f, config, "save_f");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_visible_mask, config, "save_visible_mask");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_scar_cells, config, "save_scar_cells");

        ((struct common_persistent_data *)config->persistent_data)->first_save_call = false;
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

    // TODO: change this. We dont need the current_t here
    output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, current_t);

    bool read_only_data = ((struct common_persistent_data *)config->persistent_data)->grid != NULL;

    new_vtk_unstructured_grid_from_alg_grid(&(((struct common_persistent_data *)config->persistent_data)->grid), the_grid, clip_with_plain, plain_coords,
                                            clip_with_bounds, bounds, read_only_data, save_f, save_scar_cells, NULL);

    save_vtk_unstructured_grid_as_legacy_vtk(((struct common_persistent_data *)config->persistent_data)->grid, output_dir_with_file, binary, save_f, NULL);

    if(save_visible_mask) {
        save_visibility_mask(output_dir_with_file, (((struct common_persistent_data *)config->persistent_data)->grid)->cell_visibility);
    }

    if(the_grid->adaptive) {
        free_vtk_unstructured_grid(((struct common_persistent_data *)config->persistent_data)->grid);
        ((struct common_persistent_data *)config->persistent_data)->grid = NULL;
    }

    sdsfree(output_dir_with_file);
    sdsfree(base_name);

    CALL_EXTRA_FUNCTIONS(save_mesh_fn, time_info, config, the_grid, ode_solver, purkinje_ode_solver);
}

SAVE_MESH(save_as_vtu) {

    int iteration_count = time_info->iteration;

    if(((struct common_persistent_data *)config->persistent_data)->first_save_call) {
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config, "output_dir");
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(file_prefix, config, "file_prefix");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_plain, config, "clip_with_plain");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_bounds, config, "clip_with_bounds");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(binary, config, "binary");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_pvd, config, "save_pvd");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(compress, config, "compress");
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, compression_level, config, "compression_level");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_visible_mask, config, "save_visible_mask");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_scar_cells, config, "save_scar_cells");

        if(compress)
            binary = true;

        if(!save_pvd) {
            ((struct common_persistent_data *)config->persistent_data)->first_save_call = false;
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
        add_file_to_pvd(current_t, output_dir, base_name, ((struct common_persistent_data *)config->persistent_data)->first_save_call);
        ((struct common_persistent_data *)config->persistent_data)->first_save_call = false;
    }

    bool read_only_data = ((struct common_persistent_data *)config->persistent_data)->grid != NULL;
    new_vtk_unstructured_grid_from_alg_grid(&((struct common_persistent_data *)config->persistent_data)->grid, the_grid, clip_with_plain, plain_coords,
                                            clip_with_bounds, bounds, read_only_data, save_f, save_scar_cells, NULL);

    if(compress) {
        save_vtk_unstructured_grid_as_vtu_compressed(((struct common_persistent_data *)config->persistent_data)->grid, output_dir_with_file, compression_level);
    } else {
        save_vtk_unstructured_grid_as_vtu(((struct common_persistent_data *)config->persistent_data)->grid, output_dir_with_file, binary);
    }

    if(save_visible_mask) {
        save_visibility_mask(output_dir_with_file, (((struct common_persistent_data *)config->persistent_data)->grid)->cell_visibility);
    }

    sdsfree(output_dir_with_file);
    sdsfree(base_name);

    // TODO: I do not know if we should to this here or call the end and init save functions on the adaptivity step.....
    if(the_grid->adaptive) {
        free_vtk_unstructured_grid(((struct common_persistent_data *)config->persistent_data)->grid);
        ((struct common_persistent_data *)config->persistent_data)->grid = NULL;
    }

    CALL_EXTRA_FUNCTIONS(save_mesh_fn, time_info, config, the_grid, ode_solver, purkinje_ode_solver);
}

INIT_SAVE_MESH(init_save_as_ensight) {
    if(config->persistent_data == NULL) {
        config->persistent_data = calloc(1, sizeof(struct common_persistent_data));
    }
}

SAVE_MESH(save_as_ensight) {

    struct common_persistent_data *persistent_data = (struct common_persistent_data *)config->persistent_data;

    if(the_grid == NULL) {
        log_error_and_exit("Error in save_as_ensight. No grid defined\n");
    }

    if(the_grid->num_active_cells == 0 && the_grid->purkinje == NULL) {
        log_error_and_exit("Error in save_as_ensight. No grid and/or no purkinje grid defined\n");
    }

    if(the_grid->adaptive) {
        log_error_and_exit("save_as_ensight function does not support adaptive meshes yet! Aborting\n");
    }

    if(!persistent_data->geometry_saved) {

        int print_rate = 1;

        char *mesh_format = NULL;
        GET_PARAMETER_STRING_VALUE_OR_USE_DEFAULT(mesh_format, config, "mesh_format");

        // We are getting called from save_with_activation_times
        if(mesh_format != NULL) {
            GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, print_rate, config, "mesh_print_rate");
        } else { // We are being directly called
            GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, print_rate, config, "print_rate");
        }

        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config, "output_dir");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(binary, config, "binary");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_visible_mask, config, "save_visible_mask");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_ode_state_variables, config, "save_ode_state_variables");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_purkinje, config, "save_purkinje");

        persistent_data->num_files = ((time_info->final_t / time_info->dt) / print_rate) + 1;

        sds output_dir_with_file = sdsnew(output_dir);
        output_dir_with_file = sdscat(output_dir_with_file, "/geometry.geo");

        struct ensight_grid *ensight_grid = new_ensight_grid_from_alg_grid(the_grid, false, NULL, false, NULL, false, false, save_purkinje);
        save_ensight_grid_as_ensight6_geometry(ensight_grid, output_dir_with_file, binary);

        if(save_visible_mask) {
            save_visibility_mask(output_dir_with_file, ensight_grid->parts[0].cell_visibility);
        }

        free_ensight_grid(ensight_grid);

        sdsfree(output_dir_with_file);

        output_dir_with_file = sdsnew(output_dir);
        output_dir_with_file = sdscat(output_dir_with_file, "/simulation_result.case");

        if(save_ode_state_variables) {
            persistent_data->n_state_vars = ode_solver->model_data.number_of_ode_equations - 1; // Vm is always saved
        }

        save_case_file(output_dir_with_file, persistent_data->num_files, time_info->dt, print_rate, persistent_data->n_state_vars);

        sdsfree(output_dir_with_file);
        persistent_data->geometry_saved = true;
    }

    sds output_dir_with_file = sdsnew(output_dir);
    output_dir_with_file = sdscat(output_dir_with_file, "/");

    if(persistent_data->n_digits == 0) {
        persistent_data->n_digits = log10(persistent_data->num_files * 500) + 1;
    }

    sds base_name = sdscatprintf(sdsempty(), "Vm.Esca%%0%dd", persistent_data->n_digits);

    char tmp[8192];
    sprintf(tmp, base_name, persistent_data->file_count);

    output_dir_with_file = sdscatprintf(output_dir_with_file, "/%s", tmp);

    save_en6_result_file(output_dir_with_file, the_grid, binary, save_purkinje);

    sdsfree(base_name);
    sdsfree(output_dir_with_file);

    if(persistent_data->n_state_vars) {
        size_t num_sv_entries = ode_solver->model_data.number_of_ode_equations;
        base_name = sdscatprintf(sdsempty(), "Sv%%d.Esca%%0%dd", persistent_data->n_digits);
        real *sv_cpu;

        if(ode_solver->gpu) {

#ifdef COMPILE_GPU
            sv_cpu = MALLOC_ARRAY_OF_TYPE(real, ode_solver->original_num_cells * num_sv_entries);
            // check_cuda_error(cudaMemcpy2D(sv_cpu, ode_solver->original_num_cells * sizeof(real), ode_solver->sv, ode_solver->pitch,
            //                              ode_solver->original_num_cells * sizeof(real), num_sv_entries, cudaMemcpyDeviceToHost));

            memcpy2d_device(sv_cpu, ode_solver->original_num_cells * sizeof(real), ode_solver->sv, ode_solver->pitch,
                            ode_solver->original_num_cells * sizeof(real), num_sv_entries, DEVICE_TO_HOST);
#endif
        } else {
            sv_cpu = ode_solver->sv;
        }

        for(int i = 1; i <= persistent_data->n_state_vars; i++) {

            char tmp[8192];
            sprintf(tmp, base_name, i, persistent_data->file_count);

            sds output_dir_with_file = sdsnew(output_dir);
            output_dir_with_file = sdscat(output_dir_with_file, "/");

            output_dir_with_file = sdscatprintf(output_dir_with_file, "/%s", tmp);

            save_en6_result_file_state_vars(output_dir_with_file, sv_cpu, ode_solver->original_num_cells, num_sv_entries, i, binary, ode_solver->gpu);
            sdsfree(output_dir_with_file);
        }

        sdsfree(base_name);

        if(ode_solver->gpu) {
            free(sv_cpu);
        }
    }

    persistent_data->file_count++;

    CALL_EXTRA_FUNCTIONS(save_mesh_fn, time_info, config, the_grid, ode_solver, purkinje_ode_solver);
}

END_SAVE_MESH(end_save_as_ensight) {
    free(config->persistent_data);
    config->persistent_data = NULL;
}

INIT_SAVE_MESH(init_save_with_activation_times) {

    if(config->persistent_data == NULL) {
        config->persistent_data = calloc(1, sizeof(struct common_persistent_data));

        struct common_persistent_data *cpd = (struct common_persistent_data *)config->persistent_data;

        hmdefault(cpd->cell_was_active, 0.0);
        hmdefault(cpd->last_time_v, -100.0);
        hmdefault(cpd->num_activations, 0);
        hmdefault(cpd->activation_times, NULL);
        hmdefault(cpd->apds, NULL);
        cpd->first_save_call = true;

        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, cpd->print_rate, config, "print_rate");
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, cpd->mesh_print_rate, config, "mesh_print_rate");
    }
}

END_SAVE_MESH(end_save_with_activation_times) {
    free(config->persistent_data);
    config->persistent_data = NULL;
}

SAVE_MESH(save_with_activation_times) {

    int iteration_count = time_info->iteration;
    struct common_persistent_data *cpd = (struct common_persistent_data *)config->persistent_data;

    char *mesh_format = NULL;
    GET_PARAMETER_STRING_VALUE_OR_USE_DEFAULT(mesh_format, config, "mesh_format");

    if(mesh_format && cpd->mesh_print_rate) {
        if(iteration_count % cpd->mesh_print_rate == 0) {
            if(STRINGS_EQUAL("vtk", mesh_format)) {
                save_as_vtk(time_info, config, the_grid, ode_solver, purkinje_ode_solver);
            } else if(STRINGS_EQUAL("vtu", mesh_format)) {
                save_as_vtu(time_info, config, the_grid, ode_solver, purkinje_ode_solver);
            } else if(STRINGS_EQUAL("ensight", mesh_format)) {
                save_as_ensight(time_info, config, the_grid, ode_solver, purkinje_ode_solver);
            } else if(STRINGS_EQUAL("txt", mesh_format)) {
                save_as_text_or_binary(time_info, config, the_grid, ode_solver, purkinje_ode_solver);
            }
        }
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

    struct common_persistent_data *persistent_data = (struct common_persistent_data *)config->persistent_data;

    struct cell_node *grid_cell = the_grid->first_cell;

    real_cpu center_x, center_y, center_z, dx, dy, dz;
    real_cpu v;

    FILE *act_file = fopen(output_dir_with_file, "w");

    fprintf(act_file, "%d\n", (last_t - current_t) <= dt); // rounding errors

    while(grid_cell != 0) {

        if(grid_cell->active || (grid_cell->mesh_extra_info && (FIBROTIC(grid_cell) || BORDER_ZONE(grid_cell)))) {

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
            } else {
                fprintf(act_file, "%g,%g,%g,%g,%g,%g ", center_x, center_y, center_z, dx, dy, dz);
            }

            int n_activations = 0;
            float *apds_array = NULL;
            float *activation_times_array = NULL;

            if(grid_cell->active) {
                float last_v = hmget(persistent_data->last_time_v, cell_coordinates);

                n_activations = (int)hmget(persistent_data->num_activations, cell_coordinates);
                activation_times_array = (float *)hmget(persistent_data->activation_times, cell_coordinates);
                apds_array = (float *)hmget(persistent_data->apds, cell_coordinates);

                int act_times_len = arrlen(activation_times_array);

                if(current_t == 0.0f) {
                    hmput(persistent_data->last_time_v, cell_coordinates, v);
                } else {
                    if((last_v < activation_threshold) && (v >= activation_threshold)) {

                        if(act_times_len == 0) {
                            n_activations++;
                            hmput(persistent_data->num_activations, cell_coordinates, n_activations);
                            arrput(activation_times_array, current_t);
                            float tmp = hmget(persistent_data->cell_was_active, cell_coordinates);
                            hmput(persistent_data->cell_was_active, cell_coordinates, tmp + 1);
                            hmput(persistent_data->activation_times, cell_coordinates, activation_times_array);
                        } else { // This is to avoid spikes in the middle of an Action Potential
                            float last_act_time = activation_times_array[act_times_len - 1];
                            if(current_t - last_act_time > time_threshold) {
                                n_activations++;
                                hmput(persistent_data->num_activations, cell_coordinates, n_activations);
                                arrput(activation_times_array, current_t);
                                float tmp = hmget(persistent_data->cell_was_active, cell_coordinates);
                                hmput(persistent_data->cell_was_active, cell_coordinates, tmp + 1);
                                hmput(persistent_data->activation_times, cell_coordinates, activation_times_array);
                            }
                        }
                    }

                    // CHECK APD
                    bool was_active = (hmget(persistent_data->cell_was_active, cell_coordinates) != 0.0);
                    if(was_active) {
                        if(v <= apd_threshold || (hmget(persistent_data->cell_was_active, cell_coordinates) == 2.0)) {

                            int tmp = (int)hmget(persistent_data->cell_was_active, cell_coordinates);
                            int act_time_array_len = arrlen(activation_times_array);
                            // if this in being calculated because we had a new activation before the cell achieved the rest potential,
                            //  we need to get the activation before this one
                            real_cpu last_act_time = activation_times_array[act_time_array_len - tmp];
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

            for(unsigned long i = 0; i < n_activations; i++) {
                fprintf(act_file, "%lf ", activation_times_array[i]);
            }
            fprintf(act_file, "] ");

            fprintf(act_file, "[ ");

            for(unsigned long i = 0; i < arrlen(apds_array); i++) {
                fprintf(act_file, "%lf ", apds_array[i]);
            }
            fprintf(act_file, "]\n");
        }

        grid_cell = grid_cell->next;
    }

    fclose(act_file);

    CALL_EXTRA_FUNCTIONS(save_mesh_fn, time_info, config, the_grid, ode_solver, purkinje_ode_solver);
}

SAVE_MESH(save_vm_matrix) {

    static FILE *f = NULL;

    real save_after = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, save_after, config, "start_saving_at");

    if(time_info->current_t >= save_after) {

        if(f == NULL) {
            GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config, "output_dir");
            sds output_dir_with_file = sdsnew(output_dir);
            output_dir_with_file = sdscat(output_dir_with_file, "/Vm_matrix.txt");
            f = fopen(output_dir_with_file, "w");
            sdsfree(output_dir_with_file);
        }

        uint32_t n_active = the_grid->num_active_cells;

        fprintf(f, "%e ", time_info->current_t);

        for(uint32_t i = 0; i < n_active; i++) {
            fprintf(f, "%e ", the_grid->active_cells[i]->v);
        }
    }

    fprintf(f, "\n");

    if(time_info->current_t >= time_info->final_t) {
        fclose(f);
    }
}

SAVE_MESH(no_save) {
    // Nop
}

INIT_SAVE_MESH(init_save_multiple_cell_state_variables) {
    config->persistent_data = malloc(sizeof(struct save_multiple_cell_state_variables_persistent_data));
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(uint32_t, ((struct save_multiple_cell_state_variables_persistent_data *)config->persistent_data)->num_cells,
                                                config, "num_cells");
    GET_PARAMETER_MATRIX_VALUE_OR_USE_DEFAULT(((struct save_multiple_cell_state_variables_persistent_data *)config->persistent_data)->cell_centers, config,
                                              "cell_centers", ((struct save_multiple_cell_state_variables_persistent_data *)config->persistent_data)->num_cells,
                                              3);
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(((struct save_multiple_cell_state_variables_persistent_data *)config->persistent_data)->file_name_prefix, config,
                                               "file_name_prefix");

    uint32_t num_cells = ((struct save_multiple_cell_state_variables_persistent_data *)config->persistent_data)->num_cells;
    char *file_name_prefix = ((struct save_multiple_cell_state_variables_persistent_data *)config->persistent_data)->file_name_prefix;

    ((struct save_multiple_cell_state_variables_persistent_data *)config->persistent_data)->files = MALLOC_ARRAY_OF_TYPE(FILE *, num_cells);
    ((struct save_multiple_cell_state_variables_persistent_data *)config->persistent_data)->cell_sv_positions = MALLOC_ARRAY_OF_TYPE(uint32_t, num_cells);

    for(int i = 0; i < num_cells; i++) {

        sds base_name = NULL;
        base_name = create_base_name(file_name_prefix, i, "dat");

        ((struct save_multiple_cell_state_variables_persistent_data *)config->persistent_data)->files[i] = fopen(base_name, "w");
        ((struct save_multiple_cell_state_variables_persistent_data *)config->persistent_data)->cell_sv_positions[i] = -1;
    }
}

SAVE_MESH(save_multiple_cell_state_variables) {
    struct save_multiple_cell_state_variables_persistent_data *params = ((struct save_multiple_cell_state_variables_persistent_data *)config->persistent_data);

    for(uint32_t i = 0; i < params->num_cells; i++) {
        if(params->cell_sv_positions[i] == -1) {
            if(!the_grid->adaptive) {
                FOR_EACH_CELL(the_grid) {
                    if(cell->center.x == params->cell_centers[i * 3] && cell->center.y == params->cell_centers[i * 3 + 1] &&
                       cell->center.z == params->cell_centers[i * 3 + 2]) {
                        params->cell_sv_positions[i] = cell->sv_position;
                        break;
                    }
                }
            }
        }
    }

    if(ode_solver->gpu) {
#ifdef COMPILE_GPU
        for(uint32_t k = 0; k < params->num_cells; k++) {
            real *cell_sv;

            cell_sv = MALLOC_ARRAY_OF_TYPE(real, ode_solver->model_data.number_of_ode_equations);

            // check_cuda_error(cudaMemcpy2D(cell_sv, sizeof(real), ode_solver->sv + params->cell_sv_positions[k], ode_solver->pitch, sizeof(real),
            //                                        ode_solver->model_data.number_of_ode_equations, cudaMemcpyDeviceToHost));

            memcpy2d_device(cell_sv, sizeof(real), ode_solver->sv + params->cell_sv_positions[k], ode_solver->pitch, sizeof(real),
                            ode_solver->model_data.number_of_ode_equations, DEVICE_TO_HOST);

            fprintf(params->files[k], "%lf ", time_info->current_t);
            for(int i = 0; i < ode_solver->model_data.number_of_ode_equations; i++) {
                fprintf(params->files[k], "%lf ", cell_sv[i]);
            }
            fprintf(params->files[k], "\n");

            free(cell_sv);
        }

#endif
    } else {
        for(uint32_t k = 0; k < params->num_cells; k++) {
            real *cell_sv = &ode_solver->sv[params->cell_sv_positions[k] * ode_solver->model_data.number_of_ode_equations];

            fprintf(params->files[k], "%lf ", time_info->current_t);
            for(int i = 0; i < ode_solver->model_data.number_of_ode_equations; i++) {
                fprintf(params->files[k], "%lf ", cell_sv[i]);
            }
            fprintf(params->files[k], "\n");
        }
    }
}

END_SAVE_MESH(end_save_multiple_cell_state_variables) {
    uint32_t num_cells = ((struct save_multiple_cell_state_variables_persistent_data *)config->persistent_data)->num_cells;
    free(((struct save_multiple_cell_state_variables_persistent_data *)config->persistent_data)->file_name_prefix);
    free(((struct save_multiple_cell_state_variables_persistent_data *)config->persistent_data)->cell_sv_positions);
    free(((struct save_multiple_cell_state_variables_persistent_data *)config->persistent_data)->cell_centers);
    for(uint32_t i = 0; i < num_cells; i++) {
        fclose(((struct save_multiple_cell_state_variables_persistent_data *)config->persistent_data)->files[i]);
    }
    free(((struct save_multiple_cell_state_variables_persistent_data *)config->persistent_data)->files);
    free(config->persistent_data);
    config->persistent_data = NULL;
}

SAVE_MESH(save_transmurality_as_vtk) {

    struct extra_data_for_torord *extra_data = (struct extra_data_for_torord *)ode_solver->ode_extra_data;
    real *transmurality = extra_data->transmurality;

    int iteration_count = time_info->iteration;

    if(((struct common_persistent_data *)config->persistent_data)->first_save_call) {
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config, "output_dir");
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(file_prefix, config, "file_prefix");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_plain, config, "clip_with_plain");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_bounds, config, "clip_with_bounds");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(binary, config, "binary");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_f, config, "save_f");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_visible_mask, config, "save_visible_mask");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_scar_cells, config, "save_scar_cells");

        ((struct common_persistent_data *)config->persistent_data)->first_save_call = false;
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

    // TODO: change this. We dont need the current_t here
    output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, current_t);

    bool read_only_data = ((struct common_persistent_data *)config->persistent_data)->grid != NULL;

    new_vtk_unstructured_grid_from_alg_grid(&(((struct common_persistent_data *)config->persistent_data)->grid), the_grid, clip_with_plain, plain_coords,
                                            clip_with_bounds, bounds, read_only_data, save_f, save_scar_cells, transmurality);

    save_vtk_unstructured_grid_as_legacy_vtk(((struct common_persistent_data *)config->persistent_data)->grid, output_dir_with_file, binary, save_f, NULL);

    if(save_visible_mask) {
        save_visibility_mask(output_dir_with_file, (((struct common_persistent_data *)config->persistent_data)->grid)->cell_visibility);
    }

    if(the_grid->adaptive) {
        free_vtk_unstructured_grid(((struct common_persistent_data *)config->persistent_data)->grid);
        ((struct common_persistent_data *)config->persistent_data)->grid = NULL;
    }

    sdsfree(output_dir_with_file);
    sdsfree(base_name);

    CALL_EXTRA_FUNCTIONS(save_mesh_fn, time_info, config, the_grid, ode_solver, purkinje_ode_solver);

    exit(EXIT_SUCCESS);
}
