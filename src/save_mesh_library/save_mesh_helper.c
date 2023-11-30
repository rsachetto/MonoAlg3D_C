#include "save_mesh_helper.h"
#include "../domains_library/mesh_info_data.h"

static void write_pvd_header(FILE *pvd_file) {
    fprintf(pvd_file, "<VTKFile type=\"Collection\" version=\"0.1\" compressor=\"vtkZLibDataCompressor\">\n");
    fprintf(pvd_file, "\t<Collection>\n");
    fprintf(pvd_file, "\t</Collection>\n");
    fprintf(pvd_file, "</VTKFile>");
}

sds create_base_name(char *f_prefix, int iteration_count, char *extension) {
    return sdscatprintf(sdsempty(), "%s_it_%d.%s", f_prefix, iteration_count, extension);
}

void add_file_to_pvd(real_cpu current_t, const char *output_dir, const char *base_name, bool first_call) {

    sds pvd_name = sdsnew(output_dir);
    pvd_name = sdscat(pvd_name, "/simulation_result.pvd");

    static FILE *pvd_file = NULL;
    pvd_file = fopen(pvd_name, "r+");

    if(!pvd_file) {
        pvd_file = fopen(pvd_name, "w");
        write_pvd_header(pvd_file);
    } else {
        if(first_call) {
            fclose(pvd_file);
            pvd_file = fopen(pvd_name, "w");
            write_pvd_header(pvd_file);
        }
    }

    sdsfree(pvd_name);

    fseek(pvd_file, -26, SEEK_END);

    fprintf(pvd_file, "\n\t\t<DataSet timestep=\"%lf\" group=\"\" part=\"0\" file=\"%s\"/>\n", current_t, base_name);
    fprintf(pvd_file, "\t</Collection>\n");
    fprintf(pvd_file, "</VTKFile>");
    fclose(pvd_file);
}

void calculate_purkinje_activation_time_and_apd(struct time_info *time_info, struct config *config, struct grid *the_grid, const real_cpu time_threshold,
                                                const real_cpu purkinje_activation_threshold, const real_cpu purkinje_apd_threshold) {

    real_cpu current_t = time_info->current_t;
    real_cpu last_t = time_info->final_t;
    real_cpu dt = time_info->dt;

    struct save_coupling_with_activation_times_persistent_data *persistent_data =
        (struct save_coupling_with_activation_times_persistent_data *)config->persistent_data;

    struct cell_node *purkinje_grid_cell = the_grid->purkinje->first_cell;

    real_cpu center_x, center_y, center_z;
    real_cpu v;

    while(purkinje_grid_cell != 0) {

        if(purkinje_grid_cell->active || (purkinje_grid_cell->mesh_extra_info && (FIBROTIC(purkinje_grid_cell) || BORDER_ZONE(purkinje_grid_cell)))) {

            center_x = purkinje_grid_cell->center.x;
            center_y = purkinje_grid_cell->center.y;
            center_z = purkinje_grid_cell->center.z;

            v = purkinje_grid_cell->v;

            struct point_3d cell_coordinates;
            cell_coordinates.x = center_x;
            cell_coordinates.y = center_y;
            cell_coordinates.z = center_z;

            int n_activations = 0;
            float *apds_array = NULL;
            float *activation_times_array = NULL;

            if(purkinje_grid_cell->active) {

                float last_v = hmget(persistent_data->purkinje_last_time_v, cell_coordinates);

                n_activations = (int)hmget(persistent_data->purkinje_num_activations, cell_coordinates);
                activation_times_array = (float *)hmget(persistent_data->purkinje_activation_times, cell_coordinates);
                apds_array = (float *)hmget(persistent_data->purkinje_apds, cell_coordinates);

                int act_times_len = arrlen(activation_times_array);

                if(current_t == 0.0f) {
                    hmput(persistent_data->purkinje_last_time_v, cell_coordinates, v);
                } else {
                    if((last_v < purkinje_activation_threshold) && (v >= purkinje_activation_threshold)) {

                        if(act_times_len == 0) {
                            n_activations++;
                            hmput(persistent_data->purkinje_num_activations, cell_coordinates, n_activations);
                            arrput(activation_times_array, current_t);
                            float tmp = hmget(persistent_data->purkinje_cell_was_active, cell_coordinates);
                            hmput(persistent_data->purkinje_cell_was_active, cell_coordinates, tmp + 1);
                            hmput(persistent_data->purkinje_activation_times, cell_coordinates, activation_times_array);
                        } else { // This is to avoid spikes in the middle of an Action Potential
                            float last_act_time = activation_times_array[act_times_len - 1];
                            if(current_t - last_act_time > time_threshold) {
                                n_activations++;
                                hmput(persistent_data->purkinje_num_activations, cell_coordinates, n_activations);
                                arrput(activation_times_array, current_t);
                                float tmp = hmget(persistent_data->purkinje_cell_was_active, cell_coordinates);
                                hmput(persistent_data->purkinje_cell_was_active, cell_coordinates, tmp + 1);
                                hmput(persistent_data->purkinje_activation_times, cell_coordinates, activation_times_array);
                            }
                        }
                    }

                    // CHECK APD
                    bool was_active = (hmget(persistent_data->purkinje_cell_was_active, cell_coordinates) != 0.0);
                    if(was_active) {
                        if(v <= purkinje_apd_threshold || (hmget(persistent_data->purkinje_cell_was_active, cell_coordinates) == 2.0) ||
                           (last_t - current_t) <= dt) {

                            int tmp = (int)hmget(persistent_data->purkinje_cell_was_active, cell_coordinates);
                            int act_time_array_len = arrlen(activation_times_array);
                            // if this in being calculated because we had a new activation before the cell achieved the rest potential,
                            // we need to get the activation before this one
                            real_cpu last_act_time = activation_times_array[act_time_array_len - tmp];
                            real_cpu apd = current_t - last_act_time;
                            arrput(apds_array, apd);
                            hmput(persistent_data->purkinje_apds, cell_coordinates, apds_array);
                            hmput(persistent_data->purkinje_cell_was_active, cell_coordinates, tmp - 1);
                        }
                    }

                    hmput(persistent_data->purkinje_last_time_v, cell_coordinates, v);
                }
            }
        }
        purkinje_grid_cell = purkinje_grid_cell->next;
    }
}

void write_purkinje_activation_time_maps(struct config *config, struct grid *the_grid, char *output_dir, char *file_prefix, bool clip_with_plain,
                                         bool clip_with_bounds, bool binary, bool compress, int compression_level) {

    if(((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->first_save_call) {
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config, "output_dir");
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(file_prefix, config, "file_prefix");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_plain, config, "clip_with_plain");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_bounds, config, "clip_with_bounds");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(binary, config, "binary");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(compress, config, "compress");
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, compression_level, config, "compression_level");

        if(compress)
            binary = true;
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

    write_purkinje_activation_time_for_each_pulse(config, the_grid, output_dir, plain_coords, bounds, clip_with_plain, clip_with_bounds, binary, compress,
                                                  compression_level);
}

void write_purkinje_activation_time_for_each_pulse(struct config *config, struct grid *the_grid, char *output_dir, float plain_coords[], float bounds[],
                                                   bool clip_with_plain, bool clip_with_bounds, bool binary, bool compress, int compression_level) {

    struct save_coupling_with_activation_times_persistent_data **data = (struct save_coupling_with_activation_times_persistent_data **)&config->persistent_data;

    struct cell_node **grid_cell = the_grid->purkinje->purkinje_cells;

    real_cpu center_x, center_y, center_z;

    center_x = grid_cell[0]->center.x;
    center_y = grid_cell[0]->center.y;
    center_z = grid_cell[0]->center.z;

    struct point_3d cell_coordinates;
    cell_coordinates.x = center_x;
    cell_coordinates.y = center_y;
    cell_coordinates.z = center_z;

    // Get the number of pulses using one cell of the grid
    int n_pulses = (int)hmget((*data)->purkinje_num_activations, cell_coordinates);

    // Write the activation time map for each pulse
    for(int cur_pulse = 0; cur_pulse < n_pulses; cur_pulse++) {

        sds base_name = create_base_name("purkinje_activation_time_map_pulse", cur_pulse, "vtp");
        sds output_dir_with_file = sdsnew(output_dir);
        output_dir_with_file = sdscat(output_dir_with_file, "/");
        output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, cur_pulse);

        bool read_only_data = ((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->purkinje_grid != NULL;
        new_vtk_polydata_grid_from_purkinje_grid(&((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->purkinje_grid,
                                                 the_grid->purkinje, clip_with_plain, plain_coords, clip_with_bounds, bounds, read_only_data);

        set_purkinje_vtk_values_with_activation_time_from_current_pulse(&config->persistent_data, the_grid, cur_pulse);

        if(compress) {
            save_vtk_polydata_grid_as_vtp_compressed(((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->purkinje_grid,
                                                     output_dir_with_file, compression_level);
        } else {
            save_vtk_polydata_grid_as_vtp(((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->purkinje_grid,
                                          output_dir_with_file, binary);
        }

        free_vtk_polydata_grid(((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->purkinje_grid);
        ((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->purkinje_grid = NULL;

        sdsfree(output_dir_with_file);
        sdsfree(base_name);
    }
}

void set_purkinje_vtk_values_with_activation_time_from_current_pulse(void **persistent_data, struct grid *the_grid, const int cur_pulse) {

    struct save_coupling_with_activation_times_persistent_data **data = (struct save_coupling_with_activation_times_persistent_data **)persistent_data;

    struct cell_node **purkinje_grid_cell = the_grid->purkinje->purkinje_cells;
    uint32_t num_active_cells = the_grid->purkinje->num_active_purkinje_cells;

    real_cpu center_x, center_y, center_z;

    // We need to pass through the cells in the same order as we did when the VTU grid was created
    for(uint32_t i = 0; i < num_active_cells; i++) {

        center_x = purkinje_grid_cell[i]->center.x;
        center_y = purkinje_grid_cell[i]->center.y;
        center_z = purkinje_grid_cell[i]->center.z;

        struct point_3d cell_coordinates;
        cell_coordinates.x = center_x;
        cell_coordinates.y = center_y;
        cell_coordinates.z = center_z;

        float *activation_times_array = NULL;

        if(purkinje_grid_cell[i]->active) {

            activation_times_array = (float *)hmget((*data)->purkinje_activation_times, cell_coordinates);

            if(activation_times_array) {
            // Get the activation time from the current pulse for that particular cell
            float at = activation_times_array[cur_pulse];

            // Update the scalar value from the "vtk_unstructured_grid"
            (*data)->purkinje_grid->values[i] = at;
            }
        }
    }
}

void write_purkinje_apd_map(struct config *config, struct grid *the_grid, char *output_dir, char *file_prefix, bool clip_with_plain, bool clip_with_bounds,
                            bool binary, bool compress, int compression_level) {

    real_cpu current_t = 0.0;

    if(((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->first_save_call) {
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config, "output_dir");
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(file_prefix, config, "file_prefix");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_plain, config, "clip_with_plain");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_bounds, config, "clip_with_bounds");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(binary, config, "binary");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(compress, config, "compress");
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, compression_level, config, "compression_level");

        if(compress)
            binary = true;
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

    sds base_name = create_base_name("purkinje_apd_map", 0, "vtp");
    sds output_dir_with_file = sdsnew(output_dir);
    output_dir_with_file = sdscat(output_dir_with_file, "/");
    output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, current_t);

    bool read_only_data = ((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->purkinje_grid != NULL;
    new_vtk_polydata_grid_from_purkinje_grid(&((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->purkinje_grid,
                                             the_grid->purkinje, clip_with_plain, plain_coords, clip_with_bounds, bounds, read_only_data);

    set_purkinje_vtk_values_with_mean_apd(&config->persistent_data, the_grid);

    if(compress) {
        save_vtk_polydata_grid_as_vtp_compressed(((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->purkinje_grid,
                                                 output_dir_with_file, compression_level);
    } else {
        save_vtk_polydata_grid_as_vtp(((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->purkinje_grid,
                                      output_dir_with_file, binary);
    }

    free_vtk_polydata_grid(((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->purkinje_grid);
    ((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->purkinje_grid = NULL;

    sdsfree(output_dir_with_file);
    sdsfree(base_name);
}

void set_purkinje_vtk_values_with_mean_apd(void **persistent_data, struct grid *the_grid) {
    struct save_coupling_with_activation_times_persistent_data **data = (struct save_coupling_with_activation_times_persistent_data **)persistent_data;

    struct cell_node **purkinje_grid_cell = the_grid->purkinje->purkinje_cells;
    uint32_t num_active_cells = the_grid->purkinje->num_active_purkinje_cells;

    real_cpu center_x, center_y, center_z;

    // We need to pass through the cells in the same order as we did when the VTU grid was created
    for(uint32_t i = 0; i < num_active_cells; i++) {

        center_x = purkinje_grid_cell[i]->center.x;
        center_y = purkinje_grid_cell[i]->center.y;
        center_z = purkinje_grid_cell[i]->center.z;

        struct point_3d cell_coordinates;
        cell_coordinates.x = center_x;
        cell_coordinates.y = center_y;
        cell_coordinates.z = center_z;

        float *apds_array = NULL;

        if(purkinje_grid_cell[i]->active) {

            apds_array = (float *)hmget((*data)->purkinje_apds, cell_coordinates);

            uint64_t apd_len = arrlen(apds_array);

            // Update the scalar value from the "vtk_unstructured_grid"
            (*data)->purkinje_grid->values[i] = calculate_mean(apds_array, apd_len);
        }
    }
}

void calculate_tissue_activation_time_and_apd(struct time_info *time_info, struct config *config, struct grid *the_grid, const real_cpu time_threshold,
                                              const real_cpu tissue_activation_threshold, const real_cpu tissue_apd_threshold) {

    real_cpu current_t = time_info->current_t;
    real_cpu last_t = time_info->final_t;
    real_cpu dt = time_info->dt;

    struct save_coupling_with_activation_times_persistent_data *persistent_data =
        (struct save_coupling_with_activation_times_persistent_data *)config->persistent_data;

    struct cell_node *tissue_grid_cell = the_grid->first_cell;

    real_cpu center_x, center_y, center_z;
    real_cpu v;

    while(tissue_grid_cell != 0) {

        if(tissue_grid_cell->active || (tissue_grid_cell->mesh_extra_info && (FIBROTIC(tissue_grid_cell) || BORDER_ZONE(tissue_grid_cell)))) {

            center_x = tissue_grid_cell->center.x;
            center_y = tissue_grid_cell->center.y;
            center_z = tissue_grid_cell->center.z;

            v = tissue_grid_cell->v;

            struct point_3d cell_coordinates;
            cell_coordinates.x = center_x;
            cell_coordinates.y = center_y;
            cell_coordinates.z = center_z;

            int n_activations = 0;
            float *apds_array = NULL;
            float *activation_times_array = NULL;

            if(tissue_grid_cell->active) {

                float last_v = hmget(persistent_data->tissue_last_time_v, cell_coordinates);

                n_activations = (int)hmget(persistent_data->tissue_num_activations, cell_coordinates);
                activation_times_array = (float *)hmget(persistent_data->tissue_activation_times, cell_coordinates);
                apds_array = (float *)hmget(persistent_data->tissue_apds, cell_coordinates);

                int act_times_len = arrlen(activation_times_array);

                if(current_t == 0.0f) {
                    hmput(persistent_data->tissue_last_time_v, cell_coordinates, v);
                } else {
                    if((last_v < tissue_activation_threshold) && (v >= tissue_activation_threshold)) {
                        if(act_times_len == 0) {
                            n_activations++;
                            hmput(persistent_data->tissue_num_activations, cell_coordinates, n_activations);
                            arrput(activation_times_array, current_t);
                            float tmp = hmget(persistent_data->tissue_cell_was_active, cell_coordinates);
                            hmput(persistent_data->tissue_cell_was_active, cell_coordinates, tmp + 1);
                            hmput(persistent_data->tissue_activation_times, cell_coordinates, activation_times_array);
                        } else { // This is to avoid spikes in the middle of an Action Potential
                            float last_act_time = activation_times_array[act_times_len - 1];
                            if(current_t - last_act_time > time_threshold) {
                                n_activations++;
                                hmput(persistent_data->tissue_num_activations, cell_coordinates, n_activations);
                                arrput(activation_times_array, current_t);
                                float tmp = hmget(persistent_data->tissue_cell_was_active, cell_coordinates);
                                hmput(persistent_data->tissue_cell_was_active, cell_coordinates, tmp + 1);
                                hmput(persistent_data->tissue_activation_times, cell_coordinates, activation_times_array);
                            }
                        }
                    }

                    // CHECK APD
                    bool was_active = (hmget(persistent_data->tissue_cell_was_active, cell_coordinates) != 0.0);
                    if(was_active) {
                        if(v <= tissue_apd_threshold || (hmget(persistent_data->tissue_cell_was_active, cell_coordinates) == 2.0) ||
                           (last_t - current_t) <= dt) {

                            int tmp = (int)hmget(persistent_data->tissue_cell_was_active, cell_coordinates);
                            int act_time_array_len = arrlen(activation_times_array);
                            // if this in being calculated because we had a new activation before the cell achieved the rest potential,
                            // we need to get the activation before this one
                            real_cpu last_act_time = activation_times_array[act_time_array_len - tmp];
                            real_cpu apd = current_t - last_act_time;
                            arrput(apds_array, apd);
                            hmput(persistent_data->tissue_apds, cell_coordinates, apds_array);
                            hmput(persistent_data->tissue_cell_was_active, cell_coordinates, tmp - 1);
                        }
                    }

                    hmput(persistent_data->tissue_last_time_v, cell_coordinates, v);
                }
            }
        }
        tissue_grid_cell = tissue_grid_cell->next;
    }
}

void write_tissue_apd_map(struct config *config, struct grid *the_grid, char *output_dir, char *file_prefix, bool clip_with_plain, bool clip_with_bounds,
                          bool binary, bool compress, int compression_level, bool save_f) {

    real_cpu current_t = 0.0;

    if(((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->first_save_call) {
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config, "output_dir");
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(file_prefix, config, "file_prefix");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_plain, config, "clip_with_plain");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_bounds, config, "clip_with_bounds");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(binary, config, "binary");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(compress, config, "compress");
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, compression_level, config, "compression_level");

        if(compress)
            binary = true;

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

    sds base_name = create_base_name("tissue_apd_map", 0, "vtu");
    sds output_dir_with_file = sdsnew(output_dir);
    output_dir_with_file = sdscat(output_dir_with_file, "/");
    output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, current_t);

    bool read_only_data = ((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->tissue_grid != NULL;
    new_vtk_unstructured_grid_from_alg_grid(&((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->tissue_grid, the_grid,
                                            clip_with_plain, plain_coords, clip_with_bounds, bounds, read_only_data, save_f, false, NULL);

    set_tissue_vtk_values_with_mean_apd(&config->persistent_data, the_grid);

    if(compress) {
        save_vtk_unstructured_grid_as_vtu_compressed(((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->tissue_grid,
                                                     output_dir_with_file, compression_level);
    } else {
        save_vtk_unstructured_grid_as_vtu(((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->tissue_grid,
                                          output_dir_with_file, binary);
    }

    free_vtk_unstructured_grid(((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->tissue_grid);
    ((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->tissue_grid = NULL;

    sdsfree(output_dir_with_file);
    sdsfree(base_name);
}

void write_tissue_activation_time_maps(struct config *config, struct grid *the_grid, char *output_dir, char *file_prefix, bool clip_with_plain,
                                       bool clip_with_bounds, bool binary, bool compress, int compression_level, bool save_f) {

    if(((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->first_save_call) {
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config, "output_dir");
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(file_prefix, config, "file_prefix");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_plain, config, "clip_with_plain");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_bounds, config, "clip_with_bounds");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(binary, config, "binary");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(compress, config, "compress");
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, compression_level, config, "compression_level");

        if(compress)
            binary = true;

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

    write_tissue_activation_time_for_each_pulse(config, the_grid, output_dir, plain_coords, bounds, clip_with_plain, clip_with_bounds, binary, compress,
                                                compression_level, save_f);
}

void write_tissue_activation_time_for_each_pulse(struct config *config, struct grid *the_grid, char *output_dir, float plain_coords[], float bounds[],
                                                 bool clip_with_plain, bool clip_with_bounds, bool binary, bool compress, int compression_level, bool save_f) {

    struct save_coupling_with_activation_times_persistent_data **data = (struct save_coupling_with_activation_times_persistent_data **)&config->persistent_data;

    struct cell_node **grid_cell = the_grid->active_cells;

    real_cpu center_x, center_y, center_z;

    center_x = grid_cell[0]->center.x;
    center_y = grid_cell[0]->center.y;
    center_z = grid_cell[0]->center.z;

    struct point_3d cell_coordinates;
    cell_coordinates.x = center_x;
    cell_coordinates.y = center_y;
    cell_coordinates.z = center_z;

    // Get the number of pulses using one cell of the grid
    int n_pulses = (int)hmget((*data)->tissue_num_activations, cell_coordinates);

    // Write the activation time map for each pulse
    for(int cur_pulse = 0; cur_pulse < n_pulses; cur_pulse++) {

        sds base_name = create_base_name("tissue_activation_time_map_pulse", cur_pulse, "vtu");
        sds output_dir_with_file = sdsnew(output_dir);
        output_dir_with_file = sdscat(output_dir_with_file, "/");
        output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, cur_pulse);

        bool read_only_data = ((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->tissue_grid != NULL;
        new_vtk_unstructured_grid_from_alg_grid(&((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->tissue_grid, the_grid,
                                                clip_with_plain, plain_coords, clip_with_bounds, bounds, read_only_data, save_f, false, NULL);

        set_tissue_vtk_values_with_activation_time_from_current_pulse(&config->persistent_data, the_grid, cur_pulse);

        if(compress) {
            save_vtk_unstructured_grid_as_vtu_compressed(((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->tissue_grid,
                                                         output_dir_with_file, compression_level);
        } else {
            save_vtk_unstructured_grid_as_vtu(((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->tissue_grid,
                                              output_dir_with_file, binary);
        }

        free_vtk_unstructured_grid(((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->tissue_grid);
        ((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->tissue_grid = NULL;

        sdsfree(output_dir_with_file);
        sdsfree(base_name);
    }
}

void set_tissue_vtk_values_with_activation_time_from_current_pulse(void **persistent_data, struct grid *the_grid, const int cur_pulse) {

    struct save_coupling_with_activation_times_persistent_data **data = (struct save_coupling_with_activation_times_persistent_data **)persistent_data;

    struct cell_node **grid_cell = the_grid->active_cells;
    uint32_t num_active_cells = the_grid->num_active_cells;

    real_cpu center_x, center_y, center_z;

    // We need to pass through the cells in the same order as we did when the VTU grid was created
    for(uint32_t i = 0; i < num_active_cells; i++) {

        center_x = grid_cell[i]->center.x;
        center_y = grid_cell[i]->center.y;
        center_z = grid_cell[i]->center.z;

        struct point_3d cell_coordinates;
        cell_coordinates.x = center_x;
        cell_coordinates.y = center_y;
        cell_coordinates.z = center_z;

        float *activation_times_array = NULL;

        if(grid_cell[i]->active) {

            activation_times_array = (float *)hmget((*data)->tissue_activation_times, cell_coordinates);

            // TODO: Check the case where: (activation_times_array == NULL)
            // Get the activation time from the current pulse for that particular cell
            float at;
            if(activation_times_array == NULL) {
                at = -1;
            } else {
                at = activation_times_array[cur_pulse];
            }

            // Update the scalar value from the "vtk_unstructured_grid"
            (*data)->tissue_grid->values[i] = at;
        }
    }
}

void set_tissue_vtk_values_with_mean_apd(void **persistent_data, struct grid *the_grid) {

    struct save_coupling_with_activation_times_persistent_data **data = (struct save_coupling_with_activation_times_persistent_data **)persistent_data;

    struct cell_node **grid_cell = the_grid->active_cells;
    uint32_t num_active_cells = the_grid->num_active_cells;

    real_cpu center_x, center_y, center_z;

    // We need to pass through the cells in the same order as we did when the VTU grid was created
    for(uint32_t i = 0; i < num_active_cells; i++) {

        center_x = grid_cell[i]->center.x;
        center_y = grid_cell[i]->center.y;
        center_z = grid_cell[i]->center.z;

        struct point_3d cell_coordinates;
        cell_coordinates.x = center_x;
        cell_coordinates.y = center_y;
        cell_coordinates.z = center_z;

        float *apds_array = NULL;

        if(grid_cell[i]->active) {
            apds_array = (float *)hmget((*data)->tissue_apds, cell_coordinates);
            uint64_t apd_len = arrlen(apds_array);

            // Update the scalar value from the "vtk_unstructured_grid"
            (*data)->tissue_grid->values[i] = calculate_mean(apds_array, apd_len);
        }
    }
}

// TODO: Think about how we are going to calculate the propagation velocity in the retrograde direction
void print_purkinje_propagation_velocity(struct config *config, struct grid *the_grid) {

    assert(config);
    assert(the_grid);

    struct save_coupling_with_activation_times_persistent_data *persistent_data =
        (struct save_coupling_with_activation_times_persistent_data *)config->persistent_data;
    struct cell_node **purkinje_cells = the_grid->purkinje->purkinje_cells;

    // Compute the shortest distance from the root to all other nodes
    struct graph *the_network = the_grid->purkinje->network;
    double *shortest_dist = dijkstra(the_network, 0);

    // Calculate the propagation velocity of each terminal in the Purkinje network
    struct node *n = the_network->list_nodes;
    while(n != NULL) {
        if(is_terminal(n)) {
            struct point_3d cell_coordinates;
            real_cpu center_x, center_y, center_z;

            uint32_t id = n->id;
            center_x = purkinje_cells[id]->center.x;
            center_y = purkinje_cells[id]->center.y;
            center_z = purkinje_cells[id]->center.z;

            cell_coordinates.x = center_x;
            cell_coordinates.y = center_y;
            cell_coordinates.z = center_z;

            // Get the total number of pulses
            size_t n_pulses = hmget(persistent_data->purkinje_num_activations, cell_coordinates);

            // Get the terminal cell LAT
            float *activation_times_array = (float *)hmget(persistent_data->purkinje_activation_times, cell_coordinates);

            // Calculate the stimulation period
            float period = 0.0;
            if (n_pulses > 1) {
                period = activation_times_array[1] - activation_times_array[0];
            }
            for(size_t j = 0; j < n_pulses; j++) {
                real_cpu delta_t = activation_times_array[j] - j*period;
                real_cpu delta_s = shortest_dist[id];
                real_cpu v = delta_s / delta_t;

                log_info("[Pulse %u] Purkinje cell %u || Delta_s = %g um || Delta_t = %g ms || v = %g um/mm\n", j + 1, id, delta_s, delta_t, v);
            }
        }
        n = n->next;
    }

    free(shortest_dist);
}
