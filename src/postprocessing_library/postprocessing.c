#include "../3dparty/sds/sds.h"
#include "../config/postprocessor_config.h"
#include "../config_helpers/config_helpers.h"
#include "../utils/file_utils.h"
#include "../vtk_utils/vtk_unstructured_grid.h"

#include <stdio.h>

POSTPROCESS(calculate_all_cells_activation_time) {

    char *input_dir = NULL;
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(input_dir, config, "input_dir");

    char *output_dir = NULL;
    GET_PARAMETER_STRING_VALUE_OR_USE_DEFAULT(output_dir, config, "output_dir");

    char *input_type = NULL;
    GET_PARAMETER_STRING_VALUE_OR_USE_DEFAULT(input_type, config, "input_type");

    float time_threshold = 10.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, time_threshold, config, "time_threshold");

    float activation_threshold = -30.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, activation_threshold, config, "activation_threshold");

    float apd_threshold = -83.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, apd_threshold, config, "apd_threshold");

    float dt = -1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, dt, config, "dt");

    if(input_type == NULL) {
        input_type = strdup("text");
    }

    if(strncmp(input_type, "txt", 4) == 0 || strncmp(input_type, "vtk", 3) == 0 || strncmp(input_type, "vtu", 3) == 0 || strncmp(input_type, "bin", 3) == 0) {
        printf("[INFO] Postprocessing %s files\n", input_type);
    } else {
        fprintf(stderr, "input_type = %s is not supported. Valid types are text or vtk!\n", input_type);
        exit(EXIT_FAILURE);
    }

    struct path_information input_info;

    get_path_information(input_dir, &input_info);

    if(!input_info.exists) {
        fprintf(stderr, "[ERR] Invalid directory (%s)! The input_dir parameter should be and valid directory containing the input files!\n", input_dir);
        exit(EXIT_FAILURE);
    }

    if(!output_dir) {
        output_dir = sdsnew(input_dir);
    }

    create_dir(output_dir);

    string_array files_list = list_files_from_dir(input_dir, NULL, input_type, NULL, true);

    int num_files = arrlen(files_list);

    if(num_files == 0) {
        fprintf(stderr, "[WARN] No files with extension %s found in directory %s Exiting!\n", input_type, input_dir);
        exit(EXIT_FAILURE);
    }

    struct vtk_unstructured_grid *vtk_grid;
    struct point_hash_entry *last_time_v = NULL;
    struct point_hash_entry *num_activations = NULL;
    struct point_hash_entry *cell_was_active = NULL;
    struct point_voidp_hash_entry *activation_times = NULL;
    struct point_voidp_hash_entry *apds = NULL;

    hmdefault(cell_was_active, 0.0);
    hmdefault(last_time_v, -100.0);
    hmdefault(num_activations, 0);
    hmdefault(activation_times, NULL);
    hmdefault(apds, NULL);

    int step;
    int step1;
    int step2 = 0;
    int final_step;
    float final_time;
    float current_t = 0.0;

    step1 = get_step_from_filename(files_list[0]);

    if(num_files > 1) {
        step2 = get_step_from_filename(files_list[1]);
        step = step2 - step1;
    } else {
        step = step1;
    }

    final_step = get_step_from_filename(files_list[num_files - 1]);

    if(dt == -1) {
        final_time = final_step;
        dt = step;

    } else {
        final_time = final_step * dt;
    }

    ARRAY_FOR_EACH(files_list) {

        char *file_name = files_list[i];
        char full_input_path[strlen(input_dir) + strlen(file_name) + 2];

        sprintf(full_input_path, "%s/%s", input_dir, file_name);

        vtk_grid = new_vtk_unstructured_grid_from_file(full_input_path, false);

        if(vtk_grid) {

            int64_t *cells = vtk_grid->cells;
            point3d_array points = vtk_grid->points;

            uint32_t n_active = vtk_grid->num_cells;

            int num_points = vtk_grid->points_per_cell;
            int j = num_points;

            struct point_3d cell_coordinates;

            for(uint32_t i = 0; i < n_active * num_points; i += num_points) {
                float center_x, center_y, center_z;
                float dx, dy, dz, v;

                dx = fabs((points[cells[i]].x - points[cells[i + 1]].x));
                dy = fabs((points[cells[i]].y - points[cells[i + 3]].y));
                dz = fabs((points[cells[i]].z - points[cells[i + 4]].z));

                center_x = points[cells[i]].x + dx / 2.0f;
                center_y = points[cells[i]].y + dy / 2.0f;
                center_z = points[cells[i]].z + dz / 2.0f;

                cell_coordinates.x = center_x;
                cell_coordinates.y = center_y;
                cell_coordinates.z = center_z;

                v = vtk_grid->values[j - num_points];
                j += 1;

                int n_activations = 0;
                float *apds_array = NULL;
                float *activation_times_array = NULL;

                float last_v = hmget(last_time_v, cell_coordinates);

                n_activations = (int)hmget(num_activations, cell_coordinates);
                activation_times_array = (float *)hmget(activation_times, cell_coordinates);
                apds_array = (float *)hmget(apds, cell_coordinates);

                int act_times_len = arrlen(activation_times_array);

                if(current_t == 0.0f) {
                    hmput(last_time_v, cell_coordinates, v);
                } else {
                    if((last_v < activation_threshold) && (v >= activation_threshold)) {

                        if(act_times_len == 0) {
                            n_activations++;
                            hmput(num_activations, cell_coordinates, n_activations);
                            arrput(activation_times_array, current_t);
                            float tmp = hmget(cell_was_active, cell_coordinates);
                            hmput(cell_was_active, cell_coordinates, tmp + 1);
                            hmput(activation_times, cell_coordinates, activation_times_array);
                        } else { // This is to avoid spikes in the middle of an Action Potential
                            float last_act_time = activation_times_array[act_times_len - 1];
                            if(current_t - last_act_time > time_threshold) {
                                n_activations++;
                                hmput(num_activations, cell_coordinates, n_activations);
                                arrput(activation_times_array, current_t);
                                float tmp = hmget(cell_was_active, cell_coordinates);
                                hmput(cell_was_active, cell_coordinates, tmp + 1);
                                hmput(activation_times, cell_coordinates, activation_times_array);
                            }
                        }
                    }

                    // CHECK APD
                    bool was_active = (hmget(cell_was_active, cell_coordinates) != 0.0);
                    if(was_active) {
                        if(v <= apd_threshold || (hmget(cell_was_active, cell_coordinates) == 2.0) || (final_time - current_t) <= dt) {

                            int tmp = (int)hmget(cell_was_active, cell_coordinates);
                            int act_time_array_len = arrlen(activation_times_array);
                            // if this in being calculated because we had a new activation before the cell achieved the rest potential,
                            // we need to get the activation before this one
                            real_cpu last_act_time = activation_times_array[act_time_array_len - tmp];
                            real_cpu apd = current_t - last_act_time;
                            arrput(apds_array, apd);
                            hmput(apds, cell_coordinates, apds_array);
                            hmput(cell_was_active, cell_coordinates, tmp - 1);
                        }
                    }

                    hmput(last_time_v, cell_coordinates, v);
                }
            }
        } else {
            fprintf(stderr, "[WARN] Failed to processs %s/%s\n", input_dir, files_list[i]);
        }

        current_t += dt;
    }

    FILE *act_file = stdout;

    for(int i = 0; i < hmlen(num_activations); i++) {

        int n_activations = 0;
        float *apds_array = NULL;
        float *activation_times_array = NULL;

        n_activations = num_activations[i].value;
        activation_times_array = (float *)hmget(activation_times, num_activations[i].key);
        apds_array = (float *)hmget(apds, num_activations[i].key);

        fprintf(act_file, "%lf %lf %lf ", num_activations[i].key.x, num_activations[i].key.y, num_activations[i].key.z);

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

    free(input_type);
    free(input_dir);
    free(output_dir);

    return 1;
}
