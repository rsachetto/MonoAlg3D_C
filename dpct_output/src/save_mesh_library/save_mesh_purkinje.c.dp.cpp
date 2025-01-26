//
// Created by bergolho on 30/09/20.
//

#include <sycl/sycl.hpp>
#include <dpct/dpct.hpp>
#include <stdbool.h>
#include <stdlib.h>

#include "../3dparty/sds/sds.h"
#include "../alg/grid/grid.h"
#include "../config/save_mesh_config.h"
#include "../utils/utils.h"
#include "../extra_data_library/helper_functions.h"
#include "../domains_library/mesh_info_data.h"

#include "../libraries_common/common_data_structures.h"
#include "../ensight_utils/ensight_grid.h"
#include "../vtk_utils/vtk_polydata_grid.h"
#include "../vtk_utils/vtk_unstructured_grid.h"

#include "save_mesh_helper.h"

#ifdef COMPILE_CUDA
#include "../gpu_utils/gpu_utils.h"
#endif

static char *file_prefix;
static char *file_prefix_purkinje;
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
static bool initialized = false;
static bool save_ode_state_variables = false;

struct save_as_vtk_or_vtu_persistent_data {
    struct vtk_unstructured_grid *grid;
    bool first_save_call;
};

static void save_visibility_mask(sds output_dir_with_file, ui8_array visible_cells) {
    sds output_dir_with_new_file = sdsnew(output_dir_with_file);
    output_dir_with_new_file = sdscat(output_dir_with_new_file, ".vis");
    FILE *vis = fopen(output_dir_with_new_file, "wb");
    fwrite(visible_cells, sizeof(uint8_t), arrlen(visible_cells), vis);
    sdsfree(output_dir_with_new_file);
    fclose(vis);
}

INIT_SAVE_MESH(init_save_as_vtk_or_vtp) {
    config->persistent_data = malloc(sizeof(struct save_as_vtp_persistent_data));
    ((struct save_as_vtp_persistent_data *)config->persistent_data)->grid = NULL;
    ((struct save_as_vtp_persistent_data *)config->persistent_data)->first_save_call = true;
}

END_SAVE_MESH(end_save_as_vtk_or_vtp) {
    free_vtk_polydata_grid(((struct save_as_vtp_persistent_data *)config->persistent_data)->grid);
    free(config->persistent_data);
}

SAVE_MESH(save_as_vtp_purkinje) {

    int iteration_count = time_info->iteration;

    if(((struct save_as_vtp_persistent_data *)config->persistent_data)->first_save_call) {
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config, "output_dir");
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(file_prefix, config, "file_prefix");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_plain, config, "clip_with_plain");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_bounds, config, "clip_with_bounds");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(binary, config, "binary");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_pvd, config, "save_pvd");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(compress, config, "compress");
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, compression_level, config, "compression_level");

        if(compress)
            binary = true;

        if(!save_pvd) {
            ((struct save_as_vtp_persistent_data *)config->persistent_data)->first_save_call = false;
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
    sds base_name = create_base_name(file_prefix, iteration_count, "vtp");

    real_cpu current_t = time_info->current_t;

    output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, current_t);

    if(save_pvd) {
        add_file_to_pvd(current_t, output_dir, base_name, ((struct save_as_vtp_persistent_data *)config->persistent_data)->first_save_call);
        ((struct save_as_vtp_persistent_data *)config->persistent_data)->first_save_call = false;
    }

    bool read_only_data = ((struct save_as_vtp_persistent_data *)config->persistent_data)->grid != NULL;
    new_vtk_polydata_grid_from_purkinje_grid(&((struct save_as_vtp_persistent_data *)config->persistent_data)->grid, the_grid->purkinje, clip_with_plain,
                                             plain_coords, clip_with_bounds, bounds, read_only_data);

    if(compress) {
        save_vtk_polydata_grid_as_vtp_compressed(((struct save_as_vtp_persistent_data *)config->persistent_data)->grid, output_dir_with_file,
                                                 compression_level);
    } else {
        save_vtk_polydata_grid_as_vtp(((struct save_as_vtp_persistent_data *)config->persistent_data)->grid, output_dir_with_file, binary);
    }

    sdsfree(output_dir_with_file);
    sdsfree(base_name);
}

SAVE_MESH(save_as_vtk_purkinje) {

    char *output_dir;
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config, "output_dir");

    if(!initialized) {

        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(file_prefix, config, "file_prefix");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_plain, config, "clip_with_plain");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_bounds, config, "clip_with_bounds");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(binary, config, "binary");
        initialized = true;
    }
    float plain_coords[6] = {0, 0, 0, 0, 0, 0};
    float bounds[6] = {0, 0, 0, 0, 0, 0};

    if(clip_with_plain) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[0], config, "origin_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[1], config, "origin_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[2], config, "origin_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[3], config, "normal_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[3], config, "normal_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[4], config, "normal_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[5], config, "normal_z");
    }

    if(clip_with_bounds) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[0], config, "min_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[1], config, "min_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[2], config, "min_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[3], config, "max_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[4], config, "max_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[5], config, "max_z");
    }

    sds output_dir_with_file = sdsnew(output_dir);
    output_dir_with_file = sdscat(output_dir_with_file, "/");
    sds base_name = create_base_name(file_prefix, time_info->iteration, "vtk");
    output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, time_info->current_t);

    bool read_only_data = ((struct save_as_vtp_persistent_data *)config->persistent_data)->grid != NULL;
    new_vtk_polydata_grid_from_purkinje_grid(&((struct save_as_vtp_persistent_data *)config->persistent_data)->grid, the_grid->purkinje, clip_with_plain,
                                             plain_coords, clip_with_bounds, bounds, read_only_data);

    save_vtk_polydata_grid_as_legacy_vtk(((struct save_as_vtp_persistent_data *)config->persistent_data)->grid, output_dir_with_file, binary);

    sdsfree(output_dir_with_file);
    sdsfree(base_name);
}

INIT_SAVE_MESH(init_save_tissue_as_vtk_or_vtu_purkinje_as_vtp) {
    config->persistent_data = malloc(sizeof(struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data));
    ((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *)config->persistent_data)->grid = NULL;
    ((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *)config->persistent_data)->grid_purkinje = NULL;
    ((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *)config->persistent_data)->first_save_call = true;
}

END_SAVE_MESH(end_save_tissue_as_vtk_or_vtu_purkinje_as_vtp) {
    free_vtk_polydata_grid(((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *)config->persistent_data)->grid_purkinje);
    free_vtk_unstructured_grid(((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *)config->persistent_data)->grid);
    free(config->persistent_data);
}

SAVE_MESH(save_tissue_as_vtu_purkinje_as_vtp) {

    // [TISSUE]
    int iteration_count = time_info->iteration;

    if(((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *)config->persistent_data)->first_save_call) {
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config, "output_dir");
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(file_prefix, config, "file_prefix");
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(file_prefix_purkinje, config, "file_prefix_purkinje");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_plain, config, "clip_with_plain");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_bounds, config, "clip_with_bounds");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(binary, config, "binary");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_pvd, config, "save_pvd");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(compress, config, "compress");
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, compression_level, config, "compression_level");

        if(compress)
            binary = true;

        if(!save_pvd) {
            ((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *)config->persistent_data)->first_save_call = false;
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
        add_file_to_pvd(current_t, output_dir, base_name,
                        ((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *)config->persistent_data)->first_save_call);
        ((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *)config->persistent_data)->first_save_call = false;
    }

    bool read_only_data = ((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *)config->persistent_data)->grid != NULL;
    new_vtk_unstructured_grid_from_alg_grid(&((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *)config->persistent_data)->grid, the_grid,
                                            clip_with_plain, plain_coords, clip_with_bounds, bounds, read_only_data, save_f, false, NULL);

    if(compress) {
        save_vtk_unstructured_grid_as_vtu_compressed(((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *)config->persistent_data)->grid,
                                                     output_dir_with_file, compression_level);
    } else {
        save_vtk_unstructured_grid_as_vtu(((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *)config->persistent_data)->grid, output_dir_with_file,
                                          binary);
    }

    if(the_grid->adaptive) {
        free_vtk_unstructured_grid(((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *)config->persistent_data)->grid);
        ((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *)config->persistent_data)->grid = NULL;
    }

    // [PURKINJE]

    output_dir_with_file = sdsnew(output_dir);
    output_dir_with_file = sdscat(output_dir_with_file, "/");
    base_name = create_base_name(file_prefix_purkinje, iteration_count, "vtp");

    output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, current_t);

    if(save_pvd) {
        add_file_to_pvd(current_t, output_dir, base_name,
                        ((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *)config->persistent_data)->first_save_call);
        ((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *)config->persistent_data)->first_save_call = false;
    }

    read_only_data = ((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *)config->persistent_data)->grid_purkinje != NULL;
    new_vtk_polydata_grid_from_purkinje_grid(&((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *)config->persistent_data)->grid_purkinje,
                                             the_grid->purkinje, clip_with_plain, plain_coords, clip_with_bounds, bounds, read_only_data);

    if(compress) {
        save_vtk_polydata_grid_as_vtp_compressed(((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *)config->persistent_data)->grid_purkinje,
                                                 output_dir_with_file, compression_level);
    } else {
        save_vtk_polydata_grid_as_vtp(((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *)config->persistent_data)->grid_purkinje,
                                      output_dir_with_file, binary);
    }

    sdsfree(output_dir_with_file);
    sdsfree(base_name);
}

INIT_SAVE_MESH(init_save_purkinje_with_activation_times) {

    config->persistent_data = calloc(1, sizeof(struct save_coupling_with_activation_times_persistent_data));

    hmdefault(((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->purkinje_cell_was_active, 0.0);
    hmdefault(((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->purkinje_last_time_v, -100.0);
    hmdefault(((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->purkinje_num_activations, 0);
    hmdefault(((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->purkinje_activation_times, NULL);
    hmdefault(((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->purkinje_apds, NULL);
    ((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->purkinje_grid = NULL;

    ((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->first_save_call = true;
}

END_SAVE_MESH(end_save_purkinje_with_activation_times) {

    bool save_activation_time_map = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_activation_time_map, config, "save_activation_time");

    bool save_apd_map = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_apd_map, config, "save_apd");

    bool save_purkinje_velocity = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_purkinje_velocity, config, "save_purkinje_velocity");

    if(save_activation_time_map) {
        log_info("[!] Saving activation time maps !!!!\n");
        write_purkinje_activation_time_maps(config, the_grid, output_dir, file_prefix_purkinje, clip_with_plain, clip_with_bounds, binary, compress,
                                            compression_level);
    }

    if(save_apd_map) {
        log_info("[!] Saving APD map !!!!\n");
        write_purkinje_apd_map(config, the_grid, output_dir, file_prefix_purkinje, clip_with_plain, clip_with_bounds, binary, compress, compression_level);
    }

    if(save_purkinje_velocity) {
        log_info("[!] Calculating Purkinje propagation velocity !!!!\n");
        print_purkinje_propagation_velocity(config, the_grid);
    }

    free(config->persistent_data);
}

SAVE_MESH(save_purkinje_with_activation_times) {

    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config, "output_dir");

    float time_threshold = 10.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, time_threshold, config, "time_threshold");

    float purkinje_activation_threshold = -30.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, purkinje_activation_threshold, config, "activation_threshold_purkinje");

    float purkinje_apd_threshold = -83.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, purkinje_apd_threshold, config, "apd_threshold_purkinje");

    calculate_purkinje_activation_time_and_apd(time_info, config, the_grid, time_threshold, purkinje_activation_threshold, purkinje_apd_threshold);
}

INIT_SAVE_MESH(init_save_purkinje_coupling_with_activation_times) {

    config->persistent_data = calloc(1, sizeof(struct save_coupling_with_activation_times_persistent_data));
    hmdefault(((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->tissue_cell_was_active, 0.0);
    hmdefault(((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->tissue_last_time_v, -100.0);
    hmdefault(((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->tissue_num_activations, 0);
    hmdefault(((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->tissue_activation_times, NULL);
    hmdefault(((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->tissue_apds, NULL);
    ((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->tissue_grid = NULL;

    hmdefault(((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->purkinje_cell_was_active, 0.0);
    hmdefault(((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->purkinje_last_time_v, -100.0);
    hmdefault(((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->purkinje_num_activations, 0);
    hmdefault(((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->purkinje_activation_times, NULL);
    hmdefault(((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->purkinje_apds, NULL);
    ((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->purkinje_grid = NULL;

    ((struct save_coupling_with_activation_times_persistent_data *)config->persistent_data)->first_save_call = true;
}

END_SAVE_MESH(end_save_purkinje_coupling_with_activation_times) {

    bool save_activation_time_map = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_activation_time_map, config, "save_activation_time");

    bool save_apd_map = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_apd_map, config, "save_apd");

    bool save_purkinje_velocity = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_purkinje_velocity, config, "save_purkinje_velocity");

    if(save_activation_time_map) {
        log_info("[!] Saving activation time maps !!!!\n");
        write_tissue_activation_time_maps(config, the_grid, output_dir, file_prefix, clip_with_plain, clip_with_bounds, binary, compress,
                                          compression_level, save_f);
        write_purkinje_activation_time_maps(config, the_grid, output_dir, file_prefix_purkinje, clip_with_plain, clip_with_bounds, binary,  compress,
                                            compression_level);
    }

    if(save_apd_map) {
        log_info("[!] Saving APD map !!!!\n");
        write_tissue_apd_map(config, the_grid, output_dir, file_prefix, clip_with_plain, clip_with_bounds, binary, compress, compression_level,
                             save_f);
        write_purkinje_apd_map(config, the_grid, output_dir, file_prefix_purkinje, clip_with_plain, clip_with_bounds, binary, compress,
                               compression_level);
    }

    if(save_purkinje_velocity) {
        log_info("[!] Calculating Purkinje propagation velocity !!!!\n");
        print_purkinje_propagation_velocity(config, the_grid);
    }

    free(config->persistent_data);
}

SAVE_MESH(save_purkinje_coupling_with_activation_times) {

    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config, "output_dir");

    float time_threshold = 10.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, time_threshold, config, "time_threshold");

    float tissue_activation_threshold = -30.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, tissue_activation_threshold, config, "activation_threshold_tissue");

    float tissue_apd_threshold = -83.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, tissue_apd_threshold, config, "apd_threshold_tissue");

    float purkinje_activation_threshold = -30.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, purkinje_activation_threshold, config, "activation_threshold_purkinje");

    float purkinje_apd_threshold = -83.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, purkinje_apd_threshold, config, "apd_threshold_purkinje");

    // [TISSUE]
    calculate_tissue_activation_time_and_apd(time_info, config, the_grid, time_threshold, tissue_activation_threshold, tissue_apd_threshold);

    // [PURKINJE]
    calculate_purkinje_activation_time_and_apd(time_info, config, the_grid, time_threshold, purkinje_activation_threshold, purkinje_apd_threshold);
}

INIT_SAVE_MESH(init_save_tissue_with_activation_times) {

    config->persistent_data = calloc(1, sizeof(struct save_coupling_with_activation_times_persistent_data));
    hmdefault(((struct save_coupling_with_activation_times_persistent_data*)config->persistent_data)->tissue_cell_was_active, 0.0);
    hmdefault(((struct save_coupling_with_activation_times_persistent_data*)config->persistent_data)->tissue_last_time_v, -100.0);
    hmdefault(((struct save_coupling_with_activation_times_persistent_data*)config->persistent_data)->tissue_num_activations, 0);
    hmdefault(((struct save_coupling_with_activation_times_persistent_data*)config->persistent_data)->tissue_activation_times, NULL);
    hmdefault(((struct save_coupling_with_activation_times_persistent_data*)config->persistent_data)->tissue_apds, NULL);
    ((struct save_coupling_with_activation_times_persistent_data*) config->persistent_data)->tissue_grid = NULL;

    ((struct save_coupling_with_activation_times_persistent_data*)config->persistent_data)->first_save_call = true;
}

END_SAVE_MESH(end_save_tissue_with_activation_times) {

    bool save_activation_time_map = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_activation_time_map, config, "save_activation_time");

    bool save_apd_map = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_apd_map, config, "save_apd");

    if (save_activation_time_map) {
        log_info("[!] Saving activation time maps !!!!\n");
        write_tissue_activation_time_maps(config,the_grid,output_dir,file_prefix,clip_with_plain,clip_with_bounds,binary,compress,compression_level,save_f);
    }

    if (save_apd_map) {
        log_info("[!] Saving APD map !!!!\n");
        write_tissue_apd_map(config,the_grid,output_dir,file_prefix,clip_with_plain,clip_with_bounds,binary,compress,compression_level,save_f);
    }

    free(config->persistent_data);

}

SAVE_MESH (save_tissue_with_activation_times) {

    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config, "output_dir");

    float time_threshold = 10.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, time_threshold, config, "time_threshold");

    float tissue_activation_threshold = -30.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, tissue_activation_threshold, config, "activation_threshold_tissue");

    float tissue_apd_threshold = -83.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, tissue_apd_threshold, config, "apd_threshold_tissue");

    calculate_tissue_activation_time_and_apd(time_info,config,the_grid,time_threshold,tissue_activation_threshold,tissue_apd_threshold);
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

END_SAVE_MESH(end_save_one_cell_state_variables) {
    free(((struct save_one_cell_state_variables_persistent_data *)config->persistent_data)->file_name);
    fclose(((struct save_one_cell_state_variables_persistent_data *)config->persistent_data)->file);
    free(config->persistent_data);
}

SAVE_MESH(save_one_cell_state_variables) {

    struct save_one_cell_state_variables_persistent_data *params = ((struct save_one_cell_state_variables_persistent_data *)config->persistent_data);
    real_cpu cur_time = time_info->current_t;
    real_cpu final_time = time_info->final_t;
    real_cpu dt = time_info->dt;
    real_cpu bcl = 1000.0;

    if(params->cell_sv_position == -1) {
        if(!the_grid->adaptive) {
            FOR_EACH_PURKINJE_CELL(the_grid) {
            //FOR_EACH_CELL(the_grid) {
                if(cell->center.x == params->cell_center_x && cell->center.y == params->cell_center_y && cell->center.z == params->cell_center_z) {
                    params->cell_sv_position = cell->sv_position;
                    printf("%d\n", params->cell_sv_position);
                    break;
                }
            }
        }
    }

    //if (cur_time+dt >= final_time-bcl) {
        if(ode_solver->gpu) {
    #ifdef COMPILE_CUDA
            int num_odes = ode_solver->model_data.number_of_ode_equations;
            real *cell_sv;

            cell_sv = (real *)malloc(sizeof(real) * num_odes);

            check_cuda_error(DPCT_CHECK_ERROR(dpct::dpct_memcpy(cell_sv, sizeof(real), ode_solver->sv + params->cell_sv_position, ode_solver->pitch,
                                                                sizeof(real), ode_solver->model_data.number_of_ode_equations, cudaMemcpyDeviceToHost)));

            fprintf(params->file, "%lf ", time_info->current_t);
            for(int i = 0; i < num_odes; i++) {
                fprintf(params->file, "%lf ", cell_sv[i]);
            }
            fprintf(params->file, "\n");

            free(cell_sv);
    #endif
        } else {

            int num_odes = ode_solver->model_data.number_of_ode_equations;
            real *cell_sv = &ode_solver->sv[params->cell_sv_position * num_odes];

            fprintf(params->file, "%lf ", time_info->current_t);
            for(int i = 0; i < num_odes; i++) {
                fprintf(params->file, "%lf ", cell_sv[i]);
            }
            fprintf(params->file, "\n");
        }
    //}
}

INIT_SAVE_MESH(init_save_multiple_cell_state_variables) {
    config->persistent_data = malloc(sizeof(struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data));
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(uint32_t, ((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->num_tissue_cells,
                                                config, "tissue_num_cells");
    GET_PARAMETER_MATRIX_VALUE_OR_USE_DEFAULT(((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->tissue_cell_centers,
                                                config, "tissue_cell_centers", ((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->num_tissue_cells, 3);
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->tissue_file_name_prefix, config,
                                               "tissue_file_name_prefix");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(uint32_t, ((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->num_purkinje_cells,
                                                config, "purkinje_num_cells");
    GET_PARAMETER_MATRIX_VALUE_OR_USE_DEFAULT(((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->purkinje_cell_centers,
                                                config, "purkinje_cell_centers", ((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->num_purkinje_cells, 3);
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->purkinje_file_name_prefix, config,
                                               "purkinje_file_name_prefix");

    uint32_t tissue_num_cells = ((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->num_tissue_cells;
    char *tissue_file_name_prefix = ((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->tissue_file_name_prefix;
    uint32_t *tissue_cell_sv_positions = ((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->tissue_cell_sv_positions;
    ((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->tissue_files = MALLOC_ARRAY_OF_TYPE(FILE*, tissue_num_cells);
    ((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->tissue_cell_sv_positions = MALLOC_ARRAY_OF_TYPE(uint32_t, tissue_num_cells);

    uint32_t purkinje_num_cells = ((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->num_purkinje_cells;
    char *purkinje_file_name_prefix = ((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->purkinje_file_name_prefix;
    uint32_t *purkinje_cell_sv_positions = ((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->purkinje_cell_sv_positions;
    ((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->purkinje_files = MALLOC_ARRAY_OF_TYPE(FILE*, purkinje_num_cells);
    ((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->purkinje_cell_sv_positions = MALLOC_ARRAY_OF_TYPE(uint32_t, purkinje_num_cells);

    for (int i = 0; i < tissue_num_cells; i++) {
        
        sds base_name = NULL;
        base_name = create_base_name(tissue_file_name_prefix, i, "dat");

        ((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->tissue_files[i] = fopen(base_name, "w");
        ((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->tissue_cell_sv_positions[i] = -1;
    }

    for (int i = 0; i < purkinje_num_cells; i++) {
        
        sds base_name = NULL;
        base_name = create_base_name(purkinje_file_name_prefix, i, "dat");

        ((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->purkinje_files[i] = fopen(base_name, "w");
        ((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->purkinje_cell_sv_positions[i] = -1;
    }
}

SAVE_MESH(save_multiple_cell_state_variables) {
    struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *params = ((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data);

// [DOMAIN]
    for (uint32_t i = 0; i < params->num_tissue_cells; i++) {
        if(params->tissue_cell_sv_positions[i] == -1) {
            if(!the_grid->adaptive) {
                FOR_EACH_CELL(the_grid) {    
                    if (cell->active) {
                        if(cell->center.x == params->tissue_cell_centers[i*3] && cell->center.y == params->tissue_cell_centers[i*3+1] && cell->center.z == params->tissue_cell_centers[i*3+2]) {
                            params->tissue_cell_sv_positions[i] = cell->sv_position;
                            break;
                        }
                    }
                }
            }
        }
        //printf("[tissue] Found cell %u with coordinates: (%lf %lf %lf)\n", params->tissue_cell_sv_positions[i], \
                                                                params->tissue_cell_centers[i*3], \
                                                                params->tissue_cell_centers[i*3+1], \
                                                                params->tissue_cell_centers[i*3+2]);
    }
    
    if(ode_solver->gpu) {
#ifdef COMPILE_CUDA
        for (uint32_t k = 0; k < params->num_tissue_cells; k++) {
            real *cell_sv;

            cell_sv = MALLOC_ARRAY_OF_TYPE(real, ode_solver->model_data.number_of_ode_equations);

            check_cuda_error(DPCT_CHECK_ERROR(dpct::dpct_memcpy(cell_sv, sizeof(real), ode_solver->sv + params->tissue_cell_sv_positions[k], ode_solver->pitch,
                                                                sizeof(real), ode_solver->model_data.number_of_ode_equations, cudaMemcpyDeviceToHost)));

            // Only 'time' and 'Vm'
            fprintf(params->tissue_files[k], "%g %g\n", time_info->current_t, cell_sv[0]);

            // All state-vectors
            //fprintf(params->tissue_files[k], "%lf ", time_info->current_t);
            //for(int i = 0; i < ode_solver->model_data.number_of_ode_equations; i++) {
            //    fprintf(params->tissue_files[k], "%lf ", cell_sv[i]);
            //}
            //fprintf(params->tissue_files[k], "\n");

            free(cell_sv);
        }
        
#endif
    } else {
        for (uint32_t k = 0; k < params->num_tissue_cells; k++) {
            real *cell_sv = &ode_solver->sv[params->tissue_cell_sv_positions[k] * ode_solver->model_data.number_of_ode_equations];

            // Only 'time' and 'Vm'
            fprintf(params->tissue_files[k], "%g %g\n", time_info->current_t, cell_sv[0]);

            // All state-vectors
            //fprintf(params->tissue_files[k], "%.3lf ", time_info->current_t);
            //for(int i = 0; i < ode_solver->model_data.number_of_ode_equations; i++) {
            //    fprintf(params->tissue_files[k], "%g ", cell_sv[i]);
            //}
            //fprintf(params->tissue_files[k], "\n");
        }
    }

// [PURKINJE]
    for (uint32_t i = 0; i < params->num_purkinje_cells; i++) {
        if(params->purkinje_cell_sv_positions[i] == -1) {
            FOR_EACH_PURKINJE_CELL(the_grid) {
                real_cpu dist = calc_norm(cell->center.x,cell->center.y,cell->center.z,\
                                        params->purkinje_cell_centers[i*3],params->purkinje_cell_centers[i*3+1],params->purkinje_cell_centers[i*3+2]);
                if (dist < 1e-1) {
                    params->purkinje_cell_sv_positions[i] = cell->sv_position;
                    //printf("[purkinje] Found cell %u with coordinates: (%lf %lf %lf)\n", params->purkinje_cell_sv_positions[i], \
                                                                params->purkinje_cell_centers[i*3], \
                                                                params->purkinje_cell_centers[i*3+1], \
                                                                params->purkinje_cell_centers[i*3+2]);
                    break;
                }
            }
        }
    }
    
    if(purkinje_ode_solver->gpu) {
#ifdef COMPILE_CUDA
        for (uint32_t k = 0; k < params->num_purkinje_cells; k++) {
            real *cell_sv;

            cell_sv = MALLOC_ARRAY_OF_TYPE(real, purkinje_ode_solver->model_data.number_of_ode_equations);

            check_cuda_error(DPCT_CHECK_ERROR(dpct::dpct_memcpy(cell_sv, sizeof(real), purkinje_ode_solver->sv + params->purkinje_cell_sv_positions[k],
                                                                purkinje_ode_solver->pitch, sizeof(real),
                                                                purkinje_ode_solver->model_data.number_of_ode_equations, cudaMemcpyDeviceToHost)));

            // Only 'time' and 'Vm'
            fprintf(params->purkinje_files[k], "%g %g\n", time_info->current_t, cell_sv[0]);

            // All state-vectors
            //fprintf(params->purkinje_files[k], "%.3lf ", time_info->current_t);
            //for(int i = 0; i < purkinje_ode_solver->model_data.number_of_ode_equations; i++) {
            //    fprintf(params->purkinje_files[k], "%g ", cell_sv[i]);
            //}
            //fprintf(params->purkinje_files[k], "\n");

            free(cell_sv);
        }
        
#endif
    } else {
        for (uint32_t k = 0; k < params->num_tissue_cells; k++) {
            real *cell_sv = &purkinje_ode_solver->sv[params->purkinje_cell_sv_positions[k] * purkinje_ode_solver->model_data.number_of_ode_equations];

            // Only 'time' and 'Vm'
            fprintf(params->purkinje_files[k], "%g %g\n", time_info->current_t, cell_sv[0]);

            // All state-vectos
            //fprintf(params->purkinje_files[k], "%.3lf ", time_info->current_t);
            //for(int i = 0; i < purkinje_ode_solver->model_data.number_of_ode_equations; i++) {
            //    fprintf(params->purkinje_files[k], "%g ", cell_sv[i]);
            //}
            //fprintf(params->purkinje_files[k], "\n");
        }
    }
}

END_SAVE_MESH(end_save_multiple_cell_state_variables) {
    uint32_t num_tissue_cells = ((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->num_tissue_cells;
    free(((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->tissue_file_name_prefix);
    free(((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->tissue_cell_sv_positions);
    free(((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->tissue_cell_centers);
    for (uint32_t i = 0; i < num_tissue_cells; i++) {
        fclose(((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->tissue_files[i]);
    }
    free(((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->tissue_files);

    uint32_t num_purkinje_cells = ((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->num_purkinje_cells;
    free(((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->purkinje_file_name_prefix);
    free(((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->purkinje_cell_sv_positions);
    free(((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->purkinje_cell_centers);
    for (uint32_t i = 0; i < num_purkinje_cells; i++) {
        fclose(((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->purkinje_files[i]);
    }
    free(((struct save_multiple_cell_state_variables_purkinje_coupling_persistent_data *)config->persistent_data)->purkinje_files);

    free(config->persistent_data);
}

INIT_SAVE_MESH(init_save_as_vtk_or_vtu) {
    config->persistent_data = malloc(sizeof(struct save_as_vtk_or_vtu_persistent_data));
    ((struct save_as_vtk_or_vtu_persistent_data *)config->persistent_data)->grid = NULL;
    ((struct save_as_vtk_or_vtu_persistent_data *)config->persistent_data)->first_save_call = true;
}

END_SAVE_MESH(end_save_as_vtk_or_vtu) {
    free_vtk_unstructured_grid(((struct save_as_vtk_or_vtu_persistent_data *)config->persistent_data)->grid);
    free(config->persistent_data);
}

SAVE_MESH(save_transmurality_as_vtk) {

    struct extra_data_for_torord *extra_data = (struct extra_data_for_torord*)ode_solver->ode_extra_data;
    real *transmurality = extra_data->transmurality;

    int iteration_count = time_info->iteration;

    if(((struct save_as_vtk_or_vtu_persistent_data *)config->persistent_data)->first_save_call) {
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config, "output_dir");
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(file_prefix, config, "file_prefix");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_plain, config, "clip_with_plain");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_bounds, config, "clip_with_bounds");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(binary, config, "binary");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_f, config, "save_f");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_visible_mask, config, "save_visible_mask");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_scar_cells, config, "save_scar_cells");

        ((struct save_as_vtk_or_vtu_persistent_data *)config->persistent_data)->first_save_call = false;
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

    bool read_only_data = ((struct save_as_vtk_or_vtu_persistent_data *)config->persistent_data)->grid != NULL;

    new_vtk_unstructured_grid_from_alg_grid(&(((struct save_as_vtk_or_vtu_persistent_data *)config->persistent_data)->grid), the_grid, clip_with_plain,
                                            plain_coords, clip_with_bounds, bounds, read_only_data, save_f, save_scar_cells, transmurality);

    save_vtk_unstructured_grid_as_legacy_vtk(((struct save_as_vtk_or_vtu_persistent_data *)config->persistent_data)->grid, output_dir_with_file, binary, save_f, NULL);

    if(save_visible_mask) {
        save_visibility_mask(output_dir_with_file, (((struct save_as_vtk_or_vtu_persistent_data *)config->persistent_data)->grid)->cell_visibility);
    }

    if(the_grid->adaptive) {
        free_vtk_unstructured_grid(((struct save_as_vtk_or_vtu_persistent_data *)config->persistent_data)->grid);
        ((struct save_as_vtk_or_vtu_persistent_data *)config->persistent_data)->grid = NULL;
    }

    sdsfree(output_dir_with_file);
    sdsfree(base_name);

    CALL_EXTRA_FUNCTIONS(save_mesh_fn, time_info, config, the_grid, ode_solver, purkinje_ode_solver);

    exit(EXIT_SUCCESS);
}

INIT_SAVE_MESH(init_save_multiple_cell_state_variables_purkinje_coupling_with_activation_times) {

    config->persistent_data = calloc(1, sizeof(struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data));
    
    // [PURKINJE COUPLED MULTIPLE CELLS]
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(uint32_t, ((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->num_tissue_cells,
                                                config, "tissue_num_cells");
    GET_PARAMETER_MATRIX_VALUE_OR_USE_DEFAULT(((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->tissue_cell_centers,
                                                config, "tissue_cell_centers", ((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->num_tissue_cells, 3);
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->tissue_file_name_prefix, config,
                                               "tissue_file_name_prefix");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(uint32_t, ((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->num_purkinje_cells,
                                                config, "purkinje_num_cells");
    GET_PARAMETER_MATRIX_VALUE_OR_USE_DEFAULT(((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->purkinje_cell_centers,
                                                config, "purkinje_cell_centers", ((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->num_purkinje_cells, 3);
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->purkinje_file_name_prefix, config,
                                               "purkinje_file_name_prefix");

    uint32_t tissue_num_cells = ((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->num_tissue_cells;
    char *tissue_file_name_prefix = ((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->tissue_file_name_prefix;
    uint32_t *tissue_cell_sv_positions = ((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->tissue_cell_sv_positions;
    ((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->tissue_files = MALLOC_ARRAY_OF_TYPE(FILE*, tissue_num_cells);
    ((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->tissue_cell_sv_positions = MALLOC_ARRAY_OF_TYPE(uint32_t, tissue_num_cells);

    uint32_t purkinje_num_cells = ((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->num_purkinje_cells;
    char *purkinje_file_name_prefix = ((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->purkinje_file_name_prefix;
    uint32_t *purkinje_cell_sv_positions = ((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->purkinje_cell_sv_positions;
    ((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->purkinje_files = MALLOC_ARRAY_OF_TYPE(FILE*, purkinje_num_cells);
    ((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->purkinje_cell_sv_positions = MALLOC_ARRAY_OF_TYPE(uint32_t, purkinje_num_cells);

    for (int i = 0; i < tissue_num_cells; i++) {
        
        sds base_name = NULL;
        base_name = create_base_name(tissue_file_name_prefix, i, "dat");

        ((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->tissue_files[i] = fopen(base_name, "w");
        ((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->tissue_cell_sv_positions[i] = -1;
    }

    for (int i = 0; i < purkinje_num_cells; i++) {
        
        sds base_name = NULL;
        base_name = create_base_name(purkinje_file_name_prefix, i, "dat");

        ((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->purkinje_files[i] = fopen(base_name, "w");
        ((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->purkinje_cell_sv_positions[i] = -1;
    }

    // [PURKINJE COUPLED ACTIVATION TIMES]
    hmdefault(((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->tissue_cell_was_active, 0.0);
    hmdefault(((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->tissue_last_time_v, -100.0);
    hmdefault(((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->tissue_num_activations, 0);
    hmdefault(((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->tissue_activation_times, NULL);
    hmdefault(((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->tissue_apds, NULL);
    ((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->tissue_grid = NULL;

    hmdefault(((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->purkinje_cell_was_active, 0.0);
    hmdefault(((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->purkinje_last_time_v, -100.0);
    hmdefault(((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->purkinje_num_activations, 0);
    hmdefault(((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->purkinje_activation_times, NULL);
    hmdefault(((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->purkinje_apds, NULL);
    ((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->purkinje_grid = NULL;

    ((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->first_save_call = true;
}

END_SAVE_MESH(end_save_multiple_cell_state_variables_purkinje_coupling_with_activation_times) {

    // [PURKINJE COUPLED MULTIPLE CELLS]
    uint32_t num_tissue_cells = ((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->num_tissue_cells;
    free(((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->tissue_file_name_prefix);
    free(((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->tissue_cell_sv_positions);
    free(((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->tissue_cell_centers);
    for (uint32_t i = 0; i < num_tissue_cells; i++) {
        fclose(((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->tissue_files[i]);
    }
    free(((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->tissue_files);

    uint32_t num_purkinje_cells = ((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->num_purkinje_cells;
    free(((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->purkinje_file_name_prefix);
    free(((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->purkinje_cell_sv_positions);
    free(((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->purkinje_cell_centers);
    for (uint32_t i = 0; i < num_purkinje_cells; i++) {
        fclose(((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->purkinje_files[i]);
    }
    free(((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data)->purkinje_files);

    // [PURKINJE COUPLED ACTIVATION TIMES]
    bool save_activation_time_map = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_activation_time_map, config, "save_activation_time");

    bool save_apd_map = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_apd_map, config, "save_apd");

    bool save_purkinje_velocity = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_purkinje_velocity, config, "save_purkinje_velocity");

    if(save_activation_time_map) {
        log_info("[!] Saving activation time maps !!!!\n");
        write_tissue_activation_time_maps(config, the_grid, output_dir, file_prefix, clip_with_plain, clip_with_bounds, binary, compress,
                                          compression_level, save_f);
        write_purkinje_activation_time_maps(config, the_grid, output_dir, file_prefix_purkinje, clip_with_plain, clip_with_bounds, binary,  compress,
                                            compression_level);
    }

    if(save_apd_map) {
        log_info("[!] Saving APD map !!!!\n");
        write_tissue_apd_map(config, the_grid, output_dir, file_prefix, clip_with_plain, clip_with_bounds, binary, compress, compression_level,
                             save_f);
        write_purkinje_apd_map(config, the_grid, output_dir, file_prefix_purkinje, clip_with_plain, clip_with_bounds, binary, compress,
                               compression_level);
    }

    if(save_purkinje_velocity) {
        log_info("[!] Calculating Purkinje propagation velocity !!!!\n");
        print_purkinje_propagation_velocity(config, the_grid);
    }

    free(config->persistent_data);
}

SAVE_MESH(save_multiple_cell_state_variables_purkinje_coupling_with_activation_times) {

    // [PURKINJE COUPLED ACTIVATION TIMES]
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config, "output_dir");

    float time_threshold = 10.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, time_threshold, config, "time_threshold");

    float tissue_activation_threshold = -30.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, tissue_activation_threshold, config, "activation_threshold_tissue");

    float tissue_apd_threshold = -83.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, tissue_apd_threshold, config, "apd_threshold_tissue");

    float purkinje_activation_threshold = -30.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, purkinje_activation_threshold, config, "activation_threshold_purkinje");

    float purkinje_apd_threshold = -83.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, purkinje_apd_threshold, config, "apd_threshold_purkinje");

// [TISSUE]
    real_cpu current_t = time_info->current_t;
    real_cpu last_t = time_info->final_t;
    real_cpu dt = time_info->dt;

    struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *persistent_data =
        (struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data;
    
    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    //OMP(parallel for)
    for (uint32_t i = 0; i < num_active_cells; i++) {
        if (ac[i]->active || (ac[i]->mesh_extra_info && (FIBROTIC(ac[i])) || BORDER_ZONE(ac[i]) )) {
            
            real_cpu center_x, center_y, center_z;
            real_cpu v;
            
            center_x = ac[i]->center.x;
            center_y = ac[i]->center.y;
            center_z = ac[i]->center.z;

            v = ac[i]->v;

            struct point_3d cell_coordinates;
            cell_coordinates.x = center_x;
            cell_coordinates.y = center_y;
            cell_coordinates.z = center_z;

            int n_activations = 0;
            float *apds_array = NULL;
            float *activation_times_array = NULL;

            if(ac[i]->active) {

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
    }

// [PURKINJE]
    uint32_t num_active_purkinje_cells = the_grid->purkinje->num_active_purkinje_cells;
    struct cell_node **ac_purkinje = the_grid->purkinje->purkinje_cells;

    //OMP(parallel for)
    for (uint32_t i = 0; i < num_active_purkinje_cells; i++) {
        if (ac_purkinje[i]->active || (ac_purkinje[i]->mesh_extra_info && (FIBROTIC(ac_purkinje[i]) || BORDER_ZONE(ac_purkinje[i]) )) ) {
            real_cpu center_x, center_y, center_z;
            real_cpu v;

            center_x = ac_purkinje[i]->center.x;
            center_y = ac_purkinje[i]->center.y;
            center_z = ac_purkinje[i]->center.z;

            v = ac_purkinje[i]->v;

            struct point_3d cell_coordinates;
            cell_coordinates.x = center_x;
            cell_coordinates.y = center_y;
            cell_coordinates.z = center_z;

            int n_activations = 0;
            float *apds_array = NULL;
            float *activation_times_array = NULL;

            if(ac_purkinje[i]->active) {

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
    }

    struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *params = ((struct save_multiple_cell_state_variables_purkinje_coupling_with_activation_times_data *)config->persistent_data);

    // [DOMAIN]
    for (uint32_t i = 0; i < params->num_tissue_cells; i++) {
        if(params->tissue_cell_sv_positions[i] == -1) {
            if(!the_grid->adaptive) {
                FOR_EACH_CELL(the_grid) {    
                    if (cell->active) {
                        if(cell->center.x == params->tissue_cell_centers[i*3] && cell->center.y == params->tissue_cell_centers[i*3+1] && cell->center.z == params->tissue_cell_centers[i*3+2]) {
                            params->tissue_cell_sv_positions[i] = cell->sv_position;
                            break;
                        }
                    }
                }
            }
        }
    }
    
    if(ode_solver->gpu) {
#ifdef COMPILE_CUDA
        for (uint32_t k = 0; k < params->num_tissue_cells; k++) {
            real *cell_sv;

            cell_sv = MALLOC_ARRAY_OF_TYPE(real, ode_solver->model_data.number_of_ode_equations);

            check_cuda_error(DPCT_CHECK_ERROR(dpct::dpct_memcpy(cell_sv, sizeof(real), ode_solver->sv + params->tissue_cell_sv_positions[k], ode_solver->pitch,
                                                                sizeof(real), ode_solver->model_data.number_of_ode_equations, cudaMemcpyDeviceToHost)));

            // Only 'time' and 'Vm'
            fprintf(params->tissue_files[k], "%g %g\n", time_info->current_t, cell_sv[0]);

            free(cell_sv);
        }
        
#endif
    } else {
        for (uint32_t k = 0; k < params->num_tissue_cells; k++) {
            real *cell_sv = &ode_solver->sv[params->tissue_cell_sv_positions[k] * ode_solver->model_data.number_of_ode_equations];

            // Only 'time' and 'Vm'
            fprintf(params->tissue_files[k], "%g %g\n", time_info->current_t, cell_sv[0]);
        }
    }

// [PURKINJE]
    for (uint32_t i = 0; i < params->num_purkinje_cells; i++) {
        if(params->purkinje_cell_sv_positions[i] == -1) {
            FOR_EACH_PURKINJE_CELL(the_grid) {
                real_cpu dist = calc_norm(cell->center.x,cell->center.y,cell->center.z,\
                                        params->purkinje_cell_centers[i*3],params->purkinje_cell_centers[i*3+1],params->purkinje_cell_centers[i*3+2]);
                if (dist < 1e-1) {
                    params->purkinje_cell_sv_positions[i] = cell->sv_position;
                    break;
                }
            }
        }
    }
    
    if(purkinje_ode_solver->gpu) {
#ifdef COMPILE_CUDA
        for (uint32_t k = 0; k < params->num_purkinje_cells; k++) {
            real *cell_sv;

            cell_sv = MALLOC_ARRAY_OF_TYPE(real, purkinje_ode_solver->model_data.number_of_ode_equations);

            check_cuda_error(DPCT_CHECK_ERROR(dpct::dpct_memcpy(cell_sv, sizeof(real), purkinje_ode_solver->sv + params->purkinje_cell_sv_positions[k],
                                                                purkinje_ode_solver->pitch, sizeof(real),
                                                                purkinje_ode_solver->model_data.number_of_ode_equations, cudaMemcpyDeviceToHost)));

            // Only 'time' and 'Vm'
            fprintf(params->purkinje_files[k], "%g %g\n", time_info->current_t, cell_sv[0]);

            free(cell_sv);
        }
        
#endif
    } else {
        for (uint32_t k = 0; k < params->num_tissue_cells; k++) {
            real *cell_sv = &purkinje_ode_solver->sv[params->purkinje_cell_sv_positions[k] * purkinje_ode_solver->model_data.number_of_ode_equations];

            // Only 'time' and 'Vm'
            fprintf(params->purkinje_files[k], "%g %g\n", time_info->current_t, cell_sv[0]);
        }
    }
}
