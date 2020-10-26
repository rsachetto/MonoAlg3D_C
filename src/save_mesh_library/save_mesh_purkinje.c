//
// Created by bergolho on 30/09/20.
//

#include <stdbool.h>
#include <stdlib.h>

#include "../3dparty/sds/sds.h"
#include "../alg/grid/grid.h"
#include "../config/save_mesh_config.h"
#include "../utils/utils.h"

#include "../vtk_utils/vtk_unstructured_grid.h"
#include "../vtk_utils/vtk_polydata_grid.h"
#include "../libraries_common/common_data_structures.h"

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

static bool initialized = false;


INIT_SAVE_MESH(init_save_as_vtk_or_vtp) {
    config->persistent_data = malloc(sizeof(struct save_as_vtp_persistent_data));
    ((struct save_as_vtp_persistent_data *) config->persistent_data)->grid = NULL;
    ((struct save_as_vtp_persistent_data *) config->persistent_data)->first_save_call = true;
}

END_SAVE_MESH(end_save_as_vtk_or_vtp) {
    free_vtk_polydata_grid(((struct save_as_vtp_persistent_data *) config->persistent_data)->grid);
    free(config->persistent_data);
}

SAVE_MESH(save_as_vtp_purkinje) {

    int iteration_count = time_info->iteration;

    if(((struct save_as_vtp_persistent_data *) config->persistent_data)->first_save_call) {
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config->config_data, "output_dir");
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(file_prefix, config->config_data, "file_prefix");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_plain, config->config_data, "clip_with_plain");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_bounds, config->config_data, "clip_with_bounds");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(binary, config->config_data, "binary");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_pvd, config->config_data, "save_pvd");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(compress, config->config_data, "compress");
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, compression_level, config->config_data, "compression_level");

        if(compress) binary = true;

        if(!save_pvd) {
            ((struct save_as_vtp_persistent_data *) config->persistent_data)->first_save_call = false;
        }
    }

    float plain_coords[6] = {0, 0, 0, 0, 0, 0};
    float bounds[6] = {0, 0, 0, 0, 0, 0};

    if(clip_with_plain) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[0], config->config_data, "origin_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[1], config->config_data, "origin_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[2], config->config_data, "origin_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[3], config->config_data, "normal_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[4], config->config_data, "normal_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[5], config->config_data, "normal_z");
    }

    if(clip_with_bounds) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[0], config->config_data, "min_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[1], config->config_data, "min_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[2], config->config_data, "min_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[3], config->config_data, "max_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[4], config->config_data, "max_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[5], config->config_data, "max_z");
    }


    sds output_dir_with_file = sdsnew(output_dir);
    output_dir_with_file = sdscat(output_dir_with_file, "/");
    sds base_name = create_base_name(file_prefix, iteration_count, "vtp");

    real_cpu current_t = time_info->current_t;

    output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, current_t);

    if(save_pvd) {
        add_file_to_pvd(current_t, output_dir, base_name, ((struct save_as_vtp_persistent_data *) config->persistent_data)->first_save_call);
        ((struct save_as_vtp_persistent_data *) config->persistent_data)->first_save_call = false;
    }

    bool read_only_data = ((struct save_as_vtp_persistent_data *) config->persistent_data)->grid != NULL;
    new_vtk_polydata_grid_from_purkinje_grid(&((struct save_as_vtp_persistent_data *) config->persistent_data)->grid, the_grid->purkinje, clip_with_plain, plain_coords, clip_with_bounds, bounds, read_only_data);

    if(compress) {
        save_vtk_polydata_grid_as_vtp_compressed(((struct save_as_vtp_persistent_data *) config->persistent_data)->grid, output_dir_with_file, compression_level);
    }
    else {
        save_vtk_polydata_grid_as_vtp(((struct save_as_vtp_persistent_data *) config->persistent_data)->grid, output_dir_with_file, binary);
    }

    sdsfree(output_dir_with_file);
    sdsfree(base_name);
}

SAVE_MESH(save_as_vtk_purkinje) {

    char *output_dir;
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config->config_data, "output_dir");

    if(!initialized){

        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(file_prefix, config->config_data, "file_prefix");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_plain, config->config_data, "clip_with_plain");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_bounds, config->config_data, "clip_with_bounds");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(binary, config->config_data, "binary");
        initialized = true;
    }
    float plain_coords[6] = {0, 0, 0, 0, 0, 0};
    float bounds[6] = {0, 0, 0, 0, 0, 0};

    if(clip_with_plain) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[0], config->config_data, "origin_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[1], config->config_data, "origin_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[2], config->config_data, "origin_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[3], config->config_data, "normal_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[3], config->config_data, "normal_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[4], config->config_data, "normal_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[5], config->config_data, "normal_z");
    }

    if(clip_with_bounds) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[0], config->config_data, "min_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[1], config->config_data, "min_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[2], config->config_data, "min_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[3], config->config_data, "max_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[4], config->config_data, "max_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[5], config->config_data, "max_z");
    }

    sds output_dir_with_file = sdsnew(output_dir);
    output_dir_with_file = sdscat(output_dir_with_file, "/");
    sds base_name = create_base_name(file_prefix, time_info->iteration, "vtk");
    output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, time_info->current_t);

    bool read_only_data = ((struct save_as_vtp_persistent_data *) config->persistent_data)->grid != NULL;
    new_vtk_polydata_grid_from_purkinje_grid(&((struct save_as_vtp_persistent_data *) config->persistent_data)->grid, the_grid->purkinje, clip_with_plain, plain_coords, clip_with_bounds, bounds, read_only_data);
    
    save_vtk_polydata_grid_as_legacy_vtk(((struct save_as_vtp_persistent_data *) config->persistent_data)->grid, output_dir_with_file, binary);

    sdsfree(output_dir_with_file);
    sdsfree(base_name);

}

INIT_SAVE_MESH(init_save_tissue_as_vtk_or_vtu_purkinje_as_vtp) {
    config->persistent_data = malloc(sizeof(struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data));
    ((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *) config->persistent_data)->grid = NULL;
    ((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *) config->persistent_data)->grid_purkinje = NULL;
    ((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *) config->persistent_data)->first_save_call = true;
}

END_SAVE_MESH(end_save_tissue_as_vtk_or_vtu_purkinje_as_vtp) {
    free_vtk_polydata_grid(((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *) config->persistent_data)->grid_purkinje);
    free_vtk_unstructured_grid(((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *) config->persistent_data)->grid);
    free(config->persistent_data);
}

SAVE_MESH(save_tissue_as_vtu_purkinje_as_vtp) {

// [TISSUE]
    int iteration_count = time_info->iteration;

    if(((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *) config->persistent_data)->first_save_call) {
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config->config_data, "output_dir");
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(file_prefix, config->config_data, "file_prefix");
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(file_prefix_purkinje, config->config_data, "file_prefix_purkinje");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_plain, config->config_data, "clip_with_plain");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_bounds, config->config_data, "clip_with_bounds");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(binary, config->config_data, "binary");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_pvd, config->config_data, "save_pvd");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(compress, config->config_data, "compress");
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, compression_level, config->config_data, "compression_level");

        if(compress) binary = true;

        if(!save_pvd) {
            ((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *) config->persistent_data)->first_save_call = false;
        }
    }

    float plain_coords[6] = {0, 0, 0, 0, 0, 0};
    float bounds[6] = {0, 0, 0, 0, 0, 0};

    if(clip_with_plain) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[0], config->config_data, "origin_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[1], config->config_data, "origin_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[2], config->config_data, "origin_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[3], config->config_data, "normal_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[4], config->config_data, "normal_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[5], config->config_data, "normal_z");
    }

    if(clip_with_bounds) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[0], config->config_data, "min_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[1], config->config_data, "min_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[2], config->config_data, "min_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[3], config->config_data, "max_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[4], config->config_data, "max_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[5], config->config_data, "max_z");
    }


    sds output_dir_with_file = sdsnew(output_dir);
    output_dir_with_file = sdscat(output_dir_with_file, "/");
    sds base_name = create_base_name(file_prefix, iteration_count, "vtu");

    real_cpu current_t = time_info->current_t;

    output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, current_t);

    if(save_pvd) {
        add_file_to_pvd(current_t, output_dir, base_name, ((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *) config->persistent_data)->first_save_call);
        ((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *) config->persistent_data)->first_save_call = false;
    }

    bool read_only_data = ((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *) config->persistent_data)->grid != NULL;
    new_vtk_unstructured_grid_from_alg_grid(&((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *) config->persistent_data)->grid, the_grid, clip_with_plain, plain_coords, clip_with_bounds, bounds, read_only_data, save_f);

    if(compress) {
        save_vtk_unstructured_grid_as_vtu_compressed(((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *) config->persistent_data)->grid, output_dir_with_file, compression_level);
    }
    else {
        save_vtk_unstructured_grid_as_vtu(((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *) config->persistent_data)->grid, output_dir_with_file, binary);
    }

    if(the_grid->adaptive) {
        free_vtk_unstructured_grid(((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *) config->persistent_data)->grid);
        ((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *) config->persistent_data)->grid = NULL;
    }

// [PURKINJE]

    output_dir_with_file = sdsnew(output_dir);
    output_dir_with_file = sdscat(output_dir_with_file, "/");
    base_name = create_base_name(file_prefix_purkinje, iteration_count, "vtp");

    output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, current_t);

    if(save_pvd) {
        add_file_to_pvd(current_t, output_dir, base_name, ((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *) config->persistent_data)->first_save_call);
        ((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *) config->persistent_data)->first_save_call = false;
    }

    read_only_data = ((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *) config->persistent_data)->grid_purkinje != NULL;
    new_vtk_polydata_grid_from_purkinje_grid(&((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *) config->persistent_data)->grid_purkinje, the_grid->purkinje, clip_with_plain, plain_coords, clip_with_bounds, bounds, read_only_data);

    if(compress) {
        save_vtk_polydata_grid_as_vtp_compressed(((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *) config->persistent_data)->grid_purkinje, output_dir_with_file, compression_level);
    }
    else {
        save_vtk_polydata_grid_as_vtp(((struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data *) config->persistent_data)->grid_purkinje, output_dir_with_file, binary);
    }

    sdsfree(output_dir_with_file);
    sdsfree(base_name);

}

INIT_SAVE_MESH(init_save_purkinje_with_activation_times) {

    config->persistent_data = calloc(1, sizeof(struct save_coupling_with_activation_times_persistent_data));

    hmdefault(((struct save_coupling_with_activation_times_persistent_data*)config->persistent_data)->purkinje_cell_was_active, 0.0);
    hmdefault(((struct save_coupling_with_activation_times_persistent_data*)config->persistent_data)->purkinje_last_time_v, -100.0);
    hmdefault(((struct save_coupling_with_activation_times_persistent_data*)config->persistent_data)->purkinje_num_activations, 0);
    hmdefault(((struct save_coupling_with_activation_times_persistent_data*)config->persistent_data)->purkinje_activation_times, NULL);
    hmdefault(((struct save_coupling_with_activation_times_persistent_data*)config->persistent_data)->purkinje_apds, NULL);
    ((struct save_coupling_with_activation_times_persistent_data*) config->persistent_data)->purkinje_grid = NULL;

    ((struct save_coupling_with_activation_times_persistent_data*)config->persistent_data)->first_save_call = true;
}

END_SAVE_MESH(end_save_purkinje_with_activation_times) {

    bool save_activation_time_map = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_activation_time_map, config->config_data, "save_activation_time");

    bool save_apd_map = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_apd_map, config->config_data, "save_apd");

    if (save_activation_time_map) {
        log_to_stderr_and_file("[!] Saving activation time maps !!!!\n");
        write_purkinje_activation_time_maps(config,the_grid,output_dir,\
                                    file_prefix_purkinje,clip_with_plain,clip_with_bounds,binary,save_pvd,compress,compression_level);
    }
    
    if (save_apd_map) {
        log_to_stderr_and_file("[!] Saving APD map !!!!\n");
        write_purkinje_apd_map(config,the_grid,output_dir,\
                                    file_prefix_purkinje,clip_with_plain,clip_with_bounds,binary,save_pvd,compress,compression_level);
    } 
  
    free(config->persistent_data);

}

SAVE_MESH (save_purkinje_with_activation_times) {

    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config->config_data, "output_dir");

    float time_threshold = 10.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, time_threshold, config->config_data, "time_threshold");
    
    float purkinje_activation_threshold = -30.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, purkinje_activation_threshold, config->config_data, "activation_threshold_purkinje");

    float purkinje_apd_threshold = -83.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, purkinje_apd_threshold, config->config_data, "apd_threshold_purkinje");

    calculate_purkinje_activation_time_and_apd(time_info,config,the_grid,time_threshold,purkinje_activation_threshold,purkinje_apd_threshold);

}

INIT_SAVE_MESH(init_save_purkinje_coupling_with_activation_times) {

    config->persistent_data = calloc(1, sizeof(struct save_coupling_with_activation_times_persistent_data));
    hmdefault(((struct save_coupling_with_activation_times_persistent_data*)config->persistent_data)->tissue_cell_was_active, 0.0);
    hmdefault(((struct save_coupling_with_activation_times_persistent_data*)config->persistent_data)->tissue_last_time_v, -100.0);
    hmdefault(((struct save_coupling_with_activation_times_persistent_data*)config->persistent_data)->tissue_num_activations, 0);
    hmdefault(((struct save_coupling_with_activation_times_persistent_data*)config->persistent_data)->tissue_activation_times, NULL);
    hmdefault(((struct save_coupling_with_activation_times_persistent_data*)config->persistent_data)->tissue_apds, NULL);
    ((struct save_coupling_with_activation_times_persistent_data*) config->persistent_data)->tissue_grid = NULL;

    hmdefault(((struct save_coupling_with_activation_times_persistent_data*)config->persistent_data)->purkinje_cell_was_active, 0.0);
    hmdefault(((struct save_coupling_with_activation_times_persistent_data*)config->persistent_data)->purkinje_last_time_v, -100.0);
    hmdefault(((struct save_coupling_with_activation_times_persistent_data*)config->persistent_data)->purkinje_num_activations, 0);
    hmdefault(((struct save_coupling_with_activation_times_persistent_data*)config->persistent_data)->purkinje_activation_times, NULL);
    hmdefault(((struct save_coupling_with_activation_times_persistent_data*)config->persistent_data)->purkinje_apds, NULL);
    ((struct save_coupling_with_activation_times_persistent_data*) config->persistent_data)->purkinje_grid = NULL;

    ((struct save_coupling_with_activation_times_persistent_data*)config->persistent_data)->first_save_call = true;
}

END_SAVE_MESH(end_save_purkinje_coupling_with_activation_times) {

    bool save_activation_time_map = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_activation_time_map, config->config_data, "save_activation_time");

    bool save_apd_map = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_apd_map, config->config_data, "save_apd");

    bool save_purkinje_velocity = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_purkinje_velocity, config->config_data, "save_purkinje_velocity");

    if (save_activation_time_map) {
        log_to_stderr_and_file("[!] Saving activation time maps !!!!\n");
        write_tissue_activation_time_maps(config,the_grid,output_dir,file_prefix,clip_with_plain,clip_with_bounds,binary,save_pvd,compress,compression_level,save_f);
        write_purkinje_activation_time_maps(config,the_grid,output_dir,file_prefix_purkinje,clip_with_plain,clip_with_bounds,binary,save_pvd,compress,compression_level);
    }
    
    if (save_apd_map) {
        log_to_stderr_and_file("[!] Saving APD map !!!!\n");
        write_tissue_apd_map(config,the_grid,output_dir,file_prefix,clip_with_plain,clip_with_bounds,binary,save_pvd,compress,compression_level,save_f);
        write_purkinje_apd_map(config,the_grid,output_dir,file_prefix_purkinje,clip_with_plain,clip_with_bounds,binary,save_pvd,compress,compression_level);
    } 

    if (save_purkinje_velocity) {
        log_to_stderr_and_file("[!] Calculating Purkinje propagation velocity !!!!\n");
        print_purkinje_propagation_velocity(config,the_grid);
    }
  
    free(config->persistent_data);

}

SAVE_MESH (save_purkinje_coupling_with_activation_times) {

    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config->config_data, "output_dir");

    float time_threshold = 10.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, time_threshold, config->config_data, "time_threshold");

    float tissue_activation_threshold = -30.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, tissue_activation_threshold, config->config_data, "activation_threshold_tissue");

    float tissue_apd_threshold = -83.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, tissue_apd_threshold, config->config_data, "apd_threshold_tissue");
    
    float purkinje_activation_threshold = -30.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, purkinje_activation_threshold, config->config_data, "activation_threshold_purkinje");

    float purkinje_apd_threshold = -83.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, purkinje_apd_threshold, config->config_data, "apd_threshold_purkinje");

// [TISSUE]
    calculate_tissue_activation_time_and_apd(time_info,config,the_grid,time_threshold,tissue_activation_threshold,tissue_apd_threshold);

// [PURKINJE]
    calculate_purkinje_activation_time_and_apd(time_info,config,the_grid,time_threshold,purkinje_activation_threshold,purkinje_apd_threshold);

}

INIT_SAVE_MESH(init_save_one_cell_state_variables) {
    config->persistent_data = malloc(sizeof(struct save_one_cell_state_variables_persistent_data));
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR( ((struct save_one_cell_state_variables_persistent_data *) config->persistent_data)->file_name, config->config_data, "file_name");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, ((struct save_one_cell_state_variables_persistent_data *) config->persistent_data)->cell_center_x, config->config_data, "cell_center_x");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, ((struct save_one_cell_state_variables_persistent_data *) config->persistent_data)->cell_center_y, config->config_data, "cell_center_y");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, ((struct save_one_cell_state_variables_persistent_data *) config->persistent_data)->cell_center_z, config->config_data, "cell_center_z");

    ((struct save_one_cell_state_variables_persistent_data *) config->persistent_data)->file = fopen(((struct save_one_cell_state_variables_persistent_data *) config->persistent_data)->file_name, "w");
    ((struct save_one_cell_state_variables_persistent_data *) config->persistent_data)->cell_sv_position = -1;
}

END_SAVE_MESH(end_save_one_cell_state_variables) {
    free(((struct save_one_cell_state_variables_persistent_data *) config->persistent_data)->file_name);
    fclose(((struct save_one_cell_state_variables_persistent_data *) config->persistent_data)->file);
    free(config->persistent_data);
}

SAVE_MESH(save_one_cell_state_variables) {

    struct save_one_cell_state_variables_persistent_data *params = ((struct save_one_cell_state_variables_persistent_data *) config->persistent_data);

    if(params->cell_sv_position == -1) {
        if(!the_grid->adaptive) {
            FOR_EACH_CELL(the_grid) {
                if(cell->center.x == params->cell_center_x && cell->center.y == params->cell_center_y && cell->center.z == params->cell_center_z) {
                    params->cell_sv_position = cell->sv_position;
                    printf("%d\n", params->cell_sv_position);
                    break;
                }
            }
        }
    }

    if(ode_solver->gpu) {
#ifdef COMPILE_CUDA
            int num_odes = ode_solver->model_data.number_of_ode_equations;
            real *cell_sv;

            cell_sv = (real *)malloc(sizeof(real) * num_odes);

            check_cuda_error(cudaMemcpy2D(cell_sv, sizeof(real), ode_solver->sv + params->cell_sv_position,
                                           ode_solver->pitch, sizeof(real),
                                           ode_solver->model_data.number_of_ode_equations, cudaMemcpyDeviceToHost));

            fprintf(params->file, "%lf ", time_info->current_t);
            for (int i = 0; i < num_odes; i++) {
                fprintf(params->file,"%lf ",cell_sv[i]);
            }
            fprintf(params->file, "\n");

            free(cell_sv);
#endif
    }
    else {

        int num_odes = ode_solver->model_data.number_of_ode_equations;
        real *cell_sv =  &ode_solver->sv[params->cell_sv_position * num_odes];

        fprintf(params->file, "%lf ", time_info->current_t);
        for (int i = 0; i < num_odes; i++) {
            fprintf(params->file,"%lf ",cell_sv[i]);
        }
        fprintf(params->file, "\n");
    }

}
