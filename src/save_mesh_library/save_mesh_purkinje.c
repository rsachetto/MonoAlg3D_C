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

static bool initialized = false;

void add_file_to_pvd(real_cpu current_t, const char *output_dir, const char *base_name, bool first_save_call);

static sds create_base_name(char *f_prefix, int iteration_count, char *extension) {
    return sdscatprintf(sdsempty(), "%s_it_%d.%s", f_prefix, iteration_count, extension);
}

static inline void write_pvd_header(FILE *pvd_file) {
    fprintf(pvd_file, "<VTKFile type=\"Collection\" version=\"0.1\" compressor=\"vtkZLibDataCompressor\">\n");
    fprintf(pvd_file, "\t<Collection>\n");
    fprintf(pvd_file, "\t</Collection>\n");
    fprintf(pvd_file, "</VTKFile>");
}

void add_file_to_pvd(real_cpu current_t, const char *output_dir, const char *base_name, bool first_call) {

    sds pvd_name = sdsnew(output_dir);
    pvd_name = sdscat(pvd_name, "/simulation_result.pvd");

    static FILE *pvd_file = NULL;
    pvd_file = fopen(pvd_name, "r+");

    if(!pvd_file) {
        pvd_file = fopen(pvd_name, "w");
        write_pvd_header(pvd_file);
    }
    else {
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

struct save_as_vtp_persistent_data {
    struct vtk_polydata_grid *grid;
    bool first_save_call;
};

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

static char *file_prefix_purkinje;

struct save_tissue_as_vtu_purkinje_as_vtp_persistent_data {
    struct vtk_unstructured_grid *grid;
    struct vtk_polydata_grid *grid_purkinje;
    bool first_save_call;
};

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

