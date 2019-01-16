//
// Created by sachetto on 13/10/17.
//

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>

#include "../alg/grid/grid.h"
#include "../config/save_mesh_config.h"
#include "../hash/point_hash.h"
#include "../libraries_common/config_helpers.h"
#include "../monodomain/constants.h"
#include "../string/sds.h"
#include "../utils/utils.h"
#include "vtk_unstructured_grid.h"

char *file_prefix;
bool binary = false;
bool clip_with_plain = false;
bool clip_with_bounds = false;
bool save_pvd = false;
bool compress = false;
int compression_level = 3;
static FILE *pvd_file = NULL;

static bool initialized = false;
static bool first_save_call = true;
static int count = 0;

static struct vtk_unstructured_grid *vtk_grid = NULL;

void add_file_to_pvd(double current_dt, const char *output_dir, const char *base_name);

SAVE_MESH(save_as_text_or_binary) {

    char *output_dir = config->out_dir_name;

    if(!initialized) {

        GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(file_prefix, config->config_data.config, "file_prefix");
        GET_PARAMETER_BINARY_VALUE_OR_USE_DEFAULT(binary, config->config_data.config, "binary");
        GET_PARAMETER_BINARY_VALUE_OR_USE_DEFAULT(clip_with_plain, config->config_data.config, "clip_with_plain");
        GET_PARAMETER_BINARY_VALUE_OR_USE_DEFAULT(clip_with_bounds, config->config_data.config, "clip_with_bounds");
        initialized = true;
    }

    real min_x = 0.0;
    real min_y = 0.0;
    real min_z = 0.0;
    real max_x = 0.0;
    real max_y = 0.0;
    real max_z = 0.0;

    real p0[3] = {0, 0, 0};
    real n[3] = {0, 0, 0};

    if(clip_with_plain) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, n[0], config->config_data.config, "normal_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, n[1], config->config_data.config, "normal_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, n[2], config->config_data.config, "normal_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, p0[0], config->config_data.config, "origin_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, p0[1], config->config_data.config, "origin_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, p0[2], config->config_data.config, "origin_z");
    }

    if(clip_with_bounds) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, min_x, config->config_data.config, "min_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, min_y, config->config_data.config, "min_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, min_z, config->config_data.config, "min_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, max_x, config->config_data.config, "max_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, max_y, config->config_data.config, "max_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, max_z, config->config_data.config, "max_z");
    }

    real l = sqrtf(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    real A = n[0] / l;
    real B = n[1] / l;
    real C = n[2] / l;
    real D = -(n[0] * p0[0] + n[1] * p0[1] + n[2] * p0[2]);

    double side;

    sds tmp = sdsnew(output_dir);
    tmp = sdscatprintf(tmp, "/%s_it_%lf.vtk", file_prefix, current_dt);

    FILE *output_file = fopen(tmp, "w");

    struct cell_node *grid_cell = the_grid->first_cell;

    double center_x, center_y, center_z, dx, dy, dz;
    double v;

    while(grid_cell != 0) {

        if(grid_cell->active) {

            center_x = grid_cell->center_x;
            center_y = grid_cell->center_y;
            center_z = grid_cell->center_z;

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
            dx = grid_cell->dx;
            dy = grid_cell->dy;
            dz = grid_cell->dz;

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
        }
        grid_cell = grid_cell->next;
    }

    fclose(output_file);
}

SAVE_MESH(save_as_vtk) {

    char *output_dir = config->out_dir_name;

    if(!initialized) {
        GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(file_prefix, config->config_data.config, "file_prefix");
        GET_PARAMETER_BINARY_VALUE_OR_USE_DEFAULT(clip_with_plain, config->config_data.config, "clip_with_plain");
        GET_PARAMETER_BINARY_VALUE_OR_USE_DEFAULT(clip_with_bounds, config->config_data.config, "clip_with_bounds");
        GET_PARAMETER_BINARY_VALUE_OR_USE_DEFAULT(binary, config->config_data.config, "binary");
        initialized = true;
    }
    float plain_coords[6] = {0, 0, 0, 0, 0, 0};
    float bounds[6] = {0, 0, 0, 0, 0, 0};

    if(clip_with_plain) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[0], config->config_data.config, "origin_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[1], config->config_data.config, "origin_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[2], config->config_data.config, "origin_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[3], config->config_data.config, "normal_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[4], config->config_data.config, "normal_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[5], config->config_data.config, "normal_z");
    }

    if(clip_with_bounds) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[0], config->config_data.config, "min_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[1], config->config_data.config, "min_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[2], config->config_data.config, "min_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[3], config->config_data.config, "max_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[4], config->config_data.config, "max_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[5], config->config_data.config, "max_z");
    }

    sds output_dir_with_file = sdsnew(output_dir);
    output_dir_with_file = sdscat(output_dir_with_file, "/");
    sds base_name = sdscatprintf(sdsempty(), "%s_it_%d_ms.vtu", file_prefix, count);
    count++;
    output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, current_dt);

    new_vtk_unstructured_grid_from_alg_grid(&vtk_grid, the_grid, clip_with_plain, plain_coords, clip_with_bounds, bounds, !the_grid->adaptive);

    save_vtk_unstructured_grid_as_legacy_vtk(vtk_grid, output_dir_with_file, binary);

    if(the_grid->adaptive)
        free_vtk_unstructured_grid(vtk_grid);

    sdsfree(output_dir_with_file);
}

void add_file_to_pvd(double current_dt, const char *output_dir, const char *base_name) {
    sds pvd_name = sdsnew(output_dir);
    pvd_name = sdscat(pvd_name, "/simulation_result.pvd");

    if(first_save_call) {
        pvd_file = fopen(pvd_name, "w");
        fprintf(pvd_file, "<VTKFile type=\"Collection\" version=\"0.1\" compressor=\"vtkZLibDataCompressor\">\n");
        fprintf(pvd_file, "\t<Collection>\n");
        fprintf(pvd_file, "\t</Collection>\n");
        fprintf(pvd_file, "</VTKFile>");
        first_save_call = false;

    } else {
        pvd_file = fopen(pvd_name, "r+");
    }

    sdsfree(pvd_name);

    fseek(pvd_file, -26, SEEK_END);

    fprintf(pvd_file, "\n\t\t<DataSet timestep=\"%lf\" group=\"\" part=\"0\" file=\"%s\"/>\n", current_dt, base_name);
    fprintf(pvd_file, "\t</Collection>\n");
    fprintf(pvd_file, "</VTKFile>");
    fclose(pvd_file);
}

SAVE_MESH(save_as_vtu) {

    char *output_dir = config->out_dir_name;

    if(!initialized) {
        GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(file_prefix, config->config_data.config, "file_prefix");
        GET_PARAMETER_BINARY_VALUE_OR_USE_DEFAULT(clip_with_plain, config->config_data.config, "clip_with_plain");
        GET_PARAMETER_BINARY_VALUE_OR_USE_DEFAULT(clip_with_bounds, config->config_data.config, "clip_with_bounds");
        GET_PARAMETER_BINARY_VALUE_OR_USE_DEFAULT(binary, config->config_data.config, "binary");
        GET_PARAMETER_BINARY_VALUE_OR_USE_DEFAULT(save_pvd, config->config_data.config, "save_pvd");
        GET_PARAMETER_BINARY_VALUE_OR_USE_DEFAULT(compress, config->config_data.config, "compress");
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, compression_level, config->config_data.config, "compression_level");
        if(compress) binary = true;

        initialized = true;
    }
    float plain_coords[6] = {0, 0, 0, 0, 0, 0};
    float bounds[6] = {0, 0, 0, 0, 0, 0};

    if(clip_with_plain) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[0], config->config_data.config, "origin_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[1], config->config_data.config, "origin_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[2], config->config_data.config, "origin_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[3], config->config_data.config, "normal_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[4], config->config_data.config, "normal_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, plain_coords[5], config->config_data.config, "normal_z");
    }

    if(clip_with_bounds) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[0], config->config_data.config, "min_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[1], config->config_data.config, "min_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[2], config->config_data.config, "min_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[3], config->config_data.config, "max_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[4], config->config_data.config, "max_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, bounds[5], config->config_data.config, "max_z");
    }


    sds output_dir_with_file = sdsnew(output_dir);
    output_dir_with_file = sdscat(output_dir_with_file, "/");
    sds base_name = sdscatprintf(sdsempty(), "%s_it_%d.vtu", file_prefix, count);
    count++;
    output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, current_dt);

    if(save_pvd) {
        add_file_to_pvd(current_dt, output_dir, base_name);
    }

    new_vtk_unstructured_grid_from_alg_grid(&vtk_grid, the_grid, clip_with_plain, plain_coords, clip_with_bounds, bounds, !the_grid->adaptive);

    if(compress) {
        save_vtk_unstructured_grid_as_vtu_compressed(vtk_grid, output_dir_with_file, compression_level);
    }
    else {
        save_vtk_unstructured_grid_as_vtu(vtk_grid, output_dir_with_file, binary);
    }

    if(the_grid->adaptive)
        free_vtk_unstructured_grid(vtk_grid);

    sdsfree(output_dir_with_file);
    sdsfree(base_name);

}