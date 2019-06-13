//
// Created by sachetto on 13/10/17.
//

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#include <unistd.h>
#include "../alg/grid/grid.h"
#include "../config/save_mesh_config.h"
#include "../common_types/common_types.h"
#include "../libraries_common/config_helpers.h"
#include "../monodomain/constants.h"
#include "../string/sds.h"
#include "../utils/utils.h"

#include "../vtk_utils/vtk_unstructured_grid.h"
#include "../vtk_utils/vtk_polydata_grid.h"
#include "../libraries_common/common_data_structures.h"

char *file_prefix;
bool binary = false;
bool clip_with_plain = false;
bool clip_with_bounds = false;
bool save_pvd = true;
bool compress = false;
int compression_level = 3;

static bool initialized = false;

static struct vtk_unstructured_grid *vtk_grid = NULL;

static struct vtk_polydata_grid *vtk_polydata = NULL;

void add_file_to_pvd(real_cpu current_t, const char *output_dir, const char *base_name);

static sds create_base_name(char *file_prefix, int iteration_count, char *extension) {
    return sdscatprintf(sdsempty(), "%s_it_%d.%s", file_prefix, iteration_count, extension);
}


SAVE_MESH(save_as_text_or_binary) {

    char *output_dir = config->out_dir_name;

    if(!initialized) {

        GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(file_prefix, config->config_data.config, "file_prefix");
        GET_PARAMETER_BINARY_VALUE_OR_USE_DEFAULT(binary, config->config_data.config, "binary");
        GET_PARAMETER_BINARY_VALUE_OR_USE_DEFAULT(clip_with_plain, config->config_data.config, "clip_with_plain");
        GET_PARAMETER_BINARY_VALUE_OR_USE_DEFAULT(clip_with_bounds, config->config_data.config, "clip_with_bounds");
        initialized = true;
    }

    real_cpu min_x = 0.0;
    real_cpu min_y = 0.0;
    real_cpu min_z = 0.0;
    real_cpu max_x = 0.0;
    real_cpu max_y = 0.0;
    real_cpu max_z = 0.0;

    float p0[3] = {0, 0, 0};
    float n[3] = {0, 0, 0};

    if(clip_with_plain) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, n[0], config->config_data.config, "normal_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, n[1], config->config_data.config, "normal_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, n[2], config->config_data.config, "normal_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, p0[0], config->config_data.config, "origin_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, p0[1], config->config_data.config, "origin_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, p0[2], config->config_data.config, "origin_z");
    }

    if(clip_with_bounds) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, min_x, config->config_data.config, "min_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, min_y, config->config_data.config, "min_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, min_z, config->config_data.config, "min_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, max_x, config->config_data.config, "max_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, max_y, config->config_data.config, "max_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, max_z, config->config_data.config, "max_z");
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

    sdsfree(base_name);
    sdsfree(tmp);

    struct cell_node *grid_cell = the_grid->first_cell;

    float center_x, center_y, center_z, dx, dy, dz;
    float v;

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
            dx = grid_cell->dx/2.0;
            dy = grid_cell->dy/2.0;
            dz = grid_cell->dz/2.0;

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
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[0], config->config_data.config, "origin_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[1], config->config_data.config, "origin_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[2], config->config_data.config, "origin_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[3], config->config_data.config, "normal_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[4], config->config_data.config, "normal_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[5], config->config_data.config, "normal_z");
    }

    if(clip_with_bounds) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[0], config->config_data.config, "min_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[1], config->config_data.config, "min_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[2], config->config_data.config, "min_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[3], config->config_data.config, "max_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[4], config->config_data.config, "max_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[5], config->config_data.config, "max_z");
    }

    // Write the transmembrane potential
    if  (scalar_name == 'v')
    {
        sds output_dir_with_file = sdsnew(output_dir);
        output_dir_with_file = sdscat(output_dir_with_file, "/");
        sds base_name = create_base_name(file_prefix, iteration_count, "vtk");

        //TODO: change this. We dont need the current_t here
        output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, current_t);

        new_vtk_unstructured_grid_from_alg_grid(&vtk_grid, the_grid, clip_with_plain, plain_coords, clip_with_bounds, bounds, !the_grid->adaptive,'v');

        save_vtk_unstructured_grid_as_legacy_vtk(vtk_grid, output_dir_with_file, binary);

        sdsfree(output_dir_with_file);
        sdsfree(base_name);
    }
    else if (scalar_name == 'a')
    {
        char *output_dir = config->out_dir_name;

        float plain_coords[6] = {0, 0, 0, 0, 0, 0};
        float bounds[6] = {0, 0, 0, 0, 0, 0};

        sds output_dir_with_file = sdsnew(output_dir);
        output_dir_with_file = sdscat(output_dir_with_file, "/activation-map.vtk");

        new_vtk_unstructured_grid_from_alg_grid(&vtk_grid, the_grid, clip_with_plain, plain_coords, clip_with_bounds, bounds, !the_grid->adaptive,'a');

        save_vtk_unstructured_grid_as_legacy_vtk(vtk_grid, output_dir_with_file, binary);

        sdsfree(output_dir_with_file);
    }
    else
    {
        fprintf(stderr,"[-] ERROR! Invalid scalar name!\n");
        exit(EXIT_FAILURE);
    }

    if(the_grid->adaptive)
        free_vtk_unstructured_grid(vtk_grid);

}

void add_file_to_pvd(real_cpu current_t, const char *output_dir, const char *base_name) {
    sds pvd_name = sdsnew(output_dir);
    pvd_name = sdscat(pvd_name, "/simulation_result.pvd");

    static FILE *pvd_file = NULL;
    static bool first_save_call = true;

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

    fprintf(pvd_file, "\n\t\t<DataSet timestep=\"%lf\" group=\"\" part=\"0\" file=\"%s\"/>\n", current_t, base_name);
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

        #ifndef COMPILE_ZLIB
        compress = false;
        #endif

        if(compress) binary = true;

        initialized = true;
    }
    float plain_coords[6] = {0, 0, 0, 0, 0, 0};
    float bounds[6] = {0, 0, 0, 0, 0, 0};

    if(clip_with_plain) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[0], config->config_data.config, "origin_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[1], config->config_data.config, "origin_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[2], config->config_data.config, "origin_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[3], config->config_data.config, "normal_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[4], config->config_data.config, "normal_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[5], config->config_data.config, "normal_z");
    }

    if(clip_with_bounds) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[0], config->config_data.config, "min_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[1], config->config_data.config, "min_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[2], config->config_data.config, "min_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[3], config->config_data.config, "max_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[4], config->config_data.config, "max_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[5], config->config_data.config, "max_z");
    }

    // Write transmembrane potential
    if (scalar_name == 'v')
    {
        sds output_dir_with_file = sdsnew(output_dir);
        output_dir_with_file = sdscat(output_dir_with_file, "/");
        sds base_name = create_base_name(file_prefix, iteration_count, "vtu");

        output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, current_t);

        if(save_pvd) 
        {
            add_file_to_pvd(current_t, output_dir, base_name);
        }

        new_vtk_unstructured_grid_from_alg_grid(&vtk_grid, the_grid, clip_with_plain, plain_coords, clip_with_bounds, bounds, !the_grid->adaptive,'v');

        if(compress) 
        {
            save_vtk_unstructured_grid_as_vtu_compressed(vtk_grid, output_dir_with_file, compression_level);
        }
        else 
        {
            save_vtk_unstructured_grid_as_vtu(vtk_grid, output_dir_with_file, binary);
        }

        if(the_grid->adaptive)
            free_vtk_unstructured_grid(vtk_grid);

        sdsfree(output_dir_with_file);
        sdsfree(base_name);
    }
    // Write activation map
    else if (scalar_name == 'a')
    {
        char *output_dir = config->out_dir_name;

        float plain_coords[6] = {0, 0, 0, 0, 0, 0};
        float bounds[6] = {0, 0, 0, 0, 0, 0};

        sds output_dir_with_file = sdsnew(output_dir);
        output_dir_with_file = sdscat(output_dir_with_file, "/activation-map.vtu");

        new_vtk_unstructured_grid_from_alg_grid(&vtk_grid, the_grid, clip_with_plain, plain_coords, clip_with_bounds, bounds, !the_grid->adaptive,'a');

        save_vtk_unstructured_grid_as_vtu(vtk_grid, output_dir_with_file, binary);

        sdsfree(output_dir_with_file);
    }
    else
    {
        fprintf(stderr,"[-] ERROR! Invalid scalar name!\n");
        exit(EXIT_FAILURE);
    }
}

SAVE_MESH(save_as_vtk_purkinje) {

    char *output_dir = config->out_dir_name;

    if(!initialized) 
    {
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
    sds base_name = create_base_name(file_prefix, iteration_count, "vtk");
    output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, current_t);

    new_vtk_polydata_grid_from_purkinje_grid(&vtk_polydata, the_grid,\
                                    clip_with_plain, plain_coords, clip_with_bounds, bounds,\
                                    !the_grid->adaptive);
    save_vtk_polydata_grid_as_legacy_vtk(vtk_polydata, output_dir_with_file, binary);

    if(the_grid->adaptive)
        free_vtk_polydata_grid(vtk_polydata);

    sdsfree(output_dir_with_file);
    sdsfree(base_name);

}

SAVE_MESH(save_as_vtp_purkinje) {

    char *output_dir = config->out_dir_name;

    if(!initialized) {
        GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(file_prefix, config->config_data.config, "file_prefix");
        GET_PARAMETER_BINARY_VALUE_OR_USE_DEFAULT(clip_with_plain, config->config_data.config, "clip_with_plain");
        GET_PARAMETER_BINARY_VALUE_OR_USE_DEFAULT(clip_with_bounds, config->config_data.config, "clip_with_bounds");
        GET_PARAMETER_BINARY_VALUE_OR_USE_DEFAULT(binary, config->config_data.config, "binary");
        GET_PARAMETER_BINARY_VALUE_OR_USE_DEFAULT(save_pvd, config->config_data.config, "save_pvd");
        GET_PARAMETER_BINARY_VALUE_OR_USE_DEFAULT(compress, config->config_data.config, "compress");
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, compression_level, config->config_data.config, "compression_level");

#ifndef DCOMPILE_ZLIB
        compress = false;
#endif
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
    sds base_name = create_base_name(file_prefix, iteration_count, "vtp");

    output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, current_t);

    if(save_pvd) 
    {
        add_file_to_pvd(current_t, output_dir, base_name);
    }

    new_vtk_polydata_grid_from_purkinje_grid(&vtk_polydata, the_grid,\
                                    clip_with_plain, plain_coords, clip_with_bounds, bounds,\
                                    !the_grid->adaptive);

    if(compress) 
    {
        save_vtk_polydata_grid_as_vtp_compressed(vtk_polydata, output_dir_with_file, compression_level);
    }
    else 
    {
        save_vtk_polydata_grid_as_vtp(vtk_polydata, output_dir_with_file, binary);
    }

    if(the_grid->adaptive)
        free_vtk_polydata_grid(vtk_polydata);

    sdsfree(output_dir_with_file);
    sdsfree(base_name);

}

SAVE_MESH(save_with_activation_times) {

<<<<<<< HEAD
    save_as_vtu(iteration_count, current_t, last_t, config, the_grid, 'v');
    //save_as_text_or_binary(iteration_count, current_t, last_t, config, the_grid,'v');

=======
    save_as_text_or_binary(iteration_count, current_t, last_t, dt, config, the_grid);
>>>>>>> rsachetto-master

    float time_threshold = 0.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, time_threshold, config->config_data.config, "time_threshold");

    char *output_dir = config->out_dir_name;

    sds output_dir_with_file = sdsnew(output_dir);
    output_dir_with_file = sdscat(output_dir_with_file, "/");
    sds base_name = create_base_name("activation_info", 0, "txt");
    output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, current_t);

    static struct point_hash_entry *last_time_v = NULL;

    static struct point_hash_entry*num_activations = NULL;

    static struct point_voidp_hash_entry *activation_times = NULL;

    if(last_time_v == NULL) {
        hmdefault(last_time_v, -100.0);
    }

    if(num_activations == NULL) {
        hmdefault(num_activations, 0);
    }

    if(activation_times == NULL) {
        hmdefault(activation_times, NULL);
    }

    struct cell_node *grid_cell = the_grid->first_cell;

    float center_x, center_y, center_z, dx, dy, dz;
    float v;

    FILE *act_file = fopen(output_dir_with_file, "w");

    fprintf(act_file, "%d\n", (last_t-current_t) <= dt ); //rounding errors

    while(grid_cell != 0) {

        if( grid_cell->active || ( grid_cell->mesh_extra_info && ( FIBROTIC(grid_cell) || BORDER_ZONE(grid_cell) ) ) ) {
            center_x = grid_cell->center_x;
            center_y = grid_cell->center_y;
            center_z = grid_cell->center_z;

            v = grid_cell->v;

            struct point_3d p;
            p.x = center_x;
            p.y = center_y;
            p.z = center_z;

            dx = grid_cell->dx / 2.0;
            dy = grid_cell->dy / 2.0;
            dz = grid_cell->dz / 2.0;

            fprintf(act_file, "%g,%g,%g,%g,%g,%g,%d,%d,%d ", center_x, center_y, center_z, dx, dy, dz, grid_cell->active, FIBROTIC(grid_cell), BORDER_ZONE(grid_cell));

            float last_v = hmget(last_time_v, p);

            int n_activations = (int) hmget(num_activations, p);
            float *activation_times_array = (float *) hmget(activation_times, p);

            int act_times_len = arrlen(activation_times_array);

            if (current_t == 0.0f) {
                hmput(last_time_v, p, v);
            } else {
                if ((last_v < 0.0f) && (v >= 0.0f)) {

                    if (act_times_len == 0) {
                        n_activations++;
                        hmput(num_activations, p, n_activations);
                                arrput(activation_times_array, current_t);
                        hmput(activation_times, p, activation_times_array);
                    } else {
                        float last_act_time = activation_times_array[act_times_len - 1];
                        if (current_t - last_act_time > time_threshold) {
                            n_activations++;
                            hmput(num_activations, p, n_activations);
                            arrput(activation_times_array, current_t);
                            hmput(activation_times, p, activation_times_array);
                        }

                    }
                }
                hmput(last_time_v, p, v);

            }

            fprintf(act_file, "%d [ ", n_activations);

            for (int i = 0; i < arrlen(activation_times_array); i++) {
                fprintf(act_file, "%lf ", activation_times_array[i]);
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
