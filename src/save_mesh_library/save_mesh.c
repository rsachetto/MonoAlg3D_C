//
// Created by sachetto on 13/10/17.
//

#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>

#include "../utils/utils.h"
#include "../monodomain/constants.h"
#include "../alg/grid/grid.h"
#include "../config/assembly_matrix_config.h"
#include "../libraries_common/config_helpers.h"
#include "../string/sds.h"
#include "../hash/point_hash.h"

static int count = 0;

char * file_prefix;
bool binary = false;
bool clip_with_plain = false;
bool clip_with_bounds = false;

bool initialized = false;
SAVE_MESH(save_as_text_or_binary) {

    if(!initialized) {

        GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(file_prefix, config->config_data.config, "file_prefix");
        GET_PARAMETER_BINARY_VALUE (binary, config->config_data.config, "binary");
        GET_PARAMETER_BINARY_VALUE (clip_with_plain, config->config_data.config, "clip_with_plain");
        GET_PARAMETER_BINARY_VALUE (clip_with_bounds, config->config_data.config, "clip_with_bounds");
        initialized = true;
    }

    real min_x = 0.0;
    real min_y = 0.0;
    real min_z = 0.0;
    real max_x = 0.0;
    real max_y = 0.0;
    real max_z = 0.0;

    real p0[3] = {0, 0 ,0};
    real  n[3] = {0, 0 ,0};

    if(clip_with_plain) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (real, n[0], config->config_data.config, "normal_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (real, n[1], config->config_data.config, "normal_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (real, n[2], config->config_data.config, "normal_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (real, p0[0], config->config_data.config, "origin_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (real, p0[1], config->config_data.config, "origin_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (real, p0[2], config->config_data.config, "origin_z");
    }

    if(clip_with_bounds) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (real, min_x, config->config_data.config, "min_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (real, min_y, config->config_data.config, "min_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (real, min_z, config->config_data.config, "min_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (real, max_x, config->config_data.config, "max_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (real, max_y, config->config_data.config, "max_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (real, max_z, config->config_data.config, "max_z");

    }


    real l = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    real A = n[0] / l;
    real B = n[1] / l;
    real C = n[2] / l;
    real D = -(n[0]*p0[0] + n[1]*p0[1] + n[2]*p0[2]);

    double side;


    sds tmp = sdsnew (output_dir);
    sds c = sdsfromlonglong (count);
    tmp = sdscat (tmp, "/V_t_");
    tmp = sdscat (tmp, c);

    FILE *output_file = fopen (tmp, "w");

    struct cell_node *grid_cell = the_grid->first_cell;

    double center_x, center_y, center_z, half_face;
    double v;
    bool act = false;

    while (grid_cell != 0) {

        if (grid_cell->active) {

            center_x = grid_cell->center_x;
            center_y = grid_cell->center_y;
            center_z = grid_cell->center_z;

            if(clip_with_plain) {
                side = A*center_x + B*center_y + C* center_z + D;
                if (side < 0) {
                    grid_cell = grid_cell->next;
                    continue;
                }
            }

            if(clip_with_bounds) {
                bool ignore_cell = center_x < min_x || center_x > max_x ||
                                   center_y < min_y || center_y > max_y ||
                                   center_z < min_z || center_z > max_z;

                if(ignore_cell) {
                    grid_cell = grid_cell->next;
                    continue;
                }
            }

            v = grid_cell->v;
            half_face = grid_cell->half_face_length;

            if (grid_cell->v > vm_threshold) {
                act = true;
            }

            if(binary) {
                fwrite (&center_x, sizeof(center_x), 1, output_file);
                fwrite (&center_y, sizeof(center_y), 1, output_file);
                fwrite (&center_z, sizeof(center_z), 1, output_file);
                fwrite (&half_face, sizeof(half_face), 1, output_file);
                fwrite (&v, sizeof(v), 1, output_file);
            }
            else {
                fprintf(output_file, "%g,%g,%g,%g,%g\n", center_x, center_y, center_z, half_face, v);
            }
        }
        grid_cell = grid_cell->next;
    }

    count++;
    fclose (output_file);
    sdsfree (tmp);
    sdsfree (c);
    free(file_prefix);

    return act;

}


SAVE_MESH(save_as_vtk) {

    if(!initialized) {
        GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(file_prefix, config->config_data.config, "file_prefix");
        GET_PARAMETER_BINARY_VALUE (clip_with_plain, config->config_data.config, "clip_with_plain");
        GET_PARAMETER_BINARY_VALUE (clip_with_bounds, config->config_data.config, "clip_with_bounds");
        initialized = true;
    }

    real min_x = 0.0;
    real min_y = 0.0;
    real min_z = 0.0;
    real max_x = 0.0;
    real max_y = 0.0;
    real max_z = 0.0;

     real p0[3] = {0, 0 ,0};
     real  n[3] = {0, 0 ,0};

    if(clip_with_plain) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (real, n[0], config->config_data.config, "normal_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (real, n[1], config->config_data.config, "normal_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (real, n[2], config->config_data.config, "normal_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (real, p0[0], config->config_data.config, "origin_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (real, p0[1], config->config_data.config, "origin_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (real, p0[2], config->config_data.config, "origin_z");
    }

    if(clip_with_bounds) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (real, min_x, config->config_data.config, "min_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (real, min_y, config->config_data.config, "min_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (real, min_z, config->config_data.config, "min_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (real, max_x, config->config_data.config, "max_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (real, max_y, config->config_data.config, "max_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR (real, max_z, config->config_data.config, "max_z");
    }

    sds tmp = sdsnew (output_dir);
    sds c = sdsfromlonglong (count);
    tmp = sdscat (tmp, "/V_t_");
    tmp = sdscat (tmp, c);

    tmp = sdscat (tmp, ".vtk");


    FILE *output_file = fopen (tmp, "w");

    struct cell_node *grid_cell = the_grid->first_cell;

    double center_x, center_y, center_z, half_face;
    double v;
    bool act = false;
    float *values = NULL;
    int *cells = NULL;

    struct point_hash *hash = point_hash_create();

    struct point_3d aux1;
    struct point_3d aux2;
    struct point_3d aux3;
    struct point_3d aux4;
    struct point_3d aux5;
    struct point_3d aux6;
    struct point_3d aux7;
    struct point_3d aux8;

    int id = 0;
    int num_cells = 0;

    fprintf(output_file, "# vtk DataFile Version 4.2\n");
    fprintf(output_file, "vtk output\n");
    fprintf(output_file, "ASCII\n");
    fprintf(output_file, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(output_file, "                                                                                    \n");

    real l = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    real A = n[0] / l;
    real B = n[1] / l;
    real C = n[2] / l;
    real D = -(n[0]*p0[0] + n[1]*p0[1] + n[2]*p0[2]);
    
    double side;

    while (grid_cell != 0) {

        if (grid_cell->active) {

            center_x = grid_cell->center_x;
            center_y = grid_cell->center_y;
            center_z = grid_cell->center_z;

            if(clip_with_plain) {
               side = A*center_x + B*center_y + C* center_z + D;
               if (side < 0) {
                   grid_cell = grid_cell->next;
                   continue;
               }
            }

            if(clip_with_bounds) {
                bool ignore_cell = center_x < min_x || center_x > max_x ||
                                  center_y < min_y || center_y > max_y ||
                                  center_z < min_z || center_z > max_z;

                if(ignore_cell) {
                    grid_cell = grid_cell->next;
                    continue;
                }
            }


            v = grid_cell->v;
            half_face = grid_cell->half_face_length;

            if (count > 0) { //TODO: maybe this should be a parameter
                if (grid_cell->v > vm_threshold) {
                    act = true;
                }
            } else {
                act = true;
            }

            sb_push(values, v);

            aux1.x = center_x - half_face;
            aux1.y = center_y - half_face;
            aux1.z = center_z - half_face;

            aux2.x = center_x + half_face;
            aux2.y = center_y - half_face;
            aux2.z = center_z - half_face;

            aux3.x = center_x + half_face;
            aux3.y = center_y + half_face;
            aux3.z = center_z - half_face;

            aux4.x = center_x - half_face;
            aux4.y = center_y + half_face;
            aux4.z = center_z - half_face;

            aux5.x = center_x - half_face;
            aux5.y = center_y - half_face;
            aux5.z = center_z + half_face;

            aux6.x = center_x + half_face;
            aux6.y = center_y - half_face;
            aux6.z = center_z + half_face;

            aux7.x = center_x + half_face;
            aux7.y = center_y + half_face;
            aux7.z = center_z + half_face;

            aux8.x = center_x - half_face;
            aux8.y = center_y + half_face;
            aux8.z = center_z + half_face;


            if (point_hash_search(hash, aux1) == -1) {
                point_hash_insert(hash, aux1, id);
                id++;

                fprintf(output_file, "%lf %lf %lf\n", aux1.x, aux1.y, aux1.z);
            }

            if (point_hash_search(hash, aux2) == -1) {
                point_hash_insert(hash, aux2, id);
                id++;

                fprintf(output_file, "%lf %lf %lf\n", aux2.x, aux2.y, aux2.z);
            }

            if (point_hash_search(hash, aux3) == -1) {
                point_hash_insert(hash, aux3, id);
                id++;

                fprintf(output_file, "%lf %lf %lf\n", aux3.x, aux3.y, aux3.z);
            }

            if (point_hash_search(hash, aux4) == -1) {
                point_hash_insert(hash, aux4, id);
                id++;

                fprintf(output_file, "%lf %lf %lf\n", aux4.x, aux4.y, aux4.z);
            }

            if (point_hash_search(hash, aux5) == -1) {
                point_hash_insert(hash, aux5, id);
                id++;

                fprintf(output_file, "%lf %lf %lf\n", aux5.x, aux5.y, aux5.z);
            }

            if (point_hash_search(hash, aux6) == -1) {
                point_hash_insert(hash, aux6, id);
                id++;

                fprintf(output_file, "%lf %lf %lf\n", aux6.x, aux6.y, aux6.z);
            }

            if (point_hash_search(hash, aux7) == -1) {
                point_hash_insert(hash, aux7, id);
                id++;

                fprintf(output_file, "%lf %lf %lf\n", aux7.x, aux7.y, aux7.z);
            }

            if (point_hash_search(hash, aux8) == -1) {
                point_hash_insert(hash, aux8, id);
                id++;

                fprintf(output_file, "%lf %lf %lf\n", aux8.x, aux8.y, aux8.z);
            }

            sb_push(cells, point_hash_search(hash, aux1));
            sb_push(cells, point_hash_search(hash, aux2));
            sb_push(cells, point_hash_search(hash, aux3));
            sb_push(cells, point_hash_search(hash, aux4));
            sb_push(cells, point_hash_search(hash, aux5));
            sb_push(cells, point_hash_search(hash, aux6));
            sb_push(cells, point_hash_search(hash, aux7));
            sb_push(cells, point_hash_search(hash, aux8));
            num_cells++;


        }

        grid_cell = grid_cell->next;
    }

    fprintf(output_file, "\nCELLS %d %d\n", num_cells, 9*num_cells);

    int points_per_cell = 8;
    int cell_type = 12;

    for(int i = 0; i < num_cells; i++) {
        fprintf(output_file, "%d ", points_per_cell);
        for(int j = 0; j < points_per_cell; j++) {
            fprintf(output_file, "%d ", cells[points_per_cell*i + j]);
        }
        fprintf(output_file, "\n");
    }


    fprintf(output_file, "\nCELL_TYPES %d\n", num_cells);
    for(int i = 0; i < num_cells; i++) {
        fprintf(output_file, "%d\n", cell_type);
    }

    fprintf(output_file, "\nCELL_DATA %d\n", num_cells);
    fprintf(output_file, "SCALARS Scalars_ float\n");
    fprintf(output_file, "LOOKUP_TABLE default\n");
    {
        int num_values = sb_count(values);
        for(int i = 0; i < num_values; i++) {
            fprintf(output_file, "%lf ", values[i]);
        }
    }

    fprintf(output_file, "\nMETADATA\n");
    fprintf(output_file, "INFORMATION 0\n");

    fseek(output_file, 70, SEEK_SET);

    fprintf(output_file, "POINTS %d float", id);

    point_hash_destroy(hash);
    sb_free(cells);
    sb_free(values);

    count++;
    fclose (output_file);
    sdsfree (tmp);
    sdsfree (c);

    return act;
}