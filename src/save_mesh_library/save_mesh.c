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

SAVE_MESH(save_as_text_or_binary) {

    char * file_prefix;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(file_prefix, config->config_data.config, "file_prefix");

    bool binary;
    GET_PARAMETER_BINARY_VALUE (binary, config->config_data.config, "binary");

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

    return act;

}


SAVE_MESH(save_as_vtk) {

    char * file_prefix;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(file_prefix, config->config_data.config, "file_prefix");

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

    while (grid_cell != 0) {

        if (grid_cell->active) {

            center_x = grid_cell->center_x;
            center_y = grid_cell->center_y;
            center_z = grid_cell->center_z;

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