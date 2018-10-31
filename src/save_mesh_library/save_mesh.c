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
static FILE *pvd_file = NULL;

static bool initialized = false;
static bool first_save_call = true;
static int count = 0;

void add_file_to_pvd(double current_dt, const char *output_dir, const char *base_name);

//int invert_bytes(int data) {
//    int swapped = ((data >> 24) & 0xff) |      // move byte 3 to byte 0
//                  ((data << 8) & 0xff0000) |   // move byte 1 to byte 2
//                  ((data >> 8) & 0xff00) |     // move byte 2 to byte 1
//                  ((data << 24) & 0xff000000); // byte 0 to byte 3
//    return swapped;
//}
//
//void save_binary_float(FILE *output_file, struct point_3d *p) {
//    int a = *(int *)&(p->x);
//    int swapped = invert_bytes(a);
//
//    // fwrite(&aux1, sizeof(struct point_3d), 1, output_file);
//    fwrite(&swapped, sizeof(int), 1, output_file);
//
//    a = *(int *)&(p->y);
//    swapped = invert_bytes(a);
//
//    fwrite(&swapped, sizeof(int), 1, output_file);
//
//    a = *(int *)&(p->z);
//    swapped = invert_bytes(a);
//
//    fwrite(&swapped, sizeof(int), 1, output_file);
//}

SAVE_MESH(save_as_text_or_binary) {

    char *output_dir = config->out_dir_name;

    if(!initialized) {

        GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(file_prefix, config->config_data.config, "file_prefix");
        GET_PARAMETER_BINARY_VALUE(binary, config->config_data.config, "binary");
        GET_PARAMETER_BINARY_VALUE(clip_with_plain, config->config_data.config, "clip_with_plain");
        GET_PARAMETER_BINARY_VALUE(clip_with_bounds, config->config_data.config, "clip_with_bounds");
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
    tmp = sdscatprintf(tmp, "/V_t_%lf.vtk", current_dt);

    FILE *output_file = fopen(tmp, "w");

    struct cell_node *grid_cell = the_grid->first_cell;

    double center_x, center_y, center_z, half_face;
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
            half_face = grid_cell->half_face_length;

            if(binary) {
                fwrite(&center_x, sizeof(center_x), 1, output_file);
                fwrite(&center_y, sizeof(center_y), 1, output_file);
                fwrite(&center_z, sizeof(center_z), 1, output_file);
                fwrite(&half_face, sizeof(half_face), 1, output_file);
                fwrite(&v, sizeof(v), 1, output_file);
            } else {
                fprintf(output_file, "%g,%g,%g,%g,%g\n", center_x, center_y, center_z, half_face, v);
            }
        }
        grid_cell = grid_cell->next;
    }

    fclose(output_file);
    sdsfree(tmp);
    free(file_prefix);
}

SAVE_MESH(save_as_vtk) {

    char *output_dir = config->out_dir_name;

    if(!initialized) {
        GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(file_prefix, config->config_data.config, "file_prefix");
        GET_PARAMETER_BINARY_VALUE(clip_with_plain, config->config_data.config, "clip_with_plain");
        GET_PARAMETER_BINARY_VALUE(clip_with_bounds, config->config_data.config, "clip_with_bounds");
        GET_PARAMETER_BINARY_VALUE(binary, config->config_data.config, "binary");
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
    sds base_name = sdscatprintf(sdsempty(), "V_t_%d.vtk", count);
    count++;
    output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, current_dt);

    struct vtk_unstructured_grid *vtk_grid =
            new_vtk_unstructured_grid_from_alg_grid(the_grid, clip_with_plain, plain_coords, clip_with_bounds, bounds);

    save_vtk_unstructured_grid_as_legacy_vtk(vtk_grid, output_dir_with_file, binary);

    free_vtk_unstructured_grid(vtk_grid);
    sdsfree(output_dir_with_file);

//    char *output_dir = config->out_dir_name;
//
//    if(!initialized) {
//        GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(file_prefix, config->config_data.config, "file_prefix");
//        GET_PARAMETER_BINARY_VALUE(clip_with_plain, config->config_data.config, "clip_with_plain");
//        GET_PARAMETER_BINARY_VALUE(clip_with_bounds, config->config_data.config, "clip_with_bounds");
//        GET_PARAMETER_BINARY_VALUE(binary, config->config_data.config, "binary");
//        initialized = true;
//    }
//
//    real min_x = 0.0;
//    real min_y = 0.0;
//    real min_z = 0.0;
//    real max_x = 0.0;
//    real max_y = 0.0;
//    real max_z = 0.0;
//
//    real p0[3] = {0, 0, 0};
//    real n[3] = {0, 0, 0};
//
//    if(clip_with_plain) {
//        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, n[0], config->config_data.config, "normal_x");
//        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, n[1], config->config_data.config, "normal_y");
//        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, n[2], config->config_data.config, "normal_z");
//        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, p0[0], config->config_data.config, "origin_x");
//        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, p0[1], config->config_data.config, "origin_y");
//        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, p0[2], config->config_data.config, "origin_z");
//    }
//
//    if(clip_with_bounds) {
//        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, min_x, config->config_data.config, "min_x");
//        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, min_y, config->config_data.config, "min_y");
//        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, min_z, config->config_data.config, "min_z");
//        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, max_x, config->config_data.config, "max_x");
//        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, max_y, config->config_data.config, "max_y");
//        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, max_z, config->config_data.config, "max_z");
//    }
//
//    sds output_dir_with_file = sdsnew(output_dir);
//    sds base_name = sdscatprintf(sdsempty(), "V_%d.vtk", count);
//    count++;
//    output_dir_with_file = sdscat(output_dir_with_file, "/");
//    output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, current_dt);
//
//    FILE *output_file = fopen(output_dir_with_file, "w");
//
//    struct cell_node *grid_cell = the_grid->first_cell;
//
//    float center_x, center_y, center_z, half_face;
//    double v;
//
//    float *values = NULL;
//    int *cells = NULL;
//
//    struct point_hash *hash = point_hash_create();
//
//    struct point_3d aux1;
//    struct point_3d aux2;
//    struct point_3d aux3;
//    struct point_3d aux4;
//    struct point_3d aux5;
//    struct point_3d aux6;
//    struct point_3d aux7;
//    struct point_3d aux8;
//
//    int id = 0;
//    int num_cells = 0;
//
//    fprintf(output_file, "# vtk DataFile Version 4.2\n");
//    fprintf(output_file, "vtk output\n");
//    if(binary) {
//        fprintf(output_file, "BINARY\n");
//    } else {
//        fprintf(output_file, "ASCII \n");
//    }
//
//    fprintf(output_file, "DATASET UNSTRUCTURED_GRID\n");
//    fprintf(output_file, "                                                                                    \n");
//
//    real l = sqrtf(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
//    real A = n[0] / l;
//    real B = n[1] / l;
//    real C = n[2] / l;
//    real D = -(n[0] * p0[0] + n[1] * p0[1] + n[2] * p0[2]);
//
//    double side;
//
//    while(grid_cell != 0) {
//
//        if(grid_cell->active) {
//
//            center_x = grid_cell->center_x;
//            center_y = grid_cell->center_y;
//            center_z = grid_cell->center_z;
//
//            if(clip_with_plain) {
//                side = A * center_x + B * center_y + C * center_z + D;
//                if(side < 0) {
//                    grid_cell = grid_cell->next;
//                    continue;
//                }
//            }
//
//            if(clip_with_bounds) {
//                bool ignore_cell = center_x < min_x || center_x > max_x || center_y < min_y || center_y > max_y ||
//                                   center_z < min_z || center_z > max_z;
//
//                if(ignore_cell) {
//                    grid_cell = grid_cell->next;
//                    continue;
//                }
//            }
//
//            v = grid_cell->v;
//            half_face = grid_cell->half_face_length;
//
//            sb_push(values, v);
//
//            aux1.x = center_x - half_face;
//            aux1.y = center_y - half_face;
//            aux1.z = center_z - half_face;
//
//            aux2.x = center_x + half_face;
//            aux2.y = center_y - half_face;
//            aux2.z = center_z - half_face;
//
//            aux3.x = center_x + half_face;
//            aux3.y = center_y + half_face;
//            aux3.z = center_z - half_face;
//
//            aux4.x = center_x - half_face;
//            aux4.y = center_y + half_face;
//            aux4.z = center_z - half_face;
//
//            aux5.x = center_x - half_face;
//            aux5.y = center_y - half_face;
//            aux5.z = center_z + half_face;
//
//            aux6.x = center_x + half_face;
//            aux6.y = center_y - half_face;
//            aux6.z = center_z + half_face;
//
//            aux7.x = center_x + half_face;
//            aux7.y = center_y + half_face;
//            aux7.z = center_z + half_face;
//
//            aux8.x = center_x - half_face;
//            aux8.y = center_y + half_face;
//            aux8.z = center_z + half_face;
//
//            if(point_hash_search(hash, aux1) == -1) {
//                point_hash_insert(hash, aux1, id);
//                id++;
//
//                if(binary) {
//                    save_binary_float(output_file, &aux1);
//                } else {
//                    fprintf(output_file, "%lf %lf %lf\n", aux1.x, aux1.y, aux1.z);
//                }
//            }
//
//            if(point_hash_search(hash, aux2) == -1) {
//                point_hash_insert(hash, aux2, id);
//                id++;
//
//                if(binary) {
//                    save_binary_float(output_file, &aux2);
//
//                } else {
//                    fprintf(output_file, "%lf %lf %lf\n", aux2.x, aux2.y, aux2.z);
//                }
//            }
//
//            if(point_hash_search(hash, aux3) == -1) {
//                point_hash_insert(hash, aux3, id);
//                id++;
//
//                if(binary) {
//                    save_binary_float(output_file, &aux3);
//                } else {
//                    fprintf(output_file, "%lf %lf %lf\n", aux3.x, aux3.y, aux3.z);
//                }
//            }
//
//            if(point_hash_search(hash, aux4) == -1) {
//                point_hash_insert(hash, aux4, id);
//                id++;
//
//                if(binary) {
//                    save_binary_float(output_file, &aux4);
//                } else {
//                    fprintf(output_file, "%lf %lf %lf\n", aux4.x, aux4.y, aux4.z);
//                }
//            }
//
//            if(point_hash_search(hash, aux5) == -1) {
//                point_hash_insert(hash, aux5, id);
//                id++;
//
//                if(binary) {
//                    save_binary_float(output_file, &aux5);
//
//                } else {
//                    fprintf(output_file, "%lf %lf %lf\n", aux5.x, aux5.y, aux5.z);
//                }
//            }
//
//            if(point_hash_search(hash, aux6) == -1) {
//                point_hash_insert(hash, aux6, id);
//                id++;
//
//                if(binary) {
//                    save_binary_float(output_file, &aux6);
//                } else {
//                    fprintf(output_file, "%lf %lf %lf\n", aux6.x, aux6.y, aux6.z);
//                }
//            }
//
//            if(point_hash_search(hash, aux7) == -1) {
//                point_hash_insert(hash, aux7, id);
//                id++;
//
//                if(binary) {
//                    save_binary_float(output_file, &aux7);
//                } else {
//                    fprintf(output_file, "%lf %lf %lf\n", aux7.x, aux7.y, aux7.z);
//                }
//            }
//
//            if(point_hash_search(hash, aux8) == -1) {
//                point_hash_insert(hash, aux8, id);
//                id++;
//
//                if(binary) {
//                    save_binary_float(output_file, &aux8);
//
//                } else {
//                    fprintf(output_file, "%lf %lf %lf\n", aux8.x, aux8.y, aux8.z);
//                }
//            }
//
//            sb_push(cells, point_hash_search(hash, aux1));
//            sb_push(cells, point_hash_search(hash, aux2));
//            sb_push(cells, point_hash_search(hash, aux3));
//            sb_push(cells, point_hash_search(hash, aux4));
//            sb_push(cells, point_hash_search(hash, aux5));
//            sb_push(cells, point_hash_search(hash, aux6));
//            sb_push(cells, point_hash_search(hash, aux7));
//            sb_push(cells, point_hash_search(hash, aux8));
//            num_cells++;
//        }
//
//        grid_cell = grid_cell->next;
//    }
//
//    fprintf(output_file, "\nCELLS %d %d\n", num_cells, 9 * num_cells);
//
//    int points_per_cell = 8;
//    int cell_type = 12;
//
//    int points_per_cell_swapped = invert_bytes(points_per_cell);
//    int cell_type_swapped = invert_bytes(cell_type);
//
//    for(int i = 0; i < num_cells; i++) {
//        if(binary) {
//            fwrite(&points_per_cell_swapped, sizeof(int), 1, output_file);
//        } else {
//            fprintf(output_file, "%d ", points_per_cell);
//        }
//
//        for(int j = 0; j < points_per_cell; j++) {
//            if(binary) {
//                int aux = invert_bytes(cells[points_per_cell * i + j]);
//                fwrite(&aux, sizeof(int), 1, output_file);
//            } else {
//                fprintf(output_file, "%d ", cells[points_per_cell * i + j]);
//            }
//        }
//
//        if(!binary)
//            fprintf(output_file, "\n");
//    }
//
//    fprintf(output_file, "\nCELL_TYPES %d\n", num_cells);
//    for(int i = 0; i < num_cells; i++) {
//        if(binary) {
//            fwrite(&cell_type_swapped, sizeof(int), 1, output_file);
//        } else {
//            fprintf(output_file, "%d\n", cell_type);
//        }
//    }
//
//    fprintf(output_file, "\nCELL_DATA %d\n", num_cells);
//    fprintf(output_file, "SCALARS Scalars_ float\n");
//    fprintf(output_file, "LOOKUP_TABLE default\n");
//    {
//        size_t num_values = sb_count(values);
//
//        for(int i = 0; i < num_values; i++) {
//            if(binary) {
//                int aux = invert_bytes(*((int *)&values[i]));
//                fwrite(&aux, sizeof(int), 1, output_file);
//            } else {
//                fprintf(output_file, "%lf ", values[i]);
//            }
//        }
//    }
//
//    fprintf(output_file, "\nMETADATA\n");
//    fprintf(output_file, "INFORMATION 0\n\n");
//
//    fseek(output_file, 71, SEEK_SET);
//
//    fprintf(output_file, "POINTS %d float", id);
//
//    point_hash_destroy(hash);
//    sb_free(cells);
//    sb_free(values);
//
//    fclose(output_file);
//    sdsfree(output_dir_with_file);
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
        GET_PARAMETER_BINARY_VALUE(clip_with_plain, config->config_data.config, "clip_with_plain");
        GET_PARAMETER_BINARY_VALUE(clip_with_bounds, config->config_data.config, "clip_with_bounds");
        GET_PARAMETER_BINARY_VALUE(binary, config->config_data.config, "binary");
        GET_PARAMETER_BINARY_VALUE(save_pvd, config->config_data.config, "save_pvd");
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
    sds base_name = sdscatprintf(sdsempty(), "V_t_%d.vtu", count);
    count++;
    output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, current_dt);

    if(save_pvd) {
        add_file_to_pvd(current_dt, output_dir, base_name);
    }

    struct vtk_unstructured_grid *vtk_grid =
        new_vtk_unstructured_grid_from_alg_grid(the_grid, clip_with_plain, plain_coords, clip_with_bounds, bounds);

    save_vtk_unstructured_grid_as_vtu(vtk_grid, output_dir_with_file, binary);

    free_vtk_unstructured_grid(vtk_grid);
    sdsfree(output_dir_with_file);

}