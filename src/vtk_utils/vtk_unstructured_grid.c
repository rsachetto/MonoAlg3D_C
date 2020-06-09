//
// Created by sachetto on 30/10/18.
//

#include "vtk_unstructured_grid.h"
#include "../alg/cell/cell.h"
#include "../3dparty/sds/sds.h"
#include "data_utils.h"
#include "../common_types/common_types.h"
#include "../3dparty/stb_ds.h"
#include "../utils/file_utils.h"
#include "../3dparty/xml_parser/yxml.h"

#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <ctype.h>
#include <sys/mman.h>
#include <float.h>

struct vtk_unstructured_grid *new_vtk_unstructured_grid() {
    struct vtk_unstructured_grid *grid = (struct vtk_unstructured_grid *)malloc(sizeof(struct vtk_unstructured_grid));

    grid->cells = NULL;
    grid->values = NULL;
    grid->points = NULL;

    grid->num_cells = 0;
    grid->num_points = 0;
    grid->points_per_cell = 8;
    grid->cell_type = 12;

    grid->max_v = FLT_MIN;
    grid->min_v = FLT_MAX;

    return grid;
}

void free_vtk_unstructured_grid(struct vtk_unstructured_grid *vtk_grid) {
    if(vtk_grid) {
        arrfree(vtk_grid->cells);
        arrfree(vtk_grid->values);
        arrfree(vtk_grid->points);
        free(vtk_grid);
    }
}

static inline void set_point_data(point3d_array points, struct point_3d center, struct point_3d half_face) {
    real_cpu center_x_plus  = center.x + half_face.x;
    real_cpu center_x_minus = center.x - half_face.x;

    real_cpu center_y_plus  = center.y + half_face.y;
    real_cpu center_y_minus = center.y - half_face.y;

    real_cpu center_z_plus  = center.z + half_face.z;
    real_cpu center_z_minus = center.z - half_face.z;

    points[0].x = center_x_minus;
    points[0].y = center_y_minus;
    points[0].z = center_z_minus;

    points[1].x = center_x_plus;
    points[1].y = center_y_minus;
    points[1].z = center_z_minus;

    points[2].x = center_x_plus;
    points[2].y = center_y_plus;
    points[2].z = center_z_minus;

    points[3].x = center_x_minus;
    points[3].y = center_y_plus;
    points[3].z = center_z_minus;

    points[4].x = center_x_minus;
    points[4].y = center_y_minus;
    points[4].z = center_z_plus;

    points[5].x = center_x_plus;
    points[5].y = center_y_minus;
    points[5].z = center_z_plus;

    points[6].x = center_x_plus;
    points[6].y = center_y_plus;
    points[6].z = center_z_plus;

    points[7].x = center_x_minus;
    points[7].y = center_y_plus;
    points[7].z = center_z_plus;
}

void new_vtk_unstructured_grid_from_string_with_activation_info(struct vtk_unstructured_grid **vtk_grid, char* source, size_t source_size) {

    *vtk_grid = new_vtk_unstructured_grid();

    struct point_3d center;
    struct point_3d half_face;
    float  v;

    struct point_3d points[8];

    struct point_3d point1;
    struct point_3d point2;
    struct point_3d point3;
    struct point_3d point4;
    struct point_3d point5;
    struct point_3d point6;
    struct point_3d point7;
    struct point_3d point8;

    uint32_t id = 0;
    uint32_t num_cells = 0;

    struct point_hash_entry *hash =  NULL;
    char *line = NULL;

    int active, fibrotic, bz;

    while(*source) {

        int line_end = 0;
        while (*source != '[') {
            arrput(line, *source);
            source++;
            source_size--;
            line_end++;
        }
        source++;
        source_size--;

        while (*source != '\n') {
            source++;
            source_size--;
        }
        source++;
        source_size--;

        line[line_end-1] = '\0';

        sscanf(line, "%lf,%lf,%lf,%lf,%lf,%lf,%d,%d,%d %f", &center.x, &center.y, &center.z, &half_face.x, &half_face.y, &half_face.z, &active, &fibrotic, &bz, &v);

        arrsetlen(line, 0);

        if(!active) continue;

        arrput((*vtk_grid)->values, v);
        if(v > (*vtk_grid)->max_v) (*vtk_grid)->max_v = v;
        if(v < (*vtk_grid)->min_v) (*vtk_grid)->min_v = v;

        set_point_data(points, center, half_face);

        point1 = points[0];
        point2 = points[1];
        point3 = points[2];
        point4 = points[3];
        point5 = points[4];
        point6 = points[5];
        point7 = points[6];
        point8 = points[7];

        if(hmgeti(hash, point1) == -1) {
            arrput((*vtk_grid)->points, point1);
            hmput(hash, point1, id);
            id++;
        }

        if(hmgeti(hash, point2) == -1) {
            arrput((*vtk_grid)->points, point2);
            hmput(hash, point2, id);
            id++;
        }

        if(hmgeti(hash, point3) == -1) {
            hmput(hash, point3, id);
            arrput((*vtk_grid)->points, point3);
            id++;
        }

        if(hmgeti(hash, point4) == -1) {
            hmput(hash, point4, id);
            arrput((*vtk_grid)->points, point4);
            id++;
        }

        if(hmgeti(hash, point5) == -1) {
            arrput((*vtk_grid)->points, point5);
            hmput(hash, point5, id);
            id++;
        }

        if(hmgeti(hash, point6) == -1) {
            arrput((*vtk_grid)->points, point6);
            hmput(hash, point6, id);
            id++;
        }

        if(hmgeti(hash, point7) == -1) {
            arrput((*vtk_grid)->points, point7);
            hmput(hash, point7, id);
            id++;
        }

        if(hmgeti(hash, point8) == -1) {
            arrput((*vtk_grid)->points, point8);
            hmput(hash, point8, id);
            id++;
        }

        arrput((*vtk_grid)->cells, hmget(hash, point1));
        arrput((*vtk_grid)->cells, hmget(hash, point2));
        arrput((*vtk_grid)->cells, hmget(hash, point3));
        arrput((*vtk_grid)->cells, hmget(hash, point4));
        arrput((*vtk_grid)->cells, hmget(hash, point5));
        arrput((*vtk_grid)->cells, hmget(hash, point6));
        arrput((*vtk_grid)->cells, hmget(hash, point7));
        arrput((*vtk_grid)->cells, hmget(hash, point8));
        num_cells++;
    }

    (*vtk_grid)->num_cells = num_cells;
    (*vtk_grid)->num_points = id;

    arrfree(line);
    hmfree(hash);
}

void binary_grid_error(struct vtk_unstructured_grid **vtk_grid) {
    free(*vtk_grid);
    *vtk_grid = NULL;
}

void new_vtk_unstructured_grid_from_string(struct vtk_unstructured_grid **vtk_grid, char* source, size_t source_size, bool binary, bool read_only_values) {

    if(!read_only_values) {
        *vtk_grid = new_vtk_unstructured_grid();
    }
    else {
        if(!(*vtk_grid)) {
            fprintf(stderr,
                    "new_vtk_unstructured_grid_from_string can only be called with read_only_values if the grid is already loaded\n");
            exit(EXIT_FAILURE);
        }

        arrfree((*vtk_grid)->values);
        (*vtk_grid)->values = NULL;
    }

    struct point_3d center;
    struct point_3d half_face;

    float  v;

    struct point_3d points[8];
    
    struct point_3d point1;
    struct point_3d point2;
    struct point_3d point3;
    struct point_3d point4;
    struct point_3d point5;
    struct point_3d point6;
    struct point_3d point7;
    struct point_3d point8;

    uint32_t id = 0;
    uint32_t num_cells = 0;

    struct point_hash_entry *hash =  NULL;
    char *line = NULL;

    char* source_limit = source + source_size;

    while(source_size) {

        if(!binary) {

            while (*source != '\n') {
                arrput(line, *source);
                source++;
                source_size--;
            }
            source++;
            source_size--;

            arrput(line, '\0');

            sscanf(line, "%lf,%lf,%lf,%lf,%lf,%lf,%f", &center.x, &center.y, &center.z, &half_face.x, &half_face.y, &half_face.z, &v);

            arrsetlen(line, 0);
        }
        else {

            if(source < source_limit) {
                center.x = *(float *)(source);
                source += sizeof(center.x);
                source_size -= sizeof(center.x);
            }
            else {
                binary_grid_error(vtk_grid);
                return;
            }


            if(source < source_limit) {
                center.y = *(float *)(source);
                source += sizeof(center.y);
                source_size -= sizeof(center.y);
            }
            else {
                binary_grid_error(vtk_grid);
                return;
            }

            if(source < source_limit) {
                center.z = *(float *)(source);
                source += sizeof(center.z);
                source_size -= sizeof(center.z);
            }
            else {
                binary_grid_error(vtk_grid);
                return;
            }

            if(source < source_limit) {
                half_face.x = *(float *)(source);
                source += sizeof(half_face.x);
                source_size -= sizeof(half_face.x);
            }
            else {
                binary_grid_error(vtk_grid);
                return;
            }

            if(source < source_limit) {
                half_face.y = *(float *)(source);
                source += sizeof(half_face.y);
                source_size -= sizeof(half_face.y);
            }
            else {
                binary_grid_error(vtk_grid);
                return;
            }

            if(source < source_limit) {
                half_face.z = *(float *)(source);
                source += sizeof(half_face.z);
                source_size -= sizeof(half_face.z);
            }
            else {
                binary_grid_error(vtk_grid);
                return;
            }

            if(source < source_limit) {
                v = *(float *)(source);
                source += sizeof(v);
                source_size -= sizeof(v);
            }
            else {
                binary_grid_error(vtk_grid);
                return;
            }

        }

        arrput((*vtk_grid)->values, v);

        if(v > (*vtk_grid)->max_v) (*vtk_grid)->max_v = v;
        if(v < (*vtk_grid)->min_v) (*vtk_grid)->min_v = v;

        if(read_only_values) {
            if(!binary) {
                while (*source != '\n') {
                    arrput(line, *source);
                    source++;
                }
            }
            else {
                source += sizeof(center.x);
                source_size -= sizeof(center.x);

                source += sizeof(center.y);
                source_size -= sizeof(center.y);

                source += sizeof(center.z);
                source_size -= sizeof(center.z);

                source += sizeof(half_face.x);
                source_size -= sizeof(half_face.x);

                source += sizeof(half_face.y);
                source_size -= sizeof(half_face.y);

                source += sizeof(half_face.z);
                source_size -= sizeof(half_face.z);

                source += sizeof(v);
                source_size -= sizeof(v);
            }
            continue;
        }

        set_point_data(points, center, half_face);

        point1 = points[0];
        point2 = points[1];
        point3 = points[2];
        point4 = points[3];
        point5 = points[4];
        point6 = points[5];
        point7 = points[6];
        point8 = points[7];


        if(hmgeti(hash, point1) == -1) {
            arrput((*vtk_grid)->points, point1);
            hmput(hash, point1, id);
            id++;
        }

        if(hmgeti(hash, point2) == -1) {
            arrput((*vtk_grid)->points, point2);
            hmput(hash, point2, id);
            id++;
        }

        if(hmgeti(hash, point3) == -1) {
            hmput(hash, point3, id);
            arrput((*vtk_grid)->points, point3);
            id++;
        }

        if(hmgeti(hash, point4) == -1) {
            hmput(hash, point4, id);
            arrput((*vtk_grid)->points, point4);
            id++;
        }

        if(hmgeti(hash, point5) == -1) {
            arrput((*vtk_grid)->points, point5);
            hmput(hash, point5, id);
            id++;
        }

        if(hmgeti(hash, point6) == -1) {
            arrput((*vtk_grid)->points, point6);
            hmput(hash, point6, id);
            id++;
        }

        if(hmgeti(hash, point7) == -1) {
            arrput((*vtk_grid)->points, point7);
            hmput(hash, point7, id);
            id++;
        }

        if(hmgeti(hash, point8) == -1) {
            arrput((*vtk_grid)->points, point8);
            hmput(hash, point8, id);
            id++;
        }

        arrput((*vtk_grid)->cells, hmget(hash, point1));
        arrput((*vtk_grid)->cells, hmget(hash, point2));
        arrput((*vtk_grid)->cells, hmget(hash, point3));
        arrput((*vtk_grid)->cells, hmget(hash, point4));
        arrput((*vtk_grid)->cells, hmget(hash, point5));
        arrput((*vtk_grid)->cells, hmget(hash, point6));
        arrput((*vtk_grid)->cells, hmget(hash, point7));
        arrput((*vtk_grid)->cells, hmget(hash, point8));
        num_cells++;
    }

    if(!read_only_values) {
        (*vtk_grid)->num_cells = num_cells;
        (*vtk_grid)->num_points = id;
    }
    arrfree(line);
    hmfree(hash);
}

void new_vtk_unstructured_grid_from_alg_grid(struct vtk_unstructured_grid **vtk_grid, struct grid *grid, bool clip_with_plain,
                                             float *plain_coordinates, bool clip_with_bounds,
                                             float *bounds, bool read_only_values) {
    
    if(grid == NULL) {
        return;
    }

    uint32_t num_active_cells =  grid->num_active_cells;

    if(!read_only_values) {
        *vtk_grid = new_vtk_unstructured_grid();
        arrsetcap((*vtk_grid)->values, num_active_cells);
        arrsetcap((*vtk_grid)->cells,  num_active_cells);
        arrsetcap((*vtk_grid)->points, num_active_cells);
    }
    else {
        if(!(*vtk_grid)) {
            fprintf(stderr,
                    "Function new_vtk_unstructured_grid_from_alg_grid can only be called with read_only_values if the grid is already loaded!\n");
            exit(EXIT_FAILURE);
        }

        assert(*vtk_grid);
        arrfree((*vtk_grid)->values);
        (*vtk_grid)->values = NULL;
    }

    real_cpu min_x = 0.0;
    real_cpu min_y = 0.0;
    real_cpu min_z = 0.0;
    real_cpu max_x = 0.0;
    real_cpu max_y = 0.0;
    real_cpu max_z = 0.0;

    real_cpu p0[3] = {0, 0, 0};
    real_cpu n[3] = {0, 0, 0};

    if(!plain_coordinates) {
        clip_with_plain = false;
    } else {
        p0[0] = plain_coordinates[0];
        p0[1] = plain_coordinates[1];
        p0[2] = plain_coordinates[2];

        n[0] = plain_coordinates[3];
        n[0] = plain_coordinates[4];
        n[0] = plain_coordinates[5];
    }

    if(!bounds) {
        clip_with_bounds = false;
    } else {
        min_x = bounds[0];
        min_y = bounds[1];
        min_z = bounds[2];
        max_x = bounds[3];
        max_y = bounds[4];
        max_z = bounds[5];
    }

    struct point_3d points[8];
    struct point_3d point1, point2, point3, point4, point5, point6, point7, point8;

    uint32_t id = 0;
    uint32_t num_cells = 0;

    float l = sqrtf(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    float A = n[0] / l;
    float B = n[1] / l;
    float C = n[2] / l;
    float D = -(n[0] * p0[0] + n[1] * p0[1] + n[2] * p0[2]);

    real_cpu side;
    struct point_hash_entry *hash =  NULL;

    struct point_3d half_face;
    struct point_3d center;

    struct cell_node ** grid_cell = grid->active_cells;

    for(int i = 0; i < num_active_cells; i++) {

        center = grid_cell[i]->center;

        if(clip_with_plain) {
            side = A * center.x + B * center.y + C * center.z + D;
            if(side < 0) {
                continue;
            }
        }

        if(clip_with_bounds) {
            bool ignore_cell = center.x < min_x || center.x > max_x || center.y < min_y || center.y > max_y ||
                center.z < min_z || center.z > max_z;

            if(ignore_cell) {
                continue;
            }
        }

        arrput((*vtk_grid)->values, grid_cell[i]->v);

        if(grid_cell[i]->v > (*vtk_grid)->max_v) (*vtk_grid)->max_v = grid_cell[i]->v;
        if(grid_cell[i]->v < (*vtk_grid)->min_v) (*vtk_grid)->min_v = grid_cell[i]->v;

        if(read_only_values) {
            continue;
        }

        half_face.x = grid_cell[i]->discretization.x / 2.0f;
        half_face.y = grid_cell[i]->discretization.y / 2.0f;
        half_face.z = grid_cell[i]->discretization.z / 2.0f;

        set_point_data(points, center, half_face);

        point1 = points[0];
        point2 = points[1];
        point3 = points[2];
        point4 = points[3];
        point5 = points[4];
        point6 = points[5];
        point7 = points[6];
        point8 = points[7];

        if(hmgeti(hash, point1) == -1) {
            arrput((*vtk_grid)->points, point1);
            hmput(hash, point1, id);
            id++;
        }

        if(hmgeti(hash, point2) == -1) {
            arrput((*vtk_grid)->points, point2);
            hmput(hash, point2, id);
            id++;
        }

        if(hmgeti(hash, point3) == -1) {
            hmput(hash, point3, id);
            arrput((*vtk_grid)->points, point3);
            id++;
        }

        if(hmgeti(hash, point4) == -1) {
            hmput(hash, point4, id);
            arrput((*vtk_grid)->points, point4);
            id++;
        }

        if(hmgeti(hash, point5) == -1) {
            arrput((*vtk_grid)->points, point5);
            hmput(hash, point5, id);
            id++;
        }

        if(hmgeti(hash, point6) == -1) {
            arrput((*vtk_grid)->points, point6);
            hmput(hash, point6, id);
            id++;
        }

        if(hmgeti(hash, point7) == -1) {
            arrput((*vtk_grid)->points, point7);
            hmput(hash, point7, id);
            id++;
        }

        if(hmgeti(hash, point8) == -1) {
            arrput((*vtk_grid)->points, point8);
            hmput(hash, point8, id);
            id++;
        }

        arrput((*vtk_grid)->cells, hmget(hash, point1));
        arrput((*vtk_grid)->cells, hmget(hash, point2));
        arrput((*vtk_grid)->cells, hmget(hash, point3));
        arrput((*vtk_grid)->cells, hmget(hash, point4));
        arrput((*vtk_grid)->cells, hmget(hash, point5));
        arrput((*vtk_grid)->cells, hmget(hash, point6));
        arrput((*vtk_grid)->cells, hmget(hash, point7));
        arrput((*vtk_grid)->cells, hmget(hash, point8));
        num_cells++;
    }

    if(!read_only_values) {
        (*vtk_grid)->num_cells = num_cells;
        (*vtk_grid)->num_points = id;
    }

    hmfree(hash);
}

sds create_common_vtu_header(bool compressed, int num_points, int num_cells) {

    sds header = sdsempty();

    if(compressed) {
        header = sdscat(header, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" "
                                "header_type=\"UInt64\" compressor=\"vtkZLibDataCompressor\">\n");
    } else {
        header = sdscat(
                header,
                "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
    }

    header = sdscat(header, "  <UnstructuredGrid>\n");

    header = sdscatprintf(header, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", num_points, num_cells);

    header = sdscat(header, "      <CellData Scalars=\"Scalars_\">\n");

    return header;
}

void save_vtk_unstructured_grid_as_vtu(struct vtk_unstructured_grid *vtk_grid, char *filename, bool binary) {

    size_t offset = 0;

    uint32_t num_cells = vtk_grid->num_cells;
    uint32_t num_points = vtk_grid->num_points;

    sds file_content = create_common_vtu_header(false, num_points, num_cells);

    if(binary) {
        // First offset is always 0
        file_content = sdscat(file_content,"        <DataArray type=\"Float32\" Name=\"Scalars_\" format=\"appended\" offset=\"0\">\n");


    } else {
        file_content = sdscat(file_content, "        <DataArray type=\"Float32\" Name=\"Scalars_\" format=\"ascii\">\n");
    }

    if(!binary) {
        size_t num_values = arrlenu(vtk_grid->values);

        for(int i = 0; i < num_values; i++) {
            file_content = sdscatprintf(file_content, "     %lf ", vtk_grid->values[i]);
        }
    }

    file_content = sdscat(file_content, "        </DataArray>\n");
    file_content = sdscat(file_content, "      </CellData>\n");

    offset = (num_cells * 4) + 8;

    file_content = sdscat(file_content, "      <Points>\n");

    char data_size[3] = "32";

    if(sizeof(real_cpu) == 8) {
        data_size[0] = '6';
        data_size[1] = '4';

    }

    if(binary) {
        file_content = sdscatprintf(file_content,
                                    "        <DataArray type=\"Float%s\" Name=\"Points\" "
                                    "NumberOfComponents=\"3\" format=\"appended\" offset=\"%zu\">\n",
                                    data_size, offset);

    } else {
        file_content =
                sdscatprintf(file_content,
                       "        <DataArray type=\"Float%s\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n", data_size);
    }

    if(!binary) {
        for(int i = 0; i < num_points; i++) {
            struct point_3d p = vtk_grid->points[i];
            file_content = sdscatprintf(file_content, "%lf %lf %lf\n", p.x, p.y, p.z);
        }
    }

    file_content = sdscat(file_content, "        </DataArray>\n");
    file_content = sdscat(file_content, "      </Points>\n");

    file_content = sdscat(file_content, "      <Cells>\n");

    offset += (num_points * sizeof(real_cpu) * 3) + 8; // 3 float point number for each point

    if(binary) {
        file_content = sdscatprintf(
                file_content,
                "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"appended\" offset=\"%zu\">\n", offset);

    } else {
        file_content =
                sdscat(file_content, "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    }

    int points_per_cell = vtk_grid->points_per_cell;
    int cell_type = vtk_grid->cell_type;

    if(!binary) {
        for(int i = 0; i < num_cells; i++) {
            file_content = sdscat(file_content, "     ");
            for(int j = 0; j < points_per_cell; j++) {
                file_content = sdscatprintf(file_content, "%ld ", vtk_grid->cells[points_per_cell * i + j]);
            }

            file_content = sdscat(file_content, "\n");
        }
    }

    file_content = sdscat(file_content, "        </DataArray>\n");

    offset += (num_cells * 8 * points_per_cell) + 8; // 64 bits

    if(binary) {
        file_content = sdscatprintf(
                file_content, "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"appended\" offset=\"%zu\">\n",
                offset);
    } else {
        file_content = sdscat(file_content, "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    }

    {
        int offset_local = 8;

        if(!binary) {
            for(int i = 0; i < num_cells; i++) {
                file_content = sdscat(file_content, "     ");
                file_content = sdscatprintf(file_content, "%d ", offset_local);
                offset_local += 8;
                file_content = sdscat(file_content, "\n");
            }
        }
    }

    file_content = sdscat(file_content, "        </DataArray>\n");

    offset += (num_cells * 8) + 8; // 64 bits

    if(binary) {
        file_content = sdscatprintf(
                file_content, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\"%zu\">\n",
                offset);
    } else {
        file_content = sdscat(file_content, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    }

    if(!binary) {

        for(int i = 0; i < num_cells; i++) {
            file_content = sdscatprintf(file_content, "    %d\n", cell_type);
        }
    }

    file_content = sdscat(file_content, "        </DataArray>\n");

    file_content = sdscat(file_content, "      </Cells>\n");

    file_content = sdscat(file_content, "    </Piece>\n");
    file_content = sdscat(file_content, "  </UnstructuredGrid>\n");

    size_t size_until_now = 0;

    if(binary) {
        file_content = sdscat(file_content, "  <AppendedData encoding=\"raw\">\n   _");

        size_until_now = sdslen(file_content);

        // scalars
        uint64_t block_size = sizeof(float) * num_cells;
        file_content = sdscatlen(file_content, &block_size, sizeof(uint64_t));
        file_content = sdscatlen(file_content, vtk_grid->values, (size_t)block_size);
        size_until_now += (block_size + sizeof(uint64_t));

        // Points
        block_size = sizeof(struct point_3d) * vtk_grid->num_points;
        file_content = sdscatlen(file_content, &block_size, sizeof(uint64_t));
        file_content = sdscatlen(file_content, vtk_grid->points, (size_t)block_size);
        size_until_now += (block_size + sizeof(uint64_t));

        // connectivity
        block_size = num_cells * points_per_cell * sizeof(int64_t);
        file_content = sdscatlen(file_content, &block_size, sizeof(uint64_t));
        size_until_now += sizeof(uint64_t);

        for(int i = 0; i < num_cells; i++) {
            for(int j = 0; j < points_per_cell; j++) {
                int64_t aux = (int64_t)vtk_grid->cells[points_per_cell * i + j];
                file_content = sdscatlen(file_content, &aux, sizeof(int64_t));
                size_until_now += sizeof(int64_t);
            }
        }

        // offsets
        block_size = num_cells * sizeof(int64_t);
        file_content = sdscatlen(file_content, &block_size, sizeof(uint64_t));
        size_until_now += sizeof(uint64_t);

        int64_t offset_local = 8;
        for(int i = 0; i < num_cells; i++) {
            file_content = sdscatlen(file_content, &offset_local, sizeof(int64_t));
            offset_local += 8;
            size_until_now += sizeof(int64_t);
        }

        // types
        block_size = num_cells * sizeof(uint8_t);
        file_content = sdscatlen(file_content, &block_size, sizeof(uint64_t));
        size_until_now += sizeof(uint64_t);

        for(int i = 0; i < num_cells; i++) {
            uint8_t aux = (uint8_t)cell_type;
            file_content = sdscatlen(file_content, &aux, sizeof(uint8_t));
            size_until_now += sizeof(uint8_t);
        }

        file_content = sdscat(file_content, "\n  </AppendedData>\n");
    }

    file_content = sdscat(file_content, "</VTKFile>\n");

    size_until_now += 29;

    FILE *output_file = NULL;

    if(binary) {
        output_file = fopen(filename, "wb");
        fwrite(file_content, size_until_now, 1, output_file);
    } else {
        output_file = fopen(filename, "w");
        fprintf(output_file, "%s", file_content);
    }

    sdsfree(file_content);
    fclose(output_file);
}

//TODO: for non adaptive meshes we need to compress only the values data array. The other arrays remain the same. So we need to save them somewhere
void save_vtk_unstructured_grid_as_vtu_compressed(struct vtk_unstructured_grid *vtk_grid, const char *filename, int compression_level) {

    sds first_file_part = create_common_vtu_header(true, vtk_grid->num_points, vtk_grid->num_cells);

    size_t offset = 0;

    first_file_part = sdscat(
            first_file_part, "        <DataArray type=\"Float32\" Name=\"Scalars_\" format=\"appended\" offset=\"0\" />\n");

    first_file_part = sdscat(first_file_part, "      </CellData>\n");
    first_file_part = sdscat(first_file_part, "      <Points>\n");

    sds points_array_header;

    if(sizeof(real_cpu) == 4) {
        points_array_header = sdsnew("        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" "
                                         "format=\"appended\" offset=\"%zu\" />\n");
    }
    else {
        points_array_header = sdsnew("        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" "
                                     "format=\"appended\" offset=\"%zu\" />\n");
    }

    sds points_array_header_end = sdsnew("      </Points>\n");
    points_array_header_end = sdscat(points_array_header_end, "      <Cells>\n");

    sds connectivity_array_header =
            sdsnew("        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"appended\" offset=\"%zu\"/>\n");

    sds offsets_array_header =
            sdsnew("        <DataArray type=\"Int64\" Name=\"offsets\" format=\"appended\" offset=\"%zu\"/>\n");

    sds types_array_header =
            sdsnew("        <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\"%zu\"/>\n");


    sds data_end = sdsnew("      </Cells>\n");

    data_end = sdscat(data_end, "    </Piece>\n");
    data_end = sdscat(data_end, "  </UnstructuredGrid>\n");

    sds appended_begin = sdsnew("  <AppendedData encoding=\"raw\">\n   _");

    FILE *output_file = NULL;

    output_file = fopen(filename, "wb");

    fwrite(first_file_part, sdslen(first_file_part), 1, output_file);
    sdsfree(first_file_part);

    {
        size_t data_size;

        size_t data_size_after_compression_for_values;
        size_t num_block_for_values;
        size_t block_size_uncompressed_for_values;
        size_t *block_sizes_compressed_for_values;
        size_t last_block_size_for_values;

        unsigned char *data_to_compress;

        data_size = vtk_grid->num_cells * 4; // 32 bit float for each cell

        data_to_compress = (unsigned char *)vtk_grid->values;
        unsigned char *compressed_data_for_values;

        calculate_blocks_and_compress_data(data_size, &data_size_after_compression_for_values, data_to_compress,
                                           &(compressed_data_for_values), &num_block_for_values,
                                           &block_size_uncompressed_for_values, &block_sizes_compressed_for_values,
                                           &last_block_size_for_values, compression_level);

        offset = 8 + 8 + 8 + (num_block_for_values * 8) + data_size_after_compression_for_values;

        sds tmp = sdscatprintf(sdsempty(), points_array_header, offset);
        fwrite(tmp, sdslen(tmp), 1, output_file);
        sdsfree(tmp);
        sdsfree(points_array_header);

        fwrite(points_array_header_end, sdslen(points_array_header_end), 1, output_file);
        sdsfree(points_array_header_end);

        size_t data_size_after_compression_for_points;
        size_t num_block_for_points;
        size_t block_size_uncompressed_for_points;
        size_t *block_sizes_compressed_for_points;
        size_t last_block_size_for_points;

        data_size = vtk_grid->num_points * sizeof(real_cpu) * 3; // 3 points, with 64 bit float each point

        data_to_compress = (unsigned char *)vtk_grid->points;

        unsigned char *compressed_data_for_points;

        calculate_blocks_and_compress_data(data_size, &data_size_after_compression_for_points, data_to_compress,
                                           &(compressed_data_for_points), &num_block_for_points,
                                           &block_size_uncompressed_for_points, &block_sizes_compressed_for_points,
                                           &last_block_size_for_points, compression_level);

        offset += 8 + 8 + 8 + (num_block_for_points * 8) + data_size_after_compression_for_points;

        tmp = sdscatprintf(sdsempty(), connectivity_array_header, offset);
        fwrite(tmp, sdslen(tmp), 1, output_file);
        sdsfree(tmp);
        sdsfree(connectivity_array_header);

        size_t data_size_after_compression_for_connectivity;
        size_t num_block_for_connectivity;
        size_t block_size_uncompressed_for_connectivity;
        size_t *block_sizes_compressed_for_connectivity;
        size_t last_block_size_for_connectivity;

        sds aux_data = sdsempty();
        for(int i = 0; i < vtk_grid->num_cells; i++) {
            for(int j = 0; j < vtk_grid->points_per_cell; j++) {
                int64_t aux = (int64_t)vtk_grid->cells[vtk_grid->points_per_cell * i + j];
                aux_data = sdscatlen(aux_data, &aux, sizeof(int64_t));
            }
        }

        data_size = vtk_grid->num_cells * vtk_grid->points_per_cell * sizeof(int64_t);

        data_to_compress = (unsigned char *)aux_data;
        unsigned char *compressed_data_for_connectivity;

        calculate_blocks_and_compress_data(data_size, &data_size_after_compression_for_connectivity, data_to_compress,
                                           &compressed_data_for_connectivity, &num_block_for_connectivity,
                                           &block_size_uncompressed_for_connectivity,
                                           &block_sizes_compressed_for_connectivity, &last_block_size_for_connectivity, compression_level);

        offset += 8 + 8 + 8 + (num_block_for_connectivity * 8) + data_size_after_compression_for_connectivity;

        tmp = sdscatprintf(sdsempty(), offsets_array_header, offset);
        fwrite(tmp, sdslen(tmp), 1, output_file);
        sdsfree(tmp);
        sdsfree(offsets_array_header);

        sdsfree(aux_data);

        // offsets
        size_t data_size_after_compression_for_offsets;
        size_t num_block_for_offsets;
        size_t block_size_uncompressed_for_offsets;
        size_t *block_sizes_compressed_for_offsets;
        size_t last_block_size_for_offsets;

        aux_data = sdsempty();
        int64_t offset_local = vtk_grid->points_per_cell;
        for(int i = 0; i < vtk_grid->num_cells; i++) {
            aux_data = sdscatlen(aux_data, &offset_local, sizeof(int64_t));
            offset_local += vtk_grid->points_per_cell;
        }

        data_size = vtk_grid->num_cells * sizeof(int64_t);

        data_to_compress = (unsigned char *)aux_data;
        unsigned char *compressed_data_for_offsets;

        calculate_blocks_and_compress_data(data_size, &data_size_after_compression_for_offsets, data_to_compress,
                                           &compressed_data_for_offsets, &num_block_for_offsets,
                                           &block_size_uncompressed_for_offsets, &block_sizes_compressed_for_offsets,
                                           &last_block_size_for_offsets, compression_level);

        offset += 8 + 8 + 8 + (num_block_for_offsets * 8) + data_size_after_compression_for_offsets;

        tmp = sdscatprintf(sdsempty(), types_array_header, offset);
        fwrite(tmp, sdslen(tmp), 1, output_file);
        sdsfree(tmp);
        sdsfree(types_array_header);

        sdsfree(aux_data);

        // types
        size_t data_size_after_compression_for_types;
        size_t num_block_for_types;
        size_t block_size_uncompressed_for_types;
        size_t *block_sizes_compressed_for_types;
        size_t last_block_size_for_types;

        aux_data = sdsempty();
        for(int i = 0; i < vtk_grid->num_cells; i++) {
            uint8_t aux = (uint8_t)vtk_grid->cell_type;
            aux_data = sdscatlen(aux_data, &aux, sizeof(uint8_t));
        }

        data_size = vtk_grid->num_cells * sizeof(uint8_t);

        data_to_compress = (unsigned char *)aux_data;
        unsigned char *compressed_data_for_types;

        calculate_blocks_and_compress_data(data_size, &data_size_after_compression_for_types, data_to_compress,
                                           &compressed_data_for_types, &num_block_for_types,
                                           &block_size_uncompressed_for_types, &block_sizes_compressed_for_types,
                                           &last_block_size_for_types, compression_level);

        sdsfree(aux_data);

        fwrite(data_end, sdslen(data_end), 1, output_file);
        sdsfree(data_end);

        fwrite(appended_begin, sdslen(appended_begin), 1, output_file);
        sdsfree(appended_begin);

        // Now we can save the compressed data
        //Values
        fwrite(&num_block_for_values, sizeof(uint64_t), 1, output_file);
        fwrite(&block_size_uncompressed_for_values, sizeof(uint64_t), 1, output_file);
        fwrite(&last_block_size_for_values, sizeof(uint64_t), 1, output_file);
        fwrite(block_sizes_compressed_for_values, sizeof(uint64_t), num_block_for_values, output_file);
        fwrite(compressed_data_for_values, data_size_after_compression_for_values, 1, output_file);
        free(compressed_data_for_values);
        free(block_sizes_compressed_for_values);

        //Points
        fwrite(&num_block_for_points, sizeof(uint64_t), 1, output_file);
        fwrite(&block_size_uncompressed_for_points, sizeof(uint64_t), 1, output_file);
        fwrite(&last_block_size_for_points, sizeof(uint64_t), 1, output_file);
        fwrite(block_sizes_compressed_for_points, sizeof(uint64_t), num_block_for_points, output_file);
        fwrite(compressed_data_for_points, data_size_after_compression_for_points, 1, output_file);

        free(compressed_data_for_points);
        free(block_sizes_compressed_for_points);

        //connectivity
        fwrite(&num_block_for_connectivity, sizeof(uint64_t), 1, output_file);
        fwrite(&block_size_uncompressed_for_connectivity, sizeof(uint64_t), 1, output_file);
        fwrite(&last_block_size_for_connectivity, sizeof(uint64_t), 1, output_file);
        fwrite(block_sizes_compressed_for_connectivity, sizeof(uint64_t), num_block_for_connectivity, output_file);
        fwrite(compressed_data_for_connectivity, data_size_after_compression_for_connectivity, 1, output_file);

        free(compressed_data_for_connectivity);
        free(block_sizes_compressed_for_connectivity);

        //offsets
        fwrite(&num_block_for_offsets, sizeof(uint64_t), 1, output_file);
        fwrite(&block_size_uncompressed_for_offsets, sizeof(uint64_t), 1, output_file);
        fwrite(&last_block_size_for_offsets, sizeof(uint64_t), 1, output_file);
        fwrite(block_sizes_compressed_for_offsets, sizeof(uint64_t), num_block_for_offsets, output_file);
        fwrite(compressed_data_for_offsets, data_size_after_compression_for_offsets, 1, output_file);

        free(compressed_data_for_offsets);
        free(block_sizes_compressed_for_offsets);

        //types
        fwrite(&num_block_for_types, sizeof(uint64_t), 1, output_file);
        fwrite(&block_size_uncompressed_for_types, sizeof(uint64_t), 1, output_file);
        fwrite(&last_block_size_for_types, sizeof(uint64_t), 1, output_file);
        fwrite(block_sizes_compressed_for_types, sizeof(uint64_t), num_block_for_types, output_file);
        fwrite(compressed_data_for_types, data_size_after_compression_for_types, 1, output_file);

        free(compressed_data_for_types);
        free(block_sizes_compressed_for_types);
    }

    sds appended_end = sdsnew("\n  </AppendedData>\n");
    appended_end = sdscat(appended_end, "</VTKFile>\n");

    fwrite(appended_end, sdslen(appended_end), 1, output_file);
    sdsfree(appended_end);

    fclose(output_file);
}

void save_vtk_unstructured_grid_as_legacy_vtk(struct vtk_unstructured_grid *vtk_grid, char *filename, bool binary) {

    sds file_content = sdsempty();

    file_content = sdscat(file_content, "# vtk DataFile Version 4.2\n");
    file_content = sdscat(file_content, "vtk output\n");

    if(binary) {
        file_content = sdscat(file_content, "BINARY\n");
    } else {
        file_content = sdscat(file_content, "ASCII\n");
    }

    file_content = sdscat(file_content, "DATASET UNSTRUCTURED_GRID\n");
    file_content = sdscatprintf(file_content, "POINTS %d float\n", vtk_grid->num_points);

    size_t size_until_now = sdslen(file_content);

    size_t num_points = arrlenu(vtk_grid->points);
    for(size_t i = 0; i < num_points; i++) {
        struct point_3d p = vtk_grid->points[i];
        if(binary) {
            file_content = write_binary_point(file_content, &p);
            size_until_now += 3 * sizeof(int);
        } else {
            file_content = sdscatprintf(file_content, "%lf %lf %lf\n", p.x, p.y, p.z);
        }
    }

    int num_cells = vtk_grid->num_cells;

    {
        sds tmp = sdscatprintf(sdsempty(), "\nCELLS %d %d\n", num_cells, 9 * num_cells);

        size_until_now += sdslen(tmp);

        file_content = sdscatsds(file_content, tmp);

        sdsfree(tmp);
    }

    int points_per_cell = 8;
    int cell_type = 12;

    int points_per_cell_swapped = invert_bytes(points_per_cell);
    int cell_type_swapped = invert_bytes(cell_type);

    for(int i = 0; i < num_cells; i++) {
        if(binary) {
            file_content = sdscatlen(file_content, &points_per_cell_swapped, sizeof(int));
            size_until_now += sizeof(int);

        } else {
            file_content = sdscatprintf(file_content, "%d ", points_per_cell);
        }

        for(int j = 0; j < points_per_cell; j++) {
            if(binary) {
                int aux = invert_bytes(vtk_grid->cells[points_per_cell * i + j]);
                file_content = sdscatlen(file_content, &aux, sizeof(int));
                size_until_now += sizeof(int);
            } else {
                file_content = sdscatprintf(file_content, "%ld ", vtk_grid->cells[points_per_cell * i + j]);
            }
        }

        if(!binary)
            file_content = sdscat(file_content, "\n");
    }

    {
        sds tmp = sdscatprintf(sdsempty(), "\nCELL_TYPES %d\n", num_cells);

        size_until_now += sdslen(tmp);

        file_content = sdscatsds(file_content, tmp);

        sdsfree(tmp);
    }

    for(int i = 0; i < num_cells; i++) {
        if(binary) {
            file_content = sdscatlen(file_content, &cell_type_swapped, sizeof(int));
            size_until_now += sizeof(int);
        } else {
            file_content = sdscatprintf(file_content, "%d\n", cell_type);
        }
    }

    {
        sds tmp = sdscatprintf(sdsempty(), "\nCELL_DATA %d\n", num_cells);
        tmp = sdscat(tmp, "SCALARS Scalars_ float\n");
        tmp = sdscat(tmp, "LOOKUP_TABLE default\n");

        size_until_now += sdslen(tmp);

        file_content = sdscatsds(file_content, tmp);
        sdsfree(tmp);
    }

    size_t num_values = arrlenu(vtk_grid->values);

    for(size_t i = 0; i < num_values; i++) {
        if(binary) {
            int aux = invert_bytes(*((int *)&(vtk_grid->values[i])));
            file_content = sdscatlen(file_content, &aux, sizeof(int));
            size_until_now += sizeof(int);
        } else {
            file_content = sdscatprintf(file_content, "%lf ", vtk_grid->values[i]);
        }
    }

    FILE *output_file = fopen(filename, "w");
    if(binary) {
        fwrite(file_content, size_until_now, 1, output_file);
    } else {
        fprintf(output_file, "%s", file_content);
    }

    sdsfree(file_content);
    fclose(output_file);
}

void save_vtk_unstructured_grid_as_alg_file(struct vtk_unstructured_grid *vtk_grid, char *filename, bool binary) {

    sds file_content = sdsempty();

    int64_t *cells = vtk_grid->cells;
    point3d_array points = vtk_grid->points;

    uint32_t n_active = vtk_grid->num_cells;

    int num_points = vtk_grid->points_per_cell;
    int j = num_points;

    for (uint32_t i = 0; i < n_active*num_points; i+=num_points) {

        float mesh_center_x, mesh_center_y, mesh_center_z;
        real_cpu v;
        float dx, dy, dz;

        dx = fabs((points[cells[i]].x - points[cells[i+1]].x));
        dy = fabs((points[cells[i]].y - points[cells[i+3]].y));
        dz = fabs((points[cells[i]].z - points[cells[i+4]].z));

        mesh_center_x = points[cells[i]].x + dx/2.0f;
        mesh_center_y = points[cells[i]].y + dy/2.0f;
        mesh_center_z = points[cells[i]].z + dz/2.0f;

        v = vtk_grid->values[j-num_points];
        j += 1;

        file_content = sdscatprintf(file_content, "%g,%g,%g,%g,%g,%g,%g\n", mesh_center_x, mesh_center_y, mesh_center_z, dx/2.0f, dy/2.0f, dz/2.0f, v);

    }
    FILE *output_file = fopen(filename, "w");
    if(binary) {
        fwrite(file_content, sdslen(file_content), 1, output_file);
    } else {
        fprintf(output_file, "%s", file_content);
    }

    sdsfree(file_content);
    fclose(output_file);
}

static int parse_vtk_legacy(char *source, size_t source_size, struct parser_state *state) {

    //ignoring the first two lines...
    while (*source != '\n') { source++; source_size--;}
    source++; source_size--;

    while (*source != '\n') {source++; source_size--;}
    source++; source_size--;

    char *type = NULL;
    char *data_name = NULL;

    unsigned long num_cells = 0;

    while (!isspace(*source)) {
        arrput(type, *source);
        source++; source_size--;
    }

    arrput(type, '\0');
    source++; source_size--;

    bool binary = strcasecmp(type, "BINARY") == 0;

    state->binary = binary;
    state->ascii = !binary;

    //ignoring DATASET line as we only handle UNSTRUCTURED_GRID for now....
    while (*source != '\n') { source++; source_size--;}
    source++; source_size--;

    while (source_size > 0) {

        while (!isspace(*source)) {
            arrput(data_name, *source);
            source++; source_size--;
        }

        arrput(data_name, '\0');

        source++; //skip \n or space
        source_size--;

        if (strcasecmp(data_name, POINTS) == 0) {

            while (!isspace(*source)) {
                arrput(state->number_of_points, *source);
                source++; source_size--;
            }
            arrput(state->number_of_points, '\0');

            //ignoring the rest of the line as we only save as float
            while (*source != '\n') { source++; source_size--; }
            source++; source_size--;

            if(!binary) {
                while (isspace(*source) || isdigit(*source) || *source == '-' || *source == '.') {
                    arrput(state->points_ascii, *source);
                    source++; source_size--;
                }
                arrput(state->points_ascii, '\0');
            }
            else {
                unsigned long num_points = strtoul(state->number_of_points, NULL, 10);
                for(unsigned long i = 0; i < num_points; i++) {
                    struct point_3d p;
                    read_binary_point(source, &p);

                    for(int b = 0; b < 12; b++) {
                        arrput(state->points_ascii, *((char*)(&p) + b)); //here this array will be treated as binary data
                    }

                    source += 12;
                    source_size -= 12;
                }
                source++;source_size--;
            }

            data_name[0] = '\0';
            arrsetlen(data_name, 0);
        }
        else if (strcasecmp(data_name, CELLS) == 0) {

            while (!isspace(*source)) {
                arrput(state->number_of_cells, *source);
                source++;source_size--;
            }

            arrput(state->number_of_cells, '\0');
            source++;source_size--;

            while (*source != '\n') {source++; source_size--;}
            source++; source_size--;

            if(!binary) {
                bool add_next = false;
                while (isspace(*source) || isdigit(*source) || *source == '-' || *source == '.') {
                    if (*source == '\n') {
                        add_next = false;
                        source++; source_size--;
                    }
                    if (add_next) {
                        arrput(state->cells_connectivity_ascii, *source);
                    }
                    source++; source_size--;
                    add_next = true;
                }
            }
            else {
                num_cells = strtoul(state->number_of_cells, NULL, 10);

                for(unsigned long i = 0; i < num_cells; i++) {
                    int points_per_cell = invert_bytes(*(int*)source);
                    source += 4; source_size -= 4;

                    for (int c = 0; c < points_per_cell; c++) {                 
                        uint64_t cell_point = (uint64_t )invert_bytes(*(int*)source);
                        for(int b = 0; b < 8; b++) {
                            arrput(state->cells_connectivity_ascii, *((char*)(&cell_point) + b) ); //here this array will be treated as binary data
                        }

                        source += 4; source_size -= 4;
                    }
                }
                source++; source_size--;

            }

            data_name[0] = '\0';
            arrsetlen(data_name, 0);

        }
        else if (strcasecmp(data_name, CELL_TYPES) == 0) {
            while (!isspace(*source)) {
                //arrput(state->number_of_cells, *source);
                source++; source_size--;
            }

            //arrput(state->number_of_cells, '\0');
            source++; source_size--;

            if(!binary) {
                while (isspace(*source) || isdigit(*source) || *source == '-' || *source == '.') {
                    source++; source_size--;
                }
            }
            else {
                for(unsigned long i = 0; i < num_cells; i++) {
                    //int type = *(int*)source;
                    source += 4; source_size -= 4;
                }

                source++; source_size--;

            }

            data_name[0] = '\0';
            arrsetlen(data_name, 0);

        }
        else if (strcasecmp(data_name, CELL_DATA) == 0) {
            while (!isspace(*source)) {
                //arrput(state->number_of_cells, *source);
                source++; source_size--;
            }

            //arrput(state->number_of_cells, '\0');
            source++; source_size--;

            data_name[0] = '\0';
            arrsetlen(data_name, 0);

        }
        else if (strcasecmp(data_name, SCALARS) == 0) {

            while (*source != '\n') {source++; source_size--;}
            source++; source_size--;

            data_name[0] = '\0';
            arrsetlen(data_name, 0);

        }
        else if (strcasecmp(data_name, LOOKUP_TABLE) == 0) {

            while (*source != '\n') {source++; source_size--;}
            source++; source_size--;

            if(!binary) {
                while (isspace(*source) || isdigit(*source) || *source == '-' || *source == '.') {
                    arrput(state->celldata_ascii, *source);
                    source++; source_size--;
                }
            }
            else {
                for (unsigned long c = 0; c < num_cells; c++) {

                   int raw_value = invert_bytes(*(int*)source);
                   float float_value = *(float*) (&raw_value);

                    for(int b = 0; b < 4; b++) {
                        arrput(state->celldata_ascii, *((char*)(&float_value) + b) ); //here this array will be treated as binary data
                    }

                    source += 4; source_size -= 4;
                }
                if(source_size) {
                    source++;
                    source_size--;
                }
            }

            data_name[0] = '\0';
            arrsetlen(data_name, 0);

        }
        else {
            while (*source != '\n') { source++; source_size--;}
            source++; source_size--;

            data_name[0] = '\0';
            arrsetlen(data_name, 0);
        }
    }
    return 1;
}

static int parse_vtk_xml(yxml_t *x, yxml_ret_t r, struct parser_state *state) {

    switch(r) {
        case YXML_OK:
            break;
        case YXML_ELEMSTART:
            if(strcmp(DATAARRAY, x->elem) == 0) {
                state->in_dataarray = true;
            }

            break;
        case YXML_ELEMEND:
            state->in_dataarray = false;

            if(state->ascii) {
                if (strcmp(SCALARS_NAME, state->name_value) == 0) {
                    arrput(state->celldata_ascii, '\0');
                } else if (strcmp(POINTS, state->name_value) == 0) {
                    arrput(state->points_ascii, '\0');
                }

                state->name_value[0] = '\0';
                arrsetlen(state->name_value, 0);
            }


            break;
        case YXML_ATTRSTART:
            if (strcmp(COMPRESSOR, x->attr) == 0) {
                state->compressed = true;
            }
            break;
        case YXML_ATTREND:
            if(state->in_dataarray) {
                if (strcmp(FORMAT, x->attr) == 0) {
                    arrput(state->format, '\0');
                    if(strcmp(state->format, ASCII) == 0) state->ascii = true;
                }
                else if (strcmp(NAME, x->attr) == 0) {
                    arrput(state->name_value, '\0');
                }
                else if (arrlen(state->name_value)) {
                    if (strcmp(OFFSET, x->attr) == 0) {
                        if (strcmp(SCALARS_NAME, state->name_value) == 0) {
                            arrput(state->celldata_ofsset, '\0');
                            if(!state->compressed) state->binary = true;
                        }
                        else if (strcmp(POINTS, state->name_value) == 0) {
                            arrput(state->points_ofsset, '\0');
                        }
                        else if (strcmp(CONNECTIVITY, state->name_value) == 0) {
                            arrput(state->cells_connectivity_ofsset, '\0');
                        }
                        else if (strcmp(OFFSETS, state->name_value) == 0) {
                            arrput(state->cells_offsets_ofsset, '\0');
                        }
                        else if (strcmp(TYPES, state->name_value) == 0) {
                            arrput(state->cells_types_ofsset, '\0');
                        }

                        state->name_value[0] = '\0';
                        arrsetlen(state->name_value, 0);
                    }
                }
            }

            if(strcmp(NUMBER_OF_POINTS, x->attr) == 0) {
                arrput(state->number_of_points, '\0');
            }
            else if(strcmp(NUMBER_OF_CELLS, x->attr) == 0) {
                arrput(state->number_of_cells, '\0');
            }
            else if (strcmp(ENCODING, x->attr) == 0) {
                        arrput(state->encoding_type, '\0');
            }
            else if (strcmp(HEADER_TYPE, x->attr) == 0) {
                        arrput(state->header_type, '\0');
            }
            break;
        case YXML_PICONTENT:
        case YXML_CONTENT:
            if(state->ascii) {
                if(state->in_dataarray) {

                    if(x->data[0] == '\n') x->data[0] = ' ';
                    if (strcmp(SCALARS_NAME, state->name_value) == 0) {
                        if (x->data[0] == '.' || x->data[0] == '-' || isspace(x->data[0]) || isdigit(x->data[0])) {
                            arrput(state->celldata_ascii, x->data[0]);
                        }
                    } else if (strcmp(POINTS, state->name_value) == 0) {
                        if (x->data[0] == '.' || x->data[0] == '-' || isspace(x->data[0]) || isdigit(x->data[0])) {
                            arrput(state->points_ascii, x->data[0]);
                        }
                    }
                    else if (strcmp(CONNECTIVITY, state->name_value) == 0) {
                        if (x->data[0] == '.' || x->data[0] == '-' || isspace(x->data[0]) || isdigit(x->data[0])) {
                            arrput(state->cells_connectivity_ascii, x->data[0]);
                        }
                    }
//                else if (strcmp(OFFSETS, state->name_value) == 0) {
//                    if (isdigit(x->data[0])) {
//                        arrput(state->cells_offsets_ofsset, x->data[0]);
//                    }                        }
//                else if (strcmp(TYPES, state->name_value) == 0) {
//                    if (isdigit(x->data[0])) {
//                        arrput(state->cells_types_ofsset, x->data[0]);
//                    }
//                }
                }
            }
            else {
                if (strcmp(APPENDEDDATA, x->elem) == 0) {
                    if(strcmp(state->encoding_type, "raw") == 0) {
                        //We dont have a valid XML code anymore. So we have to
                        //return here
                        return -1;
                    }
                    else if (strcmp(state->encoding_type, "base64") == 0) {
                        //base64 is a valid xml content.
                        arrput(state->base64_content, x->data[0]);
                    }
                }
            }
            break;
        case YXML_ATTRVAL:
            if(state->in_dataarray) {
                if (strcmp(NAME, x->attr) == 0) {
                    arrput(state->name_value, x->data[0]);
                } else if(arrlen(state->name_value)) {
                    if (strcmp(FORMAT, x->attr) == 0) {
                        arrput(state->format, x->data[0]);
                    }
                    if (strcmp(OFFSET, x->attr) == 0) {
                        if (strcmp(SCALARS_NAME, state->name_value) == 0) {
                            if (isdigit(x->data[0])) {
                                arrput(state->celldata_ofsset, x->data[0]);
                            }
                        }
                        else if (strcmp(POINTS, state->name_value) == 0) {
                            if (isdigit(x->data[0])) {
                                arrput(state->points_ofsset, x->data[0]);
                            }
                        }
                        else if (strcmp(CONNECTIVITY, state->name_value) == 0) {
                            if (isdigit(x->data[0])) {
                                arrput(state->cells_connectivity_ofsset, x->data[0]);
                            }
                        }
                        else if (strcmp(OFFSETS, state->name_value) == 0) {
                            if (isdigit(x->data[0])) {
                                arrput(state->cells_offsets_ofsset, x->data[0]);
                            }                        }
                        else if (strcmp(TYPES, state->name_value) == 0) {
                            if (isdigit(x->data[0])) {
                                arrput(state->cells_types_ofsset, x->data[0]);
                            }
                        }
                    }
                }
            }

            if(strcmp(NUMBER_OF_POINTS, x->attr) == 0) {
                if(isdigit(x->data[0])) {
                    arrput(state->number_of_points, x->data[0]);
                }
            }
            else if (strcmp(NUMBER_OF_CELLS, x->attr) == 0) {
                if (isdigit(x->data[0])) {
                    arrput (state->number_of_cells, x->data[0]);
                }
            }
            else if (strcmp(ENCODING, x->attr) == 0) {
                arrput (state->encoding_type, x->data[0]);
            }
            else if (strcmp(HEADER_TYPE, x->attr) == 0) {
                arrput (state->header_type, x->data[0]);
            }

            break;
        case YXML_PISTART:
            break;
        case YXML_PIEND:
            break;
        default:
           fprintf(stderr, "Error on xml %s \n", x->attr);
          //  exit(0);
    }

    return 0;
}

static void free_parser_state(struct parser_state *parser_state) {

    arrfree(parser_state->number_of_points);
    arrfree(parser_state->number_of_cells);
    arrfree(parser_state->celldata_ofsset);
    arrfree(parser_state->points_ofsset);
    arrfree(parser_state->cells_connectivity_ofsset);
    arrfree(parser_state->cells_offsets_ofsset);
    arrfree(parser_state->cells_types_ofsset);
    arrfree(parser_state->name_value);
    arrfree(parser_state->cells_connectivity_ascii);
    arrfree(parser_state->points_ascii);
    arrfree(parser_state->celldata_ascii);
    arrfree(parser_state->encoding_type);
    arrfree(parser_state->header_type);
    arrfree(parser_state->format);
    arrfree(parser_state->base64_content);
    free(parser_state);

}

//TODO: implement read only values for non-adaptive meshes
void set_vtk_grid_from_file(struct vtk_unstructured_grid **vtk_grid, const char *file_name) {

    //TODO: this whole code is really convoluted. We can do better than this mess...
    size_t size;

    struct parser_state *parser_state = NULL;

    char *tmp = read_entire_file_with_mmap(file_name, &size);

    if(!tmp) {
        *vtk_grid = NULL;
        return;
    }

    char *source = tmp;

    //TODO: Improve the file type inference
    //TRYING TO INFER THE FILE TYPE//////////////////////
    bool legacy  = (source[0] == '#');

    bool plain_text = isdigit(source[0]);

    bool xml = (source[0] == '<');

    bool activation_info = (source[0] == '0' || source[0] == '1') && (source[1] == '\n');
    ///////////////////////////////////

    size_t base64_outlen = 0;

    if(legacy || xml ) {
        parser_state =  calloc(1, sizeof(struct parser_state));
    }

    if(activation_info) {
        new_vtk_unstructured_grid_from_string_with_activation_info(vtk_grid, &source[2], size-2);
    }
    else if(xml) {
        //VTK XML file
        static char stack[8*1024];
        yxml_t *x = (yxml_t *) malloc(sizeof(yxml_t));
        yxml_init(x, stack, sizeof(stack));

        size_t bytes_read = 0;

        for (size_t i = 0; i < size; i++) {
            yxml_ret_t r = yxml_parse(x, source[i]);
            bytes_read++;
            if (parse_vtk_xml(x, r, parser_state) == -1) {
                break;
            }
        }

        source = source + bytes_read;

        free(x);
    }
    else if(legacy) {
        //VTK legacy file
        if( parse_vtk_legacy(source, size, parser_state) == -1) {
            *vtk_grid = NULL;
            return;
        }
    }
    else {
        //Simple text or binary representation
        new_vtk_unstructured_grid_from_string(vtk_grid, source, size, !plain_text, false);
    }

    if(legacy || xml ) {
        *vtk_grid = new_vtk_unstructured_grid();

        (*vtk_grid)->num_points = (uint32_t) strtoul(parser_state->number_of_points, NULL, 10);
        (*vtk_grid)->num_cells  = (uint32_t) strtoul(parser_state->number_of_cells,  NULL, 10);

        arrsetcap((*vtk_grid)->values, (*vtk_grid)->num_cells);
        arrsetcap((*vtk_grid)->points, (*vtk_grid)->num_points);
        arrsetcap((*vtk_grid)->cells,  (*vtk_grid)->num_cells * (*vtk_grid)->points_per_cell);


        if (xml && (parser_state->compressed || parser_state->binary)) {

            //TODO: change this to read the data type in the XML
            uint64_t scalars_offset_value = strtoul(parser_state->celldata_ofsset, NULL, 10);
            uint64_t points_offset_value = strtoul(parser_state->points_ofsset, NULL, 10);
            uint64_t cells_offset_value = strtoul(parser_state->cells_connectivity_ofsset, NULL, 10);

            bool is_b64 = strcmp(parser_state->encoding_type, "base64") == 0;
            bool is_raw = strcmp(parser_state->encoding_type, "raw") == 0;

            assert(is_b64 != is_raw);

            size_t b64_size = 0;
            size_t bytes_read = 0;

            if(is_raw) {
                //We ignore the \n, spaces and _ before the real data.
                //That is how VTK works!!
                while (*source != '_') source++;
                source++;
            }
            else if(is_b64) {
                //We ignore the \n, spaces and _ before the real data.
                //That is how VTK works!!
                while (*(parser_state->base64_content) != '_') stbds_arrdel(parser_state->base64_content, 0);
                stbds_arrdel(parser_state->base64_content, 0);

                b64_size = arrlenu(parser_state->base64_content) - 1;

                while (isspace(parser_state->base64_content[b64_size])) {
                    stbds_arrdel(parser_state->base64_content, b64_size);
                    b64_size--;
                }

                b64_size = arrlenu(parser_state->base64_content);
            }
            else {
                fprintf(stderr, "%s encoding not implemented yet\n", parser_state->encoding_type);
                *vtk_grid = NULL;
                return;
            }

            size_t header_size = sizeof(uint64_t);

            if(strcmp(parser_state->header_type, "UInt32") == 0 ) {
                header_size = sizeof(uint32_t);
            }

            char *raw_data = NULL;
            uint64_t num_blocks, block_size_uncompressed, last_block_size, *block_sizes_compressed;
            size_t raw_data_after_blocks_offset = 0;

            char *data_tmp =  parser_state->base64_content + scalars_offset_value;

            if(is_raw) {
                raw_data = (source + scalars_offset_value);
            }
            else if(is_b64) {
                //TODO: maybe we don't need to allocate this amount of memory
                raw_data = malloc(b64_size);
                base64_outlen = base64_decode((unsigned char*) raw_data, data_tmp, b64_size, &bytes_read);
            }

            if (parser_state->compressed) {
                raw_data_after_blocks_offset = get_block_sizes_from_compressed_vtu_file(raw_data, header_size, &num_blocks, &block_size_uncompressed,
                                                                                       &last_block_size, &block_sizes_compressed);

                if(is_b64 && (raw_data_after_blocks_offset == base64_outlen)) { //We read only the header.. need to read the data
                    b64_size -= bytes_read;
                    base64_decode((unsigned char *) raw_data + base64_outlen, data_tmp + bytes_read, b64_size, &bytes_read);
                }


                get_data_block_from_compressed_vtu_file(raw_data + raw_data_after_blocks_offset, (*vtk_grid)->values, header_size, num_blocks, block_size_uncompressed,
                                                        last_block_size, block_sizes_compressed);
            }
            else {
                get_data_block_from_uncompressed_binary_vtu_file(raw_data, (*vtk_grid)->values, header_size);
            }

            if(is_raw) {
                raw_data = (source + points_offset_value);
            }
            else if(is_b64) {
                data_tmp =  parser_state->base64_content + points_offset_value;
                b64_size = arrlen(parser_state->base64_content) - points_offset_value;

                base64_outlen = base64_decode((unsigned char*) raw_data, data_tmp, b64_size, &bytes_read);

            }

            if (parser_state->compressed) {
                raw_data_after_blocks_offset = get_block_sizes_from_compressed_vtu_file(raw_data, header_size, &num_blocks, &block_size_uncompressed,
                                                                                        &last_block_size, &block_sizes_compressed);

                if(is_b64 && (raw_data_after_blocks_offset == base64_outlen)) { //We read only the header.. need to read the data
                    b64_size -= bytes_read;
                    base64_decode((unsigned char *) raw_data + base64_outlen, data_tmp + bytes_read, b64_size,  &bytes_read);
                }

                get_data_block_from_compressed_vtu_file(raw_data + raw_data_after_blocks_offset, (*vtk_grid)->points, header_size, num_blocks, block_size_uncompressed,
                                                        last_block_size, block_sizes_compressed);
            }
            else {
                get_data_block_from_uncompressed_binary_vtu_file(raw_data, (*vtk_grid)->points, header_size);
            }

            if(is_raw) {
                raw_data = (source + cells_offset_value);
            }
            else if(is_b64) {

                data_tmp =  parser_state->base64_content + cells_offset_value;
                b64_size = arrlen(parser_state->base64_content) - cells_offset_value;

                base64_outlen = base64_decode((unsigned char*) raw_data, data_tmp, b64_size, &bytes_read);
            }
            if (parser_state->compressed) {
                raw_data_after_blocks_offset = get_block_sizes_from_compressed_vtu_file(raw_data, header_size, &num_blocks, &block_size_uncompressed,
                                                                                        &last_block_size, &block_sizes_compressed);

                if(is_b64 && (raw_data_after_blocks_offset == base64_outlen)) { //We read only the header.. need to read the data
                    b64_size -= bytes_read;
                    base64_decode((unsigned char *) raw_data + base64_outlen,
                                  data_tmp + bytes_read, b64_size,
                                  &bytes_read);
                }

                get_data_block_from_compressed_vtu_file(raw_data + raw_data_after_blocks_offset, (*vtk_grid)->cells, header_size, num_blocks, block_size_uncompressed,
                                                        last_block_size, block_sizes_compressed);
            }
            else
                get_data_block_from_uncompressed_binary_vtu_file(raw_data, (*vtk_grid)->cells, header_size);

            if(base64_outlen) {
                free(raw_data);
            }

        }
        else if (legacy && parser_state->binary) {
            memcpy((*vtk_grid)->points, parser_state->points_ascii, (*vtk_grid)->num_points * sizeof(struct point_3d));
            memcpy((*vtk_grid)->cells, parser_state->cells_connectivity_ascii,
                   (*vtk_grid)->num_cells * (*vtk_grid)->points_per_cell * sizeof(uint64_t));
            memcpy((*vtk_grid)->values, parser_state->celldata_ascii, (*vtk_grid)->num_cells * sizeof(float));
        }
        else if (parser_state->ascii) {
            char *tmp_data = parser_state->celldata_ascii;
            char *pEnd;

            while (*tmp_data != '-' && !isdigit(*tmp_data)) tmp_data++;

            for (int i = 0; i < (*vtk_grid)->num_cells; i++) {
                arrput((*vtk_grid)->values, strtof(tmp_data, &pEnd));
                tmp_data = pEnd;
            }

            tmp_data = parser_state->points_ascii;

            while (*tmp_data != '-' && !isdigit(*tmp_data)) tmp_data++;

            for (int i = 0; i < (*vtk_grid)->num_points; i++) {
                struct point_3d p = (struct point_3d) {0.0, 0.0, 0.0};
                for (int j = 0; j < (*vtk_grid)->points_per_cell; j++) {

                    switch (j) {
                        case 0:
                            p.x = strtof(tmp_data, &pEnd);
                            tmp_data = pEnd;
                            break;
                        case 1:
                            p.y = strtof(tmp_data, &pEnd);
                            tmp_data = pEnd;
                            break;
                        case 2:
                            p.z = strtof(tmp_data, &pEnd);
                            tmp_data = pEnd;
                            break;
                        default:
                            break;
                    }

                }
                arrput((*vtk_grid)->points, p);
            }

            tmp_data = parser_state->cells_connectivity_ascii;

            while (!isdigit(*tmp_data)) tmp_data++;

            for (uint32_t i = 0; i < (*vtk_grid)->num_cells * (*vtk_grid)->points_per_cell; i++) {
                arrput((*vtk_grid)->cells, strtol(tmp_data, &pEnd, 10));
                tmp_data = pEnd;
            }
        }

        free_parser_state(parser_state);

    }

    munmap(tmp, size);

}

struct vtk_unstructured_grid * new_vtk_unstructured_grid_from_file(const char *file_name) {
    struct vtk_unstructured_grid *vtk_grid = NULL;
    set_vtk_grid_from_file(&vtk_grid, file_name);
    return vtk_grid;
}
