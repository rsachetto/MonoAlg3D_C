//
// Created by sachetto on 30/10/18.
//

#include "ensight_grid.h"
#include "../3dparty/sds/sds.h"
#include "../3dparty/stb_ds.h"
#include "../utils/file_utils.h"
#include "../domains_library/mesh_info_data.h"
#include "../common_types/common_types.h"

#include <math.h>
#include <stdint.h>
#include <ctype.h>
#include <sys/mman.h>
#include <float.h>

void save_case_file(char *filename, uint64_t num_files, real_cpu dt, int print_rate) {

    FILE *case_file = fopen(filename, "w");

    fprintf(case_file, "FORMAT\n");
    fprintf(case_file, "type:   ensight\n\n");

    fprintf(case_file, "GEOMETRY\n");
    fprintf(case_file, "model:  %s\n\n", "geometry.geo");

    int n_digits = (num_files==0) ? 1 : log10(num_files) + 1;

    sds base_file_name = sdsnew("Vm.Esca");

    for(int i = 0; i < n_digits; i++) {
        base_file_name = sdscat(base_file_name, "*");
    }

    fprintf(case_file, "VARIABLE\n");
    fprintf(case_file, "scalar per element:     1    Vm  %s\n\n", base_file_name);

    sdsfree(base_file_name);

    fprintf(case_file, "TIME\n");
    fprintf(case_file, "time set: 1 Vm\n");
    fprintf(case_file, "number of steps: %zu\n", num_files);

    fprintf(case_file, "filename numbers: ");
    for(int i = 0; i < num_files; i++) {
        fprintf(case_file, "%d ", i);
    }

    fprintf(case_file, "\n");


    fprintf(case_file, "time values: ");
    for(int i = 0; i < num_files; i++) {
        fprintf(case_file, "%lf ", i*dt*print_rate);
    }

    fprintf(case_file, "\n");

    fclose(case_file);
}

struct ensight_grid * new_ensight_grid(uint32_t num_parts) {

    struct ensight_grid *grid = MALLOC_ONE_TYPE(struct ensight_grid);
    grid->parts = NULL;
    grid->points = NULL;
    grid->num_parts = num_parts;

    arrsetlen(grid->parts, num_parts);

    for(int i = 0; i < grid->num_parts; i++) {
        grid->parts[i].cells = NULL;
        grid->parts[i].cell_visibility = NULL;
    }


    //TODO: incomplete

    return grid;
}

void free_ensight_grid(struct ensight_grid *grid) {
    //TODO: incomplete
    if(grid) {
        for(int i = 0; i < grid->num_parts; i++) {
            arrfree(grid->parts[i].cells);
        }
        free(grid);
    }
}

static inline void set_point_data(struct point_3d center, struct point_3d half_face, struct ensight_grid **grid, struct point_hash_entry **hash, uint32_t *id, int part_number) {
    real_cpu center_x_plus  = center.x + half_face.x;
    real_cpu center_x_minus = center.x - half_face.x;

    real_cpu center_y_plus  = center.y + half_face.y;
    real_cpu center_y_minus = center.y - half_face.y;

    real_cpu center_z_plus  = center.z + half_face.z;
    real_cpu center_z_minus = center.z - half_face.z;

    struct point_3d points[8];

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

    struct point_3d  point1 = points[0];
    struct point_3d  point2 = points[1];
    struct point_3d  point3 = points[2];
    struct point_3d  point4 = points[3];
    struct point_3d  point5 = points[4];
    struct point_3d  point6 = points[5];
    struct point_3d  point7 = points[6];
    struct point_3d  point8 = points[7];

    int point1_idx = hmgeti(*hash, point1);
    int point2_idx = hmgeti(*hash, point2);
    int point3_idx = hmgeti(*hash, point3);
    int point4_idx = hmgeti(*hash, point4);
    int point5_idx = hmgeti(*hash, point5);
    int point6_idx = hmgeti(*hash, point6);
    int point7_idx = hmgeti(*hash, point7);
    int point8_idx = hmgeti(*hash, point8);

    if(point1_idx == -1) {
        arrput((*grid)->points, point1);
        hmput(*hash, point1, *id);
        *id += 1;
    }

    if(point2_idx == -1) {
        arrput((*grid)->points, point2);
        hmput(*hash, point2, *id);
        *id += 1;
    }

    if(point3_idx == -1) {
        hmput(*hash, point3, *id);
        arrput((*grid)->points, point3);
        *id += 1;
    }

    if(point4_idx == -1) {
        hmput(*hash, point4, *id);
        arrput((*grid)->points, point4);
        *id += 1;
    }

    if(point5_idx == -1) {
        arrput((*grid)->points, point5);
        hmput(*hash, point5, *id);
        *id += 1;
    }

    if(point6_idx == -1) {
        arrput((*grid)->points, point6);
        hmput(*hash, point6, *id);
        *id += 1;
    }

    if(point7_idx == -1) {
        arrput((*grid)->points, point7);
        hmput(*hash, point7, *id);
        *id += 1;
    }

    if(point8_idx == -1) {
        arrput((*grid)->points, point8);
        hmput(*hash, point8, *id);
        *id += 1;
    }

    arrput((*grid)->parts[part_number].cells, (point1_idx != -1) ? point1_idx : hmget(*hash, point1));
    arrput((*grid)->parts[part_number].cells, (point2_idx != -1) ? point2_idx : hmget(*hash, point2));
    arrput((*grid)->parts[part_number].cells, (point3_idx != -1) ? point3_idx : hmget(*hash, point3));
    arrput((*grid)->parts[part_number].cells, (point4_idx != -1) ? point4_idx : hmget(*hash, point4));
    arrput((*grid)->parts[part_number].cells, (point5_idx != -1) ? point5_idx : hmget(*hash, point5));
    arrput((*grid)->parts[part_number].cells, (point6_idx != -1) ? point6_idx : hmget(*hash, point6));
    arrput((*grid)->parts[part_number].cells, (point7_idx != -1) ? point7_idx : hmget(*hash, point7));
    arrput((*grid)->parts[part_number].cells, (point8_idx != -1) ? point8_idx : hmget(*hash, point8));

}

struct ensight_grid * new_ensight_grid_from_alg_grid(struct grid *grid, bool clip_with_plain,
                                                                     float *plain_coordinates, bool clip_with_bounds,
                                                                     float *bounds, bool read_fibers_f,
                                                                     bool save_fibrotic, bool save_purkinje) {

    uint32_t num_active_cells = grid->num_active_cells;

    struct ensight_grid *ensight_grid;

    if(save_purkinje) {
        uint32_t number_of_purkinje_cells = grid->purkinje->num_active_purkinje_cells;
        ensight_grid = new_ensight_grid(2);
        arrsetcap(ensight_grid->points, num_active_cells + number_of_purkinje_cells);
        arrsetcap(ensight_grid->parts[1].cells,  number_of_purkinje_cells);

        ensight_grid->parts[1].max_v = FLT_MIN;
        ensight_grid->parts[1].min_v = FLT_MAX;

    }
    else {
        ensight_grid = new_ensight_grid(1);
        arrsetcap(ensight_grid->points, num_active_cells);
    }

    arrsetcap(ensight_grid->parts[0].cell_visibility, num_active_cells);
    arrsetcap(ensight_grid->parts[0].cells,  num_active_cells);

    ensight_grid->parts[0].max_v = FLT_MIN;
    ensight_grid->parts[0].min_v = FLT_MAX;

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

    real_cpu v;

    //TODO: read purkinje grid

    FOR_EACH_CELL(grid) {

        if(!cell->active) {
            if(!save_fibrotic) {
                continue;
            }
            else if(cell->mesh_extra_info == NULL || !FIBROTIC(cell)) {
                continue;
            }

        }

        center = cell->center;
        v = cell->v;

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

        arrput(ensight_grid->parts[0].cell_visibility, cell->visible);

        if(v > ensight_grid->parts[0].max_v) ensight_grid->parts[0].max_v = v;
        if(v < ensight_grid->parts[0].min_v) ensight_grid->parts[0].min_v = v;


        half_face.x = cell->discretization.x / 2.0f;
        half_face.y = cell->discretization.y / 2.0f;
        half_face.z = cell->discretization.z / 2.0f;

        set_point_data(center, half_face, &ensight_grid, &hash, &id, 0);

        num_cells++;
    }

    ensight_grid->parts[0].num_cells = num_cells;
    ensight_grid->num_points = id;

    hmfree(hash);
    return ensight_grid;
}

void save_ensight_grid_as_ensight5_geometry(struct ensight_grid *grid, char *filename, bool binary, bool save_purkinje) {

    uint32_t num_cells = grid->parts[0].num_cells;
    uint32_t num_points = grid->num_points;

    sds file_content = sdsnew("Grid\nGrid geometry\n");
    file_content = sdscat(file_content, "node id given\n");
    file_content = sdscat(file_content, "element id given\n");
    file_content = sdscat(file_content, "coordinates\n");

    file_content = sdscatprintf(file_content, "    %8d\n", num_points);

    if(!binary) {
        for(int i = 0; i < num_points; i++) {
            struct point_3d p = grid->points[i];
            file_content = sdscatprintf(file_content, "    %8d %12.5e %12.5e %12.5e\n", i, p.x, p.y, p.z);
        }
    }

    file_content = sdscat(file_content, "part 1\nGrid\n");
    file_content = sdscatprintf(file_content, "hexa8\n");
    file_content = sdscatprintf(file_content, "%8d\n", num_cells);

    int points_per_cell = 8;
    //int cell_type = vtk_grid->cell_type;

    if(!binary) {
        for(int i = 0; i < num_cells; i++) {
            file_content = sdscatprintf(file_content, "%8d    ", i);
            for(int j = 0; j < points_per_cell; j++) {
                file_content = sdscatprintf(file_content, "%ld ", grid->parts[0].cells[points_per_cell * i + j]);
            }

            file_content = sdscat(file_content, "\n");
        }
    }

    file_content = sdscat(file_content, "\n");

    FILE *output_file = NULL;

    if(binary) {
        //output_file = fopen(filename, "wb");
        //fwrite(file_content, size_until_now, 1, output_file);
    } else {
        output_file = fopen(filename, "w");
        fprintf(output_file, "%s", file_content);
    }

    sdsfree(file_content);
    fclose(output_file);

}

/*
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

*/


struct ensight_grid * new_ensight_grid_from_file(const char *file_name) {
//    struct vtk_unstructured_grid *vtk_grid = NULL;
 //   set_vtk_grid_from_file(&vtk_grid, file_name, 6);
    return NULL;
}
