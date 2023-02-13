//
// Created by sachetto on 30/10/18.
//

#include "ensight_grid.h"
#include "../3dparty/stb_ds.h"
#include "../domains_library/mesh_info_data.h"
#include "../logger/logger.h"
#include <math.h>
#include <stdint.h>
#include <sys/mman.h>
#include <float.h>

static inline void new_line(FILE *f, bool binary) {
    if(!binary) {
        fprintf(f, "\n");
    }
}

static void write_string(char *str, FILE *f, bool binary) {

    char buffer[81];
    strncpy(buffer, str, 80);

    if(binary) {
        fwrite(buffer, sizeof(char), 80, f);
    } else {
        buffer[79] = 0; //ensight files cannot have lines longer than 79 chars
        fprintf(f, "%s", buffer);
    }
}

static inline void write_int(int i, FILE *f, bool binary) {
    if(binary) {
        fwrite(&i, sizeof(int), 1, f);
    } else {
        fprintf(f, "%10d", i);
    }
}

static inline void write_float(float n, FILE *f, bool binary) {
    if(binary) {
        fwrite(&n, sizeof(float), 1, f);
    } else {
        fprintf(f, "%12.5e", n);
    }
}

void save_case_file(char *filename, uint64_t num_files, real_cpu dt, int print_rate, int num_state_vars) {

    FILE *case_file = fopen(filename, "w");

    fprintf(case_file, "FORMAT\n");
    fprintf(case_file, "type:\tensight gold\n\n");

    fprintf(case_file, "GEOMETRY\n");
    fprintf(case_file, "model:\t%s\n\n", "geometry.geo");

    int n_digits = log10(num_files*500) + 1;

    sds base_file_name = sdsnew("Vm.Esca");

    for(int i = 0; i < n_digits; i++) {
        base_file_name = sdscat(base_file_name, "*");
    }

    fprintf(case_file, "VARIABLE\n");
    fprintf(case_file, "scalar per element:\t1\tVm\t%s\n\n", base_file_name);

    sdsfree(base_file_name);

    if(num_state_vars) {
        base_file_name = sdsnew("%s%d.Esca");

        for(int i = 0; i < n_digits; i++) {
            base_file_name = sdscat(base_file_name, "*");
        }    
    }

    for(int i = 1; i <= num_state_vars; i++) {
        sds file_name = sdscatprintf(sdsempty(), base_file_name, "Sv", i);
        fprintf(case_file, "scalar per element:\t1\tSv%d\t%s\n", i, file_name);
        sdsfree(file_name);
    }

    if(num_state_vars) {
        sdsfree(base_file_name);
    }

    fprintf(case_file, "\nTIME\n");
    fprintf(case_file, "time set: 1 Vm\n");
    fprintf(case_file, "number of steps: %zu\n", num_files);

    fprintf(case_file, "filename start number: \t0\n");
    fprintf(case_file, "filename increment: \t1");

    fprintf(case_file, "\n");

    fprintf(case_file, "time values: ");
    for(int i = 0; i < num_files; i++) {
        fprintf(case_file, "%lf ", i*dt*print_rate);
        if((i+1) % 6 == 0) {
            fprintf(case_file, "\n");
        }
    }

    fprintf(case_file, "\n");
    fclose(case_file);
}

void save_en6_result_file(char *filename, struct grid *the_grid, bool binary) {

    FILE *result_file;

    if(binary) {
        result_file = fopen(filename, "wb");
    } else {
        result_file = fopen(filename, "w");
    }

    write_string("Per element Vm", result_file, binary);
    new_line(result_file, binary);

    int part_number = 1;

    if(the_grid->num_active_cells > 0) {
        write_string("part", result_file, binary);
        new_line(result_file, binary);

        write_int(part_number, result_file, binary);
        new_line(result_file, binary);

        write_string("hexa8", result_file, binary);
        new_line(result_file, binary);

        for(int i = 0 ; i < the_grid->num_active_cells; i++) {
            float v = (float) the_grid->active_cells[i]->v;
            write_float(v, result_file, binary);
            new_line(result_file, binary);

        }
        part_number++;
    }

    if(the_grid->purkinje) {
        write_string("part", result_file, binary);
        new_line(result_file, binary);

        write_int(part_number, result_file, binary);
        new_line(result_file, binary);

        write_string("bar2", result_file, binary);
        new_line(result_file, binary);

        for(int i = 0 ; i < the_grid->purkinje->number_of_purkinje_cells - 1 ; i++) {
            float v = (float) the_grid->purkinje->purkinje_cells[i]->v;
            write_float(v, result_file, binary);
            new_line(result_file, binary);
        }
    }

    fclose(result_file);
}

void save_en6_result_file_state_vars(char *filename, real *sv_cpu, size_t num_cells, size_t num_sv_entries, int sv_entry, bool binary, bool gpu) {

    FILE *result_file;

    if(binary) {
        result_file = fopen(filename, "wb");
    } else {
        result_file = fopen(filename, "w");
    }

    if ( result_file == NULL )  {
        perror("Error occurred while opening file.\n");
        exit(1);
    }

    write_string("SV", result_file, binary);
    new_line(result_file, binary);

    int part_number = 1;

    if(num_cells > 0) {
        write_string("part", result_file, binary);
        new_line(result_file, binary);

        write_int(part_number, result_file, binary);
        new_line(result_file, binary);

        write_string("hexa8", result_file, binary);
        new_line(result_file, binary);

        for(int i = 0 ; i < num_cells; i++) {
            float value;
            if(gpu) {
                real *sv_start = sv_cpu + sv_entry*num_cells;
                value = (float) sv_start[i];
            }
            else {
                value = sv_cpu[i*num_sv_entries + sv_entry];
            }
            write_float(value, result_file, binary);
            new_line(result_file, binary);
        }
        part_number++;
    }

    fclose(result_file);

    //TODO: purkinje

}

struct ensight_grid * new_ensight_grid(uint32_t num_parts) {

    struct ensight_grid *grid = MALLOC_ONE_TYPE(struct ensight_grid);
    grid->parts = NULL;
    grid->num_parts = num_parts;

    arrsetlen(grid->parts, num_parts);

    for(int i = 0; i < grid->num_parts; i++) {
        grid->parts[i].cells = NULL;
        grid->parts[i].points = NULL;
        grid->parts[i].cell_visibility = NULL;
    }

    return grid;
}

void free_ensight_grid(struct ensight_grid *grid) {
    if(grid) {
        for(int i = 0; i < grid->num_parts; i++) {
            arrfree(grid->parts[i].cells);
            arrfree(grid->parts[i].points);
            arrfree(grid->parts[i].cell_visibility);
        }
        arrfree(grid->parts);
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
        arrput((*grid)->parts[part_number].points, point1);
        hmput(*hash, point1, *id);
        *id += 1;
    }

    if(point2_idx == -1) {
        arrput((*grid)->parts[part_number].points, point2);
        hmput(*hash, point2, *id);
        *id += 1;
    }

    if(point3_idx == -1) {
        hmput(*hash, point3, *id);
        arrput((*grid)->parts[part_number].points, point3);
        *id += 1;
    }

    if(point4_idx == -1) {
        hmput(*hash, point4, *id);
        arrput((*grid)->parts[part_number].points, point4);
        *id += 1;
    }

    if(point5_idx == -1) {
        arrput((*grid)->parts[part_number].points, point5);
        hmput(*hash, point5, *id);
        *id += 1;
    }

    if(point6_idx == -1) {
        arrput((*grid)->parts[part_number].points, point6);
        hmput(*hash, point6, *id);
        *id += 1;
    }

    if(point7_idx == -1) {
        arrput((*grid)->parts[part_number].points, point7);
        hmput(*hash, point7, *id);
        *id += 1;
    }

    if(point8_idx == -1) {
        arrput((*grid)->parts[part_number].points, point8);
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
                                                                     bool save_fibrotic) {

    struct ensight_grid *ensight_grid;

    if(grid->num_active_cells > 0 && grid->purkinje) {

        uint32_t num_active_cells = grid->num_active_cells;
        uint32_t number_of_purkinje_cells = grid->purkinje->num_active_purkinje_cells;
        ensight_grid = new_ensight_grid(2);

        arrsetcap(ensight_grid->parts[0].points, num_active_cells * 8);
        arrsetcap(ensight_grid->parts[0].cell_visibility, num_active_cells);
        arrsetcap(ensight_grid->parts[0].cells,  num_active_cells);

        arrsetcap(ensight_grid->parts[1].points, number_of_purkinje_cells * 2);
        arrsetcap(ensight_grid->parts[1].cells,  number_of_purkinje_cells);

        ensight_grid->max_v = FLT_MIN;
        ensight_grid->min_v = FLT_MAX;

    } else {
        ensight_grid = new_ensight_grid(1);

        if(grid->num_active_cells > 0) {
            uint32_t num_active_cells = grid->num_active_cells;
            arrsetcap(ensight_grid->parts[0].points, num_active_cells * 8);
            arrsetcap(ensight_grid->parts[0].cell_visibility, num_active_cells);
            arrsetcap(ensight_grid->parts[0].cells,  num_active_cells);
        } else {
            uint32_t number_of_purkinje_cells = grid->purkinje->num_active_purkinje_cells;
            arrsetcap(ensight_grid->parts[0].points, number_of_purkinje_cells * 2);
            arrsetcap(ensight_grid->parts[0].cells,  number_of_purkinje_cells);
        }

    }

    ensight_grid->max_v = FLT_MIN;
    ensight_grid->min_v = FLT_MAX;
    struct point_hash_entry *hash =  NULL;

    if(grid->num_active_cells > 0) {
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

        struct point_3d half_face;
        struct point_3d center;

        real_cpu v;

        FOR_EACH_CELL(grid) {

            if(!cell->active) {
                if(!save_fibrotic) {
                    continue;
                } else if(cell->mesh_extra_info == NULL || !FIBROTIC(cell)) {
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

            if(v > ensight_grid->max_v) ensight_grid->max_v = (float) v;
            if(v < ensight_grid->min_v) ensight_grid->min_v = (float) v;


            half_face.x = cell->discretization.x / 2.0f;
            half_face.y = cell->discretization.y / 2.0f;
            half_face.z = cell->discretization.z / 2.0f;

            set_point_data(center, half_face, &ensight_grid, &hash, &id, 0);

            num_cells++;
        }

        ensight_grid->parts[0].num_cells = num_cells;
        ensight_grid->parts[0].num_points = id;
        ensight_grid->parts[0].part_description = "Mesh";
        ensight_grid->parts[0].cell_type = "hexa8";
        ensight_grid->parts[0].points_per_cell = 8;


        hmfree(hash);
        hash = NULL;
    }

    if(grid->purkinje) {

        int part_n = 0;

        if(grid->num_active_cells > 0) {
            part_n = 1;
        }

        hmdefault(hash, -1);
        struct cell_node *grid_cell = grid->purkinje->first_cell;
        struct node *u = grid->purkinje->network->list_nodes;
        struct point_3d aux;
        int id = 0;
        int num_cells = 0;


        struct point_hash_entry *lines =  NULL;

        while (grid_cell != NULL) {

            if (grid_cell->active) {

                // Insert the point to the array of points
                aux.x = grid_cell->center.x;
                aux.y = grid_cell->center.y;
                aux.z = grid_cell->center.z;

                // Search for duplicates
                if(hmget(hash, aux) == -1) {
                    arrput(ensight_grid->parts[part_n].points, aux);
                    hmput(hash, aux, id);
                    id++;
                }

                // Insert the edge to the array of lines
                struct edge *v = u->list_edges;
                while (v != NULL) {

                    struct point_3d p = POINT3D(u->id, v->id, 0);

                    if(hmgeti(lines, p) == -1) {
                        arrput(ensight_grid->parts[part_n].cells, u->id);
                        arrput(ensight_grid->parts[part_n].cells, v->id);
                        num_cells++;
                        struct point_3d p1 = POINT3D(u->id, v->id, 0);
                        struct point_3d p2 = POINT3D(v->id, u->id, 0);

                        hmput(lines, p1, 1);
                        hmput(lines, p2, 1);

                    }
                    v = v->next;
                }

            }
            grid_cell = grid_cell->next;
            u = u->next;
        }

        ensight_grid->parts[part_n].part_description = "Purkinje";
        ensight_grid->parts[part_n].cell_type = "bar2";
        ensight_grid->parts[part_n].points_per_cell = 2;
        ensight_grid->parts[part_n].num_points = id;
        ensight_grid->parts[part_n].num_cells = num_cells;

    }

    return ensight_grid;
}

void save_ensight_grid_as_ensight6_geometry(struct ensight_grid *grid, char *filename, bool binary) {

    FILE *output_file = NULL;

    if(binary) {
        output_file = fopen(filename, "wb");
        write_string("C Binary", output_file, binary);
    }
    else {
        output_file = fopen(filename, "w");
    }

    write_string("Grid", output_file, binary);
    new_line(output_file, binary);

    write_string("Grid geometry", output_file, binary);
    new_line(output_file, binary);

    write_string("node id off", output_file, binary);
    new_line(output_file, binary);

    write_string("element id off", output_file, binary);
    new_line(output_file, binary);

    for(int part = 0; part < grid->num_parts; part++) {

        int num_cells = (int)grid->parts[part].num_cells;
        int num_points = (int)grid->parts[part].num_points;
        int points_per_cell = grid->parts[part].points_per_cell;

        write_string("part", output_file, binary);
        new_line(output_file, binary);

        int part_n = part+1;
        write_int(part_n, output_file, binary);
        new_line(output_file, binary);

        write_string(grid->parts[part].part_description, output_file, binary);
        new_line(output_file, binary);

        write_string("coordinates", output_file, binary);
        new_line(output_file, binary);

        write_int(num_points, output_file, binary);
        new_line(output_file, binary);

        for(int i = 0; i < num_points; i++) {
            struct point_3d p = grid->parts[part].points[i];
            float tmp = (float)p.x;
            write_float(tmp, output_file, binary);
            new_line(output_file, binary);
        }

        for(int i = 0; i < num_points; i++) {
            struct point_3d p = grid->parts[part].points[i];
            float tmp = (float)p.y;
            write_float(tmp, output_file, binary);
            new_line(output_file, binary);
        }

        for(int i = 0; i < num_points; i++) {
            struct point_3d p = grid->parts[part].points[i];
            float tmp = (float)p.z;
            write_float(tmp, output_file, binary);
            new_line(output_file, binary);
        }

        write_string(grid->parts[part].cell_type, output_file, binary);
        new_line(output_file, binary);

        write_int(num_cells, output_file, binary);
        new_line(output_file, binary);

        for(int i = 0; i < num_cells; i++) {
            for(int j = 0; j < points_per_cell; j++) {
                int id = (int)grid->parts[part].cells[points_per_cell * i + j] + 1;
                write_int(id, output_file, binary);
            }
            new_line(output_file, binary);
        }
    }

    fclose(output_file);

}
