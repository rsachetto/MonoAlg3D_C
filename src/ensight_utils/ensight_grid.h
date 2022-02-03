//
// Created by sachetto on 30/10/18.
//

#ifndef MONOALG3D_ENSIGHT_GRID_H
#define MONOALG3D_ENSIGHT_GRID_H

#include "../alg/grid/grid.h"
#include "../common_types/common_types.h"

struct ensight_grid {

    char *description_line1;
    char *description_line2;

    uint32_t num_points;
    point3d_array points;

    uint32_t num_parts;

    struct part {
        char *part_description;
        uint32_t num_cells;
        char *cell_type;
        int64_array cells;
        ui8_array cell_visibility;
        float min_v;
        float max_v;
    } *parts;

    /*
    f32_array values;
    real_cpu **fibers;
    ui8_array cell_visibility;
*/
};

struct ensight_grid * new_ensight_grid(uint32_t num_parts);

struct ensight_grid * new_ensight_grid_from_alg_grid(struct grid *grid, bool clip_with_plain,
                                                                     float *plain_coordinates, bool clip_with_bounds,
                                                                     float *bounds, bool read_fibers_f,
                                                                     bool save_fibrotic, bool save_purkinje);

struct ensight_grid * new_ensight_grid_from_file(const char *ensight_file_name);
struct ensight_grid * new_ensight_grid_from_string(char* source, size_t source_size);

void save_ensight_grid_as_ensight6_geometry(struct ensight_grid *grid, char *filename, bool binary, bool save_purkinje);

void save_ensight_grid_as_alg_file(struct ensight_grid *grid, char *filename, bool binary);
void free_ensight_grid(struct ensight_grid *ensight_grid);
void save_case_file(char *filename, uint64_t num_files, real_cpu dt, int print_rate);
void save_en6_result_file(char *filename, struct grid *the_grid, bool binary);

#endif //MONOALG3D_ENSIGHT_GRID_H
