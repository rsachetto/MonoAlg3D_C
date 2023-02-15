//
// Created by sachetto on 30/10/18.
//

#ifndef MONOALG3D_ENSIGHT_GRID_H
#define MONOALG3D_ENSIGHT_GRID_H

#include "../alg/grid/grid.h"
#include "../common_types/common_types.h"

struct ensight_grid {

    uint32_t num_parts;

    float min_v;
    float max_v;

    struct part {
        char *part_description;
        char *cell_type;
        int points_per_cell;
        uint32_t num_cells;
        int64_array cells;
        ui8_array cell_visibility;
        uint32_t num_points;
        point3d_array points;

       } *parts;
};

struct ensight_grid * new_ensight_grid(uint32_t num_parts);

struct ensight_grid * new_ensight_grid_from_alg_grid(struct grid *grid, bool clip_with_plain,
                                                                     float *plain_coordinates, bool clip_with_bounds,
                                                                     float *bounds, bool read_fibers_f,
                                                                     bool save_fibrotic);

void save_ensight_grid_as_ensight6_geometry(struct ensight_grid *grid, char *filename, bool binary);

void free_ensight_grid(struct ensight_grid *ensight_grid);
void save_case_file(char *filename, uint64_t num_files, real_cpu dt, int print_rate, int num_state_var);
void save_en6_result_file(char *filename, struct grid *the_grid, bool binary);
void save_en6_result_file_state_vars(char *filename, real *sv_cpu, size_t num_cells, size_t num_sv_entries, int sv_entry, bool binary, bool gpu);
#endif //MONOALG3D_ENSIGHT_GRID_H
