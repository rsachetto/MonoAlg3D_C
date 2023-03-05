//
// Created by sachetto on 30/10/18.
//

#ifndef MONOALG3D_VTK_UNSTRUCTURED_GRID_H
#define MONOALG3D_VTK_UNSTRUCTURED_GRID_H

#include "../alg/grid/grid.h"
#include "../common_types/common_types.h"

enum file_type_enum {
    VTK_LEGACY,
    VTU_XML,
    ALG_PLAIN_TEXT,
    ALG_BINARY,
    ENSIGHT_BINARY,
    ENSIGHT_ASCII,
    ACTIVATION
};

struct vtk_unstructured_grid {
    uint32_t num_points;
    uint32_t num_cells;

    //TODO: we handle only meshes with the same number of points per cell. If we need something different we will need to change this
    uint32_t points_per_cell;

    //TODO: we handle only meshes with the same cell_type. If we need something different we will need to change this
    uint8_t cell_type;

    f32_array values;
    f32_array *extra_values;
    real_cpu **fibers;
    int64_array cells;
    ui8_array cell_visibility;
    point3d_array points;

    float min_v;
    float max_v;

    f32_array min_extra_value;
    f32_array max_extra_value;

    struct point_3d average_discretization;

    //TODO: I don't know if this is the best place to put this information
    struct vtk_unstructured_grid *purkinje;

};

struct vtk_unstructured_grid *new_vtk_unstructured_grid();

void new_vtk_unstructured_grid_from_alg_grid (struct vtk_unstructured_grid **vtk_grid, struct grid *grid, bool clip_with_plain,
                                                                     float *plain_coordinates, bool clip_with_bounds,
                                                                     float *bounds, bool read_only_values, bool read_fibers_f,
                                                                     bool save_fibrotic, real *values);

void save_vtk_unstructured_grid_as_vtu(struct vtk_unstructured_grid *vtk_grid, char *filename, bool binary);
void save_vtk_unstructured_grid_as_vtu_compressed(struct vtk_unstructured_grid *vtk_grid, const char *filename, int compression_level);
void save_vtk_unstructured_grid_as_legacy_vtk(struct vtk_unstructured_grid *vtk_grid, char *filename, bool binary, bool save_f, struct string_voidp_hash_entry *extra_data_config);
void save_vtk_unstructured_grid_as_alg_file(struct vtk_unstructured_grid *vtk_grid, char *filename, bool binary);
void free_vtk_unstructured_grid(struct vtk_unstructured_grid *vtk_grid);

struct vtk_unstructured_grid * new_vtk_unstructured_grid_from_file(const char *vtu_file_name, bool calc_max_min);
struct vtk_unstructured_grid * new_vtk_unstructured_grid_from_file_with_progress(const char *file_name, bool calc_max_min, size_t *bytes_read, size_t *file_size);

void new_vtk_unstructured_grid_from_string_with_activation_info(struct vtk_unstructured_grid **vtk_grid, char* source, size_t source_size);
void set_vtk_grid_values_from_ensight_file(struct vtk_unstructured_grid *grid, const char *file_name);
void set_vtk_grid_visibility(struct vtk_unstructured_grid **vtk_grid);
void read_or_calc_visible_cells(struct vtk_unstructured_grid **vtk_grid, sds full_path);
void calc_visibility(struct vtk_unstructured_grid **vtk_grid, struct cell_hash_entry *cells, uint32_t num_cells);

#endif // MONOALG3D_VTK_UNSTRUCTURED_GRID_H
