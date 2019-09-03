//
// Created by sachetto on 30/10/18.
//

#ifndef MONOALG3D_VTK_UNSTRUCTURED_GRID_H
#define MONOALG3D_VTK_UNSTRUCTURED_GRID_H

#include "../alg/grid/grid.h"
#include "../common_types/common_types.h"

#define NUMBER_OF_POINTS "NumberOfPoints"
#define NUMBER_OF_CELLS  "NumberOfCells"
#define POINTS           "Points"
#define CELLS            "Cells"
#define CELL_TYPES       "Cell_Types"
#define CELL_DATA        "Cell_Data"
#define DATAARRAY        "DataArray"
#define SCALARS          "SCALARS"
#define OFFSET           "offset"
#define NAME             "Name"
#define SCALARS_NAME     "Scalars_"
#define OFFSETS          "offsets"
#define CONNECTIVITY     "connectivity"
#define TYPES            "types"
#define APPENDEDDATA     "AppendedData"
#define ENCODING         "encoding"
#define HEADER_TYPE      "header_type"
#define COMPRESSOR       "compressor"
#define FORMAT           "format"
#define ASCII            "ascii"
#define LOOKUP_TABLE     "LOOKUP_TABLE"


struct parser_state {
    char *number_of_points;
    char *number_of_cells;

    char *celldata_ofsset;
    char *points_ofsset;
    char *cells_connectivity_ofsset;
    char *cells_offsets_ofsset;
    char *cells_types_ofsset;
    char *name_value;
    char *format;
    char *celldata_ascii;
    char *celldata_int;
    char *points_ascii;
    char *cells_connectivity_ascii;
    char *encoding_type;
    char *header_type;
    char *base64_content;

    bool in_dataarray;
    bool compressed;
    bool binary;
    bool ascii;

};


struct vtk_unstructured_grid {
    uint32_t num_points;
    uint32_t num_cells;

    //TODO: we handle only meshes with the same number of points per cell. If we need something different we will need to change this
    uint32_t points_per_cell;

    //TODO: we handle only meshes with the same cell_type. If we need something different we will need to change this
    uint8_t cell_type;

    f32_array values;
    int64_array cells;
    point3d_array points;

};

struct vtk_unstructured_grid *new_vtk_unstructured_grid();

void new_vtk_unstructured_grid_from_alg_grid(struct vtk_unstructured_grid **vtk_grid, struct grid *grid, bool clip_with_plain,
                                             float *plain_coordinates, bool clip_with_bounds,
                                             float *bounds, bool read_only_values,\
                                             const char scalar_name);

void save_vtk_unstructured_grid_as_vtu(struct vtk_unstructured_grid *vtk_grid, char *filename, bool binary);
void save_vtk_unstructured_grid_as_vtu_compressed(struct vtk_unstructured_grid *vtk_grid, char *filename, int compression_level);
void save_vtk_unstructured_grid_as_legacy_vtk(struct vtk_unstructured_grid *vtk_grid, char *filename, bool binary);

struct vtk_unstructured_grid * new_vtk_unstructured_grid_from_vtu_file(const char *vtu_file_name);

void free_vtk_unstructured_grid(struct vtk_unstructured_grid *vtk_grid);

void new_vtk_unstructured_grid_from_string(struct vtk_unstructured_grid **vtk_grid, char* source, size_t source_size, bool binary, bool read_only_values);
void set_vtk_grid_from_file(struct vtk_unstructured_grid **vtk_grid, const char *vtu_file_name);


#endif // MONOALG3D_VTK_UNSTRUCTURED_GRID_H