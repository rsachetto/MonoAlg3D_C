//
// Created by bergolho on 25/01/19.
//

#ifndef MONOALG3D_VTK_POLYDATA_GRID_H
#define MONOALG3D_VTK_POLYDATA_GRID_H

#include "../alg/grid/grid.h"
#include "../common_types/common_types.h"


struct vtk_polydata_grid
{
    uint32_t num_points;
    uint32_t num_lines;

    f32_array values;
    point3d_array points;
    line_array lines;
};

struct vtk_polydata_grid *new_vtk_polydata_grid ();

void new_vtk_polydata_grid_from_purkinje_grid(struct vtk_polydata_grid **vtk_grid, struct grid *grid, bool clip_with_plain,
                                                                     float *plain_coordinates, bool clip_with_bounds,
                                                                     float *bounds, bool read_only_values);

void save_vtk_polydata_grid_as_vtp (struct vtk_polydata_grid *vtk_grid, char *filename, bool binary);
void save_vtk_polydata_grid_as_vtp_compressed (struct vtk_polydata_grid *vtk_grid, char *filename, int compression_level);
void save_vtk_polydata_grid_as_legacy_vtk(struct vtk_polydata_grid *vtk_grid, char *filename, bool binary);


void free_vtk_polydata_grid(struct vtk_polydata_grid *vtk_grid);

#endif