//
// Created by sachetto on 30/10/18.
//

#ifndef MONOALG3D_VTK_UNSTRUCTURED_GRID_H
#define MONOALG3D_VTK_UNSTRUCTURED_GRID_H

#include "../alg/grid/grid.h"
#include "../hash/point_hash.h"

struct vtk_unstructured_grid {
    int num_points;
    int num_cells;
    float *values;
    int *cells;
    struct point_3d *points;
};

struct vtk_unstructured_grid *new_vtk_unstructured_grid();

struct vtk_unstructured_grid *new_vtk_unstructured_grid_from_alg_grid(struct grid *grid, bool clip_with_plain,
                                                                     float *plain_coordinates, bool clip_with_bounds,
                                                                     float *bounds);

void save_vtk_unstructured_grid_as_vtu(struct vtk_unstructured_grid *vtk_grid, char *filename, bool binary);
void save_vtk_unstructured_grid_as_legacy_vtk(struct vtk_unstructured_grid *vtk_grid, char *filename, bool binary);

void free_vtk_unstructured_grid(struct vtk_unstructured_grid *vtk_grid);

#endif // MONOALG3D_VTK_UNSTRUCTURED_GRID_H
