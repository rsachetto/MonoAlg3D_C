//
// Created by bergolho on 10/09/19.
//

#ifndef MONOALG3D_GRID_PURKINJE_H
#define MONOALG3D_GRID_PURKINJE_H

#include "../cell/cell.h"
#include "../../graph/graph.h"
#include "../../common_types/common_types.h"

#include <stdlib.h>
#include <stdio.h>

struct grid_purkinje {

    struct point_3d cube_side_length;
    struct point_3d mesh_side_length;

    uint32_t number_of_purkinje_cells;      // Number of Purkinje cells of grid.
    uint32_t num_active_purkinje_cells;     // Number of active Purkinje cells of grid.

    struct cell_node **purkinje_cells;
    struct cell_node *first_cell;           // First Purkinje cell of grid.

    struct graph *network;                  // Purkinje network graph

};

struct grid_purkinje* new_grid_purkinje();

#endif //MONOALG3D_GRID_PURKINJE_H
