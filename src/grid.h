//
// Created by sachetto on 29/09/17.
//

#ifndef MONOALG3D_GRID_H
#define MONOALG3D_GRID_H

#include "cell.h"

#include <stdlib.h>
#include <stdio.h>

struct grid {

    struct cell_node *first_cell;     // First cell of grid.
    float side_length;        // Length of cube grid. Default = 1.0.
    uint64_t number_of_cells;  // Number of cells of grid.

    uint64_t original_num_cells;
    uint64_t numActiveCells;

    //TODO: @Incomplete
    /*
    vector <int> freeSVPositions;

    int *cellsToSolve;
    vector<int> refinedThisStep;

    vector<CellNode*> activeCells;
     */

    bool init_ode;
    bool parallel;
    bool gpu;

};


void initialize_grid(struct grid *the_grid);
void construct_grid(struct grid *the_grid);
void print_grid(struct grid* the_grid, FILE *output_file);
void free_grid(struct grid *the_grid);

#endif //MONOALG3D_GRID_H
