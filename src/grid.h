//
// Created by sachetto on 29/09/17.
//

#ifndef MONOALG3D_GRID_H
#define MONOALG3D_GRID_H

#include "cell.h"

struct grid {

    struct cell_node *first_cell;     // First cell of grid.
    unsigned float side_length;        // Length of cube grid. Default = 1.0.
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

    uint8_t init_ode;

};



void initialize_grid(struct grid *the_grid) {

    the_grid->first_cell = 0;
    the_grid->side_length = 1.0;
    the_grid->number_of_cells = 1
    the_grid->init_ode = false;

}

void construct_grid(struct grid *the_grid) {

    uint8_t init_ode = the_grid->init_ode;
    unsigned float side_length = the_grid->side_length;

    // Cell nodes.
    struct cell_node *front_northeast_cell,
                     *front_northwest_cell,
                     *front_southeast_cell,
                     *front_southwest_cell,
                     *back_northeast_cell,
                     *back_northwest_cell,
                     *back_southeast_cell,
                     *back_southwest_cell;

    init_cell_node(front_northeast_cell, init_ode);
    init_cell_node(front_northwest_cell, init_ode);
    init_cell_node(front_southeast_cell, init_ode);
    init_cell_node(front_southwest_cell, init_ode);
    init_cell_node(back_northeast_cell, init_ode);
    init_cell_node(back_northwest_cell, init_ode);
    init_cell_node(back_southeast_cell, init_ode);
    init_cell_node(back_southwest_cell, init_ode);


    // Transition nodes.
    struct transition_node *north_transition_node,
                           *south_transition_node,
                           *east_transition_node,
                           *west_transition_node,
                           *front_transition_node,
                           *back_transition_node;

    init_transition_node(north_transition_node);
    init_transition_node(south_transition_node;
    init_transition_node(east_transition_node);
    init_transition_node(west_transition_node);
    init_transition_node(front_transition_node);
    init_transition_node(back_transition_node);


    float half_side_length    = side_length / 2.0;
    float quarter_side_length = side_length / 4.0;

    //__________________________________________________________________________
    //              Initialization of transition nodes.
    //__________________________________________________________________________
    // East transition node.
    east_transition_node->direction = 'e';
    east_transition_node->cell_data.center_x   = half_side_length;
    east_transition_node->cell_data.center_y   = side_length;
    east_transition_node->cell_data.center_z   = half_side_length;
    east_transition_node->single_connector     = NULL;
    east_transition_node->quadruple_connector1->cell_node = front_southeast_cell;
    east_transition_node->quadruple_connector2->cell_node = back_southeast_cell;
    east_transition_node->quadruple_connector3->cell_node = back_northeast_cell;
    east_transition_node->quadruple_connector4->cell_node = front_northeast_cell;

    // North transition node.
    north_transition_node->direction = 'n';
    north_transition_node->cell_data.center_x   = half_side_length;
    north_transition_node->cell_data->center_y   = half_side_length;
    north_transition_node->cell_data->center_z   = side_length;
    north_transition_node->single_connector     = NULL;
    north_transition_node->quadruple_connector1->cell_node = front_northwest_cell;
    north_transition_node->quadruple_connector2->cell_node = front_northeast_cell;
    north_transition_node->quadruple_connector3->cell_node = back_northeast_cell;
    north_transition_node->quadruple_connector4->cell_node = back_northwest_cell;

    // West transition node.
    west_transition_node->direction = 'w';
    west_transition_node->cell_data->center_x   = half_side_length;
    west_transition_node->cell_data->center_y   = 0.0;
    west_transition_node->cell_data->center_z   = half_side_length;
    west_transition_node->single_connector     = NULL;
    west_transition_node->quadruple_connector1->cell_node = front_southwest_cell;
    west_transition_node->quadruple_connector2->cell_node = back_southwest_cell;
    west_transition_node->quadruple_connector3->cell_node = back_northwest_cell;
    west_transition_node->quadruple_connector4->cell_node = front_northwest_cell;

    // South transition node.
    south_transition_node->direction = 's';
    south_transition_node->cell_data.center_x   = half_side_length;
    south_transition_node->cell_data->center_y   = half_side_length;
    south_transition_node->cell_data->center_z   = 0.0;
    south_transition_node->single_connector     = NULL;
    south_transition_node->quadruple_connector1->cell_node = front_southwest_cell;
    south_transition_node->quadruple_connector2->cell_node = front_southeast_cell;
    south_transition_node->quadruple_connector3->cell_node = back_southeast_cell;
    south_transition_node->quadruple_connector4->cell_node = back_southwest_cell;

    // Front transition node.
    front_transition_node->direction = 'f';
    front_transition_node->cell_data->center_x   = side_length;
    front_transition_node->cell_data->center_y   = half_side_length;
    front_transition_node->cell_data->center_z   = half_side_length;
    front_transition_node->single_connector     = NULL;
    front_transition_node->quadruple_connector1->cell_node = front_southwest_cell;
    front_transition_node->quadruple_connector2->cell_node = front_southeast_cell;
    front_transition_node->quadruple_connector3->cell_node = front_northeast_cell;
    front_transition_node->quadruple_connector4->cell_node = front_northwest_cell;

    // Back transition node.
    back_transition_node->direction = 'b';
    back_transition_node->cell_data->center_x   = 0;
    back_transition_node->cell_data->center_y   = half_side_length;
    back_transition_node->cell_data->center_z   = half_side_length;
    back_transition_node->single_connector     = NULL;
    back_transition_node->quadruple_connector1->cell_node = back_southwest_cell;
    back_transition_node->quadruple_connector2->cell_node = back_southeast_cell;
    back_transition_node->quadruple_connector3->cell_node = back_northeast_cell;
    back_transition_node->quadruple_connector4->cell_node = back_northwest_cell;

    //__________________________________________________________________________
    //                      Initialization of cell nodes.
    //__________________________________________________________________________
    // front Northeast subcell initialization.
    front_northeast_cell->cell_data->face_length     = half_side_length;
    front_northeast_cell->cell_data->half_face_length = quarter_side_length;
    front_northeast_cell->cell_data->bunch_number    = 1;

    front_northeast_cell->east     = east_transition_node;
    front_northeast_cell->north    = north_transition_node;
    front_northeast_cell->west     = front_northwest_cell;
    front_northeast_cell->south    = front_southeast_cell;
    front_northeast_cell->front    = front_transition_node;
    front_northeast_cell->back     = back_northeast_cell;
    front_northeast_cell->previous = 0;
    front_northeast_cell->next     = back_northeast_cell;

    front_northeast_cell->grid_position       = 0;
    front_northeast_cell->hilbert_shape_number = 1;

    front_northeast_cell->cell_data->center_x = half_side_length + quarter_side_length;
    front_northeast_cell->cell_data->center_y = half_side_length + quarter_side_length;
    front_northeast_cell->cell_data->center_z = half_side_length + quarter_side_length;

    // back Northeast subcell initialization.
    back_northeast_cell->face_length     = half_side_length;
    back_northeast_cell->half_face_length = quarter_side_length;
    back_northeast_cell->bunch_number    = 2;

    back_northeast_cell->east     = east_transition_node;
    back_northeast_cell->north    = north_transition_node;
    back_northeast_cell->west     = back_northwest_cell;
    back_northeast_cell->south    = back_southeast_cell;
    back_northeast_cell->front    = front_northeast_cell;
    back_northeast_cell->back     = back_transition_node;
    back_northeast_cell->previous = front_northeast_cell;
    back_northeast_cell->next     = back_northwest_cell;

    back_northeast_cell->grid_position       = 1;
    back_northeast_cell->hilbert_shape_number = 2;

    back_northeast_cell->cell_data->center_x = quarter_side_length;
    back_northeast_cell->cell_data->center_y = half_side_length + quarter_side_length;
    back_northeast_cell->cell_data->center_z = half_side_length + quarter_side_length;

    // back Northwest subcell initialization.
    back_northwest_cell->face_length     = half_side_length;
    back_northwest_cell->half_face_length = quarter_side_length;
    back_northwest_cell->bunch_number    = 3;

    back_northwest_cell->east     = back_northeast_cell;
    back_northwest_cell->north    = north_transition_node;
    back_northwest_cell->west     = west_transition_node;
    back_northwest_cell->south    = back_southwest_cell;
    back_northwest_cell->front    = front_northwest_cell;
    back_northwest_cell->back     = back_transition_node;
    back_northwest_cell->previous = back_northeast_cell;
    back_northwest_cell->next     = front_northwest_cell;

    back_northwest_cell->grid_position       = 2;
    back_northwest_cell->hilbert_shape_number = 2;

    back_northwest_cell->cell_data->center_x = quarter_side_length;
    back_northwest_cell->cell_data->center_y = quarter_side_length;
    back_northwest_cell->cell_data->center_z = half_side_length + quarter_side_length;

    // front Northwest subcell initialization.
    front_northwest_cell->face_length     = half_side_length;
    front_northwest_cell->half_face_length = quarter_side_length;
    front_northwest_cell->bunch_number    = 4;

    front_northwest_cell->east     = front_northeast_cell;
    front_northwest_cell->north    = north_transition_node;
    front_northwest_cell->west     = west_transition_node;
    front_northwest_cell->south    = front_southwest_cell;
    front_northwest_cell->front    = front_transition_node;
    front_northwest_cell->back     = back_northwest_cell;
    front_northwest_cell->previous = back_northwest_cell;
    front_northwest_cell->next     = front_southwest_cell;

    front_northwest_cell->grid_position       = 3;
    front_northwest_cell->hilbert_shape_number = 3;

    front_northwest_cell->cell_data->center_x = half_side_length + quarter_side_length;
    front_northwest_cell->cell_data->center_y = quarter_side_length;
    front_northwest_cell->cell_data->center_z = half_side_length + quarter_side_length;

    // front Southwest subcell initialization.
    front_southwest_cell->face_length     = half_side_length;
    front_southwest_cell->half_face_length = quarter_side_length;
    front_southwest_cell->bunch_number    = 5;

    front_southwest_cell->east     = front_southeast_cell;
    front_southwest_cell->north    = front_northwest_cell;
    front_southwest_cell->west     = west_transition_node;
    front_southwest_cell->south    = south_transition_node;
    front_southwest_cell->front    = front_transition_node;
    front_southwest_cell->back     = back_southwest_cell;
    front_southwest_cell->previous = front_northwest_cell;
    front_southwest_cell->next     = back_southwest_cell;

    front_southwest_cell->grid_position       = 4;
    front_southwest_cell->hilbert_shape_number = 3;

    front_southwest_cell->cell_data->center_x = half_side_length + quarter_side_length;
    front_southwest_cell->cell_data->center_y = quarter_side_length;
    front_southwest_cell->cell_data->center_z = quarter_side_length;

    // back Southwest subcell initialization.
    back_southwest_cell->face_length     = half_side_length;
    back_southwest_cell->half_face_length = quarter_side_length;
    back_southwest_cell->bunch_number    = 6;

    back_southwest_cell->east     = back_southeast_cell;
    back_southwest_cell->north    = back_northwest_cell;
    back_southwest_cell->west     = west_transition_node;
    back_southwest_cell->south    = south_transition_node;
    back_southwest_cell->front    = front_southwest_cell;
    back_southwest_cell->back     = back_transition_node;
    back_southwest_cell->previous = front_southwest_cell;
    back_southwest_cell->next     = back_southeast_cell;

    back_southwest_cell->grid_position       = 5;
    back_southwest_cell->hilbert_shape_number = 4;

    back_southwest_cell->cell_data->center_x = quarter_side_length;
    back_southwest_cell->cell_data->center_y = quarter_side_length;
    back_southwest_cell->cell_data->center_z = quarter_side_length;

    // back Southeast subcell initialization.
    back_southeast_cell->face_length     = half_side_length;
    back_southeast_cell->half_face_length = quarter_side_length;
    back_southeast_cell->bunch_number    = 7;

    back_southeast_cell->east     = east_transition_node;
    back_southeast_cell->north    = back_northeast_cell;
    back_southeast_cell->west     = back_southwest_cell;
    back_southeast_cell->south    = south_transition_node;
    back_southeast_cell->front    = front_southeast_cell;
    back_southeast_cell->back     = back_transition_node;
    back_southeast_cell->previous = back_southwest_cell;
    back_southeast_cell->next     = front_southeast_cell;

    back_southeast_cell->grid_position       = 6;
    back_southeast_cell->hilbert_shape_number = 4;

    back_southeast_cell->cell_data->center_x = quarter_side_length;
    back_southeast_cell->cell_data->center_y = half_side_length + quarter_side_length;
    back_southeast_cell->cell_data->center_z = quarter_side_length;

    // front Southeast subcell initialization.->center_z
    front_southeast_cell->face_length     = half_side_length;
    front_southeast_cell->half_face_length = quarter_side_length;
    front_southeast_cell->bunch_number    = 8;

    front_southeast_cell->east     = east_transition_node;
    front_southeast_cell->north    = front_northeast_cell;
    front_southeast_cell->west     = front_southwest_cell;
    front_southeast_cell->south    = south_transition_node;
    front_southeast_cell->front    = front_transition_node;
    front_southeast_cell->back     = back_southeast_cell;
    front_southeast_cell->previous = back_southeast_cell;
    front_southeast_cell->next     = 0;

    front_southeast_cell->grid_position       = 7;
    front_southeast_cell->hilbert_shape_number = 5;

    front_southeast_cell->cell_data->center_x = half_side_length + quarter_side_length;
    front_southeast_cell->cell_data->center_y = half_side_length + quarter_side_length;
    front_southeast_cell->cell_data->center_z = quarter_side_length;

    // Grid initialization
    firstCell = front_northeast_cell;
    numberOfCells = 8;
    
}

#endif //MONOALG3D_GRID_H
