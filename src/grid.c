//
// Created by sachetto on 29/09/17.
//

#include "grid.h"

void initialize_grid(struct grid *the_grid) {

    the_grid->first_cell = 0;
    the_grid->side_length = 1.0;
    the_grid->number_of_cells = 1;
    the_grid->init_ode = false;

}

void construct_grid(struct grid *the_grid) {

    bool init_ode = the_grid->init_ode;
    float side_length = the_grid->side_length;

    size_t cell_node_size = sizeof(struct cell_node);
    size_t transition_node_size = sizeof(struct transition_node);

    // Cell nodes.
    struct cell_node
            *front_northeast_cell,
            *front_northwest_cell,
            *front_southeast_cell,
            *front_southwest_cell,
            *back_northeast_cell,
            *back_northwest_cell,
            *back_southeast_cell,
            *back_southwest_cell;


    front_northeast_cell = (struct cell_node *) malloc(cell_node_size);
    front_northwest_cell = (struct cell_node *) malloc(cell_node_size);
    front_southeast_cell = (struct cell_node *) malloc(cell_node_size);
    front_southwest_cell = (struct cell_node *) malloc(cell_node_size);
    back_northeast_cell = (struct cell_node *) malloc(cell_node_size);
    back_northwest_cell = (struct cell_node *) malloc(cell_node_size);
    back_southeast_cell = (struct cell_node *) malloc(cell_node_size);
    back_southwest_cell = (struct cell_node *) malloc(cell_node_size);


    init_cell_node(front_northeast_cell, init_ode);
    init_cell_node(front_northwest_cell, init_ode);
    init_cell_node(front_southeast_cell, init_ode);
    init_cell_node(front_southwest_cell, init_ode);
    init_cell_node(back_northeast_cell, init_ode);
    init_cell_node(back_northwest_cell, init_ode);
    init_cell_node(back_southeast_cell, init_ode);
    init_cell_node(back_southwest_cell, init_ode);


    // Transition nodes.
    struct transition_node
            *north_transition_node,
            *south_transition_node,
            *east_transition_node,
            *west_transition_node,
            *front_transition_node,
            *back_transition_node;

    north_transition_node = (struct transition_node *) malloc(transition_node_size);
    south_transition_node = (struct transition_node *) malloc(transition_node_size);
    east_transition_node = (struct transition_node *) malloc(transition_node_size);
    west_transition_node = (struct transition_node *) malloc(transition_node_size);
    front_transition_node = (struct transition_node *) malloc(transition_node_size);
    back_transition_node = (struct transition_node *) malloc(transition_node_size);

    init_transition_node(north_transition_node);
    init_transition_node(south_transition_node);
    init_transition_node(east_transition_node);
    init_transition_node(west_transition_node);
    init_transition_node(front_transition_node);
    init_transition_node(back_transition_node);


    float half_side_length = side_length / 2.0f;
    float quarter_side_length = side_length / 4.0f;
    //__________________________________________________________________________
    //              Initialization of transition nodes.
    //__________________________________________________________________________
    // East transition node.
    set_transition_node_data(east_transition_node, 'e', half_side_length, side_length, half_side_length, NULL,
                             front_southeast_cell, back_northeast_cell, front_northeast_cell, front_northwest_cell);


    set_transition_node_data(north_transition_node, 'n', half_side_length, half_side_length, side_length, NULL,
                             front_northwest_cell, front_northeast_cell, back_northeast_cell, back_northwest_cell);


    set_transition_node_data(west_transition_node, 'w', half_side_length, 0.0, half_side_length, NULL,
                             front_southwest_cell, back_southwest_cell, back_northwest_cell, front_northwest_cell);

    // South transition node.
    set_transition_node_data(south_transition_node, 's', half_side_length, half_side_length, 0.0, NULL,
                             front_southwest_cell, front_southeast_cell, back_southeast_cell, back_southwest_cell);

    // Front transition node.
    // South transition node.
    set_transition_node_data(front_transition_node, 'f', side_length, half_side_length, half_side_length, NULL,
                             front_southwest_cell, front_southeast_cell, front_northeast_cell, front_northwest_cell);

    // Back transition node.
    set_transition_node_data(back_transition_node, 'b', 0.0, half_side_length, half_side_length, NULL,
                             back_southwest_cell, back_southeast_cell, back_northeast_cell, back_northwest_cell);

    //__________________________________________________________________________
    //                      Initialization of cell nodes.
    //__________________________________________________________________________
    // front Northeast subcell initialization.
    set_cell_node_data(front_northeast_cell,
                       half_side_length,
                       quarter_side_length,
                       1,
                       east_transition_node,
                       north_transition_node,
                       front_northwest_cell,
                       front_southeast_cell,
                       front_transition_node,
                       back_northeast_cell,
                       NULL,
                       back_northeast_cell,
                       0,
                       1,
                       half_side_length + quarter_side_length,
                       half_side_length + quarter_side_length,
                       half_side_length + quarter_side_length);

    // back Northeast subcell initialization.
    set_cell_node_data(back_northeast_cell,
                       half_side_length,
                       quarter_side_length,
                       2,
                       east_transition_node,
                       north_transition_node,
                       back_northwest_cell,
                       back_southeast_cell,
                       front_northeast_cell,
                       back_transition_node,
                       front_northeast_cell,
                       back_northwest_cell,
                       1,
                       2,
                       quarter_side_length,
                       half_side_length + quarter_side_length,
                       half_side_length + quarter_side_length);


    // back Northwest subcell initialization.
    set_cell_node_data(back_northwest_cell,
                       half_side_length,
                       quarter_side_length,
                       3,
                       back_northeast_cell,
                       north_transition_node,
                       west_transition_node,
                       back_southwest_cell,
                       front_northwest_cell,
                       back_transition_node,
                       back_northeast_cell,
                       front_northwest_cell,
                       2,
                       2,
                       quarter_side_length,
                       quarter_side_length,
                       half_side_length + quarter_side_length);

    // front Northwest subcell initialization.
    set_cell_node_data(front_northwest_cell,
                       half_side_length,
                       quarter_side_length,
                       4,
                       front_northeast_cell,
                       north_transition_node,
                       west_transition_node,
                       front_southwest_cell,
                       front_transition_node,
                       back_northwest_cell,
                       back_northwest_cell,
                       front_southwest_cell,
                       3,
                       3,
                       half_side_length + quarter_side_length,
                       quarter_side_length,
                       half_side_length + quarter_side_length);

    // front Southwest subcell initialization.
    set_cell_node_data(front_southwest_cell,
                       half_side_length,
                       quarter_side_length,
                       5,
                       front_southeast_cell,
                       front_northwest_cell,
                       west_transition_node,
                       south_transition_node,
                       front_transition_node,
                       back_southwest_cell,
                       front_northwest_cell,
                       back_southwest_cell,
                       4,
                       3,
                       half_side_length + quarter_side_length,
                       quarter_side_length,
                       quarter_side_length);

    // back Southwest subcell initialization.
    set_cell_node_data(back_southwest_cell,
                       half_side_length,
                       quarter_side_length,
                       6,
                       back_southeast_cell,
                       back_northwest_cell,
                       west_transition_node,
                       south_transition_node,
                       front_southwest_cell,
                       back_transition_node,
                       front_southwest_cell,
                       back_southeast_cell,
                       5,
                       4,
                       quarter_side_length,
                       quarter_side_length,
                       quarter_side_length);


    // back Southeast subcell initialization.
    set_cell_node_data(back_southeast_cell,
                       half_side_length,
                       quarter_side_length,
                       7,
                       east_transition_node,
                       back_northeast_cell,
                       back_southwest_cell,
                       south_transition_node,
                       front_southeast_cell,
                       back_transition_node,
                       back_southwest_cell,
                       front_southeast_cell,
                       6,
                       4,
                       quarter_side_length,
                       half_side_length + quarter_side_length,
                       quarter_side_length);


    // front Southeast subcell initialization.
    set_cell_node_data(front_southeast_cell,
                       half_side_length,
                       quarter_side_length,
                       8,
                       east_transition_node,
                       front_northeast_cell,
                       front_southwest_cell,
                       south_transition_node,
                       front_transition_node,
                       back_southeast_cell,
                       back_southeast_cell,
                       NULL,
                       7,
                       5,
                       half_side_length + quarter_side_length,
                       half_side_length + quarter_side_length,
                       quarter_side_length);


    // Grid initialization
    the_grid->first_cell = front_northeast_cell;
    the_grid->number_of_cells = 8;

}

void print_grid(struct grid* the_grid, FILE * output_file) {

    //TODO:check for activity here???

    struct cell_node *grid_cell = the_grid->first_cell;

    float center_x, center_y, center_z, half_face, v;

    while( grid_cell != 0 ) {

        if(grid_cell->active) {

            center_x = grid_cell->cell_data.center_x;
            center_y = grid_cell->cell_data.center_y;
            center_z = grid_cell->cell_data.center_z;

            v = grid_cell->v;
            half_face = grid_cell->half_face_length;

            /*
            if(count > 0) {
                if (grid_cell->v > -86.0) {
                    act = true;
                }
            }
            else {
                act = true;
            }
             */

            fprintf(output_file, "%lf,%lf,%lf,%lf,%lf\n", center_x, center_y, center_z, half_face, v);

        }
        grid_cell = grid_cell->next;
    }

    //return act;
}

void order_grid_cells(struct grid *the_grid, bool update_active) {

    struct cell_node *grid_cell;
    grid_cell = the_grid->first_cell;


    if(update_active) {
        //TODO: @Incomplete
        //activeCells.clear();
        //activeCells = (CellNode**) malloc(sizeof(CellNode*)*numberOfCells);
    }

    uint64_t counter = 0;
    while( grid_cell != 0 )
    {
        if(grid_cell->active) {

            grid_cell->grid_position = counter;

            if (update_active) {
                //TODO: @Incomplete
                //activeCells[counter] = grid_cell;
                //activeCells.push_back(grid_cell);
            }

            counter++;
        }

        grid_cell = grid_cell->next;
    }

    the_grid->number_of_cells = counter;

}

void free_grid(struct grid *the_grid) {

    //TODO: @Incomplete
    struct cell_node *grid_cell = the_grid->first_cell;

    /*
    // In order to release the memory allocated for the grid, the grid is
    // derefined to level 1. Thus, the grid shape is known and each node can
    // be easily reached.
    while( numberOfCells > 8 ) {
        this->derefineAll();
    }

    // Deleting transition nodes.
    delete static_cast<TransitionNode*>(grid_cell->north);
    delete static_cast<TransitionNode*>(grid_cell->front);
    delete static_cast<TransitionNode*>(grid_cell->east);
    delete static_cast<TransitionNode*>(static_cast<CellNode*>(grid_cell->west)->west);
    delete static_cast<TransitionNode*>(static_cast<CellNode*>(grid_cell->south)->south);
    delete static_cast<TransitionNode*>(static_cast<CellNode*>(grid_cell->back)->back);
*/
    // Deleting cells nodes.
    grid_cell = grid_cell->next;
    while( grid_cell->next != 0 ) {
        free(grid_cell->previous);
        grid_cell = grid_cell->next;
    }
    free(grid_cell);

}

/**
 * Decides if the grid should be refined by traversing the whole grid, according
 * to parameters refinementLevel and refinementBound. A cell will not be refined
 * either if its refinement level  is  equal  to refinementLevel or the  highest
 * of  all  fluxes  coming  into  it  from  the  six  directions  is  less  than
 * refinementBound.
 *
 * @param refinementLevel Minimum refinement level required for the graph.
 * @param refinementBound Minimum flux required for each cell of graph.
 */
//bool refine_grid(double min_h, double refinement_bound) {
//
//        if( min_h <= 0.0 ) {
//            fprintf(stderr,"refine_grid(): Parameter min_h must be positive, passed %lf.", min_h);
//            return false;
//        }
//
//
//        struct cell_node *grid_cell,
//                *auxiliargrid_cell;
//
//        double maximumFlux;
//        bool continueRefining = true;
//        bool refinedOnce = false;
//        setFlux2();
//
//        if(gpu) {
//            refinedThisStep.clear();
//        }
//
//        while( continueRefining ) {
//            continueRefining = false;
//            grid_cell = firstCell;
//            while( grid_cell != 0 ) {
//
//                maximumFlux = grid_cell->getMaximumFlux();
//
//                //cerr << "REF: " << maximumFlux << endl;
//
//                if( ( grid_cell->canChange && grid_cell->active ) &&
//                    ( grid_cell->faceLength > minH ) &&
//                    ( maximumFlux >= refinementBound ) )
//                {
//                    //std::cout << maximumFlux << std::endl;
//                    auxiliargrid_cell = grid_cell;
//                    grid_cell = grid_cell->next;
//                    if(!gpu) {
//                        refineCell( auxiliargrid_cell );
//                    }
//                    else {
//                        refineCell2( auxiliargrid_cell );
//                    }
//                    numberOfCells += 7;
//                    continueRefining = true;
//                    refinedOnce = true;
//                }
//                else {
//                    grid_cell = grid_cell->next;
//                }
//            }
//        }
//        //ordergrid_cells();
//        return refinedOnce;
//}

//void set_grid_flux(struct grid *the_grid) {
//
//    struct cell_node *grid_cell;
//    bool parallel = the_grid->parallel;
//
//    uint64_t active_cells = the_grid->number_of_cells;
//
//    if(!parallel) {
//        grid_cell = the_grid->first_cell;
//        while ( grid_cell != 0 ) {
//            grid_cell->north_flux = 0.0;
//            grid_cell->south_flux = 0.0;
//            grid_cell->east_flux = 0.0;
//            grid_cell->west_flux  = 0.0;
//            grid_cell->front_flux = 0.0;
//            grid_cell->back_flux  = 0.0;
//
//            grid_cell = grid_cell->next;
//        }
//    }
//    else {
//#pragma omp parallel for
//        for (int i = 0; i < active_cells; i++) {
//            //TODO: @Incompelete (not parallel yet)
//            /*
//            activeCells[i]->northFlux = 0.0;
//            activeCells[i]->southFlux = 0.0;
//            activeCells[i]->eastFlux  = 0.0;
//            activeCells[i]->westFlux  = 0.0;
//            activeCells[i]->frontFlux = 0.0;
//            activeCells[i]->backFlux  = 0.0;
//             */
//
//        }
//    }
//
//    if(!parallel) {
//        grid_cell = the_grid->first_cell;
//        while ( grid_cell != 0 )	{
//            if(grid_cell -> active) {
//                grid_cell->setCellFlux2( 's' ); // Computes south flux.
//                grid_cell->setCellFlux2( 'n' ); // Computes north flux.
//                grid_cell->setCellFlux2( 'e' ); // Computes east flux.
//                grid_cell->setCellFlux2( 'w' ); // Computes west flux.
//                grid_cell->setCellFlux2( 'f' ); // Computes front flux.
//                grid_cell->setCellFlux2( 'b' ); // Computes back flux.
//            }
//            grid_cell = grid_cell->next;
//        }
//    }
//    else {
//#pragma omp parallel for
//        for (int i = 0; i < active_cells; i++) {
//            //TODO: @Incompelete (not parallel yet)
//            /*
//            activeCells[i]->setCellFlux2( 's' ); // Computes south flux.
//            activeCells[i]->setCellFlux2( 'n' ); // Computes north flux.
//            activeCells[i]->setCellFlux2( 'e' ); // Computes east flux.
//            activeCells[i]->setCellFlux2( 'w' ); // Computes west flux.
//            activeCells[i]->setCellFlux2( 'f' ); // Computes front flux.
//            activeCells[i]->setCellFlux2( 'b' ); // Computes back flux.
//             */
//
//        }
//    }
//}