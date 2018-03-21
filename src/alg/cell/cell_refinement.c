//
// Created by sachetto on 30/09/17.
//

#include "cell.h"

void refine_cell( struct cell_node *cell, uint32_t *free_sv_positions, uint32_t **refined_this_step)  {

    assert(cell);

    struct transition_node *east_transition_node,
            *north_transition_node,
            *west_transition_node,
            *south_transition_node,
            *front_transition_node,
            *back_transition_node;

    struct cell_node *front_northeast_sub_cell,
            *front_northwest_sub_cell,
            *front_southwest_sub_cell,
            *front_southeast_sub_cell,
            *back_northeast_sub_cell,
            *back_northwest_sub_cell,
            *back_southwest_sub_cell,
            *back_southeast_sub_cell;

    int number_of_hilbert_shape;

    double cell_center_x    = cell->center_x,
            cell_center_y   = cell->center_y,
            cell_center_z   = cell->center_z,
            cell_half_side    = cell->half_face_length,
            cell_quarter_side = cell->half_face_length / 2.0f;

    uint64_t old_bunch_number = cell->bunch_number;

    // Creation of the front northeast cell. This cell, which is to be refined,
    // becomes the frontNortheast cell of the new bunch.
    front_northeast_sub_cell                     = cell;
    front_northeast_sub_cell->cell_data.level    = cell->cell_data.level + (uint16_t )1;
    front_northeast_sub_cell->face_length        = cell_half_side;
    front_northeast_sub_cell->half_face_length   = cell_quarter_side;
    front_northeast_sub_cell->center_x = cell_center_x + cell_quarter_side;
    front_northeast_sub_cell->center_y = cell_center_y + cell_quarter_side;
    front_northeast_sub_cell->center_z = cell_center_z + cell_quarter_side;
    front_northeast_sub_cell->bunch_number       = old_bunch_number * 10 + 1;

    if(refined_this_step && *refined_this_step) {
        sb_push(*refined_this_step, front_northeast_sub_cell->sv_position);
    }

    // Creation of back Northeast node.
    back_northeast_sub_cell = new_cell_node();
    set_refined_cell_data(back_northeast_sub_cell,
                          front_northeast_sub_cell,
                          cell_half_side,
                          cell_quarter_side,
                          cell_center_x - cell_quarter_side,
                          cell_center_y + cell_quarter_side,
                          cell_center_z + cell_quarter_side,
                          old_bunch_number * 10 + 2, free_sv_positions, refined_this_step);


    // Creation of back Northwest node.
    back_northwest_sub_cell = new_cell_node();
    set_refined_cell_data(back_northwest_sub_cell,
                          front_northeast_sub_cell,
                          cell_half_side,
                          cell_quarter_side,
                          cell_center_x - cell_quarter_side,
                          cell_center_y - cell_quarter_side,
                          cell_center_z + cell_quarter_side,
                          old_bunch_number * 10 + 3, free_sv_positions, refined_this_step);

    // Creation of front Northwest node.
    front_northwest_sub_cell = new_cell_node();
    set_refined_cell_data(front_northwest_sub_cell,
                          front_northeast_sub_cell,
                          cell_half_side,
                          cell_quarter_side,
                          cell_center_x + cell_quarter_side,
                          cell_center_y - cell_quarter_side,
                          cell_center_z + cell_quarter_side,
                          old_bunch_number * 10 + 4, free_sv_positions, refined_this_step);


    // Creation of front Southwest node.
    front_southwest_sub_cell = new_cell_node();
    set_refined_cell_data(front_southwest_sub_cell,
                          front_northeast_sub_cell,
                          cell_half_side,
                          cell_quarter_side,
                          cell_center_x + cell_quarter_side,
                          cell_center_y - cell_quarter_side,
                          cell_center_z - cell_quarter_side,
                          old_bunch_number * 10 + 5, free_sv_positions, refined_this_step);


    // Creation of back Southwest node.
    back_southwest_sub_cell = new_cell_node();
    set_refined_cell_data(back_southwest_sub_cell,
                          front_northeast_sub_cell,
                          cell_half_side,
                          cell_quarter_side,
                          cell_center_x - cell_quarter_side,
                          cell_center_y - cell_quarter_side,
                          cell_center_z - cell_quarter_side,
                          old_bunch_number * 10 + 6, free_sv_positions, refined_this_step);



    // Creation of back Southeast node.
    back_southeast_sub_cell = new_cell_node();
    set_refined_cell_data(back_southeast_sub_cell,
                          front_northeast_sub_cell,
                          cell_half_side,
                          cell_quarter_side,
                          cell_center_x - cell_quarter_side,
                          cell_center_y + cell_quarter_side,
                          cell_center_z - cell_quarter_side,
                          old_bunch_number * 10 + 7, free_sv_positions, refined_this_step);


    // Creation of front Southeast node.
    front_southeast_sub_cell = new_cell_node();
    set_refined_cell_data(front_southeast_sub_cell,
                          front_northeast_sub_cell,
                          cell_half_side,
                          cell_quarter_side,
                          cell_center_x + cell_quarter_side,
                          cell_center_y + cell_quarter_side,
                          cell_center_z - cell_quarter_side,
                          old_bunch_number * 10 + 8, free_sv_positions, refined_this_step);

    // west transition node.
    west_transition_node = new_transition_node();
    set_refined_transition_node_data(west_transition_node, front_northeast_sub_cell, 'w');


    // north transition node.
    north_transition_node = new_transition_node();
    set_refined_transition_node_data(north_transition_node, front_northeast_sub_cell, 'n');

    // south transition node.
    south_transition_node = new_transition_node();
    set_refined_transition_node_data(south_transition_node, front_northeast_sub_cell, 's');

    // east transition node.
    east_transition_node = new_transition_node();
    set_refined_transition_node_data(east_transition_node, front_northeast_sub_cell, 'e');


    // front transition node.
    front_transition_node = new_transition_node();
    set_refined_transition_node_data(front_transition_node, front_northeast_sub_cell, 'f');

    // back transition node.
    back_transition_node = new_transition_node();
    set_refined_transition_node_data(back_transition_node, front_northeast_sub_cell, 'b');

    // Linking of new cell nodes and transition nodes.
    front_northeast_sub_cell->north = north_transition_node;
    front_northeast_sub_cell->south = front_southeast_sub_cell;
    front_northeast_sub_cell->east  = east_transition_node;
    front_northeast_sub_cell->west  = front_northwest_sub_cell;
    front_northeast_sub_cell->front = front_transition_node;
    front_northeast_sub_cell->back  = back_northeast_sub_cell;

    back_northeast_sub_cell->north = north_transition_node;
    back_northeast_sub_cell->south = back_southeast_sub_cell;
    back_northeast_sub_cell->east  = east_transition_node;
    back_northeast_sub_cell->west  = back_northwest_sub_cell;
    back_northeast_sub_cell->front = front_northeast_sub_cell;
    back_northeast_sub_cell->back  = back_transition_node;

    back_northwest_sub_cell->north = north_transition_node;
    back_northwest_sub_cell->south = back_southwest_sub_cell;
    back_northwest_sub_cell->east  = back_northeast_sub_cell;
    back_northwest_sub_cell->west  = west_transition_node;
    back_northwest_sub_cell->front = front_northwest_sub_cell;
    back_northwest_sub_cell->back  = back_transition_node;

    front_northwest_sub_cell->north = north_transition_node;
    front_northwest_sub_cell->south = front_southwest_sub_cell;
    front_northwest_sub_cell->east  = front_northeast_sub_cell;
    front_northwest_sub_cell->west  = west_transition_node;
    front_northwest_sub_cell->front = front_transition_node;
    front_northwest_sub_cell->back  = back_northwest_sub_cell;

    front_southwest_sub_cell->north = front_northwest_sub_cell;
    front_southwest_sub_cell->south = south_transition_node;
    front_southwest_sub_cell->east  = front_southeast_sub_cell;
    front_southwest_sub_cell->west  = west_transition_node;
    front_southwest_sub_cell->front = front_transition_node;
    front_southwest_sub_cell->back  = back_southwest_sub_cell;

    back_southwest_sub_cell->north = back_northwest_sub_cell;
    back_southwest_sub_cell->south = south_transition_node;
    back_southwest_sub_cell->east  = back_southeast_sub_cell;
    back_southwest_sub_cell->west  = west_transition_node;
    back_southwest_sub_cell->front = front_southwest_sub_cell;
    back_southwest_sub_cell->back  = back_transition_node;

    back_southeast_sub_cell->north = back_northeast_sub_cell;
    back_southeast_sub_cell->south = south_transition_node;
    back_southeast_sub_cell->east  = east_transition_node;
    back_southeast_sub_cell->west  = back_southwest_sub_cell;
    back_southeast_sub_cell->front = front_southeast_sub_cell;
    back_southeast_sub_cell->back  = back_transition_node;

    front_southeast_sub_cell->north = front_northeast_sub_cell;
    front_southeast_sub_cell->south = south_transition_node;
    front_southeast_sub_cell->east  = east_transition_node;
    front_southeast_sub_cell->west  = front_southwest_sub_cell;
    front_southeast_sub_cell->front = front_transition_node;
    front_southeast_sub_cell->back  = back_southeast_sub_cell;

    /* Connects the cell nodes with the transition nodes.  Quadruple  connectors
    1, 2, 3 and 4 are connected to neighbor cells  in  the  way  depicted below.

    This choice is made consistent with the one made at the function simplifyRef(),
    so that when two transition nodes of  same  level  connected  through  their
    single connectors are eliminated, the  subsequent  linking  of  corresponding
    quadruple connectors is correctly done.

                    front face           back face

                       ______              4______3
                     /|     /|            /|     /|
                   4/_|___3/ |           /_|____/ |
                   |  |___|__|          |  |___|__|
                   | /    | /           | /1   | /2
                   |/_____|/            |/_____|/
                   1      2
            ===========================================
                    west face           east face

                      3______               ______3
                     /|     /|            /|     /|
                   4/_|__ _/ |           /_|___4/ |
                   |  |___|__|          |  |___|__|
                   | /2   | /           | /    | /2
                   |/_____|/            |/_____|/
                   1                           1
            ===========================================
                    north face             south face

                      4______3              ______
                     /|     /|            /|     /|
                   1/_|___2/ |           /_|____/ |
                   |  |___|__|          |  |___|__|
                   | /    | /           | /4   | /3
                   |/_____|/            |/_____|/
                                        1      2
            ===========================================
     */

    // Front face.
    front_transition_node->quadruple_connector1 = front_southwest_sub_cell;
    front_transition_node->quadruple_connector2 = front_southeast_sub_cell;
    front_transition_node->quadruple_connector3 = front_northeast_sub_cell;
    front_transition_node->quadruple_connector4 = front_northwest_sub_cell;

    // Back face.
    back_transition_node->quadruple_connector1 = back_southwest_sub_cell;
    back_transition_node->quadruple_connector2 = back_southeast_sub_cell;
    back_transition_node->quadruple_connector3 = back_northeast_sub_cell;
    back_transition_node->quadruple_connector4 = back_northwest_sub_cell;

    // West face.
    west_transition_node->quadruple_connector1 = front_southwest_sub_cell;
    west_transition_node->quadruple_connector2 = back_southwest_sub_cell;
    west_transition_node->quadruple_connector3 = back_northwest_sub_cell;
    west_transition_node->quadruple_connector4 = front_northwest_sub_cell;

    // East face.
    east_transition_node->quadruple_connector1 = front_southeast_sub_cell;
    east_transition_node->quadruple_connector2 = back_southeast_sub_cell;
    east_transition_node->quadruple_connector3 = back_northeast_sub_cell;
    east_transition_node->quadruple_connector4 = front_northeast_sub_cell;

    // North face.
    north_transition_node->quadruple_connector1 = front_northwest_sub_cell;
    north_transition_node->quadruple_connector2 = front_northeast_sub_cell;
    north_transition_node->quadruple_connector3 = back_northeast_sub_cell;
    north_transition_node->quadruple_connector4 = back_northwest_sub_cell;

    // South face.
    south_transition_node->quadruple_connector1 = front_southwest_sub_cell;
    south_transition_node->quadruple_connector2 = front_southeast_sub_cell;
    south_transition_node->quadruple_connector3 = back_southeast_sub_cell;
    south_transition_node->quadruple_connector4 = back_southwest_sub_cell;



    // Linking bunch neighbor cells to the transition nodes just created.
    struct cell_node *neighbour_cell_node = NULL;
    struct transition_node *neighbour_transition_node = NULL;


    /*==========================================================================
                                WEST TRANSITION NODE

      Points the west neighboring cell to the transition node.

                                         t  t               B: Bunch
                                         | /                n: neighboring cell
                               n -- t -- B -- t             w: Transition node
                                       / |
                                      t  t

     *=========================================================================*/
    char node_type = ((struct basic_cell_data*)west_transition_node->single_connector)->type;
    if( node_type == 'b' ) {
        neighbour_cell_node = (struct cell_node*)(west_transition_node->single_connector);
        neighbour_cell_node->east = west_transition_node;
    }
    else if( node_type == 'w' ) {
        neighbour_transition_node = (struct transition_node*)(west_transition_node->single_connector);

        if( neighbour_transition_node->single_connector == front_northeast_sub_cell )
            neighbour_transition_node->single_connector = west_transition_node;

        else if( neighbour_transition_node->quadruple_connector1 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector1 = west_transition_node;

        else if( neighbour_transition_node->quadruple_connector2 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector2 = west_transition_node;

        else if( neighbour_transition_node->quadruple_connector3 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector3 = west_transition_node;

        else if( neighbour_transition_node->quadruple_connector4 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector4 = west_transition_node;
    }

    /*==========================================================================
                             NORTH TRANSITION NODE

      Points the north neighboring cell to the transition node.

                                      n
                                      |                  B: Bunch
                                      t  t               n: neighboring cell
                                      | /                t: Transition node
                                 t -- B -- t
                                    / |
                                   t  t

    ==========================================================================*/
    node_type = ((struct basic_cell_data*)north_transition_node->single_connector)->type;
    if( node_type == 'b' ) {
        neighbour_cell_node = (struct cell_node*)(north_transition_node->single_connector);
        neighbour_cell_node->south = north_transition_node;
    }
    else if( node_type == 'w' )	{
        neighbour_transition_node = (struct transition_node*)(north_transition_node->single_connector);

        if( neighbour_transition_node->single_connector == front_northeast_sub_cell )
            neighbour_transition_node->single_connector = north_transition_node;

        else if( neighbour_transition_node->quadruple_connector1 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector1 = north_transition_node;

        else if( neighbour_transition_node->quadruple_connector2 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector2 = north_transition_node;

        else if( neighbour_transition_node->quadruple_connector3 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector3 = north_transition_node;

        else if( neighbour_transition_node->quadruple_connector4 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector4 = north_transition_node;
    }

    /*==========================================================================
                                SOUTH TRANSITION NODE

      Points the south neighboring cell to the transition node.

                                      t  t               B: Bunch
                                      | /                n: neighboring cell
                                 t -- B -- t             w: Transition node
                                    / |
                                   t  t
                                      |
                                      n

    ==========================================================================*/
    node_type = ((struct basic_cell_data*)south_transition_node->single_connector)->type;
    if( node_type == 'b' ) {
        neighbour_cell_node = (struct cell_node*)(south_transition_node->single_connector);
        neighbour_cell_node->north = south_transition_node;
    }
    else if( node_type == 'w' )	{
        neighbour_transition_node = (struct transition_node*)(south_transition_node->single_connector);

        if( neighbour_transition_node->single_connector == front_northeast_sub_cell )
            neighbour_transition_node->single_connector = south_transition_node;

        else if( neighbour_transition_node->quadruple_connector1 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector1 = south_transition_node;

        else if( neighbour_transition_node->quadruple_connector2 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector2 = south_transition_node;

        else if( neighbour_transition_node->quadruple_connector3 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector3 = south_transition_node;

        else if( neighbour_transition_node->quadruple_connector4 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector4 = south_transition_node;
    }

    /*==========================================================================
                              EAST TRANSITION NODE

        Points the east neighboring cell to the transition node.

                                     t  t               B: Bunch
                                     | /                n: neighboring cell
                                t -- B -- t -- n        w: Transition node
                                   / |
                                  t  t

    ==========================================================================*/
    node_type = ((struct basic_cell_data*)east_transition_node->single_connector)->type;
    if( node_type == 'b' ) {
        neighbour_cell_node = (struct cell_node*)(east_transition_node->single_connector);
        neighbour_cell_node->west = east_transition_node;
    }
    else if( node_type == 'w' ) {
        neighbour_transition_node = (struct transition_node*)(east_transition_node->single_connector);

        if( neighbour_transition_node->single_connector == front_northeast_sub_cell )
            neighbour_transition_node->single_connector = east_transition_node;

        else if( neighbour_transition_node->quadruple_connector1 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector1 = east_transition_node;

        else if( neighbour_transition_node->quadruple_connector2 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector2 = east_transition_node;

        else if( neighbour_transition_node->quadruple_connector3 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector3 = east_transition_node;

        else if( neighbour_transition_node->quadruple_connector4 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector4 = east_transition_node;
    }

    /*==========================================================================
                                FRONT TRANSITION NODE

      Points the east neighboring cell to the transition node.

                                     t  t               B: Bunch
                                     | /                n: neighboring cell
                                t -- B -- t             t: TRANSITION node
                                   / |
                                  t  t
                                 /
                                n

    ==========================================================================*/
    node_type = ((struct basic_cell_data*)front_transition_node->single_connector)->type;
    if( node_type == 'b' ) {
        neighbour_cell_node = (struct cell_node*)(front_transition_node->single_connector);
        neighbour_cell_node->back = front_transition_node;
    }
    else if( node_type == 'w' ) {
        neighbour_transition_node = (struct transition_node*)(front_transition_node->single_connector);

        if( neighbour_transition_node->single_connector == front_northeast_sub_cell )
            neighbour_transition_node->single_connector = front_transition_node;

        else if( neighbour_transition_node->quadruple_connector1 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector1 = front_transition_node;

        else if( neighbour_transition_node->quadruple_connector2 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector2 = front_transition_node;

        else if( neighbour_transition_node->quadruple_connector3 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector3 = front_transition_node;

        else if( neighbour_transition_node->quadruple_connector4 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector4 = front_transition_node;
    }

    /*==========================================================================
                              BACK TRANSITION NODE

      Points the east neighboring cell to the transition node.

                                          n             B: Bunch
                                         /              n: neighboring cell
                                     t  t               t: Transition Node node
                                     | /
                                t -- B -- t
                                   / |
                                  t  t

    ==========================================================================*/
    node_type = ((struct basic_cell_data*)back_transition_node->single_connector)->type;
    if( node_type == 'b' ) {
        neighbour_cell_node = (struct cell_node*)(back_transition_node->single_connector);
        neighbour_cell_node->front = back_transition_node;
    }
    else if( node_type == 'w' ) {
        neighbour_transition_node = (struct transition_node*)(back_transition_node->single_connector);

        if( neighbour_transition_node->single_connector == front_northeast_sub_cell )
            neighbour_transition_node->single_connector = back_transition_node;

        else if( neighbour_transition_node->quadruple_connector1 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector1 = back_transition_node;

        else if( neighbour_transition_node->quadruple_connector2 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector2 = back_transition_node;

        else if( neighbour_transition_node->quadruple_connector3 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector3 = back_transition_node;

        else if( neighbour_transition_node->quadruple_connector4 == front_northeast_sub_cell )
            neighbour_transition_node->quadruple_connector4 = back_transition_node;
    }


    /*==========================================================================
                ORDERING OF CELL NODES THROUGH HILBERT'S CURVE
    ==========================================================================*/
    number_of_hilbert_shape = front_northeast_sub_cell->hilbert_shape_number;

    if( number_of_hilbert_shape == 0 )	{
        /* Shape 0
                            _______
                           /      /      b: begin
                          /     b/       e: end
                          |  ______
                          | /     /
                          |/    e/
         */

        front_northeast_sub_cell->hilbert_shape_number = 1;
        back_northeast_sub_cell->hilbert_shape_number  = 2;
        back_northwest_sub_cell->hilbert_shape_number  = 2;
        front_northwest_sub_cell->hilbert_shape_number = 3;
        front_southwest_sub_cell->hilbert_shape_number = 3;
        back_southwest_sub_cell->hilbert_shape_number  = 4;
        back_southeast_sub_cell->hilbert_shape_number  = 4;
        front_southeast_sub_cell->hilbert_shape_number = 5;

        front_southeast_sub_cell->next = front_northeast_sub_cell->next;
        front_northeast_sub_cell->next = back_northeast_sub_cell;
        back_northeast_sub_cell->next  = back_northwest_sub_cell;
        back_northwest_sub_cell->next  = front_northwest_sub_cell;
        front_northwest_sub_cell->next = front_southwest_sub_cell;
        front_southwest_sub_cell->next = back_southwest_sub_cell;
        back_southwest_sub_cell->next  = back_southeast_sub_cell;
        back_southeast_sub_cell->next  = front_southeast_sub_cell;

        front_southeast_sub_cell->previous = back_southeast_sub_cell;
        back_southeast_sub_cell->previous  = back_southwest_sub_cell;
        back_southwest_sub_cell->previous  = front_southwest_sub_cell;
        front_southwest_sub_cell->previous = front_northwest_sub_cell;
        front_northwest_sub_cell->previous = back_northwest_sub_cell;
        back_northwest_sub_cell->previous  = back_northeast_sub_cell;
        back_northeast_sub_cell->previous  = front_northeast_sub_cell;

        if( front_southeast_sub_cell->next != 0 )
            front_southeast_sub_cell->next->previous = front_southeast_sub_cell;

    }

    else if( number_of_hilbert_shape == 1 ) {
        /* Shape 1
                                       e
                               /|      |      b: begin
                             /  |   b  |      e: end
                             |  |___|__|
                             |      |
                             |______|
         */

        front_northeast_sub_cell->hilbert_shape_number = 0;
        front_southeast_sub_cell->hilbert_shape_number = 2;
        front_southwest_sub_cell->hilbert_shape_number = 2;
        front_northwest_sub_cell->hilbert_shape_number = 6;
        back_northwest_sub_cell->hilbert_shape_number  = 6;
        back_southwest_sub_cell->hilbert_shape_number  = 7;
        back_southeast_sub_cell->hilbert_shape_number  = 7;
        back_northeast_sub_cell->hilbert_shape_number  = 8;

        back_northeast_sub_cell->next  = front_northeast_sub_cell->next;
        front_northeast_sub_cell->next = front_southeast_sub_cell;
        front_southeast_sub_cell->next = front_southwest_sub_cell;
        front_southwest_sub_cell->next = front_northwest_sub_cell;
        front_northwest_sub_cell->next = back_northwest_sub_cell;
        back_northwest_sub_cell->next  = back_southwest_sub_cell;
        back_southwest_sub_cell->next  = back_southeast_sub_cell;
        back_southeast_sub_cell->next  = back_northeast_sub_cell;

        back_northeast_sub_cell->previous  = back_southeast_sub_cell;
        back_southeast_sub_cell->previous  = back_southwest_sub_cell;
        back_southwest_sub_cell->previous  = back_northwest_sub_cell;
        back_northwest_sub_cell->previous  = front_northwest_sub_cell;
        front_northwest_sub_cell->previous = front_southwest_sub_cell;
        front_southwest_sub_cell->previous = front_southeast_sub_cell;
        front_southeast_sub_cell->previous = front_northeast_sub_cell;

        if( back_northeast_sub_cell->next != 0 )
            back_northeast_sub_cell->next->previous = back_northeast_sub_cell;

    }

    else if( number_of_hilbert_shape == 2 ) {
        /* Shape 2
                               /|     /|      b: begin
                             e/ |   b/ |      e: end
                                |      |
                               /      /
                              /______/
         */

        front_northeast_sub_cell->hilbert_shape_number = 1;
        back_northeast_sub_cell->hilbert_shape_number  = 0;
        back_southeast_sub_cell->hilbert_shape_number  = 0;
        front_southeast_sub_cell->hilbert_shape_number = 9;
        front_southwest_sub_cell->hilbert_shape_number = 9;
        back_southwest_sub_cell->hilbert_shape_number  = 10;
        back_northwest_sub_cell->hilbert_shape_number  = 10;
        front_northwest_sub_cell->hilbert_shape_number = 11;

        front_northwest_sub_cell->next = front_northeast_sub_cell->next;
        front_northeast_sub_cell->next = back_northeast_sub_cell;
        back_northeast_sub_cell->next  = back_southeast_sub_cell;
        back_southeast_sub_cell->next  = front_southeast_sub_cell;
        front_southeast_sub_cell->next = front_southwest_sub_cell;
        front_southwest_sub_cell->next = back_southwest_sub_cell;
        back_southwest_sub_cell->next  = back_northwest_sub_cell;
        back_northwest_sub_cell->next  = front_northwest_sub_cell;

        front_northwest_sub_cell->previous = back_northwest_sub_cell;
        back_northwest_sub_cell->previous  = back_southwest_sub_cell;
        back_southwest_sub_cell->previous  = front_southwest_sub_cell;
        front_southwest_sub_cell->previous = front_southeast_sub_cell;
        front_southeast_sub_cell->previous = back_southeast_sub_cell;
        back_southeast_sub_cell->previous  = back_northeast_sub_cell;
        back_northeast_sub_cell->previous  = front_northeast_sub_cell;

        if( front_northwest_sub_cell->next != 0 )
            front_northwest_sub_cell->next->previous = front_northwest_sub_cell;

    }

    else if( number_of_hilbert_shape == 3 ) {
        /* Shape 3
                               /b     /|      b: begin
                              /______/ |      e: end
                                       |
                               /e     /
                              /______/
         */

        back_northwest_sub_cell->hilbert_shape_number  = 11;
        front_northwest_sub_cell->hilbert_shape_number = 7;
        front_northeast_sub_cell->hilbert_shape_number = 7;
        back_northeast_sub_cell->hilbert_shape_number  = 0;
        back_southeast_sub_cell->hilbert_shape_number  = 0;
        front_southeast_sub_cell->hilbert_shape_number = 9;
        front_southwest_sub_cell->hilbert_shape_number = 9;
        back_southwest_sub_cell->hilbert_shape_number  = 6;

        if( front_northeast_sub_cell->previous != 0 )
            front_northeast_sub_cell->previous->next = back_northwest_sub_cell;

        back_southwest_sub_cell->next  = front_northeast_sub_cell->next;
        back_northwest_sub_cell->next  = front_northwest_sub_cell;
        front_northwest_sub_cell->next = front_northeast_sub_cell;
        front_northeast_sub_cell->next = back_northeast_sub_cell;
        back_northeast_sub_cell->next  = back_southeast_sub_cell;
        back_southeast_sub_cell->next  = front_southeast_sub_cell;
        front_southeast_sub_cell->next = front_southwest_sub_cell;
        front_southwest_sub_cell->next = back_southwest_sub_cell;

        back_northwest_sub_cell->previous  = front_northeast_sub_cell->previous;
        back_southwest_sub_cell->previous  = front_southwest_sub_cell;
        front_southwest_sub_cell->previous = front_southeast_sub_cell;
        front_southeast_sub_cell->previous = back_southeast_sub_cell;
        back_southeast_sub_cell->previous  = back_northeast_sub_cell;
        back_northeast_sub_cell->previous  = front_northeast_sub_cell;
        front_northeast_sub_cell->previous = front_northwest_sub_cell;
        front_northwest_sub_cell->previous = back_northwest_sub_cell;

        if( back_southwest_sub_cell->next != 0 )
            back_southwest_sub_cell->next->previous = back_southwest_sub_cell;

    }

    else if ( number_of_hilbert_shape == 4 ) {
        /* Shape 4
                                /|     /|      b: begin
                               /_|____/ |      e: end
                                 |      |
                                /      /
                              b/     e/
         */

        front_southwest_sub_cell->hilbert_shape_number = 6;
        back_southwest_sub_cell->hilbert_shape_number  = 10;
        back_northwest_sub_cell->hilbert_shape_number  = 10;
        front_northwest_sub_cell->hilbert_shape_number = 7;
        front_northeast_sub_cell->hilbert_shape_number = 7;
        back_northeast_sub_cell->hilbert_shape_number  = 0;
        back_southeast_sub_cell->hilbert_shape_number  = 0;
        front_southeast_sub_cell->hilbert_shape_number = 5;

        if( front_northeast_sub_cell->previous != 0 )
            front_northeast_sub_cell->previous->next = front_southwest_sub_cell;

        front_southeast_sub_cell->next = front_northeast_sub_cell->next;
        front_southwest_sub_cell->next = back_southwest_sub_cell;
        back_southwest_sub_cell->next  = back_northwest_sub_cell;
        back_northwest_sub_cell->next  = front_northwest_sub_cell;
        front_northwest_sub_cell->next = front_northeast_sub_cell;
        front_northeast_sub_cell->next = back_northeast_sub_cell;
        back_northeast_sub_cell->next  = back_southeast_sub_cell;
        back_southeast_sub_cell->next  = front_southeast_sub_cell;

        front_southwest_sub_cell->previous = front_northeast_sub_cell->previous;
        front_southeast_sub_cell->previous = back_southeast_sub_cell;
        back_southeast_sub_cell->previous  = back_northeast_sub_cell;
        back_northeast_sub_cell->previous  = front_northeast_sub_cell;
        front_northeast_sub_cell->previous = front_northwest_sub_cell;
        front_northwest_sub_cell->previous = back_northwest_sub_cell;
        back_northwest_sub_cell->previous  = back_southwest_sub_cell;
        back_southwest_sub_cell->previous  = front_southwest_sub_cell;

        if ( front_southeast_sub_cell->next != 0 )
            front_southeast_sub_cell->next->previous = front_southeast_sub_cell;

    }

    else if( number_of_hilbert_shape == 5 ) {
        /* Shape 5
                                 ______
                                |      |      b: begin
                              __|___   |      e: end
                             |  |   |  |b
                             | /    |
                             |/     |e
         */

        back_southeast_sub_cell->hilbert_shape_number  = 8;
        back_northeast_sub_cell->hilbert_shape_number  = 9;
        back_northwest_sub_cell->hilbert_shape_number  = 9;
        back_southwest_sub_cell->hilbert_shape_number  = 11;
        front_southwest_sub_cell->hilbert_shape_number = 11;
        front_northwest_sub_cell->hilbert_shape_number = 4;
        front_northeast_sub_cell->hilbert_shape_number = 4;
        front_southeast_sub_cell->hilbert_shape_number = 0;

        if( front_northeast_sub_cell->previous != 0 )
            front_northeast_sub_cell->previous->next = back_southeast_sub_cell;

        front_southeast_sub_cell->next = front_northeast_sub_cell->next;
        back_southeast_sub_cell->next  = back_northeast_sub_cell;
        back_northeast_sub_cell->next  = back_northwest_sub_cell;
        back_northwest_sub_cell->next  = back_southwest_sub_cell;
        back_southwest_sub_cell->next  = front_southwest_sub_cell;
        front_southwest_sub_cell->next = front_northwest_sub_cell;
        front_northwest_sub_cell->next = front_northeast_sub_cell;
        front_northeast_sub_cell->next = front_southeast_sub_cell;

        back_southeast_sub_cell->previous  = front_northeast_sub_cell->previous;
        front_southeast_sub_cell->previous = front_northeast_sub_cell;
        front_northeast_sub_cell->previous = front_northwest_sub_cell;
        front_northwest_sub_cell->previous = front_southwest_sub_cell;
        front_southwest_sub_cell->previous = back_southwest_sub_cell;
        back_southwest_sub_cell->previous  = back_northwest_sub_cell;
        back_northwest_sub_cell->previous  = back_northeast_sub_cell;
        back_northeast_sub_cell->previous  = back_southeast_sub_cell;

        if( front_southeast_sub_cell->next != 0 )
            front_southeast_sub_cell->next->previous = front_southeast_sub_cell;

    }

    else if( number_of_hilbert_shape == 6 ) {
        /* Shape 6
                                 ______
                                |      |      b: begin
                              __|___   |      e: end
                             | e|   |  |
                             |      | /
                            b|      |/
         */

        front_southwest_sub_cell->hilbert_shape_number = 10;
        front_northwest_sub_cell->hilbert_shape_number = 4;
        front_northeast_sub_cell->hilbert_shape_number = 4;
        front_southeast_sub_cell->hilbert_shape_number = 1;
        back_southeast_sub_cell->hilbert_shape_number  = 1;
        back_northeast_sub_cell->hilbert_shape_number  = 9;
        back_northwest_sub_cell->hilbert_shape_number  = 9;
        back_southwest_sub_cell->hilbert_shape_number  = 3;

        if( front_northeast_sub_cell->previous != 0 )
            front_northeast_sub_cell->previous->next = front_southwest_sub_cell;

        back_southwest_sub_cell->next  = front_northeast_sub_cell->next;
        front_southwest_sub_cell->next = front_northwest_sub_cell;
        front_northwest_sub_cell->next = front_northeast_sub_cell;
        front_northeast_sub_cell->next = front_southeast_sub_cell;
        front_southeast_sub_cell->next = back_southeast_sub_cell;
        back_southeast_sub_cell->next  = back_northeast_sub_cell;
        back_northeast_sub_cell->next  = back_northwest_sub_cell;
        back_northwest_sub_cell->next  = back_southwest_sub_cell;

        front_southwest_sub_cell->previous = front_northeast_sub_cell->previous;
        back_southwest_sub_cell->previous  = back_northwest_sub_cell;
        back_northwest_sub_cell->previous  = back_northeast_sub_cell;
        back_northeast_sub_cell->previous  = back_southeast_sub_cell;
        back_southeast_sub_cell->previous  = front_southeast_sub_cell;
        front_southeast_sub_cell->previous = front_northeast_sub_cell;
        front_northeast_sub_cell->previous = front_northwest_sub_cell;
        front_northwest_sub_cell->previous = front_southwest_sub_cell;

        if( back_southwest_sub_cell->next != 0 )
            back_southwest_sub_cell->next->previous = back_southwest_sub_cell;

    }

    else if( number_of_hilbert_shape == 7 ) {
        /* Shape 7
                                |b     |e     b: begin
                              __|___   |      e: end
                             |  |   |  |
                             | /    | /
                             |/     |/
         */

        back_northwest_sub_cell->hilbert_shape_number  = 3;
        back_southwest_sub_cell->hilbert_shape_number  = 11;
        front_southwest_sub_cell->hilbert_shape_number = 11;
        front_northwest_sub_cell->hilbert_shape_number = 4;
        front_northeast_sub_cell->hilbert_shape_number = 4;
        front_southeast_sub_cell->hilbert_shape_number = 1;
        back_southeast_sub_cell->hilbert_shape_number  = 1;
        back_northeast_sub_cell->hilbert_shape_number  = 8;

        if( front_northeast_sub_cell->previous != 0 )
            front_northeast_sub_cell->previous->next = back_northwest_sub_cell;

        back_northeast_sub_cell->next  = front_northeast_sub_cell->next;
        back_northwest_sub_cell->next  = back_southwest_sub_cell;
        back_southwest_sub_cell->next  = front_southwest_sub_cell;
        front_southwest_sub_cell->next = front_northwest_sub_cell;
        front_northwest_sub_cell->next = front_northeast_sub_cell;
        front_northeast_sub_cell->next = front_southeast_sub_cell;
        front_southeast_sub_cell->next = back_southeast_sub_cell;
        back_southeast_sub_cell->next  = back_northeast_sub_cell;

        back_northwest_sub_cell->previous  = front_northeast_sub_cell->previous;
        back_northeast_sub_cell->previous  = back_southeast_sub_cell;
        back_southeast_sub_cell->previous  = front_southeast_sub_cell;
        front_southeast_sub_cell->previous = front_northeast_sub_cell;
        front_northeast_sub_cell->previous = front_northwest_sub_cell;
        front_northwest_sub_cell->previous = front_southwest_sub_cell;
        front_southwest_sub_cell->previous = back_southwest_sub_cell;
        back_southwest_sub_cell->previous  = back_northwest_sub_cell;

        if( back_northeast_sub_cell->next != 0 )
            back_northeast_sub_cell->next->previous = back_northeast_sub_cell;

    }

    else if( number_of_hilbert_shape == 8 ) {
        /* Shape 8
                               /|     /e      b: begin
                              /_|____/        e: end
                                |
                               /      /b
                              /______/
         */

        back_southeast_sub_cell->hilbert_shape_number  = 5;
        front_southeast_sub_cell->hilbert_shape_number = 9;
        front_southwest_sub_cell->hilbert_shape_number = 9;
        back_southwest_sub_cell->hilbert_shape_number  = 10;
        back_northwest_sub_cell->hilbert_shape_number  = 10;
        front_northwest_sub_cell->hilbert_shape_number = 7;
        front_northeast_sub_cell->hilbert_shape_number = 7;
        back_northeast_sub_cell->hilbert_shape_number  = 1;

        if( front_northeast_sub_cell->previous != 0 )
            front_northeast_sub_cell->previous->next = back_southeast_sub_cell;

        back_northeast_sub_cell->next  = front_northeast_sub_cell->next;
        back_southeast_sub_cell->next  = front_southeast_sub_cell;
        front_southeast_sub_cell->next = front_southwest_sub_cell;
        front_southwest_sub_cell->next = back_southwest_sub_cell;
        back_southwest_sub_cell->next  = back_northwest_sub_cell;
        back_northwest_sub_cell->next  = front_northwest_sub_cell;
        front_northwest_sub_cell->next = front_northeast_sub_cell;
        front_northeast_sub_cell->next = back_northeast_sub_cell;

        back_southeast_sub_cell->previous  = front_northeast_sub_cell->previous;
        back_northeast_sub_cell->previous  = front_northeast_sub_cell;
        front_northeast_sub_cell->previous = front_northwest_sub_cell;
        front_northwest_sub_cell->previous = back_northwest_sub_cell;
        back_northwest_sub_cell->previous  = back_southwest_sub_cell;
        back_southwest_sub_cell->previous  = front_southwest_sub_cell;
        front_southwest_sub_cell->previous = front_southeast_sub_cell;
        front_southeast_sub_cell->previous = back_southeast_sub_cell;

        if( back_northeast_sub_cell->next != 0 )
            back_northeast_sub_cell->next->previous = back_northeast_sub_cell;

    }

    else if( number_of_hilbert_shape == 9 ) {
        /* Shape 9
                                _______
                               /      /      b: begin
                              /      /       e: end
                              |      |
                              | /e   | /b
                              |/     |/
         */

        back_southeast_sub_cell->hilbert_shape_number  = 5;
        front_southeast_sub_cell->hilbert_shape_number = 8;
        front_northeast_sub_cell->hilbert_shape_number = 8;
        back_northeast_sub_cell->hilbert_shape_number  = 2;
        back_northwest_sub_cell->hilbert_shape_number  = 2;
        front_northwest_sub_cell->hilbert_shape_number = 3;
        front_southwest_sub_cell->hilbert_shape_number = 3;
        back_southwest_sub_cell->hilbert_shape_number  = 6;

        if( front_northeast_sub_cell->previous != 0 )
            front_northeast_sub_cell->previous->next = back_southeast_sub_cell;

        back_southwest_sub_cell->next  = front_northeast_sub_cell->next;
        back_southeast_sub_cell->next  = front_southeast_sub_cell;
        front_southeast_sub_cell->next = front_northeast_sub_cell;
        front_northeast_sub_cell->next = back_northeast_sub_cell;
        back_northeast_sub_cell->next  = back_northwest_sub_cell;
        back_northwest_sub_cell->next  = front_northwest_sub_cell;
        front_northwest_sub_cell->next = front_southwest_sub_cell;
        front_southwest_sub_cell->next = back_southwest_sub_cell;

        back_southeast_sub_cell->previous  = front_northeast_sub_cell->previous;
        back_southwest_sub_cell->previous  = front_southwest_sub_cell;
        front_southwest_sub_cell->previous = front_northwest_sub_cell;
        front_northwest_sub_cell->previous = back_northwest_sub_cell;
        back_northwest_sub_cell->previous  = back_northeast_sub_cell;
        back_northeast_sub_cell->previous  = front_northeast_sub_cell;
        front_northeast_sub_cell->previous = front_southeast_sub_cell;
        front_southeast_sub_cell->previous = back_southeast_sub_cell;

        if( back_southwest_sub_cell->next != 0 )
            back_southwest_sub_cell->next->previous = back_southwest_sub_cell;

    }

    else if( number_of_hilbert_shape == 10 ) {
        /* Shape 10
                                 _______
                                /      /      b: begin
                              e/      /       e: end
                                  ____|__
                                 /    | /
                               b/     |/
         */

        front_southwest_sub_cell->hilbert_shape_number = 6;
        back_southwest_sub_cell->hilbert_shape_number  = 4;
        back_southeast_sub_cell->hilbert_shape_number  = 4;
        front_southeast_sub_cell->hilbert_shape_number = 8;
        front_northeast_sub_cell->hilbert_shape_number = 8;
        back_northeast_sub_cell->hilbert_shape_number  = 2;
        back_northwest_sub_cell->hilbert_shape_number  = 2;
        front_northwest_sub_cell->hilbert_shape_number = 11;

        if( front_northeast_sub_cell->previous != 0 )
            front_northeast_sub_cell->previous->next = front_southwest_sub_cell;

        front_northwest_sub_cell->next = front_northeast_sub_cell->next;
        front_southwest_sub_cell->next = back_southwest_sub_cell;
        back_southwest_sub_cell->next  = back_southeast_sub_cell;
        back_southeast_sub_cell->next  = front_southeast_sub_cell;
        front_southeast_sub_cell->next = front_northeast_sub_cell;
        front_northeast_sub_cell->next = back_northeast_sub_cell;
        back_northeast_sub_cell->next  = back_northwest_sub_cell;
        back_northwest_sub_cell->next  = front_northwest_sub_cell;

        front_southwest_sub_cell->previous = front_northeast_sub_cell->previous;
        front_northwest_sub_cell->previous = back_northwest_sub_cell;
        back_northwest_sub_cell->previous  = back_northeast_sub_cell;
        back_northeast_sub_cell->previous  = front_northeast_sub_cell;
        front_northeast_sub_cell->previous = front_southeast_sub_cell;
        front_southeast_sub_cell->previous = back_southeast_sub_cell;
        back_southeast_sub_cell->previous  = back_southwest_sub_cell;
        back_southwest_sub_cell->previous  = front_southwest_sub_cell;

        if( front_northwest_sub_cell->next != 0 )
            front_northwest_sub_cell->next->previous = front_northwest_sub_cell;

    }

    else if( number_of_hilbert_shape == 11 ) {
        /* Shape 11
                                     b______
                                            |      b: begin
                                   e______  |      e: end
                                     _____|_|
                                    /     |
                                   /______|
         */

        back_northwest_sub_cell->hilbert_shape_number  = 7;
        back_northeast_sub_cell->hilbert_shape_number  = 3;
        back_southeast_sub_cell->hilbert_shape_number  = 3;
        back_southwest_sub_cell->hilbert_shape_number  = 5;
        front_southwest_sub_cell->hilbert_shape_number = 5;
        front_southeast_sub_cell->hilbert_shape_number = 10;
        front_northeast_sub_cell->hilbert_shape_number = 10;
        front_northwest_sub_cell->hilbert_shape_number = 2;

        if( front_northeast_sub_cell->previous != 0 )
            front_northeast_sub_cell->previous->next = back_northwest_sub_cell;

        front_northwest_sub_cell->next = front_northeast_sub_cell->next;
        back_northwest_sub_cell->next  = back_northeast_sub_cell;
        back_northeast_sub_cell->next  = back_southeast_sub_cell;
        back_southeast_sub_cell->next  = back_southwest_sub_cell;
        back_southwest_sub_cell->next  = front_southwest_sub_cell;
        front_southwest_sub_cell->next = front_southeast_sub_cell;
        front_southeast_sub_cell->next = front_northeast_sub_cell;
        front_northeast_sub_cell->next = front_northwest_sub_cell;

        back_northwest_sub_cell->previous  = front_northeast_sub_cell->previous;
        front_northwest_sub_cell->previous = front_northeast_sub_cell;
        front_northeast_sub_cell->previous = front_southeast_sub_cell;
        front_southeast_sub_cell->previous = front_southwest_sub_cell;
        front_southwest_sub_cell->previous = back_southwest_sub_cell;
        back_southwest_sub_cell->previous  = back_southeast_sub_cell;
        back_southeast_sub_cell->previous  = back_northeast_sub_cell;
        back_northeast_sub_cell->previous  = back_northwest_sub_cell;

        if( front_northwest_sub_cell->next != 0 )
            front_northwest_sub_cell->next->previous = front_northwest_sub_cell;

    }

    // If necessary, simplifies the graph by eliminating adjacent transition nodes
    // of same level connected through their single connectors.
    simplify_refinement( east_transition_node  );
    simplify_refinement( north_transition_node );
    simplify_refinement( west_transition_node  );
    simplify_refinement( south_transition_node );
    simplify_refinement( front_transition_node );
    simplify_refinement( back_transition_node  );
}


/**
 * Simplifies data structure eliminating adjacent transition nodes of same level.
 *
 * @param transition_node Candidate transition node to be eliminated.
 *
 */
void simplify_refinement( struct transition_node *transition_node ) {

    if( transition_node == NULL ) {
        fprintf(stderr, "simplify_refinement: Parameter transition_node is NULL. Exiting!");
        exit(10);
    }

    // Pointers used to convert the Cell in a cell node or transition node.
    struct transition_node *neighbour_transition_node;
    struct cell_node *neighbour_cell_node;

    if( transition_node->single_connector != 0 ) {

        // Both transition node and neighbor transition node must have the same
        // refinement level.
        char node_type = transition_node->cell_data.type;
        uint16_t node_level = transition_node->cell_data.level;

        uint16_t single_connector_level = ((struct basic_cell_data*)(transition_node->single_connector))->level;

        if( ( node_type == 'w') && (node_level == single_connector_level) ) {
            struct transition_node *neighbour_node = (struct transition_node*) (transition_node->single_connector);

            struct cell_node *cellNode[4];
            cellNode[0] = (struct cell_node*)(transition_node->quadruple_connector1);
            cellNode[1] = (struct cell_node*)(transition_node->quadruple_connector2);
            cellNode[2] = (struct cell_node*)(transition_node->quadruple_connector3);
            cellNode[3] = (struct cell_node*)(transition_node->quadruple_connector4);

            struct cell_node *neighborCell[4];
            neighborCell[0] = neighbour_node->quadruple_connector1;
            neighborCell[1] = neighbour_node->quadruple_connector2;
            neighborCell[2] = neighbour_node->quadruple_connector3;
            neighborCell[3] = neighbour_node->quadruple_connector4;

            char direction = transition_node->direction;
            char type;

            for( int i = 0; i < 4; i++ ) {
                switch( direction ) {
                    case 'n': { cellNode[i]->north = neighborCell[i]; break; }
                    case 's': { cellNode[i]->south = neighborCell[i]; break; }
                    case 'e': { cellNode[i]->east  = neighborCell[i]; break; }
                    case 'w': { cellNode[i]->west  = neighborCell[i]; break; }
                    case 'f': { cellNode[i]->front = neighborCell[i]; break; }
                    case 'b': { cellNode[i]->back  = neighborCell[i]; break; }
                    default: break;
                }

                type = neighborCell[i]->cell_data.type;
                switch( type ) {
                    case 'b': {
                        neighbour_cell_node = neighborCell[i];
                        switch( direction )	{
                            case 'n': { neighbour_cell_node->south = cellNode[i]; break; }
                            case 's': { neighbour_cell_node->north = cellNode[i]; break; }
                            case 'e': { neighbour_cell_node->west  = cellNode[i]; break; }
                            case 'w': { neighbour_cell_node->east  = cellNode[i]; break; }
                            case 'f': { neighbour_cell_node->back  = cellNode[i]; break; }
                            case 'b': { neighbour_cell_node->front = cellNode[i]; break; }
                            default: break;
                        }
                        break;
                    }

                    case 'w': {
                        neighbour_transition_node = (struct transition_node*)(neighborCell[i]);
                        if( neighbour_node == neighbour_transition_node->single_connector )
                            neighbour_transition_node->single_connector = cellNode[i];

                        else if( neighbour_node == neighbour_transition_node->quadruple_connector1 )
                            neighbour_transition_node->quadruple_connector1 = cellNode[i];

                        else if( neighbour_node == neighbour_transition_node->quadruple_connector2 )
                            neighbour_transition_node->quadruple_connector2 = cellNode[i];

                        else if( neighbour_node == neighbour_transition_node->quadruple_connector3 )
                            neighbour_transition_node->quadruple_connector3 = cellNode[i];

                        else if( neighbour_node == neighbour_transition_node->quadruple_connector4 )
                            neighbour_transition_node->quadruple_connector4 = cellNode[i];

                        break;
                    }

                    default: break;
                }
            }
            free(transition_node);
            free(neighbour_node);
        }
    }
}

void set_refined_cell_data(struct cell_node* the_cell, struct cell_node* other_cell,
                           double face_length, double half_face_length,
                           double center_x, double center_y, double center_z,
                           uint64_t  bunch_number, uint32_t *free_sv_positions,
                           uint32_t **refined_this_step) {


    the_cell->cell_data.level = other_cell->cell_data.level;
    the_cell->active = other_cell->active;
    the_cell->fibrotic = other_cell->fibrotic;
    the_cell->border_zone = other_cell->border_zone;
    the_cell->scar_type = other_cell->scar_type;
    the_cell->v = other_cell->v;

    the_cell->face_length = face_length;
    the_cell->half_face_length = half_face_length;
    the_cell->center_x = center_x;
    the_cell->center_y = center_y;
    the_cell->center_z = center_z;
    the_cell->bunch_number = bunch_number;

    if(free_sv_positions)
        the_cell->sv_position = sb_pop(free_sv_positions);

    if(refined_this_step && *refined_this_step)
        sb_push(*refined_this_step, the_cell->sv_position);

}

void set_refined_transition_node_data(struct transition_node *the_node, struct cell_node* other_node, char direction) {


    the_node->direction          = direction;
    the_node->cell_data.level    = other_node->cell_data.level;

    switch(direction) {
        case 'w': the_node->single_connector   = other_node->west; break;
        case 'n': the_node->single_connector   = other_node->north; break;
        case 's': the_node->single_connector   = other_node->south; break;
        case 'e': the_node->single_connector   = other_node->east; break;
        case 'f': the_node->single_connector   = other_node->front; break;
        case 'b': the_node->single_connector   = other_node->back; break;
        default:
            fprintf(stderr, "set_refined_transition_node_data() invalid direction %c Exiting!", direction);
            exit(10);
    }

}
