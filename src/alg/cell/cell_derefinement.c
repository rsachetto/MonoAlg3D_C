//
// Created by sachetto on 30/09/17.
//

#include "cell.h"

/**
* Decides if the bunch should be derefined. A bunch will not be derefined if
* the highest of all directional fluxes coming into its  six  directions  is
* greater than derefinement_bound.
*
* @param grid_cell Cell to be verified if it needs to be derefined.
* @param dererefinementBound Maximum flux allowed to this cell.
* @throw NullPointer If a null cell node is given as argument, a NullPointer
* exception is thrown.
*/
bool cell_needs_derefinement (struct cell_node *grid_cell, double derefinement_bound) {

    if (grid_cell == NULL) {
        fprintf (stderr, "cell_needs_derefinement: Parameter grid_cell is NULL.");
    }

    struct cell_node *first_cell = grid_cell;
    struct cell_node *second_cell = first_cell->next;
    struct cell_node *third_cell = second_cell->next;
    struct cell_node *fourth_cell = third_cell->next;
    struct cell_node *fifth_cell = fourth_cell->next;
    struct cell_node *sixth_cell = fifth_cell->next;
    struct cell_node *seventh_cell = sixth_cell->next;
    struct cell_node *eighth_cell = seventh_cell->next;

    double maximum1 = get_cell_maximum_flux (first_cell);
    double maximum2 = get_cell_maximum_flux (second_cell);
    double maximum3 = get_cell_maximum_flux (third_cell);
    double maximum4 = get_cell_maximum_flux (fourth_cell);
    double maximum5 = get_cell_maximum_flux (fifth_cell);
    double maximum6 = get_cell_maximum_flux (sixth_cell);
    double maximum7 = get_cell_maximum_flux (seventh_cell);
    double maximum8 = get_cell_maximum_flux (eighth_cell);

    double highest_maximum = maximum1;
    if (maximum2 > highest_maximum)
        highest_maximum = maximum2;

    if (maximum3 > highest_maximum)
        highest_maximum = maximum3;

    if (maximum4 > highest_maximum)
        highest_maximum = maximum4;

    if (maximum5 > highest_maximum)
        highest_maximum = maximum5;

    if (maximum6 > highest_maximum)
        highest_maximum = maximum6;

    if (maximum7 > highest_maximum)
        highest_maximum = maximum7;

    if (maximum8 > highest_maximum)
        highest_maximum = maximum8;

    bool derefinement_condition = false;

    if (highest_maximum <= derefinement_bound)
        derefinement_condition = true;

    return derefinement_condition;
}

void derefine_cell_bunch (struct cell_node *first_bunch_cell, uint32_t **free_sv_positions) {
    if (first_bunch_cell == 0) {
        fprintf (stderr, "derefine_cell_bunch: Parameter first_bunch_cell is NULL. Exiting!!");
        exit (10);
    }

    struct cell_node *cell_before_bunch = first_bunch_cell->previous;
    struct cell_node *cell_after_bunch =
        first_bunch_cell->next->next->next->next->next->next->next->next;

    uint16_t bunch_level = first_bunch_cell->cell_data.level;
    uint8_t hilbert_shape_number = get_father_bunch_number (first_bunch_cell);
    uint64_t bunch_number = first_bunch_cell->bunch_number;

    // New cell variable (Arithmetic mean between all cells of the bunch).
    double u = 0;

    struct cell_node *auxiliar = first_bunch_cell;
    for (int i = 0; i < 8; i++) {
        u += auxiliar->v;
        auxiliar = auxiliar->next;
    }
    u /= 8.0;

    // Front northeast node of the bunch becomes the derefined node.
    struct cell_node *new_cell = get_front_northeast_cell (first_bunch_cell);

    if(free_sv_positions && *free_sv_positions) {

        //Free Sv Map positions//////////////////////////////////
        struct cell_node *w = (struct cell_node *) new_cell->west;
        struct cell_node *b = (struct cell_node *) new_cell->back;
        struct cell_node *s = (struct cell_node *) new_cell->south;

        struct cell_node *wb = (struct cell_node *) w->back;
        struct cell_node *ws = (struct cell_node *) w->south;
        struct cell_node *sb = (struct cell_node *) s->back;

        struct cell_node *wsb = ws->back;

        sb_push(*free_sv_positions, w->sv_position);
        sb_push(*free_sv_positions, b->sv_position);
        sb_push(*free_sv_positions, s->sv_position);

        sb_push(*free_sv_positions, wb->sv_position);
        sb_push(*free_sv_positions, ws->sv_position);
        sb_push(*free_sv_positions, sb->sv_position);
        sb_push(*free_sv_positions, wsb->sv_position);
        /////////////////////////////////////////////////////////////

    }
    new_cell->previous = cell_before_bunch;
    new_cell->next = cell_after_bunch;

    if (new_cell->previous != 0)
        new_cell->previous->next = new_cell;
    if (new_cell->next != 0)
        new_cell->next->previous = new_cell;

    double aux_center_x = ((struct cell_node *)(new_cell->back))->center_x;
    double aux_center_y = ((struct cell_node *)(new_cell->west))->center_y;
    double aux_center_z = ((struct cell_node *)(new_cell->south))->center_z;

    // New geometric variables.
    new_cell->center_x = (new_cell->center_x + aux_center_x) / 2.0f;
    new_cell->center_y = (new_cell->center_y + aux_center_y) / 2.0f;
    new_cell->center_z = (new_cell->center_z + aux_center_z) / 2.0f;

    new_cell->face_length = 2.0f * new_cell->face_length;
    new_cell->half_face_length = 2.0f * new_cell->half_face_length;
    new_cell->v = u;

    new_cell->cell_data.level = bunch_level - (uint8_t)1;
    new_cell->hilbert_shape_number = hilbert_shape_number;
    new_cell->bunch_number = bunch_number / (int)10;

    struct cell_node *front_northeast_cell = new_cell;
    struct cell_node *front_southeast_cell = (struct cell_node *)(front_northeast_cell->south);
    struct cell_node *front_northwest_cell = (struct cell_node *)(front_northeast_cell->west);
    struct cell_node *front_southwest_cell = (struct cell_node *)(front_northwest_cell->south);
    struct cell_node *back_northeast_cell = (struct cell_node *)(front_northeast_cell->back);
    struct cell_node *back_southeast_cell = (struct cell_node *)(back_northeast_cell->south);
    struct cell_node *back_northwest_cell = (struct cell_node *)(back_northeast_cell->west);
    struct cell_node *back_southwest_cell = (struct cell_node *)(back_northwest_cell->south);

    // Creation of North Transition Node.
    struct transition_node *north_transition_node = new_transition_node();
    set_transition_node_data (north_transition_node, bunch_level, 'n', new_cell,
                              front_northwest_cell->north, front_northeast_cell->north,
                              back_northeast_cell->north, back_northwest_cell->north);

    // Creation of South Transition Node.
    struct transition_node *south_transition_node = new_transition_node();
    set_transition_node_data (south_transition_node, bunch_level, 's', new_cell,
                              front_southwest_cell->south, front_southeast_cell->south,
                              back_southeast_cell->south, back_southwest_cell->south);

    // Creation of East Transition Node.
    struct transition_node *east_transition_node = new_transition_node();
    set_transition_node_data (east_transition_node, bunch_level, 'e', new_cell,
                              front_southeast_cell->east, back_southeast_cell->east,
                              back_northeast_cell->east, front_northeast_cell->east);

    // Creation of West Transition Node.
    struct transition_node *west_transition_node = new_transition_node();
    set_transition_node_data (west_transition_node, bunch_level, 'w', new_cell,
                              front_southwest_cell->west, back_southwest_cell->west,
                              back_northwest_cell->west, front_northwest_cell->west);

    // Creation of Front Transition Node.
    struct transition_node *front_transition_node = new_transition_node();
    set_transition_node_data (front_transition_node, bunch_level, 'f', new_cell,
                              front_southwest_cell->front, front_southeast_cell->front,
                              front_northeast_cell->front, front_northwest_cell->front);

    // Creation of Back Transition Node.
    struct transition_node *back_transition_node = new_transition_node();
    set_transition_node_data (back_transition_node, bunch_level, 'b', new_cell,
                              back_southwest_cell->back, back_southeast_cell->back,
                              back_northeast_cell->back, back_northwest_cell->back);

    // Elimination of the seven unneeded bunch cells.
    free_cell_node (front_northwest_cell);
    free_cell_node (front_southwest_cell);
    free_cell_node (front_southeast_cell);
    free_cell_node (back_northeast_cell);
    free_cell_node (back_northwest_cell);
    free_cell_node (back_southeast_cell);
    free_cell_node (back_southwest_cell);

    // Linking of derefined cell and new transition nodes.
    new_cell->north = north_transition_node;
    new_cell->south = south_transition_node;
    new_cell->east = east_transition_node;
    new_cell->west = west_transition_node;
    new_cell->front = front_transition_node;
    new_cell->back = back_transition_node;

    // Simplification of grid Eliminating unneeded transition nodes.
    simplify_derefinement(north_transition_node);
    simplify_derefinement(south_transition_node);
    simplify_derefinement(east_transition_node);
    simplify_derefinement(west_transition_node);
    simplify_derefinement(front_transition_node);
    simplify_derefinement(back_transition_node);
}

/**
 * Gets a pointer to the frontNortheast cell of a bunch. In order to get this
 * cell, the first cell of the bunch has to be given as argument. If this cell is
 * not the first cell, the behavior of this function and its result is undefined.
 *
 * @param first_bunch_cell First cell of the bunch where the front northeast cell
 * will be returned.
 * @return The front northeast cell of a bunch.
 * @throw NullPointer If a null cell node is given as argument, a NullPointer
 * exception is thrown.
 */
struct cell_node *get_front_northeast_cell (struct cell_node *first_bunch_cell) {
    if (first_bunch_cell == NULL) {
        fprintf (stderr, "get_front_northeast_cell: Parameter first_bunch_cell is NULL.");
    }

    struct cell_node *first_cell = first_bunch_cell;
    struct cell_node *second_cell = first_cell->next;
    struct cell_node *third_cell = second_cell->next;
    struct cell_node *fourth_cell = third_cell->next;
    struct cell_node *fifth_cell = fourth_cell->next;
    struct cell_node *sixth_cell = fifth_cell->next;
    struct cell_node *seventh_cell = sixth_cell->next;
    struct cell_node *eighth_cell = seventh_cell->next;

    double coordinateSum1 = first_cell->center_x + first_cell->center_y + first_cell->center_z;
    double coordinateSum2 = second_cell->center_x + second_cell->center_y + second_cell->center_z;
    double coordinateSum3 = third_cell->center_x + third_cell->center_y + third_cell->center_z;
    double coordinateSum4 = fourth_cell->center_x + fourth_cell->center_y + fourth_cell->center_z;
    double coordinateSum5 = fifth_cell->center_x + fifth_cell->center_y + fifth_cell->center_z;
    double coordinateSum6 = sixth_cell->center_x + sixth_cell->center_y + sixth_cell->center_z;
    double coordinateSum7 =
        seventh_cell->center_x + seventh_cell->center_y + seventh_cell->center_z;
    double coordinateSum8 = eighth_cell->center_x + eighth_cell->center_y + eighth_cell->center_z;

    double maximum;
    struct cell_node *front_northeast_cell = first_cell;
    maximum = coordinateSum1;
    if (coordinateSum2 > maximum) {
        maximum = coordinateSum2;
        front_northeast_cell = second_cell;
    }
    if (coordinateSum3 > maximum) {
        maximum = coordinateSum3;
        front_northeast_cell = third_cell;
    }
    if (coordinateSum4 > maximum) {
        maximum = coordinateSum4;
        front_northeast_cell = fourth_cell;
    }
    if (coordinateSum5 > maximum) {
        maximum = coordinateSum5;
        front_northeast_cell = fifth_cell;
    }
    if (coordinateSum6 > maximum) {
        maximum = coordinateSum6;
        front_northeast_cell = sixth_cell;
    }
    if (coordinateSum7 > maximum) {
        maximum = coordinateSum7;
        front_northeast_cell = seventh_cell;
    }
    if (coordinateSum8 > maximum) {
        front_northeast_cell = eighth_cell;
    }
    return front_northeast_cell;
}

/**
 * Gets the father's bunch number of a bunch, that is, the basic hilbert shape
 * number of a bunch. It is useful mostly on derefinement procedure, where a group
 * of eight cells become a single cell with the same information of the group's
 * father. In order to get the basic hilbert shape number of a bunch, its first
 * cell  has to be given as argument. If the argument is not the first cell, the
 * value -1 is returned.
 *
 * @param first_bunch_cell First cell of the bunch where the basic hilbert shape
 * number will be returned.
 * @return The hilbert shape number of a bunch.
 * @throw NullPointer If a null cell node is given as argument, a NullPointer
 * exception is thrown.
 */
uint8_t get_father_bunch_number (struct cell_node *first_bunch_cell) {

    if (first_bunch_cell == NULL) {
        fprintf (stderr, "Grid::getFatherBunchNumber(): Parameter first_bunch_cell is NULL.");
        exit (10);
    }

    struct cell_node *firstCell = first_bunch_cell;
    struct cell_node *secondCell = firstCell->next;
    struct cell_node *thirdCell = secondCell->next;
    struct cell_node *fourthCell = thirdCell->next;
    struct cell_node *fifthCell = fourthCell->next;
    struct cell_node *sixthCell = fifthCell->next;
    struct cell_node *seventhCell = sixthCell->next;
    struct cell_node *eighthCell = seventhCell->next;

    uint8_t hilbert_shape_number = 0;
    uint8_t sum_of_hilbert_shape_numbers =
        firstCell->hilbert_shape_number + secondCell->hilbert_shape_number +
        thirdCell->hilbert_shape_number + fourthCell->hilbert_shape_number +
        fifthCell->hilbert_shape_number + sixthCell->hilbert_shape_number +
        seventhCell->hilbert_shape_number + eighthCell->hilbert_shape_number;

    switch (sum_of_hilbert_shape_numbers) {
    case 24: {
        hilbert_shape_number = 0;
        break;
    }
    case 38: {
        hilbert_shape_number = 1;
        break;
    }
    case 50: {
        hilbert_shape_number = 2;
        break;
    }
    case 49: {
        hilbert_shape_number = 3;
        break;
    }
    case 45: {
        if (secondCell->hilbert_shape_number == 10)
            hilbert_shape_number = 4;

        else if (secondCell->hilbert_shape_number == 4)
            hilbert_shape_number = 10;

        else if (secondCell->hilbert_shape_number == 3)
            hilbert_shape_number = 11;

        break;
    }
    case 56: {
        hilbert_shape_number = 5;
        break;
    }
    case 41: {
        hilbert_shape_number = 6;
        break;
    }
    case 43: {
        hilbert_shape_number = 7;
        break;
    }
    case 58: {
        hilbert_shape_number = 8;
        break;
    }
    case 37: {
        hilbert_shape_number = 9;
        break;
    }
    default: {
        hilbert_shape_number = 99;
        break;
    }
    }
    return hilbert_shape_number;
}

/**
 * Simplifies the graph by eliminating adjacent transition nodes of same  level
 * connected through their quadruple connectors. If a transition  node  is  not
 * eliminated, it's  because there is no adjacent transition node of same level,
 * then simply connects the outside to it.
 *
 * @param transition_node Candidate transition node to be eliminated.
 * @throw NullPointer If a null transition node is given as argument, a NullPointer
 * exception is thrown.
 */
void simplify_derefinement(struct transition_node *transition_node) {
    if (transition_node == NULL) {
        fprintf (stderr, "simplify_derefinement(): Parameter transition_node is NULL. Exiting!");
        exit (10);
    }

    struct cell_node *derefinedCell = (struct cell_node *)(transition_node->single_connector);
    char direction = transition_node->direction;

    struct transition_node *white_neighbor_cell;
    struct cell_node *blackNeighborCell;
    struct transition_node *neighbor_transition_node;
    void *quadrupleConnector[4];

    /* All quadruple connectors point to the  same  cell.  It  means  that  two
     * transition nodes of the same level have been connected. This connection is
     * not allowed and it is simplified through a direct connection between the two
     * cell nodes pointed by each transition node. Then, both transition nodes
     * are eliminated.
     */
    if ((transition_node->quadruple_connector1 == transition_node->quadruple_connector2) &&
        (transition_node->quadruple_connector1 == transition_node->quadruple_connector3) &&
        (transition_node->quadruple_connector1 == transition_node->quadruple_connector4)) {
        // PS: No matter which quadruple connector is chosen in the line below,
        // because each one of them points to the same cell.
        neighbor_transition_node =
            (struct transition_node *)(transition_node->quadruple_connector1);

        char neighborCellType =
            ((struct basic_cell_data *)(neighbor_transition_node->single_connector))->type;

        switch (direction) {
        case 'n': {
            derefinedCell->north = neighbor_transition_node->single_connector;
            break;
        }
        case 's': {
            derefinedCell->south = neighbor_transition_node->single_connector;
            break;
        }
        case 'e': {
            derefinedCell->east = neighbor_transition_node->single_connector;
            break;
        }
        case 'w': {
            derefinedCell->west = neighbor_transition_node->single_connector;
            break;
        }
        case 'f': {
            derefinedCell->front = neighbor_transition_node->single_connector;
            break;
        }
        case 'b': {
            derefinedCell->back = neighbor_transition_node->single_connector;
            break;
        }
        default: { break; }
        }

        if (neighborCellType == 'b') {
            blackNeighborCell = (struct cell_node *)(neighbor_transition_node->single_connector);
            switch (direction) {
            case 'n': {
                blackNeighborCell->south = derefinedCell;
                break;
            }
            case 's': {
                blackNeighborCell->north = derefinedCell;
                break;
            }
            case 'e': {
                blackNeighborCell->west = derefinedCell;
                break;
            }
            case 'w': {
                blackNeighborCell->east = derefinedCell;
                break;
            }
            case 'f': {
                blackNeighborCell->back = derefinedCell;
                break;
            }
            case 'b': {
                blackNeighborCell->front = derefinedCell;
                break;
            }
            default: { break; }
            }
        } else {
            white_neighbor_cell =
                (struct transition_node *)(neighbor_transition_node->single_connector);

            if (white_neighbor_cell->single_connector == neighbor_transition_node)
                white_neighbor_cell->single_connector = derefinedCell;

            else if (white_neighbor_cell->quadruple_connector1 == neighbor_transition_node)
                white_neighbor_cell->quadruple_connector1 = derefinedCell;

            else if (white_neighbor_cell->quadruple_connector2 == neighbor_transition_node)
                white_neighbor_cell->quadruple_connector2 = derefinedCell;

            else if (white_neighbor_cell->quadruple_connector3 == neighbor_transition_node)
                white_neighbor_cell->quadruple_connector3 = derefinedCell;

            else if (white_neighbor_cell->quadruple_connector4 == neighbor_transition_node)
                white_neighbor_cell->quadruple_connector4 = derefinedCell;
        }
        free (neighbor_transition_node);
        free (transition_node);
    }

    /* Connects outside to the transition node, if this was not deleted. That is,
     * the quadruple connectors point to different cells. */
    else {
        quadrupleConnector[0] = transition_node->quadruple_connector1;
        quadrupleConnector[1] = transition_node->quadruple_connector2;
        quadrupleConnector[2] = transition_node->quadruple_connector3;
        quadrupleConnector[3] = transition_node->quadruple_connector4;

        for (int i = 0; i < 4; i++) {
            char connector_type = ((struct basic_cell_data *)(quadrupleConnector[i]))->type;
            if (connector_type == TRANSITION_NODE_TYPE) {
                white_neighbor_cell = (struct transition_node *)(quadrupleConnector[i]);
                white_neighbor_cell->single_connector = transition_node;
            } else if (connector_type == CELL_NODE_TYPE) {
                blackNeighborCell = (struct cell_node *)(quadrupleConnector[i]);
                switch (direction) {
                case 'n': {
                    blackNeighborCell->south = transition_node;
                    break;
                }
                case 's': {
                    blackNeighborCell->north = transition_node;
                    break;
                }
                case 'e': {
                    blackNeighborCell->west = transition_node;
                    break;
                }
                case 'w': {
                    blackNeighborCell->east = transition_node;
                    break;
                }
                case 'f': {
                    blackNeighborCell->back = transition_node;
                    break;
                }
                case 'b': {
                    blackNeighborCell->front = transition_node;
                    break;
                }
                default: { break; }
                }
            }
        }
    }
}
