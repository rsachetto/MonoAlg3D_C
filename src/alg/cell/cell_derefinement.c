//
// Created by sachetto on 30/09/17.
//

#include <assert.h>

#include "../../3dparty/stb_ds.h"
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
bool cell_needs_derefinement(struct cell_node *grid_cell, real_cpu derefinement_bound) {

    if(grid_cell == NULL) {
        fprintf(stderr, "cell_needs_derefinement: Parameter grid_cell is NULL.");
        return false;
    }

    struct cell_node *first_cell = grid_cell;
    struct cell_node *second_cell = first_cell->next;
    struct cell_node *third_cell = second_cell->next;
    struct cell_node *fourth_cell = third_cell->next;
    struct cell_node *fifth_cell = fourth_cell->next;
    struct cell_node *sixth_cell = fifth_cell->next;
    struct cell_node *seventh_cell = sixth_cell->next;
    struct cell_node *eighth_cell = seventh_cell->next;

    real_cpu maximum1 = get_cell_maximum_flux(first_cell);
    real_cpu maximum2 = get_cell_maximum_flux(second_cell);
    real_cpu maximum3 = get_cell_maximum_flux(third_cell);
    real_cpu maximum4 = get_cell_maximum_flux(fourth_cell);
    real_cpu maximum5 = get_cell_maximum_flux(fifth_cell);
    real_cpu maximum6 = get_cell_maximum_flux(sixth_cell);
    real_cpu maximum7 = get_cell_maximum_flux(seventh_cell);
    real_cpu maximum8 = get_cell_maximum_flux(eighth_cell);

    real_cpu highest_maximum = maximum1;
    if(maximum2 > highest_maximum)
        highest_maximum = maximum2;

    if(maximum3 > highest_maximum)
        highest_maximum = maximum3;

    if(maximum4 > highest_maximum)
        highest_maximum = maximum4;

    if(maximum5 > highest_maximum)
        highest_maximum = maximum5;

    if(maximum6 > highest_maximum)
        highest_maximum = maximum6;

    if(maximum7 > highest_maximum)
        highest_maximum = maximum7;

    if(maximum8 > highest_maximum)
        highest_maximum = maximum8;

    bool derefinement_condition = false;

    if(highest_maximum <= derefinement_bound)
        derefinement_condition = true;

    return derefinement_condition;
}

void derefine_cell_bunch(struct cell_node *first_bunch_cell, ui32_array *free_sv_positions) {

    assert(first_bunch_cell);

    struct cell_node *cell_before_bunch = first_bunch_cell->previous;
    struct cell_node *cell_after_bunch = first_bunch_cell->next->next->next->next->next->next->next->next;

    uint8_t bunch_level = first_bunch_cell->cell_data.level;
    uint8_t hilbert_shape_number = get_father_bunch_number(first_bunch_cell);
    uint64_t bunch_number = first_bunch_cell->bunch_number;

    // New cell variable (Arithmetic mean between all cells of the bunch).
    real_cpu v_average = 0;
    real_cpu sigma_x_average = 0.0;
    real_cpu sigma_y_average = 0.0;
    real_cpu sigma_z_average = 0.0;

    struct cell_node *auxiliar = first_bunch_cell;

    for(int i = 0; i < 8; i++) {
        v_average += auxiliar->v;
        sigma_x_average += auxiliar->sigma.x;
        sigma_y_average += auxiliar->sigma.y;
        sigma_z_average += auxiliar->sigma.z;

        auxiliar = auxiliar->next;
    }

    v_average /= 8.0;

    sigma_x_average /= 8.0;
    sigma_y_average /= 8.0;
    sigma_z_average /= 8.0;

    // Front northeast node of the bunch becomes the derefined node.
    struct cell_node *new_cell = get_front_northeast_cell(first_bunch_cell);

    if(free_sv_positions && *free_sv_positions) {

        // Free Sv Map positions//////////////////////////////////
        struct cell_node *w = (struct cell_node *)new_cell->neighbours[DOWN];
        struct cell_node *b = (struct cell_node *)new_cell->neighbours[LEFT];
        struct cell_node *s = (struct cell_node *)new_cell->neighbours[BACK];

        struct cell_node *wb = (struct cell_node *)w->neighbours[LEFT];
        struct cell_node *ws = (struct cell_node *)w->neighbours[BACK];
        struct cell_node *sb = (struct cell_node *)s->neighbours[LEFT];

        struct cell_node *wsb = (struct cell_node *)ws->neighbours[LEFT];

        arrput(*free_sv_positions, w->sv_position);
        arrput(*free_sv_positions, b->sv_position);
        arrput(*free_sv_positions, s->sv_position);

        arrput(*free_sv_positions, wb->sv_position);
        arrput(*free_sv_positions, ws->sv_position);
        arrput(*free_sv_positions, sb->sv_position);
        arrput(*free_sv_positions, wsb->sv_position);
        /////////////////////////////////////////////////////////////
    }
    new_cell->previous = cell_before_bunch;
    new_cell->next = cell_after_bunch;

    if(new_cell->previous != 0)
        new_cell->previous->next = new_cell;

    if(new_cell->next != 0)
        new_cell->next->previous = new_cell;

    real_cpu aux_center_x = ((struct cell_node *)(new_cell->neighbours[LEFT]))->center.x;
    real_cpu aux_center_y = ((struct cell_node *)(new_cell->neighbours[DOWN]))->center.y;
    real_cpu aux_center_z = ((struct cell_node *)(new_cell->neighbours[BACK]))->center.z;

    // New geometric variables.
    new_cell->center.x = (new_cell->center.x + aux_center_x) / 2.0f;
    new_cell->center.y = (new_cell->center.y + aux_center_y) / 2.0f;
    new_cell->center.z = (new_cell->center.z + aux_center_z) / 2.0f;

    new_cell->discretization.x = 2.0f * new_cell->discretization.x;
    new_cell->discretization.y = 2.0f * new_cell->discretization.y;
    new_cell->discretization.z = 2.0f * new_cell->discretization.z;

    new_cell->v = v_average;

    new_cell->sigma.x = sigma_x_average;
    new_cell->sigma.y = sigma_y_average;
    new_cell->sigma.z = sigma_z_average;

    new_cell->cell_data.level = bunch_level - (uint8_t)1;
    new_cell->hilbert_shape_number = hilbert_shape_number;
    new_cell->bunch_number = bunch_number / (int)10;

    struct cell_node *front_northeast_cell = new_cell;
    struct cell_node *front_southeast_cell = (struct cell_node *)(front_northeast_cell->neighbours[BACK]);
    struct cell_node *front_northwest_cell = (struct cell_node *)(front_northeast_cell->neighbours[DOWN]);
    struct cell_node *front_southwest_cell = (struct cell_node *)(front_northwest_cell->neighbours[BACK]);
    struct cell_node *back_northeast_cell = (struct cell_node *)(front_northeast_cell->neighbours[LEFT]);
    struct cell_node *back_southeast_cell = (struct cell_node *)(back_northeast_cell->neighbours[BACK]);
    struct cell_node *back_northwest_cell = (struct cell_node *)(back_northeast_cell->neighbours[DOWN]);
    struct cell_node *back_southwest_cell = (struct cell_node *)(back_northwest_cell->neighbours[BACK]);

    // Creation of North Transition Node.
    struct transition_node *north_transition_node = new_transition_node();
    set_transition_node_data(north_transition_node, bunch_level, FRONT, new_cell, front_northwest_cell->neighbours[FRONT],
                             front_northeast_cell->neighbours[FRONT], back_northeast_cell->neighbours[FRONT], back_northwest_cell->neighbours[FRONT]);

    // Creation of South Transition Node.
    struct transition_node *south_transition_node = new_transition_node();
    set_transition_node_data(south_transition_node, bunch_level, BACK, new_cell, front_southwest_cell->neighbours[BACK], front_southeast_cell->neighbours[BACK],
                             back_southeast_cell->neighbours[BACK], back_southwest_cell->neighbours[BACK]);

    // Creation of East Transition Node.
    struct transition_node *east_transition_node = new_transition_node();
    set_transition_node_data(east_transition_node, bunch_level, TOP, new_cell, front_southeast_cell->neighbours[TOP], back_southeast_cell->neighbours[TOP],
                             back_northeast_cell->neighbours[TOP], front_northeast_cell->neighbours[TOP]);

    // Creation of West Transition Node.
    struct transition_node *west_transition_node = new_transition_node();
    set_transition_node_data(west_transition_node, bunch_level, DOWN, new_cell, front_southwest_cell->neighbours[DOWN], back_southwest_cell->neighbours[DOWN],
                             back_northwest_cell->neighbours[DOWN], front_northwest_cell->neighbours[DOWN]);

    // Creation of Front Transition Node.
    struct transition_node *front_transition_node = new_transition_node();
    set_transition_node_data(front_transition_node, bunch_level, RIGHT, new_cell, front_southwest_cell->neighbours[RIGHT],
                             front_southeast_cell->neighbours[RIGHT], front_northeast_cell->neighbours[RIGHT], front_northwest_cell->neighbours[RIGHT]);

    // Creation of Back Transition Node.
    struct transition_node *back_transition_node = new_transition_node();
    set_transition_node_data(back_transition_node, bunch_level, LEFT, new_cell, back_southwest_cell->neighbours[LEFT], back_southeast_cell->neighbours[LEFT],
                             back_northeast_cell->neighbours[LEFT], back_northwest_cell->neighbours[LEFT]);

    // Elimination of the seven unneeded bunch cells.
    free_cell_node(front_northwest_cell);
    free_cell_node(front_southwest_cell);
    free_cell_node(front_southeast_cell);
    free_cell_node(back_northeast_cell);
    free_cell_node(back_northwest_cell);
    free_cell_node(back_southeast_cell);
    free_cell_node(back_southwest_cell);

    // Linking of derefined cell and new transition nodes.
    new_cell->neighbours[FRONT] = north_transition_node;
    new_cell->neighbours[BACK] = south_transition_node;
    new_cell->neighbours[TOP] = east_transition_node;
    new_cell->neighbours[DOWN] = west_transition_node;
    new_cell->neighbours[RIGHT] = front_transition_node;
    new_cell->neighbours[LEFT] = back_transition_node;

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
struct cell_node *get_front_northeast_cell(struct cell_node *first_bunch_cell) {
    if(first_bunch_cell == NULL) {
        fprintf(stderr, "get_front_northeast_cell: Parameter first_bunch_cell is NULL.");
        return NULL;
    }

    struct cell_node *first_cell = first_bunch_cell;
    struct cell_node *second_cell = first_cell->next;
    struct cell_node *third_cell = second_cell->next;
    struct cell_node *fourth_cell = third_cell->next;
    struct cell_node *fifth_cell = fourth_cell->next;
    struct cell_node *sixth_cell = fifth_cell->next;
    struct cell_node *seventh_cell = sixth_cell->next;
    struct cell_node *eighth_cell = seventh_cell->next;

    real_cpu coordinateSum1 = first_cell->center.x + first_cell->center.y + first_cell->center.z;
    real_cpu coordinateSum2 = second_cell->center.x + second_cell->center.y + second_cell->center.z;
    real_cpu coordinateSum3 = third_cell->center.x + third_cell->center.y + third_cell->center.z;
    real_cpu coordinateSum4 = fourth_cell->center.x + fourth_cell->center.y + fourth_cell->center.z;
    real_cpu coordinateSum5 = fifth_cell->center.x + fifth_cell->center.y + fifth_cell->center.z;
    real_cpu coordinateSum6 = sixth_cell->center.x + sixth_cell->center.y + sixth_cell->center.z;
    real_cpu coordinateSum7 = seventh_cell->center.x + seventh_cell->center.y + seventh_cell->center.z;
    real_cpu coordinateSum8 = eighth_cell->center.x + eighth_cell->center.y + eighth_cell->center.z;

    real_cpu maximum;
    struct cell_node *front_northeast_cell = first_cell;
    maximum = coordinateSum1;
    if(coordinateSum2 > maximum) {
        maximum = coordinateSum2;
        front_northeast_cell = second_cell;
    }
    if(coordinateSum3 > maximum) {
        maximum = coordinateSum3;
        front_northeast_cell = third_cell;
    }
    if(coordinateSum4 > maximum) {
        maximum = coordinateSum4;
        front_northeast_cell = fourth_cell;
    }
    if(coordinateSum5 > maximum) {
        maximum = coordinateSum5;
        front_northeast_cell = fifth_cell;
    }
    if(coordinateSum6 > maximum) {
        maximum = coordinateSum6;
        front_northeast_cell = sixth_cell;
    }
    if(coordinateSum7 > maximum) {
        maximum = coordinateSum7;
        front_northeast_cell = seventh_cell;
    }
    if(coordinateSum8 > maximum) {
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
uint8_t get_father_bunch_number(struct cell_node *first_bunch_cell) {

    if(first_bunch_cell == NULL) {
        fprintf(stderr, "Grid::getFatherBunchNumber(): Parameter first_bunch_cell is NULL.");
        exit(10);
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
    uint8_t sum_of_hilbert_shape_numbers = firstCell->hilbert_shape_number + secondCell->hilbert_shape_number + thirdCell->hilbert_shape_number +
                                           fourthCell->hilbert_shape_number + fifthCell->hilbert_shape_number + sixthCell->hilbert_shape_number +
                                           seventhCell->hilbert_shape_number + eighthCell->hilbert_shape_number;

    switch(sum_of_hilbert_shape_numbers) {
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
        if(secondCell->hilbert_shape_number == 10)
            hilbert_shape_number = 4;

        else if(secondCell->hilbert_shape_number == 4)
            hilbert_shape_number = 10;

        else if(secondCell->hilbert_shape_number == 3)
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
 */
void simplify_derefinement(struct transition_node *transition_node) {

    assert(transition_node);

    struct cell_node *derefined_cell = (struct cell_node *)(transition_node->single_connector);
    enum transition_direction direction = transition_node->direction;

    struct transition_node *white_neighbor_cell;
    struct cell_node *black_neighbor_cell;
    struct transition_node *neighbor_transition_node;

    /* All quadruple connectors point to the  same  cell.  It  means  that  two
     * transition nodes of the same level have been connected. This connection is
     * not allowed and it is simplified through a direct connection between the two
     * cell nodes pointed by each transition node. Then, both transition nodes
     * are eliminated.
     */
    if((transition_node->quadruple_connector1 == transition_node->quadruple_connector2) &&
       (transition_node->quadruple_connector1 == transition_node->quadruple_connector3) &&
       (transition_node->quadruple_connector1 == transition_node->quadruple_connector4)) {
        // PS: No matter which quadruple connector is chosen in the line below,
        // because each one of them points to the same cell.
        neighbor_transition_node = (struct transition_node *)(transition_node->quadruple_connector1);

        enum cell_type neighbor_cell_type = ((struct basic_cell_data *)(neighbor_transition_node->single_connector))->type;

        if(VALID_SIMPLE_DIRECTION(direction)) {
            derefined_cell->neighbours[direction] = neighbor_transition_node->single_connector;
        } else {
            fprintf(stderr, "simplify_derefinement(). Invalid cell direction %d! Exiting...\n", direction);
            exit(10);
        }

        if(neighbor_cell_type == CELL_NODE) {
            black_neighbor_cell = (struct cell_node *)(neighbor_transition_node->single_connector);
            black_neighbor_cell->neighbours[get_inverse_direction(direction)] = derefined_cell;
        } else {
            white_neighbor_cell = (struct transition_node *)(neighbor_transition_node->single_connector);

            if(white_neighbor_cell->single_connector == neighbor_transition_node) {
                white_neighbor_cell->single_connector = derefined_cell;
            } else if(white_neighbor_cell->quadruple_connector1 == neighbor_transition_node) {
                white_neighbor_cell->quadruple_connector1 = derefined_cell;
            } else if(white_neighbor_cell->quadruple_connector2 == neighbor_transition_node) {
                white_neighbor_cell->quadruple_connector2 = derefined_cell;
            } else if(white_neighbor_cell->quadruple_connector3 == neighbor_transition_node) {
                white_neighbor_cell->quadruple_connector3 = derefined_cell;
            } else if(white_neighbor_cell->quadruple_connector4 == neighbor_transition_node) {
                white_neighbor_cell->quadruple_connector4 = derefined_cell;
            }
        }

        free(neighbor_transition_node);
        free(transition_node);
    } else {

        void *quadruple_connector[4];

        /* Connects outside to the transition node, if this was not deleted. That is,
         * the quadruple connectors point to different cells. */
        quadruple_connector[0] = transition_node->quadruple_connector1;
        quadruple_connector[1] = transition_node->quadruple_connector2;
        quadruple_connector[2] = transition_node->quadruple_connector3;
        quadruple_connector[3] = transition_node->quadruple_connector4;

        for(int i = 0; i < 4; i++) {
            enum cell_type connector_type = ((struct basic_cell_data *)(quadruple_connector[i]))->type;
            if(connector_type == TRANSITION_NODE) {
                white_neighbor_cell = (struct transition_node *)(quadruple_connector[i]);
                white_neighbor_cell->single_connector = transition_node;
            } else if(connector_type == CELL_NODE) {
                black_neighbor_cell = (struct cell_node *)(quadruple_connector[i]);
                black_neighbor_cell->neighbours[get_inverse_direction(direction)] = transition_node;
            }
        }
    }
}
