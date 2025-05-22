//
// Created by sachetto on 30/09/17.
//

#include "../../3dparty/stb_ds.h"
#include "cell.h"
#include "string.h"

#include <assert.h>

#define SET_TRANSITION_NODE(t_node, direction)                                                                                                                 \
    do {                                                                                                                                                       \
        enum cell_type node_type = ((struct basic_cell_data *)t_node->single_connector)->type;                                                                 \
        if(node_type == CELL_NODE) {                                                                                                                           \
            neighbour_cell_node = (struct cell_node *)(t_node->single_connector);                                                                              \
            neighbour_cell_node->neighbours[direction] = t_node;                                                                                               \
        } else if(node_type == TRANSITION_NODE) {                                                                                                              \
            neighbour_transition_node = (struct transition_node *)(t_node->single_connector);                                                                  \
                                                                                                                                                               \
            if(neighbour_transition_node->single_connector == right_front_top_sub_cell)                                                                        \
                neighbour_transition_node->single_connector = t_node;                                                                                          \
                                                                                                                                                               \
            else if(neighbour_transition_node->quadruple_connector1 == right_front_top_sub_cell)                                                               \
                neighbour_transition_node->quadruple_connector1 = t_node;                                                                                      \
                                                                                                                                                               \
            else if(neighbour_transition_node->quadruple_connector2 == right_front_top_sub_cell)                                                               \
                neighbour_transition_node->quadruple_connector2 = t_node;                                                                                      \
                                                                                                                                                               \
            else if(neighbour_transition_node->quadruple_connector3 == right_front_top_sub_cell)                                                               \
                neighbour_transition_node->quadruple_connector3 = t_node;                                                                                      \
                                                                                                                                                               \
            else if(neighbour_transition_node->quadruple_connector4 == right_front_top_sub_cell)                                                               \
                neighbour_transition_node->quadruple_connector4 = t_node;                                                                                      \
        }                                                                                                                                                      \
    } while(0)

void refine_cell(struct cell_node *cell, ui32_array free_sv_positions, ui32_array *refined_this_step) {

    assert(cell);

    struct transition_node *top_transition_node;
    struct transition_node *front_transition_node;
    struct transition_node *down_transition_node;
    struct transition_node *back_transition_node;
    struct transition_node *right_transition_node;
    struct transition_node *left_transition_node;

    struct cell_node *right_front_top_sub_cell, *right_front_down_sub_cell, *right_back_top_sub_cell, *right_back_down_sub_cell, *left_front_top_sub_cell,
        *left_front_down_sub_cell, *left_back_top_sub_cell, *left_back_down_sub_cell;

    uint8_t number_of_hilbert_shape;

    real_cpu cell_center_x = cell->center.x, cell_center_y = cell->center.y, cell_center_z = cell->center.z, cell_half_side_x = cell->discretization.x / 2.0f,
             cell_half_side_y = cell->discretization.y / 2.0f, cell_half_side_z = cell->discretization.z / 2.0f,
             cell_quarter_side_x = cell->discretization.x / 4.0f, cell_quarter_side_y = cell->discretization.y / 4.0f,
             cell_quarter_side_z = cell->discretization.z / 4.0f;

    uint64_t old_bunch_number = cell->bunch_number;

    // Creation of the front northeast cell.
    // This cell, which is to be refined,
    // becomes the front northeast cell of the new bunch.
    right_front_top_sub_cell = cell;
    right_front_top_sub_cell->cell_data.level = cell->cell_data.level + (uint8_t)1;
    right_front_top_sub_cell->discretization.x = cell_half_side_x;
    right_front_top_sub_cell->discretization.y = cell_half_side_y;
    right_front_top_sub_cell->discretization.z = cell_half_side_z;

    right_front_top_sub_cell->center.x = cell_center_x + cell_quarter_side_x;
    right_front_top_sub_cell->center.y = cell_center_y + cell_quarter_side_y;
    right_front_top_sub_cell->center.z = cell_center_z + cell_quarter_side_z;

    right_front_top_sub_cell->bunch_number = old_bunch_number * 10 + 1;

    if(refined_this_step && *refined_this_step) {
        arrput(*refined_this_step, right_front_top_sub_cell->sv_position);
    }

    left_front_top_sub_cell = new_cell_node();
    set_refined_cell_data(left_front_top_sub_cell, right_front_top_sub_cell, POINT3D(cell_half_side_x, cell_half_side_y, cell_half_side_z),
                          POINT3D(cell_center_x - cell_quarter_side_x, cell_center_y + cell_quarter_side_y, cell_center_z + cell_quarter_side_z),
                          old_bunch_number * 10 + 2, free_sv_positions, refined_this_step);

    left_front_down_sub_cell = new_cell_node();
    set_refined_cell_data(left_front_down_sub_cell, right_front_top_sub_cell, POINT3D(cell_half_side_x, cell_half_side_y, cell_half_side_z),
                          POINT3D(cell_center_x - cell_quarter_side_x, cell_center_y - cell_quarter_side_y, cell_center_z + cell_quarter_side_z),
                          old_bunch_number * 10 + 3, free_sv_positions, refined_this_step);

    right_front_down_sub_cell = new_cell_node();
    set_refined_cell_data(right_front_down_sub_cell, right_front_top_sub_cell, POINT3D(cell_half_side_x, cell_half_side_y, cell_half_side_z),
                          POINT3D(cell_center_x + cell_quarter_side_x, cell_center_y - cell_quarter_side_y, cell_center_z + cell_quarter_side_z),
                          old_bunch_number * 10 + 4, free_sv_positions, refined_this_step);

    right_back_down_sub_cell = new_cell_node();
    set_refined_cell_data(right_back_down_sub_cell, right_front_top_sub_cell, POINT3D(cell_half_side_x, cell_half_side_y, cell_half_side_z),
                          POINT3D(cell_center_x + cell_quarter_side_x, cell_center_y - cell_quarter_side_y, cell_center_z - cell_quarter_side_z),
                          old_bunch_number * 10 + 5, free_sv_positions, refined_this_step);

    left_back_down_sub_cell = new_cell_node();
    set_refined_cell_data(left_back_down_sub_cell, right_front_top_sub_cell, POINT3D(cell_half_side_x, cell_half_side_y, cell_half_side_z),
                          POINT3D(cell_center_x - cell_quarter_side_x, cell_center_y - cell_quarter_side_y, cell_center_z - cell_quarter_side_z),
                          old_bunch_number * 10 + 6, free_sv_positions, refined_this_step);

    left_back_top_sub_cell = new_cell_node();
    set_refined_cell_data(left_back_top_sub_cell, right_front_top_sub_cell, POINT3D(cell_half_side_x, cell_half_side_y, cell_half_side_z),
                          POINT3D(cell_center_x - cell_quarter_side_x, cell_center_y + cell_quarter_side_y, cell_center_z - cell_quarter_side_z),
                          old_bunch_number * 10 + 7, free_sv_positions, refined_this_step);

    right_back_top_sub_cell = new_cell_node();
    set_refined_cell_data(right_back_top_sub_cell, right_front_top_sub_cell, POINT3D(cell_half_side_x, cell_half_side_y, cell_half_side_z),
                          POINT3D(cell_center_x + cell_quarter_side_x, cell_center_y + cell_quarter_side_y, cell_center_z - cell_quarter_side_z),
                          old_bunch_number * 10 + 8, free_sv_positions, refined_this_step);

    front_transition_node = new_transition_node();
    set_refined_transition_node_data(front_transition_node, right_front_top_sub_cell, FRONT);

    back_transition_node = new_transition_node();
    set_refined_transition_node_data(back_transition_node, right_front_top_sub_cell, BACK);

    top_transition_node = new_transition_node();
    set_refined_transition_node_data(top_transition_node, right_front_top_sub_cell, TOP);

    down_transition_node = new_transition_node();
    set_refined_transition_node_data(down_transition_node, right_front_top_sub_cell, DOWN);

    right_transition_node = new_transition_node();
    set_refined_transition_node_data(right_transition_node, right_front_top_sub_cell, RIGHT);

    left_transition_node = new_transition_node();
    set_refined_transition_node_data(left_transition_node, right_front_top_sub_cell, LEFT);

    //    top_right_transition_node  = new_transition_node();
    //    set_refined_transition_node_data(top_right_transition_node, right_front_top_sub_cell, TOP_RIGHT);
    //
    //    top_left_transition_node  = new_transition_node();
    //    set_refined_transition_node_data(top_left_transition_node, right_front_top_sub_cell, TOP_LEFT);
    //
    //    top_front_transition_node  = new_transition_node();
    //    set_refined_transition_node_data(top_front_transition_node, right_front_top_sub_cell, TOP_FRONT);
    //
    //    top_back_transition_node  = new_transition_node();
    //    set_refined_transition_node_data(top_back_transition_node, right_front_top_sub_cell, TOP_BACK);
    //
    //    down_right_transition_node  = new_transition_node();
    //    set_refined_transition_node_data(down_right_transition_node, right_front_top_sub_cell, DOWN_RIGHT);
    //
    //    down_left_transition_node  = new_transition_node();
    //    set_refined_transition_node_data(down_left_transition_node, right_front_top_sub_cell, DOWN_LEFT);
    //
    //    down_front_transition_node  = new_transition_node();
    //    set_refined_transition_node_data(down_front_transition_node, right_front_top_sub_cell, DOWN_FRONT);
    //
    //    down_back_transition_node  = new_transition_node();
    //    set_refined_transition_node_data(down_back_transition_node, right_front_top_sub_cell, DOWN_BACK);
    //
    //    right_front_transition_node  = new_transition_node();
    //    set_refined_transition_node_data(right_front_transition_node, right_front_top_sub_cell, RIGHT_FRONT);
    //
    //    right_back_transition_node  = new_transition_node();
    //    set_refined_transition_node_data(right_back_transition_node, right_front_top_sub_cell, RIGHT_BACK);
    //
    //    left_front_transition_node  = new_transition_node();
    //    set_refined_transition_node_data(left_front_transition_node, right_front_top_sub_cell, LEFT_FRONT);
    //
    //    left_back_transition_node  = new_transition_node();
    //    set_refined_transition_node_data(left_back_transition_node, right_front_top_sub_cell, LEFT_BACK);
    //
    //    front_left_top_transition_node  = new_transition_node();
    //    set_refined_transition_node_data(front_left_top_transition_node, right_front_top_sub_cell, FRONT_LEFT_TOP);
    //
    //    front_left_down_transition_node  = new_transition_node();
    //    set_refined_transition_node_data(front_left_down_transition_node, right_front_top_sub_cell, FRONT_LEFT_DOWN);
    //
    //    front_right_top_transition_node  = new_transition_node();
    //    set_refined_transition_node_data(front_right_top_transition_node, right_front_top_sub_cell, FRONT_RIGHT_TOP);
    //
    //    front_right_down_transition_node  = new_transition_node();
    //    set_refined_transition_node_data(front_right_down_transition_node, right_front_top_sub_cell, FRONT_RIGHT_DOWN);
    //
    //    back_left_top_transition_node  = new_transition_node();
    //    set_refined_transition_node_data(back_left_top_transition_node, right_front_top_sub_cell, BACK_LEFT_TOP);
    //
    //    back_left_down_transition_node  = new_transition_node();
    //    set_refined_transition_node_data(back_left_down_transition_node, right_front_top_sub_cell, BACK_LEFT_DOWN);
    //
    //    back_right_top_transition_node  = new_transition_node();
    //    set_refined_transition_node_data(back_right_top_transition_node, right_front_top_sub_cell, BACK_RIGHT_TOP);
    //
    //    back_right_down_transition_node  = new_transition_node();
    //    set_refined_transition_node_data(back_right_down_transition_node, right_front_top_sub_cell, BACK_RIGHT_DOWN);

    // Linking of new cell nodes and transition nodes.
    right_front_top_sub_cell->neighbours[FRONT] = front_transition_node;
    right_front_top_sub_cell->neighbours[BACK] = right_back_top_sub_cell;
    right_front_top_sub_cell->neighbours[TOP] = top_transition_node;
    right_front_top_sub_cell->neighbours[DOWN] = right_front_down_sub_cell;
    right_front_top_sub_cell->neighbours[RIGHT] = right_transition_node;
    right_front_top_sub_cell->neighbours[LEFT] = left_front_top_sub_cell;

    //    right_front_top_sub_cell->neighbours[TOP]_right = top_right_transition_node;
    //    right_front_top_sub_cell->neighbours[TOP]_left = top_left_transition_node;
    //    right_front_top_sub_cell->neighbours[TOP]_front = top_front_transition_node;
    //    right_front_top_sub_cell->neighbours[TOP]_back = top_back_transition_node;
    //    right_front_top_sub_cell->neighbours[DOWN]_right = right_transition_node;
    //    right_front_top_sub_cell->neighbours[DOWN]_left = left_front_down_sub_cell;
    //    right_front_top_sub_cell->neighbours[DOWN]_front = front_transition_node;
    //    right_front_top_sub_cell->neighbours[DOWN]_back = right_back_down_sub_cell;
    //    right_front_top_sub_cell->neighbours[RIGHT]_front = right_front_transition_node ;
    //    right_front_top_sub_cell->neighbours[RIGHT]_back = right_back_transition_node;
    //    right_front_top_sub_cell->neighbours[LEFT]_front = left_front_transition_node;
    //    right_front_top_sub_cell->neighbours[LEFT]_back = left_back_top_sub_cell;
    //    right_front_top_sub_cell->neighbours[FRONT]_left_top = front_left_top_transition_node;
    //    right_front_top_sub_cell->neighbours[FRONT]_left_down = front_left_down_transition_node;
    //    right_front_top_sub_cell->neighbours[FRONT]_right_top = front_right_top_transition_node;
    //    right_front_top_sub_cell->neighbours[FRONT]_right_down = front_right_down_transition_node;
    //    right_front_top_sub_cell->neighbours[BACK]_left_top = back_left_top_transition_node;
    //    right_front_top_sub_cell->neighbours[BACK]_left_down = left_back_down_sub_cell;
    //    right_front_top_sub_cell->neighbours[BACK]_right_top = back_right_top_transition_node;
    //    right_front_top_sub_cell->neighbours[BACK]_right_down = back_right_down_transition_node;

    left_front_top_sub_cell->neighbours[FRONT] = front_transition_node;
    left_front_top_sub_cell->neighbours[BACK] = left_back_top_sub_cell;
    left_front_top_sub_cell->neighbours[TOP] = top_transition_node;
    left_front_top_sub_cell->neighbours[DOWN] = left_front_down_sub_cell;
    left_front_top_sub_cell->neighbours[RIGHT] = right_front_top_sub_cell;
    left_front_top_sub_cell->neighbours[LEFT] = left_transition_node;

    //    left_front_top_sub_cell->neighbours[TOP]_right = top_right_transition_node;
    //    left_front_top_sub_cell->neighbours[TOP]_left = top_left_transition_node;
    //    left_front_top_sub_cell->neighbours[TOP]_front = top_front_transition_node;
    //    left_front_top_sub_cell->neighbours[TOP]_back = top_back_transition_node;
    //    left_front_top_sub_cell->neighbours[DOWN]_right = right_front_down_sub_cell;
    //    left_front_top_sub_cell->neighbours[DOWN]_left = down_left_transition_node;
    //    left_front_top_sub_cell->neighbours[DOWN]_front = front_transition_node;
    //    left_front_top_sub_cell->neighbours[DOWN]_back = left_back_down_sub_cell;
    //    left_front_top_sub_cell->neighbours[RIGHT]_front = right_front_transition_node;
    //    left_front_top_sub_cell->neighbours[RIGHT]_back = right_back_top_sub_cell;
    //    left_front_top_sub_cell->neighbours[LEFT]_front = left_front_transition_node;
    //    left_front_top_sub_cell->neighbours[LEFT]_back = left_back_transition_node;
    //    left_front_top_sub_cell->neighbours[FRONT]_left_top = front_left_top_transition_node;
    //    left_front_top_sub_cell->neighbours[FRONT]_left_down = front_left_down_transition_node;
    //    left_front_top_sub_cell->neighbours[FRONT]_right_top = front_right_top_transition_node;
    //    left_front_top_sub_cell->neighbours[FRONT]_right_down = front_right_down_transition_node;
    //    left_front_top_sub_cell->neighbours[BACK]_left_top = back_left_top_transition_node;
    //    left_front_top_sub_cell->neighbours[BACK]_left_down = back_left_down_transition_node;
    //    left_front_top_sub_cell->neighbours[BACK]_right_top = back_right_top_transition_node;
    //    left_front_top_sub_cell->neighbours[BACK]_right_down = right_back_down_sub_cell;

    left_front_down_sub_cell->neighbours[FRONT] = front_transition_node;
    left_front_down_sub_cell->neighbours[BACK] = left_back_down_sub_cell;
    left_front_down_sub_cell->neighbours[TOP] = left_front_top_sub_cell;
    left_front_down_sub_cell->neighbours[DOWN] = down_transition_node;
    left_front_down_sub_cell->neighbours[RIGHT] = right_front_down_sub_cell;
    left_front_down_sub_cell->neighbours[LEFT] = left_transition_node;

    //    left_front_down_sub_cell->neighbours[TOP]_right = right_front_top_sub_cell;
    //    left_front_down_sub_cell->neighbours[TOP]_left = left_transition_node;
    //    left_front_down_sub_cell->neighbours[TOP]_front = front_transition_node;
    //    left_front_down_sub_cell->neighbours[TOP]_back = left_back_top_sub_cell;
    //    left_front_down_sub_cell->neighbours[DOWN]_right = down_right_transition_node;
    //    left_front_down_sub_cell->neighbours[DOWN]_left = down_left_transition_node;
    //    left_front_down_sub_cell->neighbours[DOWN]_front = down_front_transition_node;
    //    left_front_down_sub_cell->neighbours[DOWN]_back = down_back_transition_node;
    //    left_front_down_sub_cell->neighbours[RIGHT]_front = right_front_transition_node;
    //    left_front_down_sub_cell->neighbours[RIGHT]_back = right_back_down_sub_cell;
    //    left_front_down_sub_cell->neighbours[LEFT]_front = left_front_transition_node;
    //    left_front_down_sub_cell->neighbours[LEFT]_back = left_back_transition_node;
    //    left_front_down_sub_cell->neighbours[FRONT]_left_top = front_left_top_transition_node;
    //    left_front_down_sub_cell->neighbours[FRONT]_left_down = front_left_down_transition_node;
    //    left_front_down_sub_cell->neighbours[FRONT]_right_top = front_right_top_transition_node;
    //    left_front_down_sub_cell->neighbours[FRONT]_right_down = front_right_down_transition_node;
    //    left_front_down_sub_cell->neighbours[BACK]_left_top = back_left_top_transition_node;
    //    left_front_down_sub_cell->neighbours[BACK]_left_down = back_left_down_transition_node;
    //    left_front_down_sub_cell->neighbours[BACK]_right_top = right_back_top_sub_cell;
    //    left_front_down_sub_cell->neighbours[BACK]_right_down = back_right_down_transition_node;

    right_front_down_sub_cell->neighbours[FRONT] = front_transition_node;
    right_front_down_sub_cell->neighbours[BACK] = right_back_down_sub_cell;
    right_front_down_sub_cell->neighbours[TOP] = right_front_top_sub_cell;
    right_front_down_sub_cell->neighbours[DOWN] = down_transition_node;
    right_front_down_sub_cell->neighbours[RIGHT] = right_transition_node;
    right_front_down_sub_cell->neighbours[LEFT] = left_front_down_sub_cell;

    //    right_front_down_sub_cell->neighbours[TOP]_right = right_transition_node;
    //    right_front_down_sub_cell->neighbours[TOP]_left = left_front_top_sub_cell;
    //    right_front_down_sub_cell->neighbours[TOP]_front = front_transition_node;
    //    right_front_down_sub_cell->neighbours[TOP]_back = right_back_top_sub_cell;
    //    right_front_down_sub_cell->neighbours[DOWN]_right = down_right_transition_node;
    //    right_front_down_sub_cell->neighbours[DOWN]_left = down_transition_node;
    //    right_front_down_sub_cell->neighbours[DOWN]_front = down_front_transition_node;
    //    right_front_down_sub_cell->neighbours[DOWN]_back = down_back_transition_node;
    //    right_front_down_sub_cell->neighbours[RIGHT]_front = right_front_transition_node;
    //    right_front_down_sub_cell->neighbours[RIGHT]_back = right_back_transition_node;
    //    right_front_down_sub_cell->neighbours[LEFT]_front = left_front_transition_node;
    //    right_front_down_sub_cell->neighbours[LEFT]_back = left_back_down_sub_cell;
    //    right_front_down_sub_cell->neighbours[FRONT]_left_top = front_left_top_transition_node;
    //    right_front_down_sub_cell->neighbours[FRONT]_left_down = front_left_down_transition_node;
    //    right_front_down_sub_cell->neighbours[FRONT]_right_top = front_right_top_transition_node;
    //    right_front_down_sub_cell->neighbours[FRONT]_right_down = front_right_down_transition_node;
    //    right_front_down_sub_cell->neighbours[BACK]_left_top = left_back_top_sub_cell;
    //    right_front_down_sub_cell->neighbours[BACK]_left_down = back_left_down_transition_node;
    //    right_front_down_sub_cell->neighbours[BACK]_right_top = back_right_top_transition_node;
    //    right_front_down_sub_cell->neighbours[BACK]_right_down = back_right_down_transition_node;

    right_back_down_sub_cell->neighbours[FRONT] = right_front_down_sub_cell;
    right_back_down_sub_cell->neighbours[BACK] = back_transition_node;
    right_back_down_sub_cell->neighbours[TOP] = right_back_top_sub_cell;
    right_back_down_sub_cell->neighbours[DOWN] = down_transition_node;
    right_back_down_sub_cell->neighbours[RIGHT] = right_transition_node;
    right_back_down_sub_cell->neighbours[LEFT] = left_back_down_sub_cell;

    //    right_back_down_sub_cell->neighbours[TOP]_right = right_transition_node;
    //    right_back_down_sub_cell->neighbours[TOP]_left = left_back_top_sub_cell;
    //    right_back_down_sub_cell->neighbours[TOP]_front = right_front_top_sub_cell;
    //    right_back_down_sub_cell->neighbours[TOP]_back = back_transition_node;
    //    right_back_down_sub_cell->neighbours[DOWN]_right = down_right_transition_node;
    //    right_back_down_sub_cell->neighbours[DOWN]_left = down_transition_node;
    //    right_back_down_sub_cell->neighbours[DOWN]_front = down_front_transition_node;
    //    right_back_down_sub_cell->neighbours[DOWN]_back = down_back_transition_node;
    //    right_back_down_sub_cell->neighbours[RIGHT]_front = right_front_transition_node;
    //    right_back_down_sub_cell->neighbours[RIGHT]_back = right_back_transition_node;
    //    right_back_down_sub_cell->neighbours[LEFT]_front = left_front_down_sub_cell;
    //    right_back_down_sub_cell->neighbours[LEFT]_back = left_back_transition_node;
    //    right_back_down_sub_cell->neighbours[FRONT]_left_top = left_front_top_sub_cell;
    //    right_back_down_sub_cell->neighbours[FRONT]_left_down = front_left_down_transition_node;
    //    right_back_down_sub_cell->neighbours[FRONT]_right_top = front_right_top_transition_node;
    //    right_back_down_sub_cell->neighbours[FRONT]_right_down = front_right_down_transition_node;
    //    right_back_down_sub_cell->neighbours[BACK]_left_top = back_left_top_transition_node;
    //    right_back_down_sub_cell->neighbours[BACK]_left_down = back_left_down_transition_node;
    //    right_back_down_sub_cell->neighbours[BACK]_right_top = back_right_top_transition_node;
    //    right_back_down_sub_cell->neighbours[BACK]_right_down = back_right_down_transition_node;

    left_back_down_sub_cell->neighbours[FRONT] = left_front_down_sub_cell;
    left_back_down_sub_cell->neighbours[BACK] = back_transition_node;
    left_back_down_sub_cell->neighbours[TOP] = left_back_top_sub_cell;
    left_back_down_sub_cell->neighbours[DOWN] = down_transition_node;
    left_back_down_sub_cell->neighbours[RIGHT] = right_back_down_sub_cell;
    left_back_down_sub_cell->neighbours[LEFT] = left_transition_node;

    //    left_back_down_sub_cell->neighbours[TOP]_right = right_back_top_sub_cell;
    //    left_back_down_sub_cell->neighbours[TOP]_left = left_transition_node;
    //    left_back_down_sub_cell->neighbours[TOP]_front = left_front_top_sub_cell;
    //    left_back_down_sub_cell->neighbours[TOP]_back = back_transition_node;
    //    left_back_down_sub_cell->neighbours[DOWN]_right = down_right_transition_node;
    //    left_back_down_sub_cell->neighbours[DOWN]_left = down_left_transition_node;
    //    left_back_down_sub_cell->neighbours[DOWN]_front = down_front_transition_node;
    //    left_back_down_sub_cell->neighbours[DOWN]_back = down_back_transition_node;
    //    left_back_down_sub_cell->neighbours[RIGHT]_front = right_front_down_sub_cell;
    //    left_back_down_sub_cell->neighbours[RIGHT]_back = right_back_transition_node;
    //    left_back_down_sub_cell->neighbours[LEFT]_front = left_front_transition_node;
    //    left_back_down_sub_cell->neighbours[LEFT]_back = left_back_transition_node;
    //    left_back_down_sub_cell->neighbours[FRONT]_left_top = front_left_top_transition_node;
    //    left_back_down_sub_cell->neighbours[FRONT]_left_down = front_left_down_transition_node;
    //    left_back_down_sub_cell->neighbours[FRONT]_right_top = right_front_top_sub_cell;
    //    left_back_down_sub_cell->neighbours[FRONT]_right_down = front_right_down_transition_node;
    //    left_back_down_sub_cell->neighbours[BACK]_left_top = back_left_top_transition_node;
    //    left_back_down_sub_cell->neighbours[BACK]_left_down = back_left_down_transition_node;
    //    left_back_down_sub_cell->neighbours[BACK]_right_top = back_right_top_transition_node;
    //    left_back_down_sub_cell->neighbours[BACK]_right_down = back_right_down_transition_node;

    left_back_top_sub_cell->neighbours[FRONT] = left_front_top_sub_cell;
    left_back_top_sub_cell->neighbours[BACK] = back_transition_node;
    left_back_top_sub_cell->neighbours[TOP] = top_transition_node;
    left_back_top_sub_cell->neighbours[DOWN] = left_back_down_sub_cell;
    left_back_top_sub_cell->neighbours[RIGHT] = right_back_top_sub_cell;
    left_back_top_sub_cell->neighbours[LEFT] = left_transition_node;

    //    left_back_top_sub_cell->neighbours[TOP]_right = top_right_transition_node;
    //    left_back_top_sub_cell->neighbours[TOP]_left = top_left_transition_node;
    //    left_back_top_sub_cell->neighbours[TOP]_front = top_front_transition_node;
    //    left_back_top_sub_cell->neighbours[TOP]_back = top_back_transition_node;
    //    left_back_top_sub_cell->neighbours[DOWN]_right = right_back_down_sub_cell;
    //    left_back_top_sub_cell->neighbours[DOWN]_left = down_left_transition_node;
    //    left_back_top_sub_cell->neighbours[DOWN]_front = left_front_down_sub_cell;
    //    left_back_top_sub_cell->neighbours[DOWN]_back = back_transition_node;
    //    left_back_top_sub_cell->neighbours[RIGHT]_front = right_front_top_sub_cell;
    //    left_back_top_sub_cell->neighbours[RIGHT]_back = right_back_transition_node;
    //    left_back_top_sub_cell->neighbours[LEFT]_front = left_front_transition_node;
    //    left_back_top_sub_cell->neighbours[LEFT]_back = left_back_transition_node;
    //    left_back_top_sub_cell->neighbours[FRONT]_left_top = front_left_top_transition_node;
    //    left_back_top_sub_cell->neighbours[FRONT]_left_down = front_left_down_transition_node;
    //    left_back_top_sub_cell->neighbours[FRONT]_right_top = front_right_top_transition_node;
    //    left_back_top_sub_cell->neighbours[FRONT]_right_down = right_front_down_sub_cell;
    //    left_back_top_sub_cell->neighbours[BACK]_left_top = back_left_top_transition_node;
    //    left_back_top_sub_cell->neighbours[BACK]_left_down = back_left_down_transition_node;
    //    left_back_top_sub_cell->neighbours[BACK]_right_top = back_right_top_transition_node;
    //    left_back_top_sub_cell->neighbours[BACK]_right_down = back_right_down_transition_node;

    right_back_top_sub_cell->neighbours[FRONT] = right_front_top_sub_cell;
    right_back_top_sub_cell->neighbours[BACK] = back_transition_node;
    right_back_top_sub_cell->neighbours[TOP] = top_transition_node;
    right_back_top_sub_cell->neighbours[DOWN] = right_back_down_sub_cell;
    right_back_top_sub_cell->neighbours[RIGHT] = right_transition_node;
    right_back_top_sub_cell->neighbours[LEFT] = left_back_top_sub_cell;

    //    right_back_top_sub_cell->neighbours[TOP]_right = top_right_transition_node;
    //    right_back_top_sub_cell->neighbours[TOP]_left = top_left_transition_node;
    //    right_back_top_sub_cell->neighbours[TOP]_front = top_front_transition_node;
    //    right_back_top_sub_cell->neighbours[TOP]_back = top_back_transition_node;
    //    right_back_top_sub_cell->neighbours[DOWN]_right = right_transition_node;
    //    right_back_top_sub_cell->neighbours[DOWN]_left = left_back_down_sub_cell;
    //    right_back_top_sub_cell->neighbours[DOWN]_front = right_front_down_sub_cell;
    //    right_back_top_sub_cell->neighbours[DOWN]_back = back_transition_node;
    //    right_back_top_sub_cell->neighbours[RIGHT]_front = right_front_transition_node;
    //    right_back_top_sub_cell->neighbours[RIGHT]_back = right_back_transition_node;
    //    right_back_top_sub_cell->neighbours[LEFT]_front = left_front_top_sub_cell;
    //    right_back_top_sub_cell->neighbours[LEFT]_back = left_back_transition_node;
    //    right_back_top_sub_cell->neighbours[FRONT]_left_top = front_left_top_transition_node;
    //    right_back_top_sub_cell->neighbours[FRONT]_left_down = left_front_down_sub_cell;
    //    right_back_top_sub_cell->neighbours[FRONT]_right_top = front_right_top_transition_node;
    //    right_back_top_sub_cell->neighbours[FRONT]_right_down = front_right_down_transition_node;
    //    right_back_top_sub_cell->neighbours[BACK]_left_top = back_left_top_transition_node;
    //    right_back_top_sub_cell->neighbours[BACK]_left_down = back_left_down_transition_node;
    //    right_back_top_sub_cell->neighbours[BACK]_right_top = back_right_top_transition_node;
    //    right_back_top_sub_cell->neighbours[BACK]_right_down = back_right_down_transition_node;

    right_transition_node->quadruple_connector1 = right_back_down_sub_cell;
    right_transition_node->quadruple_connector2 = right_back_top_sub_cell;
    right_transition_node->quadruple_connector3 = right_front_top_sub_cell;
    right_transition_node->quadruple_connector4 = right_front_down_sub_cell;

    left_transition_node->quadruple_connector1 = left_back_down_sub_cell;
    left_transition_node->quadruple_connector2 = left_back_top_sub_cell;
    left_transition_node->quadruple_connector3 = left_front_top_sub_cell;
    left_transition_node->quadruple_connector4 = left_front_down_sub_cell;

    down_transition_node->quadruple_connector1 = right_back_down_sub_cell;
    down_transition_node->quadruple_connector2 = left_back_down_sub_cell;
    down_transition_node->quadruple_connector3 = left_front_down_sub_cell;
    down_transition_node->quadruple_connector4 = right_front_down_sub_cell;

    top_transition_node->quadruple_connector1 = right_back_top_sub_cell;
    top_transition_node->quadruple_connector2 = left_back_top_sub_cell;
    top_transition_node->quadruple_connector3 = left_front_top_sub_cell;
    top_transition_node->quadruple_connector4 = right_front_top_sub_cell;

    front_transition_node->quadruple_connector1 = right_front_down_sub_cell;
    front_transition_node->quadruple_connector2 = right_front_top_sub_cell;
    front_transition_node->quadruple_connector3 = left_front_top_sub_cell;
    front_transition_node->quadruple_connector4 = left_front_down_sub_cell;

    back_transition_node->quadruple_connector1 = right_back_down_sub_cell;
    back_transition_node->quadruple_connector2 = right_back_top_sub_cell;
    back_transition_node->quadruple_connector3 = left_back_top_sub_cell;
    back_transition_node->quadruple_connector4 = left_back_down_sub_cell;

    // Linking bunch neighbour cells to the transition nodes just created.
    struct cell_node *neighbour_cell_node = NULL;
    struct transition_node *neighbour_transition_node = NULL;

    SET_TRANSITION_NODE(top_transition_node, DOWN);
    SET_TRANSITION_NODE(down_transition_node, TOP);
    SET_TRANSITION_NODE(front_transition_node, BACK);
    SET_TRANSITION_NODE(back_transition_node, FRONT);
    SET_TRANSITION_NODE(right_transition_node, LEFT);
    SET_TRANSITION_NODE(left_transition_node, RIGHT);

    /*==========================================================================
                ORDERING OF CELL NODES THROUGH HILBERT'S CURVE
    ==========================================================================*/
    number_of_hilbert_shape = right_front_top_sub_cell->hilbert_shape_number;

    if(number_of_hilbert_shape == 0) {
        /* Shape 0
                            _______
                           /      /      b: begin
                          /     b/       e: end
                          |  ______
                          | /     /
                          |/    e/
         */

        right_front_top_sub_cell->hilbert_shape_number = 1;
        left_front_top_sub_cell->hilbert_shape_number = 2;
        left_front_down_sub_cell->hilbert_shape_number = 2;
        right_front_down_sub_cell->hilbert_shape_number = 3;
        right_back_down_sub_cell->hilbert_shape_number = 3;
        left_back_down_sub_cell->hilbert_shape_number = 4;
        left_back_top_sub_cell->hilbert_shape_number = 4;
        right_back_top_sub_cell->hilbert_shape_number = 5;

        right_back_top_sub_cell->next = right_front_top_sub_cell->next;
        right_front_top_sub_cell->next = left_front_top_sub_cell;
        left_front_top_sub_cell->next = left_front_down_sub_cell;
        left_front_down_sub_cell->next = right_front_down_sub_cell;
        right_front_down_sub_cell->next = right_back_down_sub_cell;
        right_back_down_sub_cell->next = left_back_down_sub_cell;
        left_back_down_sub_cell->next = left_back_top_sub_cell;
        left_back_top_sub_cell->next = right_back_top_sub_cell;

        right_back_top_sub_cell->previous = left_back_top_sub_cell;
        left_back_top_sub_cell->previous = left_back_down_sub_cell;
        left_back_down_sub_cell->previous = right_back_down_sub_cell;
        right_back_down_sub_cell->previous = right_front_down_sub_cell;
        right_front_down_sub_cell->previous = left_front_down_sub_cell;
        left_front_down_sub_cell->previous = left_front_top_sub_cell;
        left_front_top_sub_cell->previous = right_front_top_sub_cell;

        if(right_back_top_sub_cell->next != 0)
            right_back_top_sub_cell->next->previous = right_back_top_sub_cell;

    } else if(number_of_hilbert_shape == 1) {
        /* Shape 1
                                       e
                               /|      |      b: begin
                             /  |   b  |      e: end
                             |  |___|__|
                             |      |
                             |______|
         */

        right_front_top_sub_cell->hilbert_shape_number = 0;
        right_back_top_sub_cell->hilbert_shape_number = 2;
        right_back_down_sub_cell->hilbert_shape_number = 2;
        right_front_down_sub_cell->hilbert_shape_number = 6;
        left_front_down_sub_cell->hilbert_shape_number = 6;
        left_back_down_sub_cell->hilbert_shape_number = 7;
        left_back_top_sub_cell->hilbert_shape_number = 7;
        left_front_top_sub_cell->hilbert_shape_number = 8;

        left_front_top_sub_cell->next = right_front_top_sub_cell->next;
        right_front_top_sub_cell->next = right_back_top_sub_cell;
        right_back_top_sub_cell->next = right_back_down_sub_cell;
        right_back_down_sub_cell->next = right_front_down_sub_cell;
        right_front_down_sub_cell->next = left_front_down_sub_cell;
        left_front_down_sub_cell->next = left_back_down_sub_cell;
        left_back_down_sub_cell->next = left_back_top_sub_cell;
        left_back_top_sub_cell->next = left_front_top_sub_cell;

        left_front_top_sub_cell->previous = left_back_top_sub_cell;
        left_back_top_sub_cell->previous = left_back_down_sub_cell;
        left_back_down_sub_cell->previous = left_front_down_sub_cell;
        left_front_down_sub_cell->previous = right_front_down_sub_cell;
        right_front_down_sub_cell->previous = right_back_down_sub_cell;
        right_back_down_sub_cell->previous = right_back_top_sub_cell;
        right_back_top_sub_cell->previous = right_front_top_sub_cell;

        if(left_front_top_sub_cell->next != 0)
            left_front_top_sub_cell->next->previous = left_front_top_sub_cell;

    } else if(number_of_hilbert_shape == 2) {
        /* Shape 2
                               /|     /|      b: begin
                             e/ |   b/ |      e: end
                                |      |
                               /      /
                              /______/
         */

        right_front_top_sub_cell->hilbert_shape_number = 1;
        left_front_top_sub_cell->hilbert_shape_number = 0;
        left_back_top_sub_cell->hilbert_shape_number = 0;
        right_back_top_sub_cell->hilbert_shape_number = 9;
        right_back_down_sub_cell->hilbert_shape_number = 9;
        left_back_down_sub_cell->hilbert_shape_number = 10;
        left_front_down_sub_cell->hilbert_shape_number = 10;
        right_front_down_sub_cell->hilbert_shape_number = 11;

        right_front_down_sub_cell->next = right_front_top_sub_cell->next;
        right_front_top_sub_cell->next = left_front_top_sub_cell;
        left_front_top_sub_cell->next = left_back_top_sub_cell;
        left_back_top_sub_cell->next = right_back_top_sub_cell;
        right_back_top_sub_cell->next = right_back_down_sub_cell;
        right_back_down_sub_cell->next = left_back_down_sub_cell;
        left_back_down_sub_cell->next = left_front_down_sub_cell;
        left_front_down_sub_cell->next = right_front_down_sub_cell;

        right_front_down_sub_cell->previous = left_front_down_sub_cell;
        left_front_down_sub_cell->previous = left_back_down_sub_cell;
        left_back_down_sub_cell->previous = right_back_down_sub_cell;
        right_back_down_sub_cell->previous = right_back_top_sub_cell;
        right_back_top_sub_cell->previous = left_back_top_sub_cell;
        left_back_top_sub_cell->previous = left_front_top_sub_cell;
        left_front_top_sub_cell->previous = right_front_top_sub_cell;

        if(right_front_down_sub_cell->next != 0)
            right_front_down_sub_cell->next->previous = right_front_down_sub_cell;

    } else if(number_of_hilbert_shape == 3) {
        /* Shape 3
                               /b     /|      b: begin
                              /______/ |      e: end
                                       |
                               /e     /
                              /______/
         */

        left_front_down_sub_cell->hilbert_shape_number = 11;
        right_front_down_sub_cell->hilbert_shape_number = 7;
        right_front_top_sub_cell->hilbert_shape_number = 7;
        left_front_top_sub_cell->hilbert_shape_number = 0;
        left_back_top_sub_cell->hilbert_shape_number = 0;
        right_back_top_sub_cell->hilbert_shape_number = 9;
        right_back_down_sub_cell->hilbert_shape_number = 9;
        left_back_down_sub_cell->hilbert_shape_number = 6;

        if(right_front_top_sub_cell->previous != 0)
            right_front_top_sub_cell->previous->next = left_front_down_sub_cell;

        left_back_down_sub_cell->next = right_front_top_sub_cell->next;
        left_front_down_sub_cell->next = right_front_down_sub_cell;
        right_front_down_sub_cell->next = right_front_top_sub_cell;
        right_front_top_sub_cell->next = left_front_top_sub_cell;
        left_front_top_sub_cell->next = left_back_top_sub_cell;
        left_back_top_sub_cell->next = right_back_top_sub_cell;
        right_back_top_sub_cell->next = right_back_down_sub_cell;
        right_back_down_sub_cell->next = left_back_down_sub_cell;

        left_front_down_sub_cell->previous = right_front_top_sub_cell->previous;
        left_back_down_sub_cell->previous = right_back_down_sub_cell;
        right_back_down_sub_cell->previous = right_back_top_sub_cell;
        right_back_top_sub_cell->previous = left_back_top_sub_cell;
        left_back_top_sub_cell->previous = left_front_top_sub_cell;
        left_front_top_sub_cell->previous = right_front_top_sub_cell;
        right_front_top_sub_cell->previous = right_front_down_sub_cell;
        right_front_down_sub_cell->previous = left_front_down_sub_cell;

        if(left_back_down_sub_cell->next != 0)
            left_back_down_sub_cell->next->previous = left_back_down_sub_cell;

    } else if(number_of_hilbert_shape == 4) {
        /* Shape 4
                                /|     /|      b: begin
                               /_|____/ |      e: end
                                 |      |
                                /      /
                              b/     e/
         */

        right_back_down_sub_cell->hilbert_shape_number = 6;
        left_back_down_sub_cell->hilbert_shape_number = 10;
        left_front_down_sub_cell->hilbert_shape_number = 10;
        right_front_down_sub_cell->hilbert_shape_number = 7;
        right_front_top_sub_cell->hilbert_shape_number = 7;
        left_front_top_sub_cell->hilbert_shape_number = 0;
        left_back_top_sub_cell->hilbert_shape_number = 0;
        right_back_top_sub_cell->hilbert_shape_number = 5;

        if(right_front_top_sub_cell->previous != 0)
            right_front_top_sub_cell->previous->next = right_back_down_sub_cell;

        right_back_top_sub_cell->next = right_front_top_sub_cell->next;
        right_back_down_sub_cell->next = left_back_down_sub_cell;
        left_back_down_sub_cell->next = left_front_down_sub_cell;
        left_front_down_sub_cell->next = right_front_down_sub_cell;
        right_front_down_sub_cell->next = right_front_top_sub_cell;
        right_front_top_sub_cell->next = left_front_top_sub_cell;
        left_front_top_sub_cell->next = left_back_top_sub_cell;
        left_back_top_sub_cell->next = right_back_top_sub_cell;

        right_back_down_sub_cell->previous = right_front_top_sub_cell->previous;
        right_back_top_sub_cell->previous = left_back_top_sub_cell;
        left_back_top_sub_cell->previous = left_front_top_sub_cell;
        left_front_top_sub_cell->previous = right_front_top_sub_cell;
        right_front_top_sub_cell->previous = right_front_down_sub_cell;
        right_front_down_sub_cell->previous = left_front_down_sub_cell;
        left_front_down_sub_cell->previous = left_back_down_sub_cell;
        left_back_down_sub_cell->previous = right_back_down_sub_cell;

        if(right_back_top_sub_cell->next != 0)
            right_back_top_sub_cell->next->previous = right_back_top_sub_cell;

    } else if(number_of_hilbert_shape == 5) {
        /* Shape 5
                                 ______
                                |      |      b: begin
                              __|___   |      e: end
                             |  |   |  |b
                             | /    |
                             |/     |e
         */

        left_back_top_sub_cell->hilbert_shape_number = 8;
        left_front_top_sub_cell->hilbert_shape_number = 9;
        left_front_down_sub_cell->hilbert_shape_number = 9;
        left_back_down_sub_cell->hilbert_shape_number = 11;
        right_back_down_sub_cell->hilbert_shape_number = 11;
        right_front_down_sub_cell->hilbert_shape_number = 4;
        right_front_top_sub_cell->hilbert_shape_number = 4;
        right_back_top_sub_cell->hilbert_shape_number = 0;

        if(right_front_top_sub_cell->previous != 0)
            right_front_top_sub_cell->previous->next = left_back_top_sub_cell;

        right_back_top_sub_cell->next = right_front_top_sub_cell->next;
        left_back_top_sub_cell->next = left_front_top_sub_cell;
        left_front_top_sub_cell->next = left_front_down_sub_cell;
        left_front_down_sub_cell->next = left_back_down_sub_cell;
        left_back_down_sub_cell->next = right_back_down_sub_cell;
        right_back_down_sub_cell->next = right_front_down_sub_cell;
        right_front_down_sub_cell->next = right_front_top_sub_cell;
        right_front_top_sub_cell->next = right_back_top_sub_cell;

        left_back_top_sub_cell->previous = right_front_top_sub_cell->previous;
        right_back_top_sub_cell->previous = right_front_top_sub_cell;
        right_front_top_sub_cell->previous = right_front_down_sub_cell;
        right_front_down_sub_cell->previous = right_back_down_sub_cell;
        right_back_down_sub_cell->previous = left_back_down_sub_cell;
        left_back_down_sub_cell->previous = left_front_down_sub_cell;
        left_front_down_sub_cell->previous = left_front_top_sub_cell;
        left_front_top_sub_cell->previous = left_back_top_sub_cell;

        if(right_back_top_sub_cell->next != 0)
            right_back_top_sub_cell->next->previous = right_back_top_sub_cell;

    } else if(number_of_hilbert_shape == 6) {
        /* Shape 6
                                 ______
                                |      |      b: begin
                              __|___   |      e: end
                             | e|   |  |
                             |      | /
                            b|      |/
         */

        right_back_down_sub_cell->hilbert_shape_number = 10;
        right_front_down_sub_cell->hilbert_shape_number = 4;
        right_front_top_sub_cell->hilbert_shape_number = 4;
        right_back_top_sub_cell->hilbert_shape_number = 1;
        left_back_top_sub_cell->hilbert_shape_number = 1;
        left_front_top_sub_cell->hilbert_shape_number = 9;
        left_front_down_sub_cell->hilbert_shape_number = 9;
        left_back_down_sub_cell->hilbert_shape_number = 3;

        if(right_front_top_sub_cell->previous != 0)
            right_front_top_sub_cell->previous->next = right_back_down_sub_cell;

        left_back_down_sub_cell->next = right_front_top_sub_cell->next;
        right_back_down_sub_cell->next = right_front_down_sub_cell;
        right_front_down_sub_cell->next = right_front_top_sub_cell;
        right_front_top_sub_cell->next = right_back_top_sub_cell;
        right_back_top_sub_cell->next = left_back_top_sub_cell;
        left_back_top_sub_cell->next = left_front_top_sub_cell;
        left_front_top_sub_cell->next = left_front_down_sub_cell;
        left_front_down_sub_cell->next = left_back_down_sub_cell;

        right_back_down_sub_cell->previous = right_front_top_sub_cell->previous;
        left_back_down_sub_cell->previous = left_front_down_sub_cell;
        left_front_down_sub_cell->previous = left_front_top_sub_cell;
        left_front_top_sub_cell->previous = left_back_top_sub_cell;
        left_back_top_sub_cell->previous = right_back_top_sub_cell;
        right_back_top_sub_cell->previous = right_front_top_sub_cell;
        right_front_top_sub_cell->previous = right_front_down_sub_cell;
        right_front_down_sub_cell->previous = right_back_down_sub_cell;

        if(left_back_down_sub_cell->next != 0)
            left_back_down_sub_cell->next->previous = left_back_down_sub_cell;

    } else if(number_of_hilbert_shape == 7) {
        /* Shape 7
                                |b     |e     b: begin
                              __|___   |      e: end
                             |  |   |  |
                             | /    | /
                             |/     |/
         */

        left_front_down_sub_cell->hilbert_shape_number = 3;
        left_back_down_sub_cell->hilbert_shape_number = 11;
        right_back_down_sub_cell->hilbert_shape_number = 11;
        right_front_down_sub_cell->hilbert_shape_number = 4;
        right_front_top_sub_cell->hilbert_shape_number = 4;
        right_back_top_sub_cell->hilbert_shape_number = 1;
        left_back_top_sub_cell->hilbert_shape_number = 1;
        left_front_top_sub_cell->hilbert_shape_number = 8;

        if(right_front_top_sub_cell->previous != 0)
            right_front_top_sub_cell->previous->next = left_front_down_sub_cell;

        left_front_top_sub_cell->next = right_front_top_sub_cell->next;
        left_front_down_sub_cell->next = left_back_down_sub_cell;
        left_back_down_sub_cell->next = right_back_down_sub_cell;
        right_back_down_sub_cell->next = right_front_down_sub_cell;
        right_front_down_sub_cell->next = right_front_top_sub_cell;
        right_front_top_sub_cell->next = right_back_top_sub_cell;
        right_back_top_sub_cell->next = left_back_top_sub_cell;
        left_back_top_sub_cell->next = left_front_top_sub_cell;

        left_front_down_sub_cell->previous = right_front_top_sub_cell->previous;
        left_front_top_sub_cell->previous = left_back_top_sub_cell;
        left_back_top_sub_cell->previous = right_back_top_sub_cell;
        right_back_top_sub_cell->previous = right_front_top_sub_cell;
        right_front_top_sub_cell->previous = right_front_down_sub_cell;
        right_front_down_sub_cell->previous = right_back_down_sub_cell;
        right_back_down_sub_cell->previous = left_back_down_sub_cell;
        left_back_down_sub_cell->previous = left_front_down_sub_cell;

        if(left_front_top_sub_cell->next != 0)
            left_front_top_sub_cell->next->previous = left_front_top_sub_cell;

    } else if(number_of_hilbert_shape == 8) {
        /* Shape 8
                               /|     /e      b: begin
                              /_|____/        e: end
                                |
                               /      /b
                              /______/
         */

        left_back_top_sub_cell->hilbert_shape_number = 5;
        right_back_top_sub_cell->hilbert_shape_number = 9;
        right_back_down_sub_cell->hilbert_shape_number = 9;
        left_back_down_sub_cell->hilbert_shape_number = 10;
        left_front_down_sub_cell->hilbert_shape_number = 10;
        right_front_down_sub_cell->hilbert_shape_number = 7;
        right_front_top_sub_cell->hilbert_shape_number = 7;
        left_front_top_sub_cell->hilbert_shape_number = 1;

        if(right_front_top_sub_cell->previous != 0)
            right_front_top_sub_cell->previous->next = left_back_top_sub_cell;

        left_front_top_sub_cell->next = right_front_top_sub_cell->next;
        left_back_top_sub_cell->next = right_back_top_sub_cell;
        right_back_top_sub_cell->next = right_back_down_sub_cell;
        right_back_down_sub_cell->next = left_back_down_sub_cell;
        left_back_down_sub_cell->next = left_front_down_sub_cell;
        left_front_down_sub_cell->next = right_front_down_sub_cell;
        right_front_down_sub_cell->next = right_front_top_sub_cell;
        right_front_top_sub_cell->next = left_front_top_sub_cell;

        left_back_top_sub_cell->previous = right_front_top_sub_cell->previous;
        left_front_top_sub_cell->previous = right_front_top_sub_cell;
        right_front_top_sub_cell->previous = right_front_down_sub_cell;
        right_front_down_sub_cell->previous = left_front_down_sub_cell;
        left_front_down_sub_cell->previous = left_back_down_sub_cell;
        left_back_down_sub_cell->previous = right_back_down_sub_cell;
        right_back_down_sub_cell->previous = right_back_top_sub_cell;
        right_back_top_sub_cell->previous = left_back_top_sub_cell;

        if(left_front_top_sub_cell->next != 0)
            left_front_top_sub_cell->next->previous = left_front_top_sub_cell;

    } else if(number_of_hilbert_shape == 9) {
        /* Shape 9
                                _______
                               /      /      b: begin
                              /      /       e: end
                              |      |
                              | /e   | /b
                              |/     |/
         */

        left_back_top_sub_cell->hilbert_shape_number = 5;
        right_back_top_sub_cell->hilbert_shape_number = 8;
        right_front_top_sub_cell->hilbert_shape_number = 8;
        left_front_top_sub_cell->hilbert_shape_number = 2;
        left_front_down_sub_cell->hilbert_shape_number = 2;
        right_front_down_sub_cell->hilbert_shape_number = 3;
        right_back_down_sub_cell->hilbert_shape_number = 3;
        left_back_down_sub_cell->hilbert_shape_number = 6;

        if(right_front_top_sub_cell->previous != 0)
            right_front_top_sub_cell->previous->next = left_back_top_sub_cell;

        left_back_down_sub_cell->next = right_front_top_sub_cell->next;
        left_back_top_sub_cell->next = right_back_top_sub_cell;
        right_back_top_sub_cell->next = right_front_top_sub_cell;
        right_front_top_sub_cell->next = left_front_top_sub_cell;
        left_front_top_sub_cell->next = left_front_down_sub_cell;
        left_front_down_sub_cell->next = right_front_down_sub_cell;
        right_front_down_sub_cell->next = right_back_down_sub_cell;
        right_back_down_sub_cell->next = left_back_down_sub_cell;

        left_back_top_sub_cell->previous = right_front_top_sub_cell->previous;
        left_back_down_sub_cell->previous = right_back_down_sub_cell;
        right_back_down_sub_cell->previous = right_front_down_sub_cell;
        right_front_down_sub_cell->previous = left_front_down_sub_cell;
        left_front_down_sub_cell->previous = left_front_top_sub_cell;
        left_front_top_sub_cell->previous = right_front_top_sub_cell;
        right_front_top_sub_cell->previous = right_back_top_sub_cell;
        right_back_top_sub_cell->previous = left_back_top_sub_cell;

        if(left_back_down_sub_cell->next != 0)
            left_back_down_sub_cell->next->previous = left_back_down_sub_cell;

    } else if(number_of_hilbert_shape == 10) {
        /* Shape 10
                                 _______
                                /      /      b: begin
                              e/      /       e: end
                                  ____|__
                                 /    | /
                               b/     |/
         */

        right_back_down_sub_cell->hilbert_shape_number = 6;
        left_back_down_sub_cell->hilbert_shape_number = 4;
        left_back_top_sub_cell->hilbert_shape_number = 4;
        right_back_top_sub_cell->hilbert_shape_number = 8;
        right_front_top_sub_cell->hilbert_shape_number = 8;
        left_front_top_sub_cell->hilbert_shape_number = 2;
        left_front_down_sub_cell->hilbert_shape_number = 2;
        right_front_down_sub_cell->hilbert_shape_number = 11;

        if(right_front_top_sub_cell->previous != 0)
            right_front_top_sub_cell->previous->next = right_back_down_sub_cell;

        right_front_down_sub_cell->next = right_front_top_sub_cell->next;
        right_back_down_sub_cell->next = left_back_down_sub_cell;
        left_back_down_sub_cell->next = left_back_top_sub_cell;
        left_back_top_sub_cell->next = right_back_top_sub_cell;
        right_back_top_sub_cell->next = right_front_top_sub_cell;
        right_front_top_sub_cell->next = left_front_top_sub_cell;
        left_front_top_sub_cell->next = left_front_down_sub_cell;
        left_front_down_sub_cell->next = right_front_down_sub_cell;

        right_back_down_sub_cell->previous = right_front_top_sub_cell->previous;
        right_front_down_sub_cell->previous = left_front_down_sub_cell;
        left_front_down_sub_cell->previous = left_front_top_sub_cell;
        left_front_top_sub_cell->previous = right_front_top_sub_cell;
        right_front_top_sub_cell->previous = right_back_top_sub_cell;
        right_back_top_sub_cell->previous = left_back_top_sub_cell;
        left_back_top_sub_cell->previous = left_back_down_sub_cell;
        left_back_down_sub_cell->previous = right_back_down_sub_cell;

        if(right_front_down_sub_cell->next != 0)
            right_front_down_sub_cell->next->previous = right_front_down_sub_cell;

    } else if(number_of_hilbert_shape == 11) {
        /* Shape 11
                                     b______
                                            |      b: begin
                                   e______  |      e: end
                                     _____|_|
                                    /     |
                                   /______|
         */

        left_front_down_sub_cell->hilbert_shape_number = 7;
        left_front_top_sub_cell->hilbert_shape_number = 3;
        left_back_top_sub_cell->hilbert_shape_number = 3;
        left_back_down_sub_cell->hilbert_shape_number = 5;
        right_back_down_sub_cell->hilbert_shape_number = 5;
        right_back_top_sub_cell->hilbert_shape_number = 10;
        right_front_top_sub_cell->hilbert_shape_number = 10;
        right_front_down_sub_cell->hilbert_shape_number = 2;

        if(right_front_top_sub_cell->previous != 0)
            right_front_top_sub_cell->previous->next = left_front_down_sub_cell;

        right_front_down_sub_cell->next = right_front_top_sub_cell->next;
        left_front_down_sub_cell->next = left_front_top_sub_cell;
        left_front_top_sub_cell->next = left_back_top_sub_cell;
        left_back_top_sub_cell->next = left_back_down_sub_cell;
        left_back_down_sub_cell->next = right_back_down_sub_cell;
        right_back_down_sub_cell->next = right_back_top_sub_cell;
        right_back_top_sub_cell->next = right_front_top_sub_cell;
        right_front_top_sub_cell->next = right_front_down_sub_cell;

        left_front_down_sub_cell->previous = right_front_top_sub_cell->previous;
        right_front_down_sub_cell->previous = right_front_top_sub_cell;
        right_front_top_sub_cell->previous = right_back_top_sub_cell;
        right_back_top_sub_cell->previous = right_back_down_sub_cell;
        right_back_down_sub_cell->previous = left_back_down_sub_cell;
        left_back_down_sub_cell->previous = left_back_top_sub_cell;
        left_back_top_sub_cell->previous = left_front_top_sub_cell;
        left_front_top_sub_cell->previous = left_front_down_sub_cell;

        if(right_front_down_sub_cell->next != 0)
            right_front_down_sub_cell->next->previous = right_front_down_sub_cell;
    }

    // If necessary, simplifies the graph by eliminating adjacent transition nodes
    // of same level connected through their single connectors.
    simplify_refinement(top_transition_node);
    simplify_refinement(front_transition_node);
    simplify_refinement(down_transition_node);
    simplify_refinement(back_transition_node);
    simplify_refinement(right_transition_node);
    simplify_refinement(left_transition_node);
}

/**
 * Simplifies data structure eliminating adjacent transition nodes of same level.
 *
 * @param transition_node Candidate transition node to be eliminated.
 *
 */
void simplify_refinement(struct transition_node *transition_node) {

    assert(transition_node);

    if(transition_node->single_connector != 0) {

        // Both transition node and neighbor transition node must have the same
        // refinement level.
        enum cell_type node_type = transition_node->cell_data.type;
        uint16_t node_level = transition_node->cell_data.level;

        uint16_t single_connector_level = ((struct basic_cell_data *)(transition_node->single_connector))->level;

        if((node_type == TRANSITION_NODE) && (node_level == single_connector_level)) {
            struct transition_node *neighbour_node = (struct transition_node *)(transition_node->single_connector);

            struct cell_node *cell_node[4];
            cell_node[0] = (struct cell_node *)(transition_node->quadruple_connector1);
            cell_node[1] = (struct cell_node *)(transition_node->quadruple_connector2);
            cell_node[2] = (struct cell_node *)(transition_node->quadruple_connector3);
            cell_node[3] = (struct cell_node *)(transition_node->quadruple_connector4);

            struct cell_node *neighbor_cell[4];
            neighbor_cell[0] = (struct cell_node *)neighbour_node->quadruple_connector1;
            neighbor_cell[1] = (struct cell_node *)neighbour_node->quadruple_connector2;
            neighbor_cell[2] = (struct cell_node *)neighbour_node->quadruple_connector3;
            neighbor_cell[3] = (struct cell_node *)neighbour_node->quadruple_connector4;

            enum transition_direction direction = transition_node->direction;
            enum cell_type type;

            for(int i = 0; i < 4; i++) {
                switch(direction) {
                case FRONT: {
                    cell_node[i]->neighbours[FRONT] = neighbor_cell[i];
                    break;
                }
                case BACK: {
                    cell_node[i]->neighbours[BACK] = neighbor_cell[i];
                    break;
                }
                case TOP: {
                    cell_node[i]->neighbours[TOP] = neighbor_cell[i];
                    break;
                }
                case DOWN: {
                    cell_node[i]->neighbours[DOWN] = neighbor_cell[i];
                    break;
                }
                case RIGHT: {
                    cell_node[i]->neighbours[RIGHT] = neighbor_cell[i];
                    break;
                }
                case LEFT: {
                    cell_node[i]->neighbours[LEFT] = neighbor_cell[i];
                    break;
                }
                default:
                    break;
                }

                type = neighbor_cell[i]->cell_data.type;
                switch(type) {
                case CELL_NODE: {
                    struct cell_node *neighbour_cell_node = neighbor_cell[i];
                    switch(direction) {
                    case FRONT: {
                        neighbour_cell_node->neighbours[BACK] = cell_node[i];
                        break;
                    }
                    case BACK: {
                        neighbour_cell_node->neighbours[FRONT] = cell_node[i];
                        break;
                    }
                    case TOP: {
                        neighbour_cell_node->neighbours[DOWN] = cell_node[i];
                        break;
                    }
                    case DOWN: {
                        neighbour_cell_node->neighbours[TOP] = cell_node[i];
                        break;
                    }
                    case RIGHT: {
                        neighbour_cell_node->neighbours[LEFT] = cell_node[i];
                        break;
                    }
                    case LEFT: {
                        neighbour_cell_node->neighbours[RIGHT] = cell_node[i];
                        break;
                    }
                    default:
                        break;
                    }
                    break;
                }

                case TRANSITION_NODE: {
                    struct transition_node *neighbour_transition_node = (struct transition_node *)(neighbor_cell[i]);
                    if(neighbour_node == neighbour_transition_node->single_connector)
                        neighbour_transition_node->single_connector = cell_node[i];

                    else if(neighbour_node == neighbour_transition_node->quadruple_connector1)
                        neighbour_transition_node->quadruple_connector1 = cell_node[i];

                    else if(neighbour_node == neighbour_transition_node->quadruple_connector2)
                        neighbour_transition_node->quadruple_connector2 = cell_node[i];

                    else if(neighbour_node == neighbour_transition_node->quadruple_connector3)
                        neighbour_transition_node->quadruple_connector3 = cell_node[i];

                    else if(neighbour_node == neighbour_transition_node->quadruple_connector4)
                        neighbour_transition_node->quadruple_connector4 = cell_node[i];

                    break;
                }

                default:
                    break;
                }
            }
            free(transition_node);
            free(neighbour_node);
        }
    }
}

void set_refined_cell_data(struct cell_node *the_cell, struct cell_node *other_cell, struct point_3d discretization, struct point_3d center,
                           uint64_t bunch_number, ui32_array free_sv_positions, ui32_array *refined_this_step) {

    the_cell->cell_data.level = other_cell->cell_data.level;
    the_cell->active = other_cell->active;

    if(other_cell->mesh_extra_info) {
        the_cell->mesh_extra_info_size = other_cell->mesh_extra_info_size;
        the_cell->mesh_extra_info = malloc(the_cell->mesh_extra_info_size);
        memcpy(the_cell->mesh_extra_info, other_cell->mesh_extra_info, the_cell->mesh_extra_info_size);
    }

    if(other_cell->linear_system_solver_extra_info) {
        the_cell->linear_system_solver_extra_info_size = other_cell->linear_system_solver_extra_info_size;
        the_cell->linear_system_solver_extra_info = malloc(the_cell->linear_system_solver_extra_info_size);
        memcpy(the_cell->linear_system_solver_extra_info, other_cell->linear_system_solver_extra_info, the_cell->linear_system_solver_extra_info_size);
    }

    the_cell->v = other_cell->v;
    the_cell->sigma = other_cell->sigma;
    the_cell->discretization = discretization;
    the_cell->center = center;
    the_cell->original_position_in_file = other_cell->original_position_in_file;

    the_cell->bunch_number = bunch_number;

    if(free_sv_positions) {
        the_cell->sv_position = arrpop(free_sv_positions);
    }

    if(refined_this_step && *refined_this_step) {
        arrput(*refined_this_step, the_cell->sv_position);
    }
}

void set_refined_transition_node_data(struct transition_node *the_node, struct cell_node *other_node, enum transition_direction direction) {

    if(!VALID_SIMPLE_DIRECTION(direction)) {
        fprintf(stderr, "set_refined_transition_node_data() invalid direction %d. Exiting!", direction);
        exit(10);
    }

    the_node->direction = direction;
    the_node->cell_data.level = other_node->cell_data.level;
    the_node->single_connector = other_node->neighbours[direction];
}
