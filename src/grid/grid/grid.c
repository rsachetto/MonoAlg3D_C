//
// Created by sachetto on 29/09/17.
//

#include "grid.h"
#include "../../solvers/constants.h"
#include <inttypes.h>
void initialize_and_construct_grid (struct grid *the_grid, double side_length, uint8_t num_cell_neighbours) {
    initialize_grid (the_grid, side_length, num_cell_neighbours);
    construct_grid (the_grid);
}

void initialize_grid (struct grid *the_grid, double side_length, uint8_t num_cell_neighbours) {

    the_grid->first_cell = NULL;
    the_grid->side_length = side_length;
    the_grid->number_of_cells = 0;
    the_grid->init_ode = false;
    the_grid->init_ode = false;
    the_grid->active_cells = NULL;
    the_grid->num_cell_neighbours = num_cell_neighbours;
    the_grid->refined_this_step = uint32_vector_create(128);
    the_grid->free_sv_positions = uint32_vector_create(128);

}

void construct_grid (struct grid *the_grid) {

    bool init_ode = the_grid->init_ode;
    double side_length = the_grid->side_length;

    size_t cell_node_size = sizeof (struct cell_node);
    size_t transition_node_size = sizeof (struct transition_node);

    // Cell nodes.
    struct cell_node *front_northeast_cell, *front_northwest_cell, *front_southeast_cell, *front_southwest_cell,
        *back_northeast_cell, *back_northwest_cell, *back_southeast_cell, *back_southwest_cell;

    front_northeast_cell = new_cell_node ();
    front_northwest_cell = new_cell_node ();
    front_southeast_cell = new_cell_node ();
    front_southwest_cell = new_cell_node ();
    back_northeast_cell = new_cell_node ();
    back_northwest_cell = new_cell_node ();
    back_southeast_cell = new_cell_node ();
    back_southwest_cell = new_cell_node ();

    // Transition nodes.
    struct transition_node *north_transition_node, *south_transition_node, *east_transition_node, *west_transition_node,
        *front_transition_node, *back_transition_node;

    north_transition_node = new_transition_node ();
    south_transition_node = new_transition_node ();
    east_transition_node = new_transition_node ();
    west_transition_node = new_transition_node ();
    front_transition_node = new_transition_node ();
    back_transition_node = new_transition_node ();

    double half_side_length = side_length / 2.0f;
    double quarter_side_length = side_length / 4.0f;
    //__________________________________________________________________________
    //              Initialization of transition nodes.
    //__________________________________________________________________________
    // East transition node.
    set_transition_node_data (east_transition_node, 1, 'e', NULL, front_southeast_cell, back_southeast_cell,
                              back_northeast_cell, front_northeast_cell);

    // North transition node.
    set_transition_node_data (north_transition_node, 1, 'n', NULL, front_northwest_cell, front_northeast_cell,
                              back_northeast_cell, back_northwest_cell);

    // West transition node.
    set_transition_node_data (west_transition_node, 1, 'w', NULL, front_southwest_cell, back_southwest_cell,
                              back_northwest_cell, front_northwest_cell);

    // South transition node.
    set_transition_node_data (south_transition_node, 1, 's', NULL, front_southwest_cell, front_southeast_cell,
                              back_southeast_cell, back_southwest_cell);

    // Front transition node.
    set_transition_node_data (front_transition_node, 1, 'f', NULL, front_southwest_cell, front_southeast_cell,
                              front_northeast_cell, front_northwest_cell);

    // Back transition node.
    set_transition_node_data (back_transition_node, 1, 'b', NULL, back_southwest_cell, back_southeast_cell,
                              back_northeast_cell, back_northwest_cell);

    //__________________________________________________________________________
    //                      Initialization of cell nodes.
    //__________________________________________________________________________
    // front Northeast subcell initialization.
    set_cell_node_data (front_northeast_cell, half_side_length, quarter_side_length, 1, east_transition_node,
                        north_transition_node, front_northwest_cell, front_southeast_cell, front_transition_node,
                        back_northeast_cell, NULL, back_northeast_cell, 0, 1, half_side_length + quarter_side_length,
                        half_side_length + quarter_side_length, half_side_length + quarter_side_length);

    // back Northeast subcell initialization.
    set_cell_node_data (back_northeast_cell, half_side_length, quarter_side_length, 2, east_transition_node,
                        north_transition_node, back_northwest_cell, back_southeast_cell, front_northeast_cell,
                        back_transition_node, front_northeast_cell, back_northwest_cell, 1, 2, quarter_side_length,
                        half_side_length + quarter_side_length, half_side_length + quarter_side_length);

    // back Northwest subcell initialization.
    set_cell_node_data (back_northwest_cell, half_side_length, quarter_side_length, 3, back_northeast_cell,
                        north_transition_node, west_transition_node, back_southwest_cell, front_northwest_cell,
                        back_transition_node, back_northeast_cell, front_northwest_cell, 2, 2, quarter_side_length,
                        quarter_side_length, half_side_length + quarter_side_length);

    // front Northwest subcell initialization.
    set_cell_node_data (front_northwest_cell, half_side_length, quarter_side_length, 4, front_northeast_cell,
                        north_transition_node, west_transition_node, front_southwest_cell, front_transition_node,
                        back_northwest_cell, back_northwest_cell, front_southwest_cell, 3, 3,
                        half_side_length + quarter_side_length, quarter_side_length,
                        half_side_length + quarter_side_length);

    // front Southwest subcell initialization.
    set_cell_node_data (front_southwest_cell, half_side_length, quarter_side_length, 5, front_southeast_cell,
                        front_northwest_cell, west_transition_node, south_transition_node, front_transition_node,
                        back_southwest_cell, front_northwest_cell, back_southwest_cell, 4, 3,
                        half_side_length + quarter_side_length, quarter_side_length, quarter_side_length);

    // back Southwest subcell initialization.
    set_cell_node_data (back_southwest_cell, half_side_length, quarter_side_length, 6, back_southeast_cell,
                        back_northwest_cell, west_transition_node, south_transition_node, front_southwest_cell,
                        back_transition_node, front_southwest_cell, back_southeast_cell, 5, 4, quarter_side_length,
                        quarter_side_length, quarter_side_length);

    // back Southeast subcell initialization.
    set_cell_node_data (back_southeast_cell, half_side_length, quarter_side_length, 7, east_transition_node,
                        back_northeast_cell, back_southwest_cell, south_transition_node, front_southeast_cell,
                        back_transition_node, back_southwest_cell, front_southeast_cell, 6, 4, quarter_side_length,
                        half_side_length + quarter_side_length, quarter_side_length);

    // front Southeast subcell initialization.
    set_cell_node_data (front_southeast_cell, half_side_length, quarter_side_length, 8, east_transition_node,
                        front_northeast_cell, front_southwest_cell, south_transition_node, front_transition_node,
                        back_southeast_cell, back_southeast_cell, NULL, 7, 5, half_side_length + quarter_side_length,
                        half_side_length + quarter_side_length, quarter_side_length);

    // Grid initialization
    the_grid->first_cell = front_northeast_cell;
    the_grid->number_of_cells = 8;
}

void print_grid (struct grid *the_grid, FILE *output_file) {

    struct cell_node *grid_cell = the_grid->first_cell;

    double center_x, center_y, center_z, half_face;
    double v;

    while (grid_cell != 0) {

        if (grid_cell->active) {

            center_x = grid_cell->center_x;
            center_y = grid_cell->center_y;
            center_z = grid_cell->center_z;

            v = grid_cell->v;
            half_face = grid_cell->half_face_length;

            fprintf (output_file, "%lf,%lf,%lf,%lf,%lf\n", center_x, center_y, center_z, half_face, v);
        }
        grid_cell = grid_cell->next;
    }
}

bool print_grid_and_check_for_activity (struct grid *the_grid, FILE *output_file, int count) {

    struct cell_node *grid_cell = the_grid->first_cell;

    double center_x, center_y, center_z, half_face;
    double v;
    bool act = false;

    while (grid_cell != 0) {

        if (grid_cell->active) {

            center_x = grid_cell->center_x;
            center_y = grid_cell->center_y;
            center_z = grid_cell->center_z;

            v = grid_cell->v;
            half_face = grid_cell->half_face_length;

            if (count > 0) {
                if (grid_cell->v > -86.0) {
                    act = true;
                }
            } else {
                act = true;
            }

            fprintf (output_file, "%lf,%lf,%lf,%lf,%.4lf\n", center_x, center_y, center_z, half_face, v);
        }
        grid_cell = grid_cell->next;
    }

    return act;
}

void order_grid_cells (struct grid *the_grid) {

    struct cell_node *grid_cell;
    grid_cell = the_grid->first_cell;

    if (the_grid->active_cells != NULL) {
        free (the_grid->active_cells);
    }

    the_grid->active_cells = (struct cell_node **)malloc (sizeof (struct cell_node *) * the_grid->number_of_cells);

    uint64_t counter = 0;
    while (grid_cell != 0) {
        if (grid_cell->active) {

            grid_cell->grid_position = counter;
            the_grid->active_cells[counter] = grid_cell;
            counter++;
        }

        grid_cell = grid_cell->next;
    }

    the_grid->num_active_cells = counter;
}

void clean_grid (struct grid *the_grid) {


    uint64_t number_of_cells = the_grid->number_of_cells;

    // In order to release the memory allocated for the grid, the grid is
    // derefined to level 1. Thus, the grid shape is known and each node can
    // be easily reached.
    while (number_of_cells > 8) {
        derefine_all_grid (the_grid);
        number_of_cells = the_grid->number_of_cells;
    }

    struct cell_node *grid_cell = the_grid->first_cell;

    if(grid_cell) {

        // Deleting transition nodes.
        free((struct transition_node *) (grid_cell->north));
        free((struct transition_node *) (grid_cell->front));
        free((struct transition_node *) (grid_cell->east));
        free((struct transition_node *) (((struct cell_node *) (grid_cell->west))->west));
        free((struct transition_node *) (((struct cell_node *) (grid_cell->south))->south));
        free((struct transition_node *) (((struct cell_node *) (grid_cell->back))->back));

        // Deleting cells nodes.
        while (grid_cell) {

            struct cell_node *next = grid_cell->next;
            free_cell_node(grid_cell);
            grid_cell = next;

        }
    }



    if (the_grid->refined_this_step) {
        uint32_vector_clear(the_grid->refined_this_step);
    }

    if (the_grid->free_sv_positions) {
        uint32_vector_clear(the_grid->free_sv_positions);
    }

}

void clean_and_free_grid(struct grid* the_grid) {
    clean_grid(the_grid);

    if (the_grid->active_cells) {
        free (the_grid->active_cells);
    }

    if (the_grid->refined_this_step) {
        uint32_vector_clear_and_free_data(the_grid->refined_this_step);
        free (the_grid->refined_this_step);
    }

    if (the_grid->free_sv_positions) {
        uint32_vector_clear_and_free_data(the_grid->free_sv_positions);
        free (the_grid->free_sv_positions);
    }

    free (the_grid);
}

// Prints grid discretization matrix.
void print_grid_matrix (struct grid *the_grid, FILE *output_file) {
    if (!output_file) {
        fprintf (stderr, "print_grid_matrix: output_file is NULL! Open it first!");
        return;
    }

    struct cell_node *grid_cell;
    grid_cell = the_grid->first_cell;
    struct element element;
    struct element *cell_elements;
    int max_el = the_grid->num_cell_neighbours;

    while (grid_cell != 0) {
        if (grid_cell->active) {

            cell_elements = grid_cell->elements;
            element = cell_elements[0];

            fprintf(output_file, "%" PRIu64 " " "%" PRIu64 " %.15lf\n",
                    grid_cell->grid_position + 1,
                    (element.column) + 1,
                    element.value);

            int el_count = 1;

            while ((el_count < max_el) && (cell_elements[el_count].cell != NULL)) {

                element = grid_cell->elements[el_count];
                fprintf(output_file, "%" PRIu64 " " "%" PRIu64 " %.15lf\n",
                        grid_cell->grid_position + 1,
                        (element.column) + 1,
                        element.value);

                /*fprintf (output_file,
                         " %.6lf ("
                         "%" PRIu64 ","
                         "%" PRIu64 ") ",
                         element.value, grid_cell->grid_position + 1, (element.column) + 1);*/

                el_count++;
            }
            //fprintf (output_file, "\n");
        }
        grid_cell = grid_cell->next;
    }
    //fprintf (output_file, "________________________________________________________________________\n");
}


void print_grid_vector(struct grid* the_grid, FILE *output_file, char name)
{
    struct cell_node *grid_cell;
    grid_cell = the_grid->first_cell;

    while( grid_cell != 0 )
    {
        if( grid_cell->active )
        {
            if(name == 'b')
                fprintf(output_file, "%.15lf\n", grid_cell->b);
            else if (name == 'x')
                fprintf(output_file, "%.15lf\n", grid_cell->v);
        }
        grid_cell = grid_cell->next;

    }

}

double * grid_vector_to_array(struct grid* the_grid, char name, uint64_t *num_lines) {
    struct cell_node *grid_cell;
    grid_cell = the_grid->first_cell;

    *num_lines = the_grid->num_active_cells;
    double *vector = (double*) malloc(*num_lines*sizeof(double));

    while( grid_cell != 0 )
    {
        if( grid_cell->active )
        {
            if(name == 'b')
                vector[grid_cell->grid_position] = grid_cell->b;
            else if (name == 'x')
                vector[grid_cell->grid_position] = grid_cell->v;
        }
        grid_cell = grid_cell->next;

    }

    return vector;

}
