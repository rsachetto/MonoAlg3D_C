//
// Created by sachetto on 13/10/17.
//

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>

#include "../alg/grid/grid.h"
#include "../config/assembly_matrix_config.h"
#include "../libraries_common/config_helpers.h"
#include "../monodomain/constants.h"
#include "../utils/utils.h"


static inline double ALPHA(double beta, double cm, double dt, double h) {
    return (((beta * cm) / dt) * UM2_TO_CM2) * pow(h, 3.0);
}

void initialize_diagonal_elements(struct monodomain_solver *the_solver, struct grid *the_grid) {

    double alpha, h;
    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;
    double beta = the_solver->beta;
    double cm = the_solver->cm;

    double dt = the_solver->dt;

    int i;

#pragma omp parallel for private(alpha, h)
    for(i = 0; i < num_active_cells; i++) {
        h = ac[i]->face_length;
        alpha = ALPHA(beta, cm, dt, h);

        struct element element;
        element.column = ac[i]->grid_position;
        element.cell = ac[i];
        element.value = alpha;

        if(ac[i]->elements)
            sb_free(ac[i]->elements);

        ac[i]->elements = NULL;

        sb_reserve(ac[i]->elements, 7);
        sb_push(ac[i]->elements, element);
    }
}

static void fill_discretization_matrix_elements(double sigma_x, double sigma_y, double sigma_z,
                                                struct cell_node *grid_cell, void *neighbour_grid_cell,
                                                char direction) {

    uint32_t position;
    bool has_found;
    double h;

    struct transition_node *white_neighbor_cell;
    struct cell_node *black_neighbor_cell;

    double sigma_x1 = (2.0f * sigma_x * sigma_x) / (sigma_x + sigma_x);
    double sigma_x2 = (2.0f * sigma_x * sigma_x) / (sigma_x + sigma_x);
    double sigma_y1 = (2.0f * sigma_y * sigma_y) / (sigma_y + sigma_y);
    double sigma_y2 = (2.0f * sigma_y * sigma_y) / (sigma_y + sigma_y);
    double sigma_z1 = (2.0f * sigma_z * sigma_z) / (sigma_z + sigma_z);
    double sigma_z2 = (2.0f * sigma_z * sigma_z) / (sigma_z + sigma_z);

    /* When neighbour_grid_cell is a transition node, looks for the next neighbor
     * cell which is a cell node. */
    uint16_t neighbour_grid_cell_level = ((struct basic_cell_data *)(neighbour_grid_cell))->level;
    char neighbour_grid_cell_type = ((struct basic_cell_data *)(neighbour_grid_cell))->type;

    if(neighbour_grid_cell_level > grid_cell->cell_data.level) {
        if(neighbour_grid_cell_type == TRANSITION_NODE_TYPE) {
            has_found = false;
            while(!has_found) {
                if(neighbour_grid_cell_type == TRANSITION_NODE_TYPE) {
                    white_neighbor_cell = (struct transition_node *)neighbour_grid_cell;
                    if(white_neighbor_cell->single_connector == NULL) {
                        has_found = true;
                    } else {
                        neighbour_grid_cell = white_neighbor_cell->quadruple_connector1;
                        neighbour_grid_cell_type = ((struct basic_cell_data *)(neighbour_grid_cell))->type;
                    }
                } else {
                    break;
                }
            }
        }
    } else {
        if(neighbour_grid_cell_level <= grid_cell->cell_data.level &&
           (neighbour_grid_cell_type == TRANSITION_NODE_TYPE)) {
            has_found = false;
            while(!has_found) {
                if(neighbour_grid_cell_type == TRANSITION_NODE_TYPE) {
                    white_neighbor_cell = (struct transition_node *)(neighbour_grid_cell);
                    if(white_neighbor_cell->single_connector == 0) {
                        has_found = true;
                    } else {
                        neighbour_grid_cell = white_neighbor_cell->single_connector;
                        neighbour_grid_cell_type = ((struct basic_cell_data *)(neighbour_grid_cell))->type;
                    }
                } else {
                    break;
                }
            }
        }
    }

    // We care only with the interior points
    if(neighbour_grid_cell_type == CELL_NODE_TYPE) {

        black_neighbor_cell = (struct cell_node *)(neighbour_grid_cell);

        if(black_neighbor_cell->active) {

            if(black_neighbor_cell->cell_data.level > grid_cell->cell_data.level) {
                h = black_neighbor_cell->face_length;
            } else {
                h = grid_cell->face_length;
            }

            lock_cell_node(grid_cell);

            struct element *cell_elements = grid_cell->elements;
            position = black_neighbor_cell->grid_position;

            size_t max_elements = sb_count(cell_elements);
            bool insert = true;

            for(size_t i = 1; i < max_elements; i++) {
                if(cell_elements[i].column == position) {
                    insert = false;
                    break;
                }
            }

            if(insert) {

                struct element new_element;
                new_element.column = position;
                if(direction == 'n') { // Z direction
                    new_element.value = -sigma_z1 * h;
                    cell_elements[0].value += (sigma_z1 * h);
                } else if(direction == 's') { // Z direction
                    new_element.value = -sigma_z2 * h;
                    cell_elements[0].value += sigma_z2 * h;
                } else if(direction == 'e') { // Y direction
                    new_element.value = -sigma_y1 * h;
                    cell_elements[0].value += (sigma_y1 * h);
                } else if(direction == 'w') { // Y direction
                    new_element.value = -sigma_y2 * h;
                    cell_elements[0].value += (sigma_y2 * h);
                } else if(direction == 'f') { // X direction
                    new_element.value = -sigma_x1 * h;
                    cell_elements[0].value += (sigma_x1 * h);
                } else if(direction == 'b') { // X direction
                    new_element.value = -sigma_x2 * h;
                    cell_elements[0].value += (sigma_x2 * h);
                }

                new_element.cell = black_neighbor_cell;
                sb_push(grid_cell->elements, new_element);
            }
            unlock_cell_node(grid_cell);

            lock_cell_node(black_neighbor_cell);
            cell_elements = black_neighbor_cell->elements;
            position = grid_cell->grid_position;

            max_elements = sb_count(cell_elements);

            insert = true;
            for(size_t i = 1; i < max_elements; i++) {
                if(cell_elements[i].column == position) {
                    insert = false;
                    break;
                }
            }

            if(insert) {

                struct element new_element;
                new_element.column = position;
                if(direction == 'n') { // Z direction
                    new_element.value = -sigma_z1 * h;
                    cell_elements[0].value += (sigma_z1 * h);
                } else if(direction == 's') { // Z direction
                    new_element.value = -sigma_z2 * h;
                    cell_elements[0].value += (sigma_z2 * h);
                } else if(direction == 'e') { // Y direction
                    new_element.value = -sigma_y1 * h;
                    cell_elements[0].value += (sigma_y1 * h);
                } else if(direction == 'w') { // Y direction
                    new_element.value = -sigma_y2 * h;
                    cell_elements[0].value += (sigma_y2 * h);
                } else if(direction == 'f') { // X direction
                    new_element.value = -sigma_x1 * h;
                    cell_elements[0].value += (sigma_x1 * h);
                } else if(direction == 'b') { // X direction
                    new_element.value = -sigma_x2 * h;
                    cell_elements[0].value += (sigma_x2 * h);
                }

                new_element.cell = grid_cell;
                sb_push(black_neighbor_cell->elements, new_element);
            }

            unlock_cell_node(black_neighbor_cell);
        }
    }
}

ASSEMBLY_MATRIX(random_sigma_discretization_matrix) {

    printf("Assembling matrix 1 time\n");

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    int i;

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data.config, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config->config_data.config, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config->config_data.config, "sigma_z");

    srand(time(NULL));




#pragma omp parallel for
    for(i = 0; i < num_active_cells; i++) {

        float r;

        #pragma omp critical
        r = (float)rand()/(float)RAND_MAX;

        real sigma_x_new = sigma_x*r;

        #pragma omp critical
        r = (float)rand()/(float)RAND_MAX;

        real sigma_y_new= sigma_y*r;
        real sigma_z_new= sigma_z*r;


        // Computes and designates the flux due to south cells.
        fill_discretization_matrix_elements(sigma_x_new, sigma_y_new, sigma_z_new, ac[i], ac[i]->south, 's');

        // Computes and designates the flux due to north cells.
        fill_discretization_matrix_elements(sigma_x_new, sigma_y_new, sigma_z_new, ac[i], ac[i]->north, 'n');

        // Computes and designates the flux due to east cells.
        fill_discretization_matrix_elements(sigma_x_new, sigma_y_new, sigma_z_new, ac[i], ac[i]->east, 'e');

        // Computes and designates the flux due to west cells.
        fill_discretization_matrix_elements(sigma_x_new, sigma_y_new, sigma_z_new, ac[i], ac[i]->west, 'w');

        // Computes and designates the flux due to front cells.
        fill_discretization_matrix_elements(sigma_x_new, sigma_y_new, sigma_z_new, ac[i], ac[i]->front, 'f');

        // Computes and designates the flux due to back cells.
        fill_discretization_matrix_elements(sigma_x_new, sigma_y_new, sigma_z_new, ac[i], ac[i]->back, 'b');
    }
}

ASSEMBLY_MATRIX(no_fibers_assembly_matrix) {

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    int i;

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data.config, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config->config_data.config, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config->config_data.config, "sigma_z");



#pragma omp parallel for
    for(i = 0; i < num_active_cells; i++) {

        // Computes and designates the flux due to south cells.
        fill_discretization_matrix_elements(sigma_x, sigma_y, sigma_z, ac[i], ac[i]->south, 's');

        // Computes and designates the flux due to north cells.
        fill_discretization_matrix_elements(sigma_x, sigma_y, sigma_z, ac[i], ac[i]->north, 'n');

        // Computes and designates the flux due to east cells.
        fill_discretization_matrix_elements(sigma_x, sigma_y, sigma_z, ac[i], ac[i]->east, 'e');

        // Computes and designates the flux due to west cells.
        fill_discretization_matrix_elements(sigma_x, sigma_y, sigma_z, ac[i], ac[i]->west, 'w');

        // Computes and designates the flux due to front cells.
        fill_discretization_matrix_elements(sigma_x, sigma_y, sigma_z, ac[i], ac[i]->front, 'f');

        // Computes and designates the flux due to back cells.
        fill_discretization_matrix_elements(sigma_x, sigma_y, sigma_z, ac[i], ac[i]->back, 'b');
    }
}
