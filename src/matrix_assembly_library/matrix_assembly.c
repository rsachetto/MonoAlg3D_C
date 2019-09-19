//
// Created by sachetto on 13/10/17.
//

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>

#include "../alg/grid/grid.h"
#include "../config/assembly_matrix_config.h"
#include "../utils/utils.h"
#include "../single_file_libraries/stb_ds.h"
#include "../libraries_common/common_data_structures.h"

#include "../config_helpers/config_helpers.h"

INIT_ASSEMBLY_MATRIX(set_initial_conditions_fvm) {

    real_cpu alpha;
    struct cell_node **ac = the_grid->active_cells;
    uint32_t active_cells = the_grid->num_active_cells;
    real_cpu beta = the_solver->beta;
    real_cpu cm = the_solver->cm;
    real_cpu dt = the_solver->dt;
    int i;

    #pragma omp parallel for private(alpha)
    for(i = 0; i < active_cells; i++) {

        alpha = ALPHA(beta, cm, dt, ac[i]->discretization.x, ac[i]->discretization.y, ac[i]->discretization.z);
        ac[i]->v = initial_v;
        ac[i]->b = initial_v * alpha;
    }
}

static struct element fill_element(uint32_t position, char direction, real_cpu dx, real_cpu dy, real_cpu dz,\
                                   real_cpu sigma_x, real_cpu sigma_y, real_cpu sigma_z,\
                                   struct element *cell_elements);


void initialize_diagonal_elements(struct monodomain_solver *the_solver, struct grid *the_grid) {

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;
    real_cpu beta = the_solver->beta;
    real_cpu cm = the_solver->cm;
    real_cpu dt = the_solver->dt;

    uint32_t i;

    #pragma omp parallel for
    for(i = 0; i < num_active_cells; i++) {
        real_cpu alpha, dx, dy, dz;

        dx = ac[i]->discretization.x;
        dy = ac[i]->discretization.y;
        dz = ac[i]->discretization.z;

        alpha = ALPHA(beta, cm, dt, dx, dy, dz);

        struct element element;
        element.column = ac[i]->grid_position;
        element.cell = ac[i];
        element.value = alpha;

        if(ac[i]->elements)
            arrfree(ac[i]->elements);

        ac[i]->elements = NULL;

        arrsetcap(ac[i]->elements, 7);
        arrput(ac[i]->elements, element);
    }
}

struct element fill_element(uint32_t position, char direction, real_cpu dx, real_cpu dy, real_cpu dz, real_cpu sigma_x,
                            real_cpu sigma_y, real_cpu sigma_z, struct element *cell_elements) {

    real_cpu multiplier;

    struct element new_element;
    new_element.column = position;

    if(direction == 'n') { // Z direction
        multiplier = ((dx * dy) / dz);
        new_element.value = -sigma_z * multiplier;
        cell_elements[0].value += (sigma_z * multiplier);
    } else if(direction == 's') { // Z direction
        multiplier = ((dx * dy) / dz);
        new_element.value = -sigma_z * multiplier;
        cell_elements[0].value += (sigma_z * multiplier);
    } else if(direction == 'e') { // Y direction
        multiplier = ((dx * dz) / dy);
        new_element.value = -sigma_y * multiplier;
        cell_elements[0].value += (sigma_y * multiplier);
    } else if(direction == 'w') { // Y direction
        multiplier = ((dx * dz) / dy);
        new_element.value = -sigma_y * multiplier;
        cell_elements[0].value += (sigma_y * multiplier);
    } else if(direction == 'f') { // X direction
        multiplier = ((dy * dz) / dx);
        new_element.value = -sigma_x * ((dy * dz) / dx);
        cell_elements[0].value += (sigma_x * multiplier);
    } else if(direction == 'b') { // X direction
        multiplier = ((dy * dz) / dx);
        new_element.value = -sigma_x * multiplier;
        cell_elements[0].value += (sigma_x * multiplier);
    }
    return new_element;
}

static void fill_discretization_matrix_elements(struct cell_node *grid_cell, void *neighbour_grid_cell, char direction) {

    bool has_found;

    struct transition_node *white_neighbor_cell;
    struct cell_node *black_neighbor_cell;

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

            uint32_t position;
            real_cpu dx, dy, dz;

            real_cpu sigma_x1 = grid_cell->sigma.x;
            real_cpu sigma_x2 = black_neighbor_cell->sigma.x;
            real_cpu sigma_x = 0.0;
            
            if(sigma_x1 != 0.0 && sigma_x2 != 0.0) {
                sigma_x = (2.0f * sigma_x1 * sigma_x2) / (sigma_x1 + sigma_x2);
            }

            real_cpu sigma_y1 = grid_cell->sigma.y;
            real_cpu sigma_y2 = black_neighbor_cell->sigma.y;
            real_cpu sigma_y = 0.0;

            if(sigma_y1 != 0.0 && sigma_y2 != 0.0) {
                sigma_y = (2.0f * sigma_y1 * sigma_y2) / (sigma_y1 + sigma_y2);
            }

            real_cpu sigma_z1 = grid_cell->sigma.z;
            real_cpu sigma_z2 = black_neighbor_cell->sigma.z;
            real_cpu sigma_z = 0.0;

            if(sigma_z1 != 0.0 && sigma_z2 != 0.0) {
                sigma_z = (2.0f * sigma_z1 * sigma_z2) / (sigma_z1 + sigma_z2);
            }
            
            if(black_neighbor_cell->cell_data.level > grid_cell->cell_data.level) {
                dx = black_neighbor_cell->discretization.x;
                dy = black_neighbor_cell->discretization.y;
                dz = black_neighbor_cell->discretization.z;
            } else {
                dx = grid_cell->discretization.x;
                dy = grid_cell->discretization.y;
                dz = grid_cell->discretization.z;
            }

            lock_cell_node(grid_cell);

            struct element *cell_elements = grid_cell->elements;
            position = black_neighbor_cell->grid_position;

            size_t max_elements = arrlen(cell_elements);
            bool insert = true;

            for(size_t i = 1; i < max_elements; i++) {
                if(cell_elements[i].column == position) {
                    insert = false;
                    break;
                }
            }

            if(insert) {

                struct element new_element = fill_element(position, direction, dx, dy, dz, sigma_x, sigma_y, sigma_z, cell_elements);

                new_element.cell = black_neighbor_cell;
                arrput(grid_cell->elements, new_element);
            }
            unlock_cell_node(grid_cell);

            lock_cell_node(black_neighbor_cell);
            cell_elements = black_neighbor_cell->elements;
            position = grid_cell->grid_position;

            max_elements = arrlen(cell_elements);

            insert = true;
            for(size_t i = 1; i < max_elements; i++) {
                if(cell_elements[i].column == position) {
                    insert = false;
                    break;
                }
            }

            if(insert) {

                struct element new_element = fill_element(position, direction, dx, dy, dz, sigma_x, sigma_y, sigma_z, cell_elements);

                new_element.cell = grid_cell;
                arrput(black_neighbor_cell->elements, new_element);
            }

            unlock_cell_node(black_neighbor_cell);
        }
    }
}

int randRange(int n) {
    int limit;
    int r;

    limit = RAND_MAX - (RAND_MAX % n);

    while((r = rand()) >= limit)
        ;

    return r % n;
}

ASSEMBLY_MATRIX(random_sigma_discretization_matrix) {

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    int i;

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config->config_data, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config->config_data, "sigma_z");

    srand((unsigned int)time(NULL));

    real_cpu modifiers[4] = {0.0f, 0.1f, 0.5f, 1.0f};

#pragma omp parallel for
    for(i = 0; i < num_active_cells; i++) {
        real_cpu r;

        #pragma omp critical
        r = modifiers[randRange(4)];

        real sigma_x_new = sigma_x * r;
        real sigma_y_new = sigma_y * r;
        real sigma_z_new = sigma_z * r;

        ac[i]->sigma.x = sigma_x_new;
        ac[i]->sigma.y = sigma_y_new;
        ac[i]->sigma.z = sigma_z_new;

    }

#pragma omp parallel for
    for(i = 0; i < num_active_cells; i++) {

        // Computes and designates the flux due to south cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->south, 's');

        // Computes and designates the flux due to north cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->north, 'n');

        // Computes and designates the flux due to east cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->east, 'e');

        // Computes and designates the flux due to west cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->west, 'w');

        // Computes and designates the flux due to front cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->front, 'f');

        // Computes and designates the flux due to back cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->back, 'b');
    }
}

ASSEMBLY_MATRIX(source_sink_discretization_matrix) 
{

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    int i;

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config->config_data, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config->config_data, "sigma_z");

    real channel_width = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, channel_width, config->config_data, "channel_width");

    real channel_length = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, channel_length, config->config_data, "channel_length");

    bool inside;

    real side_length_x = the_grid->mesh_side_length.x;
    real side_length_y = the_grid->mesh_side_length.y;
    real side_length_z = the_grid->mesh_side_length.z;

    real region_height = (side_length_y - channel_width) / 2.0;

    #pragma omp parallel for
    for (i = 0; i < num_active_cells; i++)
    {
        real sigma_x_new;
        real sigma_y_new;
        real sigma_z_new;

        double x = ac[i]->center.x;
        double y = ac[i]->center.y;
        double z = ac[i]->center.z;

        // Check region 1
        inside = (x >= 0.0) && (x <= channel_length) && (y >= 0.0) && (y <= region_height);
        
        // Check region 2
        inside |= (x >= 0.0) && (x <= channel_length) && (y >= region_height + channel_width) && (y <= side_length_y);

        if (inside)
        {
            sigma_x_new = 0.0;
            sigma_y_new = 0.0;
            sigma_z_new = 0.0;
        }
        else
        {
            sigma_x_new = sigma_x;
            sigma_y_new = sigma_y;
            sigma_z_new = sigma_z;
        }

        ac[i]->sigma.x = sigma_x_new;
        ac[i]->sigma.y = sigma_y_new;
        ac[i]->sigma.z = sigma_z_new;

    }

    // Then, we fill the discretization matrix
    #pragma omp parallel for
    for(i = 0; i < num_active_cells; i++) 
    {

        // Computes and designates the flux due to south cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->south, 's');

        // Computes and designates the flux due to north cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->north, 'n');

        // Computes and designates the flux due to east cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->east, 'e');

        // Computes and designates the flux due to west cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->west, 'w');

        // Computes and designates the flux due to front cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->front, 'f');

        // Computes and designates the flux due to back cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->back, 'b');
    }
}

ASSEMBLY_MATRIX(source_sink_discretization_matrix_with_different_sigma) 
{

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    int i;

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config->config_data, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config->config_data, "sigma_z");

    real channel_width = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, channel_width, config->config_data, "channel_width");

    real channel_length = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, channel_length, config->config_data, "channel_length");

    real source_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, source_factor, config->config_data, "source_factor");

    real sink_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sink_factor, config->config_data, "sink_factor");

    bool inside_3, inside_4;

    real side_length_x = the_grid->mesh_side_length.x;
    real side_length_y = the_grid->mesh_side_length.y;
    real side_length_z = the_grid->mesh_side_length.z;

    real region_height = (side_length_y - channel_width) / 2.0;

    // Set the conductivities for each cell on the grid
    #pragma omp parallel for
    for(i = 0; i < num_active_cells; i++) 
    {
        real sigma_x_new = sigma_x;
        real sigma_y_new = sigma_y;
        real sigma_z_new = sigma_z;

        real x = ac[i]->center.x;
        real y = ac[i]->center.y;
        real z = ac[i]->center.z;

        // Check region 3
        inside_3 = (x >= 0.0) && (x < channel_length) && (y >= region_height) && (y <= region_height + channel_width);

        if (inside_3)
        {
            sigma_x_new = sigma_x * source_factor;
            sigma_y_new = sigma_y * source_factor;
            sigma_z_new = sigma_z * source_factor;
        }    

        // Check region 4
        inside_4 = (x >= channel_length) && (x <= side_length_x) && (y >= 0.0) && (y <= side_length_y);

        if (inside_4)
        {
            sigma_x_new = sigma_x * sink_factor;
            sigma_y_new = sigma_y * sink_factor;
            sigma_z_new = sigma_z * sink_factor;
        }

        ac[i]->sigma.x = sigma_x_new;
        ac[i]->sigma.y = sigma_y_new;
        ac[i]->sigma.z = sigma_z_new;

    }

    // Then, we fill the discretization matrix
    #pragma omp parallel for
    for(i = 0; i < num_active_cells; i++) 
    {

        // Computes and designates the flux due to south cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->south, 's');

        // Computes and designates the flux due to north cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->north, 'n');

        // Computes and designates the flux due to east cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->east, 'e');

        // Computes and designates the flux due to west cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->west, 'w');

        // Computes and designates the flux due to front cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->front, 'f');

        // Computes and designates the flux due to back cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->back, 'b');
    }

}

ASSEMBLY_MATRIX(homogeneous_sigma_assembly_matrix) {

    static bool sigma_initialized = false;

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    int i;

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config->config_data, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config->config_data, "sigma_z");

    if(!sigma_initialized) {
        #pragma omp parallel for
        for (i = 0; i < num_active_cells; i++) {
            ac[i]->sigma.x = sigma_x;
            ac[i]->sigma.y = sigma_y;
            ac[i]->sigma.z = sigma_z;
        }

        sigma_initialized = true;
    }

    #pragma omp parallel for
    for(i = 0; i < num_active_cells; i++) {

        // Computes and designates the flux due to south cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->south, 's');

        // Computes and designates the flux due to north cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->north, 'n');

        // Computes and designates the flux due to east cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->east, 'e');

        // Computes and designates the flux due to west cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->west, 'w');

        // Computes and designates the flux due to front cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->front, 'f');

        // Computes and designates the flux due to back cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->back, 'b');
    }
}

ASSEMBLY_MATRIX(homogeneous_sigma_with_a_factor_assembly_matrix) {

    static bool sigma_initialized = false;

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    int i;

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config->config_data, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config->config_data, "sigma_z");

    real sigma_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_factor, config->config_data, "sigma_factor");

    if(!sigma_initialized) {
        #pragma omp parallel for
        for (i = 0; i < num_active_cells; i++) {
            ac[i]->sigma.x = sigma_x * sigma_factor;
            ac[i]->sigma.y = sigma_y * sigma_factor;
            ac[i]->sigma.z = sigma_z * sigma_factor;
        }

        sigma_initialized = true;
    }

    #pragma omp parallel for
    for(i = 0; i < num_active_cells; i++) {

        // Computes and designates the flux due to south cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->south, 's');

        // Computes and designates the flux due to north cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->north, 'n');

        // Computes and designates the flux due to east cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->east, 'e');

        // Computes and designates the flux due to west cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->west, 'w');

        // Computes and designates the flux due to front cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->front, 'f');

        // Computes and designates the flux due to back cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->back, 'b');
    }
}

ASSEMBLY_MATRIX(fibrotic_region_with_sigma_factor_assembly_matrix) 
{

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    int i;

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config->config_data, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config->config_data, "sigma_z");

    real sigma_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_factor, config->config_data, "sigma_factor");

    real min_x = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, min_x, config->config_data, "region_min_x");

    real max_x = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, max_x, config->config_data, "region_max_x");

    real min_y = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, min_y, config->config_data, "region_min_y");

    real max_y = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, max_y, config->config_data, "region_max_y");

    real min_z = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, min_z, config->config_data, "region_min_z");

    real max_z = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, max_z, config->config_data, "region_max_z");

    bool inside;

    // Set the conductivities for each cell on the grid
    #pragma omp parallel for
    for(i = 0; i < num_active_cells; i++) 
    {
        real sigma_x_new = sigma_x;
        real sigma_y_new = sigma_y;
        real sigma_z_new = sigma_z;

        real x = ac[i]->center.x;
        real y = ac[i]->center.y;
        real z = ac[i]->center.z;

        // Check if inside the region
        inside = (x >= min_x) && (x <= max_x) &&\
                 (y >= min_y) && (y <= max_y) &&\
                 (z >= min_z) && (z <= max_z);

        if (inside)
        {
            sigma_x_new = sigma_x * sigma_factor;
            sigma_y_new = sigma_y * sigma_factor;
            sigma_z_new = sigma_z * sigma_factor;
        }    

        ac[i]->sigma.x = sigma_x_new;
        ac[i]->sigma.y = sigma_y_new;
        ac[i]->sigma.z = sigma_z_new;

    }

    // Then, we fill the discretization matrix
    #pragma omp parallel for
    for(i = 0; i < num_active_cells; i++) 
    {

        // Computes and designates the flux due to south cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->south, 's');

        // Computes and designates the flux due to north cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->north, 'n');

        // Computes and designates the flux due to east cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->east, 'e');

        // Computes and designates the flux due to west cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->west, 'w');

        // Computes and designates the flux due to front cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->front, 'f');

        // Computes and designates the flux due to back cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->back, 'b');
    }

}

ASSEMBLY_MATRIX(heterogenous_sigma_with_factor_assembly_matrix) 
{

    static bool sigma_initialized = false;

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    struct cell_node *grid_cell;

    initialize_diagonal_elements(the_solver, the_grid);

    int i;

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config->config_data, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config->config_data, "sigma_z");

    real_cpu phi = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, phi, config->config_data, "phi");

    unsigned seed = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, seed, config->config_data, "seed");

    real sigma_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_factor, config->config_data, "sigma_factor");

    print_to_stdout_and_file("Reducing conductivity from %.2lf %% of cells\n", phi * 100.0);

    // Initialize the seed for the fibrosis
    srand(seed);

    print_to_stdout_and_file("Using %u as seed\n", seed);

    if (!sigma_initialized)
    {
	    grid_cell = the_grid->first_cell;
	    while(grid_cell != 0) 
	    {

		    if(grid_cell->active) 
		    {
	    		real_cpu p = (real_cpu)(rand()) / (RAND_MAX);
	    		if (p < phi) 
			    {
				    grid_cell->sigma.x = sigma_x * sigma_factor;
				    grid_cell->sigma.y = sigma_y * sigma_factor;
				    grid_cell->sigma.z = sigma_z * sigma_factor;
	    		}
			    else
			    {
				    grid_cell->sigma.x = sigma_x;
	    			grid_cell->sigma.y = sigma_y;
	    			grid_cell->sigma.z = sigma_z;
			    }
		    }
		    grid_cell = grid_cell->next;
    	}

	    sigma_initialized = true;
    }

    #pragma omp parallel for
    for(i = 0; i < num_active_cells; i++) 
    {

        // Computes and designates the flux due to south cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->south, 's');

        // Computes and designates the flux due to north cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->north, 'n');

        // Computes and designates the flux due to east cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->east, 'e');

        // Computes and designates the flux due to west cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->west, 'w');

        // Computes and designates the flux due to front cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->front, 'f');

        // Computes and designates the flux due to back cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->back, 'b');
    }
}

// This function will read the fibrotic regions and for each cell that is inside the region we will
// reduce its conductivity value based on the 'sigma_factor'.
ASSEMBLY_MATRIX(heterogenous_sigma_with_factor_assembly_matrix_from_file) 
{

    static bool sigma_initialized = false;

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    char *fib_file = NULL;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(fib_file, config->config_data, "fibrosis_file");

    int fib_size = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(int, fib_size, config->config_data, "size");

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config->config_data, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config->config_data, "sigma_z");

    real_cpu phi = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, phi, config->config_data, "phi");

    unsigned seed = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, seed, config->config_data, "seed");

    real sigma_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_factor, config->config_data, "sigma_factor");
  
    if(!sigma_initialized) 
    {
        #pragma omp parallel for
        for (uint32_t i = 0; i < num_active_cells; i++) 
        {
            ac[i]->sigma.x = sigma_x;
            ac[i]->sigma.y = sigma_y;
            ac[i]->sigma.z = sigma_z;
        }

        sigma_initialized = true;
    }

    // Reading the fibrotic regions from the input file
    FILE *file = fopen(fib_file, "r");

    if(!file) 
    {
        printf("Error opening file %s!!\n", fib_file);
        exit(0);
    }

    real_cpu **scar_mesh = (real_cpu **)malloc(sizeof(real_cpu *) * fib_size);

    for(int i = 0; i < fib_size; i++) 
    {
        scar_mesh[i] = (real_cpu *)malloc(sizeof(real_cpu) * 7);
        if(scar_mesh[i] == NULL) 
        {
            printf("Failed to allocate memory\n");
            exit(0);
        }
    }

    uint32_t i = 0;
    while (fscanf(file, "%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &scar_mesh[i][0], &scar_mesh[i][1], &scar_mesh[i][2], &scar_mesh[i][3], &scar_mesh[i][4], &scar_mesh[i][5], &scar_mesh[i][6]) != EOF)
    {
        i++;
    }

    fclose(file);    

    uint32_t num_fibrotic_regions = i;

    // Pass through all the cells of the grid and check if its center is inside the current
    // fibrotic region
    #pragma omp parallel for
    for(int j = 0; j < num_fibrotic_regions; j++) 
    {
        
        struct cell_node *grid_cell = the_grid->first_cell;
        
        real_cpu b_center_x = scar_mesh[j][0];
        real_cpu b_center_y = scar_mesh[j][1];

        real_cpu b_h_dx = scar_mesh[j][3];
        real_cpu b_h_dy = scar_mesh[j][4];

        bool active = (bool) (scar_mesh[j][6]);

        while(grid_cell != 0) 
        {
            if (grid_cell->active)
            {
                real_cpu center_x = grid_cell->center.x;
                real_cpu center_y = grid_cell->center.y;
                real_cpu half_dx = grid_cell->discretization.x/2.0;
                real_cpu half_dy = grid_cell->discretization.y/2.0;

                struct point_3d p;
                struct point_3d q;

                p.x = b_center_y + b_h_dy;
                p.y = b_center_y - b_h_dy;

                q.x = b_center_x + b_h_dx;
                q.y = b_center_x - b_h_dx;

                // Check if the current cell is inside the fibrotic region
                if (center_x > q.y && center_x < q.x && center_y > p.y && center_y < p.x)
                {
                    if(active == 0)
                    {
                        grid_cell->sigma.x = sigma_x * sigma_factor;
                        grid_cell->sigma.y = sigma_y * sigma_factor;
                        grid_cell->sigma.z = sigma_z * sigma_factor;
                    }		
                }
            }
            
            grid_cell = grid_cell->next;
        } 
    }
		
    #pragma omp parallel for
    for(int i = 0; i < num_active_cells; i++) 
    {

        // Computes and designates the flux due to south cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->south, 's');

        // Computes and designates the flux due to north cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->north, 'n');

        // Computes and designates the flux due to east cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->east, 'e');

        // Computes and designates the flux due to west cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->west, 'w');

        // Computes and designates the flux due to front cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->front, 'f');

        // Computes and designates the flux due to back cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->back, 'b');
    }
    
    for(int k = 0; k < fib_size; k++) 
    {
        free(scar_mesh[k]);
    }

    free(scar_mesh);
}    


// This function will generate the fibrotic region file for the 120um x 120um grid by reescaling
// the original Scientific Reports 4b grid from 40000um side_length to 48000um
// rescale_factor = 1.2 
ASSEMBLY_MATRIX(heterogenous_fibrotic_region_file_write_with_input_file)
{

    static bool sigma_initialized = false;
    int num;

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    struct cell_node *grid_cell;

    initialize_diagonal_elements(the_solver, the_grid);

    char *fib_file = NULL;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(fib_file, config->config_data, "fibrosis_file");

    char *new_fib_file = NULL;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(new_fib_file, config->config_data, "rescaled_fibrosis_file");

    int fib_size = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(int, fib_size, config->config_data, "size");

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config->config_data, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config->config_data, "sigma_z");

    real sigma_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_factor, config->config_data, "sigma_factor");

    real rescale_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, rescale_factor, config->config_data, "rescale_factor");

    FILE *file = fopen(fib_file, "r");

    if(!file) 
    {
        printf("Error opening file %s!!\n", fib_file);
        exit(0);
    }

    // Read and store the original positions of the fibrotic regions
    real_cpu **scar_mesh = (real_cpu **)malloc(sizeof(real_cpu *) * fib_size);

    for(int i = 0; i < fib_size; i++) 
    {
        scar_mesh[i] = (real_cpu *)malloc(sizeof(real_cpu) * 7);
        if(scar_mesh[i] == NULL) 
        {
            printf("Failed to allocate memory\n");
            exit(0);
        }
    }

    for(int i = 0; i < fib_size; i++) 
    {
        fscanf(file, "%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &scar_mesh[i][0], &scar_mesh[i][1], &scar_mesh[i][2], &scar_mesh[i][3], &scar_mesh[i][4], &scar_mesh[i][5], &scar_mesh[i][6]);
    }

    fclose(file);  

    // Write the new fibrotic region file based on the 'rescale_factor'	
	FILE *fileW = fopen(new_fib_file, "w");
		
    if(!file) 
    {
        printf("Error opening file %s!!\n", fib_file);
        exit(0);
    }

    // Multiple the positions of each scar by a rescale factor
    for(int i = 0; i < fib_size; i++) 
    {
        scar_mesh[i][0] = scar_mesh[i][0]*rescale_factor;
        scar_mesh[i][1] = scar_mesh[i][1]*rescale_factor;
        scar_mesh[i][2] = scar_mesh[i][2]*rescale_factor;
        scar_mesh[i][3] = scar_mesh[i][3]*rescale_factor;
        scar_mesh[i][4] = scar_mesh[i][4]*rescale_factor;
        scar_mesh[i][5] = scar_mesh[i][5]*rescale_factor;
        
        fprintf(fileW, "%g,%g,%g,%g,%g,%g,%g\n", scar_mesh[i][0], scar_mesh[i][1], scar_mesh[i][2], scar_mesh[i][3], scar_mesh[i][4], scar_mesh[i][5], scar_mesh[i][6]);
    }
    
    fclose(fileW);  	
		
    for(int k = 0; k < fib_size; k++) 
    {
        free(scar_mesh[k]);
    }

    free(scar_mesh);

    // We just leave the program after this ...
    print_to_stdout_and_file("[!] Finish writing new fibrotic region file '%s'!\n",new_fib_file);
    exit(EXIT_SUCCESS);
}

// This function will generate the fibrotic region file for the 120um x 120um grid by reescaling
// the original Scientific Reports 4b grid from 40000um side_length to 48000um. The fibrotic region
// will be mapped using the same idea used on the 'domains_library' function by using a fixed seed
// for the random number generator.
// rescale_factor = 1.2
ASSEMBLY_MATRIX(heterogenous_fibrotic_region_file_write_using_seed)
{

    static bool sigma_initialized = false;
    int num;

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    struct cell_node *grid_cell;

    initialize_diagonal_elements(the_solver, the_grid);

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config->config_data, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config->config_data, "sigma_z");

    real_cpu phi = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, phi, config->config_data, "phi");

    unsigned seed = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, seed, config->config_data, "seed");

    char *new_fib_file = NULL;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(new_fib_file, config->config_data, "rescaled_fibrosis_file");

    double rescale_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real,rescale_factor, config->config_data, "rescale_factor");

    // Write the new fibrotic region file
	FILE *fileW = fopen(new_fib_file, "w+");

    grid_cell = the_grid->first_cell;

    // Initialize the random the generator with the same seed used by the original model
    srand(seed);
    while(grid_cell != 0) 
    {

        if(grid_cell->active) 
        {
            real_cpu p = (real_cpu)(rand()) / (RAND_MAX);
            if(p < phi) 
            {
                // We reescale the cell position using the 'rescale_factor'
                double center_x = grid_cell->center.x * rescale_factor;
                double center_y = grid_cell->center.y * rescale_factor;
                double center_z = grid_cell->center.z * rescale_factor;
                double dx = grid_cell->discretization.x * rescale_factor;
                double dy = grid_cell->discretization.y * rescale_factor;
                double dz = grid_cell->discretization.z * rescale_factor;
                
                // Then, we write only the fibrotic regions to the output file
                fprintf(fileW,"%g,%g,%g,%g,%g,%g,0\n",center_x,center_y,center_z,dx/2.0,dy/2.0,dz/2.0);
            }

        }
        grid_cell = grid_cell->next;
    }

    fclose(fileW);  	
		
    // We just leave the program after this ...
    print_to_stdout_and_file("[!] Finish writing fibrotic region file '%s'!\n",new_fib_file);
    exit(EXIT_SUCCESS);
}
