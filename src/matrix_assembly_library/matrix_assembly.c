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
#include "../single_file_libraries/stb_ds.h"

static struct element fill_element(uint32_t position, char direction, real_cpu dx, real_cpu dy, real_cpu dz, real_cpu sigma_x1,
                                   real_cpu sigma_x2, real_cpu sigma_y1, real_cpu sigma_y2, real_cpu sigma_z1, real_cpu sigma_z2,
                                   struct element *cell_elements);

real_cpu calculate_kappa (const real cell_length, const real h)
{
    return (pow(cell_length,4) - pow(h,4)) / (12.0*pow(cell_length,2));
}

void initialize_diagonal_elements_ddm (struct monodomain_solver *the_solver, struct grid *the_grid,\
                                    const real_cpu dx, const real_cpu dy, const real_cpu dz,\
                                    const real_cpu kappa_x, const real_cpu kappa_y, const real_cpu kappa_z,\
                                    const real_cpu sigma_x, const real_cpu sigma_y, const real_cpu sigma_z)
{
    real_cpu alpha;
    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;
//    double beta = the_solver->beta;
//    double cm = the_solver->cm;

    real_cpu dt = the_solver->dt;

    int i;

    #pragma omp parallel for private(alpha)
    for(i = 0; i < num_active_cells; i++)
    {
        alpha = ALPHA_CM(beta, cm, dt, dx, dy, dz);
        
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

struct element fill_element_ddm (uint32_t position, char direction, real_cpu dx, real_cpu dy, real_cpu dz,\
                                 real_cpu sigma_x1, real_cpu sigma_x2, real_cpu sigma_y1, real_cpu sigma_y2, real_cpu sigma_z1, real_cpu sigma_z2,\
                                 const real_cpu kappa_x, const real_cpu kappa_y, const real_cpu kappa_z,\
                                 const real_cpu dt,\
                                struct element *cell_elements) 
{

    real_cpu multiplier;

    struct element new_element;
    new_element.column = position;
    new_element.direction = direction;

    // Z direction
    if(direction == 'n') 
    { 
        multiplier = ((dx * dy) / dz);
        new_element.value = ( multiplier * (-sigma_z1 - (kappa_z / dt)) );
        cell_elements[0].value += ( multiplier * (sigma_z1 + (kappa_z / dt)) );
    } 
    // Z direction
    else if(direction == 's') 
    { 
        multiplier = ((dx * dy) / dz);
        new_element.value = ( multiplier * (-sigma_z2 - (kappa_z / dt)) );
        cell_elements[0].value += ( multiplier * (sigma_z2 + (kappa_z / dt)) );
    } 
    // Y direction
    else if(direction == 'e') 
    { 
        multiplier = ((dx * dz) / dy);
        new_element.value = ( multiplier * (-sigma_y1 - (kappa_y / dt)) );
        cell_elements[0].value += ( multiplier * (sigma_y1 + (kappa_y / dt)) );
    } 
    // Y direction
    else if(direction == 'w') 
    { 
        multiplier = ((dx * dz) / dy);
        new_element.value = ( multiplier * (-sigma_y2 - (kappa_y / dt)) );
        cell_elements[0].value += ( multiplier * (sigma_y2 + (kappa_y / dt)) );
    }
    // X direction 
    else if(direction == 'f') 
    { 
        multiplier = ((dy * dz) / dx);
        new_element.value = ( multiplier * (-sigma_x1 - (kappa_x / dt)) );
        cell_elements[0].value += ( multiplier * (sigma_x1 + (kappa_x / dt)) );
    } 
    // X direction
    else if(direction == 'b') 
    { 
        multiplier = ((dy * dz) / dx);
        new_element.value = ( multiplier * (-sigma_x2 - (kappa_x / dt)) );
        cell_elements[0].value += ( multiplier * (sigma_x2 + (kappa_x / dt)) );
    }
    
    return new_element;
}

static void fill_discretization_matrix_elements_ddm (real_cpu sigma_x, real_cpu sigma_y, real_cpu sigma_z,
                                                const real_cpu kappa_x, const real_cpu kappa_y, const real_cpu kappa_z,
                                                const real_cpu dt,
                                                struct cell_node *grid_cell, void *neighbour_grid_cell,
                                                char direction) 
{

    uint32_t position;
    bool has_found;
    real_cpu dx, dy, dz;

    struct transition_node *white_neighbor_cell;
    struct cell_node *black_neighbor_cell;

    real_cpu sigma_x1 = 0.0;
    real_cpu sigma_x2 = 0.0;

    if(sigma_x != 0.0) 
    {
        sigma_x1 = (2.0f * sigma_x * sigma_x) / (sigma_x + sigma_x);
        sigma_x2 = (2.0f * sigma_x * sigma_x) / (sigma_x + sigma_x);
    }

    real_cpu sigma_y1 = 0.0;
    real_cpu sigma_y2 = 0.0;

    if(sigma_y != 0.0) 
    {
        sigma_y1 = (2.0f * sigma_y * sigma_y) / (sigma_y + sigma_y);
        sigma_y2 = (2.0f * sigma_y * sigma_y) / (sigma_y + sigma_y);
    }

    real_cpu sigma_z1 = 0.0;
    real_cpu sigma_z2 = 0.0;

    if(sigma_z != 0.0) 
    {
        sigma_z1 = (2.0f * sigma_z * sigma_z) / (sigma_z + sigma_z);
        sigma_z2 = (2.0f * sigma_z * sigma_z) / (sigma_z + sigma_z);
    }

    /* When neighbour_grid_cell is a transition node, looks for the next neighbor
     * cell which is a cell node. */
    uint16_t neighbour_grid_cell_level = ((struct basic_cell_data *)(neighbour_grid_cell))->level;
    char neighbour_grid_cell_type = ((struct basic_cell_data *)(neighbour_grid_cell))->type;

    if(neighbour_grid_cell_level > grid_cell->cell_data.level) 
    {
        if(neighbour_grid_cell_type == TRANSITION_NODE_TYPE) 
        {
            has_found = false;
            while(!has_found) 
            {
                if(neighbour_grid_cell_type == TRANSITION_NODE_TYPE) 
                {
                    white_neighbor_cell = (struct transition_node *)neighbour_grid_cell;
                    if(white_neighbor_cell->single_connector == NULL) 
                    {
                        has_found = true;
                    } 
                    else 
                    {
                        neighbour_grid_cell = white_neighbor_cell->quadruple_connector1;
                        neighbour_grid_cell_type = ((struct basic_cell_data *)(neighbour_grid_cell))->type;
                    }
                } 
                else 
                {
                    break;
                }
            }
        }
    } 
    else 
    {
        if(neighbour_grid_cell_level <= grid_cell->cell_data.level &&
           (neighbour_grid_cell_type == TRANSITION_NODE_TYPE)) 
           {
            has_found = false;
            while(!has_found) 
            {
                if(neighbour_grid_cell_type == TRANSITION_NODE_TYPE) 
                {
                    white_neighbor_cell = (struct transition_node *)(neighbour_grid_cell);
                    if(white_neighbor_cell->single_connector == 0) 
                    {
                        has_found = true;
                    } 
                    else 
                    {
                        neighbour_grid_cell = white_neighbor_cell->single_connector;
                        neighbour_grid_cell_type = ((struct basic_cell_data *)(neighbour_grid_cell))->type;
                    }
                } 
                else 
                {
                    break;
                }
            }
        }
    }

    // We care only with the interior points
    if(neighbour_grid_cell_type == CELL_NODE_TYPE) 
    {

        black_neighbor_cell = (struct cell_node *)(neighbour_grid_cell);

        if(black_neighbor_cell->active) 
        {

            if(black_neighbor_cell->cell_data.level > grid_cell->cell_data.level) 
            {
                dx = black_neighbor_cell->dx;
                dy = black_neighbor_cell->dy;
                dz = black_neighbor_cell->dz;
            } 
            else 
            {
                dx = grid_cell->dx;
                dy = grid_cell->dy;
                dz = grid_cell->dz;
            }

            lock_cell_node(grid_cell);

            struct element *cell_elements = grid_cell->elements;
            position = black_neighbor_cell->grid_position;

            size_t max_elements = arrlen(cell_elements);
            bool insert = true;

            for(size_t i = 1; i < max_elements; i++) 
            {
                if(cell_elements[i].column == position) 
                {
                    insert = false;
                    break;
                }
            }

            if(insert) 
            {

                struct element new_element = fill_element_ddm(position, direction, dx, dy, dz,\
                                            sigma_x1, sigma_x2, sigma_y1, sigma_y2, sigma_z1, sigma_z2,\
                                            kappa_x, kappa_y, kappa_z, dt,\
                                            cell_elements);

                new_element.cell = black_neighbor_cell;
                arrput(grid_cell->elements, new_element);
            }
            unlock_cell_node(grid_cell);

            lock_cell_node(black_neighbor_cell);
            cell_elements = black_neighbor_cell->elements;
            position = grid_cell->grid_position;

            max_elements = arrlen(cell_elements);

            insert = true;
            for(size_t i = 1; i < max_elements; i++) 
            {
                if(cell_elements[i].column == position) 
                {
                    insert = false;
                    break;
                }
            }

            if(insert) 
            {

                struct element new_element = fill_element_ddm(position, direction, dx, dy, dz,\
                                            sigma_x1, sigma_x2, sigma_y1, sigma_y2, sigma_z1, sigma_z2,\
                                            kappa_x, kappa_y, kappa_z, dt,\
                                            cell_elements);

                new_element.cell = grid_cell;
                arrput(black_neighbor_cell->elements, new_element);
            }

            unlock_cell_node(black_neighbor_cell);
        }
    }
}

void initialize_diagonal_elements(struct monodomain_solver *the_solver, struct grid *the_grid) {

    real_cpu alpha, dx, dy, dz;
    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;
    real_cpu beta = the_solver->beta;
    real_cpu cm = the_solver->cm;

    real_cpu dt = the_solver->dt;

    int i;

#pragma omp parallel for private(alpha, dx, dy, dz)
    for(i = 0; i < num_active_cells; i++) {
        dx = ac[i]->dx;
        dy = ac[i]->dy;
        dz = ac[i]->dz;

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

struct element fill_element(uint32_t position, char direction, real_cpu dx, real_cpu dy, real_cpu dz, real_cpu sigma_x1,
                            real_cpu sigma_x2, real_cpu sigma_y1, real_cpu sigma_y2, real_cpu sigma_z1, real_cpu sigma_z2,
                            struct element *cell_elements) {

    real_cpu multiplier;

    struct element new_element;
    new_element.column = position;

    if(direction == 'n') { // Z direction
        multiplier = ((dx * dy) / dz);
        new_element.value = -sigma_z1 * multiplier;
        cell_elements[0].value += (sigma_z1 * multiplier);
    } else if(direction == 's') { // Z direction
        multiplier = ((dx * dy) / dz);
        new_element.value = -sigma_z2 * multiplier;
        cell_elements[0].value += (sigma_z2 * multiplier);
    } else if(direction == 'e') { // Y direction
        multiplier = ((dx * dz) / dy);
        new_element.value = -sigma_y1 * multiplier;
        cell_elements[0].value += (sigma_y1 * multiplier);
    } else if(direction == 'w') { // Y direction
        multiplier = ((dx * dz) / dy);
        new_element.value = -sigma_y2 * multiplier;
        cell_elements[0].value += (sigma_y2 * multiplier);
    } else if(direction == 'f') { // X direction
        multiplier = ((dy * dz) / dx);
        new_element.value = -sigma_x1 * ((dy * dz) / dx);
        cell_elements[0].value += (sigma_x1 * multiplier);
    } else if(direction == 'b') { // X direction
        multiplier = ((dy * dz) / dx);
        new_element.value = -sigma_x2 * multiplier;
        cell_elements[0].value += (sigma_x2 * multiplier);
    }
    return new_element;
}

static void fill_discretization_matrix_elements(real_cpu sigma_x, real_cpu sigma_y, real_cpu sigma_z,
                                                struct cell_node *grid_cell, void *neighbour_grid_cell,
                                                char direction) {

    uint32_t position;
    bool has_found;
    real_cpu dx, dy, dz;

    struct transition_node *white_neighbor_cell;
    struct cell_node *black_neighbor_cell;

    real_cpu sigma_x1 = 0.0;
    real_cpu sigma_x2 = 0.0;

    if(sigma_x != 0.0) {
        sigma_x1 = (2.0f * sigma_x * sigma_x) / (sigma_x + sigma_x);
        sigma_x2 = (2.0f * sigma_x * sigma_x) / (sigma_x + sigma_x);
    }

    real_cpu sigma_y1 = 0.0;
    real_cpu sigma_y2 = 0.0;

    if(sigma_y != 0.0) {
        sigma_y1 = (2.0f * sigma_y * sigma_y) / (sigma_y + sigma_y);
        sigma_y2 = (2.0f * sigma_y * sigma_y) / (sigma_y + sigma_y);
    }

    real_cpu sigma_z1 = 0.0;
    real_cpu sigma_z2 = 0.0;

    if(sigma_z != 0.0) {
        sigma_z1 = (2.0f * sigma_z * sigma_z) / (sigma_z + sigma_z);
        sigma_z2 = (2.0f * sigma_z * sigma_z) / (sigma_z + sigma_z);
    }

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
                dx = black_neighbor_cell->dx;
                dy = black_neighbor_cell->dy;
                dz = black_neighbor_cell->dz;
            } else {
                dx = grid_cell->dx;
                dy = grid_cell->dy;
                dz = grid_cell->dz;
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

                struct element new_element = fill_element(position, direction, dx, dy, dz, sigma_x1, sigma_x2, sigma_y1,
                                                          sigma_y2, sigma_z1, sigma_z2, cell_elements);

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

                struct element new_element = fill_element(position, direction, dx, dy, dz, sigma_x1, sigma_x2, sigma_y1,
                                                          sigma_y2, sigma_z1, sigma_z2, cell_elements);

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
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data.config, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config->config_data.config, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config->config_data.config, "sigma_z");

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

ASSEMBLY_MATRIX(source_sink_discretization_matrix) {

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

    real channel_width = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, channel_width, config->config_data.config, "channel_width");

    real channel_length = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, channel_length, config->config_data.config, "channel_length");

    bool inside;

    real side_length_x = the_grid->side_length_x;
    real side_length_y = the_grid->side_length_y;
    real side_length_z = the_grid->side_length_z;

    real region_height = (side_length_y - channel_width) / 2.0;

    #pragma omp parallel for
    for (i = 0; i < num_active_cells; i++)
    {
        real sigma_x_new;
        real sigma_y_new;
        real sigma_z_new;

        double x = ac[i]->center_x;
        double y = ac[i]->center_y;
        double z = ac[i]->center_z;

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

ASSEMBLY_MATRIX(source_sink_discretization_matrix_with_different_sigma) 
{

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

    real channel_width = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, channel_width, config->config_data.config, "channel_width");

    real channel_length = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, channel_length, config->config_data.config, "channel_length");

    real source_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, source_factor, config->config_data.config, "source_factor");

    real sink_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sink_factor, config->config_data.config, "sink_factor");

    bool inside_1, inside_2, inside_3, inside_4;

    real side_length_x = the_grid->side_length_x;
    real side_length_y = the_grid->side_length_y;
    real side_length_z = the_grid->side_length_z;

    real region_height = (side_length_y - channel_width) / 2.0;

    #pragma omp parallel for
    for (i = 0; i < num_active_cells; i++)
    {
        real sigma_x_new = sigma_x;
        real sigma_y_new = sigma_y;
        real sigma_z_new = sigma_z;

        double x = ac[i]->center_x;
        double y = ac[i]->center_y;
        double z = ac[i]->center_z;

        // Check region 1
        inside_1 = (x >= 0.0) && (x <= channel_length) && (y >= 0.0) && (y <= region_height);
        
        // Check region 2
        inside_2 = (x >= 0.0) && (x <= channel_length) && (y >= region_height + channel_width) && (y <= side_length_y);

        if (inside_1 || inside_2)
        {
            sigma_x_new = 0.0;
            sigma_y_new = 0.0;
            sigma_z_new = 0.0;
        }

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

// TODO: The kappas can be different from one cell to another
ASSEMBLY_MATRIX(ddm_assembly_matrix) 
{

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data.config, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config->config_data.config, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config->config_data.config, "sigma_z");

    real cell_length_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, cell_length_x, config->config_data.config, "cell_length_x");

    real cell_length_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, cell_length_y, config->config_data.config, "cell_length_y");

    real cell_length_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, cell_length_z, config->config_data.config, "cell_length_z");

    // TODO: The kappas can change from one volume to another ...
    // Here we are considering dx, dy, dz equal for all cells over the domain ...
    real dx = ac[0]->dx;
    real dy = ac[0]->dy;
    real dz = ac[0]->dz;

    real kappa_x = 0.0;
    kappa_x = calculate_kappa(cell_length_x,dx);
    the_solver->kappa_x = kappa_x;

    real kappa_y = 0.0;
    kappa_y = calculate_kappa(cell_length_y,dy);
    the_solver->kappa_y = kappa_y;

    real kappa_z = 0.0;
    kappa_z = calculate_kappa(cell_length_z,dz);
    the_solver->kappa_z = kappa_z;

    printf("[!] Using DDM formulation\n");
    printf("[X] Cell length = %.10lf || sigma_x = %.10lf || dx = %.10lf || kappa_x = %.10lf\n",\
            cell_length_x,sigma_x,dx,kappa_x);
    printf("[Y] Cell length = %.10lf || sigma_y = %.10lf || dy = %.10lf || kappa_y = %.10lf\n",\
            cell_length_y,sigma_y,dy,kappa_y);
    printf("[Z] Cell length = %.10lf || sigma_z = %.10lf || dz = %.10lf || kappa_z = %.10lf\n",\
            cell_length_z,sigma_z,dz,kappa_z);

    initialize_diagonal_elements_ddm(the_solver, the_grid,\
                                    dx, dy, dz,\
                                    kappa_x, kappa_y, kappa_z,
                                    sigma_x, sigma_y, sigma_z);


    int i;

    #pragma omp parallel for
    for(i = 0; i < num_active_cells; i++) 
    {

        // Computes and designates the flux due to south cells.
        fill_discretization_matrix_elements_ddm(sigma_x, sigma_y, sigma_z,\
                                                kappa_x,kappa_y,kappa_z,\
                                                the_solver->dt,\
                                                ac[i], ac[i]->south, 's');

        // Computes and designates the flux due to north cells.
        fill_discretization_matrix_elements_ddm(sigma_x, sigma_y, sigma_z,\
                                                kappa_x,kappa_y,kappa_z,\
                                                the_solver->dt,\
                                                ac[i], ac[i]->north, 'n');

        // Computes and designates the flux due to east cells.
        fill_discretization_matrix_elements_ddm(sigma_x, sigma_y, sigma_z,\
                                                kappa_x,kappa_y,kappa_z,\
                                                the_solver->dt,\
                                                ac[i], ac[i]->east, 'e');

        // Computes and designates the flux due to west cells.
        fill_discretization_matrix_elements_ddm(sigma_x, sigma_y, sigma_z,\
                                                kappa_x,kappa_y,kappa_z,\
                                                the_solver->dt,\
                                                ac[i], ac[i]->west, 'w');

        // Computes and designates the flux due to front cells.
        fill_discretization_matrix_elements_ddm(sigma_x, sigma_y, sigma_z,\
                                                kappa_x,kappa_y,kappa_z,\
                                                the_solver->dt,\
                                                ac[i], ac[i]->front, 'f');

        // Computes and designates the flux due to back cells.
        fill_discretization_matrix_elements_ddm(sigma_x, sigma_y, sigma_z,\
                                                kappa_x,kappa_y,kappa_z,\
                                                the_solver->dt,\
                                                ac[i], ac[i]->back, 'b');
    }
    
}

// *******************************************************************************************************
// Berg code
void initialize_diagonal_elements_purkinje (struct monodomain_solver *the_solver, struct grid *the_grid) 
{
    real_cpu alpha;
    real_cpu dx, dy, dz;
    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;
    struct node *n = the_grid->the_purkinje_network->list_nodes;
    real_cpu beta = the_solver->beta;
    real_cpu cm = the_solver->cm;

    real_cpu dt = the_solver->dt;

    int i;

    for (i = 0; i < num_active_cells; i++) 
    {
        dx = ac[i]->dx;
        dy = ac[i]->dy;
        dz = ac[i]->dz;

        alpha = ALPHA(beta, cm, dt, dx, dy, dz);

        struct element element;
        element.column = ac[i]->grid_position;
        element.cell = ac[i];
        element.value = alpha;

        if (ac[i]->elements != NULL) 
        {
            arrfree(ac[i]->elements);
        }

        ac[i]->elements = NULL;
        arrsetcap(ac[i]->elements,n->num_edges);
        arrput(ac[i]->elements, element);

        n = n->next;
    }       
}

// For the Purkinje fibers we only need to solve the 1D Monodomain equation
static void fill_discretization_matrix_elements_purkinje (real_cpu sigma_x, struct cell_node **grid_cells, uint32_t num_active_cells,
                                                        struct node *pk_node) 
{
    
    struct edge *e;
    struct element **cell_elements;
    real_cpu dx;

    real_cpu sigma_x1 = (2.0f * sigma_x * sigma_x) / (sigma_x + sigma_x);

    int i;

    for (i = 0; i < num_active_cells; i++, pk_node = pk_node->next)
    {
        cell_elements = &grid_cells[i]->elements;
        dx = grid_cells[i]->dx;

        e = pk_node->list_edges;

        // Do the mapping of the edges from the graph to the sparse matrix data structure ...
        while (e != NULL)
        {
            struct element new_element;

            // Neighbour elements ...
            new_element.column = e->id;
            new_element.value = -sigma_x1 * dx;
            new_element.cell = grid_cells[e->id];

            // Diagonal element ...
            cell_elements[0]->value += (sigma_x1 * dx);

            arrput(grid_cells[i]->elements,new_element);

            e = e->next;         
        }
    }
}

ASSEMBLY_MATRIX(purkinje_fibers_assembly_matrix) 
{

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;
    struct node *pk_node = the_grid->the_purkinje_network->list_nodes;

    initialize_diagonal_elements_purkinje(the_solver, the_grid);

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data.config, "sigma_x");

    fill_discretization_matrix_elements_purkinje(sigma_x,ac,num_active_cells,pk_node);

    /*
    printf("Sigma_x = %.10lf\n",sigma_x);
    printf("dx = %.10lf\n",ac[0]->dx);
    
    for (int i = 0; i < num_active_cells; i++)
    {
        printf("\nCell %d -- Diagonal = %lf\n",i,ac[i]->elements[0].value);
        int count = arrlen(ac[i]->elements);
        printf("\tElements:\n");
        for (int j = 1; j < count; j++)
            printf("\t%d -- Column = %d -- Value = %lf\n",ac[i]->elements[j].column,ac[i]->elements[j].column,ac[i]->elements[j].value);
    }

    printf("Leaving program ...\n");
    exit(EXIT_FAILURE);
    */

}
// *******************************************************************************************************
