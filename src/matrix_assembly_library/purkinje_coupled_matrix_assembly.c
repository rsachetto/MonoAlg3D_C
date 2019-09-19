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

INIT_ASSEMBLY_MATRIX(set_initial_conditions_coupled_fvm) {

    real_cpu alpha;
    
    // Tissue parameters
    struct cell_node **ac = the_grid->active_cells;
    uint32_t active_cells = the_grid->num_active_cells;

    // Purkinje parameters
    struct cell_node **ac_purkinje = the_grid->the_purkinje->purkinje_cells;
    uint32_t active_purkinje_cells = the_grid->the_purkinje->num_active_purkinje_cells;

    // Common parameters
    real_cpu beta = the_solver->beta;
    real_cpu cm = the_solver->cm;
    real_cpu dt = the_solver->dt;
    int i;

    // Tissue section
    #pragma omp parallel for private(alpha)
    for(i = 0; i < active_cells; i++) 
    {

        alpha = ALPHA(beta, cm, dt, ac[i]->discretization.x, ac[i]->discretization.y, ac[i]->discretization.z);
        ac[i]->v = initial_v;
        ac[i]->b = initial_v * alpha;
    }

    // Purkinje section
    #pragma omp parallel for private(alpha)
    for(i = 0; i < active_purkinje_cells; i++) 
    {

        alpha = ALPHA(beta, cm, dt, ac[i]->discretization.x, ac[i]->discretization.y, ac[i]->discretization.z);
        ac_purkinje[i]->v = purkinje_initial_v;
        ac_purkinje[i]->b = purkinje_initial_v * alpha;
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

void initialize_diagonal_elements_purkinje (struct monodomain_solver *the_solver, struct grid *the_grid) 
{
    real_cpu alpha;
    real_cpu dx, dy, dz;

    uint32_t num_active_cells = the_grid->the_purkinje->num_active_purkinje_cells;
    struct cell_node **ac = the_grid->the_purkinje->purkinje_cells;

    struct node *n = the_grid->the_purkinje->the_network->list_nodes;
    
    real_cpu beta = the_solver->beta;
    real_cpu cm = the_solver->cm;
    real_cpu dt = the_solver->dt;

    int i;

    for (i = 0; i < num_active_cells; i++) 
    {
        dx = ac[i]->discretization.x;
        dy = ac[i]->discretization.y;
        dz = ac[i]->discretization.z;

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
        dx = grid_cells[i]->discretization.x;

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

ASSEMBLY_MATRIX (purkinje_coupled_endocardium_assembly_matrix)
{

    // Endocardium section
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

    real sigma_purkinje = sigma_x;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real,sigma_purkinje,config->config_data,"sigma_purkinje");

    if(!sigma_initialized) 
    {
        #pragma omp parallel for
        for (i = 0; i < num_active_cells; i++) 
        {
            ac[i]->sigma.x = sigma_x;
            ac[i]->sigma.y = sigma_y;
            ac[i]->sigma.z = sigma_z;
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

    // Purkinje section
    static bool sigma_purkinje_initialized = false;

    uint32_t num_purkinje_active_cells = the_grid->the_purkinje->num_active_purkinje_cells;
    struct cell_node **ac_purkinje = the_grid->the_purkinje->purkinje_cells;

    struct node *pk_node = the_grid->the_purkinje->the_network->list_nodes;

    initialize_diagonal_elements_purkinje(the_solver, the_grid);
    
    if(!sigma_purkinje_initialized) 
    {
        #pragma omp parallel for
        for (uint32_t i = 0; i < num_purkinje_active_cells; i++) 
        {
            ac_purkinje[i]->sigma.x = sigma_purkinje;
        }

        sigma_purkinje_initialized = true;
    }

    fill_discretization_matrix_elements_purkinje(sigma_purkinje,ac_purkinje,num_purkinje_active_cells,pk_node);

    // DEBUG
    // Endocardium cells
    /*
    for (uint32_t i = 0; i < num_active_cells; i++)
    {
        struct element *cell_elements = ac[i]->elements;
        size_t max_elements = arrlen(cell_elements);

        printf("Line %u\n",i);
        printf("\tColumn = %u -- Value = %g\n",cell_elements[0].column,cell_elements[0].value);
        for (uint32_t j = 1; j < max_elements; j++)
            printf("\tColumn = %u -- Value = %g\n",cell_elements[j].column,cell_elements[j].value);
    }

    // Purkinje cells
    for (uint32_t i = 0; i < num_purkinje_active_cells; i++)
    {
        struct element *cell_elements = ac_purkinje[i]->elements;
        size_t max_elements = arrlen(cell_elements);

        printf("Line %u\n",i);
        printf("\tColumn = %u -- Value = %g\n",cell_elements[0].column,cell_elements[0].value);
        for (uint32_t j = 1; j < max_elements; j++)
            printf("\tColumn = %u -- Value = %g\n",cell_elements[j].column,cell_elements[j].value);
    }
    */
    
}