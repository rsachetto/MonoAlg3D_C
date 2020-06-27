//
// Created by bergolho on 15/06/20.
//

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>

#include "../alg/grid/grid.h"
#include "../config/assembly_matrix_config.h"
#include "../monodomain/constants.h"
#include "../utils/utils.h"
#include "../3dparty/stb_ds.h"
#include "../libraries_common/common_data_structures.h"

#include "../config_helpers/config_helpers.h"

INIT_ASSEMBLY_MATRIX(set_initial_conditions_fvm) {

    real_cpu alpha;
    struct cell_node **ac = the_grid->purkinje->purkinje_cells;
    uint32_t active_cells = the_grid->purkinje->num_active_purkinje_cells;
    real_cpu beta = the_solver->beta;
    real_cpu cm = the_solver->cm;
    real_cpu dt = the_solver->dt;
    uint32_t i;

    OMP(parallel for private(alpha))
    for(i = 0; i < active_cells; i++) 
    {

        alpha = ALPHA(beta, cm, dt, ac[i]->discretization.x, ac[i]->discretization.y, ac[i]->discretization.z);
        ac[i]->v = purkinje_initial_v;
        ac[i]->b = purkinje_initial_v * alpha;
    }
}

void initialize_diagonal_elements_purkinje (struct monodomain_solver *the_solver, struct grid *the_grid) 
{
    real_cpu alpha;
    real_cpu dx, dy, dz;
    uint32_t num_active_cells = the_grid->purkinje->num_active_purkinje_cells;
    struct cell_node **ac = the_grid->purkinje->purkinje_cells;
    struct node *n = the_grid->purkinje->network->list_nodes;
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
static void fill_discretization_matrix_elements_purkinje (real_cpu sigma_x, bool has_point_data, struct cell_node **grid_cells, uint32_t num_active_cells,
                                                        struct node *pk_node) 
{
    
    struct edge *e;
    struct element **cell_elements;

    real_cpu dx, dy, dz;
    real_cpu multiplier;

    int i;

    for (i = 0; i < num_active_cells; i++, pk_node = pk_node->next)
    {
        cell_elements = &grid_cells[i]->elements;
        dx = grid_cells[i]->discretization.x;
        dy = grid_cells[i]->discretization.x;
        dz = grid_cells[i]->discretization.x;

        multiplier = ((dy * dz) / dx);

        e = pk_node->list_edges;

        // Do the mapping of the edges from the graph to the sparse matrix data structure ...
        while (e != NULL)
        {
            struct element new_element;

            // Calculate the conductivity between the two neighboring cells
            if (has_point_data)
            {
                real_cpu sigma_x1 = pk_node->sigma;
                real_cpu sigma_x2 = e->dest->sigma;
                
                if(sigma_x1 != 0.0 && sigma_x2 != 0.0) 
                    sigma_x = (2.0f * sigma_x1 * sigma_x2) / (sigma_x1 + sigma_x2);
            }

            // Neighbour elements ...
            new_element.column = e->id;
            new_element.value = (-sigma_x * multiplier);
            new_element.cell = grid_cells[e->id];

            // Diagonal element ...
            cell_elements[0]->value += (sigma_x * multiplier);

            arrput(grid_cells[i]->elements,new_element);

            e = e->next;         
        }

    }
}

ASSEMBLY_MATRIX(purkinje_fibers_assembly_matrix) 
{
    static bool sigma_initialized = false;

    uint32_t num_active_cells = the_grid->purkinje->num_active_purkinje_cells;
    struct cell_node **ac = the_grid->purkinje->purkinje_cells;
    struct node *pk_node = the_grid->purkinje->network->list_nodes;
    bool has_point_data = the_grid->purkinje->network->has_point_data;

    initialize_diagonal_elements_purkinje(the_solver, the_grid);

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real,sigma_x, config->config_data, "sigma_purkinje");

    if(!sigma_initialized) 
    {
        // Check if the Purkinje network file has the POINT_DATA section
        if (has_point_data)
        {
            struct node *tmp = the_grid->purkinje->network->list_nodes;
            uint32_t i = 0;
            while (tmp != NULL)
            {
                // Copy the prescribed conductivity from the Purkinje network file into the ALG cell structure
                ac[i]->sigma.x = tmp->sigma;

                i++;
                tmp = tmp->next;
            }
        }
        // Otherwise, initilize the conductivity of all cells homogenously with the value from the configuration file
        else
        {
            OMP(parallel for)
            for (uint32_t i = 0; i < num_active_cells; i++) 
            {
                ac[i]->sigma.x = sigma_x;
            }
        }
        
        sigma_initialized = true;
    }

    fill_discretization_matrix_elements_purkinje(sigma_x,has_point_data,ac,num_active_cells,pk_node);

}
