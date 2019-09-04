//
// Created by bergolho on 04/09/19.
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
#include "../libraries_common/common_data_structures.h"


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
    static bool sigma_initialized = false;

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;
    struct node *pk_node = the_grid->the_purkinje_network->list_nodes;

    initialize_diagonal_elements_purkinje(the_solver, the_grid);

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data.config, "sigma_x");

    if(!sigma_initialized) 
    {
        #pragma omp parallel for
        for (uint32_t i = 0; i < num_active_cells; i++) 
        {
            ac[i]->sigma_x = sigma_x;
        }

        sigma_initialized = true;
    }

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