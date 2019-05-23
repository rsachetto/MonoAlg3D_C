//
// Created by sachetto on 29/09/17.
//

#include <assert.h>
#include <inttypes.h>
#include <float.h>

#include "grid.h"

#include "../../single_file_libraries/stb_ds.h"

struct grid *new_grid() {
    struct grid *result = (struct grid *)malloc(sizeof(struct grid));
    result->first_cell = NULL;
    result->active_cells = NULL;
    result->adaptive = false;

    result->refined_this_step = NULL;
    result->free_sv_positions = NULL;

    arrsetcap(result->refined_this_step, 128);
    arrsetcap(result->free_sv_positions, 128);

    // Purkinje
    result->the_purkinje_network = new_graph();

    return result;
}

void initialize_and_construct_grid(struct grid *the_grid, float side_length_x, float side_length_y,
                                   float side_length_z) {
    assert(the_grid);

    initialize_grid(the_grid, side_length_x, side_length_y, side_length_z);
    construct_grid(the_grid);
}

void initialize_grid(struct grid *the_grid, float side_length_x, float side_length_y, float side_length_z) {

    assert(the_grid);

    the_grid->side_length_x = side_length_x;
    the_grid->side_length_y = side_length_y;
    the_grid->side_length_z = side_length_z;
    the_grid->number_of_cells = 0;
}

void construct_grid(struct grid *the_grid) {

    assert(the_grid);

    float side_length_x = the_grid->side_length_x;
    float side_length_y = the_grid->side_length_y;
    float side_length_z = the_grid->side_length_z;

    // Cell nodes.
    struct cell_node *front_northeast_cell, *front_northwest_cell, *front_southeast_cell, *front_southwest_cell,
        *back_northeast_cell, *back_northwest_cell, *back_southeast_cell, *back_southwest_cell;

    front_northeast_cell = new_cell_node();
    front_northwest_cell = new_cell_node();
    front_southeast_cell = new_cell_node();
    front_southwest_cell = new_cell_node();
    back_northeast_cell = new_cell_node();
    back_northwest_cell = new_cell_node();
    back_southeast_cell = new_cell_node();
    back_southwest_cell = new_cell_node();

    // Transition nodes.
    struct transition_node *north_transition_node, *south_transition_node, *east_transition_node, *west_transition_node,
        *front_transition_node, *back_transition_node;

    north_transition_node = new_transition_node();
    south_transition_node = new_transition_node();
    east_transition_node = new_transition_node();
    west_transition_node = new_transition_node();
    front_transition_node = new_transition_node();
    back_transition_node = new_transition_node();

    real_cpu half_side_length_x = side_length_x / 2.0f;
    real_cpu half_side_length_y = side_length_y / 2.0f;
    real_cpu half_side_length_z = side_length_z / 2.0f;

    real_cpu quarter_side_length_x = half_side_length_x / 2.0f;
    real_cpu quarter_side_length_y = half_side_length_y / 2.0f;
    real_cpu quarter_side_length_z = half_side_length_z / 2.0f;

    //__________________________________________________________________________
    //              Initialization of transition nodes.
    //__________________________________________________________________________
    // East transition node.
    set_transition_node_data(east_transition_node, 1, 'e', NULL, front_southeast_cell, back_southeast_cell,
                             back_northeast_cell, front_northeast_cell);

    // North transition node.
    set_transition_node_data(north_transition_node, 1, 'n', NULL, front_northwest_cell, front_northeast_cell,
                             back_northeast_cell, back_northwest_cell);

    // West transition node.
    set_transition_node_data(west_transition_node, 1, 'w', NULL, front_southwest_cell, back_southwest_cell,
                             back_northwest_cell, front_northwest_cell);

    // South transition node.
    set_transition_node_data(south_transition_node, 1, 's', NULL, front_southwest_cell, front_southeast_cell,
                             back_southeast_cell, back_southwest_cell);

    // Front transition node.
    set_transition_node_data(front_transition_node, 1, 'f', NULL, front_southwest_cell, front_southeast_cell,
                             front_northeast_cell, front_northwest_cell);

    // Back transition node.
    set_transition_node_data(back_transition_node, 1, 'b', NULL, back_southwest_cell, back_southeast_cell,
                             back_northeast_cell, back_northwest_cell);

    //__________________________________________________________________________
    //                      Initialization of cell nodes.
    //__________________________________________________________________________
    // front Northeast subcell initialization.
    set_cell_node_data(front_northeast_cell, half_side_length_x, half_side_length_y, half_side_length_z, 1,
                       east_transition_node, north_transition_node, front_northwest_cell, front_southeast_cell,
                       front_transition_node, back_northeast_cell, NULL, back_northeast_cell, 0, 1,
                       half_side_length_x + quarter_side_length_x, half_side_length_y + quarter_side_length_y,
                       half_side_length_z + quarter_side_length_z);

    // back Northeast subcell initialization.
    set_cell_node_data(back_northeast_cell, half_side_length_x, half_side_length_y, half_side_length_z, 2,
                       east_transition_node, north_transition_node, back_northwest_cell, back_southeast_cell,
                       front_northeast_cell, back_transition_node, front_northeast_cell, back_northwest_cell, 1, 2,
                       quarter_side_length_x, half_side_length_y + quarter_side_length_y,
                       half_side_length_z + quarter_side_length_z);

    // back Northwest subcell initialization.
    set_cell_node_data(back_northwest_cell, half_side_length_x, half_side_length_y, half_side_length_z, 3,
                       back_northeast_cell, north_transition_node, west_transition_node, back_southwest_cell,
                       front_northwest_cell, back_transition_node, back_northeast_cell, front_northwest_cell, 2, 2,
                       quarter_side_length_x, quarter_side_length_y, half_side_length_z + quarter_side_length_z);

    // front Northwest subcell initialization.
    set_cell_node_data(front_northwest_cell, half_side_length_x, half_side_length_y, half_side_length_z, 4,
                       front_northeast_cell, north_transition_node, west_transition_node, front_southwest_cell,
                       front_transition_node, back_northwest_cell, back_northwest_cell, front_southwest_cell, 3, 3,
                       half_side_length_x + quarter_side_length_x, quarter_side_length_y,
                       half_side_length_z + quarter_side_length_z);

    // front Southwest subcell initialization.
    set_cell_node_data(front_southwest_cell, half_side_length_x, half_side_length_y, half_side_length_z, 5,
                       front_southeast_cell, front_northwest_cell, west_transition_node, south_transition_node,
                       front_transition_node, back_southwest_cell, front_northwest_cell, back_southwest_cell, 4, 3,
                       half_side_length_x + quarter_side_length_x, quarter_side_length_y, quarter_side_length_z);

    // back Southwest subcell initialization.
    set_cell_node_data(back_southwest_cell, half_side_length_x, half_side_length_y, half_side_length_z, 6,
                       back_southeast_cell, back_northwest_cell, west_transition_node, south_transition_node,
                       front_southwest_cell, back_transition_node, front_southwest_cell, back_southeast_cell, 5, 4,
                       quarter_side_length_x, quarter_side_length_y, quarter_side_length_z);

    // back Southeast subcell initialization.
    set_cell_node_data(back_southeast_cell, half_side_length_x, half_side_length_y, half_side_length_z, 7,
                       east_transition_node, back_northeast_cell, back_southwest_cell, south_transition_node,
                       front_southeast_cell, back_transition_node, back_southwest_cell, front_southeast_cell, 6, 4,
                       quarter_side_length_x, half_side_length_y + quarter_side_length_y, quarter_side_length_z);

    // front Southeast subcell initialization.
    set_cell_node_data(front_southeast_cell, half_side_length_x, half_side_length_y, half_side_length_z, 8,
                       east_transition_node, front_northeast_cell, front_southwest_cell, south_transition_node,
                       front_transition_node, back_southeast_cell, back_southeast_cell, NULL, 7, 5,
                       half_side_length_x + quarter_side_length_x, half_side_length_y + quarter_side_length_y,
                       quarter_side_length_z);

    // Grid initialization
    the_grid->first_cell = front_northeast_cell;
    the_grid->number_of_cells = 8;
}

void print_grid(struct grid *the_grid, FILE *output_file) {

    struct cell_node *grid_cell = the_grid->first_cell;

    real_cpu center_x, center_y, center_z, dx, dy, dz;
    real_cpu v;

    while(grid_cell != 0) {

        if(grid_cell->active) {

            center_x = grid_cell->center_x;
            center_y = grid_cell->center_y;
            center_z = grid_cell->center_z;

            dx = grid_cell->dx;
            dy = grid_cell->dy;
            dz = grid_cell->dz;

            v = grid_cell->v;

            fprintf(output_file, "%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", center_x, center_y, center_z, dx, dy, dz, v);
        }
        grid_cell = grid_cell->next;
    }
}

void order_grid_cells(struct grid *the_grid) {

    struct cell_node *grid_cell;
    grid_cell = the_grid->first_cell;

    // Here we allocate the maximum number of cells we will need for the whole simulation
    if(the_grid->active_cells == NULL) {
        the_grid->active_cells = (struct cell_node **)malloc(sizeof(struct cell_node *) * the_grid->number_of_cells);
    }

    uint32_t counter = 0;
    while(grid_cell != 0) {
        if(grid_cell->active) {
            grid_cell->grid_position = counter;
            the_grid->active_cells[counter] = grid_cell;
            counter++;
        }

        grid_cell = grid_cell->next;
    }

    the_grid->num_active_cells = counter;
}

void clean_grid(struct grid *the_grid) {

    assert(the_grid);

    struct cell_node *grid_cell = NULL;

    // TODO: Think about this function when the coupling happens ...
    // Delete nodes from the Purkinje network
    if (the_grid->the_purkinje_network->list_nodes != NULL) 
    {

        // First free the Purkinje mesh structure
        free_graph(the_grid->the_purkinje_network);

        grid_cell = the_grid->first_cell;

        // Then, delete the cells from the Purkinje network 
        if(grid_cell) 
        {
            while (grid_cell) 
            {

                struct cell_node *next = grid_cell->next;
                free_cell_node(grid_cell);
                grid_cell = next;

            }
        }
    }
    // Delete the tissue cells
    else
    {

        // In order to release the memory allocated for the grid, the grid is
        // derefined to level 1. Thus, the grid shape is known and each node can
        // be easily reached.

        uint32_t number_of_cells = the_grid->number_of_cells;
        while(number_of_cells > 8)
        {
            derefine_all_grid(the_grid);
            number_of_cells = the_grid->number_of_cells;
        }

        grid_cell = the_grid->first_cell;

        if(grid_cell)
        {

            // Deleting transition nodes.
            free((struct transition_node *)(grid_cell->north));
            free((struct transition_node *)(grid_cell->front));
            free((struct transition_node *)(grid_cell->east));
            free((struct transition_node *)(((struct cell_node *)(grid_cell->west))->west));
            free((struct transition_node *)(((struct cell_node *)(grid_cell->south))->south));
            free((struct transition_node *)(((struct cell_node *)(grid_cell->back))->back));

            // Deleting cells nodes.
            while(grid_cell)
            {
                struct cell_node *next = grid_cell->next;
                free_cell_node(grid_cell);
                grid_cell = next;
            }
        }
    }

    if(the_grid->refined_this_step) 
    {
        arrreset(the_grid->refined_this_step);
    }

    if(the_grid->free_sv_positions) 
    {
        arrreset(the_grid->free_sv_positions);
    }
}

void clean_and_free_grid(struct grid *the_grid) {

    assert(the_grid);

    clean_grid(the_grid);

    if(the_grid->active_cells) 
    {
        free(the_grid->active_cells);
    }

    arrfree(the_grid->refined_this_step);
    arrfree(the_grid->free_sv_positions);

    free(the_grid);
}

// Prints grid discretization matrix.
void print_grid_matrix(struct grid *the_grid, FILE *output_file) {

    assert(the_grid);
    assert(output_file);

    struct cell_node *grid_cell;
    grid_cell = the_grid->first_cell;
    struct element element;
    element_array cell_elements;

    while(grid_cell != 0) {
        if(grid_cell->active) {

            cell_elements = grid_cell->elements;
            size_t max_el = arrlen(cell_elements);

            for(size_t i = 0; i < max_el; i++) {

                element = cell_elements[i];
                if(element.cell != NULL) {
                    fprintf(output_file,
                            "%" PRIu32 " "
                            "%" PRIu32 " %.15lf\n",
                            grid_cell->grid_position + 1, (element.column) + 1, element.value);
                } else {
                    break;
                }
            }
        }
        grid_cell = grid_cell->next;
    }
}

int compare_elements (const void * a, const void * b)
{
    if ( (*(struct element*)a).column <  (*(struct element*)b).column ) return -1;
    if ( (*(struct element*)a).column == (*(struct element*)b).column ) return 0;
    if ( (*(struct element*)a).column >  (*(struct element*)b).column ) return 1;

    return 1;//Unreachable
}

void print_grid_matrix_as_octave_matrix(struct grid *the_grid, FILE *output_file) {

    assert(the_grid);
    assert(output_file);

    struct cell_node *grid_cell;
    grid_cell = the_grid->first_cell;
    struct element element;
    element_array cell_elements;

    fprintf(output_file, "# Created by Monodomain solver\n");
    fprintf(output_file, "# name: Alg_grid_matrix\n");
    fprintf(output_file, "# type: sparse matrix\n");
    fprintf(output_file, "# nnz:                                    \n");
    fprintf(output_file, "# rows: %d\n", the_grid->num_active_cells);
    fprintf(output_file, "# columns: %d\n", the_grid->num_active_cells);

    int nnz = 0;

    while(grid_cell != 0) {
        if(grid_cell->active) {

            cell_elements = grid_cell->elements;
            size_t max_el = arrlen(cell_elements);

            qsort (cell_elements, max_el, sizeof(struct element), compare_elements);

            for(size_t i = 0; i < max_el; i++) {

                element = cell_elements[i];
                if(element.cell != NULL) {
                    nnz += 1;
                    fprintf(output_file,
                            "%" PRIu32 " "
                            "%" PRIu32 " %.15lf\n",
                            grid_cell->grid_position + 1, (element.column) + 1, element.value);
                } else {
                    break;
                }
            }
        }
        grid_cell = grid_cell->next;
    }

    fseek(output_file, 84, SEEK_SET);
    fprintf(output_file, "%d", nnz);
}

void print_grid_vector(struct grid *the_grid, FILE *output_file, char name) {
    struct cell_node *grid_cell;
    grid_cell = the_grid->first_cell;

    while(grid_cell != 0) {
        if(grid_cell->active) {
            if(name == 'b')
                fprintf(output_file, "%.15lf\n", grid_cell->b);
            else if(name == 'x')
                fprintf(output_file, "%.15lf\n", grid_cell->v);
        }
        grid_cell = grid_cell->next;
    }
}

real_cpu *grid_vector_to_array(struct grid *the_grid, char name, uint32_t *num_lines) {
    struct cell_node *grid_cell;
    grid_cell = the_grid->first_cell;

    *num_lines = the_grid->num_active_cells;
    real_cpu *vector = (real_cpu *)malloc(*num_lines * sizeof(real_cpu));

    while(grid_cell != 0) {
        if(grid_cell->active) {
            if(name == 'b')
                vector[grid_cell->grid_position] = grid_cell->b;
            else if(name == 'x')
                vector[grid_cell->grid_position] = grid_cell->v;
        }
        grid_cell = grid_cell->next;
    }

    return vector;
}

void save_grid_domain(struct grid *the_grid, const char *file_name) {
    struct cell_node *grid_cell = the_grid->first_cell;
    FILE *f = fopen(file_name, "w");

    while(grid_cell != 0) {
        if(grid_cell->active) {
            fprintf(f, "%lf,%lf,%lf,%lf,%lf,%lf\n", grid_cell->center_x, grid_cell->center_y, grid_cell->center_z,
                    grid_cell->dx, grid_cell->dy, grid_cell->dz);
        }
        grid_cell = grid_cell->next;
    }
    fclose(f);
}

int get_num_refinement_steps_to_discretization(float side_len, real_cpu h) {

    int num_steps = 0;
    real_cpu aux = side_len;

    while(aux > h) {
        num_steps++;
        aux /= 2.0;
    }

    return num_steps - 1;
}

void initialize_grid_purkinje (struct grid *the_grid)
{
    assert(the_grid);

    the_grid->number_of_cells = 0;
}

void construct_grid_purkinje (struct grid *the_grid)
{
    assert(the_grid);

    // TODO: Allow dx, dy, dz to be different in the Purkinje code
    real_cpu side_length_x = the_grid->the_purkinje_network->dx;
    real_cpu side_length_y = the_grid->the_purkinje_network->dx;
    real_cpu side_length_z = the_grid->the_purkinje_network->dx;

    real_cpu half_side_length_x = side_length_x / 2.0f;
    real_cpu half_side_length_y = side_length_y / 2.0f;
    real_cpu half_side_length_z = side_length_z / 2.0f;

//    real_cpu quarter_side_length_x = half_side_length_x / 2.0f;
//    real_cpu quarter_side_length_y = half_side_length_y / 2.0f;
//    real_cpu quarter_side_length_z = half_side_length_z / 2.0f;

    int total_nodes = the_grid->the_purkinje_network->total_nodes;
    
    // Create an array of cell nodes
    struct cell_node **cells = (struct cell_node**)malloc(sizeof(struct cell_node*)*total_nodes);
    for (int i = 0; i < total_nodes; i++)
        cells[i] = new_cell_node();
    
    // Pass through the Purkinje graph and set the cell nodes.
    struct node *n = the_grid->the_purkinje_network->list_nodes;
    for (int i = 0; i < total_nodes; i++)
    {
        
        if (i == 0)
            set_cell_node_data (cells[i],half_side_length_x,half_side_length_y,half_side_length_z,\
                            0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,cells[i+1],i,0,\
                            n->x,n->y,n->z);
        else if (i == total_nodes-1)
            set_cell_node_data (cells[i],half_side_length_x,half_side_length_y,half_side_length_z,\
                        0,NULL,NULL,NULL,NULL,NULL,NULL,\
                        cells[i-1],NULL,i,0,\
                        n->x,n->y,n->z);
        else
            set_cell_node_data (cells[i],half_side_length_x,half_side_length_y,half_side_length_z,\
                        0,NULL,NULL,NULL,NULL,NULL,NULL,\
                        cells[i-1],cells[i+1],i,0,\
                        n->x,n->y,n->z);

        // Do not refine the Purkinje cells !
        cells[i]->can_change = false;

        n = n->next;
    }
    
    // Grid initialization
    the_grid->first_cell = cells[0];
    the_grid->number_of_cells = total_nodes;
    
}

void initialize_and_construct_grid_purkinje (struct grid *the_grid)
{
    assert(the_grid);

    initialize_grid_purkinje(the_grid);
    construct_grid_purkinje(the_grid);
}

void translate_mesh_to_origin(struct grid *grid) {

    real_cpu minx = FLT_MAX;
    real_cpu miny = FLT_MAX;
    real_cpu minz = FLT_MAX;

    struct cell_node *grid_cell;

    float center_x;
    float center_y;
    float center_z;

    grid_cell = grid->first_cell;

    while(grid_cell != 0) {

        center_x = grid_cell->center_x;
        center_y = grid_cell->center_y;
        center_z = grid_cell->center_z;

        if(center_x < minx){
            minx = center_x;
        }

        if(center_y < miny){
            miny = center_y;
        }

        if(center_z < minz){
            minz = center_z;
        }

        grid_cell = grid_cell->next;
    }

    grid_cell = grid->first_cell;

    while(grid_cell != 0) {

        grid_cell->center_x = grid_cell->center_x - minx + (grid_cell->dx/2.0f);
        grid_cell->center_y = grid_cell->center_y - miny + (grid_cell->dy/2.0f);
        grid_cell->center_z = grid_cell->center_z - minz + (grid_cell->dz/2.0f);

        grid_cell = grid_cell->next;
    }

}