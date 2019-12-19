//
// Created by sachetto on 29/09/17.
//

#include <assert.h>
#include <inttypes.h>
#include <float.h>

#include "grid.h"

#include "../../single_file_libraries/stb_ds.h"
#include "../../utils/file_utils.h"

struct grid *new_grid() {
    struct grid *result = (struct grid *)malloc(sizeof(struct grid));
    result->first_cell = NULL;
    result->active_cells = NULL;
    result->adaptive = false;

    result->refined_this_step = NULL;
    result->free_sv_positions = NULL;
    result->num_active_cells = result->number_of_cells = 0;

    arrsetcap(result->refined_this_step, 128);
    arrsetcap(result->free_sv_positions, 128);

    // Purkinje
    result->the_purkinje = new_grid_purkinje();

    return result;
}

void initialize_and_construct_grid(struct grid *the_grid, struct point_3d side_length) {
    assert(the_grid);

    initialize_grid(the_grid, side_length);
    construct_grid(the_grid);
}

void initialize_grid(struct grid *the_grid, struct point_3d side_length) {

    assert(the_grid);
    the_grid->cube_side_length = side_length;
    the_grid->number_of_cells = 0;
}

void construct_grid(struct grid *the_grid) {

    assert(the_grid);

    struct point_3d side_length = the_grid->cube_side_length;

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


    struct point_3d half_side_length = POINT3D(side_length.x / 2.0f, side_length.y / 2.0,side_length.z / 2.0f);
    struct point_3d quarter_side_length = POINT3D(half_side_length.x / 2.0f, half_side_length.y / 2.0,half_side_length.z / 2.0f);

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
    set_cell_node_data(front_northeast_cell, half_side_length, 1,
                       east_transition_node, north_transition_node, front_northwest_cell, front_southeast_cell,
                       front_transition_node, back_northeast_cell, NULL, back_northeast_cell, 0, 1,
                       POINT3D(half_side_length.x + quarter_side_length.x, half_side_length.y + quarter_side_length.y,
                               half_side_length.z + quarter_side_length.z), ZERO_POINT3D);

    // back Northeast subcell initialization.
    set_cell_node_data(back_northeast_cell, half_side_length, 2,
                       east_transition_node, north_transition_node, back_northwest_cell, back_southeast_cell,
                       front_northeast_cell, back_transition_node, front_northeast_cell, back_northwest_cell, 1, 2,
                       POINT3D(quarter_side_length.x, half_side_length.y + quarter_side_length.y,
                       half_side_length.z + quarter_side_length.z), ZERO_POINT3D);

    // back Northwest subcell initialization.
    set_cell_node_data(back_northwest_cell, half_side_length, 3,
                       back_northeast_cell, north_transition_node, west_transition_node, back_southwest_cell,
                       front_northwest_cell, back_transition_node, back_northeast_cell, front_northwest_cell, 2, 2,
                       POINT3D(quarter_side_length.x, quarter_side_length.y, half_side_length.z + quarter_side_length.z),
                       ZERO_POINT3D);

    // front Northwest subcell initialization.
    set_cell_node_data(front_northwest_cell, half_side_length, 4,
                       front_northeast_cell, north_transition_node, west_transition_node, front_southwest_cell,
                       front_transition_node, back_northwest_cell, back_northwest_cell, front_southwest_cell, 3, 3,
                       POINT3D(half_side_length.x + quarter_side_length.x, quarter_side_length.y,
                       half_side_length.z + quarter_side_length.z), ZERO_POINT3D);

    // front Southwest subcell initialization.
    set_cell_node_data(front_southwest_cell, half_side_length, 5,
                       front_southeast_cell, front_northwest_cell, west_transition_node, south_transition_node,
                       front_transition_node, back_southwest_cell, front_northwest_cell, back_southwest_cell, 4, 3,
                       POINT3D(half_side_length.x + quarter_side_length.x, quarter_side_length.y, quarter_side_length.z),
                       ZERO_POINT3D);

    // back Southwest subcell initialization.
    set_cell_node_data(back_southwest_cell, half_side_length, 6,
                       back_southeast_cell, back_northwest_cell, west_transition_node, south_transition_node,
                       front_southwest_cell, back_transition_node, front_southwest_cell, back_southeast_cell, 5, 4,
                       POINT3D(quarter_side_length.x, quarter_side_length.y, quarter_side_length.z),
                       ZERO_POINT3D);

    // back Southeast subcell initialization.
    set_cell_node_data(back_southeast_cell, half_side_length, 7,
                       east_transition_node, back_northeast_cell, back_southwest_cell, south_transition_node,
                       front_southeast_cell, back_transition_node, back_southwest_cell, front_southeast_cell, 6, 4,
                       POINT3D(quarter_side_length.x, half_side_length.y + quarter_side_length.y, quarter_side_length.z),
                       ZERO_POINT3D);

    // front Southeast subcell initialization.
    set_cell_node_data(front_southeast_cell, half_side_length, 8,
                       east_transition_node, front_northeast_cell, front_southwest_cell, south_transition_node,
                       front_transition_node, back_southeast_cell, back_southeast_cell, NULL, 7, 5,
                       POINT3D(half_side_length.x + quarter_side_length.x, half_side_length.y + quarter_side_length.y,
                       quarter_side_length.z), ZERO_POINT3D);

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

            center_x = grid_cell->center.x;
            center_y = grid_cell->center.y;
            center_z = grid_cell->center.z;

            dx = grid_cell->discretization.x;
            dy = grid_cell->discretization.y;
            dz = grid_cell->discretization.z;

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

    // TODO: Think about this function when the coupling happens ..
    // Delete nodes from the Purkinje network
    if (the_grid->the_purkinje->the_network->list_nodes != NULL) 
    {

        // First free the Purkinje mesh structure
        free_graph(the_grid->the_purkinje->the_network);

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

    if(the_grid->refined_this_step) {
        arrsetlen(the_grid->refined_this_step, 0);
    }

    if(the_grid->free_sv_positions) {
        arrsetlen(the_grid->free_sv_positions, 0);
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

    //free(the_grid->the_purkinje_network); // TODO: Check for leaks with Valgrind
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
            fprintf(f, "%lf,%lf,%lf,%lf,%lf,%lf\n", grid_cell->center.x, grid_cell->center.y, grid_cell->center.z,
                    grid_cell->discretization.x, grid_cell->discretization.y, grid_cell->discretization.z);
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
    assert(the_grid->the_purkinje);

    the_grid->the_purkinje->number_of_purkinje_cells = 0;

}

void construct_grid_purkinje (struct grid *the_grid)
{
    assert(the_grid);
    assert(the_grid->the_purkinje);

    struct grid_purkinje *the_purkinje = the_grid->the_purkinje;

    // TODO: Allow dx, dy, dz to be different in the Purkinje code
    real_cpu side_length_x = the_purkinje->the_network->dx;
    real_cpu side_length_y = the_purkinje->the_network->dx;
    real_cpu side_length_z = the_purkinje->the_network->dx;

    struct point_3d side_length = POINT3D(side_length_x, side_length_y, side_length_z);
    struct point_3d half_side_length = POINT3D(side_length_x / 2.0f, side_length_y / 2.0f, side_length_z / 2.0f);

    uint32_t total_purkinje_nodes = the_purkinje->the_network->total_nodes;
    
    // Create the Purkinje cells array
    struct cell_node **purkinje_cells = the_purkinje->purkinje_cells;
    purkinje_cells = (struct cell_node**)malloc(sizeof(struct cell_node*)*total_purkinje_nodes);
    for (int i = 0; i < total_purkinje_nodes; i++)
        purkinje_cells[i] = new_cell_node();
    
    // Pass through the Purkinje graph and set the cell nodes.
    struct node *n = the_purkinje->the_network->list_nodes;
    for (int i = 0; i < total_purkinje_nodes; i++)
    {
        
        if (i == 0)
            set_cell_node_data (purkinje_cells[i],side_length, 0, NULL,NULL,NULL,NULL,NULL,NULL,NULL,purkinje_cells[i+1],i,0,\
                            POINT3D(n->x,n->y,n->z), ZERO_POINT3D);
        else if (i == total_purkinje_nodes-1)
            set_cell_node_data (purkinje_cells[i],side_length,\
                        0,NULL,NULL,NULL,NULL,NULL,NULL,\
                        purkinje_cells[i-1],NULL,i,0,\
                        POINT3D(n->x,n->y,n->z), ZERO_POINT3D);
        else
            set_cell_node_data (purkinje_cells[i],side_length,\
                        0,NULL,NULL,NULL,NULL,NULL,NULL,\
                        purkinje_cells[i-1],purkinje_cells[i+1],i,0,\
                        POINT3D(n->x,n->y,n->z), ZERO_POINT3D);

        // Do not refine the Purkinje cells !
        purkinje_cells[i]->can_change = false;

        // Set the cell as active
        purkinje_cells[i]->active = true;

        n = n->next;
    }
    
    // Purkinje initialization
    the_purkinje->first_cell = purkinje_cells[0];
    the_purkinje->purkinje_cells = purkinje_cells;
    the_purkinje->number_of_purkinje_cells = total_purkinje_nodes;
    the_purkinje->num_active_purkinje_cells = total_purkinje_nodes;

    //the_grid->first_purkinje_cell = purkinje_cells[0];
    //the_grid->purkinje_cells = purkinje_cells;
    //the_grid->number_of_purkinje_cells = total_nodes;
    
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

    grid_cell = grid->first_cell;

    while(grid_cell != 0) {
        real_cpu center_x = grid_cell->center.x;
        real_cpu center_y = grid_cell->center.y;
        real_cpu center_z = grid_cell->center.z;

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

        grid_cell->translated_center.x = grid_cell->center.x - minx + (grid_cell->discretization.x/2.0f);
        grid_cell->translated_center.y = grid_cell->center.y - miny + (grid_cell->discretization.y/2.0f);
        grid_cell->translated_center.z = grid_cell->center.z - minz + (grid_cell->discretization.z/2.0f);

        grid_cell = grid_cell->next;
    }

}

static void sort_elements(struct element *cell_elements, int tam) {
    int i, j, min;
    struct element aux;
    for (i = 0; i < (tam-1); i++)
    {
        min = i;
        for (j = (i+1); j < tam; j++) {
            if(cell_elements[j].column < cell_elements[min].column)
                min = j;
        }
        if (cell_elements[i].column != cell_elements[min].column) {
            aux = cell_elements[i];
            cell_elements[i] = cell_elements[min];
            cell_elements[min] = aux;
        }
    }
}

#define for_each_cell(grid) \
    for(struct cell_node *cell = grid->first_cell; cell != NULL; cell = cell->next)

void grid_to_csr(struct grid *the_grid, real **A, int **IA, int **JA) {

    struct element element;

    arrpush(*IA, 0);

    int i = 0;
    int nnz = 0;
    size_t max_el = 0;
    int nnz_local;

    for_each_cell(the_grid) {

        bool insert = cell->active;

        if(arrlen(cell->elements) == 1 && cell->elements[0].value == 0.0) insert = false;

        if(insert) {

            if(i > 0) {
                int tmp = (*IA)[i - 1];
                arrpush(*IA, tmp + nnz_local);
            }

            nnz_local = 0;

            struct element *cell_elements = cell->elements;
            max_el = arrlen(cell_elements);

            sort_elements(cell_elements, max_el);

            for(int el = 0; el < max_el; el++) {
                element = cell_elements[el];
                if(element.value != 0) {
                    arrpush(*A, element.value);
                    arrpush(*JA, element.column);
                    nnz++;
                    nnz_local++;
                }
            }

            i++;

        }
    }

    arrpush(*IA, nnz);

}

void construct_grid_from_file(struct grid *grid, FILE *matrix_a, FILE *vector_b) {

    uint32_t n_cells;
    int num_lines_m = 0;
    int num_lines_v = 0;
    int nnz = 0;

    real_cpu **matrix = read_octave_mat_file_to_array(matrix_a, &num_lines_m, &nnz);
    real_cpu *vector = NULL;

    if(vector_b)
        vector = read_octave_vector_file_to_array(vector_b, &num_lines_v);

    initialize_and_construct_grid(grid, POINT3D(1.0, 1.0, 1.0));

    n_cells = grid->number_of_cells;
    while (n_cells < num_lines_m) {
        refine_grid(grid, 1);
        n_cells = grid->number_of_cells;
    }

    struct cell_node *cell = grid->first_cell;
    while (cell) {
        cell->active = false;
        cell = cell->next;
    }

    int item_count = 0;
    cell = grid->first_cell;
    while (item_count < num_lines_m) {
        cell->active = true;
        cell = cell->next;
        item_count++;
    }

    order_grid_cells(grid);

    cell = grid->first_cell;
    uint32_t cell_position;

    real_cpu m_value;

    for (int i = 0; i < num_lines_m; i++) {

        cell_position = cell->grid_position;
        m_value = matrix[cell_position][cell_position];
        struct element el;
        el.value = m_value;
        el.column = cell_position;
        el.cell = cell;

        arrsetcap(cell->elements, 7);
        arrput(cell->elements, el);

        for (int j = 0; j < num_lines_m; j++) {
            if (cell_position != j) {
                m_value = matrix[cell_position][j];

                if (m_value != 0.0) {
                    struct element el2;
                    el2.value = m_value;
                    el2.column = (uint32_t) j;

                    struct cell_node *aux = grid->first_cell;
                    while (aux) {
                        if (aux->grid_position == j)
                            break;
                        aux = aux->next;
                    }
                    el2.cell = aux;
                    arrput(cell->elements, el2);
                }
            }
        }

        cell = cell->next;
    }

    if(vector) {
        cell = grid->first_cell;
        for (int i = 0; i < num_lines_v; i++) {
            cell->b = vector[cell->grid_position];
            cell->v = 1.0;
            cell = cell->next;
        }
    }

    for (int i = 0; i < num_lines_m; i++) {
        free(matrix[i]);
    }

    free(matrix);
    free(vector);
}


// TODO: Include a parameter that links more cells
struct terminal* link_purkinje_to_endocardium (struct grid *the_grid)
{
    struct grid_purkinje *the_purkinje = the_grid->the_purkinje;

    struct graph *the_network = the_purkinje->the_network;

    uint32_t number_of_terminals = the_network->number_of_terminals;

    struct terminal *the_terminals = (struct terminal *)malloc(sizeof(struct terminal)*number_of_terminals);

    struct cell_node **ac_purkinje = the_purkinje->purkinje_cells;

    uint32_t j = 0;
    struct node *n = the_network->list_nodes;
    while (n != NULL)
    {
        if (n->num_edges == 1 && n->id != 0)
        {
            uint32_t n_active = the_grid->num_active_cells;
            struct cell_node **ac = the_grid->active_cells;
            
            uint32_t purkinje_index = n->id;
            struct node *purkinje_cell = n;
            
            uint32_t closest_index;
            double closest_dist = __DBL_MAX__;
            for (uint32_t i = 0; i < n_active; i++)
            {
                double dist = calc_norm(n->x,n->y,n->z,ac[i]->center.x,ac[i]->center.y,ac[i]->center.z);
                if (dist < closest_dist)
                {
                    closest_dist = dist;
                    closest_index = i;
                }
            }

            struct cell_node *endocardium_cell = ac[closest_index];
            uint32_t endocardium_index = ac[closest_index]->sv_position;
            
            the_terminals[j].endocardium_cell = endocardium_cell;
            the_terminals[j].endocardium_index = endocardium_index;
            the_terminals[j].purkinje_index = purkinje_index;
            the_terminals[j].purkinje_cell = purkinje_cell;

            // Change the position of the Purkinje terminal to be on the center of the Endocardium cell
            the_terminals[j].purkinje_cell->x = the_terminals[j].endocardium_cell->center.x;
            the_terminals[j].purkinje_cell->y = the_terminals[j].endocardium_cell->center.y;
            the_terminals[j].purkinje_cell->z = the_terminals[j].endocardium_cell->center.z;

            // Those lines are commented to prevent duplicate points ...
            //ac_purkinje[purkinje_index]->center.x = the_terminals[j].endocardium_cell->center.x;
            //ac_purkinje[purkinje_index]->center.y = the_terminals[j].endocardium_cell->center.y;
            //ac_purkinje[purkinje_index]->center.z = the_terminals[j].endocardium_cell->center.z;
            
            j++;

        }
        n = n->next;
    }

    //print_to_stdout_and_file("On 'link_purkinje_to_endocardium'\n");
    //for (uint32_t i = 0; i < number_of_terminals; i++)
    //    print_to_stdout_and_file("Terminal %u -- purkinje_index = %u -- endocardium_index = %u\n",i,the_terminals[i].purkinje_index,the_terminals[i].endocardium_index);

    return the_terminals;
}

void update_link_purkinje_to_endocardium (struct grid *the_grid, struct terminal *the_terminals)
{
    struct grid_purkinje *the_purkinje = the_grid->the_purkinje;

    struct graph *the_network = the_purkinje->the_network;

    uint32_t number_of_terminals = the_network->number_of_terminals;

    struct cell_node **ac_purkinje = the_purkinje->purkinje_cells;

    uint32_t j = 0;
    struct node *n = the_network->list_nodes;
    while (n != NULL)
    {
        if (n->num_edges == 1 && n->id != 0)
        {
            uint32_t n_active = the_grid->num_active_cells;
            struct cell_node **ac = the_grid->active_cells;
            
            uint32_t purkinje_index = n->id;
            struct node *purkinje_cell = n;
            
            uint32_t closest_index;
            double closest_dist = __DBL_MAX__;
            for (uint32_t i = 0; i < n_active; i++)
            {
                double dist = calc_norm(n->x,n->y,n->z,ac[i]->center.x,ac[i]->center.y,ac[i]->center.z);
                if (dist < closest_dist)
                {
                    closest_dist = dist;
                    closest_index = i;
                }
            }

            struct cell_node *endocardium_cell = ac[closest_index];
            uint32_t endocardium_index = ac[closest_index]->sv_position;
            
            the_terminals[j].endocardium_cell = endocardium_cell;
            the_terminals[j].endocardium_index = endocardium_index;
            //the_terminals[j].purkinje_index = purkinje_index;
            //the_terminals[j].purkinje_cell = purkinje_cell;

            // Change the position of the Purkinje terminal to be on the center of the Endocardium cell
            the_terminals[j].purkinje_cell->x = the_terminals[j].endocardium_cell->center.x;
            the_terminals[j].purkinje_cell->y = the_terminals[j].endocardium_cell->center.y;
            the_terminals[j].purkinje_cell->z = the_terminals[j].endocardium_cell->center.z;

            ac_purkinje[purkinje_index]->center.x = the_terminals[j].endocardium_cell->center.x;
            ac_purkinje[purkinje_index]->center.y = the_terminals[j].endocardium_cell->center.y;
            ac_purkinje[purkinje_index]->center.z = the_terminals[j].endocardium_cell->center.z;
            
            j++;

        }
        n = n->next;
    }

    //print_to_stdout_and_file("On 'update_link_purkinje_to_endocardium'\n");
    //for (uint32_t i = 0; i < number_of_terminals; i++)
    //    print_to_stdout_and_file("Terminal %u -- purkinje_index = %u -- endocardium_index = %u\n",i,the_terminals[i].purkinje_index,the_terminals[i].endocardium_index);
}