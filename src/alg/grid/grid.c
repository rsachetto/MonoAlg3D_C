//
// Created by sachetto on 29/09/17.
//

#include <assert.h>
#include <inttypes.h>

#include "grid.h"

#include "../../3dparty/stb_ds.h"
#include "../../utils/file_utils.h"
#include "../../utils/utils.h"

struct grid *new_grid() {
    struct grid *result = MALLOC_ONE_TYPE(struct grid);
    result->first_cell = NULL;
    result->active_cells = NULL;
    result->adaptive = false;

    result->refined_this_step = NULL;
    result->free_sv_positions = NULL;
    result->num_active_cells = result->number_of_cells = 0;

    arrsetcap(result->refined_this_step, 128);
    arrsetcap(result->free_sv_positions, 128);

    result->purkinje = NULL;

    result->extra_info = NULL;

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

    // Cell nodes. Variables names uses the center of the cube as initial reference.
    struct cell_node *right_front_top_cell, *right_front_down_cell, *right_back_top_cell, *right_back_down_cell,
        *left_front_top_cell, *left_front_down_cell, *left_back_top_cell, *left_back_down_cell;

    right_front_top_cell = new_cell_node();
    right_front_down_cell = new_cell_node();
    right_back_top_cell = new_cell_node();
    right_back_down_cell = new_cell_node();
    left_front_top_cell = new_cell_node();
    left_front_down_cell = new_cell_node();
    left_back_top_cell = new_cell_node();
    left_back_down_cell = new_cell_node();

    // Transition nodes.
    struct transition_node *front_transition_node;
    struct transition_node *back_transition_node;
    struct transition_node *top_transition_node;
    struct transition_node *down_transition_node;
    struct transition_node *right_transition_node;
    struct transition_node *left_transition_node;

    front_transition_node            = new_transition_node();
    back_transition_node             = new_transition_node();
    top_transition_node              = new_transition_node();
    down_transition_node             = new_transition_node();
    right_transition_node            = new_transition_node();
    left_transition_node             = new_transition_node();

    struct point_3d half_side_length    = POINT3D(side_length.x / 2.0f, side_length.y / 2.0, side_length.z / 2.0f);
    struct point_3d quarter_side_length = POINT3D(half_side_length.x / 2.0f, half_side_length.y / 2.0, half_side_length.z / 2.0f);

    //__________________________________________________________________________
    //              Initialization of transition nodes.
    //__________________________________________________________________________

    //Front transition node.
    set_transition_node_data(front_transition_node, 1, FRONT, NULL, right_front_down_cell, right_front_top_cell,
                             left_front_top_cell, left_front_down_cell);

    //Back transition node.
    set_transition_node_data(back_transition_node, 1, BACK, NULL, right_back_down_cell, right_back_top_cell,
                             left_back_top_cell, left_back_down_cell);
    //Top transition node.
    set_transition_node_data(top_transition_node, 1, TOP, NULL, right_back_top_cell, left_back_top_cell,
                             left_front_top_cell, right_front_top_cell);

    //Down transition node.
    set_transition_node_data(down_transition_node, 1, DOWN, NULL, right_back_down_cell, left_back_down_cell,
                             left_front_down_cell, right_front_down_cell);

    //Right transition node.
    set_transition_node_data(right_transition_node, 1, RIGHT, NULL, right_back_down_cell, right_back_top_cell,
                             right_front_top_cell, right_front_down_cell);

    //Left transition node.
    set_transition_node_data(left_transition_node, 1, LEFT, NULL, left_back_down_cell, left_back_top_cell,
                             left_front_top_cell, left_front_down_cell);

    //__________________________________________________________________________
    //                      Initialization of cell nodes.
    //__________________________________________________________________________


    void *neighbours[NUM_NEIGHBOURS];

    //right_front_top_cell neighbours
    neighbours[TOP]   = top_transition_node;
    neighbours[FRONT] = front_transition_node;
    neighbours[DOWN]  =  right_front_down_cell;
    neighbours[BACK]  = right_back_top_cell;
    neighbours[RIGHT] = right_transition_node;
    neighbours[LEFT]  = left_front_top_cell;

    set_cell_node_data(right_front_top_cell, half_side_length, 1, neighbours, NULL, left_front_top_cell, 0, 1,
                       POINT3D(half_side_length.x + quarter_side_length.x, half_side_length.y + quarter_side_length.y,
                               half_side_length.z + quarter_side_length.z));

    //left_front_top_cell neighbours
    neighbours[TOP]   = top_transition_node;
    neighbours[FRONT] = front_transition_node;
    neighbours[DOWN]  =  left_front_down_cell;
    neighbours[BACK]  = left_back_top_cell;
    neighbours[RIGHT] = right_front_top_cell;
    neighbours[LEFT]  = left_transition_node;

    set_cell_node_data(left_front_top_cell, half_side_length, 2, neighbours, right_front_top_cell, left_front_down_cell,
                       1, 2, POINT3D(quarter_side_length.x, half_side_length.y + quarter_side_length.y,
                               half_side_length.z + quarter_side_length.z));

    //left_front_down_cell neighbours
    neighbours[TOP]   = left_front_top_cell;
    neighbours[FRONT] = front_transition_node;
    neighbours[DOWN]  =  down_transition_node;
    neighbours[BACK]  = left_back_down_cell;
    neighbours[RIGHT] = right_front_down_cell;
    neighbours[LEFT]  = left_transition_node;

    set_cell_node_data(left_front_down_cell, half_side_length, 3, neighbours,
                       left_front_top_cell, right_front_down_cell, 2, 2,
                       POINT3D(quarter_side_length.x, quarter_side_length.y,
                               half_side_length.z + quarter_side_length.z));

    //right_front_down_cell neighbours
    neighbours[TOP]   = right_front_top_cell;
    neighbours[FRONT] = front_transition_node;
    neighbours[DOWN]  =  down_transition_node;
    neighbours[BACK]  = right_back_down_cell;
    neighbours[RIGHT] = right_transition_node;
    neighbours[LEFT]  = left_front_down_cell;

    set_cell_node_data(right_front_down_cell, half_side_length, 4, neighbours,
                       left_front_down_cell, right_back_down_cell, 3, 3,
                       POINT3D(half_side_length.x + quarter_side_length.x, quarter_side_length.y,
                               half_side_length.z + quarter_side_length.z));

    //right_back_down_cell neighbours
    neighbours[TOP]   = right_back_top_cell;
    neighbours[FRONT] = right_front_down_cell;
    neighbours[DOWN]  =  down_transition_node;
    neighbours[BACK]  = back_transition_node;
    neighbours[RIGHT] = right_transition_node;
    neighbours[LEFT]  = left_back_down_cell;

    set_cell_node_data(right_back_down_cell, half_side_length, 5, neighbours,
                       right_front_down_cell, left_back_down_cell, 4, 3,
                       POINT3D(half_side_length.x + quarter_side_length.x, quarter_side_length.y,
                               quarter_side_length.z));

    //left_back_down_cell neighbours
    neighbours[TOP]   = left_back_top_cell;
    neighbours[FRONT] = left_front_down_cell;
    neighbours[DOWN]  =  down_transition_node;
    neighbours[BACK]  = back_transition_node;
    neighbours[RIGHT] = right_back_down_cell;
    neighbours[LEFT]  = left_transition_node;

    set_cell_node_data(left_back_down_cell, half_side_length, 6, neighbours,
                       right_back_down_cell, left_back_top_cell, 5, 4,
                       POINT3D(quarter_side_length.x, quarter_side_length.y, quarter_side_length.z));

    //left_back_top_cell neighbours
    neighbours[TOP]   = top_transition_node;
    neighbours[FRONT] = left_front_top_cell;
    neighbours[DOWN]  =  left_back_down_cell;
    neighbours[BACK]  = back_transition_node;
    neighbours[RIGHT] = right_back_top_cell;
    neighbours[LEFT]  = left_transition_node;

    set_cell_node_data(left_back_top_cell, half_side_length, 7, neighbours, left_back_down_cell, right_back_top_cell, 6, 4,
                       POINT3D(quarter_side_length.x, half_side_length.y + quarter_side_length.y, quarter_side_length.z));

    //right_back_top_cell neighbours
    neighbours[TOP]   = top_transition_node;
    neighbours[FRONT] = right_front_top_cell;
    neighbours[DOWN]  =  right_back_down_cell;
    neighbours[BACK]  = back_transition_node;
    neighbours[RIGHT] = right_transition_node;
    neighbours[LEFT]  = left_back_top_cell;

    set_cell_node_data(right_back_top_cell, half_side_length, 8, neighbours, left_back_top_cell, NULL, 7, 5,
                       POINT3D(half_side_length.x + quarter_side_length.x, half_side_length.y + quarter_side_length.y,
                               quarter_side_length.z));

    the_grid->first_cell = right_front_top_cell;
    the_grid->number_of_cells = 8;
}

void print_grid(struct grid *the_grid, FILE *output_file) {

    real_cpu center_x, center_y, center_z, dx, dy, dz;
    real_cpu v;

    FOR_EACH_CELL(the_grid) {

        if(cell->active) {

            center_x = cell->center.x;
            center_y = cell->center.y;
            center_z = cell->center.z;

            dx = cell->discretization.x;
            dy = cell->discretization.y;
            dz = cell->discretization.z;

            v = cell->v;

            fprintf(output_file, "%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", center_x, center_y, center_z, dx, dy, dz, v);
        }
    }
}

void order_grid_cells(struct grid *the_grid) {

    // Here we allocate the maximum number of cells we will need for the whole simulation
    if(the_grid->active_cells == NULL) {
        the_grid->active_cells = MALLOC_ARRAY_OF_TYPE(struct cell_node *, the_grid->number_of_cells);
    }

    uint32_t counter = 0;
    FOR_EACH_CELL(the_grid) {
        if(cell->active) {
            cell->grid_position = counter;
            the_grid->active_cells[counter] = cell;
            cell->visible = get_visibility_mask(cell);
            counter++;
        }
    }

    the_grid->num_active_cells = counter;
}

void clean_grid(struct grid *the_grid) {

    assert(the_grid);

    struct cell_node *grid_cell = NULL;

    // TODO: Think about this function when the coupling happens ..
    // Delete nodes from the Purkinje network
    if(the_grid->purkinje && the_grid->purkinje->network->list_nodes != NULL) {

        grid_cell = the_grid->first_cell;

        // Then, delete the cells from the Purkinje network
        if(grid_cell) {
            while(grid_cell) {

                struct cell_node *next = grid_cell->next;
                free_cell_node(grid_cell);
                grid_cell = next;
            }
        }
    } else {
        // Delete the tissue cells
        // In order to release the memory allocated for the grid, the grid is
        // derefined to level 1. Thus, the grid shape is known and each node can
        // be easily reached.

        uint32_t number_of_cells = the_grid->number_of_cells;
        while(number_of_cells > 8) {
            derefine_all_grid(the_grid);
            number_of_cells = the_grid->number_of_cells;
        }

        grid_cell = the_grid->first_cell;

        if(grid_cell) {

            // Deleting transition nodes.
            free((struct transition_node *)(grid_cell->neighbours[FRONT]));
            free((struct transition_node *)(grid_cell->neighbours[RIGHT]));
            free((struct transition_node *)(grid_cell->neighbours[TOP]));
            free((struct transition_node *)(((struct cell_node *)(grid_cell->neighbours[DOWN]))->neighbours[DOWN]));
            free((struct transition_node *)(((struct cell_node *)(grid_cell->neighbours[BACK]))->neighbours[BACK]));
            free((struct transition_node *)(((struct cell_node *)(grid_cell->neighbours[LEFT]))->neighbours[LEFT]));

            // Deleting cells nodes.
            while(grid_cell) {
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

    if(the_grid->active_cells) {
        free(the_grid->active_cells);
    }

    arrfree(the_grid->refined_this_step);
    arrfree(the_grid->free_sv_positions);

    if(the_grid->purkinje) {
        free_graph(the_grid->purkinje->network);
        free(the_grid->purkinje); // TODO: Check for leaks with Valgrind
    }

    free(the_grid);
}

// Prints grid discretization matrix.
void print_grid_matrix(struct grid *the_grid, FILE *output_file) {

    assert(the_grid);
    assert(output_file);

    struct element element;
    element_array cell_elements;

    FOR_EACH_CELL(the_grid) {
        if(cell->active) {

            cell_elements = cell->elements;
            size_t max_el = arrlen(cell_elements);

            for(size_t i = 0; i < max_el; i++) {

                element = cell_elements[i];
                if(element.cell != NULL) {
                    fprintf(output_file,
                            "%" PRIu32 " "
                            "%" PRIu32 " %.15lf\n",
                            cell->grid_position + 1, (element.column) + 1, element.value);
                } else {
                    break;
                }
            }
        }
    }
}

int compare_elements(const void *a, const void *b) {
    if((*(struct element *)a).column < (*(struct element *)b).column)
        return -1;
    if((*(struct element *)a).column == (*(struct element *)b).column)
        return 0;
    if((*(struct element *)a).column > (*(struct element *)b).column)
        return 1;

    return 1; // Unreachable
}

void print_grid_matrix_as_octave_matrix(struct grid *the_grid, FILE *output_file) {

    assert(the_grid);
    assert(output_file);

    struct element element;
    element_array cell_elements;

    fprintf(output_file, "# Created by Monodomain solver\n");
    fprintf(output_file, "# name: Alg_grid_matrix\n");
    fprintf(output_file, "# type: sparse matrix\n");
    fprintf(output_file, "# nnz:                                    \n");
    fprintf(output_file, "# rows: %d\n", the_grid->num_active_cells);
    fprintf(output_file, "# columns: %d\n", the_grid->num_active_cells);

    int nnz = 0;

    FOR_EACH_CELL(the_grid) {
        if(cell->active) {

            cell_elements = cell->elements;
            size_t max_el = arrlen(cell_elements);

            qsort(cell_elements, max_el, sizeof(struct element), compare_elements);

            for(size_t i = 0; i < max_el; i++) {

                element = cell_elements[i];
                if(element.cell != NULL) {
                    nnz += 1;
                    fprintf(output_file,
                            "%" PRIu32 " "
                            "%" PRIu32 " %.15lf\n",
                            cell->grid_position + 1, (element.column) + 1, element.value);
                } else {
                    break;
                }
            }
        }
    }

    fseek(output_file, 84, SEEK_SET);
    fprintf(output_file, "%d", nnz);
}

void print_grid_vector(struct grid *the_grid, FILE *output_file, char name) {

    FOR_EACH_CELL(the_grid) {
        if(cell->active) {
            if(name == 'b')
                fprintf(output_file, "%.15lf\n", cell->b);
            else if(name == 'x')
                fprintf(output_file, "%.15lf\n", cell->v);
        }
    }
}

real_cpu *grid_vector_to_array(struct grid *the_grid, char name, uint32_t *num_lines) {

    *num_lines = the_grid->num_active_cells;
    real_cpu *vector = MALLOC_ARRAY_OF_TYPE(real_cpu, *num_lines);

    FOR_EACH_CELL(the_grid) {
        if(cell->active) {
            if(name == 'b')
                vector[cell->grid_position] = cell->b;
            else if(name == 'x')
                vector[cell->grid_position] = cell->v;
        }
    }

    return vector;
}

void save_grid_domain(struct grid *the_grid, const char *file_name) {
    FILE *f = fopen(file_name, "w");

    FOR_EACH_CELL(the_grid) {
        if(cell->active) {
            fprintf(f, "%lf,%lf,%lf,%lf,%lf,%lf\n", cell->center.x, cell->center.y, cell->center.z,
                    cell->discretization.x, cell->discretization.y, cell->discretization.z);
        }
    }
    fclose(f);
}

int get_num_refinement_steps_to_discretization(real_cpu side_len, real_cpu h) {

    int num_steps = 0;
    real_cpu aux = side_len;

    while(aux > h) {
        num_steps++;
        aux /= 2.0;
    }

    return num_steps - 1;
}

void initialize_grid_purkinje(struct grid *the_grid) {
    assert(the_grid);
    assert(the_grid->purkinje);

    the_grid->purkinje->number_of_purkinje_cells = 0;
}

void construct_grid_purkinje(struct grid *the_grid) {
    assert(the_grid);
    assert(the_grid->purkinje);

    struct grid_purkinje *the_purkinje = the_grid->purkinje;

    real_cpu side_length_x = the_purkinje->network->dx;
    real_cpu side_length_y = the_purkinje->network->dx;
    real_cpu side_length_z = the_purkinje->network->dx;

    struct point_3d side_length = POINT3D(side_length_x, side_length_y, side_length_z);

    uint32_t total_purkinje_nodes = the_purkinje->network->total_nodes;

    // Create the Purkinje cells array
    struct cell_node **purkinje_cells;
    purkinje_cells = MALLOC_ARRAY_OF_TYPE(struct cell_node *, total_purkinje_nodes);
    for(int i = 0; i < total_purkinje_nodes; i++)
        purkinje_cells[i] = new_cell_node();

    // Pass through the Purkinje graph and set the cell nodes.
    struct node *n = the_purkinje->network->list_nodes;
    void **neighbours = NULL;
    for(uint32_t i = 0; i < total_purkinje_nodes; i++) {

        if(i == 0)
            set_cell_node_data(purkinje_cells[i], side_length, 0, neighbours, NULL,
                               purkinje_cells[i + 1], i, 0, POINT3D(n->pos[0], n->pos[1], n->pos[2]));

        else if(i == total_purkinje_nodes - 1)
            set_cell_node_data(purkinje_cells[i], side_length, 0, neighbours,
                               purkinje_cells[i - 1], NULL, i, 0, POINT3D(n->pos[0], n->pos[1], n->pos[2]));
        else
            set_cell_node_data(purkinje_cells[i], side_length, 0, neighbours,
                               purkinje_cells[i - 1], purkinje_cells[i + 1], i, 0, POINT3D(n->pos[0], n->pos[1], n->pos[2]));

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
}

void initialize_and_construct_grid_purkinje(struct grid *the_grid) {
    assert(the_grid);

    initialize_grid_purkinje(the_grid);
    construct_grid_purkinje(the_grid);
}

static void sort_elements(struct element *cell_elements, size_t size) {
    int i, j, min;
    struct element aux;
    for(i = 0; i < (size - 1); i++) {
        min = i;
        for(j = (i + 1); j < size; j++) {
            if(cell_elements[j].column < cell_elements[min].column)
                min = j;
        }
        if(cell_elements[i].column != cell_elements[min].column) {
            aux = cell_elements[i];
            cell_elements[i] = cell_elements[min];
            cell_elements[min] = aux;
        }
    }
}

void grid_to_csr(struct grid *the_grid, float **A, int **IA, int **JA, bool is_purkinje) {
    grid_to_csr_for_ecg(the_grid, A, IA, JA, is_purkinje, false);
}

void grid_to_csr_for_ecg(struct grid *the_grid, float **A, int **IA, int **JA, bool is_purkinje, bool for_ecg) {

    struct element element;

    arrpush(*IA, 0);

    int nnz = 0;
    int nnz_local;

    struct cell_node *cell;

    if(is_purkinje) {
        cell = the_grid->purkinje->first_cell;
    } else {
        cell = the_grid->first_cell;
    }

    int i = 0;
    for(; cell != NULL; cell = cell->next) {

        bool insert = cell->active;

        if(arrlen(cell->elements) == 1 && cell->elements[0].value == 0.0) {
            insert = false;
        }

        if(insert) {

            if(i > 0) {
                int tmp = (*IA)[i - 1];
                arrpush(*IA, tmp + nnz_local);
            }

            nnz_local = 0;

            struct element *cell_elements = NULL;
            size_t max_el = arrlen(cell->elements);

            for(int el = 0; el < max_el; el++) {
                arrpush(cell_elements, cell->elements[el]);
            }

            sort_elements(cell_elements, max_el);

            for(int el = 0; el < max_el; el++) {
                element = cell_elements[el];
                if(for_ecg) {
                    if(element.value_ecg != 0) {
                        arrpush(*A, element.value_ecg);
                        arrpush(*JA, element.column);
                        nnz++;
                        nnz_local++;
                    }
                }
                else {
                    if(element.value != 0) {
                        arrpush(*A, element.value);
                        arrpush(*JA, element.column);
                        nnz++;
                        nnz_local++;
                    }

                }
            }

            arrfree(cell_elements);
            i++;
        }
    }

    arrpush(*IA, nnz);

}

void construct_grid_from_file(struct grid *the_grid, FILE *matrix_a, FILE *vector_b) {

    uint32_t n_cells;
    uint64_t num_lines_m = 0;
    uint64_t num_lines_v = 0;
    uint64_t nnz = 0;

    real_cpu **matrix = read_octave_mat_file_to_array(matrix_a, &num_lines_m, &nnz);
    real_cpu *vector = NULL;

    if(vector_b)
        vector = read_octave_vector_file_to_array(vector_b, &num_lines_v);

    initialize_and_construct_grid(the_grid, POINT3D(1.0, 1.0, 1.0));

    n_cells = the_grid->number_of_cells;
    while(n_cells < num_lines_m) {
        refine_grid(the_grid, 1);
        n_cells = the_grid->number_of_cells;
    }

    FOR_EACH_CELL(the_grid) {
        cell->active = false;
    }

    int item_count = 0;

    struct cell_node *cell = the_grid->first_cell;
    while(item_count < num_lines_m) {
        cell->active = true;
        cell = cell->next;
        item_count++;
    }

    order_grid_cells(the_grid);

    cell = the_grid->first_cell;
    uint32_t cell_position;

    real_cpu m_value;

    for(uint64_t i = 0; i < num_lines_m; i++) {

        cell_position = cell->grid_position;
        m_value = matrix[cell_position][cell_position];
        struct element el;
        el.value = m_value;
        el.column = cell_position;
        el.cell = cell;

        arrsetcap(cell->elements, 7);
        arrput(cell->elements, el);

        for(uint64_t j = 0; j < num_lines_m; j++) {
            if(cell_position != j) {
                m_value = matrix[cell_position][j];

                if(m_value != 0.0) {
                    struct element el2;
                    el2.value = m_value;
                    el2.column = (uint32_t)j;

                    struct cell_node *aux = the_grid->first_cell;
                    while(aux) {
                        if(aux->grid_position == j)
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
        cell = the_grid->first_cell;
        for(uint64_t i = 0; i < num_lines_v; i++) {
            cell->b = vector[cell->grid_position];
            cell->v = 1.0;
            cell = cell->next;
        }
    }

    for(uint64_t i = 0; i < num_lines_m; i++) {
        free(matrix[i]);
    }

    free(matrix);
    free(vector);
}


struct terminal *link_purkinje_to_tissue (struct grid *the_grid) {

    struct graph *the_network = the_grid->purkinje->network;

    if (the_network->has_pmj_location) {
        return link_purkinje_to_tissue_using_pmj_locations(the_grid);
    } else {
        return link_purkinje_to_tissue_default(the_grid);
    }
}

void free_terminals (struct terminal *the_terminals, const uint32_t number_of_terminals) {

    for (uint32_t i = 0; i < number_of_terminals; i++) {

        if (the_terminals[i].purkinje_cell)
            the_terminals[i].purkinje_cell = NULL;

        if (the_terminals[i].tissue_cells) {
            for (uint32_t j = 0; j < arrlen(the_terminals[i].tissue_cells); j++)
                the_terminals[i].tissue_cells[j] = NULL;
            arrfree(the_terminals[i].tissue_cells);
        }
    }
    free(the_terminals);
}

struct terminal* link_purkinje_to_tissue_default (struct grid *the_grid) {

    struct graph *the_network = the_grid->purkinje->network;

    uint32_t number_of_terminals = the_network->number_of_terminals;
    const real_cpu pmj_scale = the_network->pmj_scale;
    const real_cpu nmin_pmj = the_network->nmin_pmj;
    const real_cpu nmax_pmj = the_network->nmax_pmj;

    struct terminal *the_terminals = MALLOC_ARRAY_OF_TYPE(struct terminal, number_of_terminals);

    uint32_t j = 0;
    struct node *n = the_network->list_nodes;
    while(n != NULL) {

        if( is_terminal(n) ) {

            uint32_t n_active = the_grid->num_active_cells;
            struct cell_node **ac = the_grid->active_cells;

            // Save the current Purkinje terminal cell
            struct node *purkinje_cell = n;
            the_terminals[j].purkinje_cell = purkinje_cell;

            // All the terminals are active
            the_terminals[j].active = true;

            // Search for all the tissue cells that are within the sphere that
            // has a radius less or equal to 'pmj_scale'
            uint32_t *tissue_cells_to_link = NULL;
            real_cpu *dist_array = NULL;
            real_cpu scale = pmj_scale;
            while (arrlen(tissue_cells_to_link) < nmin_pmj) {

                for(uint32_t i = 0; i < n_active; i++) {

                    real_cpu dist = calc_norm(n->pos[0], n->pos[1], n->pos[2],\
                                            ac[i]->center.x, ac[i]->center.y, ac[i]->center.z);

                    if(dist < scale) {
                        arrput(dist_array,dist);
                        arrput(tissue_cells_to_link,i);
                    }
                }

                // Increase the 'pmj_scale' by 10%
                if (arrlen(tissue_cells_to_link) < nmin_pmj) {
                    scale *= 1.1;
                }
            }

            // QuickSort: Sort the distance array together with the indexes from the tissue cells to link
            sort_vector_by_distance(dist_array,tissue_cells_to_link,arrlen(dist_array));

            // Save the tissue cells indexes we are going to link
            // until we achieve the maximum number of points inside the PMJ region or
            // until the maximum size of the link array is reached
            uint32_t size = (arrlen(tissue_cells_to_link) < nmax_pmj) ? arrlen(tissue_cells_to_link) : nmax_pmj;
            the_terminals[j].tissue_cells = NULL;
            for (uint32_t i = 0; i < size; i++) {
                uint32_t index = tissue_cells_to_link[i];

                arrput(the_terminals[j].tissue_cells,ac[index]);
            }

            arrfree(tissue_cells_to_link);

            j++;
        }
        n = n->next;
    }

    return the_terminals;
}

struct terminal* link_purkinje_to_tissue_using_pmj_locations (struct grid *the_grid) {

    struct graph *the_network = the_grid->purkinje->network;

    uint32_t number_of_terminals = the_network->number_of_terminals;
    real_cpu nmin_pmj = the_network->nmin_pmj;

    struct terminal *the_terminals = MALLOC_ARRAY_OF_TYPE(struct terminal,  number_of_terminals);

    // Set all the terminal Purkinje cells
    uint32_t j = 0;
    struct node *n = the_network->list_nodes;
    while(n != NULL) {
        if( is_terminal(n) ) {
            // Save the current Purkinje terminal cell
            struct node *purkinje_cell = n;
            the_terminals[j].purkinje_cell = purkinje_cell;

            // Set all the terminals to inactivated initially
            the_terminals[j].active = false;

            j++;
        }
        n = n->next;
    }

    // Load the PMJ locations from the file and activate the closest terminal from the Purkinje network to each PMJ
    char *pmj_location_filename = the_network->pmj_location_filename;
    set_active_terminals(the_terminals,number_of_terminals,pmj_location_filename);

    // Calculate the distance between each tissue cell to the Purkinje terminal and store them in a table
    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;
    real_cpu *dist_array = MALLOC_ARRAY_OF_TYPE(real_cpu, number_of_terminals*n_active);
    uint32_t *tissue_ids = MALLOC_ARRAY_OF_TYPE(uint32_t, n_active);
    bool *tissue_taken = MALLOC_ARRAY_OF_TYPE(bool, n_active);

    // Grid cells loop
    for (uint32_t i = 0; i < n_active; i++) {
        
        // Purkinje terminals loop
        for (j = 0; j < number_of_terminals; j++) {
            n = the_terminals[j].purkinje_cell;
            real_cpu dist = calc_norm(n->pos[0], n->pos[1], n->pos[2], ac[i]->center.x, ac[i]->center.y, ac[i]->center.z);

            dist_array[j*n_active+i] = dist;
        }
        // Set the current tissue cell as not taken
        tissue_taken[i] = false;

        // Reset the tissue cell index array
        tissue_ids[i] = i;
    }

    for (uint32_t i = 0; i < number_of_terminals; i++) {

        if (the_terminals[i].active) {
            
            real_cpu *tissue_dist = &dist_array[i*n_active];
            nmin_pmj = the_terminals[i].purkinje_cell->nmin_pmj;
            sort_vector_by_distance(tissue_dist,tissue_ids, n_active);
            
            // Fill the 'tissue_cells_to_link' array until we achieve the 'nmin_pmjs'
            uint32_t *tissue_cells_to_link = NULL;
            j = 0;
            while (arrlen(tissue_cells_to_link) < nmin_pmj) {
                uint32_t id = tissue_ids[j];
                if (!tissue_taken[id]) {
                    arrput(tissue_cells_to_link,tissue_ids[j]);
                    tissue_taken[id] = true;
                }
                j++;
            }

            // Set the reference to the tissue cells
            the_terminals[i].tissue_cells = NULL;
            for (j = 0; j < nmin_pmj; j++) {

                uint32_t index = tissue_cells_to_link[j];
                arrput(the_terminals[i].tissue_cells,ac[index]);
            }

            // Reset the tissue cell index array
            for (j = 0; j < n_active; j++)
                tissue_ids[j] = j;
            arrfree(tissue_cells_to_link);
        } 
        else {
            the_terminals[i].tissue_cells = NULL;
        }
    }
    
    free(tissue_taken);
    free(tissue_ids);
    free(dist_array);

    //print_terminals(the_terminals,number_of_terminals);
   
    return the_terminals;
}

// TODO: Revise this function ... Think about the adaptative scenario
void update_link_purkinje_to_endocardium(struct grid *the_grid, struct terminal *the_terminals) {
/*
    struct grid_purkinje *the_purkinje = the_grid->purkinje;

    struct graph *the_network = the_purkinje->network;

    struct cell_node **ac_purkinje = the_purkinje->purkinje_cells;

    uint32_t j = 0;
    struct node *n = the_network->list_nodes;
    while(n != NULL) {
        if(n->num_edges == 1 && n->id != 0) {
            uint32_t n_active = the_grid->num_active_cells;
            struct cell_node **ac = the_grid->active_cells;

            uint32_t purkinje_index = n->id;

            uint32_t closest_index = 0;
            double closest_dist = __DBL_MAX__;

            for(uint32_t i = 0; i < n_active; i++) {
                double dist = calc_norm(n->x, n->y, n->z, ac[i]->center.x, ac[i]->center.y, ac[i]->center.z);
                if(dist < closest_dist) {
                    closest_dist = dist;
                    closest_index = i;
                }
            }

            struct cell_node *endocardium_cell = ac[closest_index];
            uint32_t endocardium_index = ac[closest_index]->sv_position;

            the_terminals[j].endocardium_cell = endocardium_cell;
            the_terminals[j].endocardium_index = endocardium_index;

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
*/
}

void set_active_terminals (struct terminal *the_terminals, const uint32_t number_of_terminals, const char filename[]) {
    // Read all the PMJ coordinates
    FILE *file = fopen(filename,"r");
    if (!file) {
        printf("Error! Reading pmj_location_file!\n");
        exit(1);
    }
    struct node *pmjs = NULL;
    real_cpu *rpmjs = NULL;
    uint32_t *nmin_pmjs = NULL;

    char str[200];
    while (fscanf(file,"%s",str) != EOF) {
        if (strcmp(str,"POINTS") == 0) break;
    }
    uint32_t num_pmjs, dummy;
    fscanf(file,"%u %s",&num_pmjs,str);
    for (uint32_t i = 0; i < num_pmjs; i++) {
        double pos[3];
        fscanf(file,"%lf %lf %lf",&pos[0],&pos[1],&pos[2]);

        struct node n;
        for (int j = 0; j < 3; j++) {
            n.pos[j] = pos[j];
        }
        arrput(pmjs,n);

    }
    fscanf(file,"%s %u %u",str,&dummy,&dummy);
    for (uint32_t i = 0; i < num_pmjs; i++) {
        fscanf(file,"%u %u",&dummy,&dummy);
    }
    // Read the POINT_DATA section, if it exist
    bool has_point_data = false;
    if (fscanf(file,"%s",str) != EOF) {
        if (strcmp(str,"POINT_DATA") == 0) {

            while (fscanf(file,"%s",str) != EOF) {
                if (strcmp(str,"float") == 0) break;
            }

            for (int i = 0; i < num_pmjs; i++) {
                real_cpu rpmj;
                fscanf(file,"%lf",&rpmj);
                arrput(rpmjs,rpmj);
            }

            while (fscanf(file,"%s",str) != EOF) {
                if (strcmp(str,"float") == 0) break;
            }

            for (int i = 0; i < num_pmjs; i++) {
                uint32_t nmin_pmj;
                fscanf(file,"%u",&nmin_pmj);
                arrput(nmin_pmjs,nmin_pmj);
            }

            has_point_data = true;
        }
    }
    fclose(file);

    // Activate only the closest terminal to each PMJ location
    for (uint32_t i = 0; i < num_pmjs; i++) {
        uint32_t min_index = 0;
        double min_dist = __DBL_MAX__;
        for (uint32_t j = 0; j < number_of_terminals; j++) {
            struct node *tmp = the_terminals[j].purkinje_cell;
            double dist = calc_norm(tmp->pos[0],tmp->pos[1],tmp->pos[2],pmjs[i].pos[0],pmjs[i].pos[1],pmjs[i].pos[2]);
            if (dist < min_dist) {
                min_dist = dist;
                min_index = j;
            }
        }
        if (has_point_data) {
            the_terminals[min_index].purkinje_cell->nmin_pmj = nmin_pmjs[i];
            the_terminals[min_index].purkinje_cell->rpmj = rpmjs[i];    
        }
        the_terminals[min_index].active = true;
    }

    arrfree(pmjs);
    if (rpmjs) arrfree(rpmjs);
    if (nmin_pmjs) arrfree(nmin_pmjs);
}

void print_terminals (struct terminal *the_terminals, const uint32_t number_of_terminals)
{
    for (uint32_t j = 0; j < number_of_terminals; j++) {
        
        if (the_terminals[j].active) {
            
            uint32_t num_tissue_cells = arrlen(the_terminals[j].tissue_cells);
            uint32_t num_border_cells = 0;
            for (uint32_t i = 0; i < num_tissue_cells; i++) {
                struct cell_node **tcells = the_terminals[j].tissue_cells;
                uint8_t visibility_flag = get_visibility_mask(tcells[i]);   // Tricky way to discover if a cell has all the neighbours
                if (visibility_flag != 0)
                    num_border_cells++;
                printf("\tTissue cell %u = (%g,%g,%g) || Visibility flag = %u\n",tcells[i]->sv_position,tcells[i]->center.x,tcells[i]->center.y,tcells[i]->center.z,visibility_flag);
            }
            printf("[Terminal %u] is linked to %u tissue cells. There are %u/%u border cells\n",j,num_tissue_cells,num_border_cells,num_tissue_cells);
        }
    }
}
