//
// Created by sachetto on 06/10/17.
//

#include "teste_cg.h"
#include <stdio.h>
#include <stdlib.h>
#include "../solvers/ode_solver.h"
#include "../utils/output_utils.h"
#include "../solvers/linear_system_solver.h"

#include <criterion/criterion.h>

struct grid* construct_grid_from_file(struct grid *grid, FILE *matrix_a, FILE *vector_b ) {

    uint64_t n_cells;
    char * line_a = NULL;
    char * line_b = NULL;
    size_t len = 0;
    ssize_t  read;
    char *sep = " ";
    int count = 0;
    int num_lines = 0;
    int nnz = 0;

    do {
        read = getline(&line_a, &len, matrix_a);
        sds *tmp = sdssplitlen(line_a, (int)strlen(line_a), sep, (int)strlen(sep), &count);
        if(count) {
            if(strcmp(tmp[1], "columns:") == 0) {
                num_lines = atoi(tmp[2]);
            }
            if(strcmp(tmp[1], "nnz:") == 0) {
                nnz = atoi(tmp[2]);
            }

        }
        sdsfreesplitres(tmp, count);
    } while(line_a[0] == '#');

    do {
        read = getline(&line_b, &len, vector_b);
    } while(line_b[0] == '#');

    int max_el = num_lines;

    cr_assert(num_lines);
    cr_assert(nnz);


    initialize_and_construct_grid(grid, 1.0, num_lines);

    n_cells = grid->number_of_cells;
    while(n_cells < num_lines) {
        refine_grid(grid, 1);
        n_cells = grid->number_of_cells;
    }

    struct cell_node *cell = grid->first_cell;
    while(cell) {
        cell->active = false;
        cell = cell->next;
    }


    int item_count = 0;
    cell = grid->first_cell;
    while(item_count < num_lines) {
        cell->active = true;
        cell = cell->next;
        item_count++;
    }


    double **matrix = (double**) malloc(num_lines*sizeof(double*));
    for(int i = 0; i < num_lines; i++) {
        matrix[i] = (double*)calloc(num_lines, sizeof(double));
    }

    double *vector = (double*) malloc(num_lines*sizeof(double));

    item_count = 0;
    cell = grid->first_cell;
    int m_line, m_column;
    double m_value;

    while(item_count < nnz) {

        sds *tmp = sdssplitlen(line_a, (int)strlen(line_a), sep, (int)strlen(sep), &count);
        //printf("%s", line_a);
        if(tmp[0][0] != '\n') {
            m_line = atoi(tmp[0]);
            m_column = atoi(tmp[1]);
            m_value = atof(tmp[2]);

            matrix[m_line-1][m_column-1] = m_value;
        }
        sdsfreesplitres(tmp, count);

        item_count++;
        read = getline(&line_a, &len, matrix_a);
    }

    item_count = 0;
    while((item_count < num_lines) && read) {
        sds *tmp = sdssplitlen(line_b, (int)strlen(line_b), sep, (int)strlen(sep), &count);

        if(tmp[0][0] != '\n') {
            vector[item_count] = atof(tmp[1]);
        }

        sdsfreesplitres(tmp, count);

        item_count++;
        read = getline(&line_b, &len, vector_b);
    }

    order_grid_cells(grid);

    cr_assert_eq(num_lines, grid->num_active_cells);

    cell = grid->first_cell;
    uint64_t cell_position;

    for(int i = 0; i < num_lines; i++) {

        cell_position = cell->grid_position;
        m_value = matrix[cell_position][cell_position];
        struct element el;
        el.value = m_value;
        el.column = cell_position;
        el.cell = cell;

        cell->elements = new_element_array(max_el);
        cell->elements[0] = el;

        int next_element = 1;

        for(int j = 0; j < num_lines; j++) {
            if(cell_position != j) {
                //we should find a better way to compare the floating points
                m_value = matrix[cell_position][j];

                if (m_value != 0.0) {
                    struct element el2;
                    el2.value = m_value;
                    el2.column = (uint64_t) j;

                    struct cell_node *aux = grid->first_cell;
                    while (aux) {
                        if (aux->grid_position == j) break;
                        aux = aux->next;
                    }
                    el2.cell = aux;
                    cell->elements[next_element] = el2;
                    next_element++;
                }
            }
        }

        cell = cell->next;
    }

    cell = grid->first_cell;
    for(int i = 0; i < num_lines; i++) {
        cell->b = vector[cell->grid_position];
        cell->v = 1.0;
        cell = cell->next;
    }

    //print_grid_matrix(grid, stdout);
    //print_grid_vector(grid, stdout, 'b');

    if(line_a) {
        free(line_a);
    }
    if(line_b) {
        free(line_b);
    }

    for(int i = 0; i < num_lines; i++) {
        free(matrix[i]);
    }

    free(matrix);
    free(vector);

}

main() {
    FILE *A = fopen("src/tests/A.txt", "r");
    FILE *B = fopen("src/tests/B.txt", "r");
    double error;

    struct grid *grid = (struct grid*)malloc(sizeof(struct grid));
    cr_assert(grid);

    construct_grid_from_file(grid, A, B);
    conjugate_gradient(grid, 200, 1e-10, false, &error);

    print_grid_vector(grid, stdout, 'x');

    clean_and_free_grid(grid);
    fclose(A);
    fclose(B);
}

/*Test(misc, cg1) {
    FILE *A = fopen("src/tests/A.txt", "r");
    FILE *B = fopen("src/tests/B.txt", "r");
    cr_assert(A);
    cr_assert(B);

    construct_grid_from_file(A, B);
    fclose(A);
}*/



void run_cg1() {




}