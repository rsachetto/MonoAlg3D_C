//
// Created by sachetto on 06/10/17.
//
#include "../monodomain/linear_system_solver.h"
#include "../monodomain/output_utils.h"
#include <criterion/criterion.h>

double *read_octave_vector_file_to_array (FILE *vec_file, int *num_lines);

double **read_octave_mat_file_to_array (FILE *matrix_file, int *num_lines, int *nnz);

void construct_grid_from_file (struct grid *grid, FILE *matrix_a, FILE *vector_b) {

    uint32_t n_cells;
    int num_lines_m = 0;
    int num_lines_v = 0;
    int nnz = 0;

    double **matrix = read_octave_mat_file_to_array (matrix_a, &num_lines_m, &nnz);
    double *vector = read_octave_vector_file_to_array (vector_b, &num_lines_v);

    cr_assert_eq (num_lines_m, num_lines_v);
    cr_assert (nnz);

    initialize_and_construct_grid (grid, 1.0);

    n_cells = grid->number_of_cells;
    while (n_cells < num_lines_m) {
        refine_grid (grid, 1);
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

    order_grid_cells (grid);

    cr_assert_eq (num_lines_m, grid->num_active_cells);

    cell = grid->first_cell;
    uint32_t cell_position;

    double m_value;

    for (int i = 0; i < num_lines_m; i++) {

        cell_position = cell->grid_position;
        m_value = matrix[cell_position][cell_position];
        struct element el;
        el.value = m_value;
        el.column = cell_position;
        el.cell = cell;

        sb_reserve(cell->elements, 7);
        sb_push(cell->elements, el);

        for (int j = 0; j < num_lines_m; j++) {
            if (cell_position != j) {
                // we should find a better way to compare the floating points
                m_value = matrix[cell_position][j];

                if (m_value != 0.0) {
                    struct element el2;
                    el2.value = m_value;
                    el2.column = (uint32_t)j;

                    struct cell_node *aux = grid->first_cell;
                    while (aux) {
                        if (aux->grid_position == j)
                            break;
                        aux = aux->next;
                    }
                    el2.cell = aux;
                    sb_push(cell->elements, el2);
                }
            }
        }

        cell = cell->next;
    }

    cell = grid->first_cell;
    for (int i = 0; i < num_lines_v; i++) {
        cell->b = vector[cell->grid_position];
        cell->v = 1.0;
        cell = cell->next;
    }

    // print_grid_matrix(grid, stdout);
    // print_grid_vector(grid, stdout, 'b');

    for (int i = 0; i < num_lines_m; i++) {
        free (matrix[i]);
    }

    free (matrix);
    free (vector);
}

double **read_octave_mat_file_to_array (FILE *matrix_file, int *num_lines, int *nnz) {
    const char *sep = " ";
    char *line_a = NULL;
    size_t len;
    int count;

    do {
        getline (&line_a, &len, matrix_file);
        sds *tmp = sdssplitlen (line_a, (int)strlen (line_a), sep, (int)strlen (sep), &count);
        if (count) {
            if (strcmp (tmp[1], "columns:") == 0) {
                (*num_lines) = atoi (tmp[2]);
            }
            if (strcmp (tmp[1], "nnz:") == 0) {
                (*nnz) = atoi (tmp[2]);
            }
        }
        sdsfreesplitres (tmp, count);
    } while ((line_a)[0] == '#');

    double **matrix = (double **)malloc (*num_lines * sizeof (double *));

    for (int i = 0; i < *num_lines; i++) {
        matrix[i] = (double *)calloc (*num_lines, sizeof (double));
    }

    int item_count = 0;
    int m_line, m_column;
    double m_value;

    while (item_count < *nnz) {

        sds *tmp = sdssplitlen (line_a, (int)strlen (line_a), sep, (int)strlen (sep), &count);
        if (tmp[0][0] != '\n') {
            m_line = atoi (tmp[0]);
            m_column = atoi (tmp[1]);
            m_value = atof (tmp[2]);

            matrix[m_line - 1][m_column - 1] = m_value;
        }
        sdsfreesplitres (tmp, count);

        item_count++;
        getline (&line_a, &len, matrix_file);
    }

    if (line_a)
        free (line_a);

    return matrix;
}

double *read_octave_vector_file_to_array (FILE *vec_file, int *num_lines) {

    ssize_t read;
    size_t len;
    char *line_b = NULL;
    int count;
    char *sep = " ";

    do {
        read = getline (&line_b, &len, vec_file);
        sds *tmp = sdssplitlen (line_b, (int)strlen (line_b), sep, (int)strlen (sep), &count);
        if (count) {
            if (strcmp (tmp[1], "rows:") == 0) {
                (*num_lines) = atoi (tmp[2]);
            }
        }
        sdsfreesplitres (tmp, count);
    } while ((line_b)[0] == '#');

    double *vector = (double *)malloc (*num_lines * sizeof (double));

    int item_count = 0;
    while ((item_count < *num_lines) && read) {
        sds *tmp = sdssplitlen (line_b, (int)strlen (line_b), sep, (int)strlen (sep), &count);

        if (tmp[0][0] != '\n') {
            vector[item_count] = atof (tmp[1]);
        }

        sdsfreesplitres (tmp, count);

        item_count++;
        read = getline (&line_b, &len, vec_file);
    }

    if (line_b)
        free (line_b);

    return vector;
}

void run_cg (bool jacobi) {
    FILE *A = fopen ("src/tests/A.txt", "r");
    FILE *B = fopen ("src/tests/B.txt", "r");
    FILE *X = fopen ("src/tests/X.txt", "r");
    double error;

    cr_assert(A);
    cr_assert(B);
    cr_assert(X);


    struct grid *grid = (struct grid *)malloc (sizeof (struct grid));
    cr_assert (grid);

    construct_grid_from_file (grid, A, B);

    int nt = 1;

#if defined(_OPENMP)
    nt = omp_get_max_threads();
#endif


    if(jacobi) {
        printf("Testing CG with jacobi preconditioner using %d treads\n", nt);
    }

    else {
        printf("Testing CG using %d treads\n", nt);
    }

    conjugate_gradient (grid, 200, 1e-16, jacobi, &error);

    int n_lines1;
    uint32_t n_lines2;

    double *x = read_octave_vector_file_to_array (X, &n_lines1);
    double *x_grid = grid_vector_to_array (grid, 'x', &n_lines2);

    cr_assert_eq (n_lines1, n_lines2);

    for (int i = 0; i < n_lines1; i++) {
        cr_assert_float_eq (x[i], x_grid[i], 1e-10);
    }

    clean_and_free_grid (grid);
    fclose (A);
    fclose (B);
}

Test (solvers, cg_jacobi) {
    run_cg (true);
}

Test (solvers, cg_no_jacobi) {
    run_cg (false);
}

Test (utils, stretchy_buffer_int) {

    int *v = NULL;

    sb_reserve(v, 1);

    cr_assert_eq(sb_count(v) ,0);
    cr_assert_eq(stb__sbm(v)  ,1);

    sb_push(v, 0);
    sb_push(v, 1);
    sb_push(v, 2);

    cr_assert_eq(sb_count(v) ,3);

    cr_assert_eq(v[0], 0);
    cr_assert_eq(v[1], 1);
    cr_assert_eq(v[2], 2);
}

Test (utils, stretchy_buffer_float) {

    float *v = NULL;

    sb_reserve(v, 1);

    cr_assert_eq(sb_count(v) ,0);
    cr_assert_eq(stb__sbm(v)  ,1);

    sb_push(v, 0);
    sb_push(v, 0.5);
    sb_push(v, 2.5);

    cr_assert_eq(sb_count(v) ,3);

    cr_assert_float_eq(v[0], 0.0, 1e-10);
    cr_assert_float_eq(v[1], 0.5, 1e-10);
    cr_assert_float_eq(v[2], 2.5, 1e-10);
}

Test (utils, stretchy_buffer_double) {

    double *v = NULL;

    sb_reserve(v, 1);

    cr_assert_eq(sb_count(v) ,0);
    cr_assert_eq(stb__sbm(v)  ,1);

            sb_push(v, 0);
            sb_push(v, 0.5);
            sb_push(v, 2.5);

    cr_assert_eq(sb_count(v) ,3);

    cr_assert_float_eq(v[0], 0.0, 1e-10);
    cr_assert_float_eq(v[1], 0.5, 1e-10);
    cr_assert_float_eq(v[2], 2.5, 1e-10);
}

Test (utils, stretchy_buffer_element) {

    struct element *v = NULL;

    sb_reserve(v, 1);

    cr_assert_eq(sb_count(v) ,0);
    cr_assert_eq(stb__sbm(v)  ,1);

    struct cell_node *c  = new_cell_node();
    struct element a = {0, 1, c};

    sb_push(v, a);

    a.column = 2;
    a.value = -2.2;
    sb_push(v, a);

    a.column = 3;
    a.value = 3.5;
    sb_push(v, a);

    cr_assert_eq(sb_count(v) ,3);

    cr_assert_eq(v[0].column, 1);
    cr_assert_float_eq(v[0].value, 0.0, 1e-10);
    cr_assert_eq(v[0].cell, c);

    cr_assert_eq(v[1].column, 2);
    cr_assert_float_eq(v[1].value, -2.2, 1e-10);
    cr_assert_eq(v[1].cell, c);

    cr_assert_eq(v[2].column, 3);
    cr_assert_float_eq(v[2].value, 3.5, 1e-10);
    cr_assert_eq(v[2].cell, c);

    struct element b = sb_pop(v);
    cr_assert_eq(sb_count(v), 2);

    cr_assert_eq(b.column, 3);

    free(c);

}