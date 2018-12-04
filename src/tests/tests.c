////
//// Created by sachetto on 06/10/17.
////
#include "../monodomain/output_utils.h"
#include <criterion/criterion.h>
#include "../alg/grid/grid.h"
#include "../config/linear_system_solver_config.h"
#include "../config/domain_config.h"
#include "../config/save_mesh_config.h"
#include <signal.h>



double *read_octave_vector_file_to_array(FILE *vec_file, int *num_lines);

double **read_octave_mat_file_to_array(FILE *matrix_file, int *num_lines, int *nnz);

double calc_mse(const double *x, const double *xapp, int n) {

    double sum_sq = 0;

    for (int i = 0; i < n; ++i) {
        double err = x[i] - xapp[i];
        sum_sq += (err * err);
    }

    return sum_sq / n;
}


void construct_grid_from_file(struct grid *grid, FILE *matrix_a, FILE *vector_b) {

    uint32_t n_cells;
    int num_lines_m = 0;
    int num_lines_v = 0;
    int nnz = 0;

    double **matrix = read_octave_mat_file_to_array(matrix_a, &num_lines_m, &nnz);
    double *vector  = read_octave_vector_file_to_array(vector_b, &num_lines_v);

    cr_assert_eq (num_lines_m, num_lines_v);
    cr_assert (nnz);

    initialize_and_construct_grid(grid, 1.0, 1.0, 1.0);

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

    for (int i = 0; i < num_lines_m; i++) {
        free(matrix[i]);
    }

    free(matrix);
    free(vector);
}

double **read_octave_mat_file_to_array(FILE *matrix_file, int *num_lines, int *nnz) {
    const char *sep = " ";
    char *line_a = NULL;
    size_t len;
    int count;

    do {
        getline(&line_a, &len, matrix_file);
        sds *tmp = sdssplitlen(line_a, (int) strlen(line_a), sep, (int) strlen(sep), &count);
        if (count) {
            if (strcmp(tmp[1], "columns:") == 0) {
                (*num_lines) = atoi(tmp[2]);
            }
            if (strcmp(tmp[1], "nnz:") == 0) {
                (*nnz) = atoi(tmp[2]);
            }
        }
        sdsfreesplitres(tmp, count);
    } while ((line_a)[0] == '#');

    double **matrix = (double **) malloc(*num_lines * sizeof(double *));

    for (int i = 0; i < *num_lines; i++) {
        matrix[i] = (double *) calloc(*num_lines, sizeof(double));
    }

    int item_count = 0;
    int m_line, m_column;
    double m_value;

    while (item_count < *nnz) {

        sds *tmp = sdssplitlen(line_a, (int) strlen(line_a), sep, (int) strlen(sep), &count);
        if (tmp[0][0] != '\n') {
            m_line = atoi(tmp[0]);
            m_column = atoi(tmp[1]);
            m_value = atof(tmp[2]);

            matrix[m_line - 1][m_column - 1] = m_value;
        }
        sdsfreesplitres(tmp, count);

        item_count++;
        getline(&line_a, &len, matrix_file);
    }

    if (line_a)
        free(line_a);

    return matrix;
}

double *read_octave_vector_file_to_array(FILE *vec_file, int *num_lines) {

    ssize_t read;
    size_t len;
    char *line_b = NULL;
    int count;
    char *sep = " ";

    do {
        read = getline(&line_b, &len, vec_file);
        sds *tmp = sdssplitlen(line_b, (int) strlen(line_b), sep, (int) strlen(sep), &count);
        if (count) {
            if (strcmp(tmp[1], "rows:") == 0) {
                (*num_lines) = atoi(tmp[2]);
            }
        }
        sdsfreesplitres(tmp, count);
    } while ((line_b)[0] == '#');

    double *vector = (double *) malloc(*num_lines * sizeof(double));

    int item_count = 0;
    while ((item_count < *num_lines) && read) {
        sds *tmp = sdssplitlen(line_b, (int) strlen(line_b), sep, (int) strlen(sep), &count);

        if (tmp[0][0] != '\n') {
            vector[item_count] = atof(tmp[1]);
        }

        sdsfreesplitres(tmp, count);

        item_count++;
        read = getline(&line_b, &len, vec_file);
    }

    if (line_b)
        free(line_b);

    return vector;
}

void test_solver(bool preconditioner, char *method_name, int nt, int version) {

    FILE *A = NULL;
    FILE *B = NULL;
    FILE *X = NULL;

    if (version == 1) {
        A = fopen("src/tests/A.txt", "r");
        B = fopen("src/tests/B.txt", "r");
        X = fopen("src/tests/X.txt", "r");
    } else if (version == 2) {
        A = fopen("src/tests/A1.txt", "r");
        B = fopen("src/tests/B1.txt", "r");
        X = fopen("src/tests/X1.txt", "r");
    }

    cr_assert(A);
    cr_assert(B);
    cr_assert(X);

    double error;

    struct grid *grid = new_grid();
    cr_assert (grid);

    construct_grid_from_file(grid, A, B);


#if defined(_OPENMP)
    omp_set_num_threads(nt);
    nt = omp_get_max_threads();
#endif

    struct linear_system_solver_config *linear_system_solver_config;

    linear_system_solver_config = new_linear_system_solver_config();
    linear_system_solver_config->config_data.function_name = method_name;

    string_hash_insert_or_overwrite(linear_system_solver_config->config_data.config, "tolerance", "1e-16");
    if (preconditioner)
        string_hash_insert_or_overwrite(linear_system_solver_config->config_data.config, "use_preconditioner", "yes");
    else
        string_hash_insert_or_overwrite(linear_system_solver_config->config_data.config, "use_preconditioner", "no");

    string_hash_insert_or_overwrite(linear_system_solver_config->config_data.config, "max_iterations", "200");

    uint32_t n_iter;

    init_linear_system_solver_functions(linear_system_solver_config);


    linear_system_solver_config->solve_linear_system(linear_system_solver_config, grid, &n_iter, &error);

    int n_lines1;
    uint32_t n_lines2;

    double *x = read_octave_vector_file_to_array(X, &n_lines1);
    double *x_grid = grid_vector_to_array(grid, 'x', &n_lines2);

    if(preconditioner)
        printf("MSE using %s with preconditioner and %d threads: %e\n", method_name, nt, calc_mse(x, x_grid, n_lines1));
    else
        printf("MSE using %s without preconditioner and %d threads:: %e\n", method_name, nt, calc_mse(x, x_grid, n_lines1));

    cr_assert_eq (n_lines1, n_lines2);

    for (int i = 0; i < n_lines1; i++) {
        cr_assert_float_eq (x[i], x_grid[i], 1e-10);
    }

    clean_and_free_grid(grid);
    fclose(A);
    fclose(B);
}

int test_cuboid_mesh(double start_dx, double start_dy, double start_dz, char* side_length_x, char* side_length_y, char* side_length_z, bool save) {

    struct grid *grid = new_grid();
    struct domain_config *domain_config;

    domain_config = new_domain_config();

    domain_config->start_dx  = start_dx;
    domain_config->start_dy  = start_dy;
    domain_config->start_dz  = start_dz;

    domain_config->config_data.function_name = strdup("initialize_grid_with_cuboid_mesh");
    domain_config->domain_name = strdup("Test cuboid");

    string_hash_insert(domain_config->config_data.config, "side_length_x", side_length_x);
    string_hash_insert(domain_config->config_data.config, "side_length_y", side_length_y);
    string_hash_insert(domain_config->config_data.config, "side_length_z", side_length_z);

    init_domain_functions(domain_config);

    int success = domain_config->set_spatial_domain(domain_config, grid);

    if(!success ) {
        return 0;
    }

    order_grid_cells(grid);

    double sx = grid->side_length_x;
    double sy = grid->side_length_y;
    double sz = grid->side_length_z;

    double nx = sx / start_dx;
    double ny = sy / start_dy;
    double nz = sz / start_dz;

    struct cell_node *cell = grid->first_cell;

    double max_x = 0;
    double max_y = 0;
    double max_z = 0;

    while(cell) {

        if(cell->active) {
            if (cell->center_x > max_x) {
                max_x = cell->center_x;
            }

            if (cell->center_y > max_y) {
                max_y = cell->center_y;
            }

            if (cell->center_z > max_z) {
                max_z = cell->center_z;
            }
        }

        cell = cell->next;
    }

    if(save) {
        struct save_mesh_config *save_mesh_config = new_save_mesh_config();

        save_mesh_config->config_data.function_name = "save_as_vtu";
        save_mesh_config->out_dir_name = "./tests_bin";
        save_mesh_config->print_rate = 1;

        sds file_prefix = sdscatprintf(sdsempty(), "test_%lf_%lf_%lf_%s_%s_%s", start_dx, start_dy, start_dz,
                                       side_length_x, side_length_y, side_length_z);
        init_save_mesh_functions(save_mesh_config);
        string_hash_insert(save_mesh_config->config_data.config, "file_prefix", file_prefix);
        string_hash_insert(save_mesh_config->config_data.config, "compress", "yes");

        save_mesh_config->save_mesh(0.0, save_mesh_config, grid);
    }

    cr_assert_float_eq(max_x+(start_dx/2.0), atof(side_length_x), 1e-16);
    cr_assert_float_eq(max_y+(start_dy/2.0), atof(side_length_y), 1e-16);
    cr_assert_float_eq(max_z+(start_dz/2.0), atof(side_length_z), 1e-16);
    cr_assert_eq(nx*ny*nz, grid->num_active_cells);

    return 1;

}

Test (mesh_load, cuboid_mesh_100_100_100_1000_1000_1000) {
    int success  = test_cuboid_mesh(100, 100, 100, "1000", "1000", "1000", false);
    assert(success);
}

Test (mesh_load, cuboid_mesh_200_100_100_1000_1000_1000) {
    int success = test_cuboid_mesh(200, 100, 100, "1000", "1000", "1000", false);
    assert(success);
}

Test (mesh_load, cuboid_mesh_100_200_100_1000_1000_1000) {

    int success = test_cuboid_mesh(100, 200, 100, "1000", "1000", "1000", false);
    assert(success);
}

Test (mesh_load, cuboid_mesh_100_100_200_1000_1000_1000) {

    int success = test_cuboid_mesh(100, 100, 200, "1000", "1000", "1000", false);
    assert(success);
}

Test (mesh_load, cuboid_mesh_100_100_100_1000_1000_2000) {
    int success = test_cuboid_mesh(100, 100, 100, "1000", "1000", "2000", false);
    assert(success);
}

Test (mesh_load, cuboid_mesh_150_150_150_1500_1500_1500) {
    int success = test_cuboid_mesh(150, 150, 150, "1500", "1500", "1500", false);
    assert(success);
}


Test (mesh_load, cuboid_mesh_150_150_150_1500_1500_3000) {
    int success  = test_cuboid_mesh(150, 150, 150, "1500", "1500", "3000", false);
    assert(success);
}

Test (mesh_load, cuboid_mesh_300_150_150_1500_1500_3000) {
    int success = test_cuboid_mesh(300, 150, 150, "1500", "1500", "3000", false);
    assert(!success);
}

Test (solvers, cg_jacobi_1t) {
    test_solver(true, "conjugate_gradient", 1, 1);
}

Test (solvers, cg_no_jacobi_1t) {
    test_solver(false, "conjugate_gradient", 1, 1);
}

Test (solvers, bcg_jacobi_1t) {
    test_solver(true, "biconjugate_gradient", 1, 1);
}

Test (solvers, jacobi_1t) {
    test_solver(false, "jacobi", 1, 1);
}

#if defined(_OPENMP)

Test (solvers, cg_jacobi_6t) {
    test_solver(true, "conjugate_gradient", 6, 1);
}

Test (solvers, cg_no_jacobi_6t) {
    test_solver(false, "conjugate_gradient", 6, 1);
}


Test (solvers, bcg_no_jacobi_6t) {
    test_solver(false, "biconjugate_gradient", 6, 1);
}

Test (solvers, jacobi_6t) {
    test_solver(false, "jacobi", 6, 1);
}

#endif


Test (solvers, cg_jacobi_1t_2) {
    test_solver(true, "conjugate_gradient", 1, 2);
}

Test (solvers, cg_no_jacobi_1t_2) {
    test_solver(false, "conjugate_gradient", 1, 2);
}

Test (solvers, bcg_jacobi_1t_2) {
    test_solver(true, "biconjugate_gradient", 1, 2);
}

Test (solvers, jacobi_1t_2) {
    test_solver(false, "jacobi", 1, 2);
}

#if defined(_OPENMP)

Test (solvers, cg_jacobi_6t_2) {
    test_solver(true, "conjugate_gradient", 6, 2);
}

Test (solvers, cg_no_jacobi_6t_2) {
    test_solver(false, "conjugate_gradient", 6, 2);
}


Test (solvers, bcg_no_jacobi_6t_2) {
    test_solver(false, "biconjugate_gradient", 6, 2);
}

Test (solvers, jacobi_6t_2) {
    test_solver(false, "jacobi", 6, 2);
}

#endif


Test (utils, stretchy_buffer_int) {

    int *v = NULL;

            sb_reserve(v, 1);

    cr_assert_eq(sb_count(v), 0);
    cr_assert_eq(stb__sbm(v), 1);

            sb_push(v, 0);
            sb_push(v, 1);
            sb_push(v, 2);

    cr_assert_eq(sb_count(v), 3);

    cr_assert_eq(v[0], 0);
    cr_assert_eq(v[1], 1);
    cr_assert_eq(v[2], 2);
}

Test (utils, stretchy_buffer_float) {

    float *v = NULL;

            sb_reserve(v, 1);

    cr_assert_eq(sb_count(v), 0);
    cr_assert_eq(stb__sbm(v), 1);

            sb_push(v, 0);
            sb_push(v, 0.5);
            sb_push(v, 2.5);

    cr_assert_eq(sb_count(v), 3);

    cr_assert_float_eq(v[0], 0.0, 1e-10);
    cr_assert_float_eq(v[1], 0.5, 1e-10);
    cr_assert_float_eq(v[2], 2.5, 1e-10);
}

Test (utils, stretchy_buffer_double) {

    double *v = NULL;

            sb_reserve(v, 1);

    cr_assert_eq(sb_count(v), 0);
    cr_assert_eq(stb__sbm(v), 1);

            sb_push(v, 0);
            sb_push(v, 0.5);
            sb_push(v, 2.5);

    cr_assert_eq(sb_count(v), 3);

    cr_assert_float_eq(v[0], 0.0, 1e-10);
    cr_assert_float_eq(v[1], 0.5, 1e-10);
    cr_assert_float_eq(v[2], 2.5, 1e-10);
}

Test (utils, stretchy_buffer_element) {

    struct element *v = NULL;

            sb_reserve(v, 1);

    cr_assert_eq(sb_count(v), 0);
    cr_assert_eq(stb__sbm(v), 1);

    struct cell_node *c = new_cell_node();
    struct element a = {0, 1, c};

            sb_push(v, a);

    a.column = 2;
    a.value = -2.2;
            sb_push(v, a);

    a.column = 3;
    a.value = 3.5;
            sb_push(v, a);

    cr_assert_eq(sb_count(v), 3);

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
