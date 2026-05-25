////
//// Created by sachetto on 06/10/17.
////

#include <criterion/criterion.h>

#include "../3dparty/ini_parser/ini.h"
#include "../3dparty/sds/sds.h"
#include "../alg/grid/grid.h"
#include "../config/linear_system_solver_config.h"
#include "../utils/file_utils.h"

#include "../3dparty/stb_ds.h"

real_cpu calc_mse(const real_cpu *x, const real_cpu *xapp, int n) {

    real_cpu sum_sq = 0;

    for(int i = 0; i < n; ++i) {
        real_cpu err = x[i] - xapp[i];
        sum_sq += (err * err);
    }

    return sum_sq / n;
}

void test_solver(bool preconditioner, char *method_name, char *init_name, char *end_name, int nt, int version) {

    FILE *A = NULL;
    FILE *B = NULL;
    FILE *X = NULL;

    if(version == 1) {
        A = fopen("tests_bin/A1.txt", "r");
        B = fopen("tests_bin/B1.txt", "r");
        X = fopen("tests_bin/X1.txt", "r");
    } else if(version == 2) {
        A = fopen("tests_bin/A2.txt", "r");
        B = fopen("tests_bin/B2.txt", "r");
        X = fopen("tests_bin/X2.txt", "r");
    } else if(version == 3) {
        A = fopen("tests_bin/A3.txt", "r");
        B = fopen("tests_bin/B3.txt", "r");
        X = fopen("tests_bin/X3.txt", "r");
    }

    cr_assert(A);
    cr_assert(B);
    cr_assert(X);

    real_cpu error;

    struct grid *grid = new_grid();
    cr_assert(grid);

    construct_grid_from_file(grid, A, B);

#if defined(_OPENMP)
    omp_set_num_threads(nt);
    nt = omp_get_max_threads();
#endif

    struct config *linear_system_solver_config;

    linear_system_solver_config = alloc_and_init_config_data();
    linear_system_solver_config->main_function_name = method_name;

    if(init_name)
        linear_system_solver_config->init_function_name = init_name;
    if(end_name)
        linear_system_solver_config->end_function_name = end_name;

    shput(linear_system_solver_config->config_data, "tolerance", "1e-16");
    if(preconditioner)
        shput(linear_system_solver_config->config_data, "use_preconditioner", "yes");
    else
        shput(linear_system_solver_config->config_data, "use_preconditioner", "no");

    shput(linear_system_solver_config->config_data, "max_iterations", "200");

    uint32_t n_iter;

    init_config_functions(linear_system_solver_config, "./shared_libs/libdefault_linear_system_solver.so", "linear_system_solver");

    struct time_info ti = ZERO_TIME_INFO;

    CALL_INIT_LINEAR_SYSTEM(linear_system_solver_config, grid, false);
    ((linear_system_solver_fn *)linear_system_solver_config->main_function)(&ti, linear_system_solver_config, grid, grid->num_active_cells, grid->active_cells,
                                                                            &n_iter, &error);
    CALL_END_LINEAR_SYSTEM(linear_system_solver_config);

    uint64_t n_lines1;
    uint32_t n_lines2;

    real_cpu *x = read_octave_vector_file_to_array(X, &n_lines1);
    real_cpu *x_grid = grid_vector_to_array(grid, 'x', &n_lines2);

    // if(preconditioner)
    //     printf("MSE using %s with preconditioner and %d threads: %e\n", method_name, nt, calc_mse(x, x_grid, n_lines1));
    // else
    //     printf("MSE using %s without preconditioner and %d threads:: %e\n", method_name, nt, calc_mse(x, x_grid, n_lines1));

    cr_assert_eq(n_lines1, n_lines2);

    for(int i = 0; i < n_lines1; i++) {
        cr_assert_float_eq(x[i], x_grid[i], 1e-3, "Found %lf, Expected %lf.", x_grid[i], x[i]);
    }

    clean_and_free_grid(grid);
    fclose(A);
    fclose(B);
}
#ifndef COMPILE_SYCL

Test(solvers, cpu_cg_jacobi_preconditioner_1t) {
    test_solver(true, "cpu_conjugate_gradient", "init_cpu_conjugate_gradient", NULL, 1, 1);
}

Test(solvers, cpu_cg_no_preconditioner_1t) {
    test_solver(false, "cpu_conjugate_gradient", "init_cpu_conjugate_gradient", NULL, 1, 1);
}

#ifdef COMPILE_CUDA

Test(solvers, gpu_cg_jacobi_preconditioner_1t) {
    test_solver(true, "gpu_conjugate_gradient", "init_gpu_conjugate_gradient", "end_gpu_conjugate_gradient", 1, 1);
}

Test(solvers, gpu_cg_no_preconditioner_1t) {
    test_solver(false, "gpu_conjugate_gradient", "init_gpu_conjugate_gradient", "end_gpu_conjugate_gradient", 1, 1);
}

Test(solvers, bcg_gpu_no_preconditioner_1t) {
    test_solver(false, "gpu_biconjugate_gradient", "init_gpu_biconjugate_gradient", NULL, 1, 1);
}
#endif

Test(solvers, bcg_jacobi_preconditioner_1t) {
    test_solver(true, "biconjugate_gradient", NULL, NULL, 1, 1);
}

Test(solvers, jacobi_1t) {
    test_solver(false, "jacobi", NULL, NULL, 1, 1);
}

#if defined(_OPENMP)

Test(solvers, cg_jacobi_6t) {
    test_solver(true, "cpu_conjugate_gradient", "init_cpu_conjugate_gradient", NULL, 6, 1);
}

Test(solvers, cg_no_jacobi_6t) {
    test_solver(false, "cpu_conjugate_gradient", "init_cpu_conjugate_gradient", NULL, 6, 1);
}

Test(solvers, bcg_no_jacobi_6t) {
    test_solver(false, "biconjugate_gradient", NULL, NULL, 6, 1);
}

Test(solvers, jacobi_6t) {
    test_solver(false, "jacobi", NULL, NULL, 6, 1);
}

#endif

// ###########################################################################################

Test(solvers, cg_jacobi_1t_2) {
    test_solver(true, "cpu_conjugate_gradient", "init_cpu_conjugate_gradient", NULL, 1, 2);
}

Test(solvers, cg_no_jacobi_1t_2) {
    test_solver(false, "cpu_conjugate_gradient", "init_cpu_conjugate_gradient", NULL, 1, 2);
}

Test(solvers, bcg_jacobi_1t_2) {
    test_solver(true, "biconjugate_gradient", NULL, NULL, 1, 2);
}

Test(solvers, jacobi_1t_2) {
    test_solver(false, "jacobi", NULL, NULL, 1, 2);
}

#if defined(_OPENMP)

Test(solvers, cg_jacobi_6t_2) {
    test_solver(true, "cpu_conjugate_gradient", "init_cpu_conjugate_gradient", NULL, 6, 2);
}

Test(solvers, cg_no_jacobi_6t_2) {
    test_solver(false, "cpu_conjugate_gradient", "init_cpu_conjugate_gradient", NULL, 6, 2);
}

Test(solvers, bcg_no_jacobi_6t_2) {
    test_solver(false, "biconjugate_gradient", NULL, NULL, 6, 2);
}

Test(solvers, jacobi_6t_2) {
    test_solver(false, "jacobi", NULL, NULL, 6, 2);
}

#endif

// #######################################################################################################

Test(solvers, cg_jacobi_1t_3) {
    test_solver(true, "cpu_conjugate_gradient", "init_cpu_conjugate_gradient", NULL, 1, 3);
}

Test(solvers, cg_no_jacobi_1t_3) {
    test_solver(false, "cpu_conjugate_gradient", "init_cpu_conjugate_gradient", NULL, 1, 3);
}

Test(solvers, bcg_jacobi_1t_3) {
    test_solver(true, "biconjugate_gradient", NULL, NULL, 1, 3);
}

Test(solvers, jacobi_1t_3) {
    test_solver(false, "jacobi", NULL, NULL, 1, 3);
}

#if defined(_OPENMP)

Test(solvers, cg_jacobi_6t_3) {
    test_solver(true, "cpu_conjugate_gradient", "init_cpu_conjugate_gradient", NULL, 6, 3);
}

Test(solvers, cg_no_jacobi_6t_3) {
    test_solver(false, "cpu_conjugate_gradient", "init_cpu_conjugate_gradient", NULL, 6, 3);
}

Test(solvers, bcg_no_jacobi_6t_3) {
    test_solver(false, "biconjugate_gradient", NULL, NULL, 6, 3);
}

#ifdef COMPILE_CUDA
Test(solvers, bcg_gpu_no_preconditioner_1t_2) {
    test_solver(false, "gpu_biconjugate_gradient", "init_gpu_biconjugate_gradient", NULL, 1, 2);
}

Test(solvers, bcg_gpu_no_preconditioner_6t_2) {
    test_solver(false, "gpu_biconjugate_gradient", "init_gpu_biconjugate_gradient", NULL, 6, 2);
}

Test(solvers, amgx_6t_2) {
    test_solver(false, "amgx", "init_amgx", NULL, 6, 2);
}

#endif

#endif

#else
Test(solvers, gpu_cg_no_preconditioner_1t) {
    test_solver(false, "sycl_conjugate_gradient", "init_sycl_conjugate_gradient", "end_sycl_conjugate_gradient", 1, 1);
}
#endif
