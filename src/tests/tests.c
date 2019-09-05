////
//// Created by sachetto on 06/10/17.
////
#include <criterion/criterion.h>
#include "../alg/grid/grid.h"
#include "../config/linear_system_solver_config.h"
#include "../config/domain_config.h"
#include "../config/save_mesh_config.h"
#include "../config/config_parser.h"
#include "../utils/file_utils.h"
#include "../ini_parser/ini.h"
#include "../string/sds.h"
#include <signal.h>

//#define STB_DS_IMPLEMENTATION
#include "../single_file_libraries/stb_ds.h"
#include "../config_helpers/config_helpers.h"


real_cpu *read_octave_vector_file_to_array(FILE *vec_file, int *num_lines);

real_cpu **read_octave_mat_file_to_array(FILE *matrix_file, int *num_lines, int *nnz);

real_cpu calc_mse(const real_cpu *x, const real_cpu *xapp, int n) {

    real_cpu sum_sq = 0;

    for (int i = 0; i < n; ++i) {
        real_cpu err = x[i] - xapp[i];
        sum_sq += (err * err);
    }

    return sum_sq / n;
}

void test_solver(bool preconditioner, char *method_name, char *init_name, char *end_name, int nt, int version) {

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

    real_cpu error;

    struct grid *grid = new_grid();
    cr_assert (grid);

    construct_grid_from_file(grid, A, B);


#if defined(_OPENMP)
    omp_set_num_threads(nt);
    nt = omp_get_max_threads();
#endif

    struct config *linear_system_solver_config;

    linear_system_solver_config = alloc_and_init_config_data();
    linear_system_solver_config->main_function_name = method_name;

    if(init_name)  linear_system_solver_config->init_function_name = init_name;
    if(end_name)  linear_system_solver_config->end_function_name = end_name;

    shput(linear_system_solver_config->config_data, "tolerance", "1e-16");
    if (preconditioner)
        shput(linear_system_solver_config->config_data, "use_preconditioner", "yes");
    else
        shput(linear_system_solver_config->config_data, "use_preconditioner", "no");

    shput(linear_system_solver_config->config_data, "max_iterations", "200");

    uint32_t n_iter;

    init_config_functions(linear_system_solver_config, "./shared_libs/libdefault_linear_system_solver.so", "linear_system_solver");


    CALL_INIT_LINEAR_SYSTEM(linear_system_solver_config, grid);
    ((linear_system_solver_fn*)linear_system_solver_config->main_function)(linear_system_solver_config, grid, &n_iter, &error);
    CALL_END_LINEAR_SYSTEM(linear_system_solver_config);

    int n_lines1;
    uint32_t n_lines2;

    real_cpu *x = read_octave_vector_file_to_array(X, &n_lines1);
    real_cpu *x_grid = grid_vector_to_array(grid, 'x', &n_lines2);

    if(preconditioner)
        printf("MSE using %s with preconditioner and %d threads: %e\n", method_name, nt, calc_mse(x, x_grid, n_lines1));
    else
        printf("MSE using %s without preconditioner and %d threads:: %e\n", method_name, nt, calc_mse(x, x_grid, n_lines1));

    cr_assert_eq (n_lines1, n_lines2);

    for (int i = 0; i < n_lines1; i++) {
        cr_assert_float_eq (x[i], x_grid[i], 1e-3, "Found %lf, Expected %lf.", x_grid[i], x[i]);
    }

    clean_and_free_grid(grid);
    fclose(A);
    fclose(B);
}

int test_perlin_mesh(char* mesh_file, char *start_dx, char* side_length_x, bool save, bool compress,  bool binary) {

    struct grid *grid = new_grid();
    struct config *domain_config;

    domain_config = alloc_and_init_config_data();

    domain_config->main_function_name = strdup("set_perlin_square_mesh");

    shput_dup_value(domain_config->config_data, "name", "Perlin mesh");

    shput(domain_config->config_data, "side_length", side_length_x);
    shput(domain_config->config_data, "mesh_file",   mesh_file);
    shput(domain_config->config_data, "start_discretization",  start_dx);

    init_config_functions(domain_config, "./shared_libs/libdefault_domains.so", "domain");

    int success = ((set_spatial_domain_fn *)domain_config->main_function)(domain_config, grid);

    if(!success ) {
        return 0;
    }

    if(save) {
        struct config *save_mesh_config = alloc_and_init_config_data();

        save_mesh_config->main_function_name = "save_as_vtu";

        shput_dup_value(save_mesh_config->config_data, "output_dir", "./tests_bin");
        shput_dup_value(save_mesh_config->config_data, "print_rate", "1");

        sds file_prefix = sdscatprintf(sdsempty(), "test_perlin");
        init_config_functions(save_mesh_config, "shared_libs/libdefault_save_mesh.so", "save_result");

        shput(save_mesh_config->config_data, "file_prefix", file_prefix);
        if(compress)
        shput(save_mesh_config->config_data, "compress", "yes");
        else if(binary)
        shput(save_mesh_config->config_data, "binary", "yes");

        ((save_mesh_fn*)save_mesh_config->main_function)(save_mesh_config, grid, 0, 0.0, 0.0, 0.0);

    }

    return 1;

}

int test_cuboid_mesh(char *start_dx, char* start_dy, char* start_dz, char* side_length_x, char* side_length_y, char* side_length_z, bool save, bool compress,  bool binary, int id) {

    struct grid *grid = new_grid();
    struct config *domain_config;

    domain_config = alloc_and_init_config_data();

    shput_dup_value(domain_config->config_data, "start_dx", start_dx);
    shput_dup_value(domain_config->config_data, "start_dy", start_dy);
    shput_dup_value(domain_config->config_data, "start_dz", start_dz);

    domain_config->main_function_name = strdup("initialize_grid_with_cuboid_mesh");
    shput_dup_value(domain_config->config_data, "name", "Test cuboid");

    shput(domain_config->config_data, "side_length_x", strdup(side_length_x));
    shput(domain_config->config_data, "side_length_y", strdup(side_length_y));
    shput(domain_config->config_data, "side_length_z", strdup(side_length_z));

     init_config_functions(domain_config, "./shared_libs/libdefault_domains.so", "domain");

    int success = ((set_spatial_domain_fn*)domain_config->main_function)(domain_config, grid);

    if(!success ) {
        clean_and_free_grid(grid);
        free_config_data(domain_config);
        return 0;
    }

    order_grid_cells(grid);

    real_cpu sx = grid->mesh_side_length.x;
    real_cpu sy = grid->mesh_side_length.y;
    real_cpu sz = grid->mesh_side_length.z;

    real_cpu start_dx_value = strtod(start_dx, NULL);
    real_cpu start_dy_value = strtod(start_dy, NULL);
    real_cpu start_dz_value = strtod(start_dz, NULL);

    real_cpu nx = sx / start_dx_value;
    real_cpu ny = sy / start_dy_value;
    real_cpu nz = sz / start_dz_value;

    struct cell_node *cell = grid->first_cell;

    real_cpu max_x = 0;
    real_cpu max_y = 0;
    real_cpu max_z = 0;

    while(cell) {

        if(cell->active) {
            if (cell->center.x > max_x) {
                max_x = cell->center.x;
            }

            if (cell->center.y > max_y) {
                max_y = cell->center.y;
            }

            if (cell->center.z > max_z) {
                max_z = cell->center.z;
            }
        }

        cell = cell->next;
    }

    if(save) {
        struct config *save_mesh_config = alloc_and_init_config_data();

        save_mesh_config->main_function_name = strdup("save_as_vtu");
        shput_dup_value(save_mesh_config->config_data, "output_dir", "./tests_bin");
        shput_dup_value(save_mesh_config->config_data, "print_rate", "1");

        sds file_prefix = sdscatprintf(sdsempty(), "test_%s_%s_%s_%s_%s_%s_%d", start_dx, start_dy, start_dz,
                                       side_length_x, side_length_y, side_length_z, id);

        init_config_functions(save_mesh_config, "shared_libs/libdefault_save_mesh.so", "save_result");

        shput(save_mesh_config->config_data, "file_prefix", strdup(file_prefix));

        if(compress) {
            shput(save_mesh_config->config_data, "compress", strdup("yes"));
        }
        else if(binary) {
            shput(save_mesh_config->config_data, "binary", strdup("yes"));
        }

        shput(save_mesh_config->config_data, "save_pvd", strdup("no"));

        ((save_mesh_fn *)save_mesh_config->main_function)(save_mesh_config, grid, 0, 0.0, 0.0, 0.0);

        free_config_data(save_mesh_config);
        sdsfree(file_prefix);

    }

    cr_assert_float_eq(max_x+(start_dx_value/2.0), atof(side_length_x), 1e-16);
    cr_assert_float_eq(max_y+(start_dy_value/2.0), atof(side_length_y), 1e-16);
    cr_assert_float_eq(max_z+(start_dz_value/2.0), atof(side_length_z), 1e-16);
    cr_assert_eq(nx*ny*nz, grid->num_active_cells);

    clean_and_free_grid(grid);
    free_config_data(domain_config);

    return 1;

}

struct user_options *load_options_from_file(char *config_file) {
    // Here we parse the config file

    struct user_options *options = new_user_options();

    if(config_file) {
        options->config_file = strdup(config_file);

        if(ini_parse(options->config_file, parse_config_file, options) < 0) {
            fprintf(stderr, "Error: Can't load the config file %s\n", options->config_file);
            return NULL;
        }

        return options;
    }

    return NULL;

}

int run_simulation_with_config(struct user_options *options, char *out_dir) {

    struct grid *the_grid;
    the_grid = new_grid();

    struct monodomain_solver *monodomain_solver;
    monodomain_solver = new_monodomain_solver();

    struct ode_solver *ode_solver;
    ode_solver = new_ode_solver();

    no_stdout = true;

    shput_dup_value(options->save_mesh_config->config_data, "output_dir", out_dir);

    char *out_dir_name = NULL;
    GET_PARAMETER_VALUE_CHAR_OR_USE_DEFAULT(out_dir_name, options->save_mesh_config->config_data, "output_dir");

    // Create the output dir and the logfile
    if(out_dir_name) {
        remove_directory(out_dir_name);
        create_dir(out_dir_name);

        sds buffer_log = sdsempty();
        buffer_log = sdscatfmt(buffer_log, "%s/outputlog.txt", out_dir_name);
        open_logfile(buffer_log);

        sdsfree(buffer_log);

    }
    else {
        return 0;
    }

    configure_ode_solver_from_options(ode_solver, options);
    configure_monodomain_solver_from_options(monodomain_solver, options);
    configure_grid_from_options(the_grid, options);

    #ifndef COMPILE_CUDA
    if(ode_solver->gpu) {
        print_to_stdout_and_file("Cuda runtime not found in this system. Fallbacking to CPU solver!!\n");
        ode_solver->gpu = false;
    }
    #endif

    int np = monodomain_solver->num_threads;

    if(np == 0)
        np = 1;

    #if defined(_OPENMP)
    omp_set_num_threads(np);
    #endif

    solve_monodomain(monodomain_solver, ode_solver, the_grid, options);

    clean_and_free_grid(the_grid);
    free_ode_solver(ode_solver);

    free(monodomain_solver);

    close_logfile();

    return 1;
}

int check_output_equals(const sds gold_output, const sds tested_output, float tol) {

    string_array files_gold = list_files_from_dir(gold_output, "V_it_");
    string_array files_tested_sim = list_files_from_dir(tested_output, "V_it_");

    cr_assert(files_gold != NULL);
    cr_assert(files_tested_sim != NULL);

    ptrdiff_t n_files_gold = arrlen(files_gold);
    ptrdiff_t n_files_tested = arrlen(files_tested_sim);

    cr_assert_eq(n_files_gold, n_files_tested);

    for(int i = 0; i < n_files_gold; i++) {
        sds full_path_gold = sdsempty();
        full_path_gold = sdscatprintf(full_path_gold, "%s/", gold_output);
        full_path_gold = sdscat(full_path_gold, files_gold[i]);

        sds full_path_tested = sdsempty();
        full_path_tested = sdscatprintf(full_path_tested, "%s/", tested_output);
        full_path_tested = sdscat(full_path_tested, files_tested_sim[i]);

        string_array lines_gold = read_lines(full_path_gold);
        string_array lines_tested = read_lines(full_path_tested);

        cr_assert(lines_gold);
        cr_assert(lines_tested);

        ptrdiff_t n_lines_gold = arrlen(lines_gold);
        ptrdiff_t n_lines_tested = arrlen(lines_tested);

        cr_assert_eq(n_lines_gold, n_lines_tested, "%s %ld lines, %s %ld lines", full_path_gold, n_lines_gold, full_path_tested, n_lines_tested );

        for(int j = 0; j < n_lines_gold; j++) {

            int count_gold;
            int count_tested;

            sds *gold_values = sdssplit(lines_gold[j], ",", &count_gold);
            sds *tested_simulation_values = sdssplit(lines_tested[j], ",", &count_tested);

            cr_assert_eq(count_gold, count_tested);

            for(int k = 0; k < count_gold-1; k++) {
                real_cpu value_gold = strtod(gold_values[k], NULL);
                real_cpu value_tested = strtod(tested_simulation_values[k], NULL);
                cr_assert_eq(value_gold, value_tested);
            }

            real_cpu value_gold = strtod(gold_values[count_gold-1], NULL);
            real_cpu value_tested = strtod(tested_simulation_values[count_gold-1], NULL);

            cr_assert_float_eq(value_gold, value_tested, tol, "Found %lf, Expected %lf (error %e) on line %d of %s", value_tested, value_gold, fabs(value_tested-value_gold), i+1, full_path_tested);
            sdsfreesplitres(gold_values, count_gold);
            sdsfreesplitres(tested_simulation_values, count_tested);
        }

        sdsfree(full_path_gold);
        sdsfree(full_path_tested);
    }

    return 1;
}

int compare_two_binary_files(FILE *fp1, FILE *fp2)
{
    char ch1, ch2;
    int flag = 0;

    while (((ch1 = fgetc(fp1)) != EOF) &&((ch2 = fgetc(fp2)) != EOF))
    {
        /*
          * character by character comparision
          * if equal then continue by comparing till the end of files
          */
        if (ch1 == ch2)
        {
            flag = 1;
            continue;
        }
            /*
              * If not equal then returns the byte position
              */
        else
        {
            fseek(fp1, -1, SEEK_CUR);
            flag = 0;
            break;
        }
    }

    if (flag == 0)
    {
       return ftell(fp1)+1;
    }
    else
    {
        return -1;
    }
}

Test (mesh_load, cuboid_mesh_100_100_100_1000_1000_1000) {
    int success  = test_cuboid_mesh("100", "100", "100", "1000", "1000", "1000", false, false, false, 0);
    cr_assert(success);
}

Test (mesh_load, cuboid_mesh_200_100_100_1000_1000_1000) {
    int success = test_cuboid_mesh("200", "100", "100", "1000", "1000", "1000", false, false, false, 0);
    cr_assert(success);
}

Test (mesh_load, cuboid_mesh_100_200_100_1000_1000_1000) {

    int success = test_cuboid_mesh("100", "200", "100", "1000", "1000", "1000", false, false, false, 0);
    cr_assert(success);
}

Test (mesh_load, cuboid_mesh_100_100_200_1000_1000_1000) {
    int success = test_cuboid_mesh("100", "100", "200", "1000", "1000", "1000", false, false, false, 0);
    cr_assert(success);

}

Test (mesh_load, cuboid_mesh_100_100_100_1000_1000_2000) {
    int success = test_cuboid_mesh("100", "100", "100", "1000", "1000", "2000", false, false, false, 0);
    cr_assert(success);
}

Test (mesh_load, cuboid_mesh_150_150_150_1500_1500_1500) {
    int success = test_cuboid_mesh("150", "150", "150", "1500", "1500", "1500", false, false, false, 0);
    cr_assert(success);
}


Test (mesh_load, cuboid_mesh_150_150_150_1500_1500_3000) {
    int success  = test_cuboid_mesh("150", "150", "150", "1500", "1500", "3000", false, false, false, 0);
    cr_assert(success);
}
Test (mesh_load, cuboid_mesh_300_150_150_1500_1500_3000) {
    int success = test_cuboid_mesh("300", "150", "150", "1500", "1500", "3000", false, false, false, 0);
    cr_assert(!success);
}

Test (mesh_load_and_check_save, cuboid_mesh_100_100_200_1000_1000_1000_check_binary) {

    int success = test_cuboid_mesh("100", "100", "200", "1000", "1000", "1000", true, false, true, 1);
    cr_assert(success);

    FILE *f1 = fopen("tests_bin/test_100.000000_100.000000_200.000000_1000_1000_1000_1_it_0.vtu", "r");
    FILE *f2 = fopen("tests_bin/gold_vtu_mesh_binary.vtu", "r");

    success = compare_two_binary_files(f1, f2);

    fclose(f1);
    fclose(f2);

    cr_assert(success == -1);
}

Test (mesh_load_and_check_save, cuboid_mesh_100_100_200_1000_1000_1000_check_compressed) {

    int success = test_cuboid_mesh("100", "100", "200", "1000", "1000", "1000", true, true, false, 2);
    cr_assert(success);

    FILE *f1 = fopen("tests_bin/test_100.000000_100.000000_200.000000_1000_1000_1000_2_it_0.vtu", "r");
    FILE *f2 = fopen("tests_bin/gold_vtu_mesh_compressed.vtu", "r");

    success = compare_two_binary_files(f1, f2);

    fclose(f1);
    fclose(f2);

    cr_assert(success == -1);
}

Test (mesh_load_and_check_save, cuboid_mesh_100_100_200_1000_1000_1000_check_plain) {

    int success = test_cuboid_mesh("100", "100", "200", "1000", "1000", "1000", true, false, false, 3);
    cr_assert(success);

    FILE *f1 = fopen("tests_bin/test_100.000000_100.000000_200.000000_1000_1000_1000_3_it_0.vtu", "r");
    FILE *f2 = fopen("tests_bin/gold_vtu_mesh.vtu", "r");

    cr_assert(f1);
    cr_assert(f2);

    success = compare_two_binary_files(f1, f2);

    fclose(f1);
    fclose(f2);

    cr_assert(success == -1);
}

Test (solvers, cpu_cg_jacobi_1t) {
    test_solver(true, "cpu_conjugate_gradient", NULL, NULL, 1, 1);
}

Test (solvers, cpu_cg_no_jacobi_1t) {
    test_solver(false, "cpu_conjugate_gradient", NULL, NULL, 1, 1);
}

Test (solvers, gpu_cg_jacobi_1t) {
    test_solver(true, "gpu_conjugate_gradient", "init_gpu_conjugate_gradient", "end_gpu_conjugate_gradient", 1, 1);
}

Test (solvers, gpu_cg_no_jacobi_1t) {
    test_solver(false, "gpu_conjugate_gradient", "init_gpu_conjugate_gradient", "end_gpu_conjugate_gradient",  1, 1);
}

Test (solvers, bcg_jacobi_1t) {
    test_solver(true, "biconjugate_gradient", NULL, NULL, 1, 1);
}

Test (solvers, jacobi_1t) {
    test_solver(false, "jacobi", NULL, NULL, 1, 1);
}

#if defined(_OPENMP)

Test (solvers, cg_jacobi_6t) {
    test_solver(true, "conjugate_gradient", NULL, NULL, 6, 1);
}

Test (solvers, cg_no_jacobi_6t) {
    test_solver(false, "conjugate_gradient", NULL, NULL, 6, 1);
}


Test (solvers, bcg_no_jacobi_6t) {
    test_solver(false, "biconjugate_gradient", NULL, NULL, 6, 1);
}

Test (solvers, jacobi_6t) {
    test_solver(false, "jacobi", NULL, NULL, 6, 1);
}

#endif


Test (solvers, cg_jacobi_1t_2) {
    test_solver(true, "conjugate_gradient", NULL, NULL, 1, 2);
}

Test (solvers, cg_no_jacobi_1t_2) {
    test_solver(false, "conjugate_gradient", NULL, NULL, 1, 2);
}

Test (solvers, bcg_jacobi_1t_2) {
    test_solver(true, "biconjugate_gradient", NULL, NULL, 1, 2);
}

Test (solvers, jacobi_1t_2) {
    test_solver(false, "jacobi", NULL, NULL, 1, 2);
}

#if defined(_OPENMP)

Test (solvers, cg_jacobi_6t_2) {
    test_solver(true, "conjugate_gradient", NULL, NULL, 6, 2);
}

Test (solvers, cg_no_jacobi_6t_2) {
    test_solver(false, "conjugate_gradient", NULL, NULL, 6, 2);
}


Test (solvers, bcg_no_jacobi_6t_2) {
    test_solver(false, "biconjugate_gradient", NULL, NULL, 6, 2);
}

Test (solvers, jacobi_6t_2) {
    test_solver(false, "jacobi", NULL, NULL, 6, 2);
}

#endif
Test(run_gold_simulation, gpu_no_adapt) {

    printf("Running simulation for testing\n");

    char *out_dir  = "tests_bin/gold_tmp_no_gpu";

    struct user_options *options = load_options_from_file("example_configs/gold_simulation_no_adapt.ini");

    int success = run_simulation_with_config(options, out_dir);
    cr_assert(success);

    sds gold_dir = sdsnew("tests_bin/gold_simulation_no_adapt_gpu/");
    sds tested_simulation_dir = sdsnew(out_dir);

    success = check_output_equals(gold_dir, tested_simulation_dir, 1e-3f);
    cr_assert(success);

    free_user_options(options);
}

Test(run_gold_simulation, gpu_no_adapt_cg_gpu) {

    printf("Running simulation for testing\n");

    struct user_options *options = load_options_from_file("example_configs/gold_simulation_no_adapt_cg_gpu.ini");

    char *out_dir  = "tests_bin/gold_tmp_gpu";
    int success = run_simulation_with_config(options, out_dir);
    cr_assert(success);

    sds gold_dir = sdsnew("tests_bin/gold_simulation_no_adapt_gpu/");
    sds tested_simulation_dir = sdsnew(out_dir);

    success = check_output_equals(gold_dir, tested_simulation_dir, 5e-2f);
    cr_assert(success);
    free_user_options(options);

}
#ifdef COMPILE_CUDA
Test(run_circle_simulation, gc_gpu_vs_cg_no_cpu) {

    char *out_dir_no_gpu_no_precond  = "tests_bin/circle_cg_no_gpu_no_precond";
    char *out_dir_no_gpu_precond  = "tests_bin/circle_cg_no_gpu_precond";

    char *out_dir_gpu_no_precond  = "tests_bin/circle_cg_gpu_no_precond";
    char *out_dir_gpu_precond  = "tests_bin/circle_cg_gpu_precond";


    struct user_options *options = load_options_from_file("example_configs/plain_mesh_with_fibrosis_and_border_zone_inside_circle_example_2cm.ini");
    options->final_time = 10.0;

    free(options->save_mesh_config->main_function_name);

    options->save_mesh_config->main_function_name = strdup("save_as_text_or_binary");
    shput_dup_value(options->save_mesh_config->config_data, "print_rate", "50");
    shput_dup_value(options->save_mesh_config->config_data, "file_prefix", "V");

    shput_dup_value(options->domain_config->config_data, strdup("seed"), "150");

    free(options->linear_system_solver_config->main_function_name);
    free(options->linear_system_solver_config->init_function_name);
    free(options->linear_system_solver_config->end_function_name);

    options->linear_system_solver_config->main_function_name = strdup("conjugate_gradient");
    options->linear_system_solver_config->init_function_name = strdup("init_conjugate_gradient");
    options->linear_system_solver_config->end_function_name = strdup("end_conjugate_gradient");

    shput_dup_value(options->linear_system_solver_config->config_data, "use_gpu", "false");
    shput_dup_value(options->linear_system_solver_config->config_data, "use_preconditioner", "false");

    int success = run_simulation_with_config(options, out_dir_no_gpu_no_precond);
    cr_assert(success);

    shput_dup_value(options->linear_system_solver_config->config_data, "use_preconditioner", "true");
    success = run_simulation_with_config(options, out_dir_no_gpu_precond);
    cr_assert(success);

    shput_dup_value(options->linear_system_solver_config->config_data, "use_gpu", "true");
    shput_dup_value(options->linear_system_solver_config->config_data, "use_preconditioner", "false");

    success = run_simulation_with_config(options, out_dir_gpu_no_precond);
    cr_assert(success);

    shput_dup_value(options->linear_system_solver_config->config_data, "use_preconditioner", "true");

    success = run_simulation_with_config(options, out_dir_gpu_precond);
    cr_assert(success);

    success = check_output_equals(out_dir_no_gpu_no_precond, out_dir_no_gpu_precond, 5e-2f);
    success &= check_output_equals(out_dir_no_gpu_no_precond, out_dir_gpu_no_precond, 5e-2f);
    success &= check_output_equals(out_dir_no_gpu_no_precond, out_dir_gpu_precond, 5e-2f);
    success &= check_output_equals(out_dir_no_gpu_precond, out_dir_gpu_no_precond, 5e-2f);
    success &= check_output_equals(out_dir_no_gpu_precond, out_dir_gpu_precond, 5e-2f);
    success &= check_output_equals(out_dir_gpu_precond, out_dir_gpu_no_precond, 5e-2f);

    cr_assert(success);



    free_user_options(options);
}
#endif


Test (utils, arr_int) {

    int *v = NULL;

    arrsetcap(v, 1);

    cr_assert_eq(arrlen(v), 0);
    cr_assert_geq(arrcap(v), 1);

    arrput(v, 0);
    arrput(v, 1);
    arrput(v, 2);

    cr_assert_eq(arrlen(v), 3);

    cr_assert_eq(v[0], 0);
    cr_assert_eq(v[1], 1);
    cr_assert_eq(v[2], 2);
}

Test (utils, arr_float) {

    float *v = NULL;

    arrsetcap(v, 1);

    cr_assert_eq(arrlen(v), 0);
    cr_assert_geq(arrcap(v), 1);

    arrput(v, 0);
    arrput(v, 0.5);
    arrput(v, 2.5);

    cr_assert_eq(arrlen(v), 3);

    cr_assert_float_eq(v[0], 0.0, 1e-10);
    cr_assert_float_eq(v[1], 0.5, 1e-10);
    cr_assert_float_eq(v[2], 2.5, 1e-10);
}

Test (utils, arr_double) {

    real_cpu *v = NULL;

    arrsetcap(v, 1);

    cr_assert_eq(arrlen(v), 0);
    cr_assert_geq(arrcap(v), 1);

    arrput(v, 0);
    arrput(v, 0.5);
    arrput(v, 2.5);

    cr_assert_eq(arrlen(v), 3);

    cr_assert_float_eq(v[0], 0.0, 1e-10);
    cr_assert_float_eq(v[1], 0.5, 1e-10);
    cr_assert_float_eq(v[2], 2.5, 1e-10);

}

Test (utils, arr_element) {

    struct element *v = NULL;

    arrsetcap(v, 1);

    cr_assert_eq(arrlen(v), 0);
    cr_assert_geq(arrcap(v), 1);

    struct cell_node *c = new_cell_node();
    struct element a = {'a', 0, 1, c};

    arrput(v, a);

    a.column = 2;
    a.value = -2.2;
    arrput(v, a);

    a.column = 3;
    a.value = 3.5;
    arrput(v, a);

    cr_assert_eq(arrlen(v), 3);

    cr_assert_eq(v[0].column, 1);
    cr_assert_float_eq(v[0].value, 0.0, 1e-10);
    cr_assert_eq(v[0].cell, c);

    cr_assert_eq(v[1].column, 2);
    cr_assert_float_eq(v[1].value, -2.2, 1e-10);
    cr_assert_eq(v[1].cell, c);

    cr_assert_eq(v[2].column, 3);
    cr_assert_float_eq(v[2].value, 3.5, 1e-10);
    cr_assert_eq(v[2].cell, c);

    struct element b = arrpop(v);
    cr_assert_eq(arrlen(v), 2);

    cr_assert_eq(b.column, 3);

    free(c);
}
