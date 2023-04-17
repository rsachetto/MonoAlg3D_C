////
//// Created by sachetto on 06/10/17.
////
#include "../3dparty/sds/sds.h"
#include "../3dparty/stb_ds.h"
#include "../utils/file_utils.h"
#include "common.h"
#include <criterion/criterion.h>

static int check_ecg_file_equals(sds gold_ecg, sds tested_ecg, float tol) {

    string_array lines_gold = read_lines(gold_ecg);
    string_array lines_tested = read_lines(tested_ecg);

    cr_assert(lines_gold);
    cr_assert(lines_tested);

    ptrdiff_t n_lines_gold = arrlen(lines_gold);
    ptrdiff_t n_lines_tested = arrlen(lines_tested);

    cr_assert_eq(n_lines_gold, n_lines_tested, "%s %ld lines, %s %ld lines", gold_ecg, n_lines_gold, tested_ecg, n_lines_tested );

    for(int j = 0; j < n_lines_gold; j++) {

        int count_gold;
        int count_tested;

        sds *gold_values = sdssplit(lines_gold[j], " ", &count_gold);
        sds *tested_simulation_values = sdssplit(lines_tested[j], " ", &count_tested);

        cr_assert_eq(count_gold, count_tested);

        for(int k = 0; k < count_gold; k++) {
            real_cpu value_gold = strtod(gold_values[k], NULL);
            real_cpu value_tested = strtod(tested_simulation_values[k], NULL);
            cr_assert_float_eq(value_gold, value_tested, tol, "Found %lf, Expected %lf (error %e) on line %d of %s when comparing with %s", value_tested, value_gold, fabs(value_tested-value_gold), j+1, tested_ecg, gold_ecg);

        }

        sdsfreesplitres(gold_values, count_gold);
        sdsfreesplitres(tested_simulation_values, count_tested);
    }

    return 1;

}

static int check_output_equals(sds gold_output, sds tested_output, float tol) {

    string_array files_gold = list_files_from_dir(gold_output, "V_it_", "txt", NULL, true);
    string_array files_tested_sim = list_files_from_dir(tested_output, "V_it_", "txt", NULL, true);

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

            cr_assert_float_eq(value_gold, value_tested, tol, "Found %lf, Expected %lf (error %e) on line %d of %s when comparing with %s", value_tested, value_gold, fabs(value_tested-value_gold), j+1, full_path_tested, full_path_gold);
            sdsfreesplitres(gold_values, count_gold);
            sdsfreesplitres(tested_simulation_values, count_tested);
        }

        sdsfree(full_path_gold);
        sdsfree(full_path_tested);
    }

    return 1;
}

static int check_output_vtp_equals(sds gold_output, sds tested_output, float tol) {

    string_array files_gold = list_files_from_dir(gold_output, "V_it_", "vtp", NULL, true);
    string_array files_tested_sim = list_files_from_dir(tested_output, "V_it_", "vtp", NULL, true);

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

        int count_gold;
        int count_tested;

        sds *gold_values = sdssplit(lines_gold[5], ",", &count_gold);
        sds *tested_simulation_values = sdssplit(lines_tested[5], ",", &count_tested);

        cr_assert_eq(count_gold, count_tested);

        for(int k = 0; k < count_gold; k++) {
            real_cpu value_gold = strtod(gold_values[k], NULL);
            real_cpu value_tested = strtod(tested_simulation_values[k], NULL);
            cr_assert_float_eq(value_gold, value_tested, tol, "Found %lf, Expected %lf (error %e) on line %d of %s when comparing with %s", value_tested, value_gold, fabs(value_tested-value_gold), k+1, full_path_tested, full_path_gold);
        }

        sdsfreesplitres(gold_values, count_gold);
        sdsfreesplitres(tested_simulation_values, count_tested);
        sdsfree(full_path_gold);
        sdsfree(full_path_tested);
    }

    return 1;
}

#ifdef COMPILE_CUDA

Test(run_circle_simulation, gc_gpu_vs_cg_no_cpu) {

    char *out_dir_no_gpu_no_precond = "tests_bin/circle_cg_no_gpu_no_precond";
    char *out_dir_no_gpu_precond  = "tests_bin/circle_cg_no_gpu_precond";

    char *out_dir_gpu_no_precond = "tests_bin/circle_cg_gpu_no_precond";
    char *out_dir_gpu_precond  = "tests_bin/circle_cg_gpu_precond";


    struct user_options *options = load_options_from_file("example_configs/plain_mesh_with_fibrosis_and_border_zone_inside_circle_example_2cm.ini");
    options->final_time = 10.0;


    free(options->save_mesh_config->main_function_name);
    options->save_mesh_config->main_function_name = strdup("save_as_text_or_binary");

    free(options->save_mesh_config->init_function_name);
    free(options->save_mesh_config->end_function_name);

    options->save_mesh_config->init_function_name = NULL;
    options->save_mesh_config->end_function_name = NULL;

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
    int success = 1;

    success = run_simulation_with_config(options, out_dir_no_gpu_no_precond);
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

Test(run_gold_simulation, mesh_simulation) {
    char *out_dir_gold_gpu = "tests_bin/gold_sim_mesh_cg_gpu_ode_gpu/";
    char *out_dir_sim_gpu  = "tests_bin/sim_mesh_cg_gpu_ode_gpu/";

    char *ecg_gold = "tests_bin/gold_sim_mesh_cg_gpu_ode_gpu/ecg_gold.txt";
    char *ecg_sim  = "tests_bin/sim_mesh_cg_gpu_ode_gpu/ecg.txt";

    struct user_options *options = load_options_from_file("tests_bin/sim_mesh_cg_gpu_ode_gpu.ini");
    free(options->save_mesh_config->main_function_name);

    options->save_mesh_config->main_function_name = strdup("save_as_text_or_binary");
    shput_dup_value(options->save_mesh_config->config_data, "print_rate", "250");
    shput_dup_value(options->save_mesh_config->config_data, "file_prefix", "V");
    shput_dup_value(options->domain_config->config_data, strdup("seed"), "1526531136");
    shput_dup_value(options->calc_ecg_config->config_data, strdup("use_gpu"), "no");

    int success = 1;
    success = run_simulation_with_config(options, out_dir_sim_gpu);
    cr_assert(success);

    success &= check_output_equals(out_dir_gold_gpu, out_dir_sim_gpu, 10e-2f);
    success &= check_ecg_file_equals(ecg_gold, ecg_sim, 10e-2f);

    cr_assert(success);

    free_user_options(options);
}

Test(run_gold_simulation, mesh_simulation_ecg_gpu) {
    char *out_dir_gold_gpu = "tests_bin/gold_sim_mesh_cg_gpu_ode_gpu/";
    char *out_dir_sim_gpu  = "tests_bin/sim_mesh_cg_gpu_ode_gpu/";

    char *ecg_gold = "tests_bin/gold_sim_mesh_cg_gpu_ode_gpu/ecg_gold.txt";
    char *ecg_sim  = "tests_bin/sim_mesh_cg_gpu_ode_gpu/ecg.txt";

    struct user_options *options = load_options_from_file("tests_bin/sim_mesh_cg_gpu_ode_gpu.ini");
    free(options->save_mesh_config->main_function_name);

    options->save_mesh_config->main_function_name = strdup("save_as_text_or_binary");
    shput_dup_value(options->save_mesh_config->config_data, "print_rate", "250");
    shput_dup_value(options->save_mesh_config->config_data, "file_prefix", "V");
    shput_dup_value(options->domain_config->config_data, strdup("seed"), "1526531136");

    int success = 1;
    success = run_simulation_with_config(options, out_dir_sim_gpu);
    cr_assert(success);

    success &= check_output_equals(out_dir_gold_gpu, out_dir_sim_gpu, 10e-2f);
    success &= check_ecg_file_equals(ecg_gold, ecg_sim, 10e-2f);

    cr_assert(success);

    free_user_options(options);
}

Test(run_gold_simulation, purkinje_simulation) {
    char *out_dir_gold = "tests_bin/gold_purkinje/";

    char *out_dir_sim_cpu  = "tests_bin/purkinje_ode_gpu/";
    char *out_dir_sim_gpu  = "tests_bin/purkinje_ode_cpu/";

    struct user_options *options = load_options_from_file("tests_bin/sim_purkinje_gold.ini");

    int success = 1;
    success = run_simulation_with_config(options, out_dir_sim_gpu);
    cr_assert(success);

    options->gpu = false;

    success = run_simulation_with_config(options, out_dir_sim_cpu);
    cr_assert(success);

    success = check_output_vtp_equals(out_dir_gold, out_dir_sim_cpu, 5e-2f);
    success &= check_output_vtp_equals(out_dir_gold, out_dir_sim_gpu, 5e-2f);

    cr_assert(success);

    free_user_options(options);
}
#endif
