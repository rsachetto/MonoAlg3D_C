////
//// Created by sachetto on 06/10/17.
////
#include <criterion/criterion.h>
#include <signal.h>

#include "../alg/grid/grid.h"
#include "../config/linear_system_solver_config.h"
#include "../config/domain_config.h"
#include "../config/save_mesh_config.h"
#include "../config/config_parser.h"
#include "../utils/file_utils.h"
#include "../3dparty/ini_parser/ini.h"
#include "../3dparty/sds/sds.h"
#include "../logger/logger.h"

#include "../3dparty/stb_ds.h"
#include "../config_helpers/config_helpers.h"


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

    shput_dup_value(options->save_mesh_config->config_data, "output_dir", out_dir);

    char *out_dir_name = NULL;
    GET_PARAMETER_STRING_VALUE_OR_USE_DEFAULT(out_dir_name, options->save_mesh_config->config_data, "output_dir");

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
        log_to_stdout_and_file("Cuda runtime not found in this system. Fallbacking to CPU solver!!\n");
        ode_solver->gpu = false;
    }
#endif

    int np = monodomain_solver->num_threads;

    if(np == 0)
        np = 1;

#if defined(_OPENMP)
    omp_set_num_threads(np);
#endif

    set_no_stdout(true);

    solve_monodomain(monodomain_solver, ode_solver, the_grid, options);

    clean_and_free_grid(the_grid);
    free_ode_solver(ode_solver);
    free(monodomain_solver);

    close_logfile();

    return 1;
}

int check_output_equals(const sds gold_output, const sds tested_output, float tol) {

    string_array files_gold = list_files_from_dir(gold_output, "V_it_", "txt", true);
    string_array files_tested_sim = list_files_from_dir(tested_output, "V_it_", "txt", true);

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

            cr_assert_float_eq(value_gold, value_tested, tol, "Found %lf, Expected %lf (error %e) on line %d of %s when comparing with %s", value_tested, value_gold, fabs(value_tested-value_gold), i+1, full_path_tested, full_path_gold);
            sdsfreesplitres(gold_values, count_gold);
            sdsfreesplitres(tested_simulation_values, count_tested);
        }

        sdsfree(full_path_gold);
        sdsfree(full_path_tested);
    }

    return 1;
}

#ifdef COMPILE_CUDA
Test(run_circle_simulation, gc_gpu_vs_cg_no_cpu) {
//int main() {

    char *out_dir_no_gpu_no_precond = "tests_bin/circle_cg_no_gpu_no_precond";
    char *out_dir_no_gpu_precond  = "tests_bin/circle_cg_no_gpu_precond";

    char *out_dir_gpu_no_precond = "tests_bin/circle_cg_gpu_no_precond";
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
#endif

