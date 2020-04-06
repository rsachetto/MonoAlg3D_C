
////
//// Created by sachetto on 06/10/17.
////

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

void test_solver(bool preconditioner, char *method_name, char *init_name, char *end_name, int nt, int version) {

    FILE *A = NULL;
    FILE *B = NULL;
    FILE *X = NULL;

    if (version == 1) {
        A = fopen("src/tests/A1.txt", "r");
        B = fopen("src/tests/B1.txt", "r");
        X = fopen("src/tests/X1.txt", "r");
    } else if (version == 2) {
        A = fopen("src/tests/A2.txt", "r");
        B = fopen("src/tests/B2.txt", "r");
        X = fopen("src/tests/X2.txt", "r");
    }  else if (version == 3) {
        A = fopen("src/tests/A3.txt", "r");
        B = fopen("src/tests/B3.txt", "r");
        X = fopen("src/tests/X3.txt", "r");
    }

    assert(A);
    assert(B);
    assert(X);

    real_cpu error;

    struct grid *grid = new_grid();
    assert (grid);

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

    struct time_info ti = ZERO_TIME_INFO;

    CALL_INIT_LINEAR_SYSTEM(linear_system_solver_config, grid);
    ((linear_system_solver_fn*)linear_system_solver_config->main_function)(&ti, linear_system_solver_config, grid, grid->num_active_cells, grid->active_cells, &n_iter, &error);
    CALL_END_LINEAR_SYSTEM(linear_system_solver_config);

    clean_and_free_grid(grid);
    fclose(A);
    fclose(B);
}

int main() {

    test_solver(false, "gpu_conjugate_gradient", "init_gpu_conjugate_gradient", "end_gpu_conjugate_gradient",  6, 3);

}