//
// Created by sachetto on 01/03/19.
//

#include "../alg/grid/grid.h"
#include "../config/domain_config.h"
#include "../config/save_mesh_config.h"
#include "../single_file_libraries/stb_ds.h"
#include "../string/sds.h"

int main() {

//    FILE *A = NULL;
//
//    A = fopen("src/tests/A_CSR_Test.txt", "r");
//
//
//    struct grid *grid = new_grid();
//
//    construct_grid_from_file(grid, A, NULL);
//    print_grid_matrix(grid, stdout);
//
//    real *MA = NULL;
//    int *IA = NULL;
//    int *JA = NULL;
//
//    grid_to_csr(grid, &MA, &IA, &JA);
//
//    for(int i = 0; i < arrlen(MA);i++) {
//        printf("%lf, ", MA[i]);
//    }
//
//    printf("\n");
//
//    for(int i = 0; i < arrlen(IA);i++) {
//        printf("%d, ", IA[i]);
//    }
//    printf("\n");
//
//    for(int i = 0; i < arrlen(JA);i++) {
//        printf("%d, ", JA[i]);
//    }
//
//    printf("\n");

    struct grid *grid = new_grid();
    grid->adaptive = true;

    initialize_and_construct_grid(grid, POINT3D(1, 1, 1));

    struct config *save_mesh_config = alloc_and_init_config_data();

    save_mesh_config->main_function_name = strdup("save_as_vtu");

    shput_dup_value(save_mesh_config->config_data, "output_dir", "./tests_bin");
    shput_dup_value(save_mesh_config->config_data, "print_rate", "1");

    init_config_functions(save_mesh_config, "./shared_libs/libdefault_save_mesh.so", "save_result");
    shput_dup_value(save_mesh_config->config_data, "file_prefix", "test_valgrind");
    shput_dup_value(save_mesh_config->config_data, "compress", "yes");

    ((save_mesh_fn*)save_mesh_config->main_function)(save_mesh_config, grid, 0, 0.0, 0.0, 0.0,'v');

    refine_grid_cell(grid, grid->first_cell);

    ((save_mesh_fn*)save_mesh_config->main_function)(save_mesh_config, grid, 1, 1.0, 1.0, 0.0,'v');

    derefine_grid_cell(grid, grid->first_cell);

    ((save_mesh_fn*)save_mesh_config->main_function)(save_mesh_config, grid, 2, 2.0, 2.0, 0.0,'v');

    free_config_data(save_mesh_config);

    clean_and_free_grid(grid);

}
