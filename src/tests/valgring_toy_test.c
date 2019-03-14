//
// Created by sachetto on 01/03/19.
//

#include "../alg/grid/grid.h"
#include "../config/domain_config.h"
#include "../config/save_mesh_config.h"
#include "../single_file_libraries/stb_ds.h"
#include "../string/sds.h"

int main() {

    struct grid *grid = new_grid();
    grid->adaptive = true;

    initialize_and_construct_grid(grid, 1, 1, 1);

    struct save_mesh_config *save_mesh_config = new_save_mesh_config();

    save_mesh_config->config_data.function_name = strdup("save_as_vtu");
    save_mesh_config->out_dir_name = strdup("./tests_bin");
    save_mesh_config->print_rate = 1;

    init_save_mesh_functions(save_mesh_config);
    shput(save_mesh_config->config_data.config, "file_prefix", strdup("test_valgrind"));
    shput(save_mesh_config->config_data.config, "compress", strdup("yes"));

    save_mesh_config->save_mesh(0, 0.0, save_mesh_config, grid);

    refine_grid_cell(grid, grid->first_cell);

    save_mesh_config->save_mesh(1, 1.0, save_mesh_config, grid);

    derefine_grid_cell(grid, grid->first_cell);

    save_mesh_config->save_mesh(2, 2.0, save_mesh_config, grid);

    free_save_mesh_config(save_mesh_config);

    clean_and_free_grid(grid);

}