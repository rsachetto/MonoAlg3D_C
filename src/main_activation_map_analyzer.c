#include "alg/grid/grid.h"
#include "ini_parser/ini.h"
#include "monodomain/monodomain_solver.h"
#include "monodomain/ode_solver.h"
#include "string/sds.h"
#include "utils/file_utils.h"
#include <string.h>

#ifdef COMPILE_OPENGL
#include "config_helpers/config_helpers.h"

#endif

int main(int argc, char **argv) {

    struct grid *the_grid;

    the_grid = new_grid();

    struct config *domain_config = alloc_and_init_config_data();

    domain_config->main_function_name = strdup("initialize_from_activation_map_file");
    shput_dup_value(domain_config->config_data, "mesh_file", strdup(argv[1]));
    shput_dup_value(domain_config->config_data, "start_dx", "100");
    shput_dup_value(domain_config->config_data, "start_dy", "100");
    shput_dup_value(domain_config->config_data, "start_dz", "100");

    shput_dup_value(domain_config->config_data, "num_layers", "1");
    shput_dup_value(domain_config->config_data, "side_length", "40000");

    init_config_functions(domain_config, "shared_libs/libdefault_domains.so", "domain");

    ((set_spatial_domain_fn*)domain_config->main_function)(domain_config, the_grid);

    if(argc == 3) {
        struct config *save_mesh_config = alloc_and_init_config_data();
        save_mesh_config->main_function_name = strdup("save_as_text_or_binary");
        shput_dup_value(save_mesh_config->config_data, "output_dir", argv[2]);
        shput_dup_value(save_mesh_config->config_data, "file_prefix", "V");
        init_config_functions(save_mesh_config, "shared_libs/libdefault_save_mesh.so", "save_mesh");
        create_dir(argv[2]);
        ((save_mesh_fn *) save_mesh_config->main_function)(save_mesh_config, the_grid, 0, 0.0, 0.0, 0.0);
    }

    return EXIT_SUCCESS;
}
