#include "3dparty/ini_parser/ini.h"
#include "3dparty/stb_ds.h"
#include "alg/cell/cell.h"
#include "alg/grid/grid.h"
#include "config/config_parser.h"
#include "config/domain_config.h"
#include "eikonal/eikonal_solver.h"
#include "logger/logger.h"
#include "utils/file_utils.h"
#include <stdlib.h>

static int compare_coordinates(const void *a, const void *b) {

    struct cell_node *coord1 = *(struct  cell_node **) a;
    struct cell_node *coord2 = *(struct  cell_node **) b;

    if (coord1->center.y != coord2->center.y) {
        return (int) (coord1->center.y - coord2->center.y);
    }

    if (coord1->center.x != coord2->center.x) {
          return (int) (coord1->center.x - coord2->center.x);
    }


    return (int) (coord1->center.z - coord2->center.z);

}

int main(int argc, char *argv[]) {

    struct eikonal_options *eikonal_options = new_eikonal_options();
    parse_eikonal_options(argc, argv, eikonal_options);
    struct grid *grid = new_grid();

    if(ini_parse(eikonal_options->config_file, parse_eikonal_config_file, eikonal_options) < 0) {
        log_error_and_exit("Error: Can't load the config file %s\n", eikonal_options->config_file);
    }

    struct config *domain_config = eikonal_options->domain_config;

    // Configure the functions and set the mesh domain
    if(domain_config) {

        init_config_functions(domain_config, "./shared_libs/libdefault_domains.so", "domain");

        print_domain_config_values(domain_config);
        log_msg(LOG_LINE_SEPARATOR);

        bool success = ((set_spatial_domain_fn *)eikonal_options->domain_config->main_function)(domain_config, grid);

        if(!success) {
            log_error_and_exit("Error configuring the tissue domain!\n");
        }

        order_grid_cells(grid);

        qsort(grid->active_cells, grid->num_active_cells, sizeof(grid->active_cells[0]), compare_coordinates);

    }

    //Translating to 0
    for(int i = 0 ; i < grid->num_active_cells; i++) {
        double new_pos_x = (grid->active_cells[i]->center.x-grid->active_cells[i]->discretization.x/2.0)/grid->active_cells[i]->discretization.x;
        double new_pos_y = (grid->active_cells[i]->center.y-grid->active_cells[i]->discretization.y/2.0)/grid->active_cells[i]->discretization.y;
        double new_pos_z = (grid->active_cells[i]->center.z-grid->active_cells[i]->discretization.z/2.0)/grid->active_cells[i]->discretization.z;

        if(new_pos_x != floor(new_pos_x) || new_pos_y != floor(new_pos_y) || new_pos_z != floor(new_pos_z)) {
            log_error_and_exit("The current version only accepts integer coordinates\n");
        }

        grid->active_cells[i]->center.x = new_pos_x;
        grid->active_cells[i]->center.y = new_pos_y;
        grid->active_cells[i]->center.z = new_pos_z;
    }

    size_t itersPerBlock = 10, type = 1;

    struct eikonal_solver *solver = new_eikonal_solver(false);
    solver->width  = (int) (grid->cube_side_length.x/grid->active_cells[0]->discretization.x);
    solver->height = (int) (grid->cube_side_length.y/grid->active_cells[0]->discretization.y);
    solver->depth  = (int) (grid->cube_side_length.z/grid->active_cells[0]->discretization.z);

    solver->solver_type = type;
    solver->iters_per_block = itersPerBlock;

    size_t *initial_seed = calloc(sizeof(size_t), 3);
    initial_seed[0] = 128;
    initial_seed[1] = 48;
    initial_seed[2] = 63;

    arrput(solver->seeds, initial_seed);
    solver->active_cells = grid->active_cells;
    solver->num_active_cells = grid->num_active_cells;

    solve_eikonal(solver);

    struct config * save_result_config = eikonal_options->save_mesh_config;

    if(save_result_config) {

        init_config_functions(save_result_config, "./shared_libs/libdefault_save_mesh.so", "save_result");

        print_save_mesh_config_values(save_result_config);
        log_msg(LOG_LINE_SEPARATOR);

        char *out_dir_name = NULL;
        GET_PARAMETER_STRING_VALUE_OR_USE_DEFAULT(out_dir_name, save_result_config, "output_dir");
        if(out_dir_name != NULL) {
            create_dir(out_dir_name);
            struct time_info ti = ZERO_TIME_INFO;
            ((save_mesh_fn *)save_result_config->main_function)(&ti, save_result_config, grid, NULL, NULL);
        } else {
            log_warn("Not output dir provided. The result will not be saved!");
        }
    }

    free(initial_seed);
    free_eikonal_solver(solver);
    free_eikonal_options(eikonal_options);

}
