#include <stdlib.h>
#include "alg/cell/cell.h"
#include "alg/grid/grid.h"
#include "3dparty/ini_parser/ini.h"
#include "config/domain_config.h"
#include "eikonal/eikonal_solver.h"
#include "3dparty/stb_ds.h"
#include "logger/logger.h"
#include "config/config_parser.h"


static int compare_coordinates(const void *a, const void *b) {

    struct cell_node *coord1 = *(struct  cell_node **) a;
    struct cell_node *coord2 = *(struct  cell_node **) b;

    if (coord1->center.y != coord2->center.y) {
        return coord1->center.y - coord2->center.y;
    }

    if (coord1->center.x != coord2->center.x) {
          return coord1->center.x - coord2->center.x;
    }


    return coord1->center.z - coord2->center.z;

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

    for(int i = 0 ; i < grid->num_active_cells; i++) {
        int new_pos_x = (grid->active_cells[i]->center.x-grid->active_cells[i]->discretization.x/2)/grid->active_cells[i]->discretization.x;
        int new_pos_y = (grid->active_cells[i]->center.y-grid->active_cells[i]->discretization.y/2)/grid->active_cells[i]->discretization.y;
        int new_pos_z = (grid->active_cells[i]->center.z-grid->active_cells[i]->discretization.z/2)/grid->active_cells[i]->discretization.z;

        grid->active_cells[i]->center.x = new_pos_x;
        grid->active_cells[i]->center.y = new_pos_y;
        grid->active_cells[i]->center.z = new_pos_z;
    }

    size_t itersPerBlock = 10, type = 1;

    struct eikonal_solver *solver = new_eikonal_solver(false);
    solver->width = grid->cube_side_length.x/grid->active_cells[0]->discretization.x;
    solver->height = grid->cube_side_length.y/grid->active_cells[0]->discretization.y;
    solver->depth = grid->cube_side_length.z/grid->active_cells[0]->discretization.z;

    solver->solver_type = type;
    solver->iters_per_block = itersPerBlock;

    size_t *initial_seed = calloc(sizeof(size_t), 3);
    initial_seed[0] = grid->active_cells[0]->center.x;
    initial_seed[1] = grid->active_cells[0]->center.y;
    initial_seed[2] = grid->active_cells[0]->center.z;

    arrput(solver->seeds, initial_seed);
    solver->active_cells = grid->active_cells;
    solver->num_active_cells = grid->num_active_cells;

    solve_eikonal(solver);
    write_alg(solver, "test.alg");

    free_eikonal_solver(solver);

}
