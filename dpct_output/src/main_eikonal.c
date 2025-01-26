#include "3dparty/ini_parser/ini.h"
#include "3dparty/stb_ds.h"
#include "alg/cell/cell.h"
#include "alg/grid/grid.h"
#include "config/config_parser.h"
#include "config/domain_config.h"
#include "eikonal/eikonal_solver.h"
#include "config/stim_config.h"
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
    struct string_voidp_hash_entry *stimuli_configs = eikonal_options->stim_configs;

    struct eikonal_solver *solver = new_eikonal_solver(false);

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

    if(stimuli_configs) {

        size_t n = shlen(stimuli_configs);
        struct time_info time_info = {0.0, 0.0, 0.0, 0};

        if(n > 0) {
            STIM_CONFIG_HASH_FOR_INIT_FUNCTIONS(stimuli_configs);
        }

        for(int i = 0; i < n; i++) {

            struct string_voidp_hash_entry e = stimuli_configs[i];
            log_info("Stimulus name: %s\n", e.key);
            struct config *tmp = (struct config *)e.value;
            print_stim_config_values(tmp);
            log_msg(LOG_LINE_SEPARATOR);
        }

        set_spatial_stim(&time_info, stimuli_configs, grid, false);

        struct config *tmp = NULL;

        uint32_t n_active = grid->num_active_cells;
        struct cell_node **ac = grid->active_cells;

        for(size_t k = 0; k < n; k++) {

            tmp = (struct config *)stimuli_configs[k].value;

            for(uint32_t i = 0; i < n_active; i++) {
                real value =  ((real *)(tmp->persistent_data))[i];
                if(value != 0.0) {
                    size_t *initial_seed = (size_t *)calloc(sizeof(size_t), 3);
                    double new_pos_x = (ac[i]->center.x-ac[i]->discretization.x/2.0)/ac[i]->discretization.x;
                    double new_pos_y = (ac[i]->center.y-ac[i]->discretization.y/2.0)/ac[i]->discretization.y;
                    double new_pos_z = (ac[i]->center.z-ac[i]->discretization.z/2.0)/ac[i]->discretization.z;

                    initial_seed[0] = new_pos_x;
                    initial_seed[1] = new_pos_y;
                    initial_seed[2] = new_pos_z;
                    arrput(solver->seeds, initial_seed);
                }
            }
        }
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

    solver->width  = (int) (grid->cube_side_length.x/grid->active_cells[0]->discretization.x);
    solver->height = (int) (grid->cube_side_length.y/grid->active_cells[0]->discretization.y);
    solver->depth  = (int) (grid->cube_side_length.z/grid->active_cells[0]->discretization.z);

    solver->solver_type = type;
    solver->iters_per_block = itersPerBlock;

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

    free_eikonal_solver(solver);
    free_eikonal_options(eikonal_options);

}
