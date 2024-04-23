#include <stdlib.h>
#include <math.h>
#include "alg/grid/grid.h"
#include "config/domain_config.h"
#include "config/save_mesh_config.h"
#include "config/stim_config.h"
#include "utils/heap.h"
#include "logger/logger.h"
#include "3dparty/ini_parser/ini.h"
#include "utils/file_utils.h"

struct heap_point_hash_entry {
    struct point_3d key;
    struct heap_point value;
};

static struct point_3d compute_wave_speed(struct point_3d conductivity_tensor) {
    //Compute wave speed based on conductivity tensor (e.g., using eigenvalues)
    //Example: For isotropic conductivity, return sqrt(1 / conductivity_tensor)
    //For anisotropic conductivity, compute wave speed based on eigenvalues
    return POINT3D(sqrt(conductivity_tensor.x), sqrt(conductivity_tensor.y), sqrt(conductivity_tensor.z));
}

//TODO: add an hash map in the heap to make the search faster
static int in_heap(struct point_distance_heap *heap, struct heap_point x, int n) {
    for(int i = 0; i < n; i++) {
        struct cell_node *cell = x.grid_cell;
        if(heap->arr[i].grid_cell->center.x == cell->center.x && heap->arr[i].grid_cell->center.y == cell->center.y && 
            heap->arr[i].grid_cell->center.z == cell->center.z) {
            return i;
        }
    }
    return -1;
}

static real euclidean_distance(struct cell_node *cell1, struct cell_node *cell2) {

    real dx = cell1->center.x - cell2->center.x;
    real dy = cell1->center.y - cell2->center.y;
    real dz = cell1->center.z - cell2->center.z;

    return sqrt(dx * dx + dy * dy + dz * dz);
}

void compute_t(struct heap_point x, void *neighbour, struct heap_point_hash_entry *known, struct point_3d condutivity_tensor, struct point_distance_heap *heap) {    

    void *neighbour_grid_cell = neighbour;
    struct cell_node *grid_cell = x.grid_cell;
    bool has_found;

    struct transition_node *white_neighbor_cell;

    uint16_t neighbour_grid_cell_level = ((struct basic_cell_data *)(neighbour_grid_cell))->level;
    enum cell_type neighbour_grid_cell_type = ((struct basic_cell_data *)(neighbour_grid_cell))->type;

    if(neighbour_grid_cell_level > grid_cell->cell_data.level) {
        if(neighbour_grid_cell_type == TRANSITION_NODE) {
            has_found = false;
            while(!has_found) {
                if(neighbour_grid_cell_type == TRANSITION_NODE) {
                    white_neighbor_cell = (struct transition_node *)neighbour_grid_cell;
                    if(white_neighbor_cell->single_connector == NULL) {
                        has_found = true;
                    } else {
                        neighbour_grid_cell = white_neighbor_cell->quadruple_connector1;
                        neighbour_grid_cell_type = ((struct basic_cell_data *)(neighbour_grid_cell))->type;
                    }
                } else {
                    break;
                }
            }
        }
    } else {
        if(neighbour_grid_cell_level <= grid_cell->cell_data.level && (neighbour_grid_cell_type == TRANSITION_NODE)) {
            has_found = false;
            while(!has_found) {
                if(neighbour_grid_cell_type == TRANSITION_NODE) {
                    white_neighbor_cell = (struct transition_node *)(neighbour_grid_cell);
                    if(white_neighbor_cell->single_connector == 0) {
                        has_found = true;
                    } else {
                        neighbour_grid_cell = white_neighbor_cell->single_connector;
                        neighbour_grid_cell_type = ((struct basic_cell_data *)(neighbour_grid_cell))->type;
                    }
                } else {
                    break;
                }
            }
        }
    }

    if(neighbour_grid_cell_type == CELL_NODE) {

        struct cell_node *neighbour_cell = (struct cell_node *)(neighbour_grid_cell);                
        
        struct heap_point tmp = hmget(known, neighbour_cell->center);

        if(!neighbour_cell->active || tmp.grid_cell != NULL) {
            return;
        }

        struct heap_point xi;
        xi.grid_cell = neighbour_cell;

        int index = in_heap(heap, xi, heap->size);        

        real neighbour_distance = INFINITY;

        real wave_speed = compute_wave_speed(condutivity_tensor).x;

        //TODO: we have to know the neighbour's direction to calculate the wave speed
        real tentative_distance = x.distance + wave_speed * euclidean_distance(grid_cell, neighbour_grid_cell);
        real time = tentative_distance/wave_speed;

        if(index != -1) {            
            neighbour_distance = heap->arr[index].distance;
        } else {
            xi.distance = tentative_distance;   
            xi.grid_cell->v = time;             
            heap_push(heap, xi);
            return;
        }

        if(tentative_distance < neighbour_distance) {       
            heap->arr[index].distance = tentative_distance;
            heap->arr[index].grid_cell->v = time;
        }
    }
}


int main(int argc, char **argv) {

    struct eikonal_options *eikonal_options = new_eikonal_options();
    parse_eikonal_options(argc, argv, eikonal_options);
    struct grid *grid = new_grid();  
    
    if (ini_parse(eikonal_options->config_file, parse_eikonal_config_file, eikonal_options) < 0) {
        fprintf(stderr, "Error: Can't load the config file %s\n", eikonal_options->config_file);
        exit(EXIT_FAILURE);
    }

    struct config *domain_config = eikonal_options->domain_config;
    struct config *save_mesh_config = eikonal_options->save_mesh_config;
    struct string_voidp_hash_entry *stimuli_configs = eikonal_options->stim_configs;

    bool save_to_file = (save_mesh_config != NULL);
    char *out_dir_name = NULL;
    
    if(save_to_file) {
        init_config_functions(save_mesh_config, "./shared_libs/libdefault_save_mesh.so", "save_result");
        GET_PARAMETER_STRING_VALUE_OR_USE_DEFAULT(out_dir_name, save_mesh_config, "output_dir");

        if(out_dir_name == NULL) {
            log_error("No output directory provided to save the results! Exiting!\n");         
            exit(EXIT_FAILURE);
        }
        
        create_dir(out_dir_name);

    } else {
        log_error("No configuration provided to save the results! Exiting!\n");
        free(out_dir_name);
        exit(EXIT_FAILURE);
    }


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
    }        

    struct point_3d condutivity_tensor = SAME_POINT3D(0.00005336);       
    
    struct heap_point *data = (struct heap_point *) malloc(grid->num_active_cells * sizeof(struct heap_point));
    struct point_distance_heap *trial = build_heap(data, 0, grid->num_active_cells);

    struct heap_point_hash_entry *known = NULL;
    struct heap_point default_entry = {NULL, -1};
    hmdefault(known, default_entry);


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

            real stim_start = 0.0;
            real stim_dur = 0.0;
            real stim_period = 0.0;

            tmp = (struct config *)stimuli_configs[k].value;
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_start, tmp, "start");
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_dur, tmp, "duration");
            GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, stim_period, tmp, "period");

            for(uint32_t i = 0; i < n_active; i++) {
                real value =  ((real *)(tmp->persistent_data))[i];
                if(value != 0.0) {
                    struct heap_point node;
                    node.grid_cell = ac[i];
                    if(stim_start == 0.0) {
                        node.distance = 0.0;
                    } else {
                        node.distance = compute_wave_speed(condutivity_tensor).x / stim_start;
                    }
                    hmput(known, ac[i]->center, node);
                }
            }
            

            // if(stim_period > 0.0) {
            //     if(time >= stim_start + stim_period) {
            //         stim_start = stim_start + stim_period;
            //         sds stim_start_char = sdscatprintf(sdsempty(), "%lf", stim_start);
            //         shput_dup_value(tmp->config_data, "start", stim_start_char);
            //         sdsfree(stim_start_char);
            //     }
            // }

            // time = cur_time;
        }

    }

    for(int i = 0; i < hmlen(known); i++) {
        struct heap_point x = known[i].value;
        for (int j = 0; j <= LEFT; j++) {
            compute_t(x, x.grid_cell->neighbours[j], known, condutivity_tensor, trial); 
        }        
    }

    while(trial->size > 0) {
        struct heap_point x = heap_pop(trial);
        hmput(known, x.grid_cell->center, x);
        for (int i = 0; i <= LEFT; i++) {
            compute_t(x, x.grid_cell->neighbours[i], known, condutivity_tensor, trial); 
        }
    }

    struct time_info ti = ZERO_TIME_INFO;
    CALL_INIT_SAVE_MESH(save_mesh_config);
    ((save_mesh_fn *)save_mesh_config->main_function)(&ti, save_mesh_config, grid, NULL, NULL);
    CALL_END_SAVE_MESH(save_mesh_config, grid);
    free_config_data(save_mesh_config);

}
