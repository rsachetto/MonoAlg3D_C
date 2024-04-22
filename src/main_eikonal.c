#include <stdlib.h>
#include <math.h>
#include "3dparty/stb_ds.h"
#include "alg/grid/grid.h"
#include "config/domain_config.h"
#include "config/save_mesh_config.h"
#include "utils/heap.h"
#include "logger/logger.h"


static int in_array(struct heap_point *arr, struct cell_node *cell) {

    int n = arrlen(arr);
    for(int i = 0; i < n; i++) {
        if(arr[i].grid_cell->center.x == cell->center.x && arr[i].grid_cell->center.y == cell->center.y && 
            arr[i].grid_cell->center.z == cell->center.z) {
            return i;
        }
    }
    return -1;
}

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

void compute_t(struct heap_point x, void *neighbour, struct heap_point *known, real wave_speed, struct point_distance_heap *heap, struct point_hash_entry **distances, struct point_hash_entry **time) {
    

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
        
        if(!neighbour_cell->active || in_array(known, neighbour_cell) != -1) {
            return;
        }

        real tentative_distance = x.distance + wave_speed * euclidean_distance(grid_cell, neighbour_grid_cell);
        real neighbour_distance = hmget(*distances, neighbour_cell->center);

        if(tentative_distance < neighbour_distance) {
            struct heap_point xi;
            xi.grid_cell = neighbour_cell;
            
            hmput(*distances, neighbour_cell->center, tentative_distance);
            hmput(*time, neighbour_cell->center, tentative_distance/wave_speed);
            
            int index = in_heap(heap, xi, heap->size);


        printf("Computing for %lf, %lf, %lf. Tentative %lf, distance %lf, index %d\n", 
                neighbour_cell->center.x, neighbour_cell->center.y, neighbour_cell->center.z,
                tentative_distance, neighbour_distance, index);

            if(index == -1) {
                xi.distance = tentative_distance;                
                heap_push(heap, xi);
            } else {
                heap->arr[index].distance = tentative_distance;
            }
        }
    }
}


int main(int argc, char **argv) {

    struct grid *grid = new_grid();   

    char *discretization = "200.0";
    char *side_length = "20000.0";
    char *num_layers = "1";
    real wave_speed = 1.0;       

    struct config *domain_config = domain_config = alloc_and_init_config_data();

    domain_config->main_function_name = strdup("initialize_grid_with_square_mesh");
    shput_dup_value(domain_config->config_data, "name", "Test custom mesh");

    shput(domain_config->config_data, "start_dx", strdup(discretization));
    shput(domain_config->config_data, "start_dy", strdup(discretization));
    shput(domain_config->config_data, "start_dz", strdup(discretization));
    shput(domain_config->config_data, "side_length", strdup(side_length));
    shput(domain_config->config_data, "num_layers", strdup(num_layers));

    init_config_functions(domain_config, "./shared_libs/libdefault_domains.so", "domain");

    int success = ((set_spatial_domain_fn*)domain_config->main_function)(domain_config, grid);

    if(!success) {
        clean_and_free_grid(grid);
        free_config_data(domain_config);
        log_error_and_exit("Error loading the domain");
    }

    order_grid_cells(grid);

    struct heap_point *data = (struct heap_point *) malloc(grid->num_active_cells * sizeof(struct heap_point));
    struct point_distance_heap *trial = build_heap(data, 0, grid->num_active_cells);

    struct point_hash_entry *distances = NULL;
    hmdefault(distances, INFINITY);

    struct point_hash_entry *times = NULL;
    hmdefault(times, INFINITY);

    struct heap_point *known = NULL;

    //TODO: this has to come from the stimulus
    struct heap_point first_node;

    for(int i = 0; i < grid->num_active_cells; i++) {
        struct cell_node *cell = grid->active_cells[i];
        if(cell->center.x == 100.0 && cell->center.y == 100.0 && cell->center.z == 100.0) {
            first_node.grid_cell = cell;
            first_node.distance = 0.0;
            arrpush(known, first_node);
            hmput(distances, cell->center, 0.0);
            hmput(times, cell->center, 0.0);
        }
        
    }
    
    printf("First cell: %f, %f, %f\n", first_node.grid_cell->center.x, first_node.grid_cell->center.y, first_node.grid_cell->center.z);

    for(int i = 0; i < arrlen(known); i++) {
        struct heap_point x = known[i];
        for (int j = 0; j <= LEFT; j++) {
            compute_t(x, x.grid_cell->neighbours[j], known, wave_speed, trial, &distances, &times); 
        }        
    }

    while(trial->size > 0) {
        struct heap_point x = heap_pop(trial);
        arrpush(known, x);
        for (int i = 0; i <= LEFT; i++) {
            compute_t(x, x.grid_cell->neighbours[i], known, wave_speed, trial, &distances, &times); 
        }
    }

    // Print times hash map
    for (int i = 0; i < arrlen(known); i++) {
        struct cell_node *cell = known[i].grid_cell;
        real time = hmget(times, cell->center);
        cell->v = time;
        printf("Time for cell (%f, %f, %f): %f\n", cell->center.x, cell->center.y, cell->center.z, time);
    }

    struct config *save_mesh_config = alloc_and_init_config_data();

    save_mesh_config->main_function_name = strdup("save_as_text_or_binary");
    shput_dup_value(save_mesh_config->config_data, "output_dir", "./test_eikonal");
    shput_dup_value(save_mesh_config->config_data, "print_rate", "1");
    init_config_functions(save_mesh_config, "shared_libs/libdefault_save_mesh.so", "save_result");

    shput(save_mesh_config->config_data, "file_prefix", strdup("AT"));

    struct time_info ti = ZERO_TIME_INFO;

    ((save_mesh_fn *)save_mesh_config->main_function)(&ti, save_mesh_config, grid, NULL, NULL);

    free_config_data(save_mesh_config);

}
