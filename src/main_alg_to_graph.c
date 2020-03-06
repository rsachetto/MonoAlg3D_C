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

int **graph = NULL;
int **cycles = NULL;
int num_edges = 0;

// Function to mark the vertex with
// different colors for different cycles
void dfs_cycle(int u, int p, int color[], int mark[], int par[], int *cyclenumber) {

    // already (completely) visited vertex.
    if(color[u] == 2) {
        return;
    }

    // seen vertex, but was not completely visited -> cycle detected.
    // backtrack based on parents to find the complete cycle.
    if(color[u] == 1) {

        *cyclenumber = *cyclenumber + 1;
        int cur = p;
        mark[cur] = *cyclenumber-1;

        // backtrack the vertex which are
        // in the current cycle thats found
        while(cur != u) {
            cur = par[cur];
            mark[cur] = *cyclenumber-1;
        }
        return;
    }
    par[u] = p;

    // partially visited.
    color[u] = 1;

    // simple dfs on graph
    for(int i = 0; i < arrlen(graph[u]); i++) {
        int v = graph[u][i];
        // if it has not been visited previously
        if(v == par[u]) {
            continue;
        }
        dfs_cycle(v, u, color, mark, par, cyclenumber);
    }

    // completely visited.
    color[u] = 2;
}

// add the edges to the graph
void addEdge(int u, int v) {
    arrpush(graph[u], v);
    num_edges++;
    arrput(cycles, NULL);
}

// Function to print the cycles
void printCycles(int edges, int mark[], int *cyclenumber) {

    // push the edges that into the
    // cycle adjacency list
    for(int i = 0; i < edges; i++) {
        if(mark[i] != 0)
            arrpush(cycles[mark[i]], i);
    }

    // print all the vertex with same cycle
    for(int i = 0; i < *cyclenumber; i++) {
        // Print the i-th cycle
        printf("Cycle Number %d: ", i + 1);
        for(int j = 0; j < arrlen(cycles[i]); j++) {
            printf("%d ", cycles[i][j]);
        }
        printf("\n");
    }
}

void print_alg_graph() {

    int n = arrlen(graph);

    for(int i = 0; i < n; i++) {
        int m = arrlen(graph[i]);
        printf("%d ",i);
        for(int j = 0; j < m; j++) {
            printf("%d ", graph[i][j]);
        }

        printf("\n");
    }

}



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

    ((set_spatial_domain_fn *)domain_config->main_function)(NULL, domain_config, the_grid);

    order_grid_cells(the_grid);

    if(argc == 3) {
        struct config *save_mesh_config = alloc_and_init_config_data();
        save_mesh_config->main_function_name = strdup("save_as_adjacency_list");
        shput_dup_value(save_mesh_config->config_data, "output_dir", argv[2]);
        shput_dup_value(save_mesh_config->config_data, "file_prefix", "V");
        init_config_functions(save_mesh_config, "shared_libs/libdefault_save_mesh.so", "save_mesh");
        create_dir(argv[2]);

        struct time_info t;
        t.iteration = 0;
        t.dt = 0.0;
        t.current_t = 0;

        ((save_mesh_fn *)save_mesh_config->main_function)(&t, save_mesh_config, the_grid);
    }


//    for(int i = 0; i < 13; i++) {
//        arrput(graph, NULL);
//    }
//
//    // add edges
//    addEdge(0, 1);
//    addEdge(1, 0);
//
//    addEdge(1, 2);
//    addEdge(2, 1);
//
//    addEdge(2, 3);
//    addEdge(3, 2);
//
//    addEdge(3, 5);
//    addEdge(5, 3);
//
//    addEdge(3, 6);
//    addEdge(6, 3);
//
//    addEdge(4, 5);
//    addEdge(5, 4);
//
//
//    addEdge(2, 4);
//    addEdge(4, 2);
//
//    addEdge(6, 7);
//    addEdge(7, 6);
//
//    addEdge(5, 9);
//    addEdge(9, 5);
//
//    addEdge(4, 8);
//    addEdge(8, 4);
//
//    addEdge(9, 10);
//    addEdge(10, 9);
//
//    addEdge(10, 11);
//    addEdge(11, 10);
//
//    addEdge(10, 12);
//    addEdge(12, 10);
//
//    addEdge(11, 12);
//    addEdge(12, 11);
//
//    struct cell_node *neighbour;
//
//    FOR_EACH_CELL(the_grid) {
//
//        if(cell->active) {
//            arrput(graph, NULL);
//
//            neighbour = get_cell_neighbour(cell, cell->north);
//            if(neighbour) {
//                addEdge(cell->grid_position, neighbour->grid_position);
//            }
//            neighbour = get_cell_neighbour(cell, cell->south);
//            if(neighbour) {
//                addEdge(cell->grid_position, neighbour->grid_position);
//            }
//
//            neighbour = get_cell_neighbour(cell, cell->west);
//            if(neighbour) {
//                addEdge(cell->grid_position, neighbour->grid_position);
//            }
//
//            neighbour = get_cell_neighbour(cell, cell->east);
//            if(neighbour) {
//                addEdge(cell->grid_position, neighbour->grid_position);
//            }
//
//            neighbour = get_cell_neighbour(cell, cell->front);
//            if(neighbour) {
//                addEdge(cell->grid_position, neighbour->grid_position);
//            }
//
//            neighbour = get_cell_neighbour(cell, cell->back);
//            if(neighbour) {
//                addEdge(cell->grid_position, neighbour->grid_position);
//            }
//        }
//    }
//
////    print_alg_graph();
////    exit(0);
//
//    int *color = calloc(num_edges, sizeof(int));
//    int *par = calloc(num_edges, sizeof(int));
//
//    // mark with unique numbers
//    int *mark = calloc(num_edges, sizeof(int));
//
//    // store the numbers of cycle
//    int cyclenumber = 0;
//
//    // call DFS to mark the cycles
//    dfs_cycle(0, 0, color, mark, par, &cyclenumber);
//
//    // function to print the cycles
//    printCycles(num_edges, mark, &cyclenumber);

    return EXIT_SUCCESS;
}
