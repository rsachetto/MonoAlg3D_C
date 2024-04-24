
#include "../common_types/common_types.h"

struct heap_point {
    struct cell_node *grid_cell;
    real distance;
    real time;
    real repolarization_time;
};

struct point_distance_heap {
    struct heap_point * arr;
    int size;
    int capacity;
};

struct point_distance_heap * build_heap(struct heap_point* data, int size, int capacity);
struct heap_point heap_pop(struct point_distance_heap* h); 
void heap_push(struct point_distance_heap *h, struct heap_point data);
