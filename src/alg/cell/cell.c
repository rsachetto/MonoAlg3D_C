//
// Created by sachetto on 29/09/17.
//

#include "cell.h"
#include <math.h>

#include "../../3dparty/stb_ds.h"

struct cell_node *new_cell_node() {
    struct cell_node *result = (struct cell_node *)malloc(sizeof(struct cell_node));
    init_cell_node(result);
    return result;
}

void init_cell_node(struct cell_node *cell_node) {

    cell_node->cell_data.type = CELL_NODE;
    cell_node->cell_data.level = 1;

    cell_node->center.x = 0.0;
    cell_node->center.y = 0.0;
    cell_node->center.z = 0.0;

    cell_node->active = true;

    cell_node->bunch_number = 0;

    cell_node->neighbours = (void **)calloc(sizeof(void *), NUM_NEIGHBOURS);

    cell_node->previous = NULL;
    cell_node->next = NULL;

    cell_node->grid_position = 0;
    cell_node->sv_position = 0;
    cell_node->hilbert_shape_number = 0;
    cell_node->discretization.x = 1.0;
    cell_node->discretization.y = 1.0;
    cell_node->discretization.z = 1.0;

    cell_node->v = 0;

    cell_node->front_flux = 0.0;
    cell_node->back_flux = 0.0;
    cell_node->top_flux = 0.0;
    cell_node->down_flux = 0.0;
    cell_node->right_flux = 0.0;
    cell_node->left_flux = 0.0;

    cell_node->b = 0.0;

    cell_node->can_change = true;
    cell_node->visited = false;
    cell_node->visible = TOP_IS_VISIBLE | RIGHT_IS_VISIBLE | DOWN_IS_VISIBLE | LEFT_IS_VISIBLE | BACK_IS_VISIBLE | FRONT_IS_VISIBLE;

    cell_node->elements = NULL;

    cell_node->linear_system_solver_extra_info = NULL;
    cell_node->mesh_extra_info = NULL;

    cell_node->original_position_in_file = -1;

#ifdef ENABLE_DDM
    cell_node->kappa.x = 0.0;
    cell_node->kappa.y = 0.0;
    cell_node->kappa.z = 0.0;
#endif

#if defined(_OPENMP)
    omp_init_lock(&(cell_node->updating));
#endif
}

void free_cell_node(struct cell_node *cell_node) {

    arrfree(cell_node->elements);
    free(cell_node->linear_system_solver_extra_info);
    free(cell_node->mesh_extra_info);
    free(cell_node->neighbours);

#if defined(_OPENMP)
    omp_destroy_lock(&(cell_node->updating));
#endif
    free(cell_node);
}

inline void lock_cell_node(struct cell_node *cell_node) {
#if defined(_OPENMP)
    omp_set_lock(&(cell_node->updating));
#endif
}

inline void unlock_cell_node(struct cell_node *cell_node) {
#if defined(_OPENMP)
    omp_unset_lock(&(cell_node->updating));
#endif
}

struct transition_node *new_transition_node() {
    struct transition_node *result = (struct transition_node *)malloc(sizeof(struct transition_node));
    init_transition_node(result);
    return result;
}

void init_transition_node(struct transition_node *transition_node) {

    transition_node->cell_data.type = TRANSITION_NODE;
    transition_node->cell_data.level = 1;

    transition_node->single_connector = NULL;
    transition_node->quadruple_connector1 = NULL;
    transition_node->quadruple_connector2 = NULL;
    transition_node->quadruple_connector3 = NULL;
    transition_node->quadruple_connector4 = NULL;
    transition_node->direction = NUM_DIRECTIONS;
}

enum transition_direction get_inverse_direction(enum transition_direction direction) {

    switch(direction) {
    case FRONT:
        return BACK;
    case BACK:
        return FRONT;
    case TOP:
        return DOWN;
    case DOWN:
        return TOP;
    case RIGHT:
        return LEFT;
    case LEFT:
        return RIGHT;
    default:
        fprintf(stderr, "get_inverse_direction(). Invalid cell direction %d! Exiting...\n", direction);
        exit(10);
    }
}

void set_transition_node_data(struct transition_node *the_transition_node, uint16_t level, enum transition_direction direction, void *single_connector,
                              void *quadruple_connector1, void *quadruple_connector2, void *quadruple_connector3, void *quadruple_connector4) {

    the_transition_node->cell_data.level = level;

    the_transition_node->direction = direction;

    the_transition_node->single_connector = single_connector;

    the_transition_node->quadruple_connector1 = quadruple_connector1;
    the_transition_node->quadruple_connector2 = quadruple_connector2;
    the_transition_node->quadruple_connector3 = quadruple_connector3;
    the_transition_node->quadruple_connector4 = quadruple_connector4;
}

void set_cell_node_data(struct cell_node *the_cell, struct point_3d discretization, uint64_t bunch_number, void **neighbours, void *previous, void *next,
                        uint32_t grid_position, uint8_t hilbert_shape_number, struct point_3d center) {

    the_cell->discretization = discretization;
    the_cell->bunch_number = bunch_number;

    if(neighbours) {
        memcpy(the_cell->neighbours, neighbours, sizeof(void *) * NUM_NEIGHBOURS);
    } else {
        the_cell->neighbours = NULL;
    }

    the_cell->previous = (struct cell_node *)previous;
    the_cell->next = (struct cell_node *)next;
    the_cell->grid_position = grid_position;
    the_cell->hilbert_shape_number = hilbert_shape_number;
    the_cell->center = center;
}

void set_cell_flux(struct cell_node *the_cell, enum transition_direction direction) {

    void *neighbour_grid_cell;
    struct transition_node *white_neighbor_cell;
    struct cell_node *black_neighbor_cell;

    if(VALID_SIMPLE_DIRECTION(direction)) {
        neighbour_grid_cell = the_cell->neighbours[direction];
    } else {
        fprintf(stderr, "set_cell_flux(). Invalid cell direction %c! Exiting...\n", direction);
        exit(10);
    }

    real_cpu least_distance_x = the_cell->discretization.x / 2.0;
    real_cpu least_distance_y = the_cell->discretization.y / 2.0;
    real_cpu least_distance_z = the_cell->discretization.z / 2.0;

    bool has_found;

    /* When neighbour_grid_cell is a transition node, looks for the next neighbor
     * cell which is a cell node. */
    uint16_t neighbour_grid_cell_level = ((struct basic_cell_data *)(neighbour_grid_cell))->level;
    uint8_t neighbour_grid_cell_type = ((struct basic_cell_data *)(neighbour_grid_cell))->type;

    if(neighbour_grid_cell_level > the_cell->cell_data.level) {
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
        if(neighbour_grid_cell_level <= the_cell->cell_data.level && (neighbour_grid_cell_type == TRANSITION_NODE)) {
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

        black_neighbor_cell = (struct cell_node *)(neighbour_grid_cell);

        if(black_neighbor_cell->discretization.x / 2.0 < least_distance_x)
            least_distance_x = black_neighbor_cell->discretization.x / 2.0;

        if(black_neighbor_cell->discretization.y / 2.0 < least_distance_y)
            least_distance_y = black_neighbor_cell->discretization.y / 2.0;

        if(black_neighbor_cell->discretization.z / 2.0 < least_distance_z)
            least_distance_z = black_neighbor_cell->discretization.z / 2.0;

        real_cpu local_flux_x = (the_cell->v - black_neighbor_cell->v) * (2.0 * least_distance_x);
        real_cpu local_flux_y = (the_cell->v - black_neighbor_cell->v) * (2.0 * least_distance_y);
        real_cpu local_flux_z = (the_cell->v - black_neighbor_cell->v) * (2.0 * least_distance_z);

        lock_cell_node(the_cell);

        switch(direction) {
        case BACK:
            if(local_flux_y > the_cell->back_flux)
                the_cell->back_flux += local_flux_y;
            break;

        case FRONT:
            if(local_flux_y > the_cell->front_flux)
                the_cell->front_flux += local_flux_y;
            break;

        case TOP:
            if(local_flux_x > the_cell->top_flux)
                the_cell->top_flux += local_flux_x;
            break;

        case DOWN:
            if(local_flux_x > the_cell->down_flux)
                the_cell->down_flux += local_flux_x;
            break;

        case RIGHT:
            if(local_flux_z > the_cell->right_flux)
                the_cell->right_flux += local_flux_z;
            break;

        case LEFT:
            if(local_flux_z > the_cell->left_flux)
                the_cell->left_flux += local_flux_z;
            break;

        default:
            break;
        }

        unlock_cell_node(the_cell);

        lock_cell_node(black_neighbor_cell);

        switch(direction) {
        case BACK:
            if(local_flux_y > black_neighbor_cell->front_flux)
                black_neighbor_cell->front_flux += local_flux_y;
            break;

        case FRONT:
            if(local_flux_y > black_neighbor_cell->back_flux)
                black_neighbor_cell->back_flux += local_flux_y;
            break;

        case TOP:
            if(local_flux_x > black_neighbor_cell->down_flux)
                black_neighbor_cell->down_flux += local_flux_x;
            break;

        case DOWN:
            if(local_flux_x > black_neighbor_cell->top_flux)
                black_neighbor_cell->top_flux += local_flux_x;
            break;

        case RIGHT:
            if(local_flux_z > black_neighbor_cell->left_flux)
                black_neighbor_cell->left_flux += local_flux_z;
            break;

        case LEFT:
            if(local_flux_z > black_neighbor_cell->right_flux)
                black_neighbor_cell->right_flux += local_flux_z;
            break;

        default:
            break;
        }

        unlock_cell_node(black_neighbor_cell);
    }
}

real_cpu get_cell_maximum_flux(struct cell_node *the_cell) {

    real_cpu maximumFlux = fabs(the_cell->top_flux);
    if(fabs(the_cell->front_flux) > maximumFlux)
        maximumFlux = fabs(the_cell->front_flux);

    if(fabs(the_cell->back_flux) > maximumFlux)
        maximumFlux = fabs(the_cell->back_flux);

    if(fabs(the_cell->down_flux) > maximumFlux)
        maximumFlux = fabs(the_cell->down_flux);

    if(fabs(the_cell->right_flux) > maximumFlux)
        maximumFlux = fabs(the_cell->right_flux);

    if(fabs(the_cell->left_flux) > maximumFlux)
        maximumFlux = fabs(the_cell->left_flux);

    return maximumFlux;
}

struct cell_node *get_cell_neighbour(struct cell_node *grid_cell, void *neighbour_grid_cell) {

    struct cell_node *black_neighbor_cell;

    uint16_t neighbour_grid_cell_level = ((struct basic_cell_data *)(neighbour_grid_cell))->level;
    enum cell_type neighbour_grid_cell_type = ((struct basic_cell_data *)(neighbour_grid_cell))->type;

    struct transition_node *white_neighbor_cell;

    if(neighbour_grid_cell_type == TRANSITION_NODE) {

        while(true) {

            if(neighbour_grid_cell_type == TRANSITION_NODE) {
                white_neighbor_cell = (struct transition_node *)neighbour_grid_cell;
                if(white_neighbor_cell->single_connector == NULL) {
                    return NULL;
                } else {

                    if(neighbour_grid_cell_level > grid_cell->cell_data.level) {
                        neighbour_grid_cell = white_neighbor_cell->quadruple_connector1;
                    } else {
                        neighbour_grid_cell = white_neighbor_cell->single_connector;
                    }

                    neighbour_grid_cell_type = ((struct basic_cell_data *)(neighbour_grid_cell))->type;
                }
            } else {

                black_neighbor_cell = (struct cell_node *)(neighbour_grid_cell);
                if(black_neighbor_cell->active) {
                    return black_neighbor_cell;
                } else {
                    return NULL;
                }
            }
        }
    } else {
        black_neighbor_cell = (struct cell_node *)(neighbour_grid_cell);
        if(black_neighbor_cell->active) {
            return black_neighbor_cell;
        } else {
            return NULL;
        }
    }
}

bool cell_has_neighbour(struct cell_node *grid_cell, void *neighbour_grid_cell) {

    struct cell_node *black_neighbor_cell = get_cell_neighbour(grid_cell, neighbour_grid_cell);

    if(black_neighbor_cell) {
        return true;
    } else {
        return false;
    }
}

struct cell_node *get_cell_neighbour_with_same_refinement_level(struct cell_node *grid_cell, enum transition_direction direction) {

    struct basic_cell_data *cell_data = NULL;

    if(VALID_SIMPLE_DIRECTION(direction)) {
        cell_data = (struct basic_cell_data *)grid_cell->neighbours[direction];
        if(cell_data->type == TRANSITION_NODE) {
            return NULL;
        } else {
            return (struct cell_node *)cell_data;
        }
    } else if(VALID_DIAGONAL_POINT_DIRECTION(direction)) {

        if(STARTS_WITH_FRONT(direction)) {
            cell_data = (struct basic_cell_data *)grid_cell->neighbours[FRONT];
        }

        if(STARTS_WITH_BACK(direction)) {
            cell_data = (struct basic_cell_data *)grid_cell->neighbours[BACK];
        }

        if(STARTS_WITH_TOP(direction)) {
            cell_data = (struct basic_cell_data *)grid_cell->neighbours[TOP];
        }

        if(STARTS_WITH_DOWN(direction)) {
            cell_data = (struct basic_cell_data *)grid_cell->neighbours[DOWN];
        }

        if(cell_data->type == TRANSITION_NODE) {
            return NULL;
        }

        struct cell_node *node = (struct cell_node *)cell_data;

        switch(direction) {
        case FRONT_TOP:
        case BACK_TOP: {
            cell_data = (struct basic_cell_data *)node->neighbours[TOP];
            break;
        }
        case FRONT_DOWN:
        case BACK_DOWN: {
            cell_data = (struct basic_cell_data *)node->neighbours[DOWN];
            break;
        }

        case FRONT_RIGHT:
        case BACK_RIGHT:
        case TOP_RIGHT:
        case DOWN_RIGHT: {
            cell_data = (struct basic_cell_data *)node->neighbours[RIGHT];
            break;
        }

        case FRONT_LEFT:
        case BACK_LEFT:
        case TOP_LEFT:

        case DOWN_LEFT: {
            cell_data = (struct basic_cell_data *)node->neighbours[LEFT];
            break;
        }

        case FRONT_TOP_RIGHT:
        case BACK_TOP_RIGHT: {
            cell_data = (struct basic_cell_data *)get_cell_neighbour_with_same_refinement_level(node, TOP_RIGHT);
            break;
        }
        case FRONT_TOP_LEFT:
        case BACK_TOP_LEFT: {
            cell_data = (struct basic_cell_data *)get_cell_neighbour_with_same_refinement_level(node, TOP_LEFT);
            break;
        }
        case FRONT_DOWN_RIGHT:
        case BACK_DOWN_RIGHT: {
            cell_data = (struct basic_cell_data *)get_cell_neighbour_with_same_refinement_level(node, DOWN_RIGHT);
            break;
        }
        case FRONT_DOWN_LEFT:
        case BACK_DOWN_LEFT: {
            cell_data = (struct basic_cell_data *)get_cell_neighbour_with_same_refinement_level(node, DOWN_LEFT);
            break;
        }
        }

        if(!cell_data || cell_data->type == TRANSITION_NODE) {
            return NULL;
        } else {
            return (struct cell_node *)cell_data;
        }

    } else {
        fprintf(stderr, "get_cell_neighbour_with_same_refinement_level() invalid direction %d. Exiting!", direction);
        exit(10);
    }
}

int find_neighbour_index(struct cell_node *grid_cell, struct cell_node *neighbour) {

    size_t max_elements = arrlen(grid_cell->elements);
    int el_index = -1;

    for(size_t i = 0; i < max_elements; i++) {
        if(grid_cell->elements[i].column == neighbour->grid_position) {
            el_index = i;
            break;
        }
    }

    return el_index;
}

uint8_t get_visibility_mask(struct cell_node *grid_cell) {

    uint8_t visibility_mask = 0;

    if(!cell_has_neighbour(grid_cell, grid_cell->neighbours[FRONT])) {
        visibility_mask |= FRONT_IS_VISIBLE;
    }

    if(!cell_has_neighbour(grid_cell, grid_cell->neighbours[BACK])) {
        visibility_mask |= BACK_IS_VISIBLE;
    }

    if(!cell_has_neighbour(grid_cell, grid_cell->neighbours[TOP])) {
        visibility_mask |= TOP_IS_VISIBLE;
    }

    if(!cell_has_neighbour(grid_cell, grid_cell->neighbours[DOWN])) {
        visibility_mask |= DOWN_IS_VISIBLE;
    }

    if(!cell_has_neighbour(grid_cell, grid_cell->neighbours[RIGHT])) {
        visibility_mask |= RIGHT_IS_VISIBLE;
    }

    if(!cell_has_neighbour(grid_cell, grid_cell->neighbours[LEFT])) {
        visibility_mask |= LEFT_IS_VISIBLE;
    }

    return visibility_mask;
}
