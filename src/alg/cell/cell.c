//
// Created by sachetto on 29/09/17.
//

#include "cell.h"
#include <math.h>

#include "../../3dparty/stb_ds.h"

void init_basic_cell_data_with_type(struct basic_cell_data *data, char type) {
    data->type = type;
    data->level = 1;
}

struct cell_node *new_cell_node() {
    struct cell_node *result = (struct cell_node *)malloc(sizeof(struct cell_node));
    init_cell_node(result);
    return result;
}

void init_cell_node(struct cell_node *cell_node) {

    init_basic_cell_data_with_type(&(cell_node->cell_data), CELL_NODE_TYPE);

    cell_node->center.x = 0.0;
    cell_node->center.y = 0.0;
    cell_node->center.z = 0.0;

    cell_node->active = true;

    cell_node->bunch_number = 0;

    cell_node->z_front = NULL;
    cell_node->z_back = NULL;
    cell_node->y_top = NULL;
    cell_node->y_down = NULL;
    cell_node->x_right = NULL;
    cell_node->x_left = NULL;

    cell_node->previous = NULL;
    cell_node->next = NULL;

    cell_node->grid_position = 0;
    cell_node->sv_position = 0;
    cell_node->hilbert_shape_number = 0;
    cell_node->discretization.x = 1.0;
    cell_node->discretization.y = 1.0;
    cell_node->discretization.z = 1.0;

    cell_node->v = 0;

    cell_node->z_front_flux = 0.0;
    cell_node->z_back_flux = 0.0;
    cell_node->y_top_flux = 0.0;
    cell_node->y_down_flux = 0.0;
    cell_node->x_right_flux = 0.0;
    cell_node->x_left_flux = 0.0;

    cell_node->b = 0.0;

    cell_node->can_change = true;
    cell_node->visited = false;

    cell_node->elements = NULL;

    cell_node->linear_system_solver_extra_info = NULL;
    cell_node->mesh_extra_info = NULL;

    cell_node->sigma.x = 0.0;
    cell_node->sigma.y = 0.0;
    cell_node->sigma.z = 0.0;

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

    init_basic_cell_data_with_type(&(transition_node->cell_data), TRANSITION_NODE_TYPE);

    transition_node->single_connector = NULL;
    transition_node->quadruple_connector1 = NULL;
    transition_node->quadruple_connector2 = NULL;
    transition_node->quadruple_connector3 = NULL;
    transition_node->quadruple_connector4 = NULL;
    transition_node->direction = 0;
}

void set_transition_node_data(struct transition_node *the_transition_node, uint16_t level, char direction,
                              void *single_connector, void *quadruple_connector1, void *quadruple_connector2,
                              void *quadruple_connector3, void *quadruple_connector4) {

    the_transition_node->cell_data.level = level;

    the_transition_node->direction = direction;

    the_transition_node->single_connector = single_connector;

    the_transition_node->quadruple_connector1 = quadruple_connector1;
    the_transition_node->quadruple_connector2 = quadruple_connector2;
    the_transition_node->quadruple_connector3 = quadruple_connector3;
    the_transition_node->quadruple_connector4 = quadruple_connector4;

}

void set_cell_node_data(struct cell_node *the_cell, struct point_3d discretization,
                        uint64_t bunch_number, void *y_top, void *z_front, void *y_down, void *z_back,
                        void *x_rigth, void *x_left, void *previous, void *next,
                        uint32_t grid_position, uint8_t hilbert_shape_number, struct point_3d center) {

    the_cell->discretization = discretization;
    the_cell->bunch_number = bunch_number;
    the_cell->y_top = y_top;
    the_cell->z_front = z_front;
    the_cell->y_down = y_down;
    the_cell->z_back = z_back;
    the_cell->x_right = x_rigth;
    the_cell->x_left = x_left;
    the_cell->previous = previous;
    the_cell->next = next;
    the_cell->grid_position = grid_position;
    the_cell->hilbert_shape_number = hilbert_shape_number;
    the_cell->center = center;

}

void set_cell_flux(struct cell_node *the_cell, char direction) {

    void *neighbour_grid_cell;
    struct transition_node *white_neighbor_cell;
    struct cell_node *black_neighbor_cell;

    switch(direction) {
    case 'n':
        neighbour_grid_cell = the_cell->z_front;
        break;

    case 's':
        neighbour_grid_cell = the_cell->z_back;
        break;

    case 'e':
        neighbour_grid_cell = the_cell->y_top;
        break;

    case 'w':
        neighbour_grid_cell = the_cell->y_down;
        break;

    case 'f':
        neighbour_grid_cell = the_cell->x_right;
        break;

    case 'b':
        neighbour_grid_cell = the_cell->x_left;
        break;
    default:
        fprintf(stderr, "Invalid cell direction %c! Exiting...", direction);
        exit(10);
    }

    real_cpu least_distance_x = the_cell->discretization.x/2.0;
    real_cpu least_distance_y = the_cell->discretization.y/2.0;
    real_cpu least_distance_z = the_cell->discretization.z/2.0;

    real_cpu local_flux_x;
    real_cpu local_flux_y;
    real_cpu local_flux_z;

    bool has_found;

    /* When neighbour_grid_cell is a transition node, looks for the next neighbor
     * cell which is a cell node. */
    uint16_t neighbour_grid_cell_level = ((struct basic_cell_data *)(neighbour_grid_cell))->level;
    char neighbour_grid_cell_type = ((struct basic_cell_data *)(neighbour_grid_cell))->type;

    if(neighbour_grid_cell_level > the_cell->cell_data.level) {
        if(neighbour_grid_cell_type == TRANSITION_NODE_TYPE) {
            has_found = false;
            while(!has_found) {
                if(neighbour_grid_cell_type == TRANSITION_NODE_TYPE) {
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
        if(neighbour_grid_cell_level <= the_cell->cell_data.level && (neighbour_grid_cell_type == 'w')) {
            has_found = false;
            while(!has_found) {
                if(neighbour_grid_cell_type == TRANSITION_NODE_TYPE) {
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

    if(neighbour_grid_cell_type == CELL_NODE_TYPE) {

        black_neighbor_cell = (struct cell_node *)(neighbour_grid_cell);

        if(black_neighbor_cell->discretization.x/2.0 < least_distance_x)
            least_distance_x = black_neighbor_cell->discretization.x / 2.0;

        if(black_neighbor_cell->discretization.y/2.0 < least_distance_y)
            least_distance_y = black_neighbor_cell->discretization.y / 2.0;

        if(black_neighbor_cell->discretization.z/2.0 < least_distance_z)
            least_distance_z = black_neighbor_cell->discretization.z / 2.0;

        local_flux_x = (the_cell->v - black_neighbor_cell->v) * (2.0 * least_distance_x);
        local_flux_y = (the_cell->v - black_neighbor_cell->v) * (2.0 * least_distance_y);
        local_flux_z = (the_cell->v - black_neighbor_cell->v) * (2.0 * least_distance_z);

        lock_cell_node(the_cell);

        switch(direction) {
            case 's':
                if(local_flux_y > the_cell->z_back_flux)
                    the_cell->z_back_flux += local_flux_y;
                break;

            case 'n':
                if(local_flux_y > the_cell->z_front_flux)
                    the_cell->z_front_flux += local_flux_y;
                break;

            case 'e':
                if(local_flux_x > the_cell->y_top_flux)
                    the_cell->y_top_flux += local_flux_x;
                break;

            case 'w':
                if(local_flux_x > the_cell->y_down_flux)
                    the_cell->y_down_flux += local_flux_x;
                break;

            case 'f':
                if(local_flux_z > the_cell->x_right_flux)
                    the_cell->x_right_flux += local_flux_z;
                break;

            case 'b':
                if(local_flux_z > the_cell->x_left_flux)
                    the_cell->x_left_flux += local_flux_z;
                break;

            default:
                break;
        }

        unlock_cell_node(the_cell);

        lock_cell_node(black_neighbor_cell);

        switch(direction) {
            case 's':
                if(local_flux_y > black_neighbor_cell->z_front_flux)
                    black_neighbor_cell->z_front_flux += local_flux_y;
                break;

            case 'n':
                if(local_flux_y > black_neighbor_cell->z_back_flux)
                    black_neighbor_cell->z_back_flux += local_flux_y;
                break;

            case 'e':
                if(local_flux_x > black_neighbor_cell->y_down_flux)
                    black_neighbor_cell->y_down_flux += local_flux_x;
                break;

            case 'w':
                if(local_flux_x > black_neighbor_cell->y_top_flux)
                    black_neighbor_cell->y_top_flux += local_flux_x;
                break;

            case 'f':
                if(local_flux_z > black_neighbor_cell->x_left_flux)
                    black_neighbor_cell->x_left_flux += local_flux_z;
                break;

            case 'b':
                if(local_flux_z > black_neighbor_cell->x_right_flux)
                    black_neighbor_cell->x_right_flux += local_flux_z;
                break;

             default:
                 break;
        }

        unlock_cell_node(black_neighbor_cell);
    }
}

real_cpu get_cell_maximum_flux(struct cell_node *the_cell) {

    real_cpu maximumFlux = fabs(the_cell->y_top_flux);
    if(fabs(the_cell->z_front_flux) > maximumFlux)
        maximumFlux = fabs(the_cell->z_front_flux);

    if(fabs(the_cell->z_back_flux) > maximumFlux)
        maximumFlux = fabs(the_cell->z_back_flux);

    if(fabs(the_cell->y_down_flux) > maximumFlux)
        maximumFlux = fabs(the_cell->y_down_flux);

    if(fabs(the_cell->x_right_flux) > maximumFlux)
        maximumFlux = fabs(the_cell->x_right_flux);

    if(fabs(the_cell->x_left_flux) > maximumFlux)
        maximumFlux = fabs(the_cell->x_left_flux);

    return maximumFlux;
}

struct cell_node *get_cell_neighbour(struct cell_node *grid_cell, void *neighbour_grid_cell) {

	struct cell_node *black_neighbor_cell;

	uint16_t neighbour_grid_cell_level = ((struct basic_cell_data *)(neighbour_grid_cell))->level;
	char neighbour_grid_cell_type = ((struct basic_cell_data *)(neighbour_grid_cell))->type;

	struct transition_node *white_neighbor_cell;

	if(neighbour_grid_cell_type == TRANSITION_NODE_TYPE) {

		while(true) {

			if(neighbour_grid_cell_type == TRANSITION_NODE_TYPE) {
				white_neighbor_cell = (struct transition_node *)neighbour_grid_cell;
				if(white_neighbor_cell->single_connector == NULL) {
					return NULL;
				}
				else {

					if(neighbour_grid_cell_level > grid_cell->cell_data.level) {
						neighbour_grid_cell = white_neighbor_cell->quadruple_connector1;
					}
					else {
						neighbour_grid_cell = white_neighbor_cell->single_connector;
					}

					neighbour_grid_cell_type = ((struct basic_cell_data *)(neighbour_grid_cell))->type;
				}
			}
			else {

				black_neighbor_cell = (struct cell_node *)(neighbour_grid_cell);
				if(black_neighbor_cell->active) {
					return black_neighbor_cell;
				}

				else {
					return NULL;
				}
			}
		}
	}
	else {
		black_neighbor_cell = (struct cell_node *)(neighbour_grid_cell);
		if(black_neighbor_cell->active) {
			return black_neighbor_cell;
		}
		else {
			return NULL;
		}
	}

}

bool cell_has_neighbour(struct cell_node *grid_cell, void *neighbour_grid_cell) {

    struct cell_node *black_neighbor_cell = get_cell_neighbour(grid_cell, neighbour_grid_cell);

    if(black_neighbor_cell) {
        return true;
    }
    else {
        return false;
    }

}

void *get_cell_neighbour_as_void(struct cell_node *grid_cell, void *neighbour_grid_cell,  char *type) {

	struct cell_node *black_neighbor_cell;

	uint16_t neighbour_grid_cell_level = ((struct basic_cell_data *)(neighbour_grid_cell))->level;
	char neighbour_grid_cell_type = ((struct basic_cell_data *)(neighbour_grid_cell))->type;

	struct transition_node *white_neighbor_cell;
	void *last_visited_node;

	if(neighbour_grid_cell_type == TRANSITION_NODE_TYPE) {

		while(true) {

			if(neighbour_grid_cell_type == TRANSITION_NODE_TYPE) {

				white_neighbor_cell = (struct transition_node *)neighbour_grid_cell;
				
				if(white_neighbor_cell->single_connector == NULL) {
					*type = 'v';					
					return NULL;
				}
				else {

					last_visited_node = neighbour_grid_cell;

					if(neighbour_grid_cell_level > grid_cell->cell_data.level) {
						neighbour_grid_cell = white_neighbor_cell->quadruple_connector1;
					}
					else {
						neighbour_grid_cell = white_neighbor_cell->single_connector;
					}

					neighbour_grid_cell_type = ((struct basic_cell_data *)(neighbour_grid_cell))->type;
				}
			}
			else {
				*type = TRANSITION_NODE_TYPE;
				return last_visited_node;
			}
		}
	}
	else {
		black_neighbor_cell = (struct cell_node *)(neighbour_grid_cell);
		if(black_neighbor_cell->active) {
			*type = CELL_NODE_TYPE;
			return black_neighbor_cell;
		}
		else {
			*type = 'v';
			return NULL;
		}
	}

}
