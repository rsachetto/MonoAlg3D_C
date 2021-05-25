//
// Created by sachetto on 30/09/17.
//

#include "../../3dparty/stb_ds.h"
#include "../../logger/logger.h"
#include "grid.h"

/**
 * Decides if the grid should be refined by traversing the whole grid, according
 * to parameters refinementLevel and refinementBound. A cell will not be refined
 * either if its refinement level  is  equal  to refinementLevel or the  highest
 * of  all  fluxes  coming  into  it  from  the  six  directions  is  less  than
 * refinementBound.
 *
 * @param min_h Minimum refinement level required for the graph.
 * @param refinement_bound Minimum flux required for each cell of graph.
 */
bool refine_grid_with_bound(struct grid *the_grid, real_cpu refinement_bound, real_cpu min_dx, real_cpu min_dy,
                            real_cpu min_dz) {

    if(min_dx <= 0.0) {
        log_error("refine_grid(): Parameter min_dx must be positive, passed %lf.", min_dx);
        return false;
    }

    if(min_dy <= 0.0) {
        log_error("refine_grid(): Parameter min_dy must be positive, passed %lf.", min_dy);
        return false;
    }

    if(min_dz <= 0.0) {
        log_error("refine_grid(): Parameter min_dz must be positive, passed %lf.", min_dz);
        return false;
    }

    struct cell_node *grid_cell, *auxiliar_grid_cell;

    real_cpu maximum_flux;
    bool continue_refining = true;
    bool refined_once = false;
    set_grid_flux(the_grid);

    uint32_t *free_sv_pos = the_grid->free_sv_positions;

    arrsetlen(the_grid->refined_this_step,0);

    while(continue_refining) {
        continue_refining = false;
        grid_cell = the_grid->first_cell;
        while(grid_cell != 0) {

            maximum_flux = get_cell_maximum_flux(grid_cell);

            if((grid_cell->can_change && grid_cell->active) && (grid_cell->discretization.x > min_dx) && (grid_cell->discretization.y > min_dy) &&
               (grid_cell->discretization.z > min_dz) && (maximum_flux >= refinement_bound)) {
                auxiliar_grid_cell = grid_cell;
                grid_cell = grid_cell->next;
                refine_cell(auxiliar_grid_cell, free_sv_pos, &(the_grid->refined_this_step));
                the_grid->number_of_cells += 7;
                continue_refining = true;
                refined_once = true;
            } else {
                grid_cell = grid_cell->next;
            }
        }
    }
    return refined_once;
}

void refine_grid(struct grid *the_grid, int num_steps) {

    if(the_grid == NULL) {
        log_error("refine_grid(): Parameter the_grid can't be null. Exiting!");
        exit(10);
    }

    struct cell_node *grid_cell, *auxiliar_grid_cell;

    for(int i = 0; i < num_steps; i++) {
        grid_cell = the_grid->first_cell;
        while(grid_cell != 0) {
            if(grid_cell->can_change && grid_cell->active) {
                auxiliar_grid_cell = grid_cell;
                grid_cell = grid_cell->next;
                refine_cell(auxiliar_grid_cell, NULL, NULL);
                the_grid->number_of_cells += 7;
            } else {
                grid_cell = grid_cell->next;
            }
        }
        log_info("Refined %d of %d (%ld cells)\n", i+1, num_steps, the_grid->number_of_cells);
    }
}

void refine_grid_with_bounds(struct grid *the_grid, int num_steps, struct point_3d min_bounds, struct point_3d max_bounds) {

    if(the_grid == NULL) {
        log_error("refine_grid(): Parameter the_grid can't be null. Exiting!");
        exit(10);
    }

    struct cell_node *grid_cell, *auxiliar_grid_cell;

    real_cpu min_x = min_bounds.x;
    real_cpu min_y = min_bounds.y;
    real_cpu min_z = min_bounds.z;

    real_cpu max_x = max_bounds.x;
    real_cpu max_y = max_bounds.y;
    real_cpu max_z = max_bounds.z;

    for(int i = 0; i < num_steps; i++) {
        grid_cell = the_grid->first_cell;
        while(grid_cell != 0) {

            real_cpu center_x = grid_cell->center.x,
                     center_y = grid_cell->center.y,
                     center_z = grid_cell->center.z;

            bool refine = center_x >= min_x && center_y >= min_y && center_z >= min_z;
            refine     &= center_x <= max_x && center_y <= max_y && center_z <= max_z;

            if(!refine) {
                refine = true;
                real_cpu x_limit;
                real_cpu y_limit;
                real_cpu z_limit;

                if(center_x < min_x) {
                    x_limit = center_x + grid_cell->discretization.x / 2;
                    refine &= min_x <= x_limit;
                }

                if (center_x > max_x) {
                    x_limit = center_x - grid_cell->discretization.x / 2;
                    refine &= max_x >= x_limit;
                }

               if(center_y < min_y) {
                    y_limit = center_y + grid_cell->discretization.y / 2;
                    refine &= min_y <= y_limit;
                }

                if (center_y > max_y) {
                    y_limit = center_y - grid_cell->discretization.y / 2;
                    refine &= max_y >= y_limit;
                }


               if(center_z < min_z) {
                    z_limit = center_z + grid_cell->discretization.z / 2;
                    refine &= min_z <= z_limit;
                }

                if (center_z > max_z) {
                    z_limit = center_z - grid_cell->discretization.z / 2;
                    refine &= max_z >= z_limit;
                }

            }

            if(grid_cell->can_change && grid_cell->active && refine) {
                auxiliar_grid_cell = grid_cell;
                grid_cell = grid_cell->next;
                refine_cell(auxiliar_grid_cell, NULL, NULL);
                the_grid->number_of_cells += 7;
            } else {
                grid_cell = grid_cell->next;
            }
        }
        log_info("Refined %d of %d (%ld cells)\n", i+1, num_steps, the_grid->number_of_cells);
    }
}

void refine_grid_cell(struct grid *the_grid, struct cell_node *grid_cell) {

    if(!grid_cell) {
        log_error("refine_grid_cell: grid_cell is NULL.\n");
        exit(10);
    }

    refine_cell(grid_cell, NULL, NULL);
    the_grid->number_of_cells += 7;
}

void set_grid_flux(struct grid *the_grid) {

    uint32_t active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    int i;

    OMP(parallel for)
    for(i = 0; i < active_cells; i++) {
        ac[i]->front_flux = 0.0;
        ac[i]->back_flux = 0.0;
        ac[i]->top_flux = 0.0;
        ac[i]->down_flux = 0.0;
        ac[i]->right_flux = 0.0;
        ac[i]->left_flux = 0.0;
    }

    OMP(parallel for)
    for(i = 0; i < active_cells; i++) {
        set_cell_flux(ac[i], BACK); // Computes south flux.
        set_cell_flux(ac[i], FRONT); // Computes north flux.
        set_cell_flux(ac[i], TOP); // Computes east flux.
        set_cell_flux(ac[i], DOWN); // Computes west flux.
        set_cell_flux(ac[i], RIGHT); // Computes front flux.
        set_cell_flux(ac[i], LEFT); // Computes back flux.
    }
}
