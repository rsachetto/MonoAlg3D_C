//
// Created by sachetto on 30/09/17.
//

#include "grid.h"
#include "../../single_file_libraries/stb_ds.h"

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
        fprintf(stderr, "refine_grid(): Parameter min_dx must be positive, passed %lf.", min_dx);
        return false;
    }

    if(min_dy <= 0.0) {
        fprintf(stderr, "refine_grid(): Parameter min_dy must be positive, passed %lf.", min_dy);
        return false;
    }

    if(min_dz <= 0.0) {
        fprintf(stderr, "refine_grid(): Parameter min_dz must be positive, passed %lf.", min_dz);
        return false;
    }

    struct cell_node *grid_cell, *auxiliar_grid_cell;

    real_cpu maximum_flux;
    bool continue_refining = true;
    bool refined_once = false;
    set_grid_flux(the_grid);

    uint32_t *free_sv_pos = the_grid->free_sv_positions;

    arrreset(the_grid->refined_this_step);

    while(continue_refining) {
        continue_refining = false;
        grid_cell = the_grid->first_cell;
        while(grid_cell != 0) {

            maximum_flux = get_cell_maximum_flux(grid_cell);

            if((grid_cell->can_change && grid_cell->active) && (grid_cell->dx > min_dx) && (grid_cell->dx > min_dy) &&
               (grid_cell->dx > min_dz) && (maximum_flux >= refinement_bound)) {
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
        fprintf(stderr, "refine_grid(): Parameter the_grid can't be null. Exiting!");
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
    }
}

void refine_grid_cell(struct grid *the_grid, struct cell_node *grid_cell) {

    if(!grid_cell) {
        fprintf(stderr, "refine_grid_cell: grid_cell is NULL.\n");
        exit(10);
    }

    refine_cell(grid_cell, NULL, NULL);
    the_grid->number_of_cells += 7;
}

void set_grid_flux(struct grid *the_grid) {

    uint32_t active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    int i;

#pragma omp parallel for
    for(i = 0; i < active_cells; i++) {
        ac[i]->north_flux = 0.0;
        ac[i]->south_flux = 0.0;
        ac[i]->east_flux = 0.0;
        ac[i]->west_flux = 0.0;
        ac[i]->front_flux = 0.0;
        ac[i]->back_flux = 0.0;
    }

#pragma omp parallel for
    for(i = 0; i < active_cells; i++) {
        set_cell_flux(ac[i], 's'); // Computes south flux.
        set_cell_flux(ac[i], 'n'); // Computes north flux.
        set_cell_flux(ac[i], 'e'); // Computes east flux.
        set_cell_flux(ac[i], 'w'); // Computes west flux.
        set_cell_flux(ac[i], 'f'); // Computes front flux.
        set_cell_flux(ac[i], 'b'); // Computes back flux.
    }
}