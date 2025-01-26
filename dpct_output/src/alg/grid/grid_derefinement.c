//
// Created by sachetto on 30/09/17.
//

#include "grid.h"

static bool can_derefine(struct cell_node *grid_cell) {

    if (grid_cell->cell_data.level > 1) {

        if ((grid_cell->next != 0) &&
                (grid_cell->next->next != 0) &&
                (grid_cell->next->next->next != 0) &&
                (grid_cell->next->next->next->next != 0) &&
                (grid_cell->next->next->next->next->next != 0) &&
                (grid_cell->next->next->next->next->next->next != 0) &&
                (grid_cell->next->next->next->next->next->next->next != 0)) {

            if ((grid_cell->cell_data.level == grid_cell->next->cell_data.level) &&
                (grid_cell->cell_data.level == grid_cell->next->next->cell_data.level) &&
                (grid_cell->cell_data.level == grid_cell->next->next->next->cell_data.level) &&
                (grid_cell->cell_data.level == grid_cell->next->next->next->next->cell_data.level) &&
                (grid_cell->cell_data.level == grid_cell->next->next->next->next->next->cell_data.level) &&
                (grid_cell->cell_data.level == grid_cell->next->next->next->next->next->next->cell_data.level) &&
                (grid_cell->cell_data.level == grid_cell->next->next->next->next->next->next->next->cell_data.level)) {

                uint64_t bunch_number1 = grid_cell->bunch_number / 10;
                uint64_t bunch_number2 = grid_cell->next->bunch_number / 10;
                uint64_t bunch_number3 = grid_cell->next->next->bunch_number / 10;
                uint64_t bunch_number4 = grid_cell->next->next->next->bunch_number / 10;
                uint64_t bunch_number5 = grid_cell->next->next->next->next->bunch_number / 10;
                uint64_t bunch_number6 = grid_cell->next->next->next->next->next->bunch_number / 10;
                uint64_t bunch_number7 = grid_cell->next->next->next->next->next->next->bunch_number / 10;
                uint64_t bunch_number8 = grid_cell->next->next->next->next->next->next->next->bunch_number / 10;

                if ((bunch_number1 == bunch_number2) &&
                    (bunch_number1 == bunch_number3) &&
                    (bunch_number1 == bunch_number4) &&
                    (bunch_number1 == bunch_number5) &&
                    (bunch_number1 == bunch_number6) &&
                    (bunch_number1 == bunch_number7) &&
                    (bunch_number1 == bunch_number8)) {

                    return true;

                }
            }
        }
    }

    return false;

}

/**
 * Decides if the grid should be derefined by traversing the whole grid,
 * according to the parameter derefinement_bound.
 *
 * @param derefinement_bound Derefinement condition.
 */
bool derefine_grid_with_bound (struct grid *the_grid, real_cpu derefinement_bound, real_cpu max_dx, real_cpu max_dy, real_cpu max_dz) {

    struct cell_node *grid_cell, *auxiliar_grid_cell;

    uint64_t bunch_number1, bunch_number2, bunch_number3, bunch_number4, bunch_number5, bunch_number6, bunch_number7,
             bunch_number8;

    bool active1, active2, active3, active4, active5, active6, active7, active8;

    bool has_been_derefined;
    bool derefined_once = false;
    set_grid_flux (the_grid);

    grid_cell = the_grid->first_cell;

    while (grid_cell != 0) {
        has_been_derefined = false;
        if (grid_cell->can_change && grid_cell->discretization.x < max_dx && grid_cell->discretization.y < max_dy && grid_cell->discretization.z < max_dz && grid_cell->active) {

            if ((grid_cell->next != 0) &&
                    (grid_cell->next->next != 0) &&
                    (grid_cell->next->next->next != 0) &&
                    (grid_cell->next->next->next->next != 0) &&
                    (grid_cell->next->next->next->next->next != 0) &&
                    (grid_cell->next->next->next->next->next->next != 0) &&
                    (grid_cell->next->next->next->next->next->next->next != 0)) {


                /* Verifies if each one of the next seven cells has  are
                 * active */
                active1 = grid_cell->active;
                active2 = grid_cell->next->active;
                active3 = grid_cell->next->next->active;
                active4 = grid_cell->next->next->next->active;
                active5 = grid_cell->next->next->next->next->active;
                active6 = grid_cell->next->next->next->next->next->active;
                active7 = grid_cell->next->next->next->next->next->next->active;
                active8 = grid_cell->next->next->next->next->next->next->next->active;

                if (active1 && active2 && active3 && active4 && active5 && active6 && active7 && active8) {

                    /* Checks if the next seven cells of the current cell can change. */
                    if ((grid_cell->next->can_change) &&
                            (grid_cell->next->next->can_change) &&
                            (grid_cell->next->next->next->can_change) &&
                            (grid_cell->next->next->next->next->can_change) &&
                            (grid_cell->next->next->next->next->next->can_change) &&
                            (grid_cell->next->next->next->next->next->next->can_change) &&
                            (grid_cell->next->next->next->next->next->next->next->can_change)) {
                        /* Verifies if each one of the next seven  cells  has  the  same
                         * refinement level. */
                        if ((grid_cell->cell_data.level == grid_cell->next->cell_data.level) &&
                                (grid_cell->cell_data.level == grid_cell->next->next->cell_data.level) &&
                                (grid_cell->cell_data.level == grid_cell->next->next->next->cell_data.level) &&
                                (grid_cell->cell_data.level == grid_cell->next->next->next->next->cell_data.level) &&
                                (grid_cell->cell_data.level == grid_cell->next->next->next->next->next->cell_data.level) &&
                                (grid_cell->cell_data.level == grid_cell->next->next->next->next->next->next->cell_data.level) &&
                                (grid_cell->cell_data.level == grid_cell->next->next->next->next->next->next->next->cell_data.level)) {
                            /* Checks if this cell and the next seven cells belong to the
                             * same bunch. */
                            bunch_number1 = grid_cell->bunch_number / 10;
                            bunch_number2 = grid_cell->next->bunch_number / 10;
                            bunch_number3 = grid_cell->next->next->bunch_number / 10;
                            bunch_number4 = grid_cell->next->next->next->bunch_number / 10;
                            bunch_number5 = grid_cell->next->next->next->next->bunch_number / 10;
                            bunch_number6 = grid_cell->next->next->next->next->next->bunch_number / 10;
                            bunch_number7 = grid_cell->next->next->next->next->next->next->bunch_number / 10;
                            bunch_number8 = grid_cell->next->next->next->next->next->next->next->bunch_number / 10;
                            if ((bunch_number1 == bunch_number2) && (bunch_number1 == bunch_number3) &&
                                    (bunch_number1 == bunch_number4) && (bunch_number1 == bunch_number5) &&
                                    (bunch_number1 == bunch_number6) && (bunch_number1 == bunch_number7) &&
                                    (bunch_number1 == bunch_number8)) {

                                /* Notice that by the ordering conferred by the  Hilbert
                                 * curve, the program always  enters  the  bunch  to  be
                                 * derefined by its first member cell.
                                 */
                                if (cell_needs_derefinement (grid_cell, derefinement_bound)) {
                                    auxiliar_grid_cell = grid_cell->next->next->next->next->next->next->next->next;
                                    derefine_cell_bunch (grid_cell, &(the_grid->free_sv_positions));

                                    the_grid->number_of_cells -= 7;
                                    grid_cell = auxiliar_grid_cell;
                                    has_been_derefined = true;
                                    derefined_once = true;
                                }
                            }
                        }
                    }
                }
            }
        }

        if (has_been_derefined == false)
            grid_cell = grid_cell->next;
    }
    return derefined_once;
}

/**
 * Derefines whole grid one level below the current refinement level.
 */
void derefine_all_grid (struct grid *the_grid) {

    if(!the_grid) return;

    struct cell_node *grid_cell, *auxiliar_grid_cell;

    grid_cell = the_grid->first_cell;

    while (grid_cell != 0) {
        bool has_been_derefined = false;
        if(can_derefine(grid_cell)) {
            auxiliar_grid_cell = grid_cell->next->next->next->next->next->next->next->next;
            derefine_cell_bunch (grid_cell, NULL);
            the_grid->number_of_cells -= 7;
            grid_cell = auxiliar_grid_cell;
            has_been_derefined = true;
        }
        if (!has_been_derefined) {
            grid_cell = grid_cell->next;
        }
    }
}

/**
 * Derefines all inactive cells of grid.
 */
void derefine_grid_inactive_cells (struct grid *the_grid) {

    struct cell_node *grid_cell, *auxiliar_grid_cell;

    bool active1, active2, active3, active4, active5, active6, active7, active8;

    grid_cell = the_grid->first_cell;

    while (grid_cell != 0) {
        bool has_been_derefined = false;

        if(can_derefine(grid_cell)) {

            /* Verifies if each one of the next seven cells has  are
             * inactive */
            active1 = grid_cell->active;
            active2 = grid_cell->next->active;
            active3 = grid_cell->next->next->active;
            active4 = grid_cell->next->next->next->active;
            active5 = grid_cell->next->next->next->next->active;
            active6 = grid_cell->next->next->next->next->next->active;
            active7 = grid_cell->next->next->next->next->next->next->active;
            active8 = grid_cell->next->next->next->next->next->next->next->active;

            if ((!active1) && (!active2) && (!active3) && (!active4) && (!active5) && (!active6) &&
                    (!active7) && (!active8)) {

                auxiliar_grid_cell = grid_cell->next->next->next->next->next->next->next->next;

                derefine_cell_bunch (grid_cell, NULL);
                the_grid->number_of_cells -= 7;
                grid_cell = auxiliar_grid_cell;
                has_been_derefined = true;
            }
        }

        if (!has_been_derefined)
            grid_cell = grid_cell->next;
    }
}
