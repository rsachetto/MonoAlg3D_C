//
// Created by bergolho on 29/09/20.
//

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../3dparty/sds/sds.h"
#include "../3dparty/stb_ds.h"
#include "../alg/grid/grid.h"
#include "../config/assembly_matrix_config.h"
#include "../libraries_common/common_data_structures.h"
#include "../utils/file_utils.h"
#include "../utils/utils.h"

INIT_ASSEMBLY_MATRIX(set_initial_conditions_coupling_fvm) {

    real_cpu alpha;

    // Tissue parameters
    struct cell_node **ac = the_grid->active_cells;
    uint32_t active_cells = the_grid->num_active_cells;

    // Purkinje parameters
    struct cell_node **ac_purkinje = the_grid->purkinje->purkinje_cells;
    uint32_t active_purkinje_cells = the_grid->purkinje->num_active_purkinje_cells;

    // Common parameters
    real_cpu beta = the_solver->beta;
    real_cpu cm = the_solver->cm;
    real_cpu dt = the_solver->dt;
    uint32_t i;

    // Tissue section
    OMP(parallel for private(alpha))
    for(i = 0; i < active_cells; i++)
    {

        alpha = ALPHA(beta, cm, dt, ac[i]->discretization.x, ac[i]->discretization.y, ac[i]->discretization.z);
        ac[i]->v = initial_v;
        ac[i]->b = initial_v * alpha;
    }

    // Purkinje section
    OMP(parallel for private(alpha))
    for(i = 0; i < active_purkinje_cells; i++) {

        alpha = ALPHA(beta, cm, dt, ac[i]->discretization.x, ac[i]->discretization.y, ac[i]->discretization.z);
        ac_purkinje[i]->v = purkinje_initial_v;
        ac_purkinje[i]->b = purkinje_initial_v * alpha;
    }
}

static struct element fill_element(uint32_t position, enum transition_direction direction, real_cpu dx, real_cpu dy, real_cpu dz,\
                                   real_cpu sigma_x, real_cpu sigma_y, real_cpu sigma_z,\
                                   struct element *cell_elements);


void initialize_diagonal_elements(struct monodomain_solver *the_solver, struct grid *the_grid) {

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;
    real_cpu beta = the_solver->beta;
    real_cpu cm = the_solver->cm;
    real_cpu dt = the_solver->dt;

    uint32_t i;

    OMP(parallel for)
    for(i = 0; i < num_active_cells; i++) {
        real_cpu alpha, dx, dy, dz;

        dx = ac[i]->discretization.x;
        dy = ac[i]->discretization.y;
        dz = ac[i]->discretization.z;

        alpha = ALPHA(beta, cm, dt, dx, dy, dz);

        struct element element;
        element.column = ac[i]->grid_position;
        element.cell = ac[i];
        element.value = alpha;

        if(ac[i]->elements)
            arrfree(ac[i]->elements);

        ac[i]->elements = NULL;

        arrsetcap(ac[i]->elements, 7);
        arrput(ac[i]->elements, element);
    }
}

struct element fill_element(uint32_t position, enum transition_direction direction, real_cpu dx, real_cpu dy, real_cpu dz, real_cpu sigma_x,
                            real_cpu sigma_y, real_cpu sigma_z, struct element *cell_elements) {

    real_cpu multiplier;

    struct element new_element;
    new_element.column = position;

    if(direction == FRONT) { // Z direction
        multiplier = ((dx * dy) / dz);
        new_element.value = -sigma_z * multiplier;
        cell_elements[0].value += (sigma_z * multiplier);
    } else if(direction == BACK) { // Z direction
        multiplier = ((dx * dy) / dz);
        new_element.value = -sigma_z * multiplier;
        cell_elements[0].value += (sigma_z * multiplier);
    } else if(direction == TOP) { // Y direction
        multiplier = ((dx * dz) / dy);
        new_element.value = -sigma_y * multiplier;
        cell_elements[0].value += (sigma_y * multiplier);
    } else if(direction == DOWN) { // Y direction
        multiplier = ((dx * dz) / dy);
        new_element.value = -sigma_y * multiplier;
        cell_elements[0].value += (sigma_y * multiplier);
    } else if(direction == RIGHT) { // X direction
        multiplier = ((dy * dz) / dx);
        new_element.value = -sigma_x * ((dy * dz) / dx);
        cell_elements[0].value += (sigma_x * multiplier);
    } else if(direction == LEFT) { // X direction
        multiplier = ((dy * dz) / dx);
        new_element.value = -sigma_x * multiplier;
        cell_elements[0].value += (sigma_x * multiplier);
    }
    return new_element;
}

static void fill_discretization_matrix_elements(struct cell_node *grid_cell, void *neighbour_grid_cell, enum transition_direction direction) {

    bool has_found;

    struct transition_node *white_neighbor_cell;
    struct cell_node *black_neighbor_cell;

    /* When neighbour_grid_cell is a transition node, looks for the next neighbor
     * cell which is a cell node. */
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
        if(neighbour_grid_cell_level <= grid_cell->cell_data.level &&
           (neighbour_grid_cell_type == TRANSITION_NODE)) {
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

    // We care only with the interior points
    if(neighbour_grid_cell_type == CELL_NODE) {

        black_neighbor_cell = (struct cell_node *)(neighbour_grid_cell);

        if(black_neighbor_cell->active) {

            uint32_t position;
            real_cpu dx, dy, dz;

            real_cpu sigma_x1 = grid_cell->sigma.x;
            real_cpu sigma_x2 = black_neighbor_cell->sigma.x;
            real_cpu sigma_x = 0.0;

            if(sigma_x1 != 0.0 && sigma_x2 != 0.0) {
                sigma_x = (2.0f * sigma_x1 * sigma_x2) / (sigma_x1 + sigma_x2);
            }

            real_cpu sigma_y1 = grid_cell->sigma.y;
            real_cpu sigma_y2 = black_neighbor_cell->sigma.y;
            real_cpu sigma_y = 0.0;

            if(sigma_y1 != 0.0 && sigma_y2 != 0.0) {
                sigma_y = (2.0f * sigma_y1 * sigma_y2) / (sigma_y1 + sigma_y2);
            }

            real_cpu sigma_z1 = grid_cell->sigma.z;
            real_cpu sigma_z2 = black_neighbor_cell->sigma.z;
            real_cpu sigma_z = 0.0;

            if(sigma_z1 != 0.0 && sigma_z2 != 0.0) {
                sigma_z = (2.0f * sigma_z1 * sigma_z2) / (sigma_z1 + sigma_z2);
            }

            if(black_neighbor_cell->cell_data.level > grid_cell->cell_data.level) {
                dx = black_neighbor_cell->discretization.x;
                dy = black_neighbor_cell->discretization.y;
                dz = black_neighbor_cell->discretization.z;
            } else {
                dx = grid_cell->discretization.x;
                dy = grid_cell->discretization.y;
                dz = grid_cell->discretization.z;
            }

            lock_cell_node(grid_cell);

            struct element *cell_elements = grid_cell->elements;
            position = black_neighbor_cell->grid_position;

            size_t max_elements = arrlen(cell_elements);
            bool insert = true;

            for(size_t i = 1; i < max_elements; i++) {
                if(cell_elements[i].column == position) {
                    insert = false;
                    break;
                }
            }

            if(insert) {

                struct element new_element = fill_element(position, direction, dx, dy, dz, sigma_x, sigma_y, sigma_z, cell_elements);

                new_element.cell = black_neighbor_cell;
                arrput(grid_cell->elements, new_element);
            }
            unlock_cell_node(grid_cell);

            lock_cell_node(black_neighbor_cell);
            cell_elements = black_neighbor_cell->elements;
            position = grid_cell->grid_position;

            max_elements = arrlen(cell_elements);

            insert = true;
            for(size_t i = 1; i < max_elements; i++) {
                if(cell_elements[i].column == position) {
                    insert = false;
                    break;
                }
            }

            if(insert) {

                struct element new_element = fill_element(position, direction, dx, dy, dz, sigma_x, sigma_y, sigma_z, cell_elements);

                new_element.cell = grid_cell;
                arrput(black_neighbor_cell->elements, new_element);
            }

            unlock_cell_node(black_neighbor_cell);
        }
    }
}

#define CALC_PARTIAL_SIGMA(dir, neighbours, __sigma__, sigma_count)                                                                                            \
    do {                                                                                                                                                       \
        if(neighbours[0] && neighbours[1]) {                                                                                                                   \
            __sigma__ += neighbours[0]->sigma.dir + neighbours[1]->sigma.dir;                                                                                  \
            sigma_count += 2;                                                                                                                                  \
        }                                                                                                                                                      \
                                                                                                                                                               \
        if(neighbours[2] && neighbours[3]) {                                                                                                                   \
            __sigma__ += neighbours[2]->sigma.dir + neighbours[3]->sigma.dir;                                                                                  \
            sigma_count += 2;                                                                                                                                  \
        }                                                                                                                                                      \
                                                                                                                                                               \
        if(neighbours[4] && neighbours[5]) {                                                                                                                   \
            __sigma__ += neighbours[4]->sigma.dir + neighbours[5]->sigma.dir;                                                                                  \
            sigma_count += 2;                                                                                                                                  \
        }                                                                                                                                                      \
                                                                                                                                                               \
        if(neighbours[6] && neighbours[7]) {                                                                                                                   \
            __sigma__ += neighbours[6]->sigma.dir + neighbours[7]->sigma.dir;                                                                                  \
            sigma_count += 2;                                                                                                                                  \
        }                                                                                                                                                      \
                                                                                                                                                               \
        if(neighbours[8] && neighbours[9]) {                                                                                                                   \
            __sigma__ += neighbours[8]->sigma.dir + neighbours[9]->sigma.dir;                                                                                  \
            sigma_count += 2;                                                                                                                                  \
        }                                                                                                                                                      \
                                                                                                                                                               \
        if(neighbours[10] && neighbours[11]) {                                                                                                                 \
            __sigma__ += neighbours[10]->sigma.dir + neighbours[11]->sigma.dir;                                                                                \
            sigma_count += 2;                                                                                                                                  \
        }                                                                                                                                                      \
                                                                                                                                                               \
        if(sigma_count) {                                                                                                                                      \
            __sigma__ /= sigma_count;                                                                                                                          \
        }                                                                                                                                                      \
    } while(0)


static void calc_sigmas(struct cell_node *cell_node, struct cell_node *neighbours[26], enum transition_direction flux_direction, real_cpu *sigma_1,
                        real_cpu *sigma_2, real_cpu *sigma_3, int *count_s1, int *count_s2, int *count_s3) {

    *sigma_1 = 0.0;
    *sigma_2 = 0.0;
    *sigma_3 = 0.0;

    struct cell_node *sigma_neighbours[12];

    if(flux_direction == RIGHT) {
        if(neighbours[RIGHT]) {
            *sigma_1 = (neighbours[RIGHT]->sigma.x + cell_node->sigma.x) / 2.0;
            *count_s1 = 1;
        }

        sigma_neighbours[0] = neighbours[TOP];
        sigma_neighbours[1] = neighbours[DOWN];

        sigma_neighbours[2] = neighbours[TOP_RIGHT];
        sigma_neighbours[3] = neighbours[DOWN_RIGHT];

        sigma_neighbours[4] = neighbours[FRONT_TOP];
        sigma_neighbours[5] = neighbours[FRONT_DOWN];

        sigma_neighbours[6] = neighbours[FRONT_TOP_RIGHT];
        sigma_neighbours[7] = neighbours[FRONT_DOWN_RIGHT];

        sigma_neighbours[8] = neighbours[BACK_TOP];
        sigma_neighbours[9] = neighbours[BACK_DOWN];

        sigma_neighbours[10] = neighbours[BACK_TOP_RIGHT];
        sigma_neighbours[11] = neighbours[BACK_DOWN_RIGHT];

        CALC_PARTIAL_SIGMA(xy, sigma_neighbours, *sigma_2, *count_s2);

        /////////////////////////////////////////////////////////////

        sigma_neighbours[0] = neighbours[FRONT];
        sigma_neighbours[1] = neighbours[BACK];

        sigma_neighbours[2] = neighbours[FRONT_RIGHT];
        sigma_neighbours[3] = neighbours[BACK_RIGHT];

        sigma_neighbours[4] = neighbours[FRONT_TOP];
        sigma_neighbours[5] = neighbours[BACK_TOP];

        sigma_neighbours[6] = neighbours[FRONT_TOP_RIGHT];
        sigma_neighbours[7] = neighbours[BACK_TOP_RIGHT];

        sigma_neighbours[8] = neighbours[FRONT_DOWN];
        sigma_neighbours[9] = neighbours[BACK_DOWN];

        sigma_neighbours[10] = neighbours[FRONT_DOWN_RIGHT];
        sigma_neighbours[11] = neighbours[BACK_DOWN_RIGHT];

        CALC_PARTIAL_SIGMA(xz, sigma_neighbours, *sigma_3, *count_s3);

    } else if(flux_direction == LEFT) {
        if(neighbours[LEFT]) {
            *sigma_1 = (neighbours[LEFT]->sigma.x + cell_node->sigma.x) / 2.0;
            *count_s1 = 1;
        }

        sigma_neighbours[0] = neighbours[TOP];
        sigma_neighbours[1] = neighbours[DOWN];

        sigma_neighbours[2] = neighbours[TOP_LEFT];
        sigma_neighbours[3] = neighbours[DOWN_LEFT];

        sigma_neighbours[4] = neighbours[FRONT_TOP];
        sigma_neighbours[5] = neighbours[FRONT_DOWN];

        sigma_neighbours[6] = neighbours[FRONT_TOP_LEFT];
        sigma_neighbours[7] = neighbours[FRONT_DOWN_LEFT];

        sigma_neighbours[8] = neighbours[BACK_TOP];
        sigma_neighbours[9] = neighbours[BACK_DOWN];

        sigma_neighbours[10] = neighbours[BACK_TOP_LEFT];
        sigma_neighbours[11] = neighbours[BACK_DOWN_LEFT];

        CALC_PARTIAL_SIGMA(xy, sigma_neighbours, *sigma_2, *count_s2);

        /////////////////////////////////////////////////////////////

        sigma_neighbours[0] = neighbours[FRONT];
        sigma_neighbours[1] = neighbours[BACK];

        sigma_neighbours[2] = neighbours[FRONT_LEFT];
        sigma_neighbours[3] = neighbours[BACK_LEFT];

        sigma_neighbours[4] = neighbours[FRONT_TOP];
        sigma_neighbours[5] = neighbours[BACK_TOP];

        sigma_neighbours[6] = neighbours[FRONT_TOP_LEFT];
        sigma_neighbours[7] = neighbours[BACK_TOP_LEFT];

        sigma_neighbours[8] = neighbours[FRONT_DOWN];
        sigma_neighbours[9] = neighbours[BACK_DOWN];

        sigma_neighbours[10] = neighbours[FRONT_DOWN_LEFT];
        sigma_neighbours[11] = neighbours[BACK_DOWN_LEFT];

        CALC_PARTIAL_SIGMA(xz, sigma_neighbours, *sigma_3, *count_s3);

    } else if(flux_direction == TOP) {
        if(neighbours[TOP]) {
            *sigma_1 = (neighbours[TOP]->sigma.y + cell_node->sigma.y) / 2.0;
            *count_s1 = 1;
        }

        sigma_neighbours[0] = neighbours[RIGHT];
        sigma_neighbours[1] = neighbours[LEFT];

        sigma_neighbours[2] = neighbours[TOP_RIGHT];
        sigma_neighbours[3] = neighbours[TOP_LEFT];

        sigma_neighbours[4] = neighbours[FRONT_RIGHT];
        sigma_neighbours[5] = neighbours[FRONT_LEFT];

        sigma_neighbours[6] = neighbours[FRONT_TOP_RIGHT];
        sigma_neighbours[7] = neighbours[FRONT_TOP_LEFT];

        sigma_neighbours[8] = neighbours[BACK_RIGHT];
        sigma_neighbours[9] = neighbours[BACK_LEFT];

        sigma_neighbours[10] = neighbours[BACK_TOP_RIGHT];
        sigma_neighbours[11] = neighbours[BACK_TOP_LEFT];

        CALC_PARTIAL_SIGMA(xy, sigma_neighbours, *sigma_2, *count_s2);

        /////////////////////////////////////////////////////////////

        sigma_neighbours[0] = neighbours[FRONT];
        sigma_neighbours[1] = neighbours[BACK];

        sigma_neighbours[2] = neighbours[FRONT_TOP];
        sigma_neighbours[3] = neighbours[BACK_TOP];

        sigma_neighbours[4] = neighbours[FRONT_RIGHT];
        sigma_neighbours[5] = neighbours[BACK_RIGHT];

        sigma_neighbours[6] = neighbours[FRONT_TOP_RIGHT];
        sigma_neighbours[7] = neighbours[BACK_TOP_RIGHT];

        sigma_neighbours[8] = neighbours[FRONT_LEFT];
        sigma_neighbours[9] = neighbours[BACK_LEFT];

        sigma_neighbours[10] = neighbours[FRONT_TOP_LEFT];
        sigma_neighbours[11] = neighbours[BACK_TOP_LEFT];

        CALC_PARTIAL_SIGMA(yz, sigma_neighbours, *sigma_3, *count_s3);

    } else if(flux_direction == DOWN) {
        if(neighbours[DOWN]) {
            *sigma_1 = (neighbours[DOWN]->sigma.y + cell_node->sigma.y) / 2.0;
            *count_s1 = 1;
        }

        sigma_neighbours[0] = neighbours[RIGHT];
        sigma_neighbours[1] = neighbours[LEFT];

        sigma_neighbours[2] = neighbours[DOWN_RIGHT];
        sigma_neighbours[3] = neighbours[DOWN_LEFT];

        sigma_neighbours[4] = neighbours[FRONT_RIGHT];
        sigma_neighbours[5] = neighbours[FRONT_LEFT];

        sigma_neighbours[6] = neighbours[FRONT_DOWN_RIGHT];
        sigma_neighbours[7] = neighbours[FRONT_DOWN_LEFT];

        sigma_neighbours[8] = neighbours[BACK_RIGHT];
        sigma_neighbours[9] = neighbours[BACK_LEFT];

        sigma_neighbours[10] = neighbours[BACK_DOWN_RIGHT];
        sigma_neighbours[11] = neighbours[BACK_DOWN_LEFT];

        CALC_PARTIAL_SIGMA(xy, sigma_neighbours, *sigma_2, *count_s2);

        /////////////////////////////////////////////////////////////

        sigma_neighbours[0] = neighbours[FRONT];
        sigma_neighbours[1] = neighbours[BACK];

        sigma_neighbours[2] = neighbours[FRONT_DOWN];
        sigma_neighbours[3] = neighbours[BACK_DOWN];

        sigma_neighbours[4] = neighbours[FRONT_RIGHT];
        sigma_neighbours[5] = neighbours[BACK_RIGHT];

        sigma_neighbours[6] = neighbours[FRONT_DOWN_RIGHT];
        sigma_neighbours[7] = neighbours[BACK_DOWN_RIGHT];

        sigma_neighbours[8] = neighbours[FRONT_LEFT];
        sigma_neighbours[9] = neighbours[BACK_LEFT];

        sigma_neighbours[10] = neighbours[FRONT_DOWN_LEFT];
        sigma_neighbours[11] = neighbours[BACK_DOWN_LEFT];

        CALC_PARTIAL_SIGMA(yz, sigma_neighbours, *sigma_3, *count_s3);

    } else if(flux_direction == FRONT) {

        if(neighbours[FRONT]) {
            *sigma_1 = (neighbours[FRONT]->sigma.z + cell_node->sigma.z) / 2.0;
            *count_s1 = 1;
        }

        sigma_neighbours[0] = neighbours[RIGHT];
        sigma_neighbours[1] = neighbours[LEFT];

        sigma_neighbours[2] = neighbours[FRONT_RIGHT];
        sigma_neighbours[3] = neighbours[FRONT_LEFT];

        sigma_neighbours[4] = neighbours[TOP_RIGHT];
        sigma_neighbours[5] = neighbours[TOP_LEFT];

        sigma_neighbours[6] = neighbours[FRONT_TOP_RIGHT];
        sigma_neighbours[7] = neighbours[FRONT_TOP_LEFT];

        sigma_neighbours[8] = neighbours[DOWN_RIGHT];
        sigma_neighbours[9] = neighbours[DOWN_LEFT];

        sigma_neighbours[10] = neighbours[FRONT_DOWN_RIGHT];
        sigma_neighbours[11] = neighbours[FRONT_DOWN_LEFT];

        CALC_PARTIAL_SIGMA(xz, sigma_neighbours, *sigma_2, *count_s2);

        /////////////////////////////////////////////////////////////

        sigma_neighbours[0] = neighbours[TOP];
        sigma_neighbours[1] = neighbours[DOWN];

        sigma_neighbours[2] = neighbours[FRONT_TOP];
        sigma_neighbours[3] = neighbours[FRONT_DOWN];

        sigma_neighbours[4] = neighbours[TOP_RIGHT];
        sigma_neighbours[5] = neighbours[DOWN_RIGHT];

        sigma_neighbours[6] = neighbours[FRONT_TOP_RIGHT];
        sigma_neighbours[7] = neighbours[FRONT_DOWN_RIGHT];

        sigma_neighbours[8] = neighbours[TOP_LEFT];
        sigma_neighbours[9] = neighbours[DOWN_LEFT];

        sigma_neighbours[10] = neighbours[FRONT_TOP_LEFT];
        sigma_neighbours[11] = neighbours[FRONT_DOWN_LEFT];

        CALC_PARTIAL_SIGMA(yz, sigma_neighbours, *sigma_3, *count_s3);

    } else if(flux_direction == BACK) {
        if(neighbours[BACK]) {
            *sigma_1 = (neighbours[BACK]->sigma.z + cell_node->sigma.z) / 2.0;
            *count_s1 = 1;
        }

        sigma_neighbours[0] = neighbours[RIGHT];
        sigma_neighbours[1] = neighbours[LEFT];

        sigma_neighbours[2] = neighbours[BACK_RIGHT];
        sigma_neighbours[3] = neighbours[BACK_LEFT];

        sigma_neighbours[4] = neighbours[TOP_RIGHT];
        sigma_neighbours[5] = neighbours[TOP_LEFT];

        sigma_neighbours[6] = neighbours[BACK_TOP_RIGHT];
        sigma_neighbours[7] = neighbours[BACK_TOP_LEFT];

        sigma_neighbours[8] = neighbours[DOWN_RIGHT];
        sigma_neighbours[9] = neighbours[DOWN_LEFT];

        sigma_neighbours[10] = neighbours[BACK_DOWN_RIGHT];
        sigma_neighbours[11] = neighbours[BACK_DOWN_LEFT];

        CALC_PARTIAL_SIGMA(xz, sigma_neighbours, *sigma_2, *count_s2);

        /////////////////////////////////////////////////////////////

        sigma_neighbours[0] = neighbours[TOP];
        sigma_neighbours[1] = neighbours[DOWN];

        sigma_neighbours[2] = neighbours[BACK_TOP];
        sigma_neighbours[3] = neighbours[BACK_DOWN];

        sigma_neighbours[4] = neighbours[TOP_RIGHT];
        sigma_neighbours[5] = neighbours[DOWN_RIGHT];

        sigma_neighbours[6] = neighbours[BACK_TOP_RIGHT];
        sigma_neighbours[7] = neighbours[BACK_DOWN_RIGHT];

        sigma_neighbours[8] = neighbours[TOP_LEFT];
        sigma_neighbours[9] = neighbours[DOWN_LEFT];

        sigma_neighbours[10] = neighbours[BACK_TOP_LEFT];
        sigma_neighbours[11] = neighbours[BACK_DOWN_LEFT];

        CALC_PARTIAL_SIGMA(yz, sigma_neighbours, *sigma_3, *count_s3);
    }
}

static inline real_cpu DIVIDE(real_cpu num, real_cpu denom) {
    if(denom != 0) {
        return num / denom;
    }
    return 0.0;
}

#define UPDATE_OR_ADD_ELEMENT(g_cell, n_cell, v)                                                                                                               \
    do {                                                                                                                                                       \
        struct element el;                                                                                                                                     \
        int el_index = find_neighbour_index(g_cell, n_cell);                                                                                                   \
        if(el_index != -1) {                                                                                                                                   \
            g_cell->elements[el_index].value += v;                                                                                                             \
        } else {                                                                                                                                               \
            el.value = v;                                                                                                                                      \
            el.cell = n_cell;                                                                                                                                  \
            el.column = n_cell->grid_position;                                                                                                                 \
            arrput(g_cell->elements, el);                                                                                                                      \
        }                                                                                                                                                      \
    } while(0)

static void fill_elements_aniso(struct cell_node *grid_cell, struct cell_node *neighbours[26]) {

    struct element *elements = grid_cell->elements;

    real_cpu dx = grid_cell->discretization.x;
    real_cpu dy = grid_cell->discretization.y;
    real_cpu dz = grid_cell->discretization.z;

    real_cpu dx_squared = dx * dx;
    real_cpu dy_squared = dy * dy;
    real_cpu dz_squared = dz * dz;

    real_cpu dx_dy = dx_squared / dy;
    real_cpu dx_dz = dx_squared / dz;

    real_cpu dy_dx = dy_squared / dx;
    real_cpu dy_dz = dy_squared / dz;

    real_cpu dz_dx = dz_squared / dx;
    real_cpu dz_dy = dz_squared / dy;

    real_cpu sigma_x_r = 0.0;
    real_cpu sigma_xy_jx_r = 0.0;
    real_cpu sigma_xz_jx_r = 0.0;
    real_cpu sigma_x_l = 0.0;
    real_cpu sigma_xy_jx_l = 0.0;
    real_cpu sigma_xz_jx_l = 0.0;
    real_cpu sigma_xy_jy_t = 0.0;
    real_cpu sigma_y_t = 0.0;
    real_cpu sigma_yz_jy_t = 0.0;
    real_cpu sigma_xy_jy_d = 0.0;
    real_cpu sigma_y_d = 0.0;
    real_cpu sigma_yz_jy_d = 0.0;
    real_cpu sigma_xz_jz_f = 0.0;
    real_cpu sigma_yz_jz_f = 0.0;
    real_cpu sigma_z_f = 0.0;
    real_cpu sigma_xz_jz_b = 0.0;
    real_cpu sigma_yz_jz_b = 0.0;
    real_cpu sigma_z_b = 0.0;

    int count_sigma_x_r = 0;
    int count_sigma_xy_jx_r = 0;
    int count_sigma_xz_jx_r = 0;
    int count_sigma_x_l = 0;
    int count_sigma_xy_jx_l = 0;
    int count_sigma_xz_jx_l = 0;
    int count_sigma_xy_jy_t = 0;
    int count_sigma_y_t = 0;
    int count_sigma_yz_jy_t = 0;
    int count_sigma_xy_jy_d = 0;
    int count_sigma_y_d = 0;
    int count_sigma_yz_jy_d = 0;
    int count_sigma_xz_jz_f = 0;
    int count_sigma_yz_jz_f = 0;
    int count_sigma_z_f = 0;
    int count_sigma_xz_jz_b = 0;
    int count_sigma_yz_jz_b = 0;
    int count_sigma_z_b = 0;

    calc_sigmas(grid_cell, neighbours, RIGHT, &sigma_x_r, &sigma_xy_jx_r, &sigma_xz_jx_r, &count_sigma_x_r, &count_sigma_xy_jx_r, &count_sigma_xz_jx_r);
    calc_sigmas(grid_cell, neighbours, LEFT,  &sigma_x_l, &sigma_xy_jx_l, &sigma_xz_jx_l, &count_sigma_x_l, &count_sigma_xy_jx_l, &count_sigma_xz_jx_l);

    calc_sigmas(grid_cell, neighbours, TOP,  &sigma_y_t, &sigma_xy_jy_t, &sigma_yz_jy_t, &count_sigma_y_t, &count_sigma_xy_jy_t, &count_sigma_yz_jy_t);
    calc_sigmas(grid_cell, neighbours, DOWN, &sigma_y_d, &sigma_xy_jy_d, &sigma_yz_jy_d, &count_sigma_y_d, &count_sigma_xy_jy_d, &count_sigma_yz_jy_d);

    calc_sigmas(grid_cell, neighbours, FRONT, &sigma_z_f, &sigma_xz_jz_f, &sigma_yz_jz_f, &count_sigma_z_f, &count_sigma_xz_jz_f, &count_sigma_yz_jz_f);
    calc_sigmas(grid_cell, neighbours, BACK,  &sigma_z_b, &sigma_xz_jz_b, &sigma_yz_jz_b, &count_sigma_z_b, &count_sigma_xz_jz_b, &count_sigma_yz_jz_b);

    //MAIN DIAGONAL
    elements[0].value += dx*sigma_x_r + dx*sigma_x_l + dy*sigma_y_t + dy*sigma_y_d + dz*sigma_z_f + dz*sigma_z_b;

    real_cpu s1, s2, s3;

    // All neighbours
    for(int direction = 0; direction < NUM_DIRECTIONS; direction++) {

        if(neighbours[direction]) {

            struct element new_element;
            new_element.value = 0.0;

            switch(direction) {
            case FRONT:
                new_element.value += - sigma_z_f*dz;

                if(neighbours[BACK]) {
                    new_element.value += - DIVIDE(sigma_yz_jy_t,count_sigma_yz_jy_t) * dy_dz
                                         + DIVIDE(sigma_yz_jy_d,count_sigma_yz_jy_d) * dy_dz
                                         - DIVIDE(sigma_xz_jx_r,count_sigma_xz_jx_r) * dx_dz
                                         + DIVIDE(sigma_xz_jx_l,count_sigma_xz_jx_l) * dx_dz;
                }

                break;
            case BACK:
                new_element.value += - sigma_z_b*dz;

                if(neighbours[FRONT]) {
                    new_element.value +=   DIVIDE(sigma_yz_jy_t,count_sigma_yz_jy_t) * dy_dz
                                         - DIVIDE(sigma_yz_jy_d,count_sigma_yz_jy_d) * dy_dz
                                         + DIVIDE(sigma_xz_jx_r,count_sigma_xz_jx_r) * dx_dz
                                         - DIVIDE(sigma_xz_jx_l,count_sigma_xz_jx_l) * dx_dz;
                }
                break;

            case TOP:
                new_element.value += - sigma_y_t*dy;

                if(neighbours[DOWN]) {
                    new_element.value += - DIVIDE(sigma_yz_jz_f,count_sigma_yz_jz_f) * dz_dy
                                         + DIVIDE(sigma_yz_jz_b,count_sigma_yz_jz_b) * dz_dy
                                         - DIVIDE(sigma_xy_jx_r,count_sigma_xy_jx_r) * dx_dy
                                         + DIVIDE(sigma_xy_jx_l,count_sigma_xy_jx_l) * dx_dy;
                }
                break;
            case DOWN:
                new_element.value += - sigma_y_d*dy;
                if(neighbours[TOP]) {
                    new_element.value +=  DIVIDE(sigma_yz_jz_f,count_sigma_yz_jz_f) * dz_dy
                                        - DIVIDE(sigma_yz_jz_b,count_sigma_yz_jz_b) * dz_dy
                                        + DIVIDE(sigma_xy_jx_r,count_sigma_xy_jx_r) * dx_dy
                                        - DIVIDE(sigma_xy_jx_l,count_sigma_xy_jx_l) * dx_dy;
                }
                break;
            case RIGHT:
                new_element.value += - sigma_x_r*dx;
                if(neighbours[LEFT]) {
                    new_element.value += - DIVIDE(sigma_xz_jz_f,count_sigma_xz_jz_f) * dz_dx
                                         + DIVIDE(sigma_xz_jz_b,count_sigma_xz_jz_b) * dz_dx
                                         - DIVIDE(sigma_xy_jy_t,count_sigma_xy_jy_t) * dy_dx
                                         + DIVIDE(sigma_xy_jy_d,count_sigma_xy_jy_d) * dy_dx;
                }
                break;
            case LEFT:
                new_element.value += - sigma_x_l*dx;
                if(neighbours[RIGHT]) {
                    new_element.value += DIVIDE(sigma_xz_jz_f,count_sigma_xz_jz_f) * dz_dx
                                       - DIVIDE(sigma_xz_jz_b,count_sigma_xz_jz_b) * dz_dx
                                       + DIVIDE(sigma_xy_jy_t,count_sigma_xy_jy_t) * dy_dx
                                       - DIVIDE(sigma_xy_jy_d,count_sigma_xy_jy_d) * dy_dx;
                }
                break;

            case FRONT_TOP:

                s1 = 0.0;
                s2 = 0.0;

                if(neighbours[FRONT_DOWN]) {
                    s1 += - DIVIDE(sigma_xy_jx_r,count_sigma_xy_jx_r) * dx_dy
                          + DIVIDE(sigma_xy_jx_l,count_sigma_xy_jx_l) * dx_dy
                          - DIVIDE(sigma_yz_jz_f,count_sigma_yz_jz_f) * dz_dy;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[FRONT_DOWN], -s1);
                }

                if(neighbours[BACK_TOP]) {
                    s2 += - DIVIDE(sigma_yz_jy_t,count_sigma_yz_jy_t) * dy_dz
                          - DIVIDE(sigma_xz_jx_r,count_sigma_xz_jx_r) * dx_dz
                          + DIVIDE(sigma_xz_jx_l,count_sigma_xz_jx_l) * dx_dz;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_TOP], -s2);

                }

                new_element.value += s1 + s2;

                break;

            case FRONT_DOWN:

                s1 = 0.0;

                if(neighbours[BACK_DOWN]) {
                    s1  += DIVIDE(sigma_yz_jy_d,count_sigma_yz_jy_d) * dy_dz
                         - DIVIDE(sigma_xz_jx_r,count_sigma_xz_jx_r) * dx_dz
                         + DIVIDE(sigma_xz_jx_l,count_sigma_xz_jx_l) * dx_dz;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[FRONT_DOWN], s1);
                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_DOWN], -s1);
                }


                break;

            case BACK_TOP:

                s1 = 0.0;

                if(neighbours[BACK_DOWN]) {

                    s1 +=   DIVIDE(sigma_yz_jz_b,count_sigma_yz_jz_b) * dz_dy
                          - DIVIDE(sigma_xy_jx_r,count_sigma_xy_jx_r) * dx_dy
                          + DIVIDE(sigma_xy_jx_l,count_sigma_xy_jx_l) * dx_dy;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_TOP] ,  s1);
                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_DOWN], -s1);

                }

                break;
            case BACK_DOWN: break; //alread handled by the above cases

            case FRONT_RIGHT:
                s1 = 0.0;
                s2 = 0.0;

                if(neighbours[BACK_RIGHT]) {
                    s1 += - DIVIDE(sigma_xz_jx_r,count_sigma_xz_jx_r) * dx_dz
                          - DIVIDE(sigma_yz_jy_t,count_sigma_yz_jy_t) * dy_dz
                          + DIVIDE(sigma_yz_jy_d,count_sigma_yz_jy_d) * dy_dz;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_RIGHT] ,  -s1);
                }

                if(neighbours[FRONT_LEFT]) {
                    s2 += - DIVIDE(sigma_xz_jz_f,count_sigma_xz_jz_f) * dz_dx
                          - DIVIDE(sigma_xy_jy_t,count_sigma_xy_jy_t) * dy_dx
                          + DIVIDE(sigma_xy_jy_d,count_sigma_xy_jy_d) * dy_dx;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[FRONT_LEFT] ,  -s2);
                }

                new_element.value += s1 + s2;
                break;

            case FRONT_LEFT:

                s1 = 0.0;

                if(neighbours[BACK_LEFT]) {

                    s1 += - DIVIDE(sigma_yz_jy_t,count_sigma_yz_jy_t) * dy_dz
                          + DIVIDE(sigma_yz_jy_d,count_sigma_yz_jy_d) * dy_dz
                          + DIVIDE(sigma_xz_jx_l,count_sigma_xz_jx_l) * dx_dz;


                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[FRONT_LEFT] ,  s1);
                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_LEFT] ,  -s1);

                }
                break;
            case BACK_RIGHT:
                s1 = 0.0;
                if(neighbours[BACK_LEFT]) {
                    s1 +=   DIVIDE(sigma_xz_jz_b,count_sigma_xz_jz_b) * dz_dx
                          - DIVIDE(sigma_xy_jy_t,count_sigma_xy_jy_t) * dy_dx
                          + DIVIDE(sigma_xy_jy_d,count_sigma_xy_jy_d) * dy_dx;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_RIGHT] ,  s1);
                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_LEFT] ,  -s1);
                }

                break;

            case BACK_LEFT: break;  //alread handled by the above cases

            case TOP_RIGHT:

                s1 = 0;
                s2 = 0;
                if(neighbours[DOWN_RIGHT]) {
                    s1 += - DIVIDE(sigma_xy_jx_r,count_sigma_xy_jx_r) * dx_dy
                          - DIVIDE(sigma_yz_jz_f,count_sigma_yz_jz_f) * dz_dy
                          + DIVIDE(sigma_yz_jz_b,count_sigma_yz_jz_b) * dz_dy;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[DOWN_RIGHT], -s1);
                }

                if(neighbours[TOP_LEFT]) {
                    s2 += - DIVIDE(sigma_xz_jz_f,count_sigma_xz_jz_f) * dz_dx
                          + DIVIDE(sigma_xz_jz_b,count_sigma_xz_jz_b) * dz_dx
                          - DIVIDE(sigma_xy_jy_t,count_sigma_xy_jy_t) * dy_dx;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[TOP_LEFT], -s2);

                }

                new_element.value += s1 + s2;

                break;
            case TOP_LEFT:

                s1 = 0;
                if(neighbours[DOWN_LEFT]) {
                    s1 += - DIVIDE(sigma_yz_jz_f,count_sigma_yz_jz_f) * dz_dy
                          + DIVIDE(sigma_yz_jz_b,count_sigma_yz_jz_b) * dz_dy
                          + DIVIDE(sigma_xy_jx_l,count_sigma_xy_jx_l) * dx_dy;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[TOP_LEFT] ,  s1);
                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[DOWN_LEFT] ,  -s1);

                }
                break;
            case DOWN_RIGHT:
                s1 = 0;
                if(neighbours[DOWN_LEFT]) {
                    s1 += - DIVIDE(sigma_xz_jz_f,count_sigma_xz_jz_f) * dz_dx
                          + DIVIDE(sigma_xz_jz_b,count_sigma_xz_jz_b) * dz_dx
                          + DIVIDE(sigma_xy_jy_d,count_sigma_xy_jy_d) * dy_dx;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[DOWN_RIGHT] ,  s1);
                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[DOWN_LEFT] ,  -s1);


                }
                break;
            case DOWN_LEFT:  break;  //alread handled by the above cases

           case FRONT_TOP_RIGHT:
                s1 = 0;
                s2 = 0;
                s3 = 0;
                if(neighbours[FRONT_DOWN_RIGHT]) {
                    s1 += - DIVIDE(sigma_xy_jx_r,count_sigma_xy_jx_r) * dx_dy
                          - DIVIDE(sigma_yz_jz_f,count_sigma_yz_jz_f) * dz_dy;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[FRONT_DOWN_RIGHT], -s1);
                }
                if(neighbours[BACK_TOP_RIGHT]) {
                    s2 += - DIVIDE(sigma_xz_jx_r,count_sigma_xz_jx_r) * dx_dz
                          - DIVIDE(sigma_yz_jy_t,count_sigma_yz_jy_t) * dy_dz;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_TOP_RIGHT], -s2);
                }
                if(neighbours[FRONT_TOP_LEFT]) {
                    s3 += - DIVIDE(sigma_xz_jz_f,count_sigma_xz_jz_f) * dz_dx
                          - DIVIDE(sigma_xy_jy_t,count_sigma_xy_jy_t) * dy_dx;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[FRONT_TOP_LEFT], -s3);
                }
                new_element.value += s1 + s2 + s3;
                break;

           case FRONT_TOP_LEFT:

                s1 = 0;
                s2 = 0;

                if(neighbours[FRONT_DOWN_LEFT]) {
                    s1 +=   DIVIDE(sigma_xy_jx_l,count_sigma_xy_jx_l) * dx_dy
                          - DIVIDE(sigma_yz_jz_f,count_sigma_yz_jz_f) * dz_dy;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[FRONT_TOP_LEFT], s1);
                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[FRONT_DOWN_LEFT], -s1);
                }

                if(neighbours[BACK_TOP_LEFT]) {
                    s2 +=  DIVIDE(sigma_xz_jx_l,count_sigma_xz_jx_l) * dx_dz
                         - DIVIDE(sigma_yz_jy_t,count_sigma_yz_jy_t) * dy_dz;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[FRONT_TOP_LEFT], s2);
                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_TOP_LEFT], -s2);
                }

                break;
           case FRONT_DOWN_RIGHT:

                s1 = 0;
                s2 = 0;

                if(neighbours[BACK_DOWN_RIGHT]) {
                    s1 += - DIVIDE(sigma_xz_jx_r,count_sigma_xz_jx_r) * dx_dz
                          + DIVIDE(sigma_yz_jy_d,count_sigma_yz_jy_d) * dy_dz;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[FRONT_DOWN_RIGHT], s1);
                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_DOWN_RIGHT], -s1);
                }

                if(neighbours[FRONT_DOWN_LEFT]) {
                    s2 +=  DIVIDE(sigma_xy_jy_d,count_sigma_xy_jy_d) * dy_dx
                         - DIVIDE(sigma_xz_jz_f,count_sigma_xz_jz_f) * dz_dx;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[FRONT_DOWN_RIGHT], s2);
                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[FRONT_DOWN_LEFT], -s2);
                }

                break;
           case FRONT_DOWN_LEFT:
                s1 = 0;
                if(neighbours[BACK_DOWN_LEFT]) {
                    s1 +=   DIVIDE(sigma_xz_jx_l,count_sigma_xz_jx_l) * dx_dz
                          + DIVIDE(sigma_yz_jy_d,count_sigma_yz_jy_d) * dy_dz;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[FRONT_DOWN_LEFT], s1);
                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_DOWN_LEFT], -s1);

                }
                break;
           case BACK_TOP_RIGHT:
                s1 = 0;
                s2 = 0;

                if(neighbours[BACK_DOWN_RIGHT]) {
                    s1 += - DIVIDE(sigma_xy_jx_r,count_sigma_xy_jx_r) * dx_dy
                          + DIVIDE(sigma_yz_jz_b,count_sigma_yz_jz_b) * dz_dy;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_TOP_RIGHT], s1);
                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_DOWN_RIGHT], -s1);

                }

                if(neighbours[BACK_TOP_LEFT]) {
                    s2 += - DIVIDE(sigma_xy_jy_t,count_sigma_xy_jy_t) * dy_dx
                          + DIVIDE(sigma_xz_jz_b,count_sigma_xz_jz_b) * dz_dx;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_TOP_RIGHT], s2);
                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_TOP_LEFT], -s2);

                }
                break;
           case BACK_TOP_LEFT:
                s1 = 0;
                if(neighbours[BACK_DOWN_LEFT]) {
                    s1 +=   DIVIDE(sigma_xy_jx_l,count_sigma_xy_jx_l) * dx_dy
                          + DIVIDE(sigma_yz_jz_b,count_sigma_yz_jz_b) * dz_dy;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_TOP_LEFT], s1);
                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_DOWN_LEFT], -s1);

                }

                break;
            case BACK_DOWN_RIGHT:
                s1 = 0;
                if(neighbours[BACK_DOWN_LEFT]) {
                    s1 +=  DIVIDE(sigma_xy_jy_d,count_sigma_xy_jy_d) * dy_dx
                         + DIVIDE(sigma_xz_jz_b,count_sigma_xz_jz_b) * dz_dx;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_DOWN_RIGHT], s1);
                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_DOWN_LEFT], -s1);

                }
               break;

            case BACK_DOWN_LEFT: break;

            default:
                break;
            }

            if(new_element.value != 0.0) {
                new_element.column = neighbours[direction]->grid_position;
                new_element.cell = neighbours[direction];
                arrput(grid_cell->elements, new_element);
            }
        }
    }
}

static void debug_cell(struct cell_node *grid_cell, struct cell_node *neighbours[26]) {

    printf("%lf, %lf, %lf - %d\n", grid_cell->center.x, grid_cell->center.y, grid_cell->center.z, grid_cell->grid_position + 1);
    printf("Sigmas %lf, %lf, %lf\n", grid_cell->sigma.x, grid_cell->sigma.y, grid_cell->sigma.z);
    printf("Main diag %.20lf\n", grid_cell->elements[0].value);

    if(neighbours[TOP])
        printf("TOP %d\n", neighbours[TOP]->grid_position + 1);
    if(neighbours[DOWN])
        printf("DOWN %d\n", neighbours[DOWN]->grid_position + 1);
    if(neighbours[FRONT])
        printf("FRONT %d\n", neighbours[FRONT]->grid_position + 1);
    if(neighbours[BACK])
        printf("BACK %d\n", neighbours[BACK]->grid_position + 1);
    if(neighbours[LEFT])
        printf("LEFT %d\n", neighbours[LEFT]->grid_position + 1);
    if(neighbours[RIGHT])
        printf("RIGHT %d\n", neighbours[RIGHT]->grid_position + 1);
    if(neighbours[FRONT_RIGHT])
        printf("FRONT_RIGHT %d\n", neighbours[FRONT_RIGHT]->grid_position + 1);
    if(neighbours[FRONT_LEFT])
        printf("FRONT_LEFT %d\n", neighbours[FRONT_LEFT]->grid_position + 1);
    if(neighbours[BACK_RIGHT])
        printf("BACK_RIGHT %d\n", neighbours[BACK_RIGHT]->grid_position + 1);
    if(neighbours[BACK_LEFT])
        printf("BACK_LEFT %d\n", neighbours[BACK_LEFT]->grid_position + 1);
    if(neighbours[TOP_RIGHT])
        printf("TOP_RIGHT %d\n", neighbours[TOP_RIGHT]->grid_position + 1);
    if(neighbours[TOP_LEFT])
        printf("TOP_LEFT %d\n", neighbours[TOP_LEFT]->grid_position + 1);
    if(neighbours[DOWN_RIGHT])
        printf("DOWN_RIGHT %d\n", neighbours[DOWN_RIGHT]->grid_position + 1);
    if(neighbours[DOWN_LEFT])
        printf("DOWN_LEFT %d\n", neighbours[DOWN_LEFT]->grid_position + 1);
    if(neighbours[FRONT_TOP])
        printf("FRONT_TOP %d\n", neighbours[FRONT_TOP]->grid_position + 1);
    if(neighbours[FRONT_DOWN])
        printf("FRONT_DOWN %d\n", neighbours[FRONT_DOWN]->grid_position + 1);
    if(neighbours[FRONT_TOP_RIGHT])
        printf("FRONT_TOP_RIGHT %d\n", neighbours[FRONT_TOP_RIGHT]->grid_position + 1);
    if(neighbours[FRONT_TOP_LEFT])
        printf("FRONT_TOP_LEFT %d\n", neighbours[FRONT_TOP_LEFT]->grid_position + 1);
    if(neighbours[FRONT_DOWN_RIGHT])
        printf("FRONT_DOWN_RIGHT %d\n", neighbours[FRONT_DOWN_RIGHT]->grid_position + 1);
    if(neighbours[FRONT_DOWN_LEFT])
        printf("FRONT_DOWN_LEFT %d\n", neighbours[FRONT_DOWN_LEFT]->grid_position + 1);
    if(neighbours[BACK_TOP])
        printf("BACK_TOP %d\n", neighbours[BACK_TOP]->grid_position + 1);
    if(neighbours[BACK_DOWN])
        printf("BACK_DOWN %d\n", neighbours[BACK_DOWN]->grid_position + 1);
    if(neighbours[BACK_TOP_RIGHT])
        printf("BACK_TOP_RIGHT %d\n", neighbours[BACK_TOP_RIGHT]->grid_position + 1);
    if(neighbours[BACK_TOP_LEFT])
        printf("BACK_TOP_LEFT %d\n", neighbours[BACK_TOP_LEFT]->grid_position + 1);
    if(neighbours[BACK_DOWN_RIGHT])
        printf("BACK_DOWN_RIGHT %d\n", neighbours[BACK_DOWN_RIGHT]->grid_position + 1);
    if(neighbours[BACK_DOWN_LEFT])
        printf("BACK_DOWN_LEFT %d\n", neighbours[BACK_DOWN_LEFT]->grid_position + 1);
}

static void fill_discretization_matrix_elements_aniso(struct cell_node *grid_cell) {

    struct cell_node *neighbours[26];
    struct cell_node *neighbour;
    int n = 0;
    for(int direction = 0; direction < NUM_DIRECTIONS; direction++) {
        neighbour = get_cell_neighbour_with_same_refinement_level(grid_cell, direction);
        if(neighbour && neighbour->active) {
            neighbours[direction] = neighbour;
            n++;

        } else {
            neighbours[direction] = NULL;
        }
    }

    fill_elements_aniso(grid_cell, neighbours);

}

void outer_product_vector_vector_t(real_cpu p[3][3], real_cpu v[3]) {
    p[0][0] = v[0] * v[0];
    p[0][1] = v[0] * v[1];
    p[0][2] = v[0] * v[2];

    p[1][0] = v[1] * v[0];
    p[1][1] = v[1] * v[1];
    p[1][2] = v[1] * v[2];

    p[2][0] = v[2] * v[0];
    p[2][1] = v[2] * v[1];
    p[2][2] = v[2] * v[2];
}

static inline void scalar_tensor(real_cpu s, real_cpu t[3][3]) {
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            t[i][j] = t[i][j] * s;
        }
    }
}

static inline void sum_tensor(real_cpu tr[3][3], real_cpu t1[3][3], real_cpu t2[3][3]) {
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            tr[i][j] = t1[i][j] + t2[i][j];
        }
    }
}

static inline void print_tensor(real_cpu t[3][3]) {
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            printf("%e ", t[i][j]);
        }
        printf("\n");
    }
}

// D(x) = sigma_l f kp f + sigma_t s kp s + sigma_n n kp n
static void calc_tensor(real_cpu D[3][3], real_cpu f[3], real_cpu s[3], real_cpu n[3], real_cpu sigma_l, real_cpu sigma_t, real_cpu sigma_n) {

    real_cpu tmp[3][3];
    real_cpu fft[3][3];
    real_cpu sst[3][3];
    real_cpu nnt[3][3];

    outer_product_vector_vector_t(fft, f);
    outer_product_vector_vector_t(sst, s);
    outer_product_vector_vector_t(nnt, n);

    scalar_tensor(sigma_l, fft);
    scalar_tensor(sigma_t, sst);
    scalar_tensor(sigma_n, nnt);

    sum_tensor(tmp, fft, sst);
    sum_tensor(D, tmp, nnt);
}

static void calc_tensor2(real_cpu D[3][3], real_cpu f[3], real_cpu sigma_l, real_cpu sigma_t) {
    // D = ((sigma_L - sigma_T) * outer_product(F, FT)) + sigma_T * ident(3);
    real_cpu fft[3][3];
    real_cpu sigma_ident[3][3];

    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            if(i == j) {
                sigma_ident[i][j] = sigma_t;
            } else {
                sigma_ident[i][j] = 0.0;
            }
        }
    }

    outer_product_vector_vector_t(fft, f);
    scalar_tensor(sigma_l - sigma_t, fft);
    sum_tensor(D, fft, sigma_ident);
}

static inline void normalize(real_cpu v[3]) {
    real_cpu m = sqrt((v[0] * v[0]) + (v[1] * v[1]) + (v[2] * v[2]));
    v[0] = v[0] / m;
    v[1] = v[1] / m;
    v[2] = v[2] / m;
}

static struct fiber_coords *read_fibers(char *fiber_file_path, bool normalize_vector) {

    FILE *fibers_file = open_file_or_exit(fiber_file_path, "r");

    struct fiber_coords *fibers = NULL;
    char *line = NULL;
    size_t len;

    while((getline(&line, &len, fibers_file)) != -1) {

        int split_count;
        sds *points = sdssplit(line, " ", &split_count);
        struct fiber_coords f_coords;

        f_coords.f[0] = strtod(points[0], NULL);
        f_coords.f[1] = strtod(points[1], NULL);
        f_coords.f[2] = strtod(points[2], NULL);

        f_coords.s[0] = strtod(points[3], NULL);
        f_coords.s[1] = strtod(points[4], NULL);
        f_coords.s[2] = strtod(points[5], NULL);

        f_coords.n[0] = strtod(points[6], NULL);
        f_coords.n[1] = strtod(points[7], NULL);
        f_coords.n[2] = strtod(points[8], NULL);

        if(normalize_vector) {
            normalize(f_coords.f);
            normalize(f_coords.s);
            normalize(f_coords.n);
        }

        arrput(fibers, f_coords);
        sdsfreesplitres(points, split_count);
    }

    free(line);

    return fibers;
}

int randRange(int n) {
    int limit;
    int r;

    limit = RAND_MAX - (RAND_MAX % n);

    while((r = rand()) >= limit)
        ;

    return r % n;
}

void initialize_diagonal_elements_purkinje (struct monodomain_solver *the_solver, struct grid *the_grid) {

    real_cpu alpha;
    real_cpu dx, dy, dz;

    uint32_t num_active_cells = the_grid->purkinje->num_active_purkinje_cells;
    struct cell_node **ac = the_grid->purkinje->purkinje_cells;

    struct node *n = the_grid->purkinje->network->list_nodes;

    real_cpu beta = the_solver->beta;
    real_cpu cm = the_solver->cm;
    real_cpu dt = the_solver->dt;

    uint32_t i;

    for (i = 0; i < num_active_cells; i++) {

        dx = ac[i]->discretization.x;
        dy = ac[i]->discretization.y;
        dz = ac[i]->discretization.z;

        alpha = ALPHA(beta, cm, dt, dx, dy, dz);

        struct element element;
        element.column = ac[i]->grid_position;
        element.cell = ac[i];
        element.value = alpha;

        if (ac[i]->elements != NULL)
            arrfree(ac[i]->elements);

        ac[i]->elements = NULL;
        arrsetcap(ac[i]->elements,n->num_edges);
        arrput(ac[i]->elements, element);

        n = n->next;
    }
}

// For the Purkinje fibers we only need to solve the 1D Monodomain equation
static void fill_discretization_matrix_elements_purkinje (bool has_point_data, real_cpu sigma_x, struct cell_node **grid_cells, uint32_t num_active_cells,
                                                        struct node *pk_node) {

    struct edge *e;
    struct element **cell_elements;

    real_cpu dx, dy, dz;
    real_cpu multiplier;

    int i;

    for (i = 0; i < num_active_cells; i++, pk_node = pk_node->next) {

        cell_elements = &grid_cells[i]->elements;
        dx = grid_cells[i]->discretization.x;
        dy = grid_cells[i]->discretization.x;
        dz = grid_cells[i]->discretization.x;

        multiplier = ((dy * dz) / dx);

        e = pk_node->list_edges;

        // Do the mapping of the edges from the graph to the sparse matrix data structure ...
        while (e != NULL) {

            struct element new_element;

            // Calculate the conductivity between the two neighboring cells
            if (has_point_data) {
                real_cpu sigma_x1 = pk_node->sigma;
                real_cpu sigma_x2 = e->dest->sigma;
                
                if(sigma_x1 != 0.0 && sigma_x2 != 0.0) 
                    sigma_x = (2.0f * sigma_x1 * sigma_x2) / (sigma_x1 + sigma_x2);
            }

            // Neighbour elements ...
            new_element.column = e->id;
            new_element.value = (-sigma_x * multiplier);
            new_element.cell = grid_cells[e->id];

            // Diagonal element ...
            cell_elements[0]->value += (sigma_x * multiplier);

            arrput(grid_cells[i]->elements,new_element);

            e = e->next;
        }
    }
}

ASSEMBLY_MATRIX (purkinje_coupling_assembly_matrix) {

// [TISSUE]
    static bool sigma_initialized = false;

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    uint32_t i;

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config, "sigma_z");

    real sigma_purkinje = sigma_x;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real,sigma_purkinje,config,"sigma_purkinje");

    if(!sigma_initialized) {

        OMP(parallel for)
        for (i = 0; i < num_active_cells; i++) {
            ac[i]->sigma.x = sigma_x;
            ac[i]->sigma.y = sigma_y;
            ac[i]->sigma.z = sigma_z;
        }

        sigma_initialized = true;
    }

    OMP(parallel for)
    for(i = 0; i < num_active_cells; i++) {

        // Computes and designates the flux due to south cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[BACK], BACK);

        // Computes and designates the flux due to north cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[FRONT], FRONT);

        // Computes and designates the flux due to east cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[TOP], TOP);

        // Computes and designates the flux due to west cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[DOWN], DOWN);

        // Computes and designates the flux due to front cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[RIGHT], RIGHT);

        // Computes and designates the flux due to back cells.
        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[LEFT], LEFT);
    }

// [PURKINJE]
    static bool sigma_purkinje_initialized = false;

    uint32_t num_purkinje_active_cells = the_grid->purkinje->num_active_purkinje_cells;
    struct cell_node **ac_purkinje = the_grid->purkinje->purkinje_cells;
    struct node *pk_node = the_grid->purkinje->network->list_nodes;
    bool has_point_data = the_grid->purkinje->network->has_point_data;

    initialize_diagonal_elements_purkinje(the_solver, the_grid);

    if(!sigma_purkinje_initialized) {
        // Check if the Purkinje network file has the POINT_DATA section
        if (has_point_data) {
            struct node *tmp = the_grid->purkinje->network->list_nodes;
            uint32_t i = 0;
            while (tmp != NULL)
            {
                // Copy the prescribed conductivity from the Purkinje network file into the ALG Purkinje cell structure
                ac_purkinje[i]->sigma.x = tmp->sigma;

                tmp = tmp->next; i++;
            }
        } 
        // Otherwise, initilize the conductivity of all cells homogenously with the value from the configuration file
        else {
            OMP(parallel for)
            for (uint32_t i = 0; i < num_active_cells; i++) {
                ac[i]->sigma.x = sigma_purkinje;
            }
        }
        sigma_purkinje_initialized = true;
    }
    fill_discretization_matrix_elements_purkinje(has_point_data,sigma_purkinje,ac_purkinje,num_purkinje_active_cells,pk_node);
}

ASSEMBLY_MATRIX (purkinje_coupling_with_anisotropic_sigma_assembly_matrix) {

// [TISSUE]
    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    //      D tensor    //
    // | sx    sxy   sxz |
    // | sxy   sy    syz |
    // | sxz   syz   sz  |
    real_cpu D[3][3];
    int i;

    real_cpu sigma_l = 0.0;
    real_cpu sigma_t = 0.0;
    real_cpu sigma_n = 0.0;
    real_cpu sigma_purkinje = 0.0;

    char *fiber_file = NULL;
    GET_PARAMETER_STRING_VALUE_OR_USE_DEFAULT(fiber_file, config, "fibers_file");

    bool fibers_in_mesh = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(fibers_in_mesh, config, "fibers_in_mesh");


    struct fiber_coords *fibers = NULL;

    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_l, config, "sigma_l");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_t, config, "sigma_t");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_n, config, "sigma_n");
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real,sigma_purkinje,config,"sigma_purkinje");

    real_cpu *f = NULL;
    real_cpu *s = NULL;
    real_cpu *n = NULL;

    if(fiber_file) {
        log_info("Loading mesh fibers\n");
        fibers = read_fibers(fiber_file, false);
    }
    else if(!fibers_in_mesh) {
        GET_PARAMETER_VECTOR_VALUE_OR_USE_DEFAULT(f, config, "f", 3);
        GET_PARAMETER_VECTOR_VALUE_OR_USE_DEFAULT(s, config, "s", 3);
        GET_PARAMETER_VECTOR_VALUE_OR_USE_DEFAULT(n, config, "n", 3);

        if(!f) {
            f = malloc(sizeof(real_cpu)*3);
            f[0] = 1.0;
            f[1] = 0.0;
            f[2] = 0.0;
        }

        if(!s) {
            s = malloc(sizeof(real_cpu)*3);
            s[0] = 0.0;
            s[1] = 1.0;
            s[2] = 0.0;
        }

        if(!n) {
            n = malloc(sizeof(real_cpu)*3);
            n[0] = 0.0;
            n[1] = 0.0;
            n[2] = 1.0;
        }

    }

    OMP(parallel for private(D))
    for(i = 0; i < num_active_cells; i++) {

        if(fibers) {
            int fiber_index = ac[i]->original_position_in_file;

            if(fiber_index == -1) {
                log_error_and_exit("fiber_index should not be -1, but it is for cell in index %d - %lf, %lf, %lf\n", i, ac[i]->center.x, ac[i]->center.y, ac[i]->center.z);
            }

            if(sigma_t == sigma_n) {
                calc_tensor2(D, fibers[fiber_index].f, sigma_l, sigma_t);
            }
            else {
                calc_tensor(D, fibers[fiber_index].f, fibers[fiber_index].s, fibers[fiber_index].n, sigma_l, sigma_t, sigma_n);
            }
            ac[i]->sigma.fibers = fibers[fiber_index];
        }
        else if(fibers_in_mesh) {
            if(sigma_t == sigma_n) {
                calc_tensor2(D, ac[i]->sigma.fibers.f, sigma_l, sigma_t);
            }
            else {
                calc_tensor(D, ac[i]->sigma.fibers.f, ac[i]->sigma.fibers.s, ac[i]->sigma.fibers.n, sigma_l, sigma_t, sigma_n);
            }

        }
        else {
            if(sigma_t == sigma_n) {
                calc_tensor2(D, f, sigma_l, sigma_t);
            }
            else {
                calc_tensor(D, f, s, n, sigma_l, sigma_t, sigma_n);
            }
        }

        ac[i]->sigma.x = D[0][0];
        ac[i]->sigma.y = D[1][1];
        ac[i]->sigma.z = D[2][2];

        ac[i]->sigma.xy = D[0][1];
        ac[i]->sigma.xz = D[0][2];
        ac[i]->sigma.yz = D[1][2];

    }

    OMP(parallel for)
    for(i = 0; i < num_active_cells; i++) {
        fill_discretization_matrix_elements_aniso(ac[i]);
    }

    free(f);
    free(s);
    free(n);

// [PURKINJE]
    static bool sigma_purkinje_initialized = false;

    uint32_t num_purkinje_active_cells = the_grid->purkinje->num_active_purkinje_cells;
    struct cell_node **ac_purkinje = the_grid->purkinje->purkinje_cells;
    struct node *pk_node = the_grid->purkinje->network->list_nodes;
    bool has_point_data = the_grid->purkinje->network->has_point_data;

    initialize_diagonal_elements_purkinje(the_solver, the_grid);

    if(!sigma_purkinje_initialized) {
        // Check if the Purkinje network file has the POINT_DATA section
        if (has_point_data) {
            struct node *tmp = the_grid->purkinje->network->list_nodes;
            uint32_t i = 0;
            while (tmp != NULL)
            {
                // Copy the prescribed conductivity from the Purkinje network file into the ALG Purkinje cell structure
                ac_purkinje[i]->sigma.x = tmp->sigma;

                tmp = tmp->next; i++;
            }
        } 
        // Otherwise, initilize the conductivity of all cells homogenously with the value from the configuration file
        else {
            OMP(parallel for)
            for (uint32_t i = 0; i < num_active_cells; i++) {
                ac[i]->sigma.x = sigma_purkinje;
            }
        }
        sigma_purkinje_initialized = true;
    }
    fill_discretization_matrix_elements_purkinje(has_point_data,sigma_purkinje,ac_purkinje,num_purkinje_active_cells,pk_node);
}
