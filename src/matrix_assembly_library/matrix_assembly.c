//
// Created by sachetto on 13/10/17.
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

INIT_ASSEMBLY_MATRIX(set_initial_conditions_fvm) {

    real_cpu alpha;

    real_cpu beta = the_solver->beta;
    real_cpu cm = the_solver->cm;
    real_cpu dt = the_solver->dt;

    struct cell_node **ac = the_grid->active_cells;
    uint32_t active_cells = the_grid->num_active_cells;

    OMP(parallel for private(alpha))
    for(uint32_t i = 0; i < active_cells; i++) {
        alpha = ALPHA(beta, cm, dt, ac[i]->discretization.x, ac[i]->discretization.y, ac[i]->discretization.z);
        ac[i]->v = initial_v;
        ac[i]->b = initial_v * alpha;
    }
}

static struct element fill_element(uint32_t position, enum transition_direction direction, real_cpu dx, real_cpu dy, real_cpu dz, real_cpu sigma_x,
                                   real_cpu sigma_y, real_cpu sigma_z, struct element *cell_elements);

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
        // alpha = 0.0;

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

struct element fill_element(uint32_t position, enum transition_direction direction, real_cpu dx, real_cpu dy, real_cpu dz, real_cpu sigma_x, real_cpu sigma_y,
                            real_cpu sigma_z, struct element *cell_elements) {

    real_cpu multiplier;

    struct element new_element;
    new_element.column = position;

    if(direction == FRONT) { // Z direction front
        multiplier = ((dx * dy) / dz);
        new_element.value = -sigma_z * multiplier;
        cell_elements[0].value += (sigma_z * multiplier);
    } else if(direction == BACK) { // Z direction back
        multiplier = ((dx * dy) / dz);
        new_element.value = -sigma_z * multiplier;
        cell_elements[0].value += (sigma_z * multiplier);
    } else if(direction == TOP) { // Y direction top
        multiplier = ((dx * dz) / dy);
        new_element.value = -sigma_y * multiplier;
        cell_elements[0].value += (sigma_y * multiplier);
    } else if(direction == DOWN) { // Y direction down
        multiplier = ((dx * dz) / dy);
        new_element.value = -sigma_y * multiplier;
        cell_elements[0].value += (sigma_y * multiplier);
    } else if(direction == RIGHT) { // X direction right
        multiplier = ((dy * dz) / dx);
        new_element.value = -sigma_x * multiplier;
        cell_elements[0].value += (sigma_x * multiplier);
    } else if(direction == LEFT) { // X direction left
        multiplier = ((dy * dz) / dx);
        new_element.value = -sigma_x * multiplier;
        cell_elements[0].value += (sigma_x * multiplier);
    }
    return new_element;
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

static int rand_range(int n) {
    int limit;
    int r;

    limit = RAND_MAX - (RAND_MAX % n);

    while((r = rand()) >= limit)
        ;

    return r % n;
}

ASSEMBLY_MATRIX(random_sigma_discretization_matrix) {

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    int i;

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config->config_data, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config->config_data, "sigma_z");

    srand((unsigned int)time(NULL));

    real_cpu modifiers[4] = {0.0f, 0.1f, 0.5f, 1.0f};

    OMP(parallel for)
    for(i = 0; i < num_active_cells; i++) {
        real_cpu r;

        OMP(critical)
        r = modifiers[rand_range(4)];

        real sigma_x_new = sigma_x * r;
        real sigma_y_new = sigma_y * r;
        real sigma_z_new = sigma_z * r;

        ac[i]->sigma.x = sigma_x_new;
        ac[i]->sigma.y = sigma_y_new;
        ac[i]->sigma.z = sigma_z_new;
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

        fill_discretization_matrix_elements(ac[i], ac[i]->neighbours[LEFT], LEFT);
    }
}

ASSEMBLY_MATRIX(source_sink_discretization_matrix) {

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    uint32_t i;

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config->config_data, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config->config_data, "sigma_z");

    real channel_width = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, channel_width, config->config_data, "channel_width");

    real channel_length = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, channel_length, config->config_data, "channel_length");

    bool inside;

    // real side_length_x = the_grid->mesh_side_length.x;
    real side_length_y = the_grid->mesh_side_length.y;
    // real side_length_z = the_grid->mesh_side_length.z;

    real region_height = (side_length_y - channel_width) / 2.0;

    OMP(parallel for)
    for(i = 0; i < num_active_cells; i++) {
        real sigma_x_new;
        real sigma_y_new;
        real sigma_z_new;

        double x = ac[i]->center.x;
        double y = ac[i]->center.y;
        //        double z = ac[i]->center.z;

        // Check region 1
        inside = (x >= 0.0) && (x <= channel_length) && (y >= 0.0) && (y <= region_height);

        // Check region 2
        inside |= (x >= 0.0) && (x <= channel_length) && (y >= region_height + channel_width) && (y <= side_length_y);

        if(inside) {
            sigma_x_new = 0.0;
            sigma_y_new = 0.0;
            sigma_z_new = 0.0;
        } else {
            sigma_x_new = sigma_x;
            sigma_y_new = sigma_y;
            sigma_z_new = sigma_z;
        }

        ac[i]->sigma.x = sigma_x_new;
        ac[i]->sigma.y = sigma_y_new;
        ac[i]->sigma.z = sigma_z_new;
    }

    // Then, we fill the discretization matrix
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
}

ASSEMBLY_MATRIX(source_sink_discretization_matrix_with_different_sigma) {

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    int i;

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config->config_data, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config->config_data, "sigma_z");

    real channel_width = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, channel_width, config->config_data, "channel_width");

    real channel_length = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, channel_length, config->config_data, "channel_length");

    real source_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, source_factor, config->config_data, "source_factor");

    real sink_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sink_factor, config->config_data, "sink_factor");

    bool inside_3, inside_4;

    real side_length_x = the_grid->mesh_side_length.x;
    real side_length_y = the_grid->mesh_side_length.y;
    // real side_length_z = the_grid->mesh_side_length.z;

    real region_height = (side_length_y - channel_width) / 2.0;

    // Set the conductivities for each cell on the grid
    OMP(parallel for)
    for(i = 0; i < num_active_cells; i++) {
        real sigma_x_new = sigma_x;
        real sigma_y_new = sigma_y;
        real sigma_z_new = sigma_z;

        real x = ac[i]->center.x;
        real y = ac[i]->center.y;
        //        real z = ac[i]->center.z;

        // Check region 3
        inside_3 = (x >= 0.0) && (x < channel_length) && (y >= region_height) && (y <= region_height + channel_width);

        if(inside_3) {
            sigma_x_new = sigma_x * source_factor;
            sigma_y_new = sigma_y * source_factor;
            sigma_z_new = sigma_z * source_factor;
        }

        // Check region 4
        inside_4 = (x >= channel_length) && (x <= side_length_x) && (y >= 0.0) && (y <= side_length_y);

        if(inside_4) {
            sigma_x_new = sigma_x * sink_factor;
            sigma_y_new = sigma_y * sink_factor;
            sigma_z_new = sigma_z * sink_factor;
        }

        ac[i]->sigma.x = sigma_x_new;
        ac[i]->sigma.y = sigma_y_new;
        ac[i]->sigma.z = sigma_z_new;
    }

    // Then, we fill the discretization matrix
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
}

ASSEMBLY_MATRIX(homogeneous_sigma_assembly_matrix) {

    static bool sigma_initialized = false;

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    int i;

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config->config_data, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config->config_data, "sigma_z");

    if(!sigma_initialized) {
        OMP(parallel for)
        for(i = 0; i < num_active_cells; i++) {
            ac[i]->sigma.x = sigma_x;
            ac[i]->sigma.y = sigma_y;
            ac[i]->sigma.z = sigma_z;
        }

        // sigma_initialized = true;
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
    sds *points;
    int split_count;

    while((getline(&line, &len, fibers_file)) != -1) {

        points = sdssplit(line, " ", &split_count);
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

ASSEMBLY_MATRIX(anisotropic_sigma_assembly_matrix) {

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

	char *fiber_file = NULL;
	GET_PARAMETER_VALUE_CHAR_OR_USE_DEFAULT(fiber_file, config->config_data, "fibers_file");

	struct fiber_coords *fibers = NULL;

	GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_l, config->config_data, "sigma_l");
	GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_t, config->config_data, "sigma_t");
	GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, sigma_n, config->config_data, "sigma_n");

	if(fiber_file) {
		log_to_stdout_and_file("Loading mesh fibers\n");
		fibers = read_fibers(fiber_file, false);
	}

	OMP(parallel for private(D))
	for(i = 0; i < num_active_cells; i++) {

		if(fibers) {
			int fiber_index = ac[i]->original_position_in_file;
			// calc_tensor(D, fibers[fiber_index].f, fibers[fiber_index].s, fibers[fiber_index].n, sigma_l, sigma_t, sigma_n);
			calc_tensor2(D, fibers[fiber_index].f, sigma_l, sigma_t);
			ac[i]->sigma.fibers = fibers[fiber_index];
		}
		else {
			//TODO: allow f to be passed from the config file
			real_cpu f[3];
			f[0] = 1;
			f[1] = 1;
			f[2] = 1;
			calc_tensor2(D, f, sigma_l, sigma_t);
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

#ifdef DEBUG_INFO
	FILE *m = fopen("m.txt", "w");
	print_grid_matrix(the_grid, m);
	fclose(m);
#endif

}

ASSEMBLY_MATRIX(homogeneous_sigma_with_a_factor_assembly_matrix) {

    static bool sigma_initialized = false;

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    int i;

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config->config_data, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config->config_data, "sigma_z");

    real sigma_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_factor, config->config_data, "sigma_factor");

    if(!sigma_initialized) {
        OMP(parallel for)
        for(i = 0; i < num_active_cells; i++) {
            ac[i]->sigma.x = sigma_x * sigma_factor;
            ac[i]->sigma.y = sigma_y * sigma_factor;
            ac[i]->sigma.z = sigma_z * sigma_factor;
        }

        // sigma_initialized = true;
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
}

ASSEMBLY_MATRIX(fibrotic_region_with_sigma_factor_assembly_matrix) {

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    int i;

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config->config_data, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config->config_data, "sigma_z");

    real sigma_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_factor, config->config_data, "sigma_factor");

    real min_x = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, min_x, config->config_data, "region_min_x");

    real max_x = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, max_x, config->config_data, "region_max_x");

    real min_y = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, min_y, config->config_data, "region_min_y");

    real max_y = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, max_y, config->config_data, "region_max_y");

    real min_z = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, min_z, config->config_data, "region_min_z");

    real max_z = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, max_z, config->config_data, "region_max_z");

    bool inside;

    // Set the conductivities for each cell on the grid
    OMP(parallel for)
    for(i = 0; i < num_active_cells; i++) {
        real sigma_x_new = sigma_x;
        real sigma_y_new = sigma_y;
        real sigma_z_new = sigma_z;

        real x = ac[i]->center.x;
        real y = ac[i]->center.y;
        real z = ac[i]->center.z;

        // Check if inside the region
        inside = (x >= min_x) && (x <= max_x) && (y >= min_y) && (y <= max_y) && (z >= min_z) && (z <= max_z);

        if(inside) {
            sigma_x_new = sigma_x * sigma_factor;
            sigma_y_new = sigma_y * sigma_factor;
            sigma_z_new = sigma_z * sigma_factor;
        }

        ac[i]->sigma.x = sigma_x_new;
        ac[i]->sigma.y = sigma_y_new;
        ac[i]->sigma.z = sigma_z_new;
    }

    // Then, we fill the discretization matrix
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
}

ASSEMBLY_MATRIX(heterogenous_sigma_with_factor_assembly_matrix) {

    static bool sigma_initialized = false;

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    int i;

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config->config_data, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config->config_data, "sigma_z");

    real_cpu phi = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, phi, config->config_data, "phi");

    unsigned seed = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, seed, config->config_data, "seed");

    real sigma_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_factor, config->config_data, "sigma_factor");

    log_to_stdout_and_file("Reducing conductivity from %.2lf %% of cells\n", phi * 100.0);

    // Initialize the seed for the fibrosis
    srand(seed);

    log_to_stdout_and_file("Using %u as seed\n", seed);

    if(!sigma_initialized) {

        FOR_EACH_CELL(the_grid) {

            if(cell->active) {
                real_cpu p = (real_cpu)(rand()) / (RAND_MAX);
                if(p < phi) {
                    cell->sigma.x = sigma_x * sigma_factor;
                    cell->sigma.y = sigma_y * sigma_factor;
                    cell->sigma.z = sigma_z * sigma_factor;
                } else {
                    cell->sigma.x = sigma_x;
                    cell->sigma.y = sigma_y;
                    cell->sigma.z = sigma_z;
                }
            }
        }

        // sigma_initialized = true;
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
}

// This function will read the fibrotic regions and for each cell that is inside the region we will
// reduce its conductivity value based on the 'sigma_factor'.
ASSEMBLY_MATRIX(heterogenous_sigma_with_factor_assembly_matrix_from_file) {

    static bool sigma_initialized = false;

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    char *fib_file = NULL;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(fib_file, config->config_data, "fibrosis_file");

    int fib_size = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(int, fib_size, config->config_data, "size");

    real sigma_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_x, config->config_data, "sigma_x");

    real sigma_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_y, config->config_data, "sigma_y");

    real sigma_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_z, config->config_data, "sigma_z");

    real sigma_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, sigma_factor, config->config_data, "sigma_factor");

    if(!sigma_initialized) {
        OMP(parallel for)
        for(uint32_t i = 0; i < num_active_cells; i++) {
            ac[i]->sigma.x = sigma_x;
            ac[i]->sigma.y = sigma_y;
            ac[i]->sigma.z = sigma_z;
        }

        // sigma_initialized = true;
    }

    // Reading the fibrotic regions from the input file
    FILE *file = fopen(fib_file, "r");

    if(!file) {
        printf("Error opening file %s!!\n", fib_file);
        exit(0);
    }

    real_cpu **scar_mesh = (real_cpu **)malloc(sizeof(real_cpu *) * fib_size);

    for(int i = 0; i < fib_size; i++) {
        scar_mesh[i] = (real_cpu *)malloc(sizeof(real_cpu) * 7);
        if(scar_mesh[i] == NULL) {
            printf("Failed to allocate memory\n");
            exit(0);
        }
    }

    uint32_t i = 0;
    while(fscanf(file, "%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &scar_mesh[i][0], &scar_mesh[i][1], &scar_mesh[i][2], &scar_mesh[i][3], &scar_mesh[i][4],
                 &scar_mesh[i][5], &scar_mesh[i][6]) != EOF) {
        i++;
    }

    fclose(file);

    uint32_t num_fibrotic_regions = i;

    // Pass through all the cells of the grid and check if its center is inside the current
    // fibrotic region
    OMP(parallel for)
    for(uint32_t j = 0; j < num_fibrotic_regions; j++) {

        struct cell_node *grid_cell = the_grid->first_cell;

        real_cpu b_center_x = scar_mesh[j][0];
        real_cpu b_center_y = scar_mesh[j][1];

        real_cpu b_h_dx = scar_mesh[j][3];
        real_cpu b_h_dy = scar_mesh[j][4];

        bool active = (bool)(scar_mesh[j][6]);

        while(grid_cell != 0) {
            if(grid_cell->active) {
                real_cpu center_x = grid_cell->center.x;
                real_cpu center_y = grid_cell->center.y;
                //               real_cpu half_dx = grid_cell->discretization.x/2.0;
                //                real_cpu half_dy = grid_cell->discretization.y/2.0;

                struct point_3d p;
                struct point_3d q;

                p.x = b_center_y + b_h_dy;
                p.y = b_center_y - b_h_dy;

                q.x = b_center_x + b_h_dx;
                q.y = b_center_x - b_h_dx;

                // Check if the current cell is inside the fibrotic region
                if(center_x > q.y && center_x < q.x && center_y > p.y && center_y < p.x) {
                    if(active == 0) {
                        grid_cell->sigma.x = sigma_x * sigma_factor;
                        grid_cell->sigma.y = sigma_y * sigma_factor;
                        grid_cell->sigma.z = sigma_z * sigma_factor;
                    }
                }
            }

            grid_cell = grid_cell->next;
        }
    }

    OMP(parallel for)
    for(int i = 0; i < num_active_cells; i++) {

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

    for(int k = 0; k < fib_size; k++) {
        free(scar_mesh[k]);
    }

    free(scar_mesh);
}

// This function will generate the fibrotic region file for the 120um x 120um grid by reescaling
// the original Scientific Reports 4b grid from 40000um side_length to 48000um
// rescale_factor = 1.2
ASSEMBLY_MATRIX(heterogenous_fibrotic_region_file_write_with_input_file) {

    initialize_diagonal_elements(the_solver, the_grid);

    char *fib_file = NULL;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(fib_file, config->config_data, "fibrosis_file");

    char *new_fib_file = NULL;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(new_fib_file, config->config_data, "rescaled_fibrosis_file");

    int fib_size = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(int, fib_size, config->config_data, "size");

    real rescale_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, rescale_factor, config->config_data, "rescale_factor");

    FILE *file = fopen(fib_file, "r");

    if(!file) {
        printf("Error opening file %s!!\n", fib_file);
        exit(0);
    }

    // Read and store the original positions of the fibrotic regions
    real_cpu **scar_mesh = (real_cpu **)malloc(sizeof(real_cpu *) * fib_size);

    for(int i = 0; i < fib_size; i++) {
        scar_mesh[i] = (real_cpu *)malloc(sizeof(real_cpu) * 7);
        if(scar_mesh[i] == NULL) {
            printf("Failed to allocate memory\n");
            exit(0);
        }
    }

    for(int i = 0; i < fib_size; i++) {
        fscanf(file, "%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", &scar_mesh[i][0], &scar_mesh[i][1], &scar_mesh[i][2], &scar_mesh[i][3], &scar_mesh[i][4],
               &scar_mesh[i][5], &scar_mesh[i][6]);
    }

    fclose(file);

    // Write the new fibrotic region file based on the 'rescale_factor'
    FILE *fileW = fopen(new_fib_file, "w");

    if(!file) {
        printf("Error opening file %s!!\n", fib_file);
        exit(0);
    }

    // Multiple the positions of each scar by a rescale factor
    for(int i = 0; i < fib_size; i++) {
        scar_mesh[i][0] = scar_mesh[i][0] * rescale_factor;
        scar_mesh[i][1] = scar_mesh[i][1] * rescale_factor;
        scar_mesh[i][2] = scar_mesh[i][2] * rescale_factor;
        scar_mesh[i][3] = scar_mesh[i][3] * rescale_factor;
        scar_mesh[i][4] = scar_mesh[i][4] * rescale_factor;
        scar_mesh[i][5] = scar_mesh[i][5] * rescale_factor;

        fprintf(fileW, "%g,%g,%g,%g,%g,%g,%g\n", scar_mesh[i][0], scar_mesh[i][1], scar_mesh[i][2], scar_mesh[i][3], scar_mesh[i][4], scar_mesh[i][5],
                scar_mesh[i][6]);
    }

    fclose(fileW);

    for(int k = 0; k < fib_size; k++) {
        free(scar_mesh[k]);
    }

    free(scar_mesh);

    // We just leave the program after this ...
    log_to_stdout_and_file("[!] Finish writing new fibrotic region file '%s'!\n", new_fib_file);
    exit(EXIT_SUCCESS);
}

// This function will generate the fibrotic region file for the 120um x 120um grid by reescaling
// the original Scientific Reports 4b grid from 40000um side_length to 48000um. The fibrotic region
// will be mapped using the same idea used on the 'domains_library' function by using a fixed seed
// for the random number generator.
// rescale_factor = 1.2
ASSEMBLY_MATRIX(heterogenous_fibrotic_region_file_write_using_seed) {

    initialize_diagonal_elements(the_solver, the_grid);

    real_cpu phi = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, phi, config->config_data, "phi");

    unsigned seed = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, seed, config->config_data, "seed");

    char *new_fib_file = NULL;
    GET_PARAMETER_VALUE_CHAR_OR_REPORT_ERROR(new_fib_file, config->config_data, "rescaled_fibrosis_file");

    double rescale_factor = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, rescale_factor, config->config_data, "rescale_factor");

    // Write the new fibrotic region file
    FILE *fileW = fopen(new_fib_file, "w+");

    // Initialize the random the generator with the same seed used by the original model
    srand(seed);
    FOR_EACH_CELL(the_grid) {

        if(cell->active) {
            real_cpu p = (real_cpu)(rand()) / (RAND_MAX);
            if(p < phi) {
                // We reescale the cell position using the 'rescale_factor'
                double center_x = cell->center.x * rescale_factor;
                double center_y = cell->center.y * rescale_factor;
                double center_z = cell->center.z * rescale_factor;
                double dx = cell->discretization.x * rescale_factor;
                double dy = cell->discretization.y * rescale_factor;
                double dz = cell->discretization.z * rescale_factor;

                // Then, we write only the fibrotic regions to the output file
                fprintf(fileW, "%g,%g,%g,%g,%g,%g,0\n", center_x, center_y, center_z, dx / 2.0, dy / 2.0, dz / 2.0);
            }
        }
    }

    fclose(fileW);

    // We just leave the program after this ...
    log_to_stdout_and_file("[!] Finish writing fibrotic region file '%s'!\n", new_fib_file);
    exit(EXIT_SUCCESS);
}
