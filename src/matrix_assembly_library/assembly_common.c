static struct element fill_element(uint32_t position, enum transition_direction direction, real_cpu dx, real_cpu dy, real_cpu dz, real_cpu sigma_x,
                                   real_cpu sigma_y, real_cpu sigma_z, struct element *cell_elements, struct cell_node *cell);

static void initialize_diagonal_elements(struct monodomain_solver *the_solver, struct grid *the_grid) {

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
        element.value_ecg = 0.0;

        if(ac[i]->elements)
            arrfree(ac[i]->elements);

        ac[i]->elements = NULL;

        arrsetcap(ac[i]->elements, 7);
        arrput(ac[i]->elements, element);
    }
}

static struct element fill_element(uint32_t position, enum transition_direction direction, real_cpu dx, real_cpu dy, real_cpu dz, real_cpu sigma_x,
                                   real_cpu sigma_y, real_cpu sigma_z, struct element *cell_elements, struct cell_node *cell) {

    real_cpu multiplier;

    struct element new_element;
    new_element.column = position;
    new_element.cell = cell;

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
        }

    } else if(flux_direction == LEFT) {
        if(neighbours[LEFT]) {
            *sigma_1 = (neighbours[LEFT]->sigma.x + cell_node->sigma.x) / 2.0;
            *count_s1 = 1;

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

        }
    } else if(flux_direction == TOP) {
        if(neighbours[TOP]) {
            *sigma_1 = (neighbours[TOP]->sigma.y + cell_node->sigma.y) / 2.0;
            *count_s1 = 1;

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

        }
    } else if(flux_direction == DOWN) {
        if(neighbours[DOWN]) {
            *sigma_1 = (neighbours[DOWN]->sigma.y + cell_node->sigma.y) / 2.0;
            *count_s1 = 1;

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

        }
    } else if(flux_direction == FRONT) {

        if(neighbours[FRONT]) {
            *sigma_1 = (neighbours[FRONT]->sigma.z + cell_node->sigma.z) / 2.0;
            *count_s1 = 1;

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

        }
    } else if(flux_direction == BACK) {
        if(neighbours[BACK]) {
            *sigma_1 = (neighbours[BACK]->sigma.z + cell_node->sigma.z) / 2.0;
            *count_s1 = 1;

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
}

static inline real_cpu DIVIDE(real_cpu num, real_cpu denom) {
    if(denom != 0) {
        return num / denom;
    }
    return 0.0;
}

#define UPDATE_OR_ADD_ELEMENT(g_cell, n_cell, v)                                                                                                               \
    do {                                                                                                                                                       \
        if(v != 0) {                                                                                                                                           \
            struct element el;                                                                                                                                 \
            int el_index = find_neighbour_index(g_cell, n_cell);                                                                                               \
            if(el_index != -1) {                                                                                                                               \
                g_cell->elements[el_index].value += v;                                                                                                         \
            } else {                                                                                                                                           \
                el.value = v;                                                                                                                                  \
                el.cell = n_cell;                                                                                                                              \
                el.column = n_cell->grid_position;                                                                                                             \
                arrput(g_cell->elements, el);                                                                                                                  \
            }                                                                                                                                                  \
        }                                                                                                                                                      \
    } while(0)

static void fill_elements_aniso(struct cell_node *grid_cell, struct cell_node *neighbours[26]) {

    struct element *elements = grid_cell->elements;

    real_cpu dx = grid_cell->discretization.x;
    real_cpu dy = grid_cell->discretization.y;
    real_cpu dz = grid_cell->discretization.z;

    real_cpu dy_times_dz_over_dx = (dy * dz) / dx;
    real_cpu dx_times_dz_over_dy = (dx * dz) / dy;
    real_cpu dx_times_dy_over_dz = (dx * dy) / dz;

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
    calc_sigmas(grid_cell, neighbours, LEFT, &sigma_x_l, &sigma_xy_jx_l, &sigma_xz_jx_l, &count_sigma_x_l, &count_sigma_xy_jx_l, &count_sigma_xz_jx_l);

    calc_sigmas(grid_cell, neighbours, TOP, &sigma_y_t, &sigma_xy_jy_t, &sigma_yz_jy_t, &count_sigma_y_t, &count_sigma_xy_jy_t, &count_sigma_yz_jy_t);
    calc_sigmas(grid_cell, neighbours, DOWN, &sigma_y_d, &sigma_xy_jy_d, &sigma_yz_jy_d, &count_sigma_y_d, &count_sigma_xy_jy_d, &count_sigma_yz_jy_d);

    calc_sigmas(grid_cell, neighbours, FRONT, &sigma_z_f, &sigma_xz_jz_f, &sigma_yz_jz_f, &count_sigma_z_f, &count_sigma_xz_jz_f, &count_sigma_yz_jz_f);
    calc_sigmas(grid_cell, neighbours, BACK, &sigma_z_b, &sigma_xz_jz_b, &sigma_yz_jz_b, &count_sigma_z_b, &count_sigma_xz_jz_b, &count_sigma_yz_jz_b);

    // MAIN DIAGONAL
    elements[0].value += dy_times_dz_over_dx * sigma_x_r + dy_times_dz_over_dx * sigma_x_l + dx_times_dz_over_dy * sigma_y_t + dx_times_dz_over_dy * sigma_y_d +
                         dx_times_dy_over_dz * sigma_z_f + dx_times_dy_over_dz * sigma_z_b;

    real_cpu s1, s2, s3;

    // All neighbours
    for(int direction = 0; direction < NUM_DIRECTIONS; direction++) {

        if(neighbours[direction]) {

            struct element new_element;
            new_element.value = 0.0;

            switch(direction) {
            case FRONT:
                new_element.value += -sigma_z_f * dx_times_dy_over_dz;

                if(neighbours[BACK]) {
                    new_element.value += -DIVIDE(sigma_yz_jy_t, count_sigma_yz_jy_t) * dx + DIVIDE(sigma_yz_jy_d, count_sigma_yz_jy_d) * dx -
                                         DIVIDE(sigma_xz_jx_r, count_sigma_xz_jx_r) * dy + DIVIDE(sigma_xz_jx_l, count_sigma_xz_jx_l) * dy;
                }

                break;
            case BACK:
                new_element.value += -sigma_z_b * dx_times_dy_over_dz;

                if(neighbours[FRONT]) {
                    new_element.value += DIVIDE(sigma_yz_jy_t, count_sigma_yz_jy_t) * dx - DIVIDE(sigma_yz_jy_d, count_sigma_yz_jy_d) * dx +
                                         DIVIDE(sigma_xz_jx_r, count_sigma_xz_jx_r) * dy - DIVIDE(sigma_xz_jx_l, count_sigma_xz_jx_l) * dy;
                }
                break;

            case TOP:
                new_element.value += -sigma_y_t * dx_times_dz_over_dy;

                if(neighbours[DOWN]) {
                    new_element.value += -DIVIDE(sigma_yz_jz_f, count_sigma_yz_jz_f) * dx + DIVIDE(sigma_yz_jz_b, count_sigma_yz_jz_b) * dx -
                                         DIVIDE(sigma_xy_jx_r, count_sigma_xy_jx_r) * dz + DIVIDE(sigma_xy_jx_l, count_sigma_xy_jx_l) * dz;
                }
                break;
            case DOWN:
                new_element.value += -sigma_y_d * dx_times_dz_over_dy;
                if(neighbours[TOP]) {
                    new_element.value += DIVIDE(sigma_yz_jz_f, count_sigma_yz_jz_f) * dx - DIVIDE(sigma_yz_jz_b, count_sigma_yz_jz_b) * dx +
                                         DIVIDE(sigma_xy_jx_r, count_sigma_xy_jx_r) * dz - DIVIDE(sigma_xy_jx_l, count_sigma_xy_jx_l) * dz;
                }
                break;
            case RIGHT:
                new_element.value += -sigma_x_r * dy_times_dz_over_dx;
                if(neighbours[LEFT]) {
                    new_element.value += -DIVIDE(sigma_xz_jz_f, count_sigma_xz_jz_f) * dy + DIVIDE(sigma_xz_jz_b, count_sigma_xz_jz_b) * dy -
                                         DIVIDE(sigma_xy_jy_t, count_sigma_xy_jy_t) * dz + DIVIDE(sigma_xy_jy_d, count_sigma_xy_jy_d) * dz;
                }
                break;
            case LEFT:
                new_element.value += -sigma_x_l * dy_times_dz_over_dx;
                if(neighbours[RIGHT]) {
                    new_element.value += DIVIDE(sigma_xz_jz_f, count_sigma_xz_jz_f) * dy - DIVIDE(sigma_xz_jz_b, count_sigma_xz_jz_b) * dy +
                                         DIVIDE(sigma_xy_jy_t, count_sigma_xy_jy_t) * dz - DIVIDE(sigma_xy_jy_d, count_sigma_xy_jy_d) * dz;
                }
                break;

            case FRONT_TOP:

                s1 = 0.0;
                s2 = 0.0;

                if(neighbours[FRONT_DOWN]) {
                    s1 += -DIVIDE(sigma_xy_jx_r, count_sigma_xy_jx_r) * dz + DIVIDE(sigma_xy_jx_l, count_sigma_xy_jx_l) * dz -
                          DIVIDE(sigma_yz_jz_f, count_sigma_yz_jz_f) * dx;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[FRONT_DOWN], -s1);
                }

                if(neighbours[BACK_TOP]) {
                    s2 += -DIVIDE(sigma_yz_jy_t, count_sigma_yz_jy_t) * dx - DIVIDE(sigma_xz_jx_r, count_sigma_xz_jx_r) * dy +
                          DIVIDE(sigma_xz_jx_l, count_sigma_xz_jx_l) * dy;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_TOP], -s2);
                }

                new_element.value += s1 + s2;

                break;

            case FRONT_DOWN:

                s1 = 0.0;

                if(neighbours[BACK_DOWN]) {
                    s1 += DIVIDE(sigma_yz_jy_d, count_sigma_yz_jy_d) * dx - DIVIDE(sigma_xz_jx_r, count_sigma_xz_jx_r) * dy +
                          DIVIDE(sigma_xz_jx_l, count_sigma_xz_jx_l) * dy;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[FRONT_DOWN], s1);
                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_DOWN], -s1);
                }

                break;

            case BACK_TOP:

                s1 = 0.0;

                if(neighbours[BACK_DOWN]) {

                    s1 += DIVIDE(sigma_yz_jz_b, count_sigma_yz_jz_b) * dx - DIVIDE(sigma_xy_jx_r, count_sigma_xy_jx_r) * dz +
                          DIVIDE(sigma_xy_jx_l, count_sigma_xy_jx_l) * dz;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_TOP], s1);
                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_DOWN], -s1);
                }

                break;
            case BACK_DOWN:
                break; // alread handled by the above cases

            case FRONT_RIGHT:
                s1 = 0.0;
                s2 = 0.0;

                if(neighbours[BACK_RIGHT]) {
                    s1 += -DIVIDE(sigma_xz_jx_r, count_sigma_xz_jx_r) * dy - DIVIDE(sigma_yz_jy_t, count_sigma_yz_jy_t) * dx +
                          DIVIDE(sigma_yz_jy_d, count_sigma_yz_jy_d) * dx;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_RIGHT], -s1);
                }

                if(neighbours[FRONT_LEFT]) {
                    s2 += -DIVIDE(sigma_xz_jz_f, count_sigma_xz_jz_f) * dy - DIVIDE(sigma_xy_jy_t, count_sigma_xy_jy_t) * dz +
                          DIVIDE(sigma_xy_jy_d, count_sigma_xy_jy_d) * dz;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[FRONT_LEFT], -s2);
                }

                new_element.value += s1 + s2;
                break;

            case FRONT_LEFT:

                s1 = 0.0;

                if(neighbours[BACK_LEFT]) {

                    s1 += -DIVIDE(sigma_yz_jy_t, count_sigma_yz_jy_t) * dx + DIVIDE(sigma_yz_jy_d, count_sigma_yz_jy_d) * dx +
                          DIVIDE(sigma_xz_jx_l, count_sigma_xz_jx_l) * dy;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[FRONT_LEFT], s1);
                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_LEFT], -s1);
                }
                break;
            case BACK_RIGHT:
                s1 = 0.0;
                if(neighbours[BACK_LEFT]) {
                    s1 += DIVIDE(sigma_xz_jz_b, count_sigma_xz_jz_b) * dy - DIVIDE(sigma_xy_jy_t, count_sigma_xy_jy_t) * dz +
                          DIVIDE(sigma_xy_jy_d, count_sigma_xy_jy_d) * dz;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_RIGHT], s1);
                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_LEFT], -s1);
                }

                break;

            case BACK_LEFT:
                break; // alread handled by the above cases

            case TOP_RIGHT:

                s1 = 0;
                s2 = 0;
                if(neighbours[DOWN_RIGHT]) {
                    s1 += -DIVIDE(sigma_xy_jx_r, count_sigma_xy_jx_r) * dz - DIVIDE(sigma_yz_jz_f, count_sigma_yz_jz_f) * dx +
                          DIVIDE(sigma_yz_jz_b, count_sigma_yz_jz_b) * dx;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[DOWN_RIGHT], -s1);
                }

                if(neighbours[TOP_LEFT]) {
                    s2 += -DIVIDE(sigma_xz_jz_f, count_sigma_xz_jz_f) * dy + DIVIDE(sigma_xz_jz_b, count_sigma_xz_jz_b) * dy -
                          DIVIDE(sigma_xy_jy_t, count_sigma_xy_jy_t) * dz;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[TOP_LEFT], -s2);
                }

                new_element.value += s1 + s2;

                break;
            case TOP_LEFT:

                s1 = 0;
                if(neighbours[DOWN_LEFT]) {
                    s1 += -DIVIDE(sigma_yz_jz_f, count_sigma_yz_jz_f) * dx + DIVIDE(sigma_yz_jz_b, count_sigma_yz_jz_b) * dx +
                          DIVIDE(sigma_xy_jx_l, count_sigma_xy_jx_l) * dz;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[TOP_LEFT], s1);
                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[DOWN_LEFT], -s1);
                }
                break;
            case DOWN_RIGHT:
                s1 = 0;
                if(neighbours[DOWN_LEFT]) {
                    s1 += -DIVIDE(sigma_xz_jz_f, count_sigma_xz_jz_f) * dy + DIVIDE(sigma_xz_jz_b, count_sigma_xz_jz_b) * dy +
                          DIVIDE(sigma_xy_jy_d, count_sigma_xy_jy_d) * dz;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[DOWN_RIGHT], s1);
                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[DOWN_LEFT], -s1);
                }
                break;
            case DOWN_LEFT:
                break; // alread handled by the above cases

            case FRONT_TOP_RIGHT:
                s1 = 0;
                s2 = 0;
                s3 = 0;
                if(neighbours[FRONT_DOWN_RIGHT]) {
                    s1 += -DIVIDE(sigma_xy_jx_r, count_sigma_xy_jx_r) * dz - DIVIDE(sigma_yz_jz_f, count_sigma_yz_jz_f) * dx;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[FRONT_DOWN_RIGHT], -s1);
                }
                if(neighbours[BACK_TOP_RIGHT]) {
                    s2 += -DIVIDE(sigma_xz_jx_r, count_sigma_xz_jx_r) * dy - DIVIDE(sigma_yz_jy_t, count_sigma_yz_jy_t) * dx;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_TOP_RIGHT], -s2);
                }
                if(neighbours[FRONT_TOP_LEFT]) {
                    s3 += -DIVIDE(sigma_xz_jz_f, count_sigma_xz_jz_f) * dy - DIVIDE(sigma_xy_jy_t, count_sigma_xy_jy_t) * dz;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[FRONT_TOP_LEFT], -s3);
                }
                new_element.value += s1 + s2 + s3;
                break;

            case FRONT_TOP_LEFT:

                s1 = 0;
                s2 = 0;

                if(neighbours[FRONT_DOWN_LEFT]) {
                    s1 += DIVIDE(sigma_xy_jx_l, count_sigma_xy_jx_l) * dz - DIVIDE(sigma_yz_jz_f, count_sigma_yz_jz_f) * dx;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[FRONT_TOP_LEFT], s1);
                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[FRONT_DOWN_LEFT], -s1);
                }

                if(neighbours[BACK_TOP_LEFT]) {
                    s2 += DIVIDE(sigma_xz_jx_l, count_sigma_xz_jx_l) * dy - DIVIDE(sigma_yz_jy_t, count_sigma_yz_jy_t) * dx;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[FRONT_TOP_LEFT], s2);
                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_TOP_LEFT], -s2);
                }

                break;
            case FRONT_DOWN_RIGHT:

                s1 = 0;
                s2 = 0;

                if(neighbours[BACK_DOWN_RIGHT]) {
                    s1 += -DIVIDE(sigma_xz_jx_r, count_sigma_xz_jx_r) * dy + DIVIDE(sigma_yz_jy_d, count_sigma_yz_jy_d) * dx;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[FRONT_DOWN_RIGHT], s1);
                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_DOWN_RIGHT], -s1);
                }

                if(neighbours[FRONT_DOWN_LEFT]) {
                    s2 += DIVIDE(sigma_xy_jy_d, count_sigma_xy_jy_d) * dz - DIVIDE(sigma_xz_jz_f, count_sigma_xz_jz_f) * dy;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[FRONT_DOWN_RIGHT], s2);
                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[FRONT_DOWN_LEFT], -s2);
                }

                break;
            case FRONT_DOWN_LEFT:
                s1 = 0;
                if(neighbours[BACK_DOWN_LEFT]) {
                    s1 += DIVIDE(sigma_xz_jx_l, count_sigma_xz_jx_l) * dy + DIVIDE(sigma_yz_jy_d, count_sigma_yz_jy_d) * dx;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[FRONT_DOWN_LEFT], s1);
                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_DOWN_LEFT], -s1);
                }
                break;
            case BACK_TOP_RIGHT:
                s1 = 0;
                s2 = 0;

                if(neighbours[BACK_DOWN_RIGHT]) {
                    s1 += -DIVIDE(sigma_xy_jx_r, count_sigma_xy_jx_r) * dz + DIVIDE(sigma_yz_jz_b, count_sigma_yz_jz_b) * dx;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_TOP_RIGHT], s1);
                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_DOWN_RIGHT], -s1);
                }

                if(neighbours[BACK_TOP_LEFT]) {
                    s2 += -DIVIDE(sigma_xy_jy_t, count_sigma_xy_jy_t) * dz + DIVIDE(sigma_xz_jz_b, count_sigma_xz_jz_b) * dy;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_TOP_RIGHT], s2);
                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_TOP_LEFT], -s2);
                }
                break;
            case BACK_TOP_LEFT:
                s1 = 0;
                if(neighbours[BACK_DOWN_LEFT]) {
                    s1 += DIVIDE(sigma_xy_jx_l, count_sigma_xy_jx_l) * dz + DIVIDE(sigma_yz_jz_b, count_sigma_yz_jz_b) * dx;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_TOP_LEFT], s1);
                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_DOWN_LEFT], -s1);
                }

                break;
            case BACK_DOWN_RIGHT:
                s1 = 0;
                if(neighbours[BACK_DOWN_LEFT]) {
                    s1 += DIVIDE(sigma_xy_jy_d, count_sigma_xy_jy_d) * dz + DIVIDE(sigma_xz_jz_b, count_sigma_xz_jz_b) * dy;

                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_DOWN_RIGHT], s1);
                    UPDATE_OR_ADD_ELEMENT(grid_cell, neighbours[BACK_DOWN_LEFT], -s1);
                }
                break;

            case BACK_DOWN_LEFT:
                break;

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
                struct element new_element = fill_element(position, direction, dx, dy, dz, sigma_x, sigma_y, sigma_z, cell_elements, black_neighbor_cell);
                //new_element.cell = black_neighbor_cell;
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
                struct element new_element = fill_element(position, direction, dx, dy, dz, sigma_x, sigma_y, sigma_z, cell_elements, grid_cell);
                //new_element.cell = grid_cell;
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

static void outer_product_vector_vector_t(real_cpu p[3][3], real_cpu v[3]) {
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
