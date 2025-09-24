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

static void calc_sigmas(struct cell_node *cell_node, struct cell_node *neighbours[26],
                        enum transition_direction flux_direction, real_cpu *sigma_1,
                        real_cpu *sigma_2, real_cpu *sigma_3, int *count_s1, int *count_s2, int *count_s3) {

    *sigma_1 = 0.0;
    *sigma_2 = 0.0;
    *sigma_3 = 0.0;

    // Initialize counts to 0
    *count_s1 = 0;
    *count_s2 = 0;
    *count_s3 = 0;

    switch(flux_direction) {
    case RIGHT:
        if(neighbours[RIGHT]) {
            // Main conductivity (x-direction)
            *sigma_1 = (neighbours[RIGHT]->sigma.x + cell_node->sigma.x) / 2.0;
            *count_s1 = 1;

            // xy cross-derivative at RIGHT face
            *sigma_2 = (neighbours[RIGHT]->sigma.xy + cell_node->sigma.xy) / 2.0;
            *count_s2 = 1;

            // xz cross-derivative at RIGHT face
            *sigma_3 = (neighbours[RIGHT]->sigma.xz + cell_node->sigma.xz) / 2.0;
            *count_s3 = 1;
        }
        break;

    case LEFT:
        if(neighbours[LEFT]) {
            *sigma_1 = (neighbours[LEFT]->sigma.x + cell_node->sigma.x) / 2.0;
            *count_s1 = 1;

            *sigma_2 = (neighbours[LEFT]->sigma.xy + cell_node->sigma.xy) / 2.0;
            *count_s2 = 1;

            *sigma_3 = (neighbours[LEFT]->sigma.xz + cell_node->sigma.xz) / 2.0;
            *count_s3 = 1;
        }
        break;

    case TOP:
        if(neighbours[TOP]) {
            *sigma_1 = (neighbours[TOP]->sigma.y + cell_node->sigma.y) / 2.0;
            *count_s1 = 1;

            *sigma_2 = (neighbours[TOP]->sigma.xy + cell_node->sigma.xy) / 2.0;
            *count_s2 = 1;

            *sigma_3 = (neighbours[TOP]->sigma.yz + cell_node->sigma.yz) / 2.0;
            *count_s3 = 1;
        }
        break;

    case DOWN:
        if(neighbours[DOWN]) {
            *sigma_1 = (neighbours[DOWN]->sigma.y + cell_node->sigma.y) / 2.0;
            *count_s1 = 1;

            *sigma_2 = (neighbours[DOWN]->sigma.xy + cell_node->sigma.xy) / 2.0;
            *count_s2 = 1;

            *sigma_3 = (neighbours[DOWN]->sigma.yz + cell_node->sigma.yz) / 2.0;
            *count_s3 = 1;
        }
        break;

    case FRONT:
        if(neighbours[FRONT]) {
            *sigma_1 = (neighbours[FRONT]->sigma.z + cell_node->sigma.z) / 2.0;
            *count_s1 = 1;

            *sigma_2 = (neighbours[FRONT]->sigma.xz + cell_node->sigma.xz) / 2.0;
            *count_s2 = 1;

            *sigma_3 = (neighbours[FRONT]->sigma.yz + cell_node->sigma.yz) / 2.0;
            *count_s3 = 1;
        }
        break;

    case BACK:
        if(neighbours[BACK]) {
            *sigma_1 = (neighbours[BACK]->sigma.z + cell_node->sigma.z) / 2.0;
            *count_s1 = 1;

            *sigma_2 = (neighbours[BACK]->sigma.xz + cell_node->sigma.xz) / 2.0;
            *count_s2 = 1;

            *sigma_3 = (neighbours[BACK]->sigma.yz + cell_node->sigma.yz) / 2.0;
            *count_s3 = 1;
        }
        break;

    default:
        break;
    }
}

static void fill_elements_aniso(struct cell_node *grid_cell, struct cell_node *neighbours[26]) {

    struct element *elements = grid_cell->elements;

    real_cpu dx = grid_cell->discretization.x;
    real_cpu dy = grid_cell->discretization.y;
    real_cpu dz = grid_cell->discretization.z;

    real_cpu dy_times_dz_over_dx = (dy * dz) / dx;
    real_cpu dx_times_dz_over_dy = (dx * dz) / dy;
    real_cpu dx_times_dy_over_dz = (dx * dy) / dz;

    // Calculate all sigma values first
    real_cpu sigma_x_r = 0.0, sigma_xy_jx_r = 0.0, sigma_xz_jx_r = 0.0;
    real_cpu sigma_x_l = 0.0, sigma_xy_jx_l = 0.0, sigma_xz_jx_l = 0.0;
    real_cpu sigma_xy_jy_t = 0.0, sigma_y_t = 0.0, sigma_yz_jy_t = 0.0;
    real_cpu sigma_xy_jy_d = 0.0, sigma_y_d = 0.0, sigma_yz_jy_d = 0.0;
    real_cpu sigma_xz_jz_f = 0.0, sigma_yz_jz_f = 0.0, sigma_z_f = 0.0;
    real_cpu sigma_xz_jz_b = 0.0, sigma_yz_jz_b = 0.0, sigma_z_b = 0.0;

    int count_sigma_x_r = 0, count_sigma_xy_jx_r = 0, count_sigma_xz_jx_r = 0;
    int count_sigma_x_l = 0, count_sigma_xy_jx_l = 0, count_sigma_xz_jx_l = 0;
    int count_sigma_xy_jy_t = 0, count_sigma_y_t = 0, count_sigma_yz_jy_t = 0;
    int count_sigma_xy_jy_d = 0, count_sigma_y_d = 0, count_sigma_yz_jy_d = 0;
    int count_sigma_xz_jz_f = 0, count_sigma_yz_jz_f = 0, count_sigma_z_f = 0;
    int count_sigma_xz_jz_b = 0, count_sigma_yz_jz_b = 0, count_sigma_z_b = 0;

    calc_sigmas(grid_cell, neighbours, RIGHT, &sigma_x_r, &sigma_xy_jx_r, &sigma_xz_jx_r, &count_sigma_x_r, &count_sigma_xy_jx_r, &count_sigma_xz_jx_r);
    calc_sigmas(grid_cell, neighbours, LEFT, &sigma_x_l, &sigma_xy_jx_l, &sigma_xz_jx_l, &count_sigma_x_l, &count_sigma_xy_jx_l, &count_sigma_xz_jx_l);
    calc_sigmas(grid_cell, neighbours, TOP, &sigma_y_t, &sigma_xy_jy_t, &sigma_yz_jy_t, &count_sigma_y_t, &count_sigma_xy_jy_t, &count_sigma_yz_jy_t);
    calc_sigmas(grid_cell, neighbours, DOWN, &sigma_y_d, &sigma_xy_jy_d, &sigma_yz_jy_d, &count_sigma_y_d, &count_sigma_xy_jy_d, &count_sigma_yz_jy_d);
    calc_sigmas(grid_cell, neighbours, FRONT, &sigma_z_f, &sigma_xz_jz_f, &sigma_yz_jz_f, &count_sigma_z_f, &count_sigma_xz_jz_f, &count_sigma_yz_jz_f);
    calc_sigmas(grid_cell, neighbours, BACK, &sigma_z_b, &sigma_xz_jz_b, &sigma_yz_jz_b, &count_sigma_z_b, &count_sigma_xz_jz_b, &count_sigma_yz_jz_b);

    // MAIN DIAGONAL
    elements[0].value += dy_times_dz_over_dx * sigma_x_r + dy_times_dz_over_dx * sigma_x_l +
                         dx_times_dz_over_dy * sigma_y_t + dx_times_dz_over_dy * sigma_y_d +
                         dx_times_dy_over_dz * sigma_z_f + dx_times_dy_over_dz * sigma_z_b;

    // Helper function to safely divide and return 0 if count is 0
    #define SAFE_DIVIDE(sigma, count) ((count) > 0 ? (sigma) / (count) : 0.0)

    // Process face neighbors first
    for(int direction = 0; direction < 6; direction++) { // Only face neighbors
        if(neighbours[direction]) {
            struct element new_element;
            new_element.value = 0.0;
            new_element.column = neighbours[direction]->grid_position;
            new_element.cell = neighbours[direction];

            switch(direction) {
            case FRONT:
                new_element.value = -sigma_z_f * dx_times_dy_over_dz;
                break;
            case BACK:
                new_element.value = -sigma_z_b * dx_times_dy_over_dz;
                break;
            case TOP:
                new_element.value = -sigma_y_t * dx_times_dz_over_dy;
                break;
            case DOWN:
                new_element.value = -sigma_y_d * dx_times_dz_over_dy;
                break;
            case RIGHT:
                new_element.value = -sigma_x_r * dy_times_dz_over_dx;
                break;
            case LEFT:
                new_element.value = -sigma_x_l * dy_times_dz_over_dx;
                break;
            }

            if(new_element.value != 0.0) {
                arrput(grid_cell->elements, new_element);
            }
        }
    }

    for(int direction = 6; direction < 18; direction++) {
        if(neighbours[direction]) {
            struct element new_element;
            new_element.value = 0.0;
            new_element.column = neighbours[direction]->grid_position;
            new_element.cell = neighbours[direction];

            switch(direction) {
            case FRONT_TOP:
                new_element.value = SAFE_DIVIDE(sigma_yz_jz_f, count_sigma_yz_jz_f) * dx +
                                   SAFE_DIVIDE(sigma_xy_jy_t, count_sigma_xy_jy_t) * dz;
                break;
            case FRONT_DOWN:
                new_element.value = -SAFE_DIVIDE(sigma_yz_jz_f, count_sigma_yz_jz_f) * dx +
                                    SAFE_DIVIDE(sigma_xy_jy_d, count_sigma_xy_jy_d) * dz;
                break;
            case BACK_TOP:
                new_element.value = -SAFE_DIVIDE(sigma_yz_jz_b, count_sigma_yz_jz_b) * dx +
                                    SAFE_DIVIDE(sigma_xy_jy_t, count_sigma_xy_jy_t) * dz;
                break;
            case BACK_DOWN:
                new_element.value = SAFE_DIVIDE(sigma_yz_jz_b, count_sigma_yz_jz_b) * dx +
                                   SAFE_DIVIDE(sigma_xy_jy_d, count_sigma_xy_jy_d) * dz;
                break;
            case FRONT_RIGHT:
                new_element.value = SAFE_DIVIDE(sigma_xz_jz_f, count_sigma_xz_jz_f) * dy +
                                   SAFE_DIVIDE(sigma_xy_jx_r, count_sigma_xy_jx_r) * dz;
                break;
            case FRONT_LEFT:
                new_element.value = -SAFE_DIVIDE(sigma_xz_jz_f, count_sigma_xz_jz_f) * dy +
                                    SAFE_DIVIDE(sigma_xy_jx_l, count_sigma_xy_jx_l) * dz;
                break;
            case BACK_RIGHT:
                new_element.value = -SAFE_DIVIDE(sigma_xz_jz_b, count_sigma_xz_jz_b) * dy +
                                    SAFE_DIVIDE(sigma_xy_jx_r, count_sigma_xy_jx_r) * dz;
                break;
            case BACK_LEFT:
                new_element.value = SAFE_DIVIDE(sigma_xz_jz_b, count_sigma_xz_jz_b) * dy +
                                   SAFE_DIVIDE(sigma_xy_jx_l, count_sigma_xy_jx_l) * dz;
                break;
            case TOP_RIGHT:
                new_element.value = SAFE_DIVIDE(sigma_yz_jy_t, count_sigma_yz_jy_t) * dx +
                                   SAFE_DIVIDE(sigma_xz_jx_r, count_sigma_xz_jx_r) * dy;
                break;
            case TOP_LEFT:
                new_element.value = -SAFE_DIVIDE(sigma_yz_jy_t, count_sigma_yz_jy_t) * dx +
                                    SAFE_DIVIDE(sigma_xz_jx_l, count_sigma_xz_jx_l) * dy;
                break;
            case DOWN_RIGHT:
                new_element.value = -SAFE_DIVIDE(sigma_yz_jy_d, count_sigma_yz_jy_d) * dx +
                                    SAFE_DIVIDE(sigma_xz_jx_r, count_sigma_xz_jx_r) * dy;
                break;
            case DOWN_LEFT:
                new_element.value = SAFE_DIVIDE(sigma_yz_jy_d, count_sigma_yz_jy_d) * dx +
                                   SAFE_DIVIDE(sigma_xz_jx_l, count_sigma_xz_jx_l) * dy;
                break;
            }

            if(new_element.value != 0.0) {
                arrput(grid_cell->elements, new_element);
            }
        }
    }

    for(int direction = 18; direction < 26; direction++) {
        if(neighbours[direction]) {
            struct element new_element;
            new_element.value = 0.0;
            new_element.column = neighbours[direction]->grid_position;
            new_element.cell = neighbours[direction];

            switch(direction) {
            case FRONT_TOP_RIGHT:
                new_element.value = SAFE_DIVIDE(sigma_yz_jz_f, count_sigma_yz_jz_f) *
                                   SAFE_DIVIDE(sigma_xy_jx_r, count_sigma_xy_jx_r) / 8.0;
                break;
            case FRONT_TOP_LEFT:
                new_element.value = SAFE_DIVIDE(sigma_yz_jz_f, count_sigma_yz_jz_f) *
                                   SAFE_DIVIDE(sigma_xy_jx_l, count_sigma_xy_jx_l) / 8.0;
                break;
            case FRONT_DOWN_RIGHT:
                new_element.value = -SAFE_DIVIDE(sigma_yz_jz_f, count_sigma_yz_jz_f) *
                                    SAFE_DIVIDE(sigma_xy_jx_r, count_sigma_xy_jx_r) / 8.0;
                break;
            case FRONT_DOWN_LEFT:
                new_element.value = -SAFE_DIVIDE(sigma_yz_jz_f, count_sigma_yz_jz_f) *
                                    SAFE_DIVIDE(sigma_xy_jx_l, count_sigma_xy_jx_l) / 8.0;
                break;
            case BACK_TOP_RIGHT:
                new_element.value = -SAFE_DIVIDE(sigma_yz_jz_b, count_sigma_yz_jz_b) *
                                    SAFE_DIVIDE(sigma_xy_jx_r, count_sigma_xy_jx_r) / 8.0;
                break;
            case BACK_TOP_LEFT:
                new_element.value = -SAFE_DIVIDE(sigma_yz_jz_b, count_sigma_yz_jz_b) *
                                    SAFE_DIVIDE(sigma_xy_jx_l, count_sigma_xy_jx_l) / 8.0;
                break;
            case BACK_DOWN_RIGHT:
                new_element.value = SAFE_DIVIDE(sigma_yz_jz_b, count_sigma_yz_jz_b) *
                                   SAFE_DIVIDE(sigma_xy_jx_r, count_sigma_xy_jx_r) / 8.0;
                break;
            case BACK_DOWN_LEFT:
                new_element.value = SAFE_DIVIDE(sigma_yz_jz_b, count_sigma_yz_jz_b) *
                                   SAFE_DIVIDE(sigma_xy_jx_l, count_sigma_xy_jx_l) / 8.0;
                break;
            }

            if(new_element.value != 0.0) {
                arrput(grid_cell->elements, new_element);
            }
        }
    }

    #undef SAFE_DIVIDE
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

// Albert`s code
static struct fiber_coords_scale *read_fibers_scale(char *fiber_file_path) {

    FILE *fibers_file = open_file_or_exit(fiber_file_path, "r");

    struct fiber_coords_scale *fibers = NULL;
    char *line = NULL;
    size_t len;

    while((getline(&line, &len, fibers_file)) != -1) {

        int split_count;
        sds *points = sdssplit(line, " ", &split_count);
        struct fiber_coords_scale f_coords;

        f_coords.f[0] = strtod(points[0], NULL);
        f_coords.f[1] = strtod(points[1], NULL);
        f_coords.f[2] = strtod(points[2], NULL);

        f_coords.s[0] = strtod(points[3], NULL);
        f_coords.s[1] = strtod(points[4], NULL);
        f_coords.s[2] = strtod(points[5], NULL);

        f_coords.n[0] = strtod(points[6], NULL);
        f_coords.n[1] = strtod(points[7], NULL);
        f_coords.n[2] = strtod(points[8], NULL);

        f_coords.x[0] = strtod(points[9], NULL);
        f_coords.x[1] = strtod(points[10], NULL);
        f_coords.x[2] = strtod(points[11], NULL);

        arrput(fibers, f_coords);
        sdsfreesplitres(points, split_count);
    }

    free(line);

    return fibers;
}
