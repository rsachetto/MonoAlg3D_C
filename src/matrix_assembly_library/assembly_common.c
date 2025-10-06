#include "../alg/cell/cell.h"

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

static void add_single_entry(struct cell_node *from, struct cell_node *to, real_cpu value) {
    if (from == NULL || to == NULL || fabs(value) < 1e-16) {
        return;
    }

    // Check if entry already exists
    struct element *elements = from->elements;
    size_t num_elements = arrlen(elements);
    bool found = false;

    for(size_t i = 1; i < num_elements; i++) {  // Skip diagonal (i=0)
        if(elements[i].column == to->grid_position) {
            elements[i].value += value;  // Add to existing entry
            found = true;
            break;
        }
    }

    // If not found, create new entry
    if(!found) {
        struct element new_element;
        new_element.column = to->grid_position;
        new_element.cell = to;
        new_element.value = value;
        arrput(from->elements, new_element);
    }

    // Always subtract from diagonal
    from->elements[0].value -= value;
}

// Helper function for harmonic mean of 4 values
static real_cpu harmonic_mean_4(real_cpu a, real_cpu b, real_cpu c, real_cpu d) {
    // Harmonic mean = 4 / (1/a + 1/b + 1/c + 1/d)
    // Handle zero/small values to avoid division by zero
    real_cpu eps = 1e-16;
    if (fabs(a) < eps) a = eps;
    if (fabs(b) < eps) b = eps;
    if (fabs(c) < eps) c = eps;
    if (fabs(d) < eps) d = eps;

    return 4.0 / (1.0/a + 1.0/b + 1.0/c + 1.0/d);
}

static real_cpu harmonic_mean_2(real_cpu a, real_cpu b) {
    // Handle zero/small values to avoid division by zero
    real_cpu eps = 1e-16;
    if (fabs(a) < eps) a = eps;
    if (fabs(b) < eps) b = eps;
    return 2.0 / (1.0/a + 1.0/b);
}

static void fill_discretization_matrix_elements_aniso(struct cell_node *cell_i) {
    struct cell_node *neighbours[26];

    // Get all neighbors
    for(int direction = 0; direction < NUM_DIRECTIONS; direction++) {
        struct cell_node *neighbour = get_cell_neighbour_with_same_refinement_level(cell_i, direction);
        if(neighbour && neighbour->active) {
            neighbours[direction] = neighbour;
        } else {
            neighbours[direction] = NULL;
        }
    }

    real_cpu dx = cell_i->discretization.x;
    real_cpu dy = cell_i->discretization.y;
    real_cpu dz = cell_i->discretization.z;

    // =================================================================
    // FACE NEIGHBORS - Direct diffusion terms (∇·(σ∇u) main diagonal terms)
    // =================================================================

    if(neighbours[RIGHT]) {
        real_cpu sigma_xx = harmonic_mean_2(cell_i->sigma.x, neighbours[RIGHT]->sigma.x);
        add_single_entry(cell_i, neighbours[RIGHT], -sigma_xx * dy * dz / dx);
    }

    if(neighbours[LEFT]) {
        real_cpu sigma_xx =harmonic_mean_2(cell_i->sigma.x, neighbours[LEFT]->sigma.x);
        add_single_entry(cell_i, neighbours[LEFT], -sigma_xx * dy * dz / dx);
    }

    if(neighbours[TOP]) {
        real_cpu sigma_yy = harmonic_mean_2(cell_i->sigma.y, neighbours[TOP]->sigma.y);
        add_single_entry(cell_i, neighbours[TOP], -sigma_yy * dx * dz / dy);
    }

    if(neighbours[DOWN]) {
        real_cpu sigma_yy = harmonic_mean_2(cell_i->sigma.y, neighbours[DOWN]->sigma.y);
        add_single_entry(cell_i, neighbours[DOWN], -sigma_yy * dx * dz / dy);
    }

    if(neighbours[FRONT]) {
        real_cpu sigma_zz = harmonic_mean_2(cell_i->sigma.z, neighbours[FRONT]->sigma.z);
        add_single_entry(cell_i, neighbours[FRONT], -sigma_zz * dx * dy / dz);
    }

    if(neighbours[BACK]) {
        real_cpu sigma_zz = harmonic_mean_2(cell_i->sigma.z, neighbours[BACK]->sigma.z);
        add_single_entry(cell_i, neighbours[BACK], -sigma_zz * dx * dy / dz);
    }

    if(neighbours[TOP_RIGHT] && neighbours[TOP] && neighbours[RIGHT]) {
        // Harmonic mean of σxy over the 4 corner points
        real_cpu sigma_xy_avg = harmonic_mean_4(cell_i->sigma.xy,
                                               neighbours[RIGHT]->sigma.xy,
                                               neighbours[TOP]->sigma.xy,
                                               neighbours[TOP_RIGHT]->sigma.xy);

        // Correct coefficient for XY cross-derivative on structured mesh
        real_cpu coeff = sigma_xy_avg * dz / (dx * dy);

        // Apply the 4-point cross-derivative stencil
        add_single_entry(cell_i, neighbours[TOP_RIGHT], coeff);
        add_single_entry(cell_i, neighbours[TOP], -coeff);
        add_single_entry(cell_i, neighbours[RIGHT], -coeff);
        // The fourth point (cell_i itself) contributes to diagonal: +coeff
        cell_i->elements[0].value += coeff;
    }

    if(neighbours[TOP_LEFT] && neighbours[TOP] && neighbours[LEFT]) {
        real_cpu sigma_xy_avg = harmonic_mean_4(cell_i->sigma.xy,
                                               neighbours[LEFT]->sigma.xy,
                                               neighbours[TOP]->sigma.xy,
                                               neighbours[TOP_LEFT]->sigma.xy);
        real_cpu coeff = sigma_xy_avg * dz / (dx * dy);

        add_single_entry(cell_i, neighbours[TOP_LEFT], -coeff);
        add_single_entry(cell_i, neighbours[TOP], coeff);
        add_single_entry(cell_i, neighbours[LEFT], coeff);
        cell_i->elements[0].value -= coeff;
    }

    if(neighbours[DOWN_RIGHT] && neighbours[DOWN] && neighbours[RIGHT]) {
        real_cpu sigma_xy_avg = harmonic_mean_4(cell_i->sigma.xy,
                                               neighbours[RIGHT]->sigma.xy,
                                               neighbours[DOWN]->sigma.xy,
                                               neighbours[DOWN_RIGHT]->sigma.xy);
        real_cpu coeff = sigma_xy_avg * dz / (dx * dy);

        add_single_entry(cell_i, neighbours[DOWN_RIGHT], -coeff);
        add_single_entry(cell_i, neighbours[DOWN], coeff);
        add_single_entry(cell_i, neighbours[RIGHT], coeff);
        cell_i->elements[0].value -= coeff;
    }

    if(neighbours[DOWN_LEFT] && neighbours[DOWN] && neighbours[LEFT]) {
        real_cpu sigma_xy_avg = harmonic_mean_4(cell_i->sigma.xy,
                                               neighbours[LEFT]->sigma.xy,
                                               neighbours[DOWN]->sigma.xy,
                                               neighbours[DOWN_LEFT]->sigma.xy);
        real_cpu coeff = sigma_xy_avg * dz / (dx * dy);

        add_single_entry(cell_i, neighbours[DOWN_LEFT], coeff);
        add_single_entry(cell_i, neighbours[DOWN], -coeff);
        add_single_entry(cell_i, neighbours[LEFT], -coeff);
        cell_i->elements[0].value += coeff;
    }

    // XZ cross-derivative
    if(neighbours[FRONT_RIGHT] && neighbours[FRONT] && neighbours[RIGHT]) {
        real_cpu sigma_xz_avg = harmonic_mean_4(cell_i->sigma.xz,
                                               neighbours[RIGHT]->sigma.xz,
                                               neighbours[FRONT]->sigma.xz,
                                               neighbours[FRONT_RIGHT]->sigma.xz);
        real_cpu coeff = sigma_xz_avg * dy / (dx * dz);

        add_single_entry(cell_i, neighbours[FRONT_RIGHT], coeff);
        add_single_entry(cell_i, neighbours[FRONT], -coeff);
        add_single_entry(cell_i, neighbours[RIGHT], -coeff);
        cell_i->elements[0].value += coeff;
    }

    if(neighbours[FRONT_LEFT] && neighbours[FRONT] && neighbours[LEFT]) {
        real_cpu sigma_xz_avg = harmonic_mean_4(cell_i->sigma.xz,
                                               neighbours[LEFT]->sigma.xz,
                                               neighbours[FRONT]->sigma.xz,
                                               neighbours[FRONT_LEFT]->sigma.xz);
        real_cpu coeff = sigma_xz_avg * dy / (dx * dz);

        add_single_entry(cell_i, neighbours[FRONT_LEFT], -coeff);
        add_single_entry(cell_i, neighbours[FRONT], coeff);
        add_single_entry(cell_i, neighbours[LEFT], coeff);
        cell_i->elements[0].value -= coeff;
    }

    if(neighbours[BACK_RIGHT] && neighbours[BACK] && neighbours[RIGHT]) {
        real_cpu sigma_xz_avg = harmonic_mean_4(cell_i->sigma.xz,
                                               neighbours[RIGHT]->sigma.xz,
                                               neighbours[BACK]->sigma.xz,
                                               neighbours[BACK_RIGHT]->sigma.xz);
        real_cpu coeff = sigma_xz_avg * dy / (dx * dz);

        add_single_entry(cell_i, neighbours[BACK_RIGHT], -coeff);
        add_single_entry(cell_i, neighbours[BACK], coeff);
        add_single_entry(cell_i, neighbours[RIGHT], coeff);
        cell_i->elements[0].value -= coeff;
    }

    if(neighbours[BACK_LEFT] && neighbours[BACK] && neighbours[LEFT]) {
        real_cpu sigma_xz_avg = harmonic_mean_4(cell_i->sigma.xz,
                                               neighbours[LEFT]->sigma.xz,
                                               neighbours[BACK]->sigma.xz,
                                               neighbours[BACK_LEFT]->sigma.xz);
        real_cpu coeff = sigma_xz_avg * dy / (dx * dz);

        add_single_entry(cell_i, neighbours[BACK_LEFT], coeff);
        add_single_entry(cell_i, neighbours[BACK], -coeff);
        add_single_entry(cell_i, neighbours[LEFT], -coeff);
        cell_i->elements[0].value += coeff;
    }

    // YZ cross-derivative
    if(neighbours[FRONT_TOP] && neighbours[FRONT] && neighbours[TOP]) {
        real_cpu sigma_yz_avg = harmonic_mean_4(cell_i->sigma.yz,
                                               neighbours[TOP]->sigma.yz,
                                               neighbours[FRONT]->sigma.yz,
                                               neighbours[FRONT_TOP]->sigma.yz);
        real_cpu coeff = sigma_yz_avg * dx / (dy * dz);

        add_single_entry(cell_i, neighbours[FRONT_TOP], coeff);
        add_single_entry(cell_i, neighbours[FRONT], -coeff);
        add_single_entry(cell_i, neighbours[TOP], -coeff);
        cell_i->elements[0].value += coeff;
    }

    if(neighbours[FRONT_DOWN] && neighbours[FRONT] && neighbours[DOWN]) {
        real_cpu sigma_yz_avg = harmonic_mean_4(cell_i->sigma.yz,
                                               neighbours[DOWN]->sigma.yz,
                                               neighbours[FRONT]->sigma.yz,
                                               neighbours[FRONT_DOWN]->sigma.yz);
        real_cpu coeff = sigma_yz_avg * dx / (dy * dz);

        add_single_entry(cell_i, neighbours[FRONT_DOWN], -coeff);
        add_single_entry(cell_i, neighbours[FRONT], coeff);
        add_single_entry(cell_i, neighbours[DOWN], coeff);
        cell_i->elements[0].value -= coeff;
    }

    if(neighbours[BACK_TOP] && neighbours[BACK] && neighbours[TOP]) {
        real_cpu sigma_yz_avg = harmonic_mean_4(cell_i->sigma.yz,
                                               neighbours[TOP]->sigma.yz,
                                               neighbours[BACK]->sigma.yz,
                                               neighbours[BACK_TOP]->sigma.yz);
        real_cpu coeff = sigma_yz_avg * dx / (dy * dz);

        add_single_entry(cell_i, neighbours[BACK_TOP], -coeff);
        add_single_entry(cell_i, neighbours[BACK], coeff);
        add_single_entry(cell_i, neighbours[TOP], coeff);
        cell_i->elements[0].value -= coeff;
    }

    if(neighbours[BACK_DOWN] && neighbours[BACK] && neighbours[DOWN]) {
        real_cpu sigma_yz_avg = harmonic_mean_4(cell_i->sigma.yz,
                                               neighbours[DOWN]->sigma.yz,
                                               neighbours[BACK]->sigma.yz,
                                               neighbours[BACK_DOWN]->sigma.yz);
        real_cpu coeff = sigma_yz_avg * dx / (dy * dz);

        add_single_entry(cell_i, neighbours[BACK_DOWN], coeff);
        add_single_entry(cell_i, neighbours[BACK], -coeff);
        add_single_entry(cell_i, neighbours[DOWN], -coeff);
        cell_i->elements[0].value += coeff;
    }
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
