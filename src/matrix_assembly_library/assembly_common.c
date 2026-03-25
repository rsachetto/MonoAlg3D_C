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

// -----------------------------------------------------------------------------
// Diffusion / obstruction helpers
// -----------------------------------------------------------------------------
// A cell is considered "diffusive" if it is active AND has non-negligible
// diagonal diffusion. This allows two equivalent ways to represent obstacles:
//   (1) cell->active = false
//   (2) keep active but set sigma.{x,y,z} = 0 (and optionally cross terms too)
//
// If you choose (2), this prevents tiny "eps" leakage across obstacle interfaces
// when computing face averages.
static inline bool is_diffusive_cell(const struct cell_node *c) {
    if(!c || !c->active) return false;

    const real_cpu SIGMA_EPS = 1e-18;
    return (fabs(c->sigma.x) + fabs(c->sigma.y) + fabs(c->sigma.z)) > SIGMA_EPS;
}

// Harmonic mean for strictly positive diagonal diffusion coefficients
// (used for sigma_xx on x-faces, sigma_yy on y-faces, sigma_zz on z-faces).
// Returns 0 if either side is ~0 -> enforces no-flux for obstacles represented
// by zero diffusion.
static inline real_cpu harmonic_mean_2_pos(real_cpu a, real_cpu b) {
    const real_cpu eps = 1e-18;
    if(a <= eps || b <= eps) return 0.0;
    return (2.0 * a * b) / (a + b);
}

// Cross terms (sigma_xy, sigma_xz, sigma_yz) may be negative depending on tensor
// orientation. Harmonic means and eps-clamping can distort the sign.
// Use a simple arithmetic mean of the four corner values.
static inline real_cpu arithmetic_mean_4(real_cpu a, real_cpu b, real_cpu c, real_cpu d) {
    return 0.25 * (a + b + c + d);
}



// -----------------------------------------------------------------------------
// Anisotropic diffusion on regular cubes (finite volume, Neumann/no-flux)
// -----------------------------------------------------------------------------
//
// IMPORTANT DESIGN CHOICES (robust with heterogeneity + random obstructions):
// 1) We discretize the operator in *finite-volume flux form* for the full tensor:
//       ∫_V ∇·(σ∇u) dV  =  -∑_faces F_face
//    with diffusive flux (out of the cell)
//       F_face = -A_face * (σ ∇u)_face · n_face
//
// 2) Neumann-only (your case): if a face has no valid neighbor (domain boundary
//    OR obstacle), we simply SKIP the flux on that face -> no-flux.
//
// 3) Obstructions can be represented either by:
//      (a) cell->active = false
//      (b) keeping active but setting σ diagonal to ~0
//    We treat both as non-diffusive via is_diffusive_cell().
//
// 4) Averaging on a face:
//    - Normal component (σ_xx on x-faces, σ_yy on y-faces, σ_zz on z-faces):
//        use 2-point harmonic mean (interface-limited, robust for jumps)
//    - Off-diagonal components on that face (σ_xy, σ_xz, σ_yz):
//        use 2-point arithmetic mean (preserves sign; harmonic is unsafe)
//
// 5) Tangential derivatives on a face:
//    We approximate (∂u/∂tangent)_face = 0.5*( (∂u/∂tangent)_P + (∂u/∂tangent)_N )
//    where each cell-centered derivative uses:
//      - centered difference if both neighbors exist
//      - otherwise a one-sided difference (still annihilates constants)
//      - otherwise 0 (if both sides missing)
//
// This yields a 27-point stencil in the interior, and gracefully reduces near
// obstacles/boundaries WITHOUT special cross-term stencils or diagonal hacks.

static inline struct cell_node *diffusive_neighbour(struct cell_node *c, enum transition_direction dir) {
    struct cell_node *n = get_cell_neighbour_with_same_refinement_level(c, dir);
    return is_diffusive_cell(n) ? n : NULL;
}

// Add a stencil contribution into ROW.
static inline void add_term(struct cell_node *row, struct cell_node *col, real_cpu value) {
    if(!row || !col || fabs(value) < 1e-16) return;
    if(col == row) return;          // diagonal is implied by add_single_entry's diag -= a_ij
    add_single_entry(row, col, value);
}

// 2-point arithmetic mean (safe for off-diagonal tensor entries which may be signed).
static inline real_cpu arithmetic_mean_2(real_cpu a, real_cpu b) {
    return 0.5 * (a + b);
}

// -----------------------------------------------------------------------------
// Cell-centered derivative helpers: add contributions of (∂u/∂x)_C, (∂u/∂y)_C, (∂u/∂z)_C
// into the ROW equation of 'row' scaled by 'scale'.
// -----------------------------------------------------------------------------

static inline void accumulate_dudx_at_cell(struct cell_node *row, struct cell_node *C, real_cpu dx, real_cpu scale) {

    struct cell_node *R = diffusive_neighbour(C, RIGHT);
    struct cell_node *L = diffusive_neighbour(C, LEFT);

    if(R && L) {
        // centered: (u_R - u_L)/(2dx)
        add_term(row, R,  scale * ( 1.0 / (2.0*dx)));
        add_term(row, L,  scale * (-1.0 / (2.0*dx)));
    }
    else if(R) {
        // forward: (u_R - u_C)/dx
        add_term(row, R,  scale * ( 1.0 / dx));
        add_term(row, C,  scale * (-1.0 / dx)); // if C==row, diagonal is implied
    }
    else if(L) {
        // backward: (u_C - u_L)/dx
        add_term(row, C,  scale * ( 1.0 / dx));
        add_term(row, L,  scale * (-1.0 / dx));
    }
    // else: no information -> 0
}

static inline void accumulate_dudy_at_cell(struct cell_node *row, struct cell_node *C, real_cpu dy, real_cpu scale) {

    struct cell_node *T = diffusive_neighbour(C, TOP);
    struct cell_node *D = diffusive_neighbour(C, DOWN);

    if(T && D) {
        // centered: (u_T - u_D)/(2dy)
        add_term(row, T,  scale * ( 1.0 / (2.0*dy)));
        add_term(row, D,  scale * (-1.0 / (2.0*dy)));
    }
    else if(T) {
        // forward: (u_T - u_C)/dy
        add_term(row, T,  scale * ( 1.0 / dy));
        add_term(row, C,  scale * (-1.0 / dy));
    }
    else if(D) {
        // backward: (u_C - u_D)/dy
        add_term(row, C,  scale * ( 1.0 / dy));
        add_term(row, D,  scale * (-1.0 / dy));
    }
}

static inline void accumulate_dudz_at_cell(struct cell_node *row, struct cell_node *C, real_cpu dz, real_cpu scale) {

    struct cell_node *F = diffusive_neighbour(C, FRONT);
    struct cell_node *B = diffusive_neighbour(C, BACK);

    if(F && B) {
        // centered: (u_F - u_B)/(2dz)
        add_term(row, F,  scale * ( 1.0 / (2.0*dz)));
        add_term(row, B,  scale * (-1.0 / (2.0*dz)));
    }
    else if(F) {
        // forward: (u_F - u_C)/dz
        add_term(row, F,  scale * ( 1.0 / dz));
        add_term(row, C,  scale * (-1.0 / dz));
    }
    else if(B) {
        // backward: (u_C - u_B)/dz
        add_term(row, C,  scale * ( 1.0 / dz));
        add_term(row, B,  scale * (-1.0 / dz));
    }
}

// -----------------------------------------------------------------------------
// Face flux assembly for a given row cell P
// -----------------------------------------------------------------------------

static inline void assemble_x_face_flux(struct cell_node *P, struct cell_node *N, bool is_right_face) {

    // x-face area
    const real_cpu dx = P->discretization.x;
    const real_cpu dy = P->discretization.y;
    const real_cpu dz = P->discretization.z;
    const real_cpu A  = dy * dz;

    // Face tensor row for x-normal face: [σxx, σxy, σxz]
    const real_cpu s_xx = harmonic_mean_2_pos(P->sigma.x,  N->sigma.x);
    const real_cpu s_xy = arithmetic_mean_2(P->sigma.xy, N->sigma.xy);
    const real_cpu s_xz = arithmetic_mean_2(P->sigma.xz, N->sigma.xz);

    // --------------------
    // Normal part: -A*σxx*(∂u/∂x)_face
    // We use (u_N - u_P)/dx on the right face, and (u_P - u_N)/dx on the left face.
    // In BOTH cases the resulting coefficient for u_N in the row of P is:
    //     -A * σxx / dx
    // which matches the isotropic FV part of your previous code.
    // --------------------
    add_term(P, N, -s_xx * A / dx);

    // --------------------
    // Cross parts: -A * n_x * [ σxy * (∂u/∂y)_face + σxz * (∂u/∂z)_face ]
    // with n_x = +1 (right face) or -1 (left face).
    // So multiplier for tangential derivatives is:
    //     mult = -A*n_x = (-A) on right, (+A) on left
    // --------------------
    const real_cpu mult = is_right_face ? (-A) : (+A);

    // (∂u/∂y)_face ≈ 0.5*( (∂u/∂y)_P + (∂u/∂y)_N )
    if(fabs(s_xy) > 0.0) {
        const real_cpu scale = 0.5 * mult * s_xy;
        accumulate_dudy_at_cell(P, P, dy, scale);
        accumulate_dudy_at_cell(P, N, dy, scale);
    }

    // (∂u/∂z)_face ≈ 0.5*( (∂u/∂z)_P + (∂u/∂z)_N )
    if(fabs(s_xz) > 0.0) {
        const real_cpu scale = 0.5 * mult * s_xz;
        accumulate_dudz_at_cell(P, P, dz, scale);
        accumulate_dudz_at_cell(P, N, dz, scale);
    }
}

static inline void assemble_y_face_flux(struct cell_node *P, struct cell_node *N, bool is_top_face) {

    const real_cpu dx = P->discretization.x;
    const real_cpu dy = P->discretization.y;
    const real_cpu dz = P->discretization.z;
    const real_cpu A  = dx * dz;

    // y-face tensor row: [σyx, σyy, σyz] with σyx=σxy (symmetric tensor)
    const real_cpu s_yy = harmonic_mean_2_pos(P->sigma.y,  N->sigma.y);
    const real_cpu s_xy = arithmetic_mean_2(P->sigma.xy, N->sigma.xy);
    const real_cpu s_yz = arithmetic_mean_2(P->sigma.yz, N->sigma.yz);

    // Normal part: -A*σyy*(∂u/∂y)_face  -> coefficient to u_N is -A*σyy/dy
    add_term(P, N, -s_yy * A / dy);

    // Tangential parts sign:
    // mult = -A*n_y = (-A) on top face (n_y=+1), (+A) on bottom face (n_y=-1)
    const real_cpu mult = is_top_face ? (-A) : (+A);

    // σyx * (∂u/∂x)_face
    if(fabs(s_xy) > 0.0) {
        const real_cpu scale = 0.5 * mult * s_xy;
        accumulate_dudx_at_cell(P, P, dx, scale);
        accumulate_dudx_at_cell(P, N, dx, scale);
    }

    // σyz * (∂u/∂z)_face
    if(fabs(s_yz) > 0.0) {
        const real_cpu scale = 0.5 * mult * s_yz;
        accumulate_dudz_at_cell(P, P, dz, scale);
        accumulate_dudz_at_cell(P, N, dz, scale);
    }
}

static inline void assemble_z_face_flux(struct cell_node *P, struct cell_node *N, bool is_front_face) {

    const real_cpu dx = P->discretization.x;
    const real_cpu dy = P->discretization.y;
    const real_cpu dz = P->discretization.z;
    const real_cpu A  = dx * dy;

    // z-face tensor row: [σzx, σzy, σzz] with σzx=σxz, σzy=σyz (symmetric tensor)
    const real_cpu s_zz = harmonic_mean_2_pos(P->sigma.z,  N->sigma.z);
    const real_cpu s_xz = arithmetic_mean_2(P->sigma.xz, N->sigma.xz);
    const real_cpu s_yz = arithmetic_mean_2(P->sigma.yz, N->sigma.yz);

    // Normal part: -A*σzz*(∂u/∂z)_face -> coefficient to u_N is -A*σzz/dz
    add_term(P, N, -s_zz * A / dz);

    // Tangential parts sign:
    // mult = -A*n_z = (-A) on front face (n_z=+1), (+A) on back face (n_z=-1)
    const real_cpu mult = is_front_face ? (-A) : (+A);

    // σzx * (∂u/∂x)_face
    if(fabs(s_xz) > 0.0) {
        const real_cpu scale = 0.5 * mult * s_xz;
        accumulate_dudx_at_cell(P, P, dx, scale);
        accumulate_dudx_at_cell(P, N, dx, scale);
    }

    // σzy * (∂u/∂y)_face
    if(fabs(s_yz) > 0.0) {
        const real_cpu scale = 0.5 * mult * s_yz;
        accumulate_dudy_at_cell(P, P, dy, scale);
        accumulate_dudy_at_cell(P, N, dy, scale);
    }
}

static void fill_discretization_matrix_elements_aniso(struct cell_node *cell_i) {

    if(!is_diffusive_cell(cell_i)) {
        return;
    }

    // For Neumann/no-flux BC and for internal obstacles:
    // - If the neighbor is NULL => missing (boundary) or obstacle => no-flux => skip.
    struct cell_node *R = diffusive_neighbour(cell_i, RIGHT);
    struct cell_node *L = diffusive_neighbour(cell_i, LEFT);
    struct cell_node *T = diffusive_neighbour(cell_i, TOP);
    struct cell_node *D = diffusive_neighbour(cell_i, DOWN);
    struct cell_node *F = diffusive_neighbour(cell_i, FRONT);
    struct cell_node *B = diffusive_neighbour(cell_i, BACK);

    if(R) assemble_x_face_flux(cell_i, R, true);
    if(L) assemble_x_face_flux(cell_i, L, false);

    if(T) assemble_y_face_flux(cell_i, T, true);
    if(D) assemble_y_face_flux(cell_i, D, false);

    if(F) assemble_z_face_flux(cell_i, F, true);
    if(B) assemble_z_face_flux(cell_i, B, false);
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
