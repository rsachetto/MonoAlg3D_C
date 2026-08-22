#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#include "src/alg/grid/grid.h"
#include "src/monodomain/monodomain_solver.h"
#include "src/utils/file_utils.h"
#include "src/utils/utils.h"
#include "src/logger/logger.h"
#include "src/domains_library/domain_helpers.h"
#include "src/matrix_assembly_library/assembly_common.c"

typedef struct {
    int n;
    double **data;
} DenseMatrix;

DenseMatrix create_matrix(int n) {
    DenseMatrix M;
    M.n = n;
    M.data = (double**) malloc(n * sizeof(double*));
    for(int i = 0; i < n; i++) {
        M.data[i] = (double*) calloc(n, sizeof(double));
    }
    return M;
}

void free_matrix(DenseMatrix M) {
    for(int i = 0; i < M.n; i++) {
        free(M.data[i]);
    }
    free(M.data);
}

// Extracts pure spatial operator (A_spatial) by subtracting the mass term ALPHA
DenseMatrix extract_spatial_matrix(struct grid *g) {
    int n = g->num_active_cells;
    DenseMatrix A = create_matrix(n);

    for(int i = 0; i < n; i++) {
        struct cell_node *c = g->active_cells[i];
        struct element *els = c->elements;
        size_t num_els = arrlen(els);

        for(size_t k = 0; k < num_els; k++) {
            int col = els[k].column;
            if(k == 0) {
                double alpha = ALPHA(0.14, 1.0, 0.02, c->discretization.x, c->discretization.y, c->discretization.z);
                A.data[i][i] = els[0].value - alpha;
            } else {
                A.data[i][col] = els[k].value;
            }
        }
    }
    return A;
}

// FIXED: Correctly computes min/max bounds per axis (x, y, z)
bool verify_symmetry(struct grid *g, DenseMatrix A, double tol) {
    double max_diff_interior = 0.0;
    double max_diff_global = 0.0;
    int n = A.n;
    double h = g->active_cells[0]->discretization.x;

    double min_x = 1e9, max_x = -1e9;
    double min_y = 1e9, max_y = -1e9;
    double min_z = 1e9, max_z = -1e9;

    for(int i = 0; i < n; i++) {
        struct cell_node *c = g->active_cells[i];
        if(c->center.x < min_x) min_x = c->center.x;
        if(c->center.x > max_x) max_x = c->center.x;
        if(c->center.y < min_y) min_y = c->center.y;
        if(c->center.y > max_y) max_y = c->center.y;
        if(c->center.z < min_z) min_z = c->center.z;
        if(c->center.z > max_z) max_z = c->center.z;
    }

    for(int i = 0; i < n; i++) {
        struct cell_node *ci = g->active_cells[i];
        bool is_interior_i = (ci->center.x > min_x + 0.5*h && ci->center.x < max_x - 0.5*h &&
                              ci->center.y > min_y + 0.5*h && ci->center.y < max_y - 0.5*h &&
                              ci->center.z > min_z + 0.5*h && ci->center.z < max_z - 0.5*h);

        for(int j = i + 1; j < n; j++) {
            struct cell_node *cj = g->active_cells[j];
            bool is_interior_j = (cj->center.x > min_x + 0.5*h && cj->center.x < max_x - 0.5*h &&
                                  cj->center.y > min_y + 0.5*h && cj->center.y < max_y - 0.5*h &&
                                  cj->center.z > min_z + 0.5*h && cj->center.z < max_z - 0.5*h);

            double diff = fabs(A.data[i][j] - A.data[j][i]);
            if(diff > max_diff_global) max_diff_global = diff;
            if(is_interior_i && is_interior_j) {
                if(diff > max_diff_interior) max_diff_interior = diff;
            }
        }
    }
    printf("   -> Assimetria Maxima no Interior |A_ij - A_ji|: %e\n", max_diff_interior);
    printf("   -> Assimetria Maxima Global (com Bordas) |A_ij - A_ji|: %e\n", max_diff_global);
    return max_diff_interior < tol;
}

bool verify_row_sum_zero(struct grid *g, double tol) {
    int n = g->num_active_cells;
    double max_row_sum = 0.0;

    for(int i = 0; i < n; i++) {
        struct cell_node *c = g->active_cells[i];
        size_t num_els = arrlen(c->elements);
        
        double row_sum = 0.0;
        for(size_t k = 0; k < num_els; k++) {
            row_sum += c->elements[k].value;
        }
        double alpha = ALPHA(0.14, 1.0, 0.02, c->discretization.x, c->discretization.y, c->discretization.z);
        double spatial_row_sum = row_sum - alpha;

        if(fabs(spatial_row_sum) > max_row_sum) {
            max_row_sum = fabs(spatial_row_sum);
        }
    }
    printf("   -> Soma Maxima das Linhas |sum_j A_ij^espacial|: %e\n", max_row_sum);
    return max_row_sum < tol;
}

double test_convergence_3d(int N, double L) {
    double h = L / N;
    
    struct monodomain_solver solver = {.beta = 0.14, .cm = 1.0, .dt = 0.02};
    struct grid *g = new_grid();
    g->adaptive = false;

    set_cuboid_domain_mesh(g, h, h, h, L, L, L);
    order_grid_cells(g);

    // FIXED: Full 3D fiber orientation f = (1/sqrt(3), 1/sqrt(3), 1/sqrt(3))
    // Generates non-zero Dxy, Dxz, Dyz!
    double inv_sq3 = 1.0 / sqrt(3.0);
    double fx = inv_sq3, fy = inv_sq3, fz = inv_sq3;
    double sigma_l = 0.001334;
    double sigma_t = 0.000176;

    for(int i = 0; i < g->num_active_cells; i++) {
        struct cell_node *c = g->active_cells[i];
        c->sigma.x = (sigma_l - sigma_t)*fx*fx + sigma_t;
        c->sigma.y = (sigma_l - sigma_t)*fy*fy + sigma_t;
        c->sigma.z = (sigma_l - sigma_t)*fz*fz + sigma_t;

        c->sigma.xy = (sigma_l - sigma_t)*fx*fy;
        c->sigma.xz = (sigma_l - sigma_t)*fx*fz;
        c->sigma.yz = (sigma_l - sigma_t)*fy*fz;
    }

    initialize_diagonal_elements(&solver, g);

    for(int i = 0; i < g->num_active_cells; i++) {
        fill_discretization_matrix_elements_aniso(g->active_cells[i]);
    }

    DenseMatrix A = extract_spatial_matrix(g);
    
    printf("\n=== VERIFICACAO NUMERICA 3D PARA MALHA N=%d (h=%.2f um) ===\n", N, h);
    bool sym = verify_symmetry(g, A, 1e-12);
    printf("1. Teste de Simetria no Interior: %s\n", sym ? "PASSED (A_ij == A_ji)" : "FAILED");

    bool row_sum = verify_row_sum_zero(g, 1e-12);
    printf("2. Teste de Conservacao (Campo Constante): %s\n", row_sum ? "PASSED (sum A_ij = 0)" : "FAILED");

    // FIXED: Full 3D manufactured solution u(x,y,z) = sin(kx*x) * sin(ky*y) * cos(kz*z)
    double kx = 2.0 * M_PI / L;
    double ky = 2.0 * M_PI / L;
    double kz = 2.0 * M_PI / L;
    
    double max_err = 0.0;
    int interior_count = 0;

    double min_x = 1e9, max_x = -1e9;
    double min_y = 1e9, max_y = -1e9;
    double min_z = 1e9, max_z = -1e9;

    for(int i = 0; i < g->num_active_cells; i++) {
        struct cell_node *c = g->active_cells[i];
        if(c->center.x < min_x) min_x = c->center.x;
        if(c->center.x > max_x) max_x = c->center.x;
        if(c->center.y < min_y) min_y = c->center.y;
        if(c->center.y > max_y) max_y = c->center.y;
        if(c->center.z < min_z) min_z = c->center.z;
        if(c->center.z > max_z) max_z = c->center.z;
    }

    for(int i = 0; i < g->num_active_cells; i++) {
        struct cell_node *c = g->active_cells[i];
        
        // Skip boundary cells for interior truncation error evaluation
        if(c->center.x <= min_x + 0.5*h || c->center.x >= max_x - 0.5*h ||
           c->center.y <= min_y + 0.5*h || c->center.y >= max_y - 0.5*h ||
           c->center.z <= min_z + 0.5*h || c->center.z >= max_z - 0.5*h) {
            continue;
        }

        double Au_spatial = 0.0;
        size_t max_el = arrlen(c->elements);
        for(size_t el = 0; el < max_el; el++) {
            int col_id = c->elements[el].column;
            struct cell_node *col_cell = g->active_cells[col_id];
            double u_val = sin(kx * col_cell->center.x) * sin(ky * col_cell->center.y) * cos(kz * col_cell->center.z);
            
            double weight = c->elements[el].value;
            if(el == 0) {
                double alpha = ALPHA(0.14, 1.0, 0.02, h, h, h);
                weight -= alpha;
            }
            Au_spatial += weight * u_val;
        }
        
        double num_div_flux = - Au_spatial / (h * h * h);

        double x = c->center.x;
        double y = c->center.y;
        double z = c->center.z;

        // Exact analytical derivatives of u(x,y,z) = sin(kx*x)*sin(ky*y)*cos(kz*z)
        double u_xx = -kx*kx * sin(kx*x) * sin(ky*y) * cos(kz*z);
        double u_yy = -ky*ky * sin(kx*x) * sin(ky*y) * cos(kz*z);
        double u_zz = -kz*kz * sin(kx*x) * sin(ky*y) * cos(kz*z);

        double u_xy =  kx*ky * cos(kx*x) * cos(ky*y) * cos(kz*z);
        double u_xz = -kx*kz * cos(kx*x) * sin(ky*y) * sin(kz*z);
        double u_yz = -ky*kz * sin(kx*x) * cos(ky*y) * sin(kz*z);

        double exact_div_flux = c->sigma.x  * u_xx + c->sigma.y  * u_yy + c->sigma.z  * u_zz
                              + 2.0 * c->sigma.xy * u_xy + 2.0 * c->sigma.xz * u_xz + 2.0 * c->sigma.yz * u_yz;

        double err = fabs(num_div_flux - exact_div_flux);
        if(err > max_err) max_err = err;
        interior_count++;
    }

    printf("3. Erro Maximo L_inf no Interior 3D (%d celulas): %e\n", interior_count, max_err);

    free_matrix(A);
    clean_and_free_grid(g);

    return max_err;
}

int main() {
    printf("=================================================================\n");
    printf(" SUITE DE VERIFICACAO NUMERICA COMPLETA 3D (FVM ANISOTROPICO)    \n");
    printf("=================================================================\n");

    double L = 2000.0; // Dominio 2000 um

    double err10 = test_convergence_3d(10, L); // h = 200 um
    double err20 = test_convergence_3d(20, L); // h = 100 um
    double err40 = test_convergence_3d(40, L); // h = 50 um

    double p1 = log2(err10 / err20);
    double p2 = log2(err20 / err40);

    printf("\n=================================================================\n");
    printf("  RESUMO DA ANALISE DE CONVERGENCIA E VERIFICACAO NUMERICA       \n");
    printf("=================================================================\n");
    printf("  Erro (N=10, h=200 um): %e\n", err10);
    printf("  Erro (N=20, h=100 um): %e -> Ordem p_1 = %.2f\n", err20, p1);
    printf("  Erro (N=40, h= 50 um): %e -> Ordem p_2 = %.2f\n", err40, p2);
    printf("=================================================================\n");

    return 0;
}
