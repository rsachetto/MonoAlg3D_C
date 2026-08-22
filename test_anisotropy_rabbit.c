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

bool verify_row_sum_zero_rabbit(struct grid *g, double tol, const char *case_name) {
    int n = g->num_active_cells;
    double max_row_sum = 0.0;
    int failed_count = 0;

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
        if(fabs(spatial_row_sum) > tol) {
            failed_count++;
        }
    }
    printf("[%s]\n", case_name);
    printf("   -> Soma Maxima das Linhas em toda a malha |sum_j A_ij^espacial|: %e\n", max_row_sum);
    printf("   -> Celulas que violaram conservacao: %d / %d\n", failed_count, n);
    return max_row_sum < tol;
}

int main() {
    printf("=================================================================\n");
    printf(" TESTE DE FIBRAS NO CORACAO DE COELHO (MESHES/RABHEART.ALG)      \n");
    printf("=================================================================\n");

    const char *mesh_file = "meshes/rabheart.alg";
    uint32_t num_volumes = 470197;
    real_cpu start_discretization = 250.0;

    struct monodomain_solver solver = {.beta = 0.14, .cm = 1.0, .dt = 0.02};
    struct grid *g = new_grid();
    g->adaptive = false;

    printf("--> Lendo malha anatômica: %s ...\n", mesh_file);
    uint32_t num_loaded = set_custom_mesh_from_file(g, (char*)mesh_file, num_volumes, start_discretization, 0, NULL);
    if(num_loaded == 0) {
        printf("ERRO ao carregar %s!\n", mesh_file);
        return 1;
    }
    order_grid_cells(g);
    printf("--> Total de volumes ativos carregados: %d\n", g->num_active_cells);

    double min_z = 1e9, max_z = -1e9;
    for(int i = 0; i < g->num_active_cells; i++) {
        struct cell_node *c = g->active_cells[i];
        if(c->center.z < min_z) min_z = c->center.z;
        if(c->center.z > max_z) max_z = c->center.z;
    }
    double height_z = max_z - min_z;

    double sigma_l = 0.001334;
    double sigma_t = 0.000176;

    // -------------------------------------------------------------
    // CASO A: Fibra 3D Uniforme Rotacionada f = (1/sqrt(3), 1/sqrt(3), 1/sqrt(3))
    // -------------------------------------------------------------
    double inv_sq3 = 1.0 / sqrt(3.0);
    for(int i = 0; i < g->num_active_cells; i++) {
        struct cell_node *c = g->active_cells[i];
        arrfree(c->elements); c->elements = NULL; // reset elements

        double fx = inv_sq3, fy = inv_sq3, fz = inv_sq3;
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
    bool test_A = verify_row_sum_zero_rabbit(g, 1e-12, "CASO A: Fibra 3D Uniforme f=(1/sqrt(3), 1/sqrt(3), 1/sqrt(3))");

    // -------------------------------------------------------------
    // CASO B: Campo de Fibras Helicoidal Variável f(x,y,z) com theta(z) de -60° a +60°
    // -------------------------------------------------------------
    for(int i = 0; i < g->num_active_cells; i++) {
        struct cell_node *c = g->active_cells[i];
        arrfree(c->elements); c->elements = NULL; // reset elements

        // Rotacao transmural de -60 deg a +60 deg ao longo da altura z do coracao
        double norm_z = (c->center.z - min_z) / height_z; // [0, 1]
        double theta = (-60.0 + norm_z * 120.0) * (M_PI / 180.0); // [-pi/3, +pi/3]

        double fx = cos(theta);
        double fy = sin(theta);
        double fz = 0.0;

        c->sigma.x = (sigma_l - sigma_t)*fx*fx + sigma_t;
        c->sigma.y = (sigma_l - sigma_t)*fy*fy + sigma_t;
        c->sigma.z = sigma_t;

        c->sigma.xy = (sigma_l - sigma_t)*fx*fy;
        c->sigma.xz = 0.0;
        c->sigma.yz = 0.0;
    }

    initialize_diagonal_elements(&solver, g);
    for(int i = 0; i < g->num_active_cells; i++) {
        fill_discretization_matrix_elements_aniso(g->active_cells[i]);
    }
    bool test_B = verify_row_sum_zero_rabbit(g, 1e-12, "CASO B: Fibras Helicoidais Variaveis no Espaco f(x,y,z)");

    printf("\n=================================================================\n");
    printf("  RESUMO DA VALIDACAO NO CORACAO DE COELHO (470k CELULAS)        \n");
    printf("=================================================================\n");
    printf("  Caso A (Fibra 3D Uniforme):             %s\n", test_A ? "PASSED (sum A_ij = 0)" : "FAILED");
    printf("  Caso B (Fibras Helicoidais Variaveis): %s\n", test_B ? "PASSED (sum A_ij = 0)" : "FAILED");
    printf("=================================================================\n");

    clean_and_free_grid(g);
    return 0;
}
