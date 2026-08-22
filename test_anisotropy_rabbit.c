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

// =====================================================================
// Fotografia esparsa da matriz montada, usada para comparar duas
// montagens diferentes sobre a MESMA malha (não-adaptativa, mesma
// ordem de células e mesma ordem de vizinhos por linha).
// =====================================================================
typedef struct {
    int n;
    size_t *counts;
    int **cols;
    double **vals;
} MatrixSnapshot;

MatrixSnapshot capture_snapshot(struct grid *g) {
    MatrixSnapshot s;
    s.n = g->num_active_cells;
    s.counts = (size_t*) malloc(s.n * sizeof(size_t));
    s.cols = (int**) malloc(s.n * sizeof(int*));
    s.vals = (double**) malloc(s.n * sizeof(double*));

    for(int i = 0; i < s.n; i++) {
        struct element *els = g->active_cells[i]->elements;
        size_t ne = arrlen(els);
        s.counts[i] = ne;
        s.cols[i] = (int*) malloc(ne * sizeof(int));
        s.vals[i] = (double*) malloc(ne * sizeof(double));
        for(size_t k = 0; k < ne; k++) {
            s.cols[i][k] = els[k].column;
            s.vals[i][k] = els[k].value;
        }
    }
    return s;
}

void free_snapshot(MatrixSnapshot s) {
    for(int i = 0; i < s.n; i++) {
        free(s.cols[i]);
        free(s.vals[i]);
    }
    free(s.cols);
    free(s.vals);
    free(s.counts);
}

// Compara duas fotografias assumindo ordem determinística de vizinhos
// por linha (mesma malha, mesmo código de montagem chamado duas vezes).
bool compare_snapshots(MatrixSnapshot a, MatrixSnapshot b, double tol, const char *case_name) {
    if(a.n != b.n) {
        printf("[%s] ERRO: numero de celulas diferente entre as duas montagens\n", case_name);
        return false;
    }

    double max_diff = 0.0;
    int mismatched_rows = 0;

    for(int i = 0; i < a.n; i++) {
        if(a.counts[i] != b.counts[i]) {
            mismatched_rows++;
            continue;
        }
        for(size_t k = 0; k < a.counts[i]; k++) {
            if(a.cols[i][k] != b.cols[i][k]) {
                mismatched_rows++;
                break;
            }
            double diff = fabs(a.vals[i][k] - b.vals[i][k]);
            if(diff > max_diff) max_diff = diff;
        }
    }

    printf("[%s]\n", case_name);
    printf("   -> Diferenca Maxima entre as duas montagens: %e\n", max_diff);
    printf("   -> Linhas com estrutura (colunas) diferente: %d / %d\n", mismatched_rows, a.n);
    return (max_diff < tol) && (mismatched_rows == 0);
}

// =====================================================================
// Testes sobre a matriz montada (sempre esparsos - nunca materializam
// uma matriz densa n x n)
// =====================================================================

bool verify_row_sum_zero(struct grid *g, double tol, const char *case_name) {
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

        if(fabs(spatial_row_sum) > max_row_sum) max_row_sum = fabs(spatial_row_sum);
        if(fabs(spatial_row_sum) > tol) failed_count++;
    }
    printf("[%s]\n", case_name);
    printf("   -> Soma Maxima das Linhas |sum_j A_ij^espacial|: %e\n", max_row_sum);
    printf("   -> Celulas que violaram conservacao: %d / %d\n", failed_count, n);
    return max_row_sum < tol;
}

// Simetria sem matriz densa: para cada aresta (i -> j) procura a aresta
// reciproca (j -> i) na lista esparsa de j e compara os valores.
bool verify_symmetry_sparse(struct grid *g, double tol, const char *case_name) {
    double max_diff = 0.0;
    int n = g->num_active_cells;

    for(int i = 0; i < n; i++) {
        struct element *els_i = g->active_cells[i]->elements;
        size_t ne_i = arrlen(els_i);

        for(size_t k = 1; k < ne_i; k++) { // pula k=0 (diagonal)
            int j = els_i[k].column;
            struct element *els_j = g->active_cells[j]->elements;
            size_t ne_j = arrlen(els_j);

            double val_ji = 0.0;
            bool found = false;
            for(size_t m = 0; m < ne_j; m++) {
                if(els_j[m].column == i) {
                    val_ji = els_j[m].value;
                    found = true;
                    break;
                }
            }
            double diff = found ? fabs(els_i[k].value - val_ji) : fabs(els_i[k].value);
            if(diff > max_diff) max_diff = diff;
        }
    }
    printf("[%s]\n", case_name);
    printf("   -> Assimetria Maxima |A_ij - A_ji| (toda a malha): %e\n", max_diff);
    return max_diff < tol;
}

// Traco do tensor de condutividade e invariante a rotacao da fibra:
// trace(D) = (sigma_l - sigma_t)*|f|^2 + 3*sigma_t = sigma_l + 2*sigma_t
// para qualquer fibra unitaria. Pega erro de normalizacao de f antes
// mesmo de chegar na montagem.
bool verify_trace_invariant(struct grid *g, double sigma_l, double sigma_t, double tol, const char *case_name) {
    double expected = sigma_l + 2.0 * sigma_t;
    double max_diff = 0.0;
    int n = g->num_active_cells;

    for(int i = 0; i < n; i++) {
        struct cell_node *c = g->active_cells[i];
        double trace = c->sigma.x + c->sigma.y + c->sigma.z;
        double diff = fabs(trace - expected);
        if(diff > max_diff) max_diff = diff;
    }
    printf("[%s]\n", case_name);
    printf("   -> Desvio Maximo do traco do tensor (esperado %e): %e\n", expected, max_diff);
    return max_diff < tol;
}

// Com fibra alinhada a um eixo cartesiano, os termos cruzados do tensor
// sao exatamente zero, entao nenhuma celula deveria ter na sua lista de
// elementos um vizinho deslocado em mais de um eixo simultaneamente
// (stencil de 7 pontos puro, sem "vizinho diagonal").
bool verify_axis_aligned_no_cross_terms(struct grid *g, const char *case_name) {
    int violations = 0;
    int n = g->num_active_cells;

    for(int i = 0; i < n; i++) {
        struct cell_node *ci = g->active_cells[i];
        struct element *els = ci->elements;
        size_t ne = arrlen(els);

        for(size_t k = 1; k < ne; k++) {
            struct cell_node *cj = g->active_cells[els[k].column];
            int diff_axes = (fabs(cj->center.x - ci->center.x) > 1e-6) +
                             (fabs(cj->center.y - ci->center.y) > 1e-6) +
                             (fabs(cj->center.z - ci->center.z) > 1e-6);
            if(diff_axes > 1) violations++;
        }
    }
    printf("[%s]\n", case_name);
    printf("   -> Vizinhos fora-de-eixo com fibra alinhada: %d\n", violations);
    return violations == 0;
}

// Reporta os 'top_k' pares (i,j) com maior assimetria |A_ij - A_ji|,
// junto com o discretization.{x,y,z} de cada celula - permite checar se
// a assimetria esta correlacionada com celulas vizinhas de tamanhos
// diferentes (malha localmente refinada / nao-cubica).
typedef struct { int i, j; double diff, val_ij, val_ji; } SymOffender;

void report_worst_symmetry_offenders(struct grid *g, int top_k, const char *case_name) {
    SymOffender *worst = (SymOffender*) calloc(top_k, sizeof(SymOffender));
    int n = g->num_active_cells;

    for(int i = 0; i < n; i++) {
        struct element *els_i = g->active_cells[i]->elements;
        size_t ne_i = arrlen(els_i);
        for(size_t k = 1; k < ne_i; k++) {
            int j = els_i[k].column;
            if(j < i) continue; // conta cada par uma unica vez

            struct element *els_j = g->active_cells[j]->elements;
            size_t ne_j = arrlen(els_j);
            double val_ji = 0.0;
            bool found = false;
            for(size_t m = 0; m < ne_j; m++) {
                if(els_j[m].column == i) { val_ji = els_j[m].value; found = true; break; }
            }
            double diff = found ? fabs(els_i[k].value - val_ji) : fabs(els_i[k].value);

            if(diff > worst[top_k-1].diff) {
                worst[top_k-1] = (SymOffender){i, j, diff, els_i[k].value, val_ji};
                for(int p = top_k-1; p > 0 && worst[p].diff > worst[p-1].diff; p--) {
                    SymOffender tmp = worst[p]; worst[p] = worst[p-1]; worst[p-1] = tmp;
                }
            }
        }
    }

    printf("[%s] Top %d pares mais assimetricos:\n", case_name, top_k);
    for(int r = 0; r < top_k && worst[r].diff > 0; r++) {
        struct cell_node *ci = g->active_cells[worst[r].i];
        struct cell_node *cj = g->active_cells[worst[r].j];
        int diff_axes = (fabs(cj->center.x - ci->center.x) > 1e-6) +
                         (fabs(cj->center.y - ci->center.y) > 1e-6) +
                         (fabs(cj->center.z - ci->center.z) > 1e-6);
        bool same_size = fabs(ci->discretization.x - cj->discretization.x) < 1e-6 &&
                          fabs(ci->discretization.y - cj->discretization.y) < 1e-6 &&
                          fabs(ci->discretization.z - cj->discretization.z) < 1e-6;
        printf("   #%d i=%d j=%d A_ij=%e A_ji=%e diff=%e eixos_dif=%d mesmo_tamanho=%s\n",
               r+1, worst[r].i, worst[r].j, worst[r].val_ij, worst[r].val_ji,
               worst[r].diff, diff_axes, same_size ? "sim" : "NAO");
        printf("       h_i=(%.2f,%.2f,%.2f)  h_j=(%.2f,%.2f,%.2f)\n",
               ci->discretization.x, ci->discretization.y, ci->discretization.z,
               cj->discretization.x, cj->discretization.y, cj->discretization.z);
    }
    free(worst);
}

// Conta elementos fora da diagonal positivos: um FVM anisotropico bem
// comportado deve ter A_ij <= 0 para vizinhos (fluxo). Anisotropia forte
// + fibra bem tiltada pode quebrar essa propriedade (perda de M-matrix).
// Nao falha o teste sozinho - so da visibilidade.
int count_positive_off_diagonal(struct grid *g, const char *case_name) {
    int count = 0;
    int n = g->num_active_cells;

    for(int i = 0; i < n; i++) {
        struct element *els = g->active_cells[i]->elements;
        size_t ne = arrlen(els);
        for(size_t k = 1; k < ne; k++) {
            if(els[k].value > 1e-12) count++;
        }
    }
    printf("[%s]\n", case_name);
    printf("   -> Elementos fora da diagonal positivos (possivel perda de M-matrix): %d\n", count);
    return count;
}

// =====================================================================
// Preparacao da malha / campos de fibra
// =====================================================================

void reset_elements(struct grid *g) {
    for(int i = 0; i < g->num_active_cells; i++) {
        struct cell_node *c = g->active_cells[i];
        arrfree(c->elements);
        c->elements = NULL;
    }
}

void assemble(struct grid *g, struct monodomain_solver *solver) {
    initialize_diagonal_elements(solver, g);
    for(int i = 0; i < g->num_active_cells; i++) {
        fill_discretization_matrix_elements_aniso(g->active_cells[i]);
    }
}

void set_uniform_fiber(struct grid *g, double fx, double fy, double fz, double sigma_l, double sigma_t) {
    for(int i = 0; i < g->num_active_cells; i++) {
        struct cell_node *c = g->active_cells[i];
        c->sigma.x = (sigma_l - sigma_t)*fx*fx + sigma_t;
        c->sigma.y = (sigma_l - sigma_t)*fy*fy + sigma_t;
        c->sigma.z = (sigma_l - sigma_t)*fz*fz + sigma_t;
        c->sigma.xy = (sigma_l - sigma_t)*fx*fy;
        c->sigma.xz = (sigma_l - sigma_t)*fx*fz;
        c->sigma.yz = (sigma_l - sigma_t)*fy*fz;
    }
}

// Rotacao transmural de -60 a +60 graus ao longo da altura z do coracao,
// fibra sempre no plano xy (fz = 0).
void set_helicoidal_fiber(struct grid *g, double min_z, double height_z, double sigma_l, double sigma_t) {
    for(int i = 0; i < g->num_active_cells; i++) {
        struct cell_node *c = g->active_cells[i];
        double norm_z = (c->center.z - min_z) / height_z; // [0, 1]
        double theta = (-60.0 + norm_z * 120.0) * (M_PI / 180.0); // [-pi/3, +pi/3]
        double fx = cos(theta), fy = sin(theta), fz = 0.0;

        c->sigma.x = (sigma_l - sigma_t)*fx*fx + sigma_t;
        c->sigma.y = (sigma_l - sigma_t)*fy*fy + sigma_t;
        c->sigma.z = sigma_t;
        c->sigma.xy = (sigma_l - sigma_t)*fx*fy;
        c->sigma.xz = 0.0;
        c->sigma.yz = 0.0;
    }
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

    printf("--> Lendo malha anatomica: %s ...\n", mesh_file);
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
    double tol = 1e-12;

    bool test_A, test_B, test_C, test_D, test_E;

    // -----------------------------------------------------------------
    // CASO A: Fibra 3D Uniforme Rotacionada f = (1,1,1)/sqrt(3)
    // -----------------------------------------------------------------
    double inv_sq3 = 1.0 / sqrt(3.0);
    reset_elements(g);
    set_uniform_fiber(g, inv_sq3, inv_sq3, inv_sq3, sigma_l, sigma_t);
    verify_trace_invariant(g, sigma_l, sigma_t, tol, "CASO A - traco do tensor");
    assemble(g, &solver);
    test_A = verify_row_sum_zero(g, tol, "CASO A: Fibra 3D Uniforme f=(1,1,1)/sqrt(3) - conservacao");
    test_A = verify_symmetry_sparse(g, tol, "CASO A: Fibra 3D Uniforme - simetria") && test_A;
    report_worst_symmetry_offenders(g, 15, "CASO A: diagnostico de assimetria");
    count_positive_off_diagonal(g, "CASO A: sinais fora da diagonal");
    MatrixSnapshot snap_A = capture_snapshot(g);

    // -----------------------------------------------------------------
    // CASO D: Invariancia sob f -> -f (mesma fibra 3D, sinal invertido)
    // -----------------------------------------------------------------
    reset_elements(g);
    set_uniform_fiber(g, -inv_sq3, -inv_sq3, -inv_sq3, sigma_l, sigma_t);
    assemble(g, &solver);
    MatrixSnapshot snap_D = capture_snapshot(g);
    test_D = compare_snapshots(snap_A, snap_D, tol, "CASO D: Invariancia f -> -f (deve ser identico ao Caso A)");
    free_snapshot(snap_A);
    free_snapshot(snap_D);

    // -----------------------------------------------------------------
    // CASO B: Campo de Fibras Helicoidal Variavel f(x,y,z), theta(z) de -60 a +60
    // -----------------------------------------------------------------
    reset_elements(g);
    set_helicoidal_fiber(g, min_z, height_z, sigma_l, sigma_t);
    verify_trace_invariant(g, sigma_l, sigma_t, tol, "CASO B - traco do tensor");
    assemble(g, &solver);
    test_B = verify_row_sum_zero(g, tol, "CASO B: Fibras Helicoidais Variaveis - conservacao");
    test_B = verify_symmetry_sparse(g, tol, "CASO B: Fibras Helicoidais Variaveis - simetria") && test_B;
    report_worst_symmetry_offenders(g, 15, "CASO B: diagnostico de assimetria");
    count_positive_off_diagonal(g, "CASO B: sinais fora da diagonal");

    // -----------------------------------------------------------------
    // CASO C: Fibra alinhada a eixo f = (1,0,0), meio anisotropico
    //         (termos cruzados devem ser exatamente zero na montagem)
    // -----------------------------------------------------------------
    reset_elements(g);
    set_uniform_fiber(g, 1.0, 0.0, 0.0, sigma_l, sigma_t);
    assemble(g, &solver);
    test_C = verify_row_sum_zero(g, tol, "CASO C: Fibra axis-aligned f=(1,0,0) - conservacao");
    test_C = verify_symmetry_sparse(g, tol, "CASO C: Fibra axis-aligned - simetria") && test_C;
    test_C = verify_axis_aligned_no_cross_terms(g, "CASO C: ausencia de termos cruzados") && test_C;

    // -----------------------------------------------------------------
    // CASO F (referencia) + CASO E: Invariancia sob meio isotropico.
    // Com sigma_l == sigma_t o tensor colapsa em sigma_t*I para
    // QUALQUER fibra, entao a montagem com o campo helicoidal do Caso B
    // deve ser IDENTICA a montagem com fibra fixa axis-aligned.
    // -----------------------------------------------------------------
    reset_elements(g);
    set_uniform_fiber(g, 1.0, 0.0, 0.0, sigma_t, sigma_t); // sigma_l = sigma_t
    assemble(g, &solver);
    MatrixSnapshot snap_F = capture_snapshot(g);

    reset_elements(g);
    set_helicoidal_fiber(g, min_z, height_z, sigma_t, sigma_t); // sigma_l = sigma_t
    assemble(g, &solver);
    MatrixSnapshot snap_E = capture_snapshot(g);

    test_E = compare_snapshots(snap_F, snap_E, tol,
        "CASO E: Invariancia sob meio isotropico (helicoidal deve ser identico ao axis-aligned)");
    free_snapshot(snap_F);
    free_snapshot(snap_E);

    printf("\n=================================================================\n");
    printf("  RESUMO DA VALIDACAO NO CORACAO DE COELHO (%d CELULAS)         \n", g->num_active_cells);
    printf("=================================================================\n");
    printf("  Caso A (Fibra 3D Uniforme, conservacao+simetria): %s\n", test_A ? "PASSED" : "FAILED");
    printf("  Caso B (Fibras Helicoidais, conservacao+simetria): %s\n", test_B ? "PASSED" : "FAILED");
    printf("  Caso C (Fibra axis-aligned, sem termos cruzados):  %s\n", test_C ? "PASSED" : "FAILED");
    printf("  Caso D (Invariancia f -> -f):                      %s\n", test_D ? "PASSED" : "FAILED");
    printf("  Caso E (Invariancia sob meio isotropico):          %s\n", test_E ? "PASSED" : "FAILED");
    printf("=================================================================\n");

    clean_and_free_grid(g);
    return (test_A && test_B && test_C && test_D && test_E) ? 0 : 1;
}
