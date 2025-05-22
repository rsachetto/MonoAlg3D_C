//
// Created by sachetto on 04/10/17.
//
#include "../3dparty/stb_ds.h"
#include "../config/linear_system_solver_config.h"
#include "../config_helpers/config_helpers.h"
#include "../libraries_common/common_data_structures.h"
#include "../logger/logger.h"

// TODO: remove these global variables
bool jacobi_initialized = false;
bool bcg_initialized = false;
static bool use_preconditioner = false;
static int max_its = 200;
static real_cpu tol = 1e-16;

#ifdef COMPILE_CUDA
#include "../gpu_utils/gpu_utils.h"
#include <cublas_v2.h>
#include <cusparse_v2.h>

#if CUBLAS_VER_MAJOR < 11
#pragma message("gpu linear system solver in is using file gpu_solvers_cublas_10.c")
#include "gpu_solvers_cublas_10.c"
#elif CUBLAS_VER_MAJOR == 11
#pragma message("gpu linear system solver in is using file gpu_solvers_cublas_11.c")
#include "gpu_solvers_cublas_11.c"
#else
#pragma message("gpu linear system solver in is using file gpu_solvers_cublas_12.c")
#include "gpu_solvers_cublas_12.c"
#endif

#ifdef AMGX
#include "amgx_c.h"
#endif

#elif defined(COMPILE_SYCL)
#pragma message("gpu linear system solver in is using file gpu_solvers_sycl.c")
// #include "gpu_solvers_sycl.c"
#endif // COMPILE_CUDA

INIT_LINEAR_SYSTEM(init_cpu_conjugate_gradient) {
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, tol, config, "tolerance");
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(use_preconditioner, config, "use_preconditioner");
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, max_its, config, "max_iterations");
}

END_LINEAR_SYSTEM(end_cpu_conjugate_gradient) {
}

SOLVE_LINEAR_SYSTEM(cpu_conjugate_gradient) {

    real_cpu rTr, pTAp, alpha, beta, precision = tol, rTz, r1Tz1;

    *error = 1.0;
    *number_of_iterations = 1;

    //__________________________________________________________________________
    // Computes int_vector A*x, residue r = b - Ax, scalar rTr = r^T * r and
    // sets initial search direction p.

    rTr = 0.0;
    rTz = 0.0;

    uint32_t i;

    OMP(parallel for reduction(+:rTr,rTz))
    for(i = 0; i < num_active_cells; i++) {

        if(CG_INFO(active_cells[i]) == NULL) {
            INITIALIZE_CONJUGATE_GRADIENT_INFO(active_cells[i]);
        }

        struct element *cell_elements = active_cells[i]->elements;
        active_cells[i]->Ax = 0.0;
        size_t max_el = arrlen(cell_elements);

        for(size_t el = 0; el < max_el; el++) {
            active_cells[i]->Ax += cell_elements[el].value * cell_elements[el].cell->v;
        }

        CG_R(active_cells[i]) = active_cells[i]->b - active_cells[i]->Ax;
        if(use_preconditioner) {
            real_cpu value = cell_elements[0].value;
            if(value == 0.0)
                value = 1.0;
            CG_Z(active_cells[i]) = (1.0 / value) * CG_R(active_cells[i]); // preconditioner
            rTz += CG_R(active_cells[i]) * CG_Z(active_cells[i]);
            CG_P(active_cells[i]) = CG_Z(active_cells[i]);
        } else {
            CG_P(active_cells[i]) = CG_R(active_cells[i]);
        }

        real_cpu r = CG_R(active_cells[i]);

        rTr += r * r;
    }

    *error = rTr;

    //__________________________________________________________________________
    // Conjugate gradient iterations.
    if(*error >= precision) {
        real_cpu r1Tr1;
        while(*number_of_iterations < max_its) {
            //__________________________________________________________________
            // Computes Ap and pTAp. Uses Ax to store Ap.
            pTAp = 0.0;

            OMP(parallel for reduction(+ : pTAp))
            for(i = 0; i < num_active_cells; i++) {

                active_cells[i]->Ax = 0.0;
                struct element *cell_elements = active_cells[i]->elements;

                size_t max_el = arrlen(cell_elements);
                for(size_t el = 0; el < max_el; el++) {
                    active_cells[i]->Ax += cell_elements[el].value * CG_P(cell_elements[el].cell);
                }

                pTAp += CG_P(active_cells[i]) * active_cells[i]->Ax;
            }

            //__________________________________________________________________
            // Computes alpha.
            if(use_preconditioner) {
                alpha = rTz / pTAp;
            } else {
                alpha = rTr / pTAp;
            }
            //__________________________________________________________________

            r1Tr1 = 0.0;
            r1Tz1 = 0.0;

            // Computes new value of solution: u = u + alpha*p.
            OMP(parallel for reduction (+:r1Tr1,r1Tz1))
            for(i = 0; i < num_active_cells; i++) {

                active_cells[i]->v += alpha * CG_P(active_cells[i]);

                CG_R(active_cells[i]) -= alpha * active_cells[i]->Ax;

                real_cpu r = CG_R(active_cells[i]);

                if(use_preconditioner) {
                    real_cpu value = active_cells[i]->elements[0].value;
                    if(value == 0.0)
                        value = 1.0;
                    CG_Z(active_cells[i]) = (1.0 / value) * r;
                    r1Tz1 += CG_Z(active_cells[i]) * r;
                }
                r1Tr1 += r * r;
            }
            //__________________________________________________________________
            // Computes beta.
            if(use_preconditioner) {
                beta = r1Tz1 / rTz;
            } else {
                beta = r1Tr1 / rTr;
            }

            *error = r1Tr1;

            *number_of_iterations = *number_of_iterations + 1;
            if(*error <= precision) {
                break;
            }
            //__________________________________________________________________
            // Computes int_vector p1 = r1 + beta*p and uses it to upgrade p.
            OMP(parallel for)
            for(i = 0; i < num_active_cells; i++) {
                if(use_preconditioner) {
                    CG_P1(active_cells[i]) = CG_Z(active_cells[i]) + beta * CG_P(active_cells[i]);
                } else {
                    CG_P1(active_cells[i]) = CG_R(active_cells[i]) + beta * CG_P(active_cells[i]);
                }
                CG_P(active_cells[i]) = CG_P1(active_cells[i]);
            }

            rTz = r1Tz1;
            rTr = r1Tr1;
        }

    } // end of conjugate gradient iterations.

} // end conjugateGradient() function.

SOLVE_LINEAR_SYSTEM(conjugate_gradient) {

    bool gpu = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(gpu, config, "use_gpu");

    if(gpu) {
#ifdef COMPILE_CUDA
        gpu_conjugate_gradient(time_info, config, the_grid, num_active_cells, active_cells, number_of_iterations, error);
#else
        log_warn("Cuda runtime not found in this system. Fallbacking to CPU solver!!\n");
        cpu_conjugate_gradient(time_info, config, the_grid, num_active_cells, active_cells, number_of_iterations, error);
#endif
    } else {
        cpu_conjugate_gradient(time_info, config, the_grid, num_active_cells, active_cells, number_of_iterations, error);
    }
}

INIT_LINEAR_SYSTEM(init_conjugate_gradient) {
    bool gpu = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(gpu, config, "use_gpu");

    if(gpu) {
#ifdef COMPILE_CUDA
        init_gpu_conjugate_gradient(config, the_grid, is_purkinje);
#else
        init_cpu_conjugate_gradient(config, the_grid, is_purkinje);
#endif
    } else {
        init_cpu_conjugate_gradient(config, the_grid, is_purkinje);
    }
}

END_LINEAR_SYSTEM(end_conjugate_gradient) {

    bool gpu = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(gpu, config, "use_gpu");

    if(gpu) {
#ifdef COMPILE_CUDA
        end_gpu_conjugate_gradient(config);
#else
        end_cpu_conjugate_gradient(config);
#endif
    } else {
        end_cpu_conjugate_gradient(config);
    }
}

// Berg's code
SOLVE_LINEAR_SYSTEM(jacobi) {

    if(!jacobi_initialized) {
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, tol, config, "tolerance");
        max_its = 500;
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, max_its, config, "max_iterations");
        jacobi_initialized = true;
    }

    real_cpu precision = tol;

    *error = 1.0;
    *number_of_iterations = 1;

    struct element element;

    if(*error >= precision) {
        //__________________________________________________________________________
        // Jacobi iterations.
        while(*number_of_iterations < max_its) {

            OMP(parallel for private (element))
            for(uint32_t i = 0; i < num_active_cells; i++) {

                real_cpu sigma = 0.0;

                if(JACOBI_INFO(active_cells[i]) == NULL) {
                    INITIALIZE_JACOBI_INFO(active_cells[i]);
                }

                struct element *cell_elements = active_cells[i]->elements;

                size_t max_el = arrlen(cell_elements);

                // Do not take the diagonal element
                for(size_t el = 1; el < max_el; el++) {
                    element = cell_elements[el];
                    sigma += element.value * element.cell->v;
                }

                real_cpu value = cell_elements[0].value;
                JACOBI_X_AUX(active_cells[i]) = (1.0 / value) * (active_cells[i]->b - sigma);
            }
            real_cpu residue = 0.0;
            real_cpu sum;

            OMP(parallel for private (element,sum) reduction (+:residue))
            for(uint32_t i = 0; i < num_active_cells; i++) {

                struct element *cell_elements = active_cells[i]->elements;

                size_t max_el = arrlen(cell_elements);

                // Do not take the diagonal element
                sum = 0.0;
                for(size_t el = 0; el < max_el; el++) {
                    element = cell_elements[el];
                    sum += element.value * JACOBI_X_AUX(element.cell);
                }

                active_cells[i]->v = JACOBI_X_AUX(active_cells[i]);
                residue += pow(active_cells[i]->b - sum, 2);
            }
            // The error is norm of the residue
            residue = sqrt(residue);
            *error = residue;

            *number_of_iterations = *number_of_iterations + 1;
            if(*error <= precision)
                break;
        }
    }
}

//// Berg's code
SOLVE_LINEAR_SYSTEM(biconjugate_gradient) {

    if(!bcg_initialized) {
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, tol, config, "tolerance");

        char *preconditioner_char = NULL;
        GET_PARAMETER_STRING_VALUE_OR_USE_DEFAULT(preconditioner_char, config, "use_preconditioner");
        if(preconditioner_char != NULL) {
            use_preconditioner = ((strcmp(preconditioner_char, "yes") == 0) || (strcmp(preconditioner_char, "true") == 0));
        }

        max_its = 100;
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, max_its, config, "max_iterations");
        bcg_initialized = true;
    }

    real_cpu rTr, pTAp, alpha, beta, precision = tol, rTz, r1Tz1;

    *error = 1.0;
    *number_of_iterations = 1;

    struct element element;
    //__________________________________________________________________________
    // Zero all entries on the int_vector x*A
    // And initialize the second guess vector x_aux
    OMP(parallel for)
    for(uint32_t i = 0; i < num_active_cells; i++) {

        if(BCG_INFO(active_cells[i]) == NULL) {
            INITIALIZE_BICONJUGATE_GRADIENT_INFO(active_cells[i]);
        }

        BCG_XA(active_cells[i]) = 0.0;
        BCG_X_AUX(active_cells[i]) = active_cells[i]->v;
    }

    //__________________________________________________________________________
    // Computes int_vector A*x, x*A
    // xA must be fully calculated to start doing anything over the r_aux vector
    OMP(parallel for private (element))
    for(uint32_t i = 0; i < num_active_cells; i++) {
        struct element *cell_elements = active_cells[i]->elements;
        active_cells[i]->Ax = 0.0;

        size_t max_el = arrlen(cell_elements);

        for(size_t el = 0; el < max_el; el++) {
            element = cell_elements[el];
            uint32_t col = element.column;
            active_cells[i]->Ax += element.value * element.cell->v;

            OMP(critical)
            BCG_XA(active_cells[col]) += element.value * BCG_X_AUX(active_cells[i]);
        }
    }

    rTr = 0.0;
    rTz = 0.0;

    //__________________________________________________________________________
    // Computes residues r, r_aux
    // scalar rTr = r^T * r_aux and
    // sets initial search directions p and p_aux.
    OMP(parallel for private (element) reduction(+:rTr,rTz))
    for(uint32_t i = 0; i < num_active_cells; i++) {
        struct element *cell_elements = active_cells[i]->elements;

        BCG_R(active_cells[i]) = active_cells[i]->b - active_cells[i]->Ax;
        BCG_R_AUX(active_cells[i]) = active_cells[i]->b - BCG_XA(active_cells[i]);

        if(use_preconditioner) {
            real_cpu value = cell_elements[0].value;
            if(value == 0.0)
                value = 1.0;
            BCG_Z(active_cells[i]) = (1.0 / value) * BCG_R(active_cells[i]); // preconditioner
            BCG_Z_AUX(active_cells[i]) = (1.0 / value) * BCG_R_AUX(active_cells[i]);
            rTz += BCG_R_AUX(active_cells[i]) * BCG_Z(active_cells[i]);
            BCG_P(active_cells[i]) = BCG_Z(active_cells[i]);
            BCG_P_AUX(active_cells[i]) = BCG_Z_AUX(active_cells[i]);
        } else {
            BCG_P(active_cells[i]) = BCG_R(active_cells[i]);
            BCG_P_AUX(active_cells[i]) = BCG_R_AUX(active_cells[i]);
        }
        rTr += BCG_R_AUX(active_cells[i]) * BCG_R(active_cells[i]);
    }

    *error = rTr;

    //__________________________________________________________________________
    // Biconjugate gradient iterations.
    if(*error >= precision) {
        real_cpu r1Tr1;
        while(*number_of_iterations < max_its) {
            //__________________________________________________________________
            // Computes Ap, pA and pTAp. Uses Ax to store Ap and xA to store pA
            pTAp = 0.0;

            OMP(parallel for)
            for(uint32_t i = 0; i < num_active_cells; i++)
                BCG_XA(active_cells[i]) = 0.0;

            OMP(parallel for private(element) reduction(+ : pTAp))
            for(uint32_t i = 0; i < num_active_cells; i++) {
                active_cells[i]->Ax = 0.0;
                struct element *cell_elements = active_cells[i]->elements;

                size_t max_el = arrlen(cell_elements);
                for(size_t el = 0; el < max_el; el++) {
                    element = cell_elements[el];
                    uint32_t col = element.column;
                    active_cells[i]->Ax += element.value * BCG_P(element.cell);

                    OMP(critical)
                    BCG_XA(active_cells[col]) += element.value * BCG_P_AUX(active_cells[i]);
                }

                pTAp += BCG_P_AUX(active_cells[i]) * active_cells[i]->Ax;
            }

            //__________________________________________________________________
            // Computes alpha.
            if(use_preconditioner) {
                alpha = rTz / pTAp;
            } else {
                alpha = rTr / pTAp;
            }
            //__________________________________________________________________

            r1Tr1 = 0.0;
            r1Tz1 = 0.0;

            // Computes new value of solution: u = u + alpha*p.
            //                                 u_aux = u_aux + alpha*p_aux
            OMP(parallel for reduction (+:r1Tr1,r1Tz1))
            for(uint32_t i = 0; i < num_active_cells; i++) {

                active_cells[i]->v += alpha * BCG_P(active_cells[i]);
                BCG_X_AUX(active_cells[i]) += alpha * BCG_P_AUX(active_cells[i]);

                BCG_R(active_cells[i]) -= alpha * active_cells[i]->Ax;
                BCG_R_AUX(active_cells[i]) -= alpha * BCG_XA(active_cells[i]);

                if(use_preconditioner) {
                    real_cpu value = active_cells[i]->elements[0].value;
                    if(value == 0.0)
                        value = 1.0;
                    BCG_Z(active_cells[i]) = (1.0 / value) * BCG_R(active_cells[i]);
                    BCG_Z_AUX(active_cells[i]) = (1.0 / value) * BCG_R_AUX(active_cells[i]);
                    r1Tz1 += BCG_Z(active_cells[i]) * BCG_R_AUX(active_cells[i]);
                }

                r1Tr1 += BCG_R(active_cells[i]) * BCG_R_AUX(active_cells[i]);
            }
            //__________________________________________________________________
            // Computes beta.
            if(use_preconditioner) {
                beta = r1Tz1 / rTz;
            } else {
                beta = r1Tr1 / rTr;
            }

            *error = r1Tr1;
            *number_of_iterations = *number_of_iterations + 1;
            if(*error <= precision) {
                break;
            }

            //__________________________________________________________________
            // Computes int_vector p1 = r1 + beta*p and uses it to upgrade p.
            OMP(parallel for)
            for(uint32_t i = 0; i < num_active_cells; i++) {
                if(use_preconditioner) {
                    BCG_P1(active_cells[i]) = BCG_Z(active_cells[i]) + beta * BCG_P(active_cells[i]);
                    BCG_P1_AUX(active_cells[i]) = BCG_Z_AUX(active_cells[i]) + beta * BCG_P_AUX(active_cells[i]);
                } else {
                    BCG_P1(active_cells[i]) = BCG_R(active_cells[i]) + beta * BCG_P(active_cells[i]);
                    BCG_P1_AUX(active_cells[i]) = BCG_R_AUX(active_cells[i]) + beta * BCG_P_AUX(active_cells[i]);
                }
                BCG_P(active_cells[i]) = BCG_P1(active_cells[i]);
                BCG_P_AUX(active_cells[i]) = BCG_P1_AUX(active_cells[i]);
            }

            rTz = r1Tz1;
            rTr = r1Tr1;
        }

    } // end of biconjugate gradient iterations.

} // end biconjugateGradient() function.

#ifdef AMGX

/* print callback (could be customized) */
void print_callback(const char *msg, int length) {
    // printf("%s", msg);
}

char *get_amgx_config(const char *config_name, int max_iters, float tol) {

    char cfg[2048];

#define STR_EQUALS(a, b) strcmp(a, b) == 0

    // Generate JSON string based on the configName using if-else statements
    if(STR_EQUALS(config_name, "AMG_CLASSICAL_PMIS")) {
        sprintf(cfg,
                "{\"config_version\":2,\"determinism_flag\":1,\"solver\":{\"print_grid_stats\":0,\"obtain_timings\":0,\"solver\":\"GMRES\",\"print_solve_"
                "stats\":0,\"preconditioner\":{\"interpolator\":\"D2\",\"solver\":\"AMG\",\"cycle\":\"V\",\"smoother\":{\"relaxation_factor\":1,\"scope\":"
                "\"jacobi\",\"solver\":\"JACOBI_L1\"},\"presweeps\":2,\"postsweeps\":2,\"selector\":\"PMIS\",\"coarsest_sweeps\":2,\"coarse_solver\":"
                "\"NOSOLVER\",\"max_iters\":1,\"max_row_sum\":0.9,\"min_coarse_rows\":2,\"scope\":\"amg_solver\",\"max_levels\":24,\"print_grid_stats\":1,"
                "\"aggressive_levels\":1,\"interp_max_elements\":4},\"max_iters\":%d,\"store_res_history\":1,\"monitor_residual\":1,\"gmres_n_restart\":10,"
                "\"convergence\":\"RELATIVE_INI_CORE\",\"tolerance\":%g,\"norm\":\"L2\"}}",
                max_iters, tol);

    } else if(STR_EQUALS(config_name, "GMRES_AMG_D2")) {
        sprintf(cfg,
                "{\"config_version\":2,\"determinism_flag\":1,\"exception_handling\":1,\"solver\":{\"print_grid_stats\":0,\"solver\":\"GMRES\",\"print_solve_"
                "stats\":0,\"obtain_timings\":0,\"preconditioner\":{\"interpolator\":\"D2\",\"print_grid_stats\":0,\"solver\":\"AMG\",\"smoother\":\"JACOBI_"
                "L1\",\"presweeps\":2,\"selector\":\"PMIS\",\"coarsest_sweeps\":2,\"coarse_solver\":\"NOSOLVER\",\"max_iters\":1,\"interp_max_elements\":4,"
                "\"min_coarse_rows\":2,\"scope\":\"amg_solver\",\"max_levels\":24,\"cycle\":\"V\",\"postsweeps\":2},\"max_iters\":%d,\"store_res_history\":1,"
                "\"monitor_residual\":1,\"gmres_n_restart\":10,\"convergence\":\"ABSOLUTE\",\"tolerance\":%g,\"norm\":\"L2\"}}",
                max_iters, tol);

    } else if(STR_EQUALS(config_name, "PBICGSTAB_NOPREC")) {
        sprintf(cfg,
                "{\"config_version\":2,\"solver\":{\"preconditioner\":{\"scope\":\"amg_solver\",\"solver\":\"NOSOLVER\"},\"use_scalar_norm\":1,\"store_res_"
                "history\":1,\"solver\":\"PBICGSTAB\",\"print_solve_stats\":0,\"obtain_timings\":0,\"max_iters\":%d,\"monitor_residual\":1,\"convergence\":"
                "\"ABSOLUTE\",\"scope\":\"main\",\"tolerance\":%g,\"norm\":\"L2\"}}",
                max_iters, tol);

    } else if(STR_EQUALS(config_name, "AMG_CLASSICAL_AGGRESSIVE_CHEB_L1_TRUNC")) {
        sprintf(cfg,
                "{\"config_version\": 2, \"determinism_flag\": 1, \"solver\": {\"scope\": \"main\", \"solver\": \"PCG\", \"store_res_history\": 1, "
                "\"print_solve_stats\": 0, \"obtain_timings\": 0, \"preconditioner\": {\"print_grid_stats\": 0, \"scope\": \"amg_solver\", \"interpolator\": "
                "\"D2\", \"solver\": \"AMG\", \"max_levels\": 24, \"selector\": \"PMIS\", \"cycle\": \"V\", \"presweeps\": 0, \"postsweeps\": 3, "
                "\"coarsest_sweeps\": 2, \"min_coarse_rows\": 2, \"coarse_solver\": \"NOSOLVER\", \"max_iters\": 1, \"max_row_sum\": 0.9, "
                "\"strength_threshold\": 0.25, \"error_scaling\": 3, \"print_grid_stats\": 1, \"aggressive_levels\": 1, \"interp_max_elements\": 4, "
                "\"smoother\": {\"relaxation_factor\": 0.91, \"scope\": \"jacobi\", \"solver\": \"CHEBYSHEV\", \"max_iters\": 1, \"preconditioner\": "
                "{\"solver\": \"JACOBI_L1\", \"max_iters\": 1}, \"chebyshev_polynomial_order\": 2, \"chebyshev_lambda_estimate_mode\": 2}}, \"max_iters\": %d, "
                "\"monitor_residual\": 1, \"convergence\": \"ABSOLUTE\", \"tolerance\": %g, \"norm\": \"L2\"}}",
                max_iters, tol);
    } else if(STR_EQUALS(config_name, "FGMRES_CLASSICAL_AGGRESSIVE_HMIS")) {
        sprintf(
            cfg,
            "{\"config_version\": 2, \"solver\": {\"print_grid_stats\": 0, \"solver\": \"FGMRES\", \"print_solve_stats\": 0, \"obtain_timings\": 0, "
            "\"preconditioner\": {\"interpolator\": \"D2\", \"solver\": \"AMG\", \"print_grid_stats\": 1, \"aggressive_levels\": 1, \"interp_max_elements\": "
            "4, \"smoother\": {\"relaxation_factor\": 1, \"scope\": \"jacobi\", \"solver\": \"JACOBI_L1\"}, \"presweeps\": 2, \"selector\": \"HMIS\", "
            "\"coarsest_sweeps\": 2, \"coarse_solver\": \"NOSOLVER\", \"max_iters\": 1, \"max_row_sum\": 0.9, \"strength_threshold\": 0.25, "
            "\"min_coarse_rows\": 2, \"scope\": \"amg_solver\", \"max_levels\": 24, \"cycle\": \"V\", \"postsweeps\": 2}, \"max_iters\": %d, "
            "\"store_res_history\": 1, \"monitor_residual\": 1, \"gmres_n_restart\": 100, \"convergence\": \"ABSOLUTE\", \"tolerance\": %g, \"norm\": \"L2\"}}",
            max_iters, tol);
    } else if(STR_EQUALS(config_name, "FGMRES_AGGREGATION")) {
        sprintf(cfg,
                "{\"config_version\": 2, \"solver\": {\"preconditioner\": {\"error_scaling\": 0, \"print_grid_stats\": 0, \"max_uncolored_percentage\": 0.05, "
                "\"algorithm\": \"AGGREGATION\", \"solver\": \"AMG\", \"smoother\": \"MULTICOLOR_DILU\", \"presweeps\": 0, \"selector\": \"SIZE_2\", "
                "\"coarse_solver\": \"DENSE_LU_SOLVER\", \"max_iters\": 1, \"postsweeps\": 3, \"min_coarse_rows\": 32, \"relaxation_factor\": 0.75, \"scope\": "
                "\"amg\", \"max_levels\": 100, \"matrix_coloring_scheme\": \"PARALLEL_GREEDY\", \"cycle\": \"V\"}, \"use_scalar_norm\": 1, \"solver\": "
                "\"FGMRES\", \"print_solve_stats\": 0, \"obtain_timings\": 0, \"max_iters\": %d, \"monitor_residual\": 0, \"gmres_n_restart\": 10, "
                "\"convergence\": \"ABSOLUTE\", \"scope\": \"main\", \"store_res_history\": 1, \"monitor_residual\": 1, \"tolerance\": %g, \"norm\": \"L2\"}}",
                max_iters, tol);
    } else if(STR_EQUALS(config_name, "PCG_NOPREC")) {
        sprintf(
            cfg,
            "{\"config_version\":2,\"solver\":{\"store_res_history\": 1, "
            "\"preconditioner\":{\"scope\":\"amg\",\"solver\":\"NOSOLVER\"},\"use_scalar_norm\":1,\"solver\":\"PCG\",\"print_solve_stats\":0,\"obtain_"
            "timings\":0,\"monitor_residual\":1,\"convergence\":\"RELATIVE_INI_CORE\",\"scope\":\"main\",\"max_iters\":%d, \"tolerance\":%e,\"norm\":\"L2\"}}",
            max_iters, tol);
    } else if(STR_EQUALS(config_name, "PCG_CLASSICAL_F_JACOBI")) {
        sprintf(cfg,
                "{\"config_version\":2,\"solver\":{\"preconditioner\":{\"print_grid_stats\":0,\"print_vis_data\":0,\"solver\":\"AMG\",\"smoother\":{\"scope\":"
                "\"jacobi\",\"solver\":\"BLOCK_JACOBI\",\"monitor_residual\":0,\"print_solve_stats\":0},\"print_solve_stats\":0,\"presweeps\":1,\"max_iters\":"
                "1,\"interpolator\":\"D2\",\"monitor_residual\":0,\"store_res_history\":0,\"scope\":\"amg\",\"max_levels\":100,\"cycle\":\"F\",\"postsweeps\":"
                "1},\"solver\":\"PCG\",\"print_solve_stats\":0,\"obtain_timings\":0,\"max_iters\":%d,\"store_res_history\": "
                "1,\"monitor_residual\":1,\"convergence\":\"ABSOLUTE\",\"scope\":\"main\",\"tolerance\":%e,\"norm\":\"L2\"}}",
                max_iters, tol);
    } else if(STR_EQUALS(config_name, "PCG_DILU")) {
        sprintf(cfg,
                "{\"config_version\":2,\"solver\":{\"preconditioner\":{\"scope\":\"precond\",\"solver\":\"MULTICOLOR_DILU\"},\"solver\":\"CG\",\"print_solve_"
                "stats\":0,\"obtain_timings\":0,\"max_iters\":%d,\"store_res_history\":1, "
                "\"monitor_residual\":1,\"scope\":\"main\",\"tolerance\":%e,\"norm\":\"L2\"}}",
                max_iters, tol);
    } else {
        log_error("Error: Unknown configuration name.\n");
    }

    return strdup(cfg);
}

struct amgx_persistent_data {
    int N, nz;
    // library handles
    AMGX_config_handle cfg;
    AMGX_resources_handle rsrc;
    AMGX_matrix_handle A;
    AMGX_vector_handle b, x;
    AMGX_solver_handle solver;
    // status handling
    AMGX_SOLVE_STATUS status;
};

INIT_LINEAR_SYSTEM(init_amgx) {

    struct amgx_persistent_data *persistent_data = CALLOC_ONE_TYPE(struct amgx_persistent_data);

    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, tol, config, "tolerance");
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, max_its, config, "max_iterations");

    char *amgx_config = strdup("AMG_CLASSICAL_PMIS");
    GET_PARAMETER_STRING_VALUE_OR_USE_DEFAULT(amgx_config, config, "amgx_config");

    /* init */
    AMGX_SAFE_CALL(AMGX_initialize());
    /* system */
    AMGX_SAFE_CALL(AMGX_register_print_callback(&print_callback));
    AMGX_SAFE_CALL(AMGX_install_signal_handler());

    char *tmp = get_amgx_config(amgx_config, max_its, tol);
    printf("%s\n", tmp);
    AMGX_SAFE_CALL(AMGX_config_create(&(persistent_data->cfg), tmp));
    free(tmp);

    AMGX_Mode mode = AMGX_mode_dFFI;

    AMGX_resources_create_simple(&(persistent_data->rsrc), persistent_data->cfg);
    AMGX_matrix_create(&(persistent_data->A), persistent_data->rsrc, mode);
    AMGX_vector_create(&(persistent_data->x), persistent_data->rsrc, mode);
    AMGX_vector_create(&(persistent_data->b), persistent_data->rsrc, mode);
    AMGX_solver_create(&(persistent_data->solver), persistent_data->rsrc, mode, persistent_data->cfg);

    int_array I = NULL, J = NULL;
    f32_array val = NULL;

    uint32_t num_active_cells;
    struct cell_node **active_cells = NULL;

    if(is_purkinje) {
        grid_to_csr(the_grid, &val, &I, &J, true);
        num_active_cells = the_grid->purkinje->num_active_purkinje_cells;
        active_cells = the_grid->purkinje->purkinje_cells;
    } else {
        grid_to_csr(the_grid, &val, &I, &J, false);
        num_active_cells = the_grid->num_active_cells;
        active_cells = the_grid->active_cells;
    }

    int nz = arrlen(val);
    int N = num_active_cells;

    persistent_data->nz = nz;
    persistent_data->N = N;

    AMGX_SAFE_CALL(AMGX_matrix_upload_all(persistent_data->A, N, nz, 1, 1, I, J, val, NULL));

    float *rhs = (float *)malloc(sizeof(float) * N);

    OMP(parallel for)
    for(uint32_t i = 0; i < num_active_cells; i++) {
        rhs[i] = active_cells[i]->b;
    }

    AMGX_SAFE_CALL(AMGX_vector_upload(persistent_data->x, N, 1, rhs));

    // Set up the AMGX solver
    AMGX_SAFE_CALL(AMGX_solver_setup(persistent_data->solver, persistent_data->A));

    config->persistent_data = persistent_data;
    free(rhs);
    arrfree(I);
    arrfree(J);
    arrfree(val);
    free(amgx_config);
}

SOLVE_LINEAR_SYSTEM(amgx) {

    struct amgx_persistent_data *persistent_data = (struct amgx_persistent_data *)config->persistent_data;

    if(!persistent_data) {
        log_error_and_exit(
            "The amgx solver needs to be initialized before being called. Add a init_function in the [linear_system_solver] section of the .ini file\n");
    }

    float *rhs; // Vector B
    rhs = (float *)malloc(sizeof(float) * num_active_cells);

    OMP(parallel for)
    for(uint32_t i = 0; i < num_active_cells; i++) {
        rhs[i] = active_cells[i]->b;
    }

    int N = persistent_data->N;
    AMGX_SAFE_CALL(AMGX_vector_upload(persistent_data->b, N, 1, rhs));

    // Solve the system using AMGX
    AMGX_SAFE_CALL(AMGX_solver_solve(persistent_data->solver, persistent_data->b, persistent_data->x));

    AMGX_SAFE_CALL(AMGX_vector_download(persistent_data->x, rhs));

    int its;
    double e = 0;
    AMGX_SAFE_CALL(AMGX_solver_get_iterations_number(persistent_data->solver, &its)); // alt: AMGX_solver_get_iterations_number(AMGX_solver_handle obj, int *n)
    AMGX_SAFE_CALL(AMGX_solver_get_iteration_residual(persistent_data->solver, its - 1, 0,
                                                      &e)); // AMGX_solver_get_iteration_residual(AMGX_solver_handle obj, int iter,int idx, double *res);

    *error = e;
    *number_of_iterations = (uint32_t)its;

    OMP(parallel for)
    for(uint32_t i = 0; i < num_active_cells; i++) {
        active_cells[i]->v = rhs[i];
    }

    free(rhs);
}

END_LINEAR_SYSTEM(end_amgx) {
    struct amgx_persistent_data *persistent_data = (struct amgx_persistent_data *)config->persistent_data;

    // Cleanup and finalize AMGX
    AMGX_SAFE_CALL(AMGX_solver_destroy(persistent_data->solver));
    AMGX_SAFE_CALL(AMGX_vector_destroy(persistent_data->x));
    AMGX_SAFE_CALL(AMGX_vector_destroy(persistent_data->b));
    AMGX_SAFE_CALL(AMGX_matrix_destroy(persistent_data->A));
    AMGX_SAFE_CALL(AMGX_resources_destroy(persistent_data->rsrc));
    AMGX_SAFE_CALL(AMGX_config_destroy(persistent_data->cfg));
    AMGX_SAFE_CALL(AMGX_finalize());
}

#endif
