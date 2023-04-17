//
// Created by sachetto on 04/10/17.
//
#include "../3dparty/stb_ds.h"
#include "../config/linear_system_solver_config.h"
#include "../config_helpers/config_helpers.h"
#include "../libraries_common/common_data_structures.h"
#include "../logger/logger.h"

//TODO: remove these global variables
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
        #pragma message ("gpu linear system solver in is using file gpu_solvers_cublas_10.c" )
        #include "gpu_solvers_cublas_10.c"
    #elif CUBLAS_VER_MAJOR == 11
        #pragma message ("gpu linear system solver in is using file gpu_solvers_cublas_11.c" )
        #include "gpu_solvers_cublas_11.c"
    #else
        #pragma message ("gpu linear system solver in is using file gpu_solvers_cublas_12.c" )
        #include "gpu_solvers_cublas_12.c"
    #endif
#endif //COMPILE_CUDA

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


