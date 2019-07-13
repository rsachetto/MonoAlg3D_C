//
// Created by sachetto on 04/10/17.
//

#include "../config/linear_system_solver_config.h"
#include "../libraries_common/config_helpers.h"
#include "../libraries_common/common_data_structures.h"

#include "../single_file_libraries/stb_ds.h"
#include "../models_library/model_gpu_utils.h"

bool cg_cpu_initialized = false;
bool cg_gpu_initialized = false;
bool jacobi_initialized = false;
bool bcg_initialized = false;
bool use_jacobi;
int max_its = 50;
real_cpu tol = 1e-16;

#include <cusparse_v2.h>
#include <cublas_v2.h>


#ifdef COMPILE_CUDA
SOLVE_LINEAR_SYSTEM(gpu_conjugate_gradient) {

    /* Conjugate gradient without preconditioning.
       ------------------------------------------
       Follows the description by Golub & Van Loan, "Matrix Computations 3rd ed.", Section 10.2.6  */


    if(!cg_gpu_initialized) {
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, tol, config->config_data.config, "tolerance");
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, max_its, config->config_data.config, "max_iterations");
        cg_gpu_initialized = true;
    }

    int M = 0, N = 0, nz = 0;

    static int *I = NULL, *J = NULL;
    static real *val = NULL;

    real a, b, na, r0, r1;

    static int *d_col, *d_row;
    static real *d_val, *d_x;

    real dot;
    static real *d_r, *d_p, *d_Ax;

    int k;
    real alpha, beta, alpham1;

    real *rhs; //Vector B

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node** ac = the_grid->active_cells;

    rhs = (real*) malloc(sizeof(real)*num_active_cells);

    #pragma omp parallel for
    for (int i = 0; i < num_active_cells; i++) {
        rhs[i] = ac[i]->b;
    }

    /* Get handle to the CUBLAS context */
    static cublasHandle_t cublasHandle = 0;
    cublasStatus_t cublasStatus;

    /* Get handle to the CUSPARSE context */
    static cusparseHandle_t cusparseHandle = 0;
    cusparseStatus_t cusparseStatus;

    static cusparseMatDescr_t descr = 0;

    N = M = num_active_cells;

    if(val == NULL) {

        cublasStatus = cublasCreate(&cublasHandle);
        check_cuda_error((cudaError_t)cublasStatus);

        cusparseStatus = cusparseCreate(&cusparseHandle);
        check_cuda_error((cudaError_t)cusparseStatus);

        cusparseStatus = cusparseCreateMatDescr(&descr);
        check_cuda_error((cudaError_t)cusparseStatus);

        cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL);
        cusparseSetMatIndexBase(descr,CUSPARSE_INDEX_BASE_ZERO);


        grid_to_csr(the_grid, &val, &I, &J);

        nz = arrlen(val);

        check_cuda_error(cudaMalloc((void **) &d_col, nz * sizeof(int)));
        check_cuda_error(cudaMalloc((void **) &d_row, (N + 1) * sizeof(int)));
        check_cuda_error(cudaMalloc((void **) &d_val, nz * sizeof(float)));
        check_cuda_error(cudaMalloc((void **) &d_x, N * sizeof(float)));
        check_cuda_error(cudaMalloc((void **) &d_r, N * sizeof(float)));
        check_cuda_error(cudaMalloc((void **) &d_p, N * sizeof(float)));
        check_cuda_error(cudaMalloc((void **) &d_Ax, N * sizeof(float)));

        cudaMemcpy(d_col, J, nz * sizeof(int), cudaMemcpyHostToDevice); //JA
        cudaMemcpy(d_row, I, (N + 1) * sizeof(int), cudaMemcpyHostToDevice); //IA
        cudaMemcpy(d_val, val, nz * sizeof(float), cudaMemcpyHostToDevice); //A
        cudaMemcpy(d_x, rhs, N * sizeof(float), cudaMemcpyHostToDevice); //Result
    }
    else {
        nz = arrlen(val);
    }

    cudaMemcpy(d_r, rhs, N * sizeof(float), cudaMemcpyHostToDevice); //B

    alpha = 1.0;
    alpham1 = -1.0;
    beta = 0.0;
    r0 = 0.;

    cusparseScsrmv(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &alpha, descr, d_val, d_row, d_col, d_x, &beta, d_Ax);

    cublasSaxpy(cublasHandle, N, &alpham1, d_Ax, 1, d_r, 1);
    cublasStatus = cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);

    k = 1;

    while (r1 > tol*tol && k <= max_its)
    {
        if (k > 1)
        {
            b = r1 / r0;
            cublasStatus = cublasSscal(cublasHandle, N, &b, d_p, 1);
            cublasStatus = cublasSaxpy(cublasHandle, N, &alpha, d_r, 1, d_p, 1);
        }
        else
        {
            cublasStatus = cublasScopy(cublasHandle, N, d_r, 1, d_p, 1);
        }

        cusparseScsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &alpha, descr, d_val, d_row, d_col, d_p, &beta, d_Ax);
        cublasStatus = cublasSdot(cublasHandle, N, d_p, 1, d_Ax, 1, &dot);
        a = r1 / dot;

        cublasStatus = cublasSaxpy(cublasHandle, N, &a, d_p, 1, d_x, 1);
        na = -a;
        cublasStatus = cublasSaxpy(cublasHandle, N, &na, d_Ax, 1, d_r, 1);

        r0 = r1;
        cublasStatus = cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
        cudaDeviceSynchronize();
        k++;
    }

    cudaMemcpy(rhs, d_x, N*sizeof(float), cudaMemcpyDeviceToHost);

    *number_of_iterations = k-1;
    *error = sqrt(r1);

    #pragma omp parallel for
    for (int i = 0; i < num_active_cells; i++) {
        ac[i]->v = rhs[i];
    }

    free(rhs);

}
#endif

SOLVE_LINEAR_SYSTEM(conjugate_gradient) {

    if(!cg_cpu_initialized) {
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, tol, config->config_data.config, "tolerance");
        GET_PARAMETER_BINARY_VALUE_OR_USE_DEFAULT(use_jacobi, config->config_data.config, "use_preconditioner");
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, max_its, config->config_data.config, "max_iterations");
        cg_cpu_initialized = true;
    }


    real_cpu  rTr,
            r1Tr1,
            pTAp,
            alpha,
            beta,
            precision = tol,
            rTz,
            r1Tz1;


    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node** ac = the_grid->active_cells;

    *error = 1.0;
    *number_of_iterations = 1;

    //__________________________________________________________________________
    //Computes int_vector A*x, residue r = b - Ax, scalar rTr = r^T * r and
    //sets initial search direction p.

    rTr = 0.0;
    rTz = 0.0;

    struct element element;
    int i;

    #pragma omp parallel for private (element) reduction(+:rTr,rTz)
    for (i = 0; i < num_active_cells; i++) {

        if(CG_INFO(ac[i]) == NULL) {
            INITIALIZE_CONJUGATE_GRADIENT_INFO(ac[i]);
        }

        struct element *cell_elements = ac[i]->elements;
        ac[i]->Ax = 0.0;

        size_t max_el = arrlen(cell_elements);

        for(int el = 0; el < max_el; el++) {
            element = cell_elements[el];
            ac[i]->Ax += element.value * element.cell->v;
        }

        CG_R(ac[i]) = ac[i]->b - ac[i]->Ax;
        if(use_jacobi) {
            real_cpu value = cell_elements[0].value;
            if(value == 0.0) value = 1.0;
            CG_Z(ac[i]) = (1.0/value) * CG_R(ac[i]); // preconditioner
            rTz += CG_R(ac[i]) * CG_Z(ac[i]);
            CG_P(ac[i]) = CG_Z(ac[i]);
        }
        else {
            CG_P(ac[i]) = CG_R(ac[i]);
        }

        real_cpu r = CG_R(ac[i]);

        rTr += r * r;

    }

    *error = rTr;


    //__________________________________________________________________________
    //Conjugate gradient iterations.
    if( *error >= precision ) {
        while( *number_of_iterations < max_its ) {
            //__________________________________________________________________
            // Computes Ap and pTAp. Uses Ax to store Ap.
            pTAp = 0.0;

            #pragma omp parallel for private(element) reduction(+ : pTAp)
            for (i = 0; i < num_active_cells; i++) {

                ac[i]->Ax = 0.0;
                struct element *cell_elements = ac[i]->elements;

                size_t max_el = arrlen(cell_elements);
                for(int el = 0; el < max_el; el++) {
                    element = cell_elements[el];
                    ac[i]->Ax += element.value * CG_P(element.cell);
                }

                pTAp += CG_P(ac[i]) * ac[i]->Ax;
            }

            //__________________________________________________________________
            // Computes alpha.
            if(use_jacobi) {
                alpha = rTz/pTAp;
            }
            else {
                alpha = rTr/pTAp;
            }
            //__________________________________________________________________


            r1Tr1 = 0.0;
            r1Tz1 = 0.0;

            // Computes new value of solution: u = u + alpha*p.
            #pragma omp parallel for reduction (+:r1Tr1,r1Tz1)
            for (i = 0; i < num_active_cells; i++) {
                ac[i]->v += alpha * CG_P(ac[i]);

                CG_R(ac[i]) -= alpha * ac[i]->Ax;

                if(use_jacobi) {
                    real_cpu value = ac[i]->elements[0].value;
                    if(value == 0.0) value = 1.0;
                    CG_Z(ac[i]) = (1.0/value) * CG_R(ac[i]);
                    r1Tz1 += CG_Z(ac[i]) * CG_R(ac[i]);
                }

                real_cpu r = CG_R(ac[i]);

                r1Tr1 += r * r;
            }
            //__________________________________________________________________
            //Computes beta.
            if(use_jacobi) {
                beta = r1Tz1/rTz;
            }
            else {
                beta = r1Tr1/rTr;
            }

            *error = r1Tr1;

            *number_of_iterations = *number_of_iterations + 1;
            if( *error <= precision ) {
                break;
            }
            //__________________________________________________________________
            //Computes int_vector p1 = r1 + beta*p and uses it to upgrade p.
            #pragma omp parallel for
            for (i = 0; i < num_active_cells; i++) {
                if(use_jacobi) {
                    CG_P1(ac[i]) = CG_Z(ac[i]) + beta * CG_P(ac[i]);
                }
                else {
                    CG_P1(ac[i]) = CG_R(ac[i]) + beta * CG_P(ac[i]);
                }
                CG_P(ac[i]) = CG_P1(ac[i]);
            }

            rTz = r1Tz1;
            rTr = r1Tr1;

        }

    }//end of conjugate gradient iterations.

}//end conjugateGradient() function.

// Berg's code
SOLVE_LINEAR_SYSTEM(jacobi) {


    if(!jacobi_initialized) {
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, tol, config->config_data.config, "tolerance");
        max_its = 500;
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, max_its, config->config_data.config, "max_iterations");
        jacobi_initialized = true;
    }


    real_cpu  sigma,
            precision = tol;

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node** ac = the_grid->active_cells;

    *error = 1.0;
    *number_of_iterations = 1;

    struct element element;
    int i;

    if (*error >= precision)
    {
        //__________________________________________________________________________
        //Jacobi iterations.
        while (*number_of_iterations < max_its)
        {
            #pragma omp parallel for private (element,sigma)
            for (i = 0; i < num_active_cells; i++)
            {
                if(JACOBI_INFO(ac[i]) == NULL) {
                    INITIALIZE_JACOBI_INFO(ac[i]);
                }

                struct element *cell_elements = ac[i]->elements;
                sigma = 0.0;

                size_t max_el = arrlen(cell_elements);

                // Do not take the diagonal element
                for(int el = 1; el < max_el; el++)
                {
                    element = cell_elements[el];
                    sigma += element.value * element.cell->v;
                }

                real_cpu value = cell_elements[0].value;
                JACOBI_X_AUX(ac[i]) = (1.0/value)*(ac[i]->b - sigma);
            }
            real_cpu residue = 0.0;
            real_cpu sum;
            #pragma omp parallel for private (element,sum) reduction (+:residue)
            for (i = 0; i < num_active_cells; i++)
            {
                struct element *cell_elements = ac[i]->elements;

                size_t max_el = arrlen(cell_elements);

                // Do not take the diagonal element
                sum = 0.0;
                for(int el = 0; el < max_el; el++)
                {
                    element = cell_elements[el];
                    sum += element.value * JACOBI_X_AUX(element.cell);
                }

                ac[i]->v = JACOBI_X_AUX(ac[i]);
                residue += pow(ac[i]->b - sum,2);
            }
            // The error is norm of the residue
            residue = sqrt(residue);
            *error = residue;

            *number_of_iterations = *number_of_iterations + 1;
            if( *error <= precision )
                break;
        }
    }
}

//// Berg's code
SOLVE_LINEAR_SYSTEM(biconjugate_gradient)
{


    if(!bcg_initialized) {
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, tol, config->config_data.config, "tolerance");

        char *preconditioner_char = NULL;
        GET_PARAMETER_VALUE_CHAR_OR_USE_DEFAULT(preconditioner_char, config->config_data.config, "use_preconditioner");
        if (preconditioner_char != NULL)
        {
            use_jacobi = ((strcmp (preconditioner_char, "yes") == 0) || (strcmp (preconditioner_char, "true") == 0));
        }

        max_its = 100;
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, max_its, config->config_data.config, "max_iterations");
        bcg_initialized = true;
    }


    real_cpu  rTr,
            r1Tr1,
            pTAp,
            alpha,
            beta,
            precision = tol,
            rTz,
            r1Tz1;


    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node** ac = the_grid->active_cells;

    *error = 1.0;
    *number_of_iterations = 1;

    struct element element;
    int i;

    //__________________________________________________________________________
    // Zero all entries on the int_vector x*A
    // And initialize the second guess vector x_aux
    #pragma omp parallel for
    for (i = 0; i < num_active_cells; i++)
    {

        if(BCG_INFO(ac[i]) == NULL) {
            INITIALIZE_BICONJUGATE_GRADIENT_INFO(ac[i]);
        }

        BCG_XA(ac[i]) = 0.0;
        BCG_X_AUX(ac[i]) = ac[i]->v;
    }


    //__________________________________________________________________________
    //Computes int_vector A*x, x*A
    //xA must be fully calculated to start doing anything over the r_aux vector
    #pragma omp parallel for private (element)
    for (i = 0; i < num_active_cells; i++)
    {
        struct element *cell_elements = ac[i]->elements;
        ac[i]->Ax = 0.0;

        size_t max_el = arrlen(cell_elements);

        for(int el = 0; el < max_el; el++)
        {
            element = cell_elements[el];
            uint32_t col = element.column;
            ac[i]->Ax += element.value * element.cell->v;

            #pragma omp critical
            BCG_XA(ac[col]) += element.value * BCG_X_AUX(ac[i]);
        }
    }

    rTr = 0.0;
    rTz = 0.0;

    //__________________________________________________________________________
    //Computes residues r, r_aux
    //scalar rTr = r^T * r_aux and
    //sets initial search directions p and p_aux.
    #pragma omp parallel for private (element) reduction(+:rTr,rTz)
    for (i = 0; i < num_active_cells; i++)
    {
        struct element *cell_elements = ac[i]->elements;

        BCG_R(ac[i]) = ac[i]->b - ac[i]->Ax;
        BCG_R_AUX(ac[i]) = ac[i]->b - BCG_XA(ac[i]);

        if(use_jacobi)
        {
            real_cpu value = cell_elements[0].value;
            if(value == 0.0) value = 1.0;
            BCG_Z(ac[i]) = (1.0/value) * BCG_R(ac[i]); // preconditioner
            BCG_Z_AUX(ac[i]) = (1.0/value) * BCG_R_AUX(ac[i]);
            rTz += BCG_R_AUX(ac[i]) * BCG_Z(ac[i]);
            BCG_P(ac[i]) = BCG_Z(ac[i]);
            BCG_P_AUX(ac[i]) = BCG_Z_AUX(ac[i]);
        }
        else
        {
            BCG_P(ac[i]) = BCG_R(ac[i]);
            BCG_P_AUX(ac[i])= BCG_R_AUX(ac[i]);
        }
        rTr += BCG_R_AUX(ac[i]) * BCG_R(ac[i]);
    }

    *error = rTr;

    //__________________________________________________________________________
    //Biconjugate gradient iterations.
    if( *error >= precision )
    {
        while( *number_of_iterations < max_its )
        {
            //__________________________________________________________________
            // Computes Ap, pA and pTAp. Uses Ax to store Ap and xA to store pA
            pTAp = 0.0;

            #pragma omp parallel for
            for (i = 0; i < num_active_cells; i++)
                BCG_XA(ac[i]) = 0.0;

            #pragma omp parallel for private(element) reduction(+ : pTAp)
            for (i = 0; i < num_active_cells; i++)
            {
                ac[i]->Ax = 0.0;
                struct element *cell_elements = ac[i]->elements;

                size_t max_el = arrlen(cell_elements);
                for(int el = 0; el < max_el; el++)
                {
                    element = cell_elements[el];
                    uint32_t col = element.column;
                    ac[i]->Ax += element.value * BCG_P(element.cell);

                    #pragma omp critical
                    BCG_XA(ac[col]) += element.value * BCG_P_AUX(ac[i]);
                }

                pTAp += BCG_P_AUX(ac[i]) * ac[i]->Ax;
            }

            //__________________________________________________________________
            // Computes alpha.
            if(use_jacobi)
            {
                alpha = rTz/pTAp;
            }
            else
            {
                alpha = rTr/pTAp;
            }
            //__________________________________________________________________

            r1Tr1 = 0.0;
            r1Tz1 = 0.0;

            // Computes new value of solution: u = u + alpha*p.
            //                                 u_aux = u_aux + alpha*p_aux
            #pragma omp parallel for reduction (+:r1Tr1,r1Tz1)
            for (i = 0; i < num_active_cells; i++)
            {
                ac[i]->v += alpha * BCG_P(ac[i]);
                BCG_X_AUX(ac[i]) += alpha * BCG_P_AUX(ac[i]);

                BCG_R(ac[i]) -= alpha * ac[i]->Ax;
                BCG_R_AUX(ac[i]) -= alpha * BCG_XA(ac[i]);

                if(use_jacobi)
                {
                    real_cpu value = ac[i]->elements[0].value;
                    if(value == 0.0) value = 1.0;
                    BCG_Z(ac[i]) = (1.0/value) * BCG_R(ac[i]);
                    BCG_Z_AUX(ac[i]) = (1.0/value) * BCG_R_AUX(ac[i]);
                    r1Tz1 += BCG_Z(ac[i]) * BCG_R_AUX(ac[i]);
                }

                r1Tr1 += BCG_R(ac[i]) * BCG_R_AUX(ac[i]);
            }
            //__________________________________________________________________
            //Computes beta.
            if(use_jacobi)
            {
                beta = r1Tz1/rTz;
            }
            else
            {
                beta = r1Tr1/rTr;
            }

            *error = r1Tr1;
            *number_of_iterations = *number_of_iterations + 1;
            if( *error <= precision )
            {
                break;
            }

            //__________________________________________________________________
            //Computes int_vector p1 = r1 + beta*p and uses it to upgrade p.
            #pragma omp parallel for
            for (i = 0; i < num_active_cells; i++)
            {
                if(use_jacobi)
                {
                    BCG_P1(ac[i]) = BCG_Z(ac[i]) + beta * BCG_P(ac[i]);
                    BCG_P1_AUX(ac[i]) = BCG_Z_AUX(ac[i]) + beta * BCG_P_AUX(ac[i]);
                }
                else
                {
                    BCG_P1(ac[i]) = BCG_R(ac[i]) + beta * BCG_P(ac[i]);
                    BCG_P1_AUX(ac[i]) = BCG_R_AUX(ac[i]) + beta * BCG_P_AUX(ac[i]);
                }
                BCG_P(ac[i]) = BCG_P1(ac[i]);
                BCG_P_AUX(ac[i]) = BCG_P1_AUX(ac[i]);
            }

            rTz = r1Tz1;
            rTr = r1Tr1;

        }

    }//end of biconjugate gradient iterations.

}//end biconjugateGradient() function.
