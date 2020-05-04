//
// Created by sachetto on 04/10/17.
//
#include "../config/linear_system_solver_config.h"
#include "../config_helpers/config_helpers.h"
#include "../libraries_common/common_data_structures.h"
#include "../utils/file_utils.h"

#include "../3dparty/stb_ds.h"

bool jacobi_initialized = false;
bool bcg_initialized = false;
static bool use_preconditioner = false;
static int max_its = 50;
static real_cpu tol = 1e-16;

#ifdef COMPILE_CUDA
#include "../models_library/model_gpu_utils.h"
#include <cusparse_v2.h>
#include <cublas_v2.h>

static int *d_col, *d_row;
static float *d_val, *d_x;

static float *d_r, *d_p, *d_Ax;

static int N = 0, nz = 0;

/* Get handle to the CUBLAS context */
static cublasHandle_t cublasHandle = 0;
static cublasStatus_t cublasStatus;

/* Get handle to the CUSPARSE context */
static cusparseHandle_t cusparseHandle = 0;
static cusparseStatus_t cusparseStatus;

cusparseMatDescr_t descr = 0;

static int nzILU0;
static float *d_valsILU0;
static float *d_zm1, *d_zm2, *d_rm2;
static float *d_y;

static cusparseSolveAnalysisInfo_t infoA = 0;
static cusparseSolveAnalysisInfo_t info_u;
static cusparseMatDescr_t descrL = 0;
static cusparseMatDescr_t descrU = 0;

INIT_LINEAR_SYSTEM(init_gpu_conjugate_gradient) {

    int_array I = NULL, J = NULL;
    f32_array val = NULL;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, tol, config->config_data, "tolerance");
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, max_its, config->config_data, "max_iterations");
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(use_preconditioner, config->config_data, "use_preconditioner");

    check_cuda_error((cudaError_t)cublasCreate(&cublasHandle));

    check_cuda_error((cudaError_t)cusparseCreate(&cusparseHandle));

    check_cuda_error((cudaError_t)cusparseCreateMatDescr(&descr));

    cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descr,CUSPARSE_INDEX_BASE_ZERO);

    grid_to_csr(the_grid, &val, &I, &J);

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node** active_cells = the_grid->active_cells;

    nz = arrlen(val);
    N = num_active_cells;

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
    float *rhs = (float*) malloc(sizeof(float)*num_active_cells);

    OMP(parallel for)
    for (uint32_t i = 0; i < num_active_cells; i++) {
        rhs[i] = active_cells[i]->b;
    }

    check_cuda_error(cudaMemcpy(d_x, rhs, N * sizeof(float), cudaMemcpyHostToDevice)); //Result

    if(use_preconditioner) {
        nzILU0 = 2*N-1;
        check_cuda_error(cudaMalloc((void **)&d_valsILU0, nz*sizeof(float)));
        check_cuda_error(cudaMalloc((void **)&d_zm1, (N)*sizeof(float)));
        check_cuda_error(cudaMalloc((void **)&d_zm2, (N)*sizeof(float)));
        check_cuda_error(cudaMalloc((void **)&d_rm2, (N)*sizeof(float)));
        check_cuda_error(cudaMalloc((void **)&d_y, N*sizeof(float)));


        cusparseStatus = cusparseCreateSolveAnalysisInfo(&infoA);
        check_cuda_error((cudaError_t)cusparseStatus);

        /* Perform the analysis for the Non-Transpose case */
        cusparseStatus = cusparseScsrsv_analysis(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                                 N, nz, descr, d_val, d_row, d_col, infoA);

        check_cuda_error((cudaError_t)cusparseStatus);

        /* Copy A data to ILU0 vals as input*/
        cudaMemcpy(d_valsILU0, d_val, nz*sizeof(float), cudaMemcpyDeviceToDevice);

        /* generate the Incomplete LU factor H for the matrix A using cudsparseScsrilu0 */
        cusparseStatus = cusparseScsrilu0(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, descr, d_valsILU0, d_row, d_col, infoA);

        check_cuda_error((cudaError_t)cusparseStatus);

        cusparseCreateSolveAnalysisInfo(&info_u);

        cusparseStatus = cusparseCreateMatDescr(&descrL);
        cusparseSetMatType(descrL,CUSPARSE_MATRIX_TYPE_GENERAL);
        cusparseSetMatIndexBase(descrL,CUSPARSE_INDEX_BASE_ZERO);
        cusparseSetMatFillMode(descrL, CUSPARSE_FILL_MODE_LOWER);
        cusparseSetMatDiagType(descrL, CUSPARSE_DIAG_TYPE_UNIT);

        cusparseStatus = cusparseCreateMatDescr(&descrU);
        cusparseSetMatType(descrU,CUSPARSE_MATRIX_TYPE_GENERAL);
        cusparseSetMatIndexBase(descrU,CUSPARSE_INDEX_BASE_ZERO);
        cusparseSetMatFillMode(descrU, CUSPARSE_FILL_MODE_UPPER);
        cusparseSetMatDiagType(descrU, CUSPARSE_DIAG_TYPE_NON_UNIT);
        cusparseStatus = cusparseScsrsv_analysis(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, nz, descrU, d_val, d_row, d_col, info_u);

    }

    free(rhs);
    arrfree(I);
    arrfree(J);
    arrfree(val);
}

END_LINEAR_SYSTEM(end_gpu_conjugate_gradient) {

    check_cuda_error( (cudaError_t)cusparseDestroy(cusparseHandle) );
    check_cuda_error( (cudaError_t)cublasDestroy(cublasHandle) );
    check_cuda_error( (cudaError_t)cusparseDestroyMatDescr(descr));


    check_cuda_error(cudaFree(d_col));
    check_cuda_error(cudaFree(d_row));
    check_cuda_error(cudaFree(d_val));
    check_cuda_error(cudaFree(d_x));
    check_cuda_error(cudaFree(d_r));
    check_cuda_error(cudaFree(d_p));
    check_cuda_error(cudaFree(d_Ax));
    check_cuda_error(cudaFree(d_y));

	if(use_preconditioner) {
	    /* Destroy parameters */
	    cusparseDestroySolveAnalysisInfo(infoA);
	    cusparseDestroySolveAnalysisInfo(info_u);
		check_cuda_error( (cudaError_t)cusparseDestroyMatDescr(descrL));
	    check_cuda_error( (cudaError_t)cusparseDestroyMatDescr(descrU));
		check_cuda_error(cudaFree(d_valsILU0));
		check_cuda_error(cudaFree(d_zm1));
		check_cuda_error(cudaFree(d_zm2));
		check_cuda_error(cudaFree(d_rm2));
	}
}

SOLVE_LINEAR_SYSTEM(gpu_conjugate_gradient) {

    /* Conjugate gradient without preconditioning.
       ------------------------------------------
       Follows the description by Golub & Van Loan, "Matrix Computations 3rd ed.", Section 10.2.6
    */
    float dot;

    float a, b, na, r0, r1;

    int k;
    float alpha, alpham1, beta;

    float *rhs; //Vector B
    rhs = (float*) malloc(sizeof(float)*num_active_cells);

    OMP(parallel for)
    for (uint32_t i = 0; i < num_active_cells; i++) {
        rhs[i] = active_cells[i]->b;
    }

    cudaMemcpy(d_r, rhs, N * sizeof(float), cudaMemcpyHostToDevice); //B

    alpha = 1.0;
    alpham1 = -1.0;
    beta = 0.0;
    r0 = 0.;

    float numerator, denominator;

    cusparseScsrmv(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &alpha, descr, d_val, d_row, d_col, d_x, &beta, d_Ax);

    cublasSaxpy(cublasHandle, N, &alpham1, d_Ax, 1, d_r, 1);
    cublasStatus = cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);

    k = 1;

    while (r1 >= tol && k <= max_its)
    {

        if(use_preconditioner) {
            // Forward Solve, we can re-use infoA since the sparsity pattern of A matches that of L
            cusparseStatus = cusparseScsrsv_solve(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, &alpha, descrL,
                                                  d_valsILU0, d_row, d_col, infoA, d_r, d_y);
            check_cuda_error((cudaError_t)cusparseStatus);

            // Back Substitution
            cusparseStatus = cusparseScsrsv_solve(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, &alpha, descrU,
                                                  d_valsILU0, d_row, d_col, info_u, d_y, d_zm1);
            check_cuda_error((cudaError_t)cusparseStatus);

        }

        if (k > 1)
        {
            if(use_preconditioner) {
                cublasSdot(cublasHandle, N, d_r, 1, d_zm1, 1, &numerator);
                cublasSdot(cublasHandle, N, d_rm2, 1, d_zm2, 1, &denominator);
                b = numerator/denominator;
                cublasSscal(cublasHandle, N, &b, d_p, 1);
                cublasSaxpy(cublasHandle, N, &alpha, d_zm1, 1, d_p, 1) ;

            }
            else {
                b = r1 / r0;
                cublasStatus = cublasSscal(cublasHandle, N, &b, d_p, 1);
                cublasStatus = cublasSaxpy(cublasHandle, N, &alpha, d_r, 1, d_p, 1);
            }
        }
        else
        {
            if(use_preconditioner) {
                cublasScopy(cublasHandle, N, d_zm1, 1, d_p, 1);
            }
            else {
                cublasStatus = cublasScopy(cublasHandle, N, d_r, 1, d_p, 1);
            }
        }

        if(use_preconditioner) {
            cusparseScsrmv(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nzILU0, &alpha, descrU, d_val, d_row, d_col, d_p, &beta, d_Ax);
            cublasSdot(cublasHandle, N, d_r, 1, d_zm1, 1, &numerator);
            cublasSdot(cublasHandle, N, d_p, 1, d_Ax, 1, &denominator);

            a = numerator / denominator;

            cublasSaxpy(cublasHandle, N, &a, d_p, 1, d_x, 1);
            cublasScopy(cublasHandle, N, d_r, 1, d_rm2, 1);
            cublasScopy(cublasHandle, N, d_zm1, 1, d_zm2, 1);
            na = -a;
            cublasSaxpy(cublasHandle, N, &na, d_Ax, 1, d_r, 1);

            r0 = r1;
            cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);


        }
        else {
            cusparseScsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &alpha, descr, d_val, d_row,
                           d_col, d_p, &beta, d_Ax);


            cublasStatus = cublasSdot(cublasHandle, N, d_p, 1, d_Ax, 1, &dot);
            a = r1 / dot;

            cublasStatus = cublasSaxpy(cublasHandle, N, &a, d_p, 1, d_x, 1);
            na = -a;
            cublasStatus = cublasSaxpy(cublasHandle, N, &na, d_Ax, 1, d_r, 1);

            r0 = r1;
            cublasStatus = cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
        }

        cudaDeviceSynchronize();
        k++;
    }

    cudaMemcpy(rhs, d_x, N*sizeof(float), cudaMemcpyDeviceToHost);

    *number_of_iterations = k-1;
    *error = r1;

    OMP(parallel for)
    for (uint32_t i = 0; i < num_active_cells; i++) {
        active_cells[i]->v = rhs[i];
    }

    free(rhs);
}

#endif

INIT_LINEAR_SYSTEM(init_cpu_conjugate_gradient) {
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, tol, config->config_data, "tolerance");
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(use_preconditioner, config->config_data, "use_preconditioner");
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, max_its, config->config_data, "max_iterations");
}

END_LINEAR_SYSTEM(end_cpu_conjugate_gradient) {
}

SOLVE_LINEAR_SYSTEM(cpu_conjugate_gradient) {

    real_cpu  rTr,
             r1Tr1,
             pTAp,
             alpha,
             beta,
             precision = tol,
             rTz,
             r1Tz1;


    *error = 1.0;
    *number_of_iterations = 1;

    //__________________________________________________________________________
    //Computes int_vector A*x, residue r = b - Ax, scalar rTr = r^T * r and
    //sets initial search direction p.

    rTr = 0.0;
    rTz = 0.0;

    uint32_t i;

    OMP(parallel for reduction(+:rTr,rTz))
    for (i = 0; i < num_active_cells; i++) {

        if(CG_INFO(active_cells[i]) == NULL) {
            INITIALIZE_CONJUGATE_GRADIENT_INFO(active_cells[i]);
        }

        struct element *cell_elements = active_cells[i]->elements;
        active_cells[i]->Ax = 0.0;
        size_t max_el = arrlen(cell_elements);

        for(size_t el = 0; el < max_el; el++) {
            active_cells[i]->Ax +=  cell_elements[el].value *  cell_elements[el].cell->v;
        }

        CG_R(active_cells[i]) = active_cells[i]->b - active_cells[i]->Ax;
        if(use_preconditioner) {
            real_cpu value = cell_elements[0].value;
            if(value == 0.0) value = 1.0;
            CG_Z(active_cells[i]) = (1.0/value) * CG_R(active_cells[i]); // preconditioner
            rTz += CG_R(active_cells[i]) * CG_Z(active_cells[i]);
            CG_P(active_cells[i]) = CG_Z(active_cells[i]);
        }
        else {
            CG_P(active_cells[i]) = CG_R(active_cells[i]);
        }

        real_cpu r = CG_R(active_cells[i]);

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

            OMP(parallel for reduction(+ : pTAp))
            for (i = 0; i < num_active_cells; i++) {

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
                alpha = rTz/pTAp;
            }
            else {
                alpha = rTr/pTAp;
            }
            //__________________________________________________________________

            r1Tr1 = 0.0;
            r1Tz1 = 0.0;

            // Computes new value of solution: u = u + alpha*p.
            OMP(parallel for reduction (+:r1Tr1,r1Tz1))
            for (i = 0; i < num_active_cells; i++) {

                active_cells[i]->v += alpha * CG_P(active_cells[i]);

                CG_R(active_cells[i]) -= alpha * active_cells[i]->Ax;

                real_cpu r = CG_R(active_cells[i]);

                if(use_preconditioner) {
                    real_cpu value = active_cells[i]->elements[0].value;
                    if(value == 0.0) value = 1.0;
                    CG_Z(active_cells[i]) = (1.0/value) * r;
                    r1Tz1 += CG_Z(active_cells[i]) * r;
                }
                r1Tr1 += r * r;
            }
            //__________________________________________________________________
            //Computes beta.
            if(use_preconditioner) {
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
            OMP(parallel for)
            for (i = 0; i < num_active_cells; i++) {
                if(use_preconditioner) {
                    CG_P1(active_cells[i]) = CG_Z(active_cells[i]) + beta * CG_P(active_cells[i]);
                }
                else {
                    CG_P1(active_cells[i]) = CG_R(active_cells[i]) + beta * CG_P(active_cells[i]);
                }
                CG_P(active_cells[i]) = CG_P1(active_cells[i]);
            }

            rTz = r1Tz1;
            rTr = r1Tr1;

        }

    }//end of conjugate gradient iterations.

}//end conjugateGradient() function.


SOLVE_LINEAR_SYSTEM(conjugate_gradient) {

    bool gpu = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(gpu, config->config_data, "use_gpu");

    if(gpu) {
    #ifdef COMPILE_CUDA
        gpu_conjugate_gradient(time_info, config, the_grid, num_active_cells, active_cells, number_of_iterations, error);
    #else
        log_to_stdout_and_file("Cuda runtime not found in this system. Fallbacking to CPU solver!!\n");
        cpu_conjugate_gradient(time_info, config, the_grid, num_active_cells, active_cells, number_of_iterations, error);
    #endif
    }
    else {
        cpu_conjugate_gradient(time_info, config, the_grid, num_active_cells, active_cells, number_of_iterations, error);
    }
}

INIT_LINEAR_SYSTEM(init_conjugate_gradient) {
    bool gpu = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(gpu, config->config_data, "use_gpu");

    if(gpu) {
#ifdef COMPILE_CUDA
        init_gpu_conjugate_gradient(config, the_grid);
#else
        init_cpu_conjugate_gradient(config, the_grid);
#endif
    }
    else {
        init_cpu_conjugate_gradient(config, the_grid);
    }

}

END_LINEAR_SYSTEM(end_conjugate_gradient) {

    bool gpu = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(gpu, config->config_data, "use_gpu");

    if(gpu) {
#ifdef COMPILE_CUDA
        end_gpu_conjugate_gradient(config);
#else
        end_cpu_conjugate_gradient(config);
#endif
    }
    else {
        end_cpu_conjugate_gradient(config);
    }

}

// Berg's code
SOLVE_LINEAR_SYSTEM(jacobi) {

    if(!jacobi_initialized) {
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, tol, config->config_data, "tolerance");
        max_its = 500;
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, max_its, config->config_data, "max_iterations");
        jacobi_initialized = true;
    }


    real_cpu  sigma,
            precision = tol;

    *error = 1.0;
    *number_of_iterations = 1;

    struct element element;

    if (*error >= precision)
    {
        //__________________________________________________________________________
        //Jacobi iterations.
        while (*number_of_iterations < max_its) {

            OMP(parallel for private (element,sigma))
            for (uint32_t i = 0; i < num_active_cells; i++) {
                if(JACOBI_INFO(active_cells[i]) == NULL) {
                    INITIALIZE_JACOBI_INFO(active_cells[i]);
                }

                struct element *cell_elements = active_cells[i]->elements;
                sigma = 0.0;

                size_t max_el = arrlen(cell_elements);

                // Do not take the diagonal element
                for(size_t el = 1; el < max_el; el++) {
                    element = cell_elements[el];
                    sigma += element.value * element.cell->v;
                }

                real_cpu value = cell_elements[0].value;
                JACOBI_X_AUX(active_cells[i]) = (1.0/value)*(active_cells[i]->b - sigma);
            }
            real_cpu residue = 0.0;
            real_cpu sum;

            OMP(parallel for private (element,sum) reduction (+:residue))
            for (uint32_t i = 0; i < num_active_cells; i++) {

                struct element *cell_elements = active_cells[i]->elements;

                size_t max_el = arrlen(cell_elements);

                // Do not take the diagonal element
                sum = 0.0;
                for(size_t el = 0; el < max_el; el++)
                {
                    element = cell_elements[el];
                    sum += element.value * JACOBI_X_AUX(element.cell);
                }

                active_cells[i]->v = JACOBI_X_AUX(active_cells[i]);
                residue += pow(active_cells[i]->b - sum,2);
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
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, tol, config->config_data, "tolerance");

        char *preconditioner_char = NULL;
        GET_PARAMETER_VALUE_CHAR_OR_USE_DEFAULT(preconditioner_char, config->config_data, "use_preconditioner");
        if (preconditioner_char != NULL)
        {
            use_preconditioner = ((strcmp (preconditioner_char, "yes") == 0) || (strcmp (preconditioner_char, "true") == 0));
        }

        max_its = 100;
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, max_its, config->config_data, "max_iterations");
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


    *error = 1.0;
    *number_of_iterations = 1;

    struct element element;
    //__________________________________________________________________________
    // Zero all entries on the int_vector x*A
    // And initialize the second guess vector x_aux
    OMP(parallel for)
    for (uint32_t i = 0; i < num_active_cells; i++)
    {

        if(BCG_INFO(active_cells[i]) == NULL) {
            INITIALIZE_BICONJUGATE_GRADIENT_INFO(active_cells[i]);
        }

        BCG_XA(active_cells[i]) = 0.0;
        BCG_X_AUX(active_cells[i]) = active_cells[i]->v;
    }


    //__________________________________________________________________________
    //Computes int_vector A*x, x*A
    //xA must be fully calculated to start doing anything over the r_aux vector
    OMP(parallel for private (element))
    for (uint32_t i = 0; i < num_active_cells; i++)
    {
        struct element *cell_elements = active_cells[i]->elements;
        active_cells[i]->Ax = 0.0;

        size_t max_el = arrlen(cell_elements);

        for(size_t el = 0; el < max_el; el++)
        {
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
    //Computes residues r, r_aux
    //scalar rTr = r^T * r_aux and
    //sets initial search directions p and p_aux.
    OMP(parallel for private (element) reduction(+:rTr,rTz))
    for (uint32_t i = 0; i < num_active_cells; i++)
    {
        struct element *cell_elements = active_cells[i]->elements;

        BCG_R(active_cells[i]) = active_cells[i]->b - active_cells[i]->Ax;
        BCG_R_AUX(active_cells[i]) = active_cells[i]->b - BCG_XA(active_cells[i]);

        if(use_preconditioner)
        {
            real_cpu value = cell_elements[0].value;
            if(value == 0.0) value = 1.0;
            BCG_Z(active_cells[i]) = (1.0/value) * BCG_R(active_cells[i]); // preconditioner
            BCG_Z_AUX(active_cells[i]) = (1.0/value) * BCG_R_AUX(active_cells[i]);
            rTz += BCG_R_AUX(active_cells[i]) * BCG_Z(active_cells[i]);
            BCG_P(active_cells[i]) = BCG_Z(active_cells[i]);
            BCG_P_AUX(active_cells[i]) = BCG_Z_AUX(active_cells[i]);
        }
        else
        {
            BCG_P(active_cells[i]) = BCG_R(active_cells[i]);
            BCG_P_AUX(active_cells[i])= BCG_R_AUX(active_cells[i]);
        }
        rTr += BCG_R_AUX(active_cells[i]) * BCG_R(active_cells[i]);
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

            OMP(parallel for)
            for (uint32_t i = 0; i < num_active_cells; i++)
                BCG_XA(active_cells[i]) = 0.0;

            OMP(parallel for private(element) reduction(+ : pTAp))
            for (uint32_t i = 0; i < num_active_cells; i++)
            {
                active_cells[i]->Ax = 0.0;
                struct element *cell_elements = active_cells[i]->elements;

                size_t max_el = arrlen(cell_elements);
                for(size_t el = 0; el < max_el; el++)
                {
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
            if(use_preconditioner)
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
            OMP(parallel for reduction (+:r1Tr1,r1Tz1))
            for (uint32_t i = 0; i < num_active_cells; i++) {

                active_cells[i]->v += alpha * BCG_P(active_cells[i]);
                BCG_X_AUX(active_cells[i]) += alpha * BCG_P_AUX(active_cells[i]);

                BCG_R(active_cells[i]) -= alpha * active_cells[i]->Ax;
                BCG_R_AUX(active_cells[i]) -= alpha * BCG_XA(active_cells[i]);

                if(use_preconditioner)
                {
                    real_cpu value = active_cells[i]->elements[0].value;
                    if(value == 0.0) value = 1.0;
                    BCG_Z(active_cells[i]) = (1.0/value) * BCG_R(active_cells[i]);
                    BCG_Z_AUX(active_cells[i]) = (1.0/value) * BCG_R_AUX(active_cells[i]);
                    r1Tz1 += BCG_Z(active_cells[i]) * BCG_R_AUX(active_cells[i]);
                }

                r1Tr1 += BCG_R(active_cells[i]) * BCG_R_AUX(active_cells[i]);
            }
            //__________________________________________________________________
            //Computes beta.
            if(use_preconditioner)
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
            OMP(parallel for)
            for (uint32_t i = 0; i < num_active_cells; i++)
            {
                if(use_preconditioner)
                {
                    BCG_P1(active_cells[i]) = BCG_Z(active_cells[i]) + beta * BCG_P(active_cells[i]);
                    BCG_P1_AUX(active_cells[i]) = BCG_Z_AUX(active_cells[i]) + beta * BCG_P_AUX(active_cells[i]);
                }
                else
                {
                    BCG_P1(active_cells[i]) = BCG_R(active_cells[i]) + beta * BCG_P(active_cells[i]);
                    BCG_P1_AUX(active_cells[i]) = BCG_R_AUX(active_cells[i]) + beta * BCG_P_AUX(active_cells[i]);
                }
                BCG_P(active_cells[i]) = BCG_P1(active_cells[i]);
                BCG_P_AUX(active_cells[i]) = BCG_P1_AUX(active_cells[i]);
            }

            rTz = r1Tz1;
            rTr = r1Tr1;

        }

    }//end of biconjugate gradient iterations.

}//end biconjugateGradient() function.
