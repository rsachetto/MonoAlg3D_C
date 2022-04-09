//
// Created by sachetto on 21/08/2020.
//

struct gpu_persistent_data {
    int *d_col, *d_row;
    float *d_val, *d_x;

    float *d_r, *d_p, *d_Ax, *d_rw, *d_v, *d_t;
    int N, nz;

    /* Get handle to the CUBLAS context */
    cublasHandle_t cublasHandle;
    cublasStatus_t cublasStatus;

    /* Get handle to the CUSPARSE context */
    cusparseHandle_t cusparseHandle;
    cusparseStatus_t cusparseStatus;

    cusparseMatDescr_t descr;

    int nzILU0;
    float *d_valsILU0;
    float *d_zm1, *d_zm2, *d_rm2;
    float *d_y;

    cusparseSolveAnalysisInfo_t infoA;
    cusparseSolveAnalysisInfo_t info_u;
    cusparseMatDescr_t descrL;
    cusparseMatDescr_t descrU;
};

 INIT_LINEAR_SYSTEM(init_gpu_conjugate_gradient) {

     if(the_grid->adaptive) {
         log_error_and_exit("The gpu conjugate gradient does not support mesh adaptivity. Aborting!\n");
     }

    struct gpu_persistent_data *persistent_data = CALLOC_ONE_TYPE(struct gpu_persistent_data);

    int_array I = NULL, J = NULL;
    f32_array val = NULL;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, tol, config, "tolerance");
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, max_its, config, "max_iterations");
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(use_preconditioner, config, "use_preconditioner");

    check_cuda_error((cudaError_t)cublasCreate(&(persistent_data->cublasHandle)));

    check_cuda_error((cudaError_t)cusparseCreate(&(persistent_data->cusparseHandle)));

    check_cuda_error((cudaError_t)cusparseCreateMatDescr(&(persistent_data->descr)));

    cusparseSetMatType(persistent_data->descr, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(persistent_data->descr, CUSPARSE_INDEX_BASE_ZERO);

    uint32_t num_active_cells;
    struct cell_node **active_cells = NULL;

    if (is_purkinje) {
        grid_to_csr(the_grid, &val, &I, &J, true);
        num_active_cells = the_grid->purkinje->num_active_purkinje_cells;
        active_cells = the_grid->purkinje->purkinje_cells;
    }
    else {
        grid_to_csr(the_grid, &val, &I, &J, false);
        num_active_cells = the_grid->num_active_cells;
        active_cells = the_grid->active_cells;
    }

    int N = 0, nz = 0;

    nz = arrlen(val);
    N = num_active_cells;

    persistent_data->N = N;
    persistent_data->nz = nz;

    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_col), nz * sizeof(int)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_row), (N + 1) * sizeof(int)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_val), nz * sizeof(float)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_x)  , N * sizeof(float)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_r) , N * sizeof(float)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_p) , N * sizeof(float)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_Ax), N * sizeof(float)));

    cudaMemcpy(persistent_data->d_col, J, nz * sizeof(int), cudaMemcpyHostToDevice);      // JA
    cudaMemcpy(persistent_data->d_row, I, (N + 1) * sizeof(int), cudaMemcpyHostToDevice); // IA
    cudaMemcpy(persistent_data->d_val, val, nz * sizeof(float), cudaMemcpyHostToDevice);  // A
    float *rhs = (float *)malloc(sizeof(float) * num_active_cells);

    OMP(parallel for)
    for(uint32_t i = 0; i < num_active_cells; i++) {
        rhs[i] = active_cells[i]->b;
    }

    check_cuda_error(cudaMemcpy(persistent_data->d_x, rhs, N * sizeof(float), cudaMemcpyHostToDevice)); // Result
    if(use_preconditioner) {
        persistent_data->nzILU0 = 2 * N - 1;
        check_cuda_error(cudaMalloc((void **)&(persistent_data->d_valsILU0), nz * sizeof(float)));
        check_cuda_error(cudaMalloc((void **)&(persistent_data->d_zm1), (N) * sizeof(float)));
        check_cuda_error(cudaMalloc((void **)&(persistent_data->d_zm2), (N) * sizeof(float)));
        check_cuda_error(cudaMalloc((void **)&(persistent_data->d_rm2), (N) * sizeof(float)));
        check_cuda_error(cudaMalloc((void **)&(persistent_data->d_y), N * sizeof(float)));

        persistent_data->cusparseStatus = cusparseCreateSolveAnalysisInfo(&(persistent_data->infoA));
        check_cuda_error((cudaError_t)(persistent_data->cusparseStatus));

        /* Perform the analysis for the Non-Transpose case */
        persistent_data->cusparseStatus = cusparseScsrsv_analysis(persistent_data->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, nz, persistent_data->descr, persistent_data->d_val, persistent_data->d_row, persistent_data->d_col, persistent_data->infoA);

        check_cuda_error((cudaError_t)(persistent_data->cusparseStatus));

        /* Copy A data to ILU0 vals as input*/
        cudaMemcpy(persistent_data->d_valsILU0, persistent_data->d_val, persistent_data->nz * sizeof(float), cudaMemcpyDeviceToDevice);

        /* generate the Incomplete LU factor H for the matrix A using cudsparseScsrilu0 */
        persistent_data->cusparseStatus = cusparseScsrilu0(persistent_data->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, persistent_data->N, persistent_data->descr, persistent_data->d_valsILU0, persistent_data->d_row, persistent_data->d_col, persistent_data->infoA);

        check_cuda_error((cudaError_t)(persistent_data->cusparseStatus));

        cusparseCreateSolveAnalysisInfo(&(persistent_data->info_u));

        persistent_data->cusparseStatus = cusparseCreateMatDescr(&(persistent_data->descrL));
        cusparseSetMatType(persistent_data->descrL, CUSPARSE_MATRIX_TYPE_GENERAL);
        cusparseSetMatIndexBase(persistent_data->descrL, CUSPARSE_INDEX_BASE_ZERO);
        cusparseSetMatFillMode(persistent_data->descrL, CUSPARSE_FILL_MODE_LOWER);
        cusparseSetMatDiagType(persistent_data->descrL, CUSPARSE_DIAG_TYPE_UNIT);

        persistent_data->cusparseStatus = cusparseCreateMatDescr(&(persistent_data->descrU));
        cusparseSetMatType(persistent_data->descrU, CUSPARSE_MATRIX_TYPE_GENERAL);
        cusparseSetMatIndexBase(persistent_data->descrU, CUSPARSE_INDEX_BASE_ZERO);
        cusparseSetMatFillMode(persistent_data->descrU, CUSPARSE_FILL_MODE_UPPER);
        cusparseSetMatDiagType(persistent_data->descrU, CUSPARSE_DIAG_TYPE_NON_UNIT);
        persistent_data->cusparseStatus = cusparseScsrsv_analysis(persistent_data->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, nz,persistent_data->descrU, persistent_data->d_val, persistent_data->d_row, persistent_data->d_col, persistent_data->info_u);
    }

    config->persistent_data = persistent_data;

    free(rhs);
    arrfree(I);
    arrfree(J);
    arrfree(val);
}

END_LINEAR_SYSTEM(end_gpu_conjugate_gradient) {

    struct gpu_persistent_data *persistent_data = (struct gpu_persistent_data *)config->persistent_data;

    if(!persistent_data) return;

    check_cuda_error((cudaError_t)cusparseDestroy(persistent_data->cusparseHandle));
    check_cuda_error((cudaError_t)cublasDestroy(persistent_data->cublasHandle));
    check_cuda_error((cudaError_t)cusparseDestroyMatDescr(persistent_data->descr));

    check_cuda_error(cudaFree(persistent_data->d_col));
    check_cuda_error(cudaFree(persistent_data->d_row));
    check_cuda_error(cudaFree(persistent_data->d_val));
    check_cuda_error(cudaFree(persistent_data->d_x));
    check_cuda_error(cudaFree(persistent_data->d_r));
    check_cuda_error(cudaFree(persistent_data->d_p));
    check_cuda_error(cudaFree(persistent_data->d_Ax));
    check_cuda_error(cudaFree(persistent_data->d_y));

    if(use_preconditioner) {
        /* Destroy parameters */
        cusparseDestroySolveAnalysisInfo(persistent_data->infoA);
        cusparseDestroySolveAnalysisInfo(persistent_data->info_u);
        check_cuda_error((cudaError_t)cusparseDestroyMatDescr(persistent_data->descrL));
        check_cuda_error((cudaError_t)cusparseDestroyMatDescr(persistent_data->descrU));
        check_cuda_error(cudaFree(persistent_data->d_valsILU0));
        check_cuda_error(cudaFree(persistent_data->d_zm1));
        check_cuda_error(cudaFree(persistent_data->d_zm2));
        check_cuda_error(cudaFree(persistent_data->d_rm2));
    }

    free(persistent_data);
}

SOLVE_LINEAR_SYSTEM(gpu_conjugate_gradient) {

    /* Conjugate gradient without preconditioning.
       ------------------------------------------
       Follows the description by Golub & Van Loan, "Matrix Computations 3rd ed.", Section 10.2.6
    */

    struct gpu_persistent_data *persistent_data = (struct gpu_persistent_data*)config->persistent_data;

    if(!persistent_data) {
        log_error_and_exit("[ERROR] The gpu_conjugate_gradient solver needs to be initialized before being called. Add a init_function in the [linear_system_solver] section of the .ini file\n");
    }

    float dot;

    float a, b, na, r0, r1;

    int k;
    float alpha, alpham1, beta;

    int N = persistent_data->N;
    int nz = persistent_data->nz;

    float *rhs; // Vector B
    rhs = (float *)malloc(sizeof(float) * num_active_cells);

    OMP(parallel for)
    for(uint32_t i = 0; i < num_active_cells; i++) {
        rhs[i] = active_cells[i]->b;
    }

    cudaMemcpy(persistent_data->d_r, rhs, N * sizeof(float), cudaMemcpyHostToDevice); // B

    alpha = 1.0;
    alpham1 = -1.0;
    beta = 0.0;
    r0 = 0.;

    float numerator, denominator;


    cusparseScsrmv(persistent_data->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &alpha, persistent_data->descr, persistent_data->d_val, persistent_data->d_row, persistent_data->d_col, persistent_data->d_x, &beta, persistent_data->d_Ax);

    cublasSaxpy(persistent_data->cublasHandle, N, &alpham1, persistent_data->d_Ax, 1, persistent_data->d_r, 1);
    persistent_data->cublasStatus = cublasSdot(persistent_data->cublasHandle, N, persistent_data->d_r, 1, persistent_data->d_r, 1, &r1);

    k = 1;

    while(r1 >= tol && k <= max_its) {

        if(use_preconditioner) {
            // Forward Solve, we can re-use infoA since the sparsity pattern of A matches that of L
            persistent_data->cusparseStatus =
                cusparseScsrsv_solve(persistent_data->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N,
                                     &alpha, persistent_data->descrL, persistent_data->d_valsILU0, persistent_data->d_row, persistent_data->d_col,
                                     persistent_data->infoA, persistent_data->d_r, persistent_data->d_y);
            check_cuda_error((cudaError_t)persistent_data->cusparseStatus);

            // Back Substitution
            persistent_data->cusparseStatus =
                cusparseScsrsv_solve(persistent_data->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, &alpha, persistent_data->descrU,
                                     persistent_data->d_valsILU0, persistent_data->d_row, persistent_data->d_col, persistent_data->info_u, persistent_data->d_y,
                                     persistent_data->d_zm1);

            check_cuda_error((cudaError_t)persistent_data->cusparseStatus);
        }

        if(k > 1) {
            if(use_preconditioner) {
                cublasSdot(persistent_data->cublasHandle, N, persistent_data->d_r, 1,   persistent_data->d_zm1, 1, &numerator);
                cublasSdot(persistent_data->cublasHandle, N, persistent_data->d_rm2, 1, persistent_data->d_zm2, 1, &denominator);
                b = numerator / denominator;
                cublasSscal(persistent_data->cublasHandle, N, &b, persistent_data->d_p, 1);
                cublasSaxpy(persistent_data->cublasHandle, N, &alpha, persistent_data->d_zm1, 1, persistent_data->d_p, 1);

            } else {
                b = r1 / r0;
                persistent_data->cublasStatus = cublasSscal(persistent_data->cublasHandle, N, &b, persistent_data->d_p, 1);
                persistent_data->cublasStatus = cublasSaxpy(persistent_data->cublasHandle, N, &alpha, persistent_data->d_r, 1, persistent_data->d_p, 1);
            }
        } else {
            if(use_preconditioner) {
                cublasScopy(persistent_data->cublasHandle, N, persistent_data->d_zm1, 1, persistent_data->d_p, 1);
            } else {
                persistent_data->cublasStatus = cublasScopy(persistent_data->cublasHandle, N, persistent_data->d_r, 1, persistent_data->d_p, 1);
            }
        }

        if(use_preconditioner) {
            cusparseScsrmv(persistent_data->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, persistent_data->nzILU0, &alpha, persistent_data->descrU, persistent_data->d_val, persistent_data->d_row, persistent_data->d_col, persistent_data->d_p, &beta, persistent_data->d_Ax);
            cublasSdot(persistent_data->cublasHandle, N, persistent_data->d_r, 1, persistent_data->d_zm1, 1, &numerator);
            cublasSdot(persistent_data->cublasHandle, N, persistent_data->d_p, 1, persistent_data->d_Ax, 1, &denominator);

            a = numerator / denominator;

            cublasSaxpy(persistent_data->cublasHandle, N, &a, persistent_data->d_p, 1, persistent_data->d_x, 1);
            cublasScopy(persistent_data->cublasHandle, N, persistent_data->d_r, 1, persistent_data->d_rm2, 1);
            cublasScopy(persistent_data->cublasHandle, N,persistent_data->d_zm1, 1, persistent_data->d_zm2, 1);
            na = -a;
            cublasSaxpy(persistent_data->cublasHandle, N, &na, persistent_data->d_Ax, 1, persistent_data->d_r, 1);

            r0 = r1;
            cublasSdot(persistent_data->cublasHandle, N, persistent_data->d_r, 1, persistent_data->d_r, 1, &r1);

        } else {
            cusparseScsrmv(persistent_data->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &alpha, persistent_data->descr, persistent_data->d_val,
                           persistent_data->d_row, persistent_data->d_col, persistent_data->d_p, &beta, persistent_data->d_Ax);

            persistent_data->cublasStatus = cublasSdot(persistent_data->cublasHandle, N, persistent_data->d_p, 1, persistent_data->d_Ax, 1, &dot);
            a = r1 / dot;

            persistent_data->cublasStatus = cublasSaxpy(persistent_data->cublasHandle, N, &a, persistent_data->d_p, 1, persistent_data->d_x, 1);
            na = -a;
            persistent_data->cublasStatus = cublasSaxpy(persistent_data->cublasHandle, N, &na, persistent_data->d_Ax, 1, persistent_data->d_r, 1);

            r0 = r1;
            persistent_data->cublasStatus = cublasSdot(persistent_data->cublasHandle, N, persistent_data->d_r, 1, persistent_data->d_r, 1, &r1);
        }

        cudaDeviceSynchronize();
        k++;
    }

    cudaMemcpy(rhs, persistent_data->d_x, N * sizeof(float), cudaMemcpyDeviceToHost);

    *number_of_iterations = k - 1;
    *error = r1;

    OMP(parallel for)
    for(uint32_t i = 0; i < num_active_cells; i++) {
        active_cells[i]->v = rhs[i];
    }

    free(rhs);
}

INIT_LINEAR_SYSTEM(init_gpu_biconjugate_gradient) {

    int_array I = NULL, J = NULL;
    f32_array val = NULL;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, tol, config, "tolerance");
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, max_its, config, "max_iterations");

    struct gpu_persistent_data *persistent_data = CALLOC_ONE_TYPE(struct gpu_persistent_data);

    check_cuda_error((cudaError_t)cublasCreate(&(persistent_data->cublasHandle)));

    check_cuda_error((cudaError_t)cusparseCreate(&(persistent_data->cusparseHandle)));

    check_cuda_error((cudaError_t)cusparseCreateMatDescr(&(persistent_data->descr)));

    cusparseSetMatType(persistent_data->descr, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(persistent_data->descr, CUSPARSE_INDEX_BASE_ZERO);

    grid_to_csr(the_grid, &val, &I, &J, false);

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **active_cells = the_grid->active_cells;

    int nz = arrlen(val);
    int N = num_active_cells;

    persistent_data->N = N;
    persistent_data->nz = nz;

    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_col), nz * sizeof(int)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_row), (N + 1) * sizeof(int)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_val), nz * sizeof(float)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_x), N * sizeof(float)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_r), N * sizeof(float)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_rw), N * sizeof(float)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_p), N * sizeof(float)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_v), N * sizeof(float)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_t), N * sizeof(float)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_Ax), N * sizeof(float)));

    cudaMemcpy(persistent_data->d_col, J, nz * sizeof(int), cudaMemcpyHostToDevice);      // JA
    cudaMemcpy(persistent_data->d_row, I, (N + 1) * sizeof(int), cudaMemcpyHostToDevice); // IA
    cudaMemcpy(persistent_data->d_val, val, nz * sizeof(float), cudaMemcpyHostToDevice);  // A

    check_cuda_error(cudaMemset(persistent_data->d_r, 0, N * sizeof(float)));
    check_cuda_error(cudaMemset(persistent_data->d_rw,0, N * sizeof(float)));
    check_cuda_error(cudaMemset(persistent_data->d_p, 0, N * sizeof(float)));
    check_cuda_error(cudaMemset(persistent_data->d_t, 0, N * sizeof(float)));
    check_cuda_error(cudaMemset(persistent_data->d_v, 0, N * sizeof(float)));

    float *rhs = (float *)malloc(sizeof(float) * num_active_cells);

    OMP(parallel for)
    for(uint32_t i = 0; i < num_active_cells; i++) {
        rhs[i] = active_cells[i]->b;
    }

    check_cuda_error(cudaMemcpy(persistent_data->d_x, rhs, N * sizeof(float), cudaMemcpyHostToDevice)); // Result

    free(rhs);
    arrfree(I);
    arrfree(J);

    config->persistent_data = persistent_data;
}

SOLVE_LINEAR_SYSTEM(gpu_biconjugate_gradient) {

    struct gpu_persistent_data *persistent_data = (struct gpu_persistent_data*) config->persistent_data;

    if(!persistent_data) {
        log_error_and_exit("[ERROR] The gpu_biconjugate_gradient solver needs to be initialized before being called. Add a init_function in the [linear_system_solver] section of the .ini file\n");
    }

    float rho, rhop, beta, alpha, negalpha, omega, negomega, temp, temp2;
    float nrmr, nrmr0;
    rho = 0.0f;
    float zero = 0.0f;
    float one = 1.0f;
    float mone = -1.0f;
    int i = 0;

    float *rhs; // Vector B
    rhs = (float *)malloc(sizeof(float) * num_active_cells);

    OMP(parallel for)
    for(uint32_t i = 0; i < num_active_cells; i++) {
        rhs[i] = active_cells[i]->b;
    }

    int N = persistent_data->N;
    int nz = persistent_data->nz;

    cudaMemcpy(persistent_data->d_r, rhs, N * sizeof(float), cudaMemcpyHostToDevice); // B

    // compute initial residual r0=b-Ax0 (using initial guess in x)
    check_cuda_error(
        (cudaError_t)cusparseScsrmv(persistent_data->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &one, persistent_data->descr, persistent_data->d_val, persistent_data->d_row, persistent_data->d_col, persistent_data->d_x, &zero, persistent_data->d_Ax));
    check_cuda_error((cudaError_t)cublasSaxpy(persistent_data->cublasHandle, N, &mone, persistent_data->d_Ax, 1, persistent_data->d_r, 1));

    // copy residual r into r^{\hat} and p
    check_cuda_error((cudaError_t)cublasScopy(persistent_data->cublasHandle, N, persistent_data->d_r, 1, persistent_data->d_rw, 1));
    check_cuda_error((cudaError_t)cublasScopy(persistent_data->cublasHandle, N, persistent_data->d_r, 1, persistent_data->d_p, 1));
    check_cuda_error((cudaError_t)cublasSnrm2(persistent_data->cublasHandle, N, persistent_data->d_r, 1, &nrmr0));


    for(i = 0; i < max_its;) {
        rhop = rho;
        check_cuda_error((cudaError_t)cublasSdot(persistent_data->cublasHandle, N, persistent_data->d_rw, 1, persistent_data->d_r, 1, &rho));

        if(i > 0) {
            beta = (rho / rhop) * (alpha / omega);
            negomega = -omega;
            check_cuda_error((cudaError_t)cublasSaxpy(persistent_data->cublasHandle, N, &negomega, persistent_data->d_v, 1, persistent_data->d_p, 1));
            check_cuda_error((cudaError_t)cublasSscal(persistent_data->cublasHandle, N, &beta, persistent_data->d_p, 1));
            check_cuda_error((cudaError_t)cublasSaxpy(persistent_data->cublasHandle, N, &one, persistent_data->d_r, 1, persistent_data->d_p, 1));
        }

        // matrix-vector multiplication
        check_cuda_error(
            (cudaError_t)cusparseScsrmv(persistent_data->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &one, persistent_data->descr, persistent_data->d_val, persistent_data->d_row, persistent_data->d_col, persistent_data->d_p, &zero, persistent_data->d_v));
        check_cuda_error((cudaError_t)cublasSdot(persistent_data->cublasHandle, N, persistent_data->d_rw, 1, persistent_data->d_v, 1, &temp));
        alpha = rho / temp;
        negalpha = -(alpha);
        check_cuda_error((cudaError_t)cublasSaxpy(persistent_data->cublasHandle, N, &negalpha, persistent_data->d_v, 1, persistent_data->d_r, 1));
        check_cuda_error((cudaError_t)cublasSaxpy(persistent_data->cublasHandle, N, &alpha, persistent_data->d_p, 1, persistent_data->d_x, 1));
        check_cuda_error((cudaError_t)cublasSnrm2(persistent_data->cublasHandle, N, persistent_data->d_r, 1, &nrmr));

        if(nrmr < tol * nrmr0) {
            break;
        }

        // matrix-vector multiplication
        check_cuda_error((cudaError_t)cusparseScsrmv(persistent_data->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &one, persistent_data->descr, persistent_data->d_val, persistent_data->d_row, persistent_data->d_col, persistent_data->d_r, &zero, persistent_data->d_t));

        check_cuda_error((cudaError_t)cublasSdot(persistent_data->cublasHandle, N, persistent_data->d_t, 1, persistent_data->d_r, 1, &temp));
        check_cuda_error((cudaError_t)cublasSdot(persistent_data->cublasHandle, N, persistent_data->d_t, 1, persistent_data->d_t, 1, &temp2));
        omega = temp / temp2;
        negomega = -(omega);
        check_cuda_error((cudaError_t)cublasSaxpy(persistent_data->cublasHandle, N, &omega, persistent_data->d_r, 1, persistent_data->d_x, 1));
        check_cuda_error((cudaError_t)cublasSaxpy(persistent_data->cublasHandle, N, &negomega,persistent_data-> d_t, 1, persistent_data->d_r, 1));
        check_cuda_error((cudaError_t)cublasSnrm2(persistent_data->cublasHandle, N, persistent_data->d_r, 1, &nrmr));

        if(nrmr < tol * nrmr0) {
            i++;
            break;
        }
        i++;
    }

    cudaMemcpy(rhs, persistent_data->d_x, N * sizeof(float), cudaMemcpyDeviceToHost);

    *number_of_iterations = i - 1;
    *error = nrmr;

    OMP(parallel for)
    for(uint32_t i = 0; i < num_active_cells; i++) {
        active_cells[i]->v = rhs[i];
    }

    free(rhs);
}
