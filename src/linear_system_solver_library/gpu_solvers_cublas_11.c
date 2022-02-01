//
// Created by sachetto on 21/08/2020.
//

struct gpu_persistent_data {
    int *d_col, *d_row;
    float *d_val, *d_x;

    float *d_r, *d_p, *d_Ax, *d_rw, *d_v, *d_t, *d_omega;

    int N, nz;

    /* Get handle to the CUBLAS context */
    cublasHandle_t cublasHandle;

    /* Get handle to the CUSPARSE context */
    cusparseHandle_t cusparseHandle;

    cusparseSpMatDescr_t matA;
    cusparseDnVecDescr_t vecx;
    cusparseDnVecDescr_t vecp;
    cusparseDnVecDescr_t vecv;
    cusparseDnVecDescr_t vecr;
    cusparseDnVecDescr_t vect;
    cusparseDnVecDescr_t vecAx;
    cusparseDnVecDescr_t vecomega;

    float *d_valsILU0;
    float *d_zm1, *d_zm2, *d_rm2;
    float *d_y;

    csrilu02Info_t infoILU;
    cusparseMatDescr_t descrL;
    csrsv2Info_t infoL;

    cusparseMatDescr_t descrU;
    csrsv2Info_t infoU;

    cusparseMatDescr_t descr;
    size_t bufferSize;
    void *buffer;
};

const float floatone = 1.0f;
const float floatzero = 0.0f;

INIT_LINEAR_SYSTEM(init_gpu_conjugate_gradient) {

    if(the_grid->adaptive) {
        log_error_and_exit("The gpu conjugate gradient does not support mesh adaptivity. Aborting!\n");
    }

    struct gpu_persistent_data *persistent_data = CALLOC_ONE_TYPE(struct gpu_persistent_data);

    size_t tmp;
    int stmp;

    int_array I = NULL, J = NULL;
    f32_array val = NULL;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, tol, config, "tolerance");
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, max_its, config, "max_iterations");
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(use_preconditioner, config, "use_preconditioner");

    check_cublas_error(cublasCreate(&(persistent_data->cublasHandle)));
    check_cublas_error(cusparseCreate(&(persistent_data->cusparseHandle)));

    uint32_t num_active_cells;
    struct cell_node **active_cells = NULL;

    if (is_purkinje) {
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

    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_col), nz * sizeof(int)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_row), (N + 1) * sizeof(int)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_val), nz * sizeof(float)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_x), N * sizeof(float)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_r), N * sizeof(float)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_p), N * sizeof(float)));
    check_cuda_error(cudaMalloc((void **)&(persistent_data->d_Ax), N * sizeof(float)));

    /* Wrap raw data into cuSPARSE generic API objects */
    check_cuda_error(cusparseCreateCsr(&(persistent_data->matA), N, N, nz, persistent_data->d_row, persistent_data->d_col, persistent_data->d_val, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F));
    check_cuda_error(cusparseCreateDnVec(&(persistent_data->vecx), N,  persistent_data->d_x, CUDA_R_32F));
    check_cuda_error(cusparseCreateDnVec(&(persistent_data->vecp), N,  persistent_data->d_p, CUDA_R_32F));
    check_cuda_error(cusparseCreateDnVec(&(persistent_data->vecAx), N, persistent_data->d_Ax, CUDA_R_32F));


    cudaMemcpy(persistent_data->d_col, J, nz * sizeof(int), cudaMemcpyHostToDevice);      // JA
    cudaMemcpy(persistent_data->d_row, I, (N + 1) * sizeof(int), cudaMemcpyHostToDevice); // IA
    cudaMemcpy(persistent_data->d_val, val, nz * sizeof(float), cudaMemcpyHostToDevice);  // A
    float *rhs = (float *)malloc(sizeof(float) * num_active_cells);

    OMP(parallel for)
    for(uint32_t i = 0; i < num_active_cells; i++) {
        rhs[i] = active_cells[i]->b;
    }

    check_cuda_error(cudaMemcpy(persistent_data->d_x, rhs, N * sizeof(float), cudaMemcpyHostToDevice)); // Result

    float alpha = 1.0f;
    float beta = 0.0f;

    check_cuda_error(cusparseSpMV_bufferSize(persistent_data->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, persistent_data->matA, persistent_data->vecx, &beta, persistent_data->vecAx, CUDA_R_32F, CUSPARSE_MV_ALG_DEFAULT, &(persistent_data->bufferSize)));
    check_cuda_error(cudaMalloc(&(persistent_data->buffer), persistent_data->bufferSize));

    if(use_preconditioner) {
        check_cuda_error(cudaMalloc((void **)&(persistent_data->d_valsILU0), nz * sizeof(float)));
        check_cuda_error(cudaMalloc((void **)&(persistent_data->d_zm1), (N) * sizeof(float)));
        check_cuda_error(cudaMalloc((void **)&(persistent_data->d_zm2), (N) * sizeof(float)));
        check_cuda_error(cudaMalloc((void **)&(persistent_data->d_rm2), (N) * sizeof(float)));
        check_cuda_error(cudaMalloc((void **)&(persistent_data->d_omega), N * sizeof(float)));
        check_cuda_error(cudaMalloc((void **)&(persistent_data->d_y), N * sizeof(float)));

        check_cuda_error(cusparseCreateDnVec(&(persistent_data->vecomega), N,persistent_data->d_omega, CUDA_R_32F));

        check_cuda_error(cusparseCreateCsrilu02Info(&(persistent_data->infoILU)));

        /* Create L factor descriptor and triangular solve info */
        check_cuda_error(cusparseCreateMatDescr(&(persistent_data->descrL)));
        check_cuda_error(cusparseSetMatType(persistent_data->descrL, CUSPARSE_MATRIX_TYPE_GENERAL));
        check_cuda_error(cusparseSetMatIndexBase(persistent_data->descrL, CUSPARSE_INDEX_BASE_ZERO));
        check_cuda_error(cusparseSetMatFillMode(persistent_data->descrL, CUSPARSE_FILL_MODE_LOWER));
        check_cuda_error(cusparseSetMatDiagType(persistent_data->descrL, CUSPARSE_DIAG_TYPE_UNIT));
        check_cuda_error(cusparseCreateCsrsv2Info(&(persistent_data->infoL)));

        /* Create U factor descriptor and triangular solve info */
        check_cuda_error(cusparseCreateMatDescr(&(persistent_data->descrU)));
        check_cuda_error(cusparseSetMatType(persistent_data->descrU, CUSPARSE_MATRIX_TYPE_GENERAL));
        check_cuda_error(cusparseSetMatIndexBase(persistent_data->descrU, CUSPARSE_INDEX_BASE_ZERO));
        check_cuda_error(cusparseSetMatFillMode(persistent_data->descrU, CUSPARSE_FILL_MODE_UPPER));
        check_cuda_error(cusparseSetMatDiagType(persistent_data->descrU, CUSPARSE_DIAG_TYPE_NON_UNIT));
        check_cuda_error(cusparseCreateCsrsv2Info(&(persistent_data->infoU)));

        check_cuda_error(cusparseCreateMatDescr(&(persistent_data->descr)));
        check_cuda_error(cusparseSetMatType(persistent_data->descr, CUSPARSE_MATRIX_TYPE_GENERAL));
        check_cuda_error(cusparseSetMatIndexBase(persistent_data->descr, CUSPARSE_INDEX_BASE_ZERO));

        check_cuda_error(cusparseSpMV_bufferSize(persistent_data->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &floatone, persistent_data->matA, persistent_data->vecp, &floatzero, persistent_data->vecomega, CUDA_R_32F, CUSPARSE_MV_ALG_DEFAULT, &tmp));
        if (tmp > persistent_data->bufferSize) {
            persistent_data->bufferSize = tmp;
        }
        check_cuda_error(cusparseScsrilu02_bufferSize(persistent_data->cusparseHandle, N, nz, persistent_data->descr, persistent_data->d_val, persistent_data->d_row, persistent_data->d_col, persistent_data->infoILU, &stmp));
        if (stmp > persistent_data->bufferSize) {
            persistent_data->bufferSize = stmp;
        }
        check_cuda_error(cusparseScsrsv2_bufferSize(persistent_data->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, nz, persistent_data->descrL, persistent_data->d_val, persistent_data->d_row, persistent_data->d_col, persistent_data->infoL, &stmp));
        if (stmp > persistent_data->bufferSize) {
            persistent_data->bufferSize = stmp;
        }
        check_cuda_error(cusparseScsrsv2_bufferSize(persistent_data->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, nz, persistent_data->descrU, persistent_data->d_val, persistent_data->d_row, persistent_data->d_col, persistent_data->infoU, &stmp));
        if (stmp > persistent_data->bufferSize) {
            persistent_data->bufferSize = stmp;
        }
        check_cuda_error(cudaMalloc(&(persistent_data->buffer), persistent_data->bufferSize));

        /* Perform analysis for ILU(0) */
        check_cuda_error(cusparseScsrilu02_analysis(persistent_data->cusparseHandle, N, nz, persistent_data->descr, persistent_data->d_val, persistent_data->d_row, persistent_data->d_col, persistent_data->infoILU, CUSPARSE_SOLVE_POLICY_USE_LEVEL, persistent_data->buffer));

        /* Copy A data to ILU(0) values as input*/
        check_cuda_error(cudaMemcpy(persistent_data->d_valsILU0, persistent_data->d_val, nz*sizeof(float), cudaMemcpyDeviceToDevice));

        /* generate the ILU(0) factors */
        check_cuda_error(cusparseScsrilu02(persistent_data->cusparseHandle, N, nz, persistent_data->descr, persistent_data->d_valsILU0, persistent_data->d_row, persistent_data->d_col, persistent_data->infoILU, CUSPARSE_SOLVE_POLICY_USE_LEVEL, persistent_data->buffer));
        /* perform triangular solve analysis */
        check_cuda_error(cusparseScsrsv2_analysis(persistent_data->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, nz, persistent_data->descrL, persistent_data->d_valsILU0, persistent_data->d_row, persistent_data->d_col, persistent_data->infoL, CUSPARSE_SOLVE_POLICY_USE_LEVEL, persistent_data->buffer));
        check_cuda_error(cusparseScsrsv2_analysis(persistent_data->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, nz, persistent_data->descrU, persistent_data->d_valsILU0, persistent_data->d_row, persistent_data->d_col, persistent_data->infoU, CUSPARSE_SOLVE_POLICY_USE_LEVEL, persistent_data->buffer));
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

    if (persistent_data->matA)  { check_cuda_error(cusparseDestroySpMat(persistent_data->matA)); }
    if (persistent_data->vecx)  { check_cuda_error(cusparseDestroyDnVec(persistent_data->vecx)); }
    if (persistent_data->vecAx) { check_cuda_error(cusparseDestroyDnVec(persistent_data->vecAx)); }
    if (persistent_data->vecp)  { check_cuda_error(cusparseDestroyDnVec(persistent_data->vecp)); }

    check_cuda_error(cudaFree(persistent_data->d_col));
    check_cuda_error(cudaFree(persistent_data->d_row));
    check_cuda_error(cudaFree(persistent_data->d_val));
    check_cuda_error(cudaFree(persistent_data->d_x));
    check_cuda_error(cudaFree(persistent_data->d_r));
    check_cuda_error(cudaFree(persistent_data->d_p));
    check_cuda_error(cudaFree(persistent_data->d_Ax));
    check_cuda_error(cudaFree(persistent_data->d_y));

    if(use_preconditioner) {
        /* Destroy descriptors */
        check_cuda_error(cusparseDestroyCsrsv2Info(persistent_data->infoU));
        check_cuda_error(cusparseDestroyCsrsv2Info(persistent_data->infoL));
        check_cuda_error(cusparseDestroyCsrilu02Info(persistent_data->infoILU));
        check_cuda_error(cusparseDestroyMatDescr(persistent_data->descrL));
        check_cuda_error(cusparseDestroyMatDescr(persistent_data->descrU));
        check_cuda_error(cusparseDestroyMatDescr(persistent_data->descr));
        check_cuda_error(cusparseDestroyDnVec(persistent_data->vecomega));

        check_cuda_error(cudaFree(persistent_data->d_valsILU0));
        check_cuda_error(cudaFree(persistent_data->d_zm1));
        check_cuda_error(cudaFree(persistent_data->d_zm2));
        check_cuda_error(cudaFree(persistent_data->d_rm2));
        check_cuda_error(cudaFree(persistent_data->d_omega));
    }

    free(persistent_data);
}

SOLVE_LINEAR_SYSTEM(gpu_conjugate_gradient) {

    /* Conjugate gradient.
       ------------------------------------------
       Follows the description by Golub & Van Loan, "Matrix Computations 3rd ed.", Section 10.2.6
    */
    struct gpu_persistent_data *persistent_data = (struct gpu_persistent_data *)config->persistent_data;

    if(!persistent_data) {
        log_error_and_exit("The gpu_conjugate_gradient solver needs to be initialized before being called. Add a init_function in the [linear_system_solver] section of the .ini file!\n");
    }

    float dot;

    float a, b, na, r0, r1;

    int k;
    float alpha, alpham1, beta;

    float *rhs; // Vector B
    rhs = (float *)malloc(sizeof(float) * num_active_cells);

    OMP(parallel for)
    for(uint32_t i = 0; i < num_active_cells; i++) {
        rhs[i] = active_cells[i]->b;
    }

    int N = persistent_data->N;
    int nz = persistent_data->nz;

    cudaMemcpy(persistent_data->d_r, rhs, N * sizeof(float), cudaMemcpyHostToDevice); // B

    alpha = 1.0f;
    alpham1 = -1.0f;
    beta = 0.0f;
    r0 = 0.0f;

    float numerator, denominator;

    check_cuda_error(cusparseSpMV(persistent_data->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, persistent_data->matA, persistent_data->vecx, &beta, persistent_data->vecAx, CUDA_R_32F, CUSPARSE_MV_ALG_DEFAULT, persistent_data->buffer));

    cublasSaxpy(persistent_data->cublasHandle, N, &alpham1, persistent_data->d_Ax, 1, persistent_data->d_r, 1);
    cublasSdot(persistent_data->cublasHandle, N, persistent_data->d_r, 1, persistent_data->d_r, 1, &r1);

    k = 1;

    while(r1 >= tol && k <= max_its) {

        if(use_preconditioner) {
            // preconditioner application: d_zm1 = U^-1 L^-1 d_r
            check_cuda_error(cusparseScsrsv2_solve(persistent_data->cusparseHandle,
                                                  CUSPARSE_OPERATION_NON_TRANSPOSE, N, nz, &floatone, persistent_data->descrL,
                                                 persistent_data->d_valsILU0, persistent_data->d_row,persistent_data-> d_col, persistent_data->infoL, persistent_data->d_r, persistent_data->d_y,
                                                  CUSPARSE_SOLVE_POLICY_USE_LEVEL, persistent_data->buffer));
            check_cuda_error(cusparseScsrsv2_solve(persistent_data->cusparseHandle,
                                                  CUSPARSE_OPERATION_NON_TRANSPOSE, N, nz, &floatone, persistent_data->descrU,
                                                  persistent_data->d_valsILU0, persistent_data->d_row, persistent_data->d_col, persistent_data->infoU, persistent_data->d_y, persistent_data->d_zm1,
                                                  CUSPARSE_SOLVE_POLICY_USE_LEVEL, persistent_data->buffer));
        }

        if(k > 1) {
            if(use_preconditioner) {
                cublasSdot(persistent_data->cublasHandle, N, persistent_data->d_r, 1, persistent_data->d_zm1, 1, &numerator);
                cublasSdot(persistent_data->cublasHandle, N, persistent_data->d_rm2, 1, persistent_data->d_zm2, 1, &denominator);
                b = numerator / denominator;
                cublasSscal(persistent_data->cublasHandle, N, &b, persistent_data->d_p, 1);
                cublasSaxpy(persistent_data->cublasHandle, N, &alpha, persistent_data->d_zm1, 1, persistent_data->d_p, 1);

            } else {
                b = r1 / r0;
                cublasSscal(persistent_data->cublasHandle, N, &b, persistent_data->d_p, 1);
                cublasSaxpy(persistent_data->cublasHandle, N, &alpha, persistent_data->d_r, 1, persistent_data->d_p, 1);
            }
        } else {
            if(use_preconditioner) {
                cublasScopy(persistent_data->cublasHandle, N, persistent_data->d_zm1, 1, persistent_data->d_p, 1);
            } else {
                cublasScopy(persistent_data->cublasHandle, N, persistent_data->d_r, 1, persistent_data->d_p, 1);
            }
        }

        if(use_preconditioner) {
            check_cuda_error(cusparseSpMV(persistent_data->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &floatone, persistent_data->matA, persistent_data->vecp, &floatzero, persistent_data->vecomega, CUDA_R_32F, CUSPARSE_MV_ALG_DEFAULT, persistent_data->buffer));
            cublasSdot(persistent_data->cublasHandle, N, persistent_data->d_r, 1, persistent_data->d_zm1, 1, &numerator);
            cublasSdot(persistent_data->cublasHandle, N, persistent_data->d_p, 1, persistent_data->d_omega, 1, &denominator);

            a = numerator / denominator;

            cublasSaxpy(persistent_data->cublasHandle, N, &a, persistent_data->d_p, 1, persistent_data->d_x, 1);
            cublasScopy(persistent_data->cublasHandle, N, persistent_data->d_r, 1, persistent_data->d_rm2, 1);
            cublasScopy(persistent_data->cublasHandle, N, persistent_data->d_zm1, 1, persistent_data->d_zm2, 1);
            na = -a;
            cublasSaxpy(persistent_data->cublasHandle, N, &na, persistent_data->d_omega, 1, persistent_data->d_r, 1);

            r0 = r1;
            cublasSdot(persistent_data->cublasHandle, N, persistent_data->d_r, 1, persistent_data->d_r, 1, &r1);

        } else {
            check_cuda_error(cusparseSpMV(persistent_data->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, persistent_data->matA, persistent_data->vecp, &beta, persistent_data->vecAx, CUDA_R_32F, CUSPARSE_MV_ALG_DEFAULT, persistent_data->buffer));
            cublasSdot(persistent_data->cublasHandle, N, persistent_data->d_p, 1, persistent_data->d_Ax, 1, &dot);

            cublasSdot(persistent_data->cublasHandle, N, persistent_data->d_p, 1, persistent_data->d_Ax, 1, &dot);
            a = r1 / dot;

            cublasSaxpy(persistent_data->cublasHandle, N, &a, persistent_data->d_p, 1, persistent_data->d_x, 1);
            na = -a;
            cublasSaxpy(persistent_data->cublasHandle, N, &na, persistent_data->d_Ax, 1, persistent_data->d_r, 1);

            r0 = r1;
            cublasSdot(persistent_data->cublasHandle, N, persistent_data->d_r, 1, persistent_data->d_r, 1, &r1);
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

    struct gpu_persistent_data *persistent_data = CALLOC_ONE_TYPE(struct gpu_persistent_data);

    int_array I = NULL, J = NULL;
    f32_array val = NULL;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, tol, config, "tolerance");
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, max_its, config, "max_iterations");

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

    float alpha = 1.0f;
    float beta = 0.0f;

    check_cuda_error(cusparseCreateCsr(&(persistent_data->matA), N, N, nz, persistent_data->d_row, persistent_data->d_col, persistent_data->d_val, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F));
    check_cuda_error(cusparseCreateDnVec(&(persistent_data->vecx), N, persistent_data->d_x, CUDA_R_32F));
    check_cuda_error(cusparseCreateDnVec(&(persistent_data->vecp), N, persistent_data->d_p, CUDA_R_32F));
    check_cuda_error(cusparseCreateDnVec(&(persistent_data->vecv), N, persistent_data->d_v, CUDA_R_32F));
    check_cuda_error(cusparseCreateDnVec(&(persistent_data->vecr), N, persistent_data->d_r, CUDA_R_32F));
    check_cuda_error(cusparseCreateDnVec(&(persistent_data->vect), N, persistent_data->d_t, CUDA_R_32F));
    check_cuda_error(cusparseCreateDnVec(&(persistent_data->vecAx), N, persistent_data->d_Ax, CUDA_R_32F));

    check_cuda_error(cusparseSpMV_bufferSize(persistent_data->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, persistent_data->matA, persistent_data->vecx, &beta, persistent_data->vecAx, CUDA_R_32F, CUSPARSE_MV_ALG_DEFAULT, &(persistent_data->bufferSize)));
    check_cuda_error(cudaMalloc(&(persistent_data->buffer), persistent_data->bufferSize));

    free(rhs);
    arrfree(I);
    arrfree(J);

    config->persistent_data = persistent_data;
}

SOLVE_LINEAR_SYSTEM(gpu_biconjugate_gradient) {

    struct gpu_persistent_data *persistent_data = (struct gpu_persistent_data *)config->persistent_data;

    if(!persistent_data) {
        log_error_and_exit("The gpu_biconjugate_gradient solver needs to be initialized before being called. Add a init_function in the [linear_system_solver] section of the .ini file\n");
    }

    float rho, rhop, beta, alpha, negalpha, omega, negomega, temp, temp2;
    float nrmr, nrmr0;
    rho = 0.0f;
    float zero = 0.0f;
    float one = 1.0f;
    float mone = -1.0f;

    float *rhs; // Vector B
    rhs = (float *)malloc(sizeof(float) * num_active_cells);

    OMP(parallel for)
    for(uint32_t i = 0; i < num_active_cells; i++) {
        rhs[i] = active_cells[i]->b;
    }

    int N = persistent_data->N;

    check_cuda_error(cudaMemcpy(persistent_data->d_r, rhs, N * sizeof(float), cudaMemcpyHostToDevice)); // B

    // compute initial residual r0=b-Ax0 (using initial guess in x)
    check_cublas_error(cusparseSpMV(persistent_data->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &one, persistent_data->matA, persistent_data->vecx, &zero, persistent_data->vecAx, CUDA_R_32F, CUSPARSE_MV_ALG_DEFAULT, persistent_data->buffer));

    check_cublas_error(cublasSaxpy(persistent_data->cublasHandle, N, &mone, persistent_data->d_Ax, 1, persistent_data->d_r, 1));

    // copy residual r into r^{\hat} and p
    check_cublas_error(cublasScopy(persistent_data->cublasHandle, N, persistent_data->d_r, 1, persistent_data->d_rw, 1));
    check_cublas_error(cublasScopy(persistent_data->cublasHandle, N, persistent_data->d_r, 1, persistent_data->d_p, 1));
    check_cublas_error(cublasSnrm2(persistent_data->cublasHandle, N, persistent_data->d_r, 1, &nrmr0));

    int i;
    for(i = 0; i < max_its;) {
        rhop = rho;
        check_cuda_error((cudaError_t)cublasSdot(persistent_data->cublasHandle, N, persistent_data->d_rw, 1, persistent_data->d_r, 1, &rho));

        if(i > 0) {
            beta = (rho / rhop) * (alpha / omega);
            negomega = -omega;
            check_cublas_error(cublasSaxpy(persistent_data->cublasHandle, N, &negomega, persistent_data->d_v, 1,persistent_data-> d_p, 1));
            check_cublas_error(cublasSscal(persistent_data->cublasHandle, N, &beta, persistent_data->d_p, 1));
            check_cublas_error(cublasSaxpy(persistent_data->cublasHandle, N, &one, persistent_data->d_r, 1, persistent_data->d_p, 1));
        }

        // matrix-vector multiplication
        check_cublas_error(cusparseSpMV(persistent_data->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &one, persistent_data->matA, persistent_data->vecp, &zero, persistent_data->vecv, CUDA_R_32F, CUSPARSE_MV_ALG_DEFAULT, persistent_data->buffer));
        check_cublas_error(cublasSdot(persistent_data->cublasHandle, N, persistent_data->d_rw, 1, persistent_data->d_v, 1, &temp));
        alpha = rho / temp;
        negalpha = -(alpha);
        check_cublas_error(cublasSaxpy(persistent_data->cublasHandle, N, &negalpha, persistent_data->d_v, 1, persistent_data->d_r, 1));
        check_cublas_error(cublasSaxpy(persistent_data->cublasHandle, N, &alpha, persistent_data->d_p, 1, persistent_data->d_x, 1));
        check_cublas_error(cublasSnrm2(persistent_data->cublasHandle, N, persistent_data->d_r, 1, &nrmr));

        if(nrmr < tol * nrmr0) {
            break;
        }

        // matrix-vector multiplication
        check_cublas_error(cusparseSpMV(persistent_data->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &one, persistent_data->matA, persistent_data->vecr, &zero, persistent_data->vect, CUDA_R_32F, CUSPARSE_MV_ALG_DEFAULT, persistent_data->buffer));

        check_cublas_error(cublasSdot(persistent_data->cublasHandle, N, persistent_data->d_t, 1, persistent_data->d_r, 1, &temp));
        check_cublas_error(cublasSdot(persistent_data->cublasHandle, N, persistent_data->d_t, 1, persistent_data->d_t, 1, &temp2));
        omega = temp / temp2;
        negomega = -(omega);
        check_cublas_error(cublasSaxpy(persistent_data->cublasHandle, N, &omega,    persistent_data->d_r, 1, persistent_data->d_x, 1));
        check_cublas_error(cublasSaxpy(persistent_data->cublasHandle, N, &negomega, persistent_data->d_t, 1, persistent_data->d_r, 1));

        check_cublas_error(cublasSnrm2(persistent_data->cublasHandle, N, persistent_data->d_r, 1, &nrmr));

        if(nrmr < tol * nrmr0) {
            i++;
            break;
        }
        i++;
    }

    check_cuda_error(cudaMemcpy(rhs, persistent_data->d_x, N * sizeof(float), cudaMemcpyDeviceToHost));

    *number_of_iterations = i - 1;
    *error = nrmr;

    OMP(parallel for)
    for(uint32_t j = 0; j < num_active_cells; j++) {
        active_cells[j]->v = rhs[j];
    }

    free(rhs);
}

END_LINEAR_SYSTEM(end_gpu_biconjugate_gradient) {

    struct gpu_persistent_data *persistent_data = (struct gpu_persistent_data *)config->persistent_data;

    if(!persistent_data) return;

    check_cuda_error((cudaError_t)cusparseDestroy(persistent_data->cusparseHandle));
    check_cuda_error((cudaError_t)cublasDestroy(persistent_data->cublasHandle));

    if (persistent_data->matA)  { check_cuda_error(cusparseDestroySpMat(persistent_data->matA)); }
    if (persistent_data->vecx)  { check_cuda_error(cusparseDestroyDnVec(persistent_data->vecx)); }
    if (persistent_data->vecAx) { check_cuda_error(cusparseDestroyDnVec(persistent_data->vecAx)); }
    if (persistent_data->vecp)  { check_cuda_error(cusparseDestroyDnVec(persistent_data->vecp)); }
    if (persistent_data->vecv)  { check_cuda_error(cusparseDestroyDnVec(persistent_data->vecv)); }
    if (persistent_data->vecr)  { check_cuda_error(cusparseDestroyDnVec(persistent_data->vecr)); }
    if (persistent_data->vect)  { check_cuda_error(cusparseDestroyDnVec(persistent_data->vect)); }

    check_cuda_error(cudaFree(persistent_data->d_col));
    check_cuda_error(cudaFree(persistent_data->d_row));
    check_cuda_error(cudaFree(persistent_data->d_val));
    check_cuda_error(cudaFree(persistent_data->d_x));
    check_cuda_error(cudaFree(persistent_data->d_r));
    check_cuda_error(cudaFree(persistent_data->d_rw));
    check_cuda_error(cudaFree(persistent_data->d_p));
    check_cuda_error(cudaFree(persistent_data->d_t));
    check_cuda_error(cudaFree(persistent_data->d_v));
    check_cuda_error(cudaFree(persistent_data->d_Ax));
    check_cuda_error(cudaFree(persistent_data->d_y));

    free(persistent_data);
}
