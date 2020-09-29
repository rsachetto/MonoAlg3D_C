//
// Created by sachetto on 21/08/2020.
//
static int *d_col, *d_row;
static float *d_val, *d_x;

static float *d_r, *d_p, *d_Ax, *d_rw, *d_v, *d_t, *d_omega;

static int N = 0, nz = 0;

/* Get handle to the CUBLAS context */
static cublasHandle_t cublasHandle = 0;
static cublasStatus_t cublasStatus;

/* Get handle to the CUSPARSE context */
static cusparseHandle_t cusparseHandle = 0;

static cusparseSpMatDescr_t matA = NULL;
static cusparseDnVecDescr_t vecx = NULL;
static cusparseDnVecDescr_t vecp = NULL;
static cusparseDnVecDescr_t vecv = NULL;
static cusparseDnVecDescr_t vecr = NULL;
static cusparseDnVecDescr_t vect = NULL;
static cusparseDnVecDescr_t vecAx = NULL;
static cusparseDnVecDescr_t vecomega = NULL;

static int nzILU0;
static float *d_valsILU0;
static float *d_zm1, *d_zm2, *d_rm2;
static float *d_y;

static csrilu02Info_t infoILU = NULL;
cusparseMatDescr_t descrL = NULL;
csrsv2Info_t infoL = NULL;

cusparseMatDescr_t descrU = NULL;
csrsv2Info_t infoU = NULL;

cusparseMatDescr_t descr = NULL;

static size_t bufferSize = 0;
static void *buffer = NULL;
static size_t tmp = 0;
static int stmp = 0;

const float floatone = 1.0;
const float floatzero = 0.0;

INIT_LINEAR_SYSTEM(init_gpu_conjugate_gradient) {

    int_array I = NULL, J = NULL;
    f32_array val = NULL;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, tol, config->config_data, "tolerance");
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, max_its, config->config_data, "max_iterations");
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(use_preconditioner, config->config_data, "use_preconditioner");

    check_cublas_error(cublasCreate(&cublasHandle));
    check_cublas_error(cusparseCreate(&cusparseHandle));

    grid_to_csr(the_grid, &val, &I, &J);

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **active_cells = the_grid->active_cells;

    nz = arrlen(val);
    N = num_active_cells;

    check_cuda_error(cudaMalloc((void **)&d_col, nz * sizeof(int)));
    check_cuda_error(cudaMalloc((void **)&d_row, (N + 1) * sizeof(int)));
    check_cuda_error(cudaMalloc((void **)&d_val, nz * sizeof(float)));
    check_cuda_error(cudaMalloc((void **)&d_x, N * sizeof(float)));
    check_cuda_error(cudaMalloc((void **)&d_r, N * sizeof(float)));
    check_cuda_error(cudaMalloc((void **)&d_p, N * sizeof(float)));
    check_cuda_error(cudaMalloc((void **)&d_Ax, N * sizeof(float)));

    /* Wrap raw data into cuSPARSE generic API objects */
    check_cuda_error(cusparseCreateCsr(&matA, N, N, nz, d_row, d_col, d_val, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F));
    check_cuda_error(cusparseCreateDnVec(&vecx, N, d_x, CUDA_R_32F));
    check_cuda_error(cusparseCreateDnVec(&vecp, N, d_p, CUDA_R_32F));
    check_cuda_error(cusparseCreateDnVec(&vecAx, N, d_Ax, CUDA_R_32F));


    cudaMemcpy(d_col, J, nz * sizeof(int), cudaMemcpyHostToDevice);      // JA
    cudaMemcpy(d_row, I, (N + 1) * sizeof(int), cudaMemcpyHostToDevice); // IA
    cudaMemcpy(d_val, val, nz * sizeof(float), cudaMemcpyHostToDevice);  // A
    float *rhs = (float *)malloc(sizeof(float) * num_active_cells);

    OMP(parallel for)
    for(uint32_t i = 0; i < num_active_cells; i++) {
        rhs[i] = active_cells[i]->b;
    }

    check_cuda_error(cudaMemcpy(d_x, rhs, N * sizeof(float), cudaMemcpyHostToDevice)); // Result

    float alpha = 1.0;
    float beta = 0.0;

    check_cuda_error(cusparseSpMV_bufferSize(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecx, &beta, vecAx, CUDA_R_32F, CUSPARSE_MV_ALG_DEFAULT, &bufferSize));
    check_cuda_error(cudaMalloc(&buffer, bufferSize));

    if(use_preconditioner) {
        nzILU0 = 2 * N - 1;
        check_cuda_error(cudaMalloc((void **)&d_valsILU0, nz * sizeof(float)));
        check_cuda_error(cudaMalloc((void **)&d_zm1, (N) * sizeof(float)));
        check_cuda_error(cudaMalloc((void **)&d_zm2, (N) * sizeof(float)));
        check_cuda_error(cudaMalloc((void **)&d_rm2, (N) * sizeof(float)));
        check_cuda_error(cudaMalloc((void **)&d_omega, N * sizeof(float)));
        check_cuda_error(cudaMalloc((void **)&d_y, N * sizeof(float)));

        check_cuda_error(cusparseCreateDnVec(&vecomega, N, d_omega, CUDA_R_32F));

        check_cuda_error(cusparseCreateCsrilu02Info(&infoILU));

        /* Create L factor descriptor and triangular solve info */
        check_cuda_error(cusparseCreateMatDescr(&descrL));
        check_cuda_error(cusparseSetMatType(descrL, CUSPARSE_MATRIX_TYPE_GENERAL));
        check_cuda_error(cusparseSetMatIndexBase(descrL, CUSPARSE_INDEX_BASE_ZERO));
        check_cuda_error(cusparseSetMatFillMode(descrL, CUSPARSE_FILL_MODE_LOWER));
        check_cuda_error(cusparseSetMatDiagType(descrL, CUSPARSE_DIAG_TYPE_UNIT));
        check_cuda_error(cusparseCreateCsrsv2Info(&infoL));

        /* Create U factor descriptor and triangular solve info */
        check_cuda_error(cusparseCreateMatDescr(&descrU));
        check_cuda_error(cusparseSetMatType(descrU, CUSPARSE_MATRIX_TYPE_GENERAL));
        check_cuda_error(cusparseSetMatIndexBase(descrU, CUSPARSE_INDEX_BASE_ZERO));
        check_cuda_error(cusparseSetMatFillMode(descrU, CUSPARSE_FILL_MODE_UPPER));
        check_cuda_error(cusparseSetMatDiagType(descrU, CUSPARSE_DIAG_TYPE_NON_UNIT));
        check_cuda_error(cusparseCreateCsrsv2Info(&infoU));

        check_cuda_error(cusparseCreateMatDescr(&descr));
        check_cuda_error(cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL));
        check_cuda_error(cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO));

        check_cuda_error(cusparseSpMV_bufferSize(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &floatone, matA, vecp, &floatzero, vecomega, CUDA_R_32F, CUSPARSE_MV_ALG_DEFAULT, &tmp));
        if (tmp > bufferSize) {
            bufferSize = stmp;
        }
        check_cuda_error(cusparseScsrilu02_bufferSize(cusparseHandle, N, nz, descr, d_val, d_row, d_col, infoILU, &stmp));
        if (stmp > bufferSize) {
            bufferSize = stmp;
        }
        check_cuda_error(cusparseScsrsv2_bufferSize(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, nz, descrL, d_val, d_row, d_col, infoL, &stmp));
        if (stmp > bufferSize) {
            bufferSize = stmp;
        }
        check_cuda_error(cusparseScsrsv2_bufferSize(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, nz, descrU, d_val,d_row, d_col, infoU, &stmp));
        if (stmp > bufferSize) {
            bufferSize = stmp;
        }
        check_cuda_error(cudaMalloc(&buffer, bufferSize));

        /* Perform analysis for ILU(0) */
        check_cuda_error(cusparseScsrilu02_analysis(cusparseHandle, N, nz, descr, d_val, d_row, d_col, infoILU, CUSPARSE_SOLVE_POLICY_USE_LEVEL, buffer));

        /* Copy A data to ILU(0) values as input*/
        check_cuda_error(cudaMemcpy(d_valsILU0, d_val, nz*sizeof(float), cudaMemcpyDeviceToDevice));

        /* generate the ILU(0) factors */
        check_cuda_error(cusparseScsrilu02(cusparseHandle, N, nz, descr, d_valsILU0, d_row, d_col, infoILU, CUSPARSE_SOLVE_POLICY_USE_LEVEL, buffer));
        /* perform triangular solve analysis */
        check_cuda_error(cusparseScsrsv2_analysis(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, nz, descrL, d_valsILU0, d_row, d_col, infoL, CUSPARSE_SOLVE_POLICY_USE_LEVEL, buffer));
        check_cuda_error(cusparseScsrsv2_analysis(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, nz, descrU, d_valsILU0, d_row, d_col, infoU, CUSPARSE_SOLVE_POLICY_USE_LEVEL, buffer));
    }

    free(rhs);
    arrfree(I);
    arrfree(J);
    arrfree(val);
}

END_LINEAR_SYSTEM(end_gpu_conjugate_gradient) {

    check_cuda_error((cudaError_t)cusparseDestroy(cusparseHandle));
    check_cuda_error((cudaError_t)cublasDestroy(cublasHandle));

    if (matA)  { check_cuda_error(cusparseDestroySpMat(matA)); }
    if (vecx)  { check_cuda_error(cusparseDestroyDnVec(vecx)); }
    if (vecAx) { check_cuda_error(cusparseDestroyDnVec(vecAx)); }
    if (vecp)  { check_cuda_error(cusparseDestroyDnVec(vecp)); }

    check_cuda_error(cudaFree(d_col));
    check_cuda_error(cudaFree(d_row));
    check_cuda_error(cudaFree(d_val));
    check_cuda_error(cudaFree(d_x));
    check_cuda_error(cudaFree(d_r));
    check_cuda_error(cudaFree(d_p));
    check_cuda_error(cudaFree(d_Ax));
    check_cuda_error(cudaFree(d_y));

    if(use_preconditioner) {
        /* Destroy descriptors */
        check_cuda_error(cusparseDestroyCsrsv2Info(infoU));
        check_cuda_error(cusparseDestroyCsrsv2Info(infoL));
        check_cuda_error(cusparseDestroyCsrilu02Info(infoILU));
        check_cuda_error(cusparseDestroyMatDescr(descrL));
        check_cuda_error(cusparseDestroyMatDescr(descrU));
        check_cuda_error(cusparseDestroyMatDescr(descr));
        check_cuda_error(cusparseDestroyDnVec(vecomega));

        check_cuda_error(cudaFree(d_valsILU0));
        check_cuda_error(cudaFree(d_zm1));
        check_cuda_error(cudaFree(d_zm2));
        check_cuda_error(cudaFree(d_rm2));
        check_cuda_error(cudaFree(d_omega));
    }
}

SOLVE_LINEAR_SYSTEM(gpu_conjugate_gradient) {

    /* Conjugate gradient.
       ------------------------------------------
       Follows the description by Golub & Van Loan, "Matrix Computations 3rd ed.", Section 10.2.6
    */
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

    cudaMemcpy(d_r, rhs, N * sizeof(float), cudaMemcpyHostToDevice); // B

    alpha = 1.0;
    alpham1 = -1.0;
    beta = 0.0;
    r0 = 0.;

    float numerator, denominator;

    check_cuda_error(cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecx, &beta, vecAx, CUDA_R_32F, CUSPARSE_MV_ALG_DEFAULT, buffer));

    cublasSaxpy(cublasHandle, N, &alpham1, d_Ax, 1, d_r, 1);
    cublasStatus = cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);

    k = 1;

    while(r1 >= tol && k <= max_its) {

        if(use_preconditioner) {
            // preconditioner application: d_zm1 = U^-1 L^-1 d_r
            check_cuda_error(cusparseScsrsv2_solve(cusparseHandle,
                                                  CUSPARSE_OPERATION_NON_TRANSPOSE, N, nz, &floatone, descrL,
                                                  d_valsILU0, d_row, d_col, infoL, d_r, d_y,
                                                  CUSPARSE_SOLVE_POLICY_USE_LEVEL, buffer));
            check_cuda_error(cusparseScsrsv2_solve(cusparseHandle,
                                                  CUSPARSE_OPERATION_NON_TRANSPOSE, N, nz, &floatone, descrU,
                                                  d_valsILU0, d_row, d_col, infoU, d_y, d_zm1,
                                                  CUSPARSE_SOLVE_POLICY_USE_LEVEL, buffer));
        }

        if(k > 1) {
            if(use_preconditioner) {
                cublasSdot(cublasHandle, N, d_r, 1, d_zm1, 1, &numerator);
                cublasSdot(cublasHandle, N, d_rm2, 1, d_zm2, 1, &denominator);
                b = numerator / denominator;
                cublasSscal(cublasHandle, N, &b, d_p, 1);
                cublasSaxpy(cublasHandle, N, &alpha, d_zm1, 1, d_p, 1);

            } else {
                b = r1 / r0;
                cublasStatus = cublasSscal(cublasHandle, N, &b, d_p, 1);
                cublasStatus = cublasSaxpy(cublasHandle, N, &alpha, d_r, 1, d_p, 1);
            }
        } else {
            if(use_preconditioner) {
                cublasScopy(cublasHandle, N, d_zm1, 1, d_p, 1);
            } else {
                cublasStatus = cublasScopy(cublasHandle, N, d_r, 1, d_p, 1);
            }
        }

        if(use_preconditioner) {
            check_cuda_error(cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &floatone, matA, vecp, &floatzero, vecomega, CUDA_R_32F, CUSPARSE_MV_ALG_DEFAULT, buffer));
            cublasSdot(cublasHandle, N, d_r, 1, d_zm1, 1, &numerator);
            cublasSdot(cublasHandle, N, d_p, 1, d_omega, 1, &denominator);

            a = numerator / denominator;

            cublasSaxpy(cublasHandle, N, &a, d_p, 1, d_x, 1);
            cublasScopy(cublasHandle, N, d_r, 1, d_rm2, 1);
            cublasScopy(cublasHandle, N, d_zm1, 1, d_zm2, 1);
            na = -a;
            cublasSaxpy(cublasHandle, N, &na, d_omega, 1, d_r, 1);

            r0 = r1;
            cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);

        } else {
            check_cuda_error(cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecp, &beta, vecAx, CUDA_R_32F, CUSPARSE_MV_ALG_DEFAULT, buffer));
            cublasStatus = cublasSdot(cublasHandle, N, d_p, 1, d_Ax, 1, &dot);

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

    cudaMemcpy(rhs, d_x, N * sizeof(float), cudaMemcpyDeviceToHost);

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
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, tol, config->config_data, "tolerance");
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, max_its, config->config_data, "max_iterations");
    // GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(use_preconditioner, config->config_data, "use_preconditioner");

    check_cuda_error((cudaError_t)cublasCreate(&cublasHandle));

    check_cuda_error((cudaError_t)cusparseCreate(&cusparseHandle));

    check_cuda_error((cudaError_t)cusparseCreateMatDescr(&descr));

    cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);

    grid_to_csr(the_grid, &val, &I, &J);

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **active_cells = the_grid->active_cells;

    nz = arrlen(val);
    N = num_active_cells;

    check_cuda_error(cudaMalloc((void **)&d_col, nz * sizeof(int)));
    check_cuda_error(cudaMalloc((void **)&d_row, (N + 1) * sizeof(int)));
    check_cuda_error(cudaMalloc((void **)&d_val, nz * sizeof(float)));
    check_cuda_error(cudaMalloc((void **)&d_x, N * sizeof(float)));
    check_cuda_error(cudaMalloc((void **)&d_r, N * sizeof(float)));
    check_cuda_error(cudaMalloc((void **)&d_rw, N * sizeof(float)));
    check_cuda_error(cudaMalloc((void **)&d_p, N * sizeof(float)));
    check_cuda_error(cudaMalloc((void **)&d_v, N * sizeof(float)));
    check_cuda_error(cudaMalloc((void **)&d_t, N * sizeof(float)));
    check_cuda_error(cudaMalloc((void **)&d_Ax, N * sizeof(float)));

    cudaMemcpy(d_col, J, nz * sizeof(int), cudaMemcpyHostToDevice);      // JA
    cudaMemcpy(d_row, I, (N + 1) * sizeof(int), cudaMemcpyHostToDevice); // IA
    cudaMemcpy(d_val, val, nz * sizeof(float), cudaMemcpyHostToDevice);  // A

    check_cuda_error(cudaMemset(d_r, 0, N * sizeof(float)));
    check_cuda_error(cudaMemset(d_rw,0, N * sizeof(float)));
    check_cuda_error(cudaMemset(d_p, 0, N * sizeof(float)));
    check_cuda_error(cudaMemset(d_t, 0, N * sizeof(float)));
    check_cuda_error(cudaMemset(d_v, 0, N * sizeof(float)));

    float *rhs = (float *)malloc(sizeof(float) * num_active_cells);

    OMP(parallel for)
    for(uint32_t i = 0; i < num_active_cells; i++) {
        rhs[i] = active_cells[i]->b;
    }

    check_cuda_error(cudaMemcpy(d_x, rhs, N * sizeof(float), cudaMemcpyHostToDevice)); // Result

    float alpha = 1.0;
    float beta = 0.0;

 	check_cuda_error(cusparseCreateCsr(&matA, N, N, nz, d_row, d_col, d_val, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F));
    check_cuda_error(cusparseCreateDnVec(&vecx, N, d_x, CUDA_R_32F));
    check_cuda_error(cusparseCreateDnVec(&vecp, N, d_p, CUDA_R_32F));
    check_cuda_error(cusparseCreateDnVec(&vecv, N, d_v, CUDA_R_32F));
    check_cuda_error(cusparseCreateDnVec(&vecr, N, d_r, CUDA_R_32F));
    check_cuda_error(cusparseCreateDnVec(&vect, N, d_t, CUDA_R_32F));
    check_cuda_error(cusparseCreateDnVec(&vecAx, N, d_Ax, CUDA_R_32F));

    check_cuda_error(cusparseSpMV_bufferSize(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecx, &beta, vecAx, CUDA_R_32F, CUSPARSE_MV_ALG_DEFAULT, &bufferSize));
    check_cuda_error(cudaMalloc(&buffer, bufferSize));

    free(rhs);
    arrfree(I);
    arrfree(J);

}

SOLVE_LINEAR_SYSTEM(gpu_biconjugate_gradient) {

    float rho, rhop, beta, alpha, negalpha, omega, negomega, temp, temp2;
    float nrmr, nrmr0;
    rho = 0.0;
    float zero = 0.0;
    float one = 1.0;
    float mone = -1.0;
    int i = 0;

    float *rhs; // Vector B
    rhs = (float *)malloc(sizeof(float) * num_active_cells);

    OMP(parallel for)
    for(uint32_t i = 0; i < num_active_cells; i++) {
        rhs[i] = active_cells[i]->b;
    }

    check_cuda_error(cudaMemcpy(d_r, rhs, N * sizeof(float), cudaMemcpyHostToDevice)); // B

    // compute initial residual r0=b-Ax0 (using initial guess in x)
	check_cublas_error(cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &one, matA, vecx, &zero, vecAx, CUDA_R_32F, CUSPARSE_MV_ALG_DEFAULT, buffer));
    
	check_cublas_error(cublasSaxpy(cublasHandle, N, &mone, d_Ax, 1, d_r, 1));

    // copy residual r into r^{\hat} and p
    check_cublas_error(cublasScopy(cublasHandle, N, d_r, 1, d_rw, 1));
    check_cublas_error(cublasScopy(cublasHandle, N, d_r, 1, d_p, 1));
    check_cublas_error(cublasSnrm2(cublasHandle, N, d_r, 1, &nrmr0));

    for(i = 0; i < max_its;) {
        rhop = rho;
        check_cuda_error((cudaError_t)cublasSdot(cublasHandle, N, d_rw, 1, d_r, 1, &rho));

        if(i > 0) {
            beta = (rho / rhop) * (alpha / omega);
            negomega = -omega;
            check_cublas_error(cublasSaxpy(cublasHandle, N, &negomega, d_v, 1, d_p, 1));
            check_cublas_error(cublasSscal(cublasHandle, N, &beta, d_p, 1));
            check_cublas_error(cublasSaxpy(cublasHandle, N, &one, d_r, 1, d_p, 1));
        }

        // matrix-vector multiplication
		check_cublas_error(cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &one, matA, vecp, &zero, vecv, CUDA_R_32F, CUSPARSE_MV_ALG_DEFAULT, buffer));
        check_cublas_error(cublasSdot(cublasHandle, N, d_rw, 1, d_v, 1, &temp));
        alpha = rho / temp;
        negalpha = -(alpha);
        check_cublas_error(cublasSaxpy(cublasHandle, N, &negalpha, d_v, 1, d_r, 1));
        check_cublas_error(cublasSaxpy(cublasHandle, N, &alpha, d_p, 1, d_x, 1));
        check_cublas_error(cublasSnrm2(cublasHandle, N, d_r, 1, &nrmr));

        if(nrmr < tol * nrmr0) {
            break;
        }

        // matrix-vector multiplication
		check_cublas_error(cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &one, matA, vecr, &zero, vect, CUDA_R_32F, CUSPARSE_MV_ALG_DEFAULT, buffer));

        check_cublas_error(cublasSdot(cublasHandle, N, d_t, 1, d_r, 1, &temp));
        check_cublas_error(cublasSdot(cublasHandle, N, d_t, 1, d_t, 1, &temp2));
        omega = temp / temp2;
        negomega = -(omega);
        check_cublas_error(cublasSaxpy(cublasHandle, N, &omega, d_r, 1, d_x, 1));
        check_cublas_error(cublasSaxpy(cublasHandle, N, &negomega, d_t, 1, d_r, 1));

        check_cublas_error(cublasSnrm2(cublasHandle, N, d_r, 1, &nrmr));

        if(nrmr < tol * nrmr0) {
            i++;
            break;
        }
        i++;
    }

    check_cuda_error(cudaMemcpy(rhs, d_x, N * sizeof(float), cudaMemcpyDeviceToHost));

    *number_of_iterations = i - 1;
    *error = nrmr;

    OMP(parallel for)
    for(uint32_t i = 0; i < num_active_cells; i++) {
        active_cells[i]->v = rhs[i];
    }

    free(rhs);
}

END_LINEAR_SYSTEM(end_gpu_biconjugate_gradient) {

    check_cuda_error((cudaError_t)cusparseDestroy(cusparseHandle));
    check_cuda_error((cudaError_t)cublasDestroy(cublasHandle));

    if (matA)  { check_cuda_error(cusparseDestroySpMat(matA)); }
    if (vecx)  { check_cuda_error(cusparseDestroyDnVec(vecx)); }
    if (vecAx) { check_cuda_error(cusparseDestroyDnVec(vecAx)); }
    if (vecp)  { check_cuda_error(cusparseDestroyDnVec(vecp)); }
    if (vecv)  { check_cuda_error(cusparseDestroyDnVec(vecv)); }
    if (vecr)  { check_cuda_error(cusparseDestroyDnVec(vecr)); }
    if (vect)  { check_cuda_error(cusparseDestroyDnVec(vect)); }

    check_cuda_error(cudaFree(d_col));
    check_cuda_error(cudaFree(d_row));
    check_cuda_error(cudaFree(d_val));
    check_cuda_error(cudaFree(d_x));
    check_cuda_error(cudaFree(d_r));
    check_cuda_error(cudaFree(d_rw));
    check_cuda_error(cudaFree(d_p));
    check_cuda_error(cudaFree(d_t));
    check_cuda_error(cudaFree(d_v));
    check_cuda_error(cudaFree(d_Ax));
    check_cuda_error(cudaFree(d_y));
}
