//
// Created by sachetto on 21/08/2020.
//
static int *d_col, *d_row;
static float *d_val, *d_x;

static float *d_r, *d_p, *d_Ax, *d_rw, *d_v, *d_t;

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

	cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);

	//grid_to_csr(the_grid, &val, &I, &J);

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

	cudaMemcpy(d_col, J, nz * sizeof(int), cudaMemcpyHostToDevice);      // JA
	cudaMemcpy(d_row, I, (N + 1) * sizeof(int), cudaMemcpyHostToDevice); // IA
	cudaMemcpy(d_val, val, nz * sizeof(float), cudaMemcpyHostToDevice);  // A
	float *rhs = (float *)malloc(sizeof(float) * num_active_cells);

	OMP(parallel for)
		for(uint32_t i = 0; i < num_active_cells; i++) {
			rhs[i] = active_cells[i]->b;
		}

	check_cuda_error(cudaMemcpy(d_x, rhs, N * sizeof(float), cudaMemcpyHostToDevice)); // Result

	if(use_preconditioner) {
		nzILU0 = 2 * N - 1;
		check_cuda_error(cudaMalloc((void **)&d_valsILU0, nz * sizeof(float)));
		check_cuda_error(cudaMalloc((void **)&d_zm1, (N) * sizeof(float)));
		check_cuda_error(cudaMalloc((void **)&d_zm2, (N) * sizeof(float)));
		check_cuda_error(cudaMalloc((void **)&d_rm2, (N) * sizeof(float)));
		check_cuda_error(cudaMalloc((void **)&d_y, N * sizeof(float)));

		cusparseStatus = cusparseCreateSolveAnalysisInfo(&infoA);
		check_cuda_error((cudaError_t)cusparseStatus);

		/* Perform the analysis for the Non-Transpose case */
		cusparseStatus = cusparseScsrsv_analysis(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, nz, descr, d_val, d_row, d_col, infoA);

		check_cuda_error((cudaError_t)cusparseStatus);

		/* Copy A data to ILU0 vals as input*/
		cudaMemcpy(d_valsILU0, d_val, nz * sizeof(float), cudaMemcpyDeviceToDevice);

		/* generate the Incomplete LU factor H for the matrix A using cudsparseScsrilu0 */
		cusparseStatus = cusparseScsrilu0(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, descr, d_valsILU0, d_row, d_col, infoA);

		check_cuda_error((cudaError_t)cusparseStatus);

		cusparseCreateSolveAnalysisInfo(&info_u);

		cusparseStatus = cusparseCreateMatDescr(&descrL);
		cusparseSetMatType(descrL, CUSPARSE_MATRIX_TYPE_GENERAL);
		cusparseSetMatIndexBase(descrL, CUSPARSE_INDEX_BASE_ZERO);
		cusparseSetMatFillMode(descrL, CUSPARSE_FILL_MODE_LOWER);
		cusparseSetMatDiagType(descrL, CUSPARSE_DIAG_TYPE_UNIT);

		cusparseStatus = cusparseCreateMatDescr(&descrU);
		cusparseSetMatType(descrU, CUSPARSE_MATRIX_TYPE_GENERAL);
		cusparseSetMatIndexBase(descrU, CUSPARSE_INDEX_BASE_ZERO);
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

    check_cuda_error((cudaError_t)cusparseDestroy(cusparseHandle));
    check_cuda_error((cudaError_t)cublasDestroy(cublasHandle));
    check_cuda_error((cudaError_t)cusparseDestroyMatDescr(descr));

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
        check_cuda_error((cudaError_t)cusparseDestroyMatDescr(descrL));
        check_cuda_error((cudaError_t)cusparseDestroyMatDescr(descrU));
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

    cusparseScsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &alpha, descr, d_val, d_row, d_col, d_x, &beta, d_Ax);

    cublasSaxpy(cublasHandle, N, &alpham1, d_Ax, 1, d_r, 1);
    cublasStatus = cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);

    k = 1;

    while(r1 >= tol && k <= max_its) {

        if(use_preconditioner) {
            // Forward Solve, we can re-use infoA since the sparsity pattern of A matches that of L
            cusparseStatus =
                cusparseScsrsv_solve(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, &alpha, descrL, d_valsILU0, d_row, d_col, infoA, d_r, d_y);
            check_cuda_error((cudaError_t)cusparseStatus);

            // Back Substitution
            cusparseStatus =
                cusparseScsrsv_solve(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, &alpha, descrU, d_valsILU0, d_row, d_col, info_u, d_y, d_zm1);
            check_cuda_error((cudaError_t)cusparseStatus);
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
            cusparseScsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nzILU0, &alpha, descrU, d_val, d_row, d_col, d_p, &beta, d_Ax);
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

        } else {
            cusparseScsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &alpha, descr, d_val, d_row, d_col, d_p, &beta, d_Ax);

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

    grid_to_csr(the_grid, &val, &I, &J, false);

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

    cudaMemcpy(d_r, rhs, N * sizeof(float), cudaMemcpyHostToDevice); // B

    // compute initial residual r0=b-Ax0 (using initial guess in x)
    check_cuda_error(
        (cudaError_t)cusparseScsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &one, descr, d_val, d_row, d_col, d_x, &zero, d_Ax));
    check_cuda_error((cudaError_t)cublasSaxpy(cublasHandle, N, &mone, d_Ax, 1, d_r, 1));

    // copy residual r into r^{\hat} and p
    check_cuda_error((cudaError_t)cublasScopy(cublasHandle, N, d_r, 1, d_rw, 1));
    check_cuda_error((cudaError_t)cublasScopy(cublasHandle, N, d_r, 1, d_p, 1));
    check_cuda_error((cudaError_t)cublasSnrm2(cublasHandle, N, d_r, 1, &nrmr0));


    for(i = 0; i < max_its;) {
        rhop = rho;
        check_cuda_error((cudaError_t)cublasSdot(cublasHandle, N, d_rw, 1, d_r, 1, &rho));

        if(i > 0) {
            beta = (rho / rhop) * (alpha / omega);
            negomega = -omega;
            check_cuda_error((cudaError_t)cublasSaxpy(cublasHandle, N, &negomega, d_v, 1, d_p, 1));
            check_cuda_error((cudaError_t)cublasSscal(cublasHandle, N, &beta, d_p, 1));
            check_cuda_error((cudaError_t)cublasSaxpy(cublasHandle, N, &one, d_r, 1, d_p, 1));
        }

        // matrix-vector multiplication
        check_cuda_error(
            (cudaError_t)cusparseScsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &one, descr, d_val, d_row, d_col, d_p, &zero, d_v));
        check_cuda_error((cudaError_t)cublasSdot(cublasHandle, N, d_rw, 1, d_v, 1, &temp));
        alpha = rho / temp;
        negalpha = -(alpha);
        check_cuda_error((cudaError_t)cublasSaxpy(cublasHandle, N, &negalpha, d_v, 1, d_r, 1));
        check_cuda_error((cudaError_t)cublasSaxpy(cublasHandle, N, &alpha, d_p, 1, d_x, 1));
        check_cuda_error((cudaError_t)cublasSnrm2(cublasHandle, N, d_r, 1, &nrmr));

        if(nrmr < tol * nrmr0) {
            break;
        }

        // matrix-vector multiplication
        check_cuda_error((cudaError_t)cusparseScsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &one, descr, d_val, d_row, d_col, d_r, &zero, d_t));

        check_cuda_error((cudaError_t)cublasSdot(cublasHandle, N, d_t, 1, d_r, 1, &temp));
        check_cuda_error((cudaError_t)cublasSdot(cublasHandle, N, d_t, 1, d_t, 1, &temp2));
        omega = temp / temp2;
        negomega = -(omega);
        check_cuda_error((cudaError_t)cublasSaxpy(cublasHandle, N, &omega, d_r, 1, d_x, 1));
        check_cuda_error((cudaError_t)cublasSaxpy(cublasHandle, N, &negomega, d_t, 1, d_r, 1));

        check_cuda_error((cudaError_t)cublasSnrm2(cublasHandle, N, d_r, 1, &nrmr));

        if(nrmr < tol * nrmr0) {
            i++;
            break;
        }
        i++;
    }

    cudaMemcpy(rhs, d_x, N * sizeof(float), cudaMemcpyDeviceToHost);

    *number_of_iterations = i - 1;
    *error = nrmr;

    OMP(parallel for)
    for(uint32_t i = 0; i < num_active_cells; i++) {
        active_cells[i]->v = rhs[i];
    }

    free(rhs);
}
