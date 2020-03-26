SOLVE_LINEAR_SYSTEM(gpu_conjugate_gradient) {

        /* Conjugate gradient without preconditioning.
           ------------------------------------------
           Follows the description by Golub & Van Loan, "Matrix Computations 3rd ed.", Section 10.2.6
        */
        real dot;

        real a, b, na, r0, r1;

        int k;
        real alpha, beta, alpham1;

        real *rhs; //Vector B

        uint32_t num_active_cells = the_grid->num_active_cells;
        struct cell_node** ac = the_grid->active_cells;

        rhs = (real*) malloc(sizeof(real)*num_active_cells);

        OMP(parallel for)
        for (uint32_t i = 0; i < num_active_cells; i++) {
            rhs[i] = ac[i]->b;
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

        OMP(parallel for)
        for (uint32_t i = 0; i < num_active_cells; i++) {
            ac[i]->v = rhs[i];
        }

        free(rhs);

}