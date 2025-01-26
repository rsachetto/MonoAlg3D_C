#include <sycl/sycl.hpp>
#include <dpct/dpct.hpp>
#include <dpct/blas_utils.hpp>
//
// Created by sachetto on 21/08/2020.
//

struct gpu_persistent_data {
    int *d_col, *d_row;
    float *d_val, *d_x;

    float *d_r, *d_p, *d_Ax, *d_rw, *d_v, *d_t, *d_omega;

    int N, nz;

    /* Get handle to the CUBLAS context */
    dpct::blas::descriptor_ptr cublasHandle;

    /* Get handle to the CUSPARSE context */
    dpct::sparse::descriptor_ptr cusparseHandle;

    dpct::sparse::sparse_matrix_desc_t matA;
    std::shared_ptr<dpct::sparse::dense_vector_desc> vecx;
    std::shared_ptr<dpct::sparse::dense_vector_desc> vecp;
    std::shared_ptr<dpct::sparse::dense_vector_desc> vecv;
    std::shared_ptr<dpct::sparse::dense_vector_desc> vecr;
    std::shared_ptr<dpct::sparse::dense_vector_desc> vect;
    std::shared_ptr<dpct::sparse::dense_vector_desc> vecAx;
    std::shared_ptr<dpct::sparse::dense_vector_desc> vecomega;

    float *d_valsILU0;
    float *d_zm1, *d_zm2, *d_rm2;
    float *d_y;
    
    int                 bufferSizeLU;
    size_t              bufferSizeMV, bufferSizeL, bufferSizeU;
    void*               d_bufferLU, *d_bufferMV,  *d_bufferL, *d_bufferU;
    int spsvDescrL, spsvDescrU;
    std::shared_ptr<dpct::sparse::matrix_info> matLU;
    csrilu02Info_t      infoILU;

    dpct::sparse::sparse_matrix_desc_t matM_lower, matM_upper;
};

const float floatone = 1.0f;
const float floatzero = 0.0f;

INIT_LINEAR_SYSTEM(init_gpu_conjugate_gradient) {

    if(the_grid->adaptive) {
        log_error_and_exit("The gpu conjugate gradient does not support mesh adaptivity. Aborting!\n");
    }

    struct gpu_persistent_data *persistent_data = CALLOC_ONE_TYPE(struct gpu_persistent_data);

    int_array I = NULL, J = NULL;
    f32_array val = NULL;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, tol, config, "tolerance");
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, max_its, config, "max_iterations");

    check_cublas_error(DPCT_CHECK_ERROR(persistent_data->cublasHandle = new dpct::blas::descriptor()));
    check_cublas_error(DPCT_CHECK_ERROR(persistent_data->cusparseHandle = new dpct::sparse::descriptor()));

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

    check_cuda_error(DPCT_CHECK_ERROR(persistent_data->d_col = sycl::malloc_device<int>(nz, dpct::get_in_order_queue())));
    check_cuda_error(DPCT_CHECK_ERROR(persistent_data->d_row = sycl::malloc_device<int>((N + 1), dpct::get_in_order_queue())));
    check_cuda_error(DPCT_CHECK_ERROR(persistent_data->d_val = sycl::malloc_device<float>(nz, dpct::get_in_order_queue())));
    check_cuda_error(DPCT_CHECK_ERROR(persistent_data->d_x = sycl::malloc_device<float>(N, dpct::get_in_order_queue())));
    check_cuda_error(DPCT_CHECK_ERROR(persistent_data->d_r = sycl::malloc_device<float>(N, dpct::get_in_order_queue())));
    check_cuda_error(DPCT_CHECK_ERROR(persistent_data->d_p = sycl::malloc_device<float>(N, dpct::get_in_order_queue())));
    check_cuda_error(DPCT_CHECK_ERROR(persistent_data->d_Ax = sycl::malloc_device<float>(N, dpct::get_in_order_queue())));

    /* Wrap raw data into cuSPARSE generic API objects */
    check_cuda_error(DPCT_CHECK_ERROR(persistent_data->matA = std::make_shared<dpct::sparse::sparse_matrix_desc>(
                                          N, N, nz, persistent_data->d_row, persistent_data->d_col, persistent_data->d_val, dpct::library_data_t::real_int32,
                                          dpct::library_data_t::real_int32, oneapi::mkl::index_base::zero, dpct::library_data_t::real_float,
                                          dpct::sparse::matrix_format::csr)));
    check_cuda_error(
        DPCT_CHECK_ERROR(persistent_data->vecx = std::make_shared<dpct::sparse::dense_vector_desc>(N, persistent_data->d_x, dpct::library_data_t::real_float)));
    check_cuda_error(
        DPCT_CHECK_ERROR(persistent_data->vecp = std::make_shared<dpct::sparse::dense_vector_desc>(N, persistent_data->d_p, dpct::library_data_t::real_float)));
    check_cuda_error(DPCT_CHECK_ERROR(persistent_data->vecAx =
                                          std::make_shared<dpct::sparse::dense_vector_desc>(N, persistent_data->d_Ax, dpct::library_data_t::real_float)));

    dpct::get_in_order_queue().memcpy(persistent_data->d_col, J, nz * sizeof(int));            // JA
    dpct::get_in_order_queue().memcpy(persistent_data->d_row, I, (N + 1) * sizeof(int));       // IA
    dpct::get_in_order_queue().memcpy(persistent_data->d_val, val, nz * sizeof(float)).wait(); // A
    float *rhs = (float *)malloc(sizeof(float) * num_active_cells);

    OMP(parallel for)
    for(uint32_t i = 0; i < num_active_cells; i++) {
        rhs[i] = active_cells[i]->b;
    }

    check_cuda_error(DPCT_CHECK_ERROR(dpct::get_in_order_queue().memcpy(persistent_data->d_x, rhs, N * sizeof(float)).wait())); // Result

    float alpha = 1.0f;
    float beta = 0.0f;

    check_cuda_error(DPCT_CHECK_ERROR(persistent_data->bufferSizeMV = 0));
    check_cuda_error(DPCT_CHECK_ERROR(persistent_data->d_bufferMV = (void *)sycl::malloc_device(persistent_data->bufferSizeMV, dpct::get_in_order_queue())));

    config->persistent_data = persistent_data;
    free(rhs);
    arrfree(I);
    arrfree(J);
    arrfree(val);
}

END_LINEAR_SYSTEM(end_gpu_conjugate_gradient) {

    struct gpu_persistent_data *persistent_data = (struct gpu_persistent_data *)config->persistent_data;

    if(!persistent_data) return;

    check_cuda_error((dpct::err0)DPCT_CHECK_ERROR(delete(persistent_data->cusparseHandle)));
    check_cuda_error((dpct::err0)DPCT_CHECK_ERROR(delete(persistent_data->cublasHandle)));

    if(persistent_data->matA) { check_cuda_error(DPCT_CHECK_ERROR((persistent_data->matA).reset())); }
    if(persistent_data->vecx) { check_cuda_error(DPCT_CHECK_ERROR((persistent_data->vecx).reset())); }
    if(persistent_data->vecAx) { check_cuda_error(DPCT_CHECK_ERROR((persistent_data->vecAx).reset())); }
    if(persistent_data->vecp) { check_cuda_error(DPCT_CHECK_ERROR((persistent_data->vecp).reset())); }

    check_cuda_error(DPCT_CHECK_ERROR(dpct::dpct_free(persistent_data->d_col, dpct::get_in_order_queue())));
    check_cuda_error(DPCT_CHECK_ERROR(dpct::dpct_free(persistent_data->d_row, dpct::get_in_order_queue())));
    check_cuda_error(DPCT_CHECK_ERROR(dpct::dpct_free(persistent_data->d_val, dpct::get_in_order_queue())));
    check_cuda_error(DPCT_CHECK_ERROR(dpct::dpct_free(persistent_data->d_x, dpct::get_in_order_queue())));
    check_cuda_error(DPCT_CHECK_ERROR(dpct::dpct_free(persistent_data->d_r, dpct::get_in_order_queue())));
    check_cuda_error(DPCT_CHECK_ERROR(dpct::dpct_free(persistent_data->d_p, dpct::get_in_order_queue())));
    check_cuda_error(DPCT_CHECK_ERROR(dpct::dpct_free(persistent_data->d_Ax, dpct::get_in_order_queue())));
    check_cuda_error(DPCT_CHECK_ERROR(dpct::dpct_free(persistent_data->d_y, dpct::get_in_order_queue())));

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

    dpct::get_in_order_queue().memcpy(persistent_data->d_r, rhs, N * sizeof(float)).wait(); // B

    alpha = 1.0f;
    alpham1 = -1.0f;
    beta = 0.0f;
    r0 = 0.0f;

    check_cuda_error(
        DPCT_CHECK_ERROR(dpct::sparse::spmv((persistent_data->cusparseHandle)->get_queue(), oneapi::mkl::transpose::nontrans, &alpha, persistent_data->matA,
                                            persistent_data->vecx, &beta, persistent_data->vecAx, dpct::library_data_t::real_float)));

    oneapi::mkl::blas::column_major::axpy((persistent_data->cublasHandle)->get_queue(), N, alpham1, persistent_data->d_Ax, 1, persistent_data->d_r, 1);
    /*
    DPCT1034:3: Migrated API does not return an error code. 0 is returned in the lambda. You may need to rewrite this code.
    */
    [&]() {
    dpct::blas::wrapper_float_out res_wrapper_ct6((persistent_data->cublasHandle)->get_queue(), &r1);
    oneapi::mkl::blas::column_major::dot((persistent_data->cublasHandle)->get_queue(), N, persistent_data->d_r, 1, persistent_data->d_r, 1,
                                         res_wrapper_ct6.get_ptr());
    return 0;
    }();

    k = 1;

    while(r1 >= tol && k <= max_its) {

        if(k > 1) {
            b = r1 / r0;
            oneapi::mkl::blas::column_major::scal((persistent_data->cublasHandle)->get_queue(), N, b, persistent_data->d_p, 1);
            oneapi::mkl::blas::column_major::axpy((persistent_data->cublasHandle)->get_queue(), N, alpha, persistent_data->d_r, 1, persistent_data->d_p, 1);

        } else {
            oneapi::mkl::blas::column_major::copy((persistent_data->cublasHandle)->get_queue(), N, persistent_data->d_r, 1, persistent_data->d_p, 1);
        }

        check_cuda_error(
            DPCT_CHECK_ERROR(dpct::sparse::spmv((persistent_data->cusparseHandle)->get_queue(), oneapi::mkl::transpose::nontrans, &alpha, persistent_data->matA,
                                                persistent_data->vecp, &beta, persistent_data->vecAx, dpct::library_data_t::real_float)));

        /*
        DPCT1034:4: Migrated API does not return an error code. 0 is returned in the lambda. You may need to rewrite this code.
        */
        [&]() {
        dpct::blas::wrapper_float_out res_wrapper_ct6((persistent_data->cublasHandle)->get_queue(), &dot);
        oneapi::mkl::blas::column_major::dot((persistent_data->cublasHandle)->get_queue(), N, persistent_data->d_p, 1, persistent_data->d_Ax, 1,
                                             res_wrapper_ct6.get_ptr());
        return 0;
        }();

        a = r1 / dot;

        oneapi::mkl::blas::column_major::axpy((persistent_data->cublasHandle)->get_queue(), N, a, persistent_data->d_p, 1, persistent_data->d_x, 1);
        na = -a;
        oneapi::mkl::blas::column_major::axpy((persistent_data->cublasHandle)->get_queue(), N, na, persistent_data->d_Ax, 1, persistent_data->d_r, 1);

        r0 = r1;
        /*
        DPCT1034:5: Migrated API does not return an error code. 0 is returned in the lambda. You may need to rewrite this code.
        */
        [&]() {
        dpct::blas::wrapper_float_out res_wrapper_ct6((persistent_data->cublasHandle)->get_queue(), &r1);
        oneapi::mkl::blas::column_major::dot((persistent_data->cublasHandle)->get_queue(), N, persistent_data->d_r, 1, persistent_data->d_r, 1,
                                             res_wrapper_ct6.get_ptr());
        return 0;
        }();

        dpct::get_current_device().queues_wait_and_throw();
        k++;
    }

    dpct::get_in_order_queue().memcpy(rhs, persistent_data->d_x, N * sizeof(float)).wait();

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

    check_cuda_error((dpct::err0)DPCT_CHECK_ERROR(persistent_data->cublasHandle = new dpct::blas::descriptor()));
    check_cuda_error((dpct::err0)DPCT_CHECK_ERROR(persistent_data->cusparseHandle = new dpct::sparse::descriptor()));

    grid_to_csr(the_grid, &val, &I, &J, false);

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **active_cells = the_grid->active_cells;

    int nz = arrlen(val);
    int N = num_active_cells;

    persistent_data->N = N;
    persistent_data->nz = nz;

    check_cuda_error(DPCT_CHECK_ERROR(persistent_data->d_col = sycl::malloc_device<int>(nz, dpct::get_in_order_queue())));
    check_cuda_error(DPCT_CHECK_ERROR(persistent_data->d_row = sycl::malloc_device<int>((N + 1), dpct::get_in_order_queue())));
    check_cuda_error(DPCT_CHECK_ERROR(persistent_data->d_val = sycl::malloc_device<float>(nz, dpct::get_in_order_queue())));
    check_cuda_error(DPCT_CHECK_ERROR(persistent_data->d_x = sycl::malloc_device<float>(N, dpct::get_in_order_queue())));
    check_cuda_error(DPCT_CHECK_ERROR(persistent_data->d_r = sycl::malloc_device<float>(N, dpct::get_in_order_queue())));
    check_cuda_error(DPCT_CHECK_ERROR(persistent_data->d_rw = sycl::malloc_device<float>(N, dpct::get_in_order_queue())));
    check_cuda_error(DPCT_CHECK_ERROR(persistent_data->d_p = sycl::malloc_device<float>(N, dpct::get_in_order_queue())));
    check_cuda_error(DPCT_CHECK_ERROR(persistent_data->d_v = sycl::malloc_device<float>(N, dpct::get_in_order_queue())));
    check_cuda_error(DPCT_CHECK_ERROR(persistent_data->d_t = sycl::malloc_device<float>(N, dpct::get_in_order_queue())));
    check_cuda_error(DPCT_CHECK_ERROR(persistent_data->d_Ax = sycl::malloc_device<float>(N, dpct::get_in_order_queue())));

    dpct::get_in_order_queue().memcpy(persistent_data->d_col, J, nz * sizeof(int));            // JA
    dpct::get_in_order_queue().memcpy(persistent_data->d_row, I, (N + 1) * sizeof(int));       // IA
    dpct::get_in_order_queue().memcpy(persistent_data->d_val, val, nz * sizeof(float)).wait(); // A

    check_cuda_error(DPCT_CHECK_ERROR(dpct::get_in_order_queue().memset(persistent_data->d_r, 0, N * sizeof(float)).wait()));
    check_cuda_error(DPCT_CHECK_ERROR(dpct::get_in_order_queue().memset(persistent_data->d_rw, 0, N * sizeof(float)).wait()));
    check_cuda_error(DPCT_CHECK_ERROR(dpct::get_in_order_queue().memset(persistent_data->d_p, 0, N * sizeof(float)).wait()));
    check_cuda_error(DPCT_CHECK_ERROR(dpct::get_in_order_queue().memset(persistent_data->d_t, 0, N * sizeof(float)).wait()));
    check_cuda_error(DPCT_CHECK_ERROR(dpct::get_in_order_queue().memset(persistent_data->d_v, 0, N * sizeof(float)).wait()));

    float *rhs = (float *)malloc(sizeof(float) * num_active_cells);

    OMP(parallel for)
    for(uint32_t i = 0; i < num_active_cells; i++) {
        rhs[i] = active_cells[i]->b;
    }

    check_cuda_error(DPCT_CHECK_ERROR(dpct::get_in_order_queue().memcpy(persistent_data->d_x, rhs, N * sizeof(float)).wait())); // Result

    float alpha = 1.0f;
    float beta = 0.0f;

    check_cuda_error(DPCT_CHECK_ERROR(persistent_data->matA = std::make_shared<dpct::sparse::sparse_matrix_desc>(
                                          N, N, nz, persistent_data->d_row, persistent_data->d_col, persistent_data->d_val, dpct::library_data_t::real_int32,
                                          dpct::library_data_t::real_int32, oneapi::mkl::index_base::zero, dpct::library_data_t::real_float,
                                          dpct::sparse::matrix_format::csr)));

    check_cuda_error(
        DPCT_CHECK_ERROR(persistent_data->vecx = std::make_shared<dpct::sparse::dense_vector_desc>(N, persistent_data->d_x, dpct::library_data_t::real_float)));
    check_cuda_error(
        DPCT_CHECK_ERROR(persistent_data->vecp = std::make_shared<dpct::sparse::dense_vector_desc>(N, persistent_data->d_p, dpct::library_data_t::real_float)));
    check_cuda_error(
        DPCT_CHECK_ERROR(persistent_data->vecv = std::make_shared<dpct::sparse::dense_vector_desc>(N, persistent_data->d_v, dpct::library_data_t::real_float)));
    check_cuda_error(
        DPCT_CHECK_ERROR(persistent_data->vecr = std::make_shared<dpct::sparse::dense_vector_desc>(N, persistent_data->d_r, dpct::library_data_t::real_float)));
    check_cuda_error(
        DPCT_CHECK_ERROR(persistent_data->vect = std::make_shared<dpct::sparse::dense_vector_desc>(N, persistent_data->d_t, dpct::library_data_t::real_float)));
    check_cuda_error(DPCT_CHECK_ERROR(persistent_data->vecAx =
                                          std::make_shared<dpct::sparse::dense_vector_desc>(N, persistent_data->d_Ax, dpct::library_data_t::real_float)));

    check_cuda_error(DPCT_CHECK_ERROR(persistent_data->bufferSizeMV = 0));

    check_cuda_error(DPCT_CHECK_ERROR(persistent_data->d_bufferMV = (void *)sycl::malloc_device(persistent_data->bufferSizeMV, dpct::get_in_order_queue())));

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

    check_cuda_error(DPCT_CHECK_ERROR(dpct::get_in_order_queue().memcpy(persistent_data->d_r, rhs, N * sizeof(float)).wait())); // B

    // compute initial residual r0=b-Ax0 (using initial guess in x)
    check_cublas_error(
        DPCT_CHECK_ERROR(dpct::sparse::spmv((persistent_data->cusparseHandle)->get_queue(), oneapi::mkl::transpose::nontrans, &one, persistent_data->matA,
                                            persistent_data->vecx, &zero, persistent_data->vecAx, dpct::library_data_t::real_float)));

    check_cublas_error(DPCT_CHECK_ERROR(
        oneapi::mkl::blas::column_major::axpy((persistent_data->cublasHandle)->get_queue(), N, mone, persistent_data->d_Ax, 1, persistent_data->d_r, 1)));

    // copy residual r into r^{\hat} and p
    check_cublas_error(DPCT_CHECK_ERROR(
        oneapi::mkl::blas::column_major::copy((persistent_data->cublasHandle)->get_queue(), N, persistent_data->d_r, 1, persistent_data->d_rw, 1)));
    check_cublas_error(DPCT_CHECK_ERROR(
        oneapi::mkl::blas::column_major::copy((persistent_data->cublasHandle)->get_queue(), N, persistent_data->d_r, 1, persistent_data->d_p, 1)));
    /*
    DPCT1034:6: Migrated API does not return an error code. 0 is returned in the lambda. You may need to rewrite this code.
    */
    check_cublas_error([&]() {
    dpct::blas::wrapper_float_out res_wrapper_ct4((persistent_data->cublasHandle)->get_queue(), &nrmr0);
    oneapi::mkl::blas::column_major::nrm2((persistent_data->cublasHandle)->get_queue(), N, persistent_data->d_r, 1, res_wrapper_ct4.get_ptr());
    return 0;
    }());

    int i;
    for(i = 0; i < max_its;) {
        rhop = rho;
        /*
        DPCT1034:7: Migrated API does not return an error code. 0 is returned in the lambda. You may need to rewrite this code.
        */
        check_cuda_error((dpct::err0)[&]() {
        dpct::blas::wrapper_float_out res_wrapper_ct6((persistent_data->cublasHandle)->get_queue(), &rho);
        oneapi::mkl::blas::column_major::dot((persistent_data->cublasHandle)->get_queue(), N, persistent_data->d_rw, 1, persistent_data->d_r, 1,
                                             res_wrapper_ct6.get_ptr());
        return 0;
        }());

        if(i > 0) {
            beta = (rho / rhop) * (alpha / omega);
            negomega = -omega;
            check_cublas_error(DPCT_CHECK_ERROR(oneapi::mkl::blas::column_major::axpy((persistent_data->cublasHandle)->get_queue(), N, negomega,
                                                                                      persistent_data->d_v, 1, persistent_data->d_p, 1)));
            check_cublas_error(
                DPCT_CHECK_ERROR(oneapi::mkl::blas::column_major::scal((persistent_data->cublasHandle)->get_queue(), N, beta, persistent_data->d_p, 1)));
            check_cublas_error(DPCT_CHECK_ERROR(
                oneapi::mkl::blas::column_major::axpy((persistent_data->cublasHandle)->get_queue(), N, one, persistent_data->d_r, 1, persistent_data->d_p, 1)));
        }

        // matrix-vector multiplication
        check_cublas_error(
            DPCT_CHECK_ERROR(dpct::sparse::spmv((persistent_data->cusparseHandle)->get_queue(), oneapi::mkl::transpose::nontrans, &one, persistent_data->matA,
                                                persistent_data->vecp, &zero, persistent_data->vecv, dpct::library_data_t::real_float)));

        /*
        DPCT1034:8: Migrated API does not return an error code. 0 is returned in the lambda. You may need to rewrite this code.
        */
        check_cublas_error([&]() {
        dpct::blas::wrapper_float_out res_wrapper_ct6((persistent_data->cublasHandle)->get_queue(), &temp);
        oneapi::mkl::blas::column_major::dot((persistent_data->cublasHandle)->get_queue(), N, persistent_data->d_rw, 1, persistent_data->d_v, 1,
                                             res_wrapper_ct6.get_ptr());
        return 0;
        }());
        alpha = rho / temp;
        negalpha = -(alpha);
        check_cublas_error(DPCT_CHECK_ERROR(oneapi::mkl::blas::column_major::axpy((persistent_data->cublasHandle)->get_queue(), N, negalpha,
                                                                                  persistent_data->d_v, 1, persistent_data->d_r, 1)));
        check_cublas_error(DPCT_CHECK_ERROR(
            oneapi::mkl::blas::column_major::axpy((persistent_data->cublasHandle)->get_queue(), N, alpha, persistent_data->d_p, 1, persistent_data->d_x, 1)));
        /*
        DPCT1034:9: Migrated API does not return an error code. 0 is returned in the lambda. You may need to rewrite this code.
        */
        check_cublas_error([&]() {
        dpct::blas::wrapper_float_out res_wrapper_ct4((persistent_data->cublasHandle)->get_queue(), &nrmr);
        oneapi::mkl::blas::column_major::nrm2((persistent_data->cublasHandle)->get_queue(), N, persistent_data->d_r, 1, res_wrapper_ct4.get_ptr());
        return 0;
        }());

        if(nrmr < tol * nrmr0) {
            break;
        }

        // matrix-vector multiplication
        check_cublas_error(
            DPCT_CHECK_ERROR(dpct::sparse::spmv((persistent_data->cusparseHandle)->get_queue(), oneapi::mkl::transpose::nontrans, &one, persistent_data->matA,
                                                persistent_data->vecr, &zero, persistent_data->vect, dpct::library_data_t::real_float)));

        /*
        DPCT1034:10: Migrated API does not return an error code. 0 is returned in the lambda. You may need to rewrite this code.
        */
        check_cublas_error([&]() {
        dpct::blas::wrapper_float_out res_wrapper_ct6((persistent_data->cublasHandle)->get_queue(), &temp);
        oneapi::mkl::blas::column_major::dot((persistent_data->cublasHandle)->get_queue(), N, persistent_data->d_t, 1, persistent_data->d_r, 1,
                                             res_wrapper_ct6.get_ptr());
        return 0;
        }());
        /*
        DPCT1034:11: Migrated API does not return an error code. 0 is returned in the lambda. You may need to rewrite this code.
        */
        check_cublas_error([&]() {
        dpct::blas::wrapper_float_out res_wrapper_ct6((persistent_data->cublasHandle)->get_queue(), &temp2);
        oneapi::mkl::blas::column_major::dot((persistent_data->cublasHandle)->get_queue(), N, persistent_data->d_t, 1, persistent_data->d_t, 1,
                                             res_wrapper_ct6.get_ptr());
        return 0;
        }());
        omega = temp / temp2;
        negomega = -(omega);
        check_cublas_error(DPCT_CHECK_ERROR(
            oneapi::mkl::blas::column_major::axpy((persistent_data->cublasHandle)->get_queue(), N, omega, persistent_data->d_r, 1, persistent_data->d_x, 1)));
        check_cublas_error(DPCT_CHECK_ERROR(oneapi::mkl::blas::column_major::axpy((persistent_data->cublasHandle)->get_queue(), N, negomega,
                                                                                  persistent_data->d_t, 1, persistent_data->d_r, 1)));

        /*
        DPCT1034:12: Migrated API does not return an error code. 0 is returned in the lambda. You may need to rewrite this code.
        */
        check_cublas_error([&]() {
        dpct::blas::wrapper_float_out res_wrapper_ct4((persistent_data->cublasHandle)->get_queue(), &nrmr);
        oneapi::mkl::blas::column_major::nrm2((persistent_data->cublasHandle)->get_queue(), N, persistent_data->d_r, 1, res_wrapper_ct4.get_ptr());
        return 0;
        }());

        if(nrmr < tol * nrmr0) {
            i++;
            break;
        }
        i++;
    }

    check_cuda_error(DPCT_CHECK_ERROR(dpct::get_in_order_queue().memcpy(rhs, persistent_data->d_x, N * sizeof(float)).wait()));

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

    check_cuda_error((dpct::err0)DPCT_CHECK_ERROR(delete(persistent_data->cusparseHandle)));
    check_cuda_error((dpct::err0)DPCT_CHECK_ERROR(delete(persistent_data->cublasHandle)));

    if(persistent_data->matA) { check_cuda_error(DPCT_CHECK_ERROR((persistent_data->matA).reset())); }
    if(persistent_data->vecx) { check_cuda_error(DPCT_CHECK_ERROR((persistent_data->vecx).reset())); }
    if(persistent_data->vecAx) { check_cuda_error(DPCT_CHECK_ERROR((persistent_data->vecAx).reset())); }
    if(persistent_data->vecp) { check_cuda_error(DPCT_CHECK_ERROR((persistent_data->vecp).reset())); }
    if(persistent_data->vecv) { check_cuda_error(DPCT_CHECK_ERROR((persistent_data->vecv).reset())); }
    if(persistent_data->vecr) { check_cuda_error(DPCT_CHECK_ERROR((persistent_data->vecr).reset())); }
    if(persistent_data->vect) { check_cuda_error(DPCT_CHECK_ERROR((persistent_data->vect).reset())); }

    check_cuda_error(DPCT_CHECK_ERROR(dpct::dpct_free(persistent_data->d_col, dpct::get_in_order_queue())));
    check_cuda_error(DPCT_CHECK_ERROR(dpct::dpct_free(persistent_data->d_row, dpct::get_in_order_queue())));
    check_cuda_error(DPCT_CHECK_ERROR(dpct::dpct_free(persistent_data->d_val, dpct::get_in_order_queue())));
    check_cuda_error(DPCT_CHECK_ERROR(dpct::dpct_free(persistent_data->d_x, dpct::get_in_order_queue())));
    check_cuda_error(DPCT_CHECK_ERROR(dpct::dpct_free(persistent_data->d_r, dpct::get_in_order_queue())));
    check_cuda_error(DPCT_CHECK_ERROR(dpct::dpct_free(persistent_data->d_rw, dpct::get_in_order_queue())));
    check_cuda_error(DPCT_CHECK_ERROR(dpct::dpct_free(persistent_data->d_p, dpct::get_in_order_queue())));
    check_cuda_error(DPCT_CHECK_ERROR(dpct::dpct_free(persistent_data->d_t, dpct::get_in_order_queue())));
    check_cuda_error(DPCT_CHECK_ERROR(dpct::dpct_free(persistent_data->d_v, dpct::get_in_order_queue())));
    check_cuda_error(DPCT_CHECK_ERROR(dpct::dpct_free(persistent_data->d_Ax, dpct::get_in_order_queue())));
    check_cuda_error(DPCT_CHECK_ERROR(dpct::dpct_free(persistent_data->d_y, dpct::get_in_order_queue())));

    free(persistent_data);
}
