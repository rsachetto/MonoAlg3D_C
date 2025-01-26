//
// Created by sachetto on 28/11/24.
//

#include <sycl/sycl.hpp>
#include <dpct/dpct.hpp>
#include "accel_utils.h"
#include "gpu_utils.h"
#include <dpct/blas_utils.hpp>

#include <dpct/sparse_utils.hpp>

#include <dpct/lib_common_utils.hpp>

#ifdef COMPILE_CUDA
#elif defined(COMPILE_SYCL)
#include <sycl/sycl.hpp>
#include <dpct/dpct.hpp>
#include <dpct/sparse_utils.hpp>
#include <dpct/blas_utils.hpp>
#endif

extern "C" void malloc_device(void **ptr, size_t n) {

#ifdef COMPILE_CUDA
    check_cuda_error(cudaMalloc(ptr, n));
#elif defined(COMPILE_SYCL)
    DPCT_CHECK_ERROR(*ptr = sycl::malloc_device(n, dpct::get_in_order_queue()));
#endif

}

extern "C" void free_device(void *ptr) {
#ifdef COMPILE_CUDA
    check_cuda_error(cudaFree(ptr));
#elif defined(COMPILE_SYCL)
    DPCT_CHECK_ERROR(dpct::dpct_free(ptr, dpct::get_in_order_queue()));
#endif
}

extern "C" void memcpy_device(void *dest, const void *src, size_t n, enum copy_direction kind) {

    if(kind == HOST_TO_DEVICE)  {
    #ifdef COMPILE_CUDA
        check_cuda_error(cudaMemcpy(dest, src, n, cudaMemcpyHostToDevice));
    #elif defined(COMPILE_SYCL)
        sycl::device dev_ct1;
        sycl::queue q_ct1(dev_ct1, sycl::property_list{sycl::property::queue::in_order()});
        q_ct1.memcpy(dest, src, n).wait();
    #endif
    } else if(kind == DEVICE_TO_HOST) {
    #ifdef COMPILE_CUDA
        check_cuda_error(cudaMemcpy(dest, src, n, cudaMemcpyDeviceToHost));
    #elif defined(COMPILE_SYCL)
        dpct::device_ext &dev_ct1 = dpct::get_current_device();
        sycl::queue &q_ct1 = dev_ct1.default_queue();
        q_ct1.memcpy(dest, src, n).wait();
    #endif
    }
}

extern "C" void create_sparse_handle(void **handle) {
#ifdef COMPILE_CUDA
    check_cublas_error(cusparseCreate((cusparseHandle_t *) *handle));
#elif defined(COMPILE_SYCL)
    DPCT_CHECK_ERROR(*handle = new dpct::sparse::descriptor());
#endif
}

extern "C" void create_blas_handle(void **handle) {
#ifdef COMPILE_CUDA
    check_cublas_error(cublasCreate((cublasHandle_t *) *handle));
#elif defined(COMPILE_SYCL)
    DPCT_CHECK_ERROR(*handle = new dpct::blas::descriptor());
#endif
}

#ifdef COMPILE_CUDA
static inline dpct::library_data_t get_index_type(enum sparse_index_type type) {
    switch(type) {
    case INDEX_INT32:
        return dpct::library_data_t::real_int32;
    case INDEX_INT64:
        return dpct::library_data_t::real_int64;
    }
    return dpct::library_data_t::real_int64;
}
#elif defined(COMPILE_SYCL)
static inline dpct::library_data_t get_index_type(enum sparse_index_type type) {
    switch(type) {
    case INDEX_INT32:
        return dpct::library_data_t::real_int32;
    case INDEX_INT64:
        return dpct::library_data_t::real_int64;
    }
    return dpct::library_data_t::real_int64;
}
#endif

#ifdef COMPILE_CUDA
static inline dpct::library_data_t get_data_type(enum data_type type) {
    switch(type) {
    case REAL_FLOAT:
        return 0;
    case REAL_DOUBLE:
        return 1;
    }

    return 1;
}
#elif defined(COMPILE_SYCL)
static inline dpct::library_data_t get_data_type(enum data_type type) {
    switch(type) {
    case REAL_FLOAT:
        return dpct::library_data_t::real_float;
    case REAL_DOUBLE:
        return dpct::library_data_t::real_double;
    }

    return dpct::library_data_t::real_double;
}
#endif

#ifdef COMPILE_CUDA
static inline oneapi::mkl::index_base get_index_base_type(enum sparse_index_base type) {
    switch(type) {
    case INDEX_BASE_ZERO:
        return oneapi::mkl::index_base::zero;
    }

    return oneapi::mkl::index_base::zero;
}
#elif defined(COMPILE_SYCL)
static inline oneapi::mkl::index_base  get_index_base_type(enum sparse_index_base type) {
    switch(type) {
    case INDEX_BASE_ZERO:
        return oneapi::mkl::index_base::zero;
    }
}
#endif

extern "C" void sparse_create_scr(void *mat, int64_t rows, int64_t cols, int64_t nnz,
                  void* row_ptr,
                  void* cols_ind_ptr,
                  void* vals_ptr,
                  enum sparse_index_type csr_row_offsets_type,
                  enum sparse_index_type   csr_col_ind_type,
                  enum sparse_index_base   idx_base,
                  enum data_type          value_type) {


#ifdef COMPILE_CUDA
    cusparseIndexType_t row_index_type = get_index_type(csr_row_offsets_type);
    cusparseIndexType_t col_index_type = get_index_type(csr_col_ind_type);
    cudaDataType l_value_type = get_data_type(value_type);
    cusparseIndexBase_t l_idx_base = get_index_base_type(idx_base);
    check_cuda_error(cusparseCreateCsr((cusparseSpMatDescr_t*) mat, rows, cols, nnz, row_ptr, cols_ind_ptr,
                                       vals_ptr, row_index_type, col_index_type, l_idx_base, l_value_type));
#elif defined(COMPILE_SYCL)
    dpct::library_data_t row_index_type = get_index_type(csr_row_offsets_type);
    dpct::library_data_t col_index_type = get_index_type(csr_col_ind_type);
    dpct::library_data_t l_value_type = get_data_type(value_type);
    oneapi::mkl::index_base l_idx_base = get_index_base_type(idx_base);
    DPCT_CHECK_ERROR(mat = new dpct::sparse::sparse_matrix_desc(rows, cols, nnz, row_ptr, cols_ind_ptr, vals_ptr, row_index_type,
                                col_index_type, l_idx_base, l_value_type, dpct::sparse::matrix_format::csr));
#endif
}