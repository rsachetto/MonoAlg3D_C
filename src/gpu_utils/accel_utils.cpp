//
// Created by sachetto on 28/11/24.
//

#include "accel_utils.h"
#include "gpu_utils.h"

#ifdef COMPILE_CUDA
#include <cublas_v2.h>
#include <cusparse_v2.h>
#endif

extern "C" void malloc_device(void **ptr, size_t n) {

#ifdef COMPILE_CUDA
    check_cuda_error(cudaMalloc(ptr, n));
#elif defined(COMPILE_SYCL)
    DPCT_CHECK_ERROR(ptr = sycl::malloc_device(n, dpct::get_in_order_queue()));
#endif

}

extern "C" void free_device(void *ptr) {
#ifdef COMPILE_CUDA
    check_cuda_error(cudaFree(ptr));
#elif defined(COMPILE_SYCL)
    DPCT_CHECK_ERROR(dpct::dpct_free(persistent_data->d_col, dpct::get_in_order_queue()));
#endif
}

extern "C" void memcpy_device(void *dest, const void *src, size_t n, copy_direction kind) {

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

extern "C" void create_sparse_handle(void *handle) {
#ifdef COMPILE_CUDA
    check_cublas_error(cusparseCreate((cusparseHandle_t *)handle));
#elif defined(COMPILE_SYCL)
    DPCT_CHECK_ERROR(handle = new dpct::sparse::descriptor();
#endif
}

extern "C" void create_blas_handle(void *handle) {
#ifdef COMPILE_CUDA
    check_cublas_error(cublasCreate((cublasHandle_t *)handle));
#elif defined(COMPILE_SYCL)
    DPCT_CHECK_ERROR(handle = new dpct::blas::descriptor();
#endif
}

extern "C" void sparse_create_scr(void *mat, int64_t rows, int64_t cols, int64_t nnz,
                  void* csrRowOffsets,
                  void* csrColInd,
                  void* csrValues,
                  cusparseIndexType_t   csrRowOffsetsType,
                  cusparseIndexType_t   csrColIndType,
                  cusparseIndexBase_t   idxBase,
                  cudaDataType          valueType) {
#ifdef COMPILE_CUDA
    check_cuda_error(cusparseCreateCsr(&(PSEUDO_BIDOMAIN_DATA->matA), N, N, nz, PSEUDO_BIDOMAIN_DATA->d_row, PSEUDO_BIDOMAIN_DATA->d_col,
                                       PSEUDO_BIDOMAIN_DATA->d_val, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUBLAS_SIZE));
#elif defined(COMPILE_SYCL)
    DPCT_CHECK_ERROR(mat = new dpct::sparse::sparse_matrix_desc(
                                rows, cols, nnz, csrRowOffsets,csrColInd, csrValues, dpct::library_data_t::real_int32,
                                dpct::library_data_t::real_int32, oneapi::mkl::index_base::zero, CUBLAS_SIZE, dpct::sparse::matrix_format::csr));
#endif
}