//
// Created by sachetto on 28/11/24.
//

#ifndef ACCEL_UTILS_H
#define ACCEL_UTILS_H

#include <stddef.h>
#include <stdint.h>
#include "../common_types/common_types.h"

#ifdef COMPILE_CUDA
#include <cublas_v2.h>
#include <cusparse_v2.h>
#elif defined(COMPILE_SYCL)
#include <sycl/sycl.hpp>
#include <dpct/dpct.hpp>
#include <dpct/sparse_utils.hpp>
#include <dpct/blas_utils.hpp>
#endif

enum copy_direction {
    HOST_TO_DEVICE,
    DEVICE_TO_HOST,
};

enum sparse_index_type {
    INDEX_INT32,
    INDEX_INT64,
};

enum sparse_index_base {
    INDEX_BASE_ZERO,
};

enum data_type {
    REAL_FLOAT,
    REAL_DOUBLE,
};

#ifdef __cplusplus
extern "C" {
#endif
void malloc_device(void **ptr, size_t n);
void free_device(void *ptr);
void memcpy_device(void *dest, const void *src, size_t n, enum copy_direction kind);
void create_sparse_handle(void **handle);
void create_blas_handle(void **handle);
void sparse_create_scr(void *mat, int64_t rows, int64_t cols, int64_t nnz,
                       void* row_ptr,
                       void* cols_ind_ptr,
                       void* vals_ptr,
                       enum sparse_index_type csr_row_offsets_type,
                       enum sparse_index_type csr_col_ind_type,
                       enum sparse_index_base idx_base,
                       enum data_type         value_type);

void create_dense_vector(void** descr, int64_t size, void *values, enum data_type valueType);

#ifdef COMPILE_CUDA
void sparse_spmv(cusparseHandle_t          handle,
                                cusparseOperation_t       opA,
                                const void*               alpha,
                                cusparseConstSpMatDescr_t matA,
                                cusparseConstDnVecDescr_t vecX,
                                const void*               beta,
                                cusparseDnVecDescr_t      vecY,
                                enum data_type              computeType,
                                void*                     externalBuffer);
#elif defined(COMPILE_SYCL)
void sparse_spmv(dpct::sparse::descriptor_ptr handle, oneapi::mkl::transpose opA,
                            const void *alpha, dpct::sparse::sparse_matrix_desc_t matA,
                            std::shared_ptr<dpct::sparse::dense_vector_desc> vecX, const void *beta,
                            std::shared_ptr<dpct::sparse::dense_vector_desc> vecY,
                            enum data_type computeType, void *externalBuffer);
#endif

void blas_dot(void *handle, int n, real *x, int incx, real *y, int incy, real *result);

#ifdef __cplusplus
}
#endif
#endif //ACCEL_UTILS_H
