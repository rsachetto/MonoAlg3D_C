//
// Created by sachetto on 28/11/24.
//

#ifndef ACCEL_UTILS_H
#define ACCEL_UTILS_H

#include <stddef.h>
#include <stdint.h>

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
#ifdef __cplusplus
}
#endif
#endif //ACCEL_UTILS_H
