//
// Created by sachetto on 30/04/25.
//

#ifndef ACCEL_UTILS_H
#define ACCEL_UTILS_H

#include "../ode_solver/ode_solver.h"
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

enum copy_direction { HOST_TO_DEVICE, DEVICE_TO_HOST, DEVICE_TO_DEVICE };

#ifdef __cplusplus
extern "C" {
#endif
void *create_sparse_handle(int num_rows, int num_cols, int *d_row, int *d_col, float *d_val, void *v_q);
void destroy_sparse_handle(void *v_handle, void *v_q);
char *get_device_info(struct ode_solver *the_ode_solver, bool linear_solver_on_gpu);
void malloc_device(void **ptr, size_t n);
void free_device(void *ptr);
void memcpy_device(void *dest, const void *src, size_t n, enum copy_direction kind);
void memcpy2d_device(void *dest, size_t pitch_dest, const void *src, size_t pitch_src, size_t w, size_t h, enum copy_direction kind);
void *get_sycl_queue();
void mkl_dot(int n, float *x, int incx, float *y, int incy, float *result, void *v_q);
void mkl_scal(int n, float alpha, float *x, int incx, void *v_q);
void mkl_nrm2(int n, float *x, int incx, float *res, void *v_q);
void mkl_axpy(int n, float alpha, float *x, int incx, float *y, int incy, void *v_q);
void mkl_spmv(void *v_handle, float *in_vec, float *out_vec, void *v_q);
void mkl_copy(int n, float *x, int incx, float *y, int incy, void *v_q);
#ifdef __cplusplus
}
#endif

#endif // ACCEL_UTILS_H
