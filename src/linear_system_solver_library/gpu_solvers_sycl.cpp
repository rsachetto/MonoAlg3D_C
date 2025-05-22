#include <oneapi/mkl.hpp>
#include <dpct/dpct.hpp>
#include "../config/linear_system_solver_config.h"
#include "../config_helpers/config_helpers.h"
#include "../libraries_common/common_data_structures.h"
#include "../logger/logger.h"
#include "sycl/usm.hpp"

using namespace sycl;
using namespace oneapi;

static bool use_preconditioner = false;
static int max_its = 200;
static real_cpu tol = 1e-16;

struct gpu_persistent_data {
    int *d_row, *d_col;
    float *d_val, *d_x;
    float *d_r, *d_p, *d_Ax, *d_rw, *d_v, *d_t;

    int N, nz;
    float tol;
    int max_its;
};

static void gemv_csr_sycl(queue &q, int M, int N, float alpha, const float *values, const int *col_idx, const int *row_ptr, const float *x, float beta,
                          float *y) {

    q.submit([&](handler &h) {
         h.parallel_for(range<1>(M), [=](id<1> i) {
             float sum = 0.0f;
             int row_start = row_ptr[i];
             int row_end = row_ptr[i + 1];
             for(int j = row_start; j < row_end; ++j) {
                 sum += values[j] * x[col_idx[j]];
             }
             y[i] = alpha * sum + beta * y[i];
         });
     }).wait();
}

inline void sycl_axpy(queue &q, int n, float alpha, float *x, int incx, float *y, int incy) {
    q.submit([&](handler &h) {
         h.parallel_for(range<1>(n), [=](id<1> i) {
             int idx_x = i[0] * incx;
             int idx_y = i[0] * incy;
             y[idx_y] += alpha * x[idx_x];
         });
     }).wait();
}

inline void sycl_dot(queue &q, int n, float *x, int incx, float *y, int incy, float *result) {
    float *partial = malloc_shared<float>(n, q);
    q.submit([&](handler &h) {
         h.parallel_for(range<1>(n), [=](id<1> i) {
             int ix = i[0] * incx;
             int iy = i[0] * incy;
             partial[i[0]] = x[ix] * y[iy];
         });
     }).wait();
    float sum = 0.0f;
    for(int i = 0; i < n; i++)
        sum += partial[i];
    *result = sum;
    free(partial, q);
}

inline void sycl_scal(queue &q, int n, float alpha, float *x, int incx) {
    q.submit([&](handler &h) {
         h.parallel_for(range<1>(n), [=](id<1> i) {
             int idx = i[0] * incx;
             x[idx] *= alpha;
         });
     }).wait();
}

inline void sycl_copy(queue &q, int n, float *x, int incx, float *y, int incy) {
    q.submit([&](handler &h) {
         h.parallel_for(range<1>(n), [=](id<1> i) {
             int ix = i[0] * incx;
             int iy = i[0] * incy;
             y[iy] = x[ix];
         });
     }).wait();
}

inline void sycl_nrm2(queue &q, int n, float *x, int incx, float *res) {
    float *partial = malloc_shared<float>(n, q);
    q.submit([&](handler &h) {
         h.parallel_for(range<1>(n), [=](id<1> i) {
             int idx = i[0] * incx;
             partial[i[0]] = x[idx] * x[idx];
         });
     }).wait();
    float sum = 0.0f;
    for(int i = 0; i < n; i++)
        sum += partial[i];
    *res = std::sqrt(sum);
    free(partial, q);
}

extern "C" INIT_LINEAR_SYSTEM(init_sycl_conjugate_gradient) {

    struct gpu_persistent_data *persistent_data = CALLOC_ONE_TYPE(struct gpu_persistent_data);

    dpct::device_ext &dev = dpct::get_current_device();
    sycl::queue &q = dev.default_queue();

    int_array I = NULL, J = NULL;
    f32_array val = NULL;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, tol, config, "tolerance");
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, max_its, config, "max_iterations");
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(use_preconditioner, config, "use_preconditioner");

    uint32_t num_active_cells;
    struct cell_node **active_cells = NULL;

    if(is_purkinje) {
        grid_to_csr(the_grid, &val, &I, &J, true);
        num_active_cells = the_grid->purkinje->num_active_purkinje_cells;
        active_cells = the_grid->purkinje->purkinje_cells;
    } else {
        grid_to_csr(the_grid, &val, &I, &J, false);
        num_active_cells = the_grid->num_active_cells;
        active_cells = the_grid->active_cells;
    }

    int N = 0, nz = 0;

    nz = arrlen(val);
    N = num_active_cells;

    persistent_data->N = N;
    persistent_data->nz = nz;
    persistent_data->tol = tol;
    persistent_data->max_its = max_its;

    persistent_data->d_row = (int *)sycl::malloc_device((N + 1) * sizeof(int), q);
    persistent_data->d_col = (int *)sycl::malloc_device(nz * sizeof(int), q);
    persistent_data->d_val = (float *)sycl::malloc_device(nz * sizeof(float), q);
    persistent_data->d_x = (float *)sycl::malloc_device(N * sizeof(float), q);
    persistent_data->d_r = (float *)sycl::malloc_device(N * sizeof(float), q);
    persistent_data->d_p = (float *)sycl::malloc_device(N * sizeof(float), q);
    persistent_data->d_Ax = (float *)sycl::malloc_device(N * sizeof(float), q);

    float *rhs = (float *)malloc(sizeof(float) * num_active_cells);

    OMP(parallel for)
    for(uint32_t i = 0; i < num_active_cells; i++) {
        rhs[i] = active_cells[i]->b;
    }

    q.memcpy(persistent_data->d_row, I, (N + 1) * sizeof(int)).wait();
    q.memcpy(persistent_data->d_col, J, nz * sizeof(int)).wait();
    q.memcpy(persistent_data->d_val, val, nz * sizeof(float)).wait();
    q.memcpy(persistent_data->d_x, rhs, N * sizeof(float)).wait();
    config->persistent_data = persistent_data;
}

extern "C" SOLVE_LINEAR_SYSTEM(sycl_conjugate_gradient) {

    struct gpu_persistent_data *persistent_data = (struct gpu_persistent_data *)config->persistent_data;

    dpct::device_ext &dev = dpct::get_current_device();
    sycl::queue &q = dev.default_queue();

    int N = persistent_data->N;

    float *rhs; // Vector B
    rhs = (float *)malloc(sizeof(float) * num_active_cells);

    OMP(parallel for)
    for(uint32_t i = 0; i < num_active_cells; i++) {
        rhs[i] = active_cells[i]->b;
    }

    q.memcpy(persistent_data->d_r, rhs, N * sizeof(float)).wait();

    float r0 = 0.0f, r1 = 0.0f, dot, a, b, na;

    gemv_csr_sycl(q, N, N, 1.0, persistent_data->d_val, persistent_data->d_col, persistent_data->d_row, persistent_data->d_x, 0.0, persistent_data->d_Ax);
    sycl_axpy(q, N, -1.0, persistent_data->d_Ax, 1, persistent_data->d_r, 1);
    sycl_dot(q, N, persistent_data->d_r, 1, persistent_data->d_r, 1, &r1);

    int k = 1;
    while(r1 >= persistent_data->tol && k <= persistent_data->max_its) {
        if(k > 1) {
            b = r1 / r0;
            sycl_scal(q, N, b, persistent_data->d_p, 1);
            sycl_axpy(q, N, 1.0f, persistent_data->d_r, 1, persistent_data->d_p, 1);
        } else {
            sycl_copy(q, N, persistent_data->d_r, 1, persistent_data->d_p, 1);
        }

        gemv_csr_sycl(q, N, N, 1.0, persistent_data->d_val, persistent_data->d_col, persistent_data->d_row, persistent_data->d_p, 0.0, persistent_data->d_Ax);
        sycl_dot(q, N, persistent_data->d_p, 1, persistent_data->d_Ax, 1, &dot);

        a = r1 / dot;
        sycl_axpy(q, N, a, persistent_data->d_p, 1, persistent_data->d_x, 1);
        na = -a;
        sycl_axpy(q, N, na, persistent_data->d_Ax, 1, persistent_data->d_r, 1);

        r0 = r1;
        sycl_dot(q, N, persistent_data->d_r, 1, persistent_data->d_r, 1, &r1);
        k++;
    }

    q.memcpy(rhs, persistent_data->d_x, N * sizeof(float)).wait();

    *number_of_iterations = k - 1;
    *error = r1;

    OMP(parallel for)
    for(uint32_t i = 0; i < num_active_cells; i++) {
        active_cells[i]->v = rhs[i];
    }

    free(rhs);
}

extern "C" END_LINEAR_SYSTEM(end_sycl_conjugate_gradient) {

    struct gpu_persistent_data *persistent_data = (struct gpu_persistent_data *)config->persistent_data;
    dpct::device_ext &dev = dpct::get_current_device();
    sycl::queue &q = dev.default_queue();

    if(!persistent_data)
        return;

    dpct::dpct_free(persistent_data->d_row, q);
    dpct::dpct_free(persistent_data->d_col, q);
    dpct::dpct_free(persistent_data->d_val, q);
    dpct::dpct_free(persistent_data->d_x, q);
    dpct::dpct_free(persistent_data->d_r, q);
    dpct::dpct_free(persistent_data->d_p, q);
    dpct::dpct_free(persistent_data->d_Ax, q);
}
