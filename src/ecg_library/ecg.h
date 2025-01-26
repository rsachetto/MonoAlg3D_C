#include "../libraries_common/common_data_structures.h"
#include <stdio.h>

#ifdef COMPILE_CUDA
#include <cublas_v2.h>
#include <cusparse_v2.h>
#endif

struct pseudo_bidomain_persistent_data {
    real *distances;
    real *beta_im;
    struct point_3d *leads;
    FILE *output_file;
    real_cpu scale_factor;
    uint32_t n_leads;
    uint32_t diff_curr_rate;
    real_cpu diff_curr_max_time;

#if defined(COMPILE_CUDA) || defined(COMPILE_SYCL)
#ifdef COMPILE_CUDA
    cusparseHandle_t sparseHandle;
    cublasHandle_t blasHandle;
#else
    dpct::sparse::descriptor_ptr sparseHandle;
    dpct::blas::descriptor_ptr blasHandle;
#endif
    int *d_col, *d_row, nz;
    real *d_distances;
    real *d_volumes;
    real *volumes;
    real *tmp_data;
    real *d_val;
    real *beta_im_cpu;
    size_t bufferSize;
    void *buffer;
#if defined(COMPILE_CUDA) && CUBLAS_VER_MAJOR <= 10
    cusparseMatDescr_t descr;
    real *local_sv;
#elif defined(COMPILE_CUDA)
    cusparseSpMatDescr_t matA;
    cusparseDnVecDescr_t vec_vm;
    cusparseDnVecDescr_t vec_beta_im;
#elif defined(COMPILE_SYCL)
    dpct::sparse::sparse_matrix_desc_t matA;
    std::shared_ptr<dpct::sparse::dense_vector_desc> vec_vm;
    std::shared_ptr<dpct::sparse::dense_vector_desc> vec_beta_im;
#endif

#endif


};

#define PSEUDO_BIDOMAIN_DATA ((struct pseudo_bidomain_persistent_data *)config->persistent_data)
