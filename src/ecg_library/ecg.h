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

#ifdef COMPILE_CUDA
    cusparseHandle_t cusparseHandle;
    cublasHandle_t cublasHandle;
    int *d_col, *d_row, nz;
    real *d_distances;
    real *d_volumes;
    real *volumes;
    real *tmp_data;
    real *d_val;

    size_t bufferSize;
    void *buffer;

#if CUBLAS_VER_MAJOR <= 10
    cusparseMatDescr_t descr;
    real *local_sv;
#else
    cusparseSpMatDescr_t matA;
    cusparseDnVecDescr_t vec_vm;
    cusparseDnVecDescr_t vec_beta_im;
#endif

#endif
};

#define PSEUDO_BIDOMAIN_DATA ((struct pseudo_bidomain_persistent_data *)config->persistent_data)
