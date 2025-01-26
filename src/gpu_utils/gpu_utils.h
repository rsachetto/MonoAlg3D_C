#ifndef MONOALG3D_GPU_UTILS_H
#define MONOALG3D_GPU_UTILS_H

#ifdef COMPILE_CUDA
#include <cuda_runtime.h>
#endif
#include <stdlib.h>
#include <stdio.h>

#include "../common_types/common_types.h"

#define BLOCK_SIZE 32

#ifdef COMPILE_CUDA
#define check_cuda_error(ans) do{cuda_assert((ans), #ans, __FILE__, __LINE__, "cuda");} while(0)
#define check_cublas_error(ans) do{cuda_assert((ans), #ans, __FILE__, __LINE__, "cublas");} while(0)
#endif

#ifdef __cplusplus
extern "C" {
#endif
void gpu_vec_div_vec(real *vec1, real *vec2, real* res, size_t n);
#ifdef COMPILE_CUDA
void cuda_assert(int code, char const *const func, const char *const file, int const line, const char *api);
void uuid_print (cudaUUID_t a);
#endif
#ifdef __cplusplus
}
#endif
#endif //MONOALG3D_GPU_UTILS_H
