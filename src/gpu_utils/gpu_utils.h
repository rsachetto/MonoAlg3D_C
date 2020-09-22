#ifndef MONOALG3D_GPU_UTILS_H
#define MONOALG3D_GPU_UTILS_H

#include "cuda_runtime.h"

/*! \brief Macros/inlines to assist CLion to parse Cuda files (*.cu, *.cuh) */
#ifdef __JETBRAINS_IDE__
#define __CUDACC__ 1
#define __host__
#define __device__
#define __global__
#define __forceinline__
#define __shared__
inline void __syncthreads() {}
inline void __threadfence_block() {}
template<class T> inline T __clz(const T val) { return val; }
struct __cuda_fake_struct { int x; int y; int z; };
extern __cuda_fake_struct blockDim;
extern __cuda_fake_struct threadIdx;
extern __cuda_fake_struct blockIdx;
#endif

#include <stdlib.h>
#include <stdio.h>
#include "cuda_runtime.h"

#define BLOCK_SIZE 32

#define check_cuda_error(ans) do{cuda_assert((ans), #ans, __FILE__, __LINE__, "cuda");} while(0)
#define check_cublas_error(ans) do{cuda_assert((ans), #ans, __FILE__, __LINE__, "cublas");} while(0)

#ifdef __cplusplus
extern "C" {
#endif
void cuda_assert(int code, char const *const func, const char *const file, int const line, const char *api);
#ifdef __cplusplus
}
#endif
#endif //MONOALG3D_GPU_UTILS_H
