#ifndef MONOALG3D_GPU_UTILS_H
#define MONOALG3D_GPU_UTILS_H


#include "cuda_runtime.h"

#define check_cuda_errors(ans) do{cuda_assert((ans), __FILE__, __LINE__);} while(0)
void cuda_assert(cudaError_t code, const char *file, int line);

#endif //MONOALG3D_GPU_UTILS_H