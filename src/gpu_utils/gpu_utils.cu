#include "gpu_utils.h"

__global__ void kernel_gpu_vec_div_vec(real *vec1, real *vec2, real *vec3, size_t n) {
    
    unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;

    if(i < n) {
        vec3[i] = vec1[i]/vec2[i];
    }
}

extern "C" void gpu_vec_div_vec(real *vec1, real *vec2, real *res, size_t n) {
    const int GRID  = (n + BLOCK_SIZE - 1)/BLOCK_SIZE;
    kernel_gpu_vec_div_vec<<<GRID, BLOCK_SIZE>>>(vec1, vec2, res, n);
    cudaDeviceSynchronize();
}


__global__ void kernel_update_monodomain(real *vm, float *b, real alpha, size_t n) {

    unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;

    if(i < n) {
        b[i] = vm[i] * alpha;
    }
}

extern "C" void gpu_update_monodomain(real *vm, float *b, real alpha, size_t n) {
    const int GRID  = (n + BLOCK_SIZE - 1)/BLOCK_SIZE;
    kernel_update_monodomain<<<GRID, BLOCK_SIZE>>>(vm, b, alpha, n);
    cudaDeviceSynchronize();
}
