#include "gpu_utils.h"

__global__ void gpu_ecg_integral_kernel(const real *beta_im, const real* distances, const real *volumes, int n, real *result);
__global__ void kernel_gpu_vec_div_vec(real *vec1, real *vec2, real *vec3, size_t n);

    extern "C" void gpu_vec_div_vec(real *vec1, real *vec2, real *res, size_t n) {
    const int GRID  = (n + BLOCK_SIZE - 1)/BLOCK_SIZE;
    kernel_gpu_vec_div_vec<<<GRID, BLOCK_SIZE>>>(vec1, vec2, res, n);
    cudaDeviceSynchronize();
}

__global__ void kernel_gpu_vec_div_vec(real *vec1, real *vec2, real *vec3, size_t n) {
    
    unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;

    if(i < n) {
        vec3[i] = vec1[i]/vec2[i];
    }
}

