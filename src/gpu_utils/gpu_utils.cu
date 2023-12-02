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

// Adapted from: https://stackoverflow.com/questions/68823023/set-cuda-device-by-uuid
void uuid_print (cudaUUID_t a){
  printf("GPU UUID ");
  int r[5][2] = {{0,4}, {4,6}, {6,8}, {8,10}, {10,16}};
  for (int i = 0; i < 5; i++) {
      printf("-");
      for (int j = r[i][0]; j < r[i][1]; j++) {
          printf("%02x", (unsigned)(unsigned char)a.bytes[j]);
      }
  }
  printf("\n");
}