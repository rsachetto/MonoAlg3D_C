SET(CUDA_NVCC_FLAGS "-Xcompiler -fopenmp" CACHE STRING "nvcc flags" FORCE)
#HW capability 1.3
SET (CUDA_VERBOSE_BUILD ON CACHE BOOL "nvcc verbose" FORCE)

set(GPU_FILE gpu_common.cu gpu_common.h)

SET(LIB_TYPE STATIC) #set the lib type
CUDA_ADD_LIBRARY(CudaLib ${LIB_TYPE} ${GPU_FILE})