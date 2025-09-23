GPU_UTILS_SOURCE_FILES="accel_utils.cpp"
GPU_UTILS_HEADER_FILES="accel_utils.h"
EXTRA_C_FLAGS=""

    if [ -n "$USE_SYCL" ]; then
      EXTRA_C_FLAGS="-fsycl -fsycl-targets=nvptx64-nvidia-cuda,x86_64,spir64"
      IS_SYCL="1"
    elif [ -n "$CUDA_FOUND" ]; then
        GPU_UTILS_EXTRA_LIB_PATH=$CUDA_LIBRARY_PATH
        GPU_UTILS_DYNAMIC_LIBS="c cudart"
        GPU_UTILS_SOURCE_FILES="$GPU_UTILS_SOURCE_FILES gpu_utils.c gpu_utils.cu"
    fi

COMPILE_SHARED_LIB "gpu_utils" "$GPU_UTILS_SOURCE_FILES" "$GPU_UTILS_HEADER_FILES" "" "$GPU_UTILS_DYNAMIC_LIBS" "$GPU_UTILS_EXTRA_LIB_PATH" "$EXTRA_C_FLAGS" "$CUDA_FOUND" "$IS_SYCL"

