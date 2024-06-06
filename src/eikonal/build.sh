if [ -n "$CUDA_FOUND" ]; then

    EIKONAL_SOURCE_FILES="eikonal_solver.c cuda_fim.cu cuda_fim_kernel.cu"
    EIKONAL_HEADER_FILES="eikonal_solver.h common_def.h cuda_fim.h cuda_fim_kernel.h"
    EIKONAL_EXTRA_LIB_PATH=$CUDA_LIBRARY_PATH
    EIKONAL_DYNAMIC_LIBS="c cudart"

    COMPILE_SHARED_LIB "eikonal_solver" "$EIKONAL_SOURCE_FILES" "$EIKONAL_HEADER_FILES" "" "$EIKONAL_DYNAMIC_LIBS" "$EIKONAL_EXTRA_LIB_PATH" "" "$CUDA_FOUND"

fi
