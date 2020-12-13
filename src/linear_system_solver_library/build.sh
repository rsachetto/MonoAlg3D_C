if [ -n "$CUDA_FOUND" ]; then
    EXTRA_CUDA_LIBS="cudart cublas cusparse"
fi

CHECK_CUSTOM_FILE

COMPILE_SHARED_LIB "default_linear_system_solver" "linear_system_solver.c ${CUSTOM_FILE}" "" "alg config_helpers utils tinyexpr gpu_utils" "$EXTRA_CUDA_LIBS" "$CUDA_LIBRARY_PATH"
