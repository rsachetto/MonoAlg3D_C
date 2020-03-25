if [ -n "$CUDA_FOUND" ]; then
    EXTRA_CUDA_LIBS="cudart cublas cusparse"
fi

COMPILE_SHARED_LIB "default_linear_system_solver" "linear_system_solver.c" "" "alg config_helpers utils tinyexpr" "$EXTRA_CUDA_LIBS" "$CUDA_LIBRARY_PATH"
