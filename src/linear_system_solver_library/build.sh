if [ -n "$USE_SYCL" ]; then
      SYCL_SOURCES="gpu_solvers_sycl.cpp"
      EXTRA_C_FLAGS="$EXTRA_C_FLAGS -fsycl -fsycl-targets=nvptx64-nvidia-cuda,x86_64,spir64"
      IS_SYCL="1"

      COMPILE_SHARED_LIB "default_linear_system_solver" "${SYCL_SOURCES}" "" "alg config_helpers utils tinyexpr" "$EXTRA_CUDA_LIBS" "" "$EXTRA_C_FLAGS" "" "$IS_SYCL"

elif [ -n "$CUDA_FOUND" ]; then
    EXTRA_CUDA_LIBS="cudart cublas cusparse"

    if [ -n "$AMGX_FOUND" ]; then
        EXTRA_CUDA_LIBS="$EXTRA_CUDA_LIBS $AMGX_LIBRARY"
    fi

    CHECK_CUSTOM_FILE

    COMPILE_SHARED_LIB "default_linear_system_solver" "linear_system_solver.c ${CUSTOM_FILE}" "gpu_solvers_cublas_12.c gpu_solvers_cublas_11.c gpu_solvers_cublas_10.c" "alg config_helpers utils tinyexpr" "$EXTRA_CUDA_LIBS" "$CUDA_LIBRARY_PATH $CUDA_MATH_LIBRARY_PATH $AMGX_LIBRARY_PATH" "$EXTRA_C_FLAGS"

fi


