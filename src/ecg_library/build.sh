if [ -n "$CUDA_FOUND" ]; then
    EXTRA_CUDA_LIBS="cudart cublas cusparse"
fi


CHECK_CUSTOM_FILE
COMPILE_SHARED_LIB "default_calc_ecg" "ecg.c ${CUSTOM_FILE}" "" "config_helpers vtk_utils utils alg sds tinyexpr" "m $EXTRA_CUDA_LIBS" "$CUDA_LIBRARY_PATH"
