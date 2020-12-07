if [ -n "$CUDA_FOUND" ]; then
    EXTRA_CUDA_LIBS="cudart"
fi

CHECK_CUSTOM_FILE "custom_save_state_functions.c"

COMPILE_SHARED_LIB "default_save_state" "save_state.c ${CUSTOM_FILE}" "" "alg config_helpers utils sds tinyexpr" "$EXTRA_CUDA_LIBS" "$CUDA_LIBRARY_PATH"
