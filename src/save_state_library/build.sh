if [ -n "$CUDA_FOUND" ]; then
    EXTRA_CUDA_LIBS="cudart"
fi

COMPILE_SHARED_LIB "default_save_state" "save_state.c" "" "alg config_helpers utils" "$EXTRA_CUDA_LIBS" "$CUDA_LIBRARY_PATH"
