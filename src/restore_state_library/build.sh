RESTORE_STATIC_DEPS="alg config_helpers utils sds tinyexpr"

if [ -n "$CUDA_FOUND" ]; then
    RESTORE_STATIC_DEPS="$RESTORE_STATIC_DEPS"
    EXTRA_CUDA_LIBS="c cudart "
fi

CHECK_CUSTOM_FILE

COMPILE_SHARED_LIB "default_restore_state" "restore_state.c ${CUSTOM_FILE}" "" "${RESTORE_STATIC_DEPS}" "$EXTRA_CUDA_LIBS" "$CUDA_LIBRARY_PATH"
