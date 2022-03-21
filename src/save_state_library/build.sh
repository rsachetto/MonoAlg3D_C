SAVE_STATE_STATIC_DEPS="alg config_helpers utils sds tinyexpr"

if [ -n "$CUDA_FOUND" ]; then
	SAVE_STATE_STATIC_DEPS="$SAVE_STATE_STATIC_DEPS"
    EXTRA_CUDA_LIBS="cudart"
fi

CHECK_CUSTOM_FILE

COMPILE_SHARED_LIB "default_save_state" "save_state.c ${CUSTOM_FILE}" "" "${SAVE_STATE_STATIC_DEPS}" "$EXTRA_CUDA_LIBS" "$CUDA_LIBRARY_PATH"
