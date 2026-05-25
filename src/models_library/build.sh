COMPILE_MODEL_LIB () {
    local LIB_NAME=$1
    local MODEL_FILE_CPU=$2
    local MODEL_FILE_GPU=$3
    local MODEL_FILE_SYCL=$4
    local COMMON_HEADERS="../model_common.h ../default_solvers.c $5"
    local EXTRA_C_FLAGS=$6
    local EXTRA_STATIC_LIBS=$7

    local MODEL_SOURCES="$MODEL_FILE_CPU"    

	local MODELS_STATIC_DEPS="config_helpers utils tinyexpr $EXTRA_STATIC_LIBS"

    if [ -n "$USE_SYCL" ] && [ -n "$MODEL_FILE_SYCL" ]; then
      MODEL_SOURCES="$MODEL_FILE_SYCL"
      EXTRA_C_FLAGS="$EXTRA_C_FLAGS -fsycl -fsycl-targets=nvptx64-nvidia-cuda,x86_64,spir64"
      IS_SYCL="1"
    elif [ -n "$CUDA_FOUND" ]; then
        MODELS_EXTRA_LIB_PATH=$CUDA_LIBRARY_PATH
        MODELS_DYNAMIC_LIBS="c cudart"

        MODEL_SOURCES="$MODEL_SOURCES $MODEL_FILE_GPU"
        COMMON_HEADERS="$COMMON_HEADERS ../default_solvers.cu"
    fi
    
    COMPILE_SHARED_LIB "$LIB_NAME" "$MODEL_SOURCES" "$COMMON_HEADERS" "$MODELS_STATIC_DEPS" "$MODELS_DYNAMIC_LIBS" "$MODELS_EXTRA_LIB_PATH" "$EXTRA_C_FLAGS" "$CUDA_FOUND" "$IS_SYCL"
}

for d in */ ; do
	if [ -f "${d}/build.sh" ]; then
		ADD_SUBDIRECTORY "$d"
	fi
done
