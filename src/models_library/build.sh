COMPILE_MODEL_LIB () {
    local LIB_NAME=$1
    local MODEL_FILE_CPU=$2
    local MODEL_FILE_GPU=$3
    local COMMON_HEADERS="../model_common.h ../default_solvers.c $4"
    local EXTRA_C_FLAGS=$5
    local EXTRA_STATIC_LIBS=$6

    local MODEL_SOURCES="$MODEL_FILE_CPU"    

	local MODELS_STATIC_DEPS="config_helpers utils tinyexpr $EXTRA_STATIC_LIBS"
 
    if [ -n "$CUDA_FOUND" ]; then
        MODELS_EXTRA_LIB_PATH=$CUDA_LIBRARY_PATH
        MODELS_DYNAMIC_LIBS="c cudart"

        MODEL_SOURCES="$MODEL_SOURCES $MODEL_FILE_GPU"
        COMMON_HEADERS="$COMMON_HEADERS ../default_solvers.cu"
    fi
    
    COMPILE_SHARED_LIB "$LIB_NAME" "$MODEL_SOURCES" "$COMMON_HEADERS" "$MODELS_STATIC_DEPS" "$MODELS_DYNAMIC_LIBS" "$MODELS_EXTRA_LIB_PATH" "$EXTRA_C_FLAGS" "$CUDA_FOUND"
}

for d in */ ; do
	if [ -f "${d}/build.sh" ]; then
		ADD_SUBDIRECTORY "$d"
	fi
done
