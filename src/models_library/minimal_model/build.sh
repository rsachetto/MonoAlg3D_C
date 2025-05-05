############## minimal_model ##############################
MODEL_FILE_CPU="minimal_model.c"
MODEL_FILE_GPU="minimal_model.cu"
MODEL_FILE_SYCL=""
COMMON_HEADERS="minimal_model.h"

COMPILE_MODEL_LIB "minimal_model" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$MODEL_FILE_SYCL" "$COMMON_HEADERS"
