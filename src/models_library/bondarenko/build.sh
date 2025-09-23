################ BONDARENKO ##############################
MODEL_FILE_CPU="bondarenko_2004.c"
MODEL_FILE_GPU="bondarenko_2004_GPU.cu"
MODEL_FILE_SYCL=""
COMMON_HEADERS="bondarenko_2004.h"

COMPILE_MODEL_LIB "bondarenko_2004" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$MODEL_FILE_SYCL" "$COMMON_HEADERS"
###########################################################