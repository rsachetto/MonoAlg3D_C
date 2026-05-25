############## NOBLE 1962 ##############################
MODEL_FILE_CPU="noble_1962.c"
MODEL_FILE_GPU="noble_1962.cu"
MODEL_FILE_SYCL=""
COMMON_HEADERS="noble_1962.h"

COMPILE_MODEL_LIB "noble_1962" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$MODEL_FILE_SYCL" "$COMMON_HEADERS"
#########################################################