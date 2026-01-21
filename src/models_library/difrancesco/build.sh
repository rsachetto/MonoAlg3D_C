############## DI & FRANCESCO 1985 ##############################
MODEL_FILE_CPU="difrancesco.c"
MODEL_FILE_GPU="difrancesco.cu"
MODEL_FILE_SYCL=""
COMMON_HEADERS="difrancesco.h"

COMPILE_MODEL_LIB "difrancesco" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$MODEL_FILE_SYCL" "$COMMON_HEADERS"
#########################################################