############## MITCHELL SHAEFFER 2003 ##############################
MODEL_FILE_CPU="mitchell_shaeffer_2003.c"
MODEL_FILE_GPU="mitchell_shaeffer_2003.cu"
MODEL_FILE_SYCL=""
COMMON_HEADERS="mitchell_shaeffer_2003.h"

COMPILE_MODEL_LIB "mitchell_shaeffer_2003" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$MODEL_FILE_SYCL" "$COMMON_HEADERS"
#########################################################