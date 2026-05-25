############## CRN ##############################
MODEL_FILE_CPU="courtemanche_ramirez_nattel_1998.c"
MODEL_FILE_GPU="courtemanche_ramirez_nattel_1998.cu"
MODEL_FILE_SYCL=""
COMMON_HEADERS="courtemanche_ramirez_nattel_1998.h"

COMPILE_MODEL_LIB "courtemanche_ramirez_nattel_1998" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$MODEL_FILE_SYCL" "$COMMON_HEADERS"
##########################################################

############## CRN_RL ##############################
MODEL_FILE_CPU="courtemanche_ramirez_nattel_1998_RL.c"
MODEL_FILE_GPU="courtemanche_ramirez_nattel_1998_RL.cu"
MODEL_FILE_SYCL=""
COMMON_HEADERS="courtemanche_ramirez_nattel_1998_RL.h"

COMPILE_MODEL_LIB "courtemanche_ramirez_nattel_1998_RL" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$MODEL_FILE_SYCL" "$COMMON_HEADERS"
##########################################################

