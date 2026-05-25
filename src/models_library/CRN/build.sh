############### CRN2017 TASK1 and IKACh  #####################################################
MODEL_FILE_CPU="CRN2017_TASK1_IKACh.c"
MODEL_FILE_GPU="CRN2017_TASK1_IKACh.cu"
MODEL_FILE_SYCL=""
COMMON_HEADERS="CRN2017_TASK1_IKACh.h"


COMPILE_MODEL_LIB "CRN2017_TASK1_IKACh" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$MODEL_FILE_SYCL" "$COMMON_HEADERS"
##############################################################################################
