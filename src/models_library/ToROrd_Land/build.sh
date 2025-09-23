############## ToRORd Land Mixed ENDO_MID_EPI ##############################
MODEL_FILE_CPU="ToRORd_Land_mixed_endo_mid_epi.c"
MODEL_FILE_GPU="ToRORd_Land_mixed_endo_mid_epi.cu"
MODEL_FILE_SYCL=""
COMMON_HEADERS="ToRORd_Land_mixed_endo_mid_epi.h"

COMPILE_MODEL_LIB "ToRORd_Land_mixed_endo_mid_epi" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$MODEL_FILE_SYCL" "$COMMON_HEADERS"
