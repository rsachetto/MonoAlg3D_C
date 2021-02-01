############### ToRORd (not working) ##############################
MODEL_FILE_CPU="ToRORd_fkatp_endo.c"
MODEL_FILE_GPU="ToRORd_fkatp_endo.cu"
COMMON_HEADERS="ToRORd_fkatp_endo.h"
#
COMPILE_MODEL_LIB "ToRORd_fkatp_endo" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$COMMON_HEADERS"

############## ToRORd (working) ##############################
MODEL_FILE_CPU="ToRORd_fkatp_endo_2019.c"
MODEL_FILE_GPU="ToRORd_fkatp_endo_2019.cu"
COMMON_HEADERS="ToRORd_fkatp_endo_2019.h"

COMPILE_MODEL_LIB "ToRORd_fkatp_endo_2019" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$COMMON_HEADERS"
