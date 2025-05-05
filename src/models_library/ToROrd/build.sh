############### ToRORd ##############################
MODEL_FILE_CPU="ToRORd_fkatp_endo.c"
MODEL_FILE_GPU="ToRORd_fkatp_endo.cu"
MODEL_FILE_SYCL=""
COMMON_HEADERS="ToRORd_fkatp_endo.h"
#
COMPILE_MODEL_LIB "ToRORd_fkatp_endo" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$MODEL_FILE_SYCL" "$COMMON_HEADERS"

############## ToRORd fkatp Mixed ENDO_MID_EPI ##############################
MODEL_FILE_CPU="ToRORd_fkatp_mixed_endo_mid_epi.c"
MODEL_FILE_GPU="ToRORd_fkatp_mixed_endo_mid_epi.cu"
MODEL_FILE_SYCL=""
COMMON_HEADERS="ToRORd_fkatp_mixed_endo_mid_epi.h"

COMPILE_MODEL_LIB "ToRORd_fkatp_mixed_endo_mid_epi" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$MODEL_FILE_SYCL" "$COMMON_HEADERS"

############## ToRORd dynCl Mixed ENDO_MID_EPI ##############################
MODEL_FILE_CPU="ToRORd_dynCl_mixed_endo_mid_epi.c"
MODEL_FILE_GPU="ToRORd_dynCl_mixed_endo_mid_epi.cu"
MODEL_FILE_SYCL=""
COMMON_HEADERS="ToRORd_dynCl_mixed_endo_mid_epi.h"

COMPILE_MODEL_LIB "ToRORd_dynCl_mixed_endo_mid_epi" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$MODEL_FILE_SYCL" "$COMMON_HEADERS"

############## ToRORd Land Mixed ENDO_MID_EPI ##############################
MODEL_FILE_CPU="ToRORd_Land_mixed_endo_mid_epi.c"
MODEL_FILE_GPU="ToRORd_Land_mixed_endo_mid_epi.cu"
MODEL_FILE_SYCL=""
COMMON_HEADERS="ToRORd_Land_mixed_endo_mid_epi.h"

COMPILE_MODEL_LIB "ToRORd_Land_mixed_endo_mid_epi" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$MODEL_FILE_SYCL" "$COMMON_HEADERS"

############## ToRORd fkatp GKsGKrtjca adjusted Mixed ENDO_MID_EPI ##############################
MODEL_FILE_CPU="ToRORd_fkatp_mixed_endo_mid_epi_GKsGKrtjca_adjustments.c"
MODEL_FILE_GPU="ToRORd_fkatp_mixed_endo_mid_epi_GKsGKrtjca_adjustments.cu"
MODEL_FILE_SYCL=""
COMMON_HEADERS="ToRORd_fkatp_mixed_endo_mid_epi_GKsGKrtjca_adjustments.h"

COMPILE_MODEL_LIB "ToRORd_fkatp_mixed_endo_mid_epi_GKsGKrtjca_adjustments" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$MODEL_FILE_SYCL" "$COMMON_HEADERS"

