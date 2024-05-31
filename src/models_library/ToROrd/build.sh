############### ToRORd ##############################
MODEL_FILE_CPU="ToRORd_fkatp_endo.c"
MODEL_FILE_GPU="ToRORd_fkatp_endo.cu"
COMMON_HEADERS="ToRORd_fkatp_endo.h"
#
COMPILE_MODEL_LIB "ToRORd_fkatp_endo" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$COMMON_HEADERS"

############## ToRORd fkatp Mixed ENDO_MID_EPI ##############################
MODEL_FILE_CPU="ToRORd_fkatp_mixed_endo_mid_epi.c"
MODEL_FILE_GPU="ToRORd_fkatp_mixed_endo_mid_epi.cu"
COMMON_HEADERS="ToRORd_fkatp_mixed_endo_mid_epi.h"

COMPILE_MODEL_LIB "ToRORd_fkatp_mixed_endo_mid_epi" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$COMMON_HEADERS"

############## ToRORd dynCl Mixed ENDO_MID_EPI ##############################
MODEL_FILE_CPU="ToRORd_dynCl_mixed_endo_mid_epi.c"
MODEL_FILE_GPU="ToRORd_dynCl_mixed_endo_mid_epi.cu"
COMMON_HEADERS="ToRORd_dynCl_mixed_endo_mid_epi.h"

COMPILE_MODEL_LIB "ToRORd_dynCl_mixed_endo_mid_epi" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$COMMON_HEADERS"

############## ToRORd Land Mixed ENDO_MID_EPI ##############################
MODEL_FILE_CPU="ToRORd_Land_mixed_endo_mid_epi.c"
MODEL_FILE_GPU="ToRORd_Land_mixed_endo_mid_epi.cu"
COMMON_HEADERS="ToRORd_Land_mixed_endo_mid_epi.h"

COMPILE_MODEL_LIB "ToRORd_Land_mixed_endo_mid_epi" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$COMMON_HEADERS"

############## ToRORd Land Iks BaseToApex Mixed ENDO_MID_EPI ##############################
MODEL_FILE_CPU="ToRORd_Land_Iks_baseToApex_mixed_endo_mid_epi.c"
MODEL_FILE_GPU="ToRORd_Land_Iks_baseToApex_mixed_endo_mid_epi.cu"
COMMON_HEADERS="ToRORd_Land_Iks_baseToApex_mixed_endo_mid_epi.h"

COMPILE_MODEL_LIB "ToRORd_Land_Iks_baseToApex_mixed_endo_mid_epi" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$COMMON_HEADERS"

############## ToRORd fkatp Mixed ENDO_MID_EPI with GKs GKr Tjca adjustments for T wave  ##############################
MODEL_FILE_CPU="ToRORd_fkatp_mixed_endo_mid_epi_GKsGKrtjca_adjustments.c"
MODEL_FILE_GPU="ToRORd_fkatp_mixed_endo_mid_epi_GKsGKrtjca_adjustments.cu"
COMMON_HEADERS="ToRORd_fkatp_mixed_endo_mid_epi_GKsGKrtjca_adjustments.h"

COMPILE_MODEL_LIB "ToRORd_fkatp_mixed_endo_mid_epi_GKsGKrtjca_adjustments" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$COMMON_HEADERS"


