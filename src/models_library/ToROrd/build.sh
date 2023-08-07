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

############## Paci_ToRORd dynCl PhiCaL IKCa Mixed Apicobasal InfarctRemod ##############################
#MODEL_FILE_CPU="Paci_ToRORd_dynCl_PhiCaL_IKCa_mixed_apicobasal_infarctionRemod_RL.c"
#MODEL_FILE_GPU="Paci_ToRORd_dynCl_PhiCaL_IKCa_mixed_apicobasal_infarctionRemod_RL.cu"
#COMMON_HEADERS="Paci_ToRORd_dynCl_PhiCaL_IKCa_mixed_apicobasal_infarctionRemod_RL.h"

#COMPILE_MODEL_LIB "Paci_ToRORd_dynCl_PhiCaL_IKCa_mixed_apicobasal_infarctionRemod_RL" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$COMMON_HEADERS"
