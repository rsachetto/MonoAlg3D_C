############## TEN TUSCHER 2006 ##############################
MODEL_FILE_CPU="ten_tusscher_2006_RS_CPU.c"
MODEL_FILE_GPU="ten_tusscher_2006_RS_GPU.cu"
MODEL_FILE_SYCL=""
COMMON_HEADERS="ten_tusscher_2006.h"

COMPILE_MODEL_LIB "ten_tusscher_2006" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$MODEL_FILE_SYCL" "$COMMON_HEADERS"

##########################################################

############## TEN TUSCHER 2004 ##############################
MODEL_FILE_CPU="ten_tusscher_2004_RS_CPU.c"
MODEL_FILE_GPU="ten_tusscher_2004_RS_GPU.cu"
COMMON_HEADERS="ten_tusscher_2004.h"

COMPILE_MODEL_LIB "ten_tusscher_2004_endo" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$MODEL_FILE_SYCL" "$COMMON_HEADERS"
##########################################################

############## TEN TUSCHER 3 ENDO ##############################
MODEL_FILE_CPU="ten_tusscher_3_RS_CPU.c"
MODEL_FILE_GPU="ten_tusscher_3_RS_GPU.cu"
COMMON_HEADERS="ten_tusscher_3_RS.h"
COMPILE_MODEL_LIB "ten_tusscher_3_endo" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$MODEL_FILE_SYCL" "$COMMON_HEADERS" "-DENDO"
##########################################################

############## TEN TUSCHER 3 EPI ##############################
MODEL_FILE_CPU="ten_tusscher_3_RS_CPU.c"
MODEL_FILE_GPU="ten_tusscher_3_RS_GPU.cu"
COMMON_HEADERS="ten_tusscher_3_RS.h"
COMPILE_MODEL_LIB "ten_tusscher_3_epi" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$MODEL_FILE_SYCL" "$COMMON_HEADERS" "-DEPI"
##########################################################

############## ELNAZ TENTUSSCHER ORIGINAL 2004 ##############################
MODEL_FILE_CPU="ten_tusscher_2004_RS_CPU_epi.c"
MODEL_FILE_GPU="ten_tusscher_2004_RS_GPU_epi.cu"
COMMON_HEADERS="ten_tusscher_2004_epi.h"

COMPILE_MODEL_LIB "ten_tusscher_2004_epi" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$MODEL_FILE_SYCL" "$COMMON_HEADERS"
############################################################################################

############## TENTUSSCHER MIXED MYO EPI #####################################################
MODEL_FILE_CPU="mixed_tentusscher_myo_epi_2004.c"
MODEL_FILE_GPU="mixed_tentusscher_myo_epi_2004.cu"
COMMON_HEADERS="mixed_tentusscher_myo_epi_2004.h"

COMPILE_MODEL_LIB "mixed_tentusscher_myo_epi_2004" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$MODEL_FILE_SYCL" "$COMMON_HEADERS"
############################################################################################

############## TENTUSSCHER MIXED ENDO MID EPI #####################################################
MODEL_FILE_CPU="ten_tusscher_2004_mixed_endo_mid_epi.c"
MODEL_FILE_GPU="ten_tusscher_2004_mixed_endo_mid_epi.cu"
COMMON_HEADERS="ten_tusscher_2004_mixed_endo_mid_epi.h"

COMPILE_MODEL_LIB "ten_tusscher_2004_mixed_endo_mid_epi" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$MODEL_FILE_SYCL" "$COMMON_HEADERS"
############################################################################################

############## TENTUSSCHER 3 MIXED ENDO MID EPI #####################################################
MODEL_FILE_CPU="ten_tusscher_tt3_mixed_endo_mid_epi.c"
MODEL_FILE_GPU="ten_tusscher_tt3_mixed_endo_mid_epi.cu"
COMMON_HEADERS="ten_tusscher_tt3_mixed_endo_mid_epi.h"

COMPILE_MODEL_LIB "ten_tusscher_tt3_mixed_endo_mid_epi" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$MODEL_FILE_SYCL" "$COMMON_HEADERS"
############################################################################################

