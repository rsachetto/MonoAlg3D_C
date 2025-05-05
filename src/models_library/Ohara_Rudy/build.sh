############### OHARA_RUDY 2000 ##############################
MODEL_FILE_CPU="Ohara_Rudy_2011.c"
MODEL_FILE_GPU="Ohara_Rudy_2011.cu"
MODEL_FILE_SYCL=""
COMMON_HEADERS="Ohara_Rudy_2011.h"

COMPILE_MODEL_LIB "ohara_rudy_endo_2011" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$MODEL_FILE_SYCL" "$COMMON_HEADERS"
##########################################################