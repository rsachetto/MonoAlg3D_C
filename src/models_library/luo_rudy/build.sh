############### FABERRUDY 2000 ##############################
MODEL_FILE_CPU="luo_rudy_1991.c"
MODEL_FILE_GPU="luo_rudy_1991.cu"
COMMON_HEADERS="luo_rudy_1991.h"

COMPILE_MODEL_LIB "luo_rudy_1991" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$COMMON_HEADERS"
##########################################################