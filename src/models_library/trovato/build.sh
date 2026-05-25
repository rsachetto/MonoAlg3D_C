############## TROVATO_2019 ##############################
MODEL_FILE_CPU="trovato_2019.c"
MODEL_FILE_GPU="trovato_2019.cu"
MODEL_FILE_SYCL=""
COMMON_HEADERS="trovato_2019.h"

COMPILE_MODEL_LIB "trovato_2019" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$MODEL_FILE_SYCL" "$COMMON_HEADERS"
##########################################################

############## TROVATO_2020 ##############################
MODEL_FILE_CPU="trovato_2020.c"
MODEL_FILE_GPU="trovato_2020.cu"
MODEL_FILE_SYCL=""
COMMON_HEADERS="trovato_2020.h"

COMPILE_MODEL_LIB "trovato_2020" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$MODEL_FILE_SYCL" "$COMMON_HEADERS"
##########################################################