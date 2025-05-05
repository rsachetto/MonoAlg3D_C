############## Maleckar 2008 ##########################################################
MODEL_FILE_CPU="Maleckar2008.c"
MODEL_FILE_GPU="Maleckar2008.cu"
MODEL_FILE_SYCL=""
COMMON_HEADERS="Maleckar2008.h"

COMPILE_MODEL_LIB "Maleckar2008" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$MODEL_FILE_SYCL" "$COMMON_HEADERS"
######################################################################################

############## Maleckar 2009 ##########################################################
MODEL_FILE_CPU="Maleckar2009.c"
MODEL_FILE_GPU="Maleckar2009.cu"
MODEL_FILE_SYCL=""
COMMON_HEADERS="Maleckar2009.h"

COMPILE_MODEL_LIB "Maleckar2009" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$MODEL_FILE_SYCL" "$COMMON_HEADERS"
######################################################################################
