############## ARPF 2009 ##############################
MODEL_FILE_CPU="stewart_aslanidi_noble_2009.c"
MODEL_FILE_GPU="stewart_aslanidi_noble_2009.cu"
MODEL_FILE_SYCL=""
COMMON_HEADERS="stewart_aslanidi_noble_2009.h"

COMPILE_MODEL_LIB "stewart_aslanidi_noble_2009" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$MODEL_FILE_SYCL" "$COMMON_HEADERS"
#########################################################