############## Tworld ##############################
MODEL_FILE_CPU="Tworld.c"
MODEL_FILE_GPU="Tworld.cu"
COMMON_HEADERS="Tworld.h"

COMPILE_MODEL_LIB "Tworld" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$COMMON_HEADERS"
