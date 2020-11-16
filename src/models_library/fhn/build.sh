############## FHN 1961 MOD ##############################
MODEL_FILE_CPU="fhn_mod.c"
MODEL_FILE_GPU="fhn_mod.cu"
COMMON_HEADERS="fhn_mod.h"

COMPILE_MODEL_LIB "fhn_mod" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$COMMON_HEADERS"
########################################################