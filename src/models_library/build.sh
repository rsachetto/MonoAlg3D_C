COMPILE_MODEL_LIB () {
    local LIB_NAME=$1
    local MODEL_FILE_CPU=$2
    local MODEL_FILE_GPU=$3
    local COMMON_HEADERS="model_common.h $4"
    local EXTRA_C_FLAGS=$5
    local EXTRA_STATIC_LIBS=$6

    local MODEL_SOURCES="$MODEL_FILE_CPU"    

    local MODELS_STATIC_DEPS="config_helpers utils tinyexpr $EXTRA_STATIC_LIBS"
 
    if [ -n "$CUDA_FOUND" ]; then
        MODELS_EXTRA_LIB_PATH=$CUDA_LIBRARY_PATH
        MODELS_DYNAMIC_LIBS="c cuda cudart"

        MODEL_SOURCES="$MODEL_SOURCES $MODEL_FILE_GPU"
        COMMON_HEADERS="$COMMON_HEADERS model_gpu_utils.h"
    fi
    
    COMPILE_SHARED_LIB "$LIB_NAME" "$MODEL_SOURCES" "$COMMON_HEADERS" "$MODELS_STATIC_DEPS" "$MODELS_DYNAMIC_LIBS" "$MODELS_EXTRA_LIB_PATH" "$EXTRA_C_FLAGS" "$CUDA_FOUND"

}

############## TEN TUSCHER 2006 ##############################
MODEL_FILE_CPU="ten_tusscher_2006_RS_CPU.c"
MODEL_FILE_GPU="ten_tusscher_2006_RS_GPU.cu"
COMMON_HEADERS="ten_tusscher_2006.h"

COMPILE_MODEL_LIB "ten_tusscher_2006" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$COMMON_HEADERS"

##########################################################

############## TEN TUSCHER 2004 ##############################
MODEL_FILE_CPU="ten_tusscher_2004_RS_CPU.c"
MODEL_FILE_GPU="ten_tusscher_2004_RS_GPU.cu"
COMMON_HEADERS="ten_tusscher_2004.h"

COMPILE_MODEL_LIB "ten_tusscher_2004_endo" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$COMMON_HEADERS"
##########################################################

############## TEN TUSCHER 3 ENDO ##############################
MODEL_FILE_CPU="ten_tusscher_3_RS_CPU.c"
MODEL_FILE_GPU="ten_tusscher_3_RS_GPU.cu"
COMMON_HEADERS="ten_tusscher_3_RS.h"
COMPILE_MODEL_LIB "ten_tusscher_3_endo" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$COMMON_HEADERS" "-DENDO"
##########################################################

############## TEN TUSCHER 3 EPI ##############################
MODEL_FILE_CPU="ten_tusscher_3_RS_CPU.c"
MODEL_FILE_GPU="ten_tusscher_3_RS_GPU.cu"
COMMON_HEADERS="ten_tusscher_3_RS.h"
COMPILE_MODEL_LIB "ten_tusscher_3_epi" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$COMMON_HEADERS" "-DEPI"
##########################################################

############### BONDARENKO ##############################
MODEL_FILE_CPU="bondarenko_2004.c"
MODEL_FILE_GPU="bondarenko_2004_GPU.cu"
COMMON_HEADERS="bondarenko_2004.h"

COMPILE_MODEL_LIB "bondarenko_2004" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$COMMON_HEADERS" "-I$CUDA)"
##########################################################

############### FABERRUDY 2000 ##############################
MODEL_FILE_CPU="luo_rudy_1991.c"
MODEL_FILE_GPU="luo_rudy_1991.cu"
COMMON_HEADERS="luo_rudy_1991.h"

COMPILE_MODEL_LIB "luo_rudy_1991" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$COMMON_HEADERS"
##########################################################

############### OHARA_RUDY 2000 ##############################
MODEL_FILE_CPU="Ohara_Rudy_2011.c"
MODEL_FILE_GPU="Ohara_Rudy_2011.cu"
COMMON_HEADERS="Ohara_Rudy_2011.h"

COMPILE_MODEL_LIB "ohara_rudy_endo_2011" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$COMMON_HEADERS"

##########################################################

############## FHN 1961 MOD ##############################
MODEL_FILE_CPU="fhn_mod.c"
MODEL_FILE_GPU="fhn_mod.cu"
COMMON_HEADERS="fhn_mod.h"

COMPILE_MODEL_LIB "fhn_mod" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$COMMON_HEADERS"
#########################################################

############## MITCHELL SHAEFFER 2003 ##############################
MODEL_FILE_CPU="mitchell_shaeffer_2003.c"
MODEL_FILE_GPU="mitchell_shaeffer_2003.cu"
COMMON_HEADERS="mitchell_shaeffer_2003.h"

COMPILE_MODEL_LIB "mitchell_shaeffer_2003" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$COMMON_HEADERS"

#########################################################

############## ELNAZ TENTUSSCHER ORIGINAL 2004 ##############################
MODEL_FILE_CPU="ten_tusscher_2004_RS_CPU_epi.c"
MODEL_FILE_GPU="ten_tusscher_2004_RS_GPU_epi.cu"
COMMON_HEADERS="ten_tusscher_2004_epi.h"

COMPILE_MODEL_LIB "ten_tusscher_2004_epi" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$COMMON_HEADERS"
############################################################################################

############## NOBLE 1962 ##############################
MODEL_FILE_CPU="noble_1962.c"
MODEL_FILE_GPU="noble_1962.cu"
COMMON_HEADERS="noble_1962.h"

COMPILE_MODEL_LIB "noble_1962" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$COMMON_HEADERS"
#########################################################

############## DI & FRANCESCO 1985 ##############################
MODEL_FILE_CPU="difrancesco.c"
MODEL_FILE_GPU="difrancesco.cu"
COMMON_HEADERS="difrancesco.h"

COMPILE_MODEL_LIB "difrancesco" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$COMMON_HEADERS"
#########################################################

############## ARPF 2009 ##############################
MODEL_FILE_CPU="stewart_aslanidi_noble_2009.c"
MODEL_FILE_GPU="stewart_aslanidi_noble_2009.cu"
COMMON_HEADERS="stewart_aslanidi_noble_2009.h"

COMPILE_MODEL_LIB "stewart_aslanidi_noble_2009" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$COMMON_HEADERS"
#########################################################

############## LIRUDY 2011 ##############################
MODEL_FILE_CPU="li_rudy_2011.c"
MODEL_FILE_GPU="li_rudy_2011.cu"
COMMON_HEADERS="li_rudy_2011.h"

COMPILE_MODEL_LIB "li_rudy_2011" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$COMMON_HEADERS"

############## ToTOTd ##############################
MODEL_FILE_CPU="ToRORd_fkatp_endo.c"
MODEL_FILE_GPU="ToRORd_fkatp_endo.cu"
COMMON_HEADERS="ToRORd_fkatp_endo.h"

COMPILE_MODEL_LIB "ToRORd_fkatp_endo" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$COMMON_HEADERS"