SET(LIB_TYPE MODULE) #set the lib type

############## TEN TUSCHER 2006 ##############################
set(MODEL_FILE_CPU "ten_tusscher_2006_RS_CPU.c")
set(MODEL_FILE_GPU "ten_tusscher_2006_RS_GPU.cu")
set(COMMON_HEADERS model_common.h ten_tusscher_2006.h)

if(CUDA_FOUND)
    CUDA_ADD_LIBRARY(ten_tusscher_2006 ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
    ADD_LIBRARY(ten_tusscher_2006 ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)

#########################################################

############## TEN TUSCHER 2004 ##############################
set(MODEL_FILE_CPU "ten_tusscher_2004_RS_CPU.c")
set(MODEL_FILE_GPU "ten_tusscher_2004_RS_GPU.cu")
set(COMMON_HEADERS model_common.h ten_tusscher_2004.h)

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(ten_tusscher_2004_endo ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(ten_tusscher_2004_endo ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)

#########################################################

############## BONDARENKO ##############################
set(MODEL_FILE_CPU "bondarenko_2004.c")
set(MODEL_FILE_GPU "bondarenko_2004_GPU.cu")
set(COMMON_HEADERS model_common.h bondarenko_2004.h)


if(CUDA_FOUND)
    CUDA_ADD_LIBRARY(bondarenko_2004 ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
    ADD_LIBRARY(bondarenko_2004 ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#########################################################

############## TEN TUSCHER 3 ENDO ##############################
set(MODEL_FILE_CPU "ten_tusscher_3_ENDO_RS_CPU.c")
set(MODEL_FILE_GPU "ten_tusscher_3_ENDO_RS_GPU.cu")

if(CUDA_FOUND)
    CUDA_ADD_LIBRARY(ten_tusscher_3_endo ${LIB_TYPE} model_gpu_utils.h model_common.h ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
    ADD_LIBRARY(ten_tusscher_3_endo ${LIB_TYPE} model_common.h ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#########################################################

############## TEN TUSCHER 3 EPI ##############################
set(MODEL_FILE_CPU "ten_tusscher_3_EPI_RS_CPU.c")
set(MODEL_FILE_GPU "ten_tusscher_3_EPI_RS_GPU.cu")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(ten_tusscher_3_epi ${LIB_TYPE} model_gpu_utils.h model_common.h ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(ten_tusscher_3_epi ${LIB_TYPE} model_common.h ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#########################################################

############## TEN TUSCHER 3 FIG 4B ##############################
set(MODEL_FILE_CPU "ten_tusscher_3_Fig4b.c")
set(MODEL_FILE_GPU "ten_tusscher_3_Fig_4b.cu")
set (LIB_NAME "ten_tusscher_3_fig4b")

if(CUDA_FOUND)
    CUDA_ADD_LIBRARY( ${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h model_common.h ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
    ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_common.h ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#########################################################

############## FABERRUDY 2000 ##############################
set(MODEL_FILE_CPU "luo_rudy_1991.c")
set(MODEL_FILE_GPU "luo_rudy_1991.cu")
set(COMMON_HEADERS model_common.h luo_rudy_1991.h)

if(CUDA_FOUND)
    CUDA_ADD_LIBRARY(luo_rudy_1991 ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
    ADD_LIBRARY(luo_rudy_1991 ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#########################################################

############## OHARA_RUDY 2000 ##############################
set(MODEL_FILE_CPU "Ohara_Rudy_2011.c")
set(MODEL_FILE_GPU "Ohara_Rudy_2011.cu")
set(COMMON_HEADERS model_common.h Ohara_Rudy_2011.h)

set(LIB_NAME "ohara_rudy_endo_2011")

if(CUDA_FOUND)
    CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
    ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#########################################################

############## FHN 1961 MOD ##############################
set(MODEL_FILE_CPU "fhn_mod.c")
set(MODEL_FILE_GPU "fhn_mod.cu")
set(COMMON_HEADERS model_common.h fhn_mod.h)

set(LIB_NAME "fhn_mod")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#########################################################

############## MITCHELL SHAEFFER 2003 ##############################
set(MODEL_FILE_CPU "mitchell_shaeffer_2003.c")
set(MODEL_FILE_GPU "mitchell_shaeffer_2003.cu")
set(COMMON_HEADERS model_common.h mitchell_shaeffer_2003.h)

set(LIB_NAME "mitchell_shaeffer_2003")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#########################################################

############## MIXED FHN + MITCHELL ##############################
set(MODEL_FILE_CPU "mixed_fhn_mod_mitchell.c")
set(MODEL_FILE_GPU "mixed_fhn_mod_mitchell.cu")
set(COMMON_HEADERS model_common.h mixed_fhn_mod_mitchell.h)

set(LIB_NAME "mixed_fhn_mod_mitchell")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#########################################################

############## MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM ##############################
set(MODEL_FILE_CPU "mixed_tentusscher_myo_epi_2004.c")
set(MODEL_FILE_GPU "mixed_tentusscher_myo_epi_2004.cu")
set(COMMON_HEADERS model_common.h mixed_tentusscher_myo_epi_2004.h)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
############################################################################################

############## ELNAZ TENTUSSCHER ORIGINAL 2004 ##############################
set(MODEL_FILE_CPU "ten_tusscher_2004_RS_CPU_epi.c")
set(MODEL_FILE_GPU "ten_tusscher_2004_RS_GPU_epi.cu")
set(COMMON_HEADERS model_common.h ten_tusscher_2004_epi.h)

set(LIB_NAME "ten_tusscher_2004_epi")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
############################################################################################

############## ELNAZ TENTUSSCHER 2004 FOR SCENARIO 1 ##############################
set(MODEL_FILE_CPU "ten_tusscher_2004_RS_CPU_epi_S1_tr6_pop69.c")
set(MODEL_FILE_GPU "ten_tusscher_2004_RS_GPU_epi_S1_tr6_pop69.cu")
set(COMMON_HEADERS model_common.h ten_tusscher_2004_epi_S1_tr6_pop69.h)

set(LIB_NAME "ten_tusscher_2004_epi_S1_tr6_pop69")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
############################################################################################

############## ELNAZ TENTUSSCHER 2004 FOR SCENARIO 2 ##############################
set(MODEL_FILE_CPU "ten_tusscher_2004_RS_CPU_epi_S2_tr4_pop53.c")
set(MODEL_FILE_GPU "ten_tusscher_2004_RS_GPU_epi_S2_tr4_pop53.cu")
set(COMMON_HEADERS model_common.h ten_tusscher_2004_epi_S2_tr4_pop53.h)

set(LIB_NAME "ten_tusscher_2004_epi_S2_tr4_pop53")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
############################################################################################

############## ELNAZ TENTUSSCHER 2004 FOR SCENARIO 3 ##############################
set(MODEL_FILE_CPU "ten_tusscher_2004_RS_CPU_epi_S3_tr2_pop41.c")
set(MODEL_FILE_GPU "ten_tusscher_2004_RS_GPU_epi_S3_tr2_pop41.cu")
set(COMMON_HEADERS model_common.h ten_tusscher_2004_epi_S3_tr2_pop41.h)

set(LIB_NAME "ten_tusscher_2004_epi_S3_tr2_pop41")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
############################################################################################

############## NOBLE 1962 ##############################
set(MODEL_FILE_CPU "noble_1962.c")
set(MODEL_FILE_GPU "noble_1962.cu")
set(COMMON_HEADERS model_common.h noble_1962.h)

if(CUDA_FOUND)
    CUDA_ADD_LIBRARY(noble_1962 ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
        IF(MSVC)
                set(CUDA_LIBRARIES PUBLIC ${CUDA_LIBRARIES})
                TARGET_LINK_LIBRARIES(noble_1962 PUBLIC utils)
        endif()
else()
    ADD_LIBRARY(noble_1962 ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#########################################################

############## DI & FRANCESCO 1985 ##############################
set(MODEL_FILE_CPU "difrancesco.c")
set(MODEL_FILE_GPU "difrancesco.cu")
set(COMMON_HEADERS model_common.h difrancesco.h)

if(CUDA_FOUND)
    CUDA_ADD_LIBRARY(difrancesco ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
        IF(MSVC)
                set(CUDA_LIBRARIES PUBLIC ${CUDA_LIBRARIES})
                TARGET_LINK_LIBRARIES(difrancesco PUBLIC utils)
        endif()
else()
    ADD_LIBRARY(difrancesco ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#########################################################

############## ARPF 2009 ##############################
set(MODEL_FILE_CPU "stewart_aslanidi_noble_2009.c")
set(MODEL_FILE_GPU "stewart_aslanidi_noble_2009.cu")
set(COMMON_HEADERS model_common.h stewart_aslanidi_noble_2009.h)

if(CUDA_FOUND)
    CUDA_ADD_LIBRARY(stewart_aslanidi_noble_2009 ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
        IF(MSVC)
                set(CUDA_LIBRARIES PUBLIC ${CUDA_LIBRARIES})
                TARGET_LINK_LIBRARIES(stewart_aslanidi_noble_2009 PUBLIC utils)
        endif()
else()
    ADD_LIBRARY(stewart_aslanidi_noble_2009 ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#########################################################

############## LIRUDY 2011 ##############################
set(MODEL_FILE_CPU "li_rudy_2011.c")
set(MODEL_FILE_GPU "li_rudy_2011.cu")
set(COMMON_HEADERS model_common.h li_rudy_2011.h)

if(CUDA_FOUND)
    CUDA_ADD_LIBRARY(li_rudy_2011 ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
        IF(MSVC)
                set(CUDA_LIBRARIES PUBLIC ${CUDA_LIBRARIES})
                TARGET_LINK_LIBRARIES(li_rudy_2011 PUBLIC utils)
        endif()
else()
    ADD_LIBRARY(li_rudy_2011 ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#########################################################
