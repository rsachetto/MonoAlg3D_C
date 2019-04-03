SET(LIB_TYPE MODULE) #set the lib type

############## TEN TUSCHER 2006 ##############################
set(MODEL_FILE_CPU "ten_tusscher_2006_RS_CPU.c")
set(MODEL_FILE_GPU "ten_tusscher_2006_RS_GPU.cu")
set(COMMON_HEADERS model_common.h ten_tusscher_2006.h)

if(CUDA_FOUND)
    CUDA_ADD_LIBRARY(ten_tusscher_2006 ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
	IF(MSVC)
		set(CUDA_LIBRARIES PUBLIC ${CUDA_LIBRARIES})
		TARGET_LINK_LIBRARIES(ten_tusscher_2006 PUBLIC utils)	
	ENDIF()
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
	IF(MSVC)
		set(CUDA_LIBRARIES PUBLIC ${CUDA_LIBRARIES})
		TARGET_LINK_LIBRARIES(ten_tusscher_2004_endo PUBLIC utils)
	ENDIF()
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
	IF(MSVC)
		set(CUDA_LIBRARIES PUBLIC ${CUDA_LIBRARIES})
		TARGET_LINK_LIBRARIES(bondarenko_2004 PUBLIC utils)
	ENDIF()
else()
    ADD_LIBRARY(bondarenko_2004 ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#########################################################

############## TEN TUSCHER 3 ##############################
set(MODEL_FILE_CPU "ten_tusscher_3_RS_CPU.c")
set(MODEL_FILE_GPU "ten_tusscher_3_RS_GPU.cu")

if(CUDA_FOUND)
    CUDA_ADD_LIBRARY(ten_tusscher_3 ${LIB_TYPE} model_gpu_utils.h model_common.h ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
	IF(MSVC)
		set(CUDA_LIBRARIES PUBLIC ${CUDA_LIBRARIES})
		TARGET_LINK_LIBRARIES(ten_tusscher_3 PUBLIC utils)
	endif()
else()
    ADD_LIBRARY(ten_tusscher_3 ${LIB_TYPE} model_common.h ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#########################################################

############## FABERRUDY 2000 ##############################
set(MODEL_FILE_CPU "luo_rudy_1991.c")
set(MODEL_FILE_GPU "luo_rudy_1991.cu")
set(COMMON_HEADERS model_common.h luo_rudy_1991.h)

if(CUDA_FOUND)
    CUDA_ADD_LIBRARY(luo_rudy_1991 ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
	IF(MSVC)
		set(CUDA_LIBRARIES PUBLIC ${CUDA_LIBRARIES})
		TARGET_LINK_LIBRARIES(luo_rudy_1991 PUBLIC utils)
	endif()
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
    IF(MSVC)
        set(CUDA_LIBRARIES PUBLIC ${CUDA_LIBRARIES})
        TARGET_LINK_LIBRARIES(${LIB_NAME} PUBLIC utils)
    endif()
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
	IF(MSVC)
		set(CUDA_LIBRARIES PUBLIC ${CUDA_LIBRARIES})
		TARGET_LINK_LIBRARIES(${LIB_NAME} PUBLIC utils)
	endif()
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
	IF(MSVC)
		set(CUDA_LIBRARIES PUBLIC ${CUDA_LIBRARIES})
		TARGET_LINK_LIBRARIES(${LIB_NAME} PUBLIC utils)
	endif()
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#########################################################

