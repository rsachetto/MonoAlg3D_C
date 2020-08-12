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

############## TEN TUSCHER 2004 ENDO ##############################
set(MODEL_FILE_CPU "ten_tusscher_2004_RS_CPU.c")
set(MODEL_FILE_GPU "ten_tusscher_2004_RS_GPU.cu")
set(COMMON_HEADERS model_common.h ten_tusscher_2004.h)

if(CUDA_FOUND)
	ADD_DEFINITIONS(-DENDO)
	CUDA_ADD_LIBRARY(ten_tusscher_2004_endo ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
	REMOVE_DEFINITIONS(-DENDO)
else()
	ADD_DEFINITIONS(-DENDO)
	ADD_LIBRARY(ten_tusscher_2004_endo ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
	REMOVE_DEFINITIONS(-DENDO)
endif(CUDA_FOUND)
#########################################################

############## TEN TUSCHER 2004 EPI ##############################
set(MODEL_FILE_CPU "ten_tusscher_2004_RS_CPU.c")
set(MODEL_FILE_GPU "ten_tusscher_2004_RS_GPU.cu")
set(COMMON_HEADERS model_common.h ten_tusscher_2004.h)

if(CUDA_FOUND)
	ADD_DEFINITIONS(-DEPI)
	CUDA_ADD_LIBRARY(ten_tusscher_2004_epi ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
	REMOVE_DEFINITIONS(-DEPI)
else()
	ADD_DEFINITIONS(-DEPI)
	ADD_LIBRARY(ten_tusscher_2004_epi ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
	REMOVE_DEFINITIONS(-DEPI)
endif(CUDA_FOUND)

#########################################################

############## TEN TUSCHER 2004 MCELL ##############################
set(MODEL_FILE_CPU "ten_tusscher_2004_RS_CPU.c")
set(MODEL_FILE_GPU "ten_tusscher_2004_RS_GPU.cu")
set(COMMON_HEADERS model_common.h ten_tusscher_2004.h)

if(CUDA_FOUND)
	ADD_DEFINITIONS(-DMCELL)
	CUDA_ADD_LIBRARY(ten_tusscher_2004_myo ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
	REMOVE_DEFINITIONS(-DMCELL)
else()
	ADD_DEFINITIONS(-DMCELL)
	ADD_LIBRARY(ten_tusscher_2004_myo ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
	REMOVE_DEFINITIONS(-DMCELL)
endif(CUDA_FOUND)

#########################################################

############## TEN TUSCHER 3 ENDO ##############################
set(MODEL_FILE_CPU "ten_tusscher_3_RS_CPU.c")
set(MODEL_FILE_GPU "ten_tusscher_3_RS_GPU.cu")
set(COMMON_HEADERS model_common.h ten_tusscher_3_RS.h)

if(CUDA_FOUND)
	ADD_DEFINITIONS(-DENDO)
	CUDA_ADD_LIBRARY(ten_tusscher_3_endo ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
	REMOVE_DEFINITIONS(-DENDO)
else()
	ADD_DEFINITIONS(-DENDO)
	ADD_LIBRARY(ten_tusscher_3_endo ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
	REMOVE_DEFINITIONS(-DENDO)
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

############## ToRORd (Working version) ##############################
set(MODEL_FILE_CPU "new_torord_fkatp_endo.c")
set(MODEL_FILE_GPU "new_torord_fkatp_endo.cu")
set(COMMON_HEADERS model_common.h new_torord_fkatp_endo.h)

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(ToRORd_fkatp_endo ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(ToRORd_fkatp_endo ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
################################################################################

############## CRN ##############################
set(MODEL_FILE_CPU "courtemanche_ramirez_nattel_1998.c")
set(MODEL_FILE_GPU "courtemanche_ramirez_nattel_1998.cu")
set(COMMON_HEADERS model_common.h courtemanche_ramirez_nattel_1998.h)

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(courtemanche_ramirez_nattel_1998 ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(courtemanche_ramirez_nattel_1998 ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
################################################################################

############## TROVATO ##############################
set(MODEL_FILE_CPU "trovato_2019.c")
set(MODEL_FILE_GPU "trovato_2019.cu")
set(COMMON_HEADERS model_common.h trovato_2019.h)

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(trovato_2019 ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(trovato_2019 ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
################################################################################


#########################################################################################################################################################
# MIXED MODELS SECTION
#########################################################################################################################################################

############## MIXED FHN + MITCHELL ##############################
set(MODEL_FILE_CPU "00_Mixed_Models/fhn_mitchell.c")
set(MODEL_FILE_GPU "00_Mixed_Models/fhn_mitchell.cu")
set(COMMON_HEADERS model_common.h 00_Mixed_Models/fhn_mitchell.h)

set(LIB_NAME "mixed_fhn_mitchell")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#########################################################

############## TT3: MCELL + EPI ##############################
set(MODEL_FILE_CPU "00_Mixed_Models/tt2004_myo_epi.c")
set(MODEL_FILE_GPU "00_Mixed_Models/tt2004_myo_epi.cu")
set(COMMON_HEADERS model_common.h 00_Mixed_Models/tt2004_myo_epi.h)

set(LIB_NAME "mixed_tt2004_myo_epi")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#########################################################

#########################################################################################################################################################
# PEDRO SECTION
#########################################################################################################################################################

#[[
############## PEDRO TENTUSSCHER ENDO SENSIBILITY ##############################
set(MODEL_FILE_CPU "ten_tusscher_3_sensibility.c")
set(MODEL_FILE_GPU "ten_tusscher_3_sensibility.cu")
set(COMMON_HEADERS model_common.h ten_tusscher_3_sensibility.h)

if(CUDA_FOUND)
    CUDA_ADD_LIBRARY(ten_tusscher_3_sensibility ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
        IF(MSVC)
                set(CUDA_LIBRARIES PUBLIC ${CUDA_LIBRARIES})
                TARGET_LINK_LIBRARIES(ten_tusscher_3_sensibility PUBLIC utils)
        endif()
else()
    ADD_LIBRARY(ten_tusscher_3_sensibility ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
################################################################################

############## PEDRO TENTUSSCHER TWO REGIONS ###################################
set(MODEL_FILE_CPU "ten_tusscher_3_two_region.c")
set(MODEL_FILE_GPU "ten_tusscher_3_two_region.cu")
set(COMMON_HEADERS model_common.h ten_tusscher_3_two_region.h)

if(CUDA_FOUND)
    CUDA_ADD_LIBRARY(ten_tusscher_3_two_region ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
        IF(MSVC)
                set(CUDA_LIBRARIES PUBLIC ${CUDA_LIBRARIES})
                TARGET_LINK_LIBRARIES(ten_tusscher_3_two_region PUBLIC utils)
        endif()
else()
    ADD_LIBRARY(ten_tusscher_3_two_region ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
################################################################################

############## PEDRO TENTUSSCHER ENDO PERCOLATION ##############################
set(MODEL_FILE_CPU "ten_tusscher_3_percolation_ENDO_CPU.c")
set(MODEL_FILE_GPU "ten_tusscher_3_percolation_ENDO_GPU.cu")
set(COMMON_HEADERS model_common.h ten_tusscher_3_percolation.h)

if(CUDA_FOUND)
    CUDA_ADD_LIBRARY(ten_tusscher_3_percolation ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
        IF(MSVC)
                set(CUDA_LIBRARIES PUBLIC ${CUDA_LIBRARIES})
                TARGET_LINK_LIBRARIES(ten_tusscher_3_percolation PUBLIC utils)
        endif()
else()
    ADD_LIBRARY(ten_tusscher_3_percolation ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
################################################################################
]]
#########################################################################################################################################################
# ELNAZ SECTION
#########################################################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 1 ##############################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S1.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S1.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S1.h)

set(LIB_NAME "tentusscher_epi_2004_S1")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 1 ##############################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S1.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S1.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S1.h)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S1")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 2.1 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_1.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_1.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_1.h)

set(LIB_NAME "tentusscher_epi_2004_S2_1")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 2.2 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_2.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_2.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_2.h)

set(LIB_NAME "tentusscher_epi_2004_S2_2")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 2.3 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_3.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_3.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_3.h)

set(LIB_NAME "tentusscher_epi_2004_S2_3")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 2.4 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_4.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_4.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_4.h)

set(LIB_NAME "tentusscher_epi_2004_S2_4")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 2.5 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_5.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_5.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_5.h)

set(LIB_NAME "tentusscher_epi_2004_S2_5")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 2.6 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_6.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_6.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_6.h)

set(LIB_NAME "tentusscher_epi_2004_S2_6")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 2.7 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_7.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_7.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_7.h)

set(LIB_NAME "tentusscher_epi_2004_S2_7")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 2.8 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_8.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_8.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_8.h)

set(LIB_NAME "tentusscher_epi_2004_S2_8")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 2.9 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_9.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_9.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_9.h)

set(LIB_NAME "tentusscher_epi_2004_S2_9")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 2.10 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_10.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_10.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_10.h)

set(LIB_NAME "tentusscher_epi_2004_S2_10")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 2.11 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_11.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_11.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_11.h)

set(LIB_NAME "tentusscher_epi_2004_S2_11")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 2.12 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_12.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_12.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_12.h)

set(LIB_NAME "tentusscher_epi_2004_S2_12")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 2.13 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_13.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_13.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_13.h)

set(LIB_NAME "tentusscher_epi_2004_S2_13")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 2.14 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_14.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_14.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_14.h)

set(LIB_NAME "tentusscher_epi_2004_S2_14")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 2.15 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_15.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_15.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_15.h)

set(LIB_NAME "tentusscher_epi_2004_S2_15")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 2.16 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_16.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_16.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_16.h)

set(LIB_NAME "tentusscher_epi_2004_S2_16")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 2.17 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_17.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_17.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_17.h)

set(LIB_NAME "tentusscher_epi_2004_S2_17")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 2.18 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_18.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_18.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_18.h)

set(LIB_NAME "tentusscher_epi_2004_S2_18")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 2.19 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_19.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_19.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_19.h)

set(LIB_NAME "tentusscher_epi_2004_S2_19")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 2.20 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_20.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_20.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S2_20.h)

set(LIB_NAME "tentusscher_epi_2004_S2_20")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 3.1 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_1.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_1.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_1.h)

set(LIB_NAME "tentusscher_epi_2004_S3_1")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 3.2 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_2.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_2.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_2.h)

set(LIB_NAME "tentusscher_epi_2004_S3_2")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 3.3 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_3.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_3.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_3.h)

set(LIB_NAME "tentusscher_epi_2004_S3_3")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 3.4 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_4.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_4.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_4.h)

set(LIB_NAME "tentusscher_epi_2004_S3_4")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 3.5 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_5.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_5.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_5.h)

set(LIB_NAME "tentusscher_epi_2004_S3_5")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 3.6 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_6.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_6.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_6.h)

set(LIB_NAME "tentusscher_epi_2004_S3_6")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 3.7 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_7.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_7.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_7.h)

set(LIB_NAME "tentusscher_epi_2004_S3_7")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 3.8 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_8.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_8.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_8.h)

set(LIB_NAME "tentusscher_epi_2004_S3_8")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 3.9 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_9.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_9.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_9.h)

set(LIB_NAME "tentusscher_epi_2004_S3_9")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 3.10 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_10.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_10.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_10.h)

set(LIB_NAME "tentusscher_epi_2004_S3_10")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 3.11 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_11.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_11.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_11.h)

set(LIB_NAME "tentusscher_epi_2004_S3_11")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 3.12 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_12.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_12.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_12.h)

set(LIB_NAME "tentusscher_epi_2004_S3_12")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 3.13 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_13.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_13.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_13.h)

set(LIB_NAME "tentusscher_epi_2004_S3_13")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 3.14 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_14.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_14.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_14.h)

set(LIB_NAME "tentusscher_epi_2004_S3_14")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 3.15 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_15.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_15.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_15.h)

set(LIB_NAME "tentusscher_epi_2004_S3_15")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 3.16 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_16.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_16.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_16.h)

set(LIB_NAME "tentusscher_epi_2004_S3_16")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 3.17 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_17.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_17.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_17.h)

set(LIB_NAME "tentusscher_epi_2004_S3_17")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 3.18 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_18.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_18.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_18.h)

set(LIB_NAME "tentusscher_epi_2004_S3_18")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 3.19 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_19.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_19.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_19.h)

set(LIB_NAME "tentusscher_epi_2004_S3_19")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ TENTUSSCHER 2004 EPICARDIUM FOR SCENARIO 3.20 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_20.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_20.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/01_Single_Models/ten_tusscher_2004_epi_S3_20.h)

set(LIB_NAME "tentusscher_epi_2004_S3_20")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 2.1 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_1.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_1.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_1.h)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S2_1")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 2.2 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_2.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_2.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_2)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S2_2")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 2.3 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_3.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_3.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_3)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S2_3")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 2.4 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_4.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_4.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_4)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S2_4")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 2.5 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_5.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_5.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_5)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S2_5")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 2.6 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_6.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_6.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_6)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S2_6")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 2.7 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_7.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_7.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_7)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S2_7")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 2.8 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_8.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_8.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_8)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S2_8")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 2.9 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_9.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_9.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_9)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S2_9")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 2.10 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_10.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_10.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_10)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S2_10")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 2.11 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_11.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_11.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_11)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S2_11")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 2.12 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_12.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_12.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_12)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S2_12")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 2.13 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_13.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_13.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_13)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S2_13")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 2.14 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_14.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_14.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_14)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S2_14")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 2.15 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_15.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_15.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_15)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S2_15")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 2.16 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_16.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_16.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_16)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S2_16")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 2.17 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_17.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_17.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_17)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S2_17")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 2.18 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_18.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_18.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_18)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S2_18")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 2.19 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_19.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_19.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_19)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S2_19")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 2.20 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_20.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_20.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S2_20)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S2_20")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 3.1 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_1.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_1.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_1)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S3_1")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 3.2 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_2.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_2.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_2)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S3_2")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 3.3 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_3.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_3.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_3)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S3_3")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 3.4 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_4.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_4.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_4)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S3_4")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 3.5 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_5.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_5.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_5)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S3_5")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 3.6 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_6.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_6.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_6)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S3_6")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 3.7 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_7.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_7.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_7)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S3_7")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 3.8 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_8.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_8.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_8)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S3_8")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 3.9 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_9.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_9.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_9)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S3_9")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 3.10 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_10.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_10.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_10)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S3_10")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 3.11 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_11.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_11.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_11)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S3_11")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 3.12 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_12.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_12.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_12)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S3_12")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 3.13 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_13.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_13.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_13)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S3_13")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 3.14 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_14.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_14.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_14)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S3_14")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 3.15 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_15.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_15.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_15)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S3_15")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 3.16 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_16.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_16.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_16)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S3_16")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 3.17 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_17.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_17.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_17)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S3_17")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 3.18 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_18.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_18.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_18)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S3_18")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 3.19 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_19.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_19.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_19)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S3_19")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################

############## ELNAZ MIXED TENTUSSCHER 2004 MYOCARDIUM + EPICARDIUM FOR SCENARIO 3.20 ##########################################
set(MODEL_FILE_CPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_20.c")
set(MODEL_FILE_GPU "01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_20.cu")
set(COMMON_HEADERS model_common.h 01_Elnaz_Models/02_Mixed_Models/mixed_tentusscher_myo_epi_2004_S3_20)

set(LIB_NAME "mixed_tentusscher_myo_epi_2004_S3_20")

if(CUDA_FOUND)
	CUDA_ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} model_gpu_utils.h ${COMMON_HEADERS} ${MODEL_FILE_GPU} ${MODEL_FILE_CPU})
else()
	ADD_LIBRARY(${LIB_NAME} ${LIB_TYPE} ${COMMON_HEADERS} ${MODEL_FILE_CPU})
endif(CUDA_FOUND)
#################################################################################################################
