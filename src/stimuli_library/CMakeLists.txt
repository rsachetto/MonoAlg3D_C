SET(LIB_TYPE SHARED) #set the lib type
SET(COMMON_SRC ../libraries_common/config_helpers.c ../libraries_common/config_helpers.c)



############## DEFAULT STIM LIBRARY ##############################
ADD_LIBRARY(default_stimuli ${LIB_TYPE} ${COMMON_SRC} stimuli.c)

IF(MSVC)
	TARGET_LINK_LIBRARIES(default_stimuli hashes utils)
ENDIF()

#########################################################

