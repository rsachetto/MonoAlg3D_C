SET(LIB_TYPE MODULE) #set the lib type
SET(COMMON_SRC ../libraries_common/config_helpers.c ../libraries_common/config_helpers.c)

############## DEFAULT STIM LIBRARY 2006 ##############################
ADD_LIBRARY(default_domains ${LIB_TYPE} ${COMMON_SRC} grid_domain.c domain_helpers.c domain_helpers.h)

IF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
	TARGET_LINK_LIBRARIES(default_domains hashes utils alg)
ENDIF()

#########################################################

